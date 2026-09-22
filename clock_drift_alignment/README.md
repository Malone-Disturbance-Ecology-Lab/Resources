# Temperature-dependent clock-drift alignment

`align_clock_drift_temperature.R`

Aligns eddy-covariance **flux** data with **met/biomet** data when the two
streams are logged on separate instruments whose clocks drift apart — and where
the drift *rate* changes with temperature.

## The problem

Independent datalogger clocks use uncompensated quartz crystals whose frequency
follows a parabolic temperature response (`rate ≈ a·(T − T0)² + b`, turnover
`T0 ≈ 25 °C`). Because a field enclosure swings diurnally and seasonally, the
drift rate oscillates and the **accumulated** offset is the time-integral of a
temperature-dependent rate. A single linear (or piecewise-linear) offset fit
therefore leaves structured, temperature-correlated residuals.

This script estimates the offset from data, fits the *rate* as a function of the
logger's panel temperature, integrates it to a continuous `offset(t)`, and
reindexes one stream onto the other.

## Before you run it (Step 0)

Rule out convention artefacts first — they masquerade as drift:

- **Period label** — start / middle / end of the averaging interval. Declared
  per stream via `flux_label` / `met_label` and normalized to a common instant.
- **Timezone / DST** — this script assumes timestamps are in **UTC**
  (no DST); parse them accordingly (`ymd_hms(x, tz = "UTC")`).
- **Fixed integer-hour offsets** set at logger install.

## Method

1. `normalize_period_label()` — put both streams on a common (midpoint) instant.
2. `estimate_offsets()` — windowed cross-correlation of a shared signal on a fine
   grid (sub-averaging-interval lags via parabolic peak refinement) → anchor
   offsets `o_i(t_i)`, weighted by correlation. Windows are multi-day so the
   drift signal does not alias with the reference variable's diurnal cycle.
3. `fit_temperature_drift()` — fit `rate(T) = a·(T−T0)² + b` by integrating over
   the corrected logger's panel temperature (`PTemp`) between anchors.
4. `reconstruct_offset()` — reuse the exact fitted integral over the full `PTemp`
   record → continuous `offset(t)` (fills gaps in anchor coverage).
5. `align_streams()` — apply the offset and reindex onto the reference grid.
6. `plot_drift_diagnostics()` — offset timeline, rate-vs-temperature curve, and
   the **residual-vs-PTemp** panel (the key check: it should be flat).

## Usage

```r
source("align_clock_drift_temperature.R")

flux <- read.csv("flux_30min.csv");  met <- read.csv("biomet_30min.csv")
flux$TIMESTAMP <- ymd_hms(flux$TIMESTAMP, tz = "UTC")
met$TIMESTAMP  <- ymd_hms(met$TIMESTAMP,  tz = "UTC")

res <- align_clock_drift(
  flux, met,
  reference = "solar_blend",   # see "Reference options" below
  correct   = "met",           # which stream to drift-correct
  # map YOUR dataframe's columns:
  flux_cols = flux_columns(time = "TIMESTAMP", shared = "Ts", sw = "SW_IN"),
  met_cols  = met_columns (time = "TIMESTAMP", shared = "Ta", sw = "SW_IN", ptemp = "PTemp"),
  flux_label = "end", met_label = "start", interval_min = 30,
  lat = 25.4, lon = -80.6,     # site coords (solar anchoring)
  window_days = 7, step_days = 2, min_corr = 0.7)

res$fit$par            # a, b, T0, s0
summary(res$residual_offsets)   # post-correction residuals (sec) ≈ 0 if good
head(res$merged)       # reference stream + drift-corrected "<var>_corr" columns
```

## Setting column names from your dataframe

Point the pipeline at whatever columns your data uses via two mapping helpers —
each a named list with keys `time`, `shared`, `sw`, `ptemp`:

```r
flux_cols = flux_columns(time = "TIMESTAMP", shared = "Ts", sw = "SW_IN")
met_cols  = met_columns (time = "TIMESTAMP", shared = "Ta", sw = "SW_IN", ptemp = "PTemp")
```

- `time` — timestamp column (POSIXct); `shared` — signal for cross-correlation
  (sonic vs. air temperature); `sw` — incoming shortwave (only for solar
  anchoring; set `NA` if unused); `ptemp` — logger panel temperature, **required
  on whichever stream is being corrected**.
- Every mapped name is validated against the dataframe you pass; a missing column
  raises an error that lists the columns actually present.
- To discover names automatically, call `suggest_columns(df)` — it matches common
  EC/AmeriFlux aliases (case-insensitive) and returns a mapping list for you to
  review and edit before use.

## Reference options (`reference =`)

You choose the trusted time base:

- `"flux"` — trust the flux clock; correct the other stream to it.
- `"met"` — trust the met clock; correct the other stream to it.
- `"solar"` — pin the corrected stream to true solar time via its `SW_in`
  vs. modelled potential radiation (absolute, not just relative).
- `"shared"` — relative alignment between the two streams (no absolute pin).
- `"custom"` — supply your own `ref_time` and `ref_val` (e.g. a GPS-synced
  signal).
- `"solar_blend"` — solar-primary with shared-signal gap-fill in windows where
  the solar anchor is weak (night-heavy, overcast). *Note:* the shared fallback
  assumes the other stream is ~on true time; solar windows are preferred.

## Outputs

`align_clock_drift()` returns a list: `anchors`, `fit` (parameters + rate
function + anchor residuals), `corrected` (the corrected stream with
`offset_sec` and `time_corr`), `merged` (reference frame + `<var>_corr`
columns), and `residual_offsets` (post-correction verification). Diagnostic
PNGs are written to `outdir` (default `drift_diag/`).

## Caveats

- You cannot resolve sub-averaging-interval drift from aggregated products; fit
  the seasonal term and don't chase the diurnal wiggle if its amplitude is small
  relative to the averaging interval.
- Don't confuse this with the sonic-vs-IRGA transport lag (handled by covariance
  maximization on high-frequency data) — that is a physical delay, not clock
  drift.
- The durable fix going forward is hardware: GPS/NTP time sync, or logging both
  streams on one clock.

## Dependencies

Base R + `lubridate`; `ggplot2` for the diagnostic plots only.

## Gap handling and fit checks

- `max_gap_sec = 3600` limits the time between interpolation endpoints for
  anchoring and output alignment. Explicit missing signal values remain missing;
  gaps from absent rows longer than this limit also remain missing. Set this to
  match your sampling interval and acceptable interpolation span.
- `max_ptemp_gap_sec = 7200` permits linear temperature interpolation only between
  valid readings at most two hours apart. Missing temperatures at either record
  edge, larger gaps, and duplicate timestamps stop the fit with an explanation.
  Supply temperature coverage or fit separate covered records. The fitted integral
  is stored and reused for reconstruction, so fitting and correction cannot apply
  different gap policies. Reconstruction is restricted to that fitted time range;
  refit if the temperature data change.
- Fitting requires at least five distinct valid anchors and a full-rank integrated
  temperature design. Constant temperature cannot identify the four parameters.
  Overlapping windows still produce correlated anchors; the count/rank checks do
  not establish statistical independence or quantify parameter uncertainty.
- A failed optimizer stops the pipeline before a correction is applied. A linear
  fit of the integrated polynomial supplies starting values for the nonlinear fit.
- Short records return an empty anchor table; the orchestrator explains when too
  few anchors are available. Solar blending handles empty solar or shared tables.
- `reference` and `correct` cannot name the same logger. For `reference="custom"`,
  supply UTC midpoint `ref_time`; verification uses that same custom reference.
  `merged` still uses the other logger's grid. Solar-blend verification remains
  solar-only and may have no residual anchors when shortwave data are unavailable.
- Default output variables are numeric columns. Explicitly selected nonnumeric
  columns raise an error; metadata are not interpolated automatically. Numeric
  quality flags should be excluded with `out_vars` if interpolation is unsuitable.

Run the synthetic regression checks from this directory:

```sh
Rscript test_clock_drift.R align_clock_drift_temperature.R
```

These checks do not replace validation against field data.
