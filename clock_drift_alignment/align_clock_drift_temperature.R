# =============================================================================
# Temperature-dependent clock-drift alignment for EC flux vs. met/biomet data
# -----------------------------------------------------------------------------
# Two independent logger clocks drift apart, and the drift RATE depends on
# logger temperature (uncompensated 32.768 kHz crystal: rate ~ a*(T - T0)^2 + b).
# The accumulated offset is the time-integral of that temperature-dependent rate,
# so a single linear fit leaves structured, temperature-correlated residuals.
#
# WHAT YOU CONTROL (set once in align_clock_drift(), the orchestrator):
#   * columns   : map YOUR dataframe's columns via flux_cols / met_cols
#                 (flux_columns()/met_columns(); suggest_columns() auto-detects)
#   * reference : which time base to trust  -- "flux" | "met" | "solar" |
#                 "shared" | "custom"  (you provide the reference)
#   * correct   : which stream gets corrected -- "met" | "flux"
#   * period labels : start / mid / end of averaging period, PER STREAM
#   * timestamps assumed UTC (no DST); solar geometry uses UTC directly.
#
# Pipeline:
#   0. normalize_period_label()  put both streams on a common (midpoint) instant
#   1. estimate_offsets()        windowed cross-correlation -> anchor offsets
#   2. fit_temperature_drift()   fit rate(T)=a*(T-T0)^2+b by integrating over the
#                                CORRECTED logger's panel temp (PTemp)
#   3. reconstruct_offset()      integrate fitted rate over full PTemp -> offset(t)
#   4. align_streams()           apply offset, reindex onto the reference grid
#   5. plot_drift_diagnostics()  offset(t), rate(T), residual-vs-T
#
# Convention: offset = seconds ADDED to the CORRECTED stream's timestamps to
# align them with the reference time base. Dependencies: base R + lubridate,
# plus ggplot2 for diagnostics only.
# =============================================================================

suppressPackageStartupMessages(library(lubridate))

# Bounded interpolation: preserve explicit missing observations and never bridge
# an outage longer than max_gap_sec. Exact valid observations remain available.
.gap_approx <- function(x, y, xout, max_gap_sec = 3600) {
  if (length(x) != length(y) || !is.numeric(y))
    stop("Interpolation requires equally sized timestamps and numeric values.")
  if (length(max_gap_sec) != 1L || !is.finite(max_gap_sec) || max_gap_sec <= 0)
    stop("max_gap_sec must be finite and positive.")
  x <- as.numeric(x); xo <- as.numeric(xout)
  keep <- is.finite(x); x <- x[keep]; y <- y[keep]
  ord <- order(x); x <- x[ord]; y <- y[ord]
  if (anyDuplicated(x)) stop("Duplicate timestamps must be resolved before interpolation.")
  y[!is.finite(y)] <- NA_real_
  ans <- rep(NA_real_, length(xo))
  if (!length(x)) return(ans)
  exact <- match(xo, x); hit <- !is.na(exact)
  ans[hit] <- y[exact[hit]]
  if (length(x) < 2L) return(ans)
  idx <- which(is.finite(xo) & !hit & xo > x[1] & xo < tail(x, 1))
  j <- findInterval(xo[idx], x)
  good <- x[j + 1L] - x[j] <= max_gap_sec
  idx <- idx[good]; j <- j[good]
  ans[idx] <- y[j] + (y[j + 1L] - y[j]) *
    (xo[idx] - x[j]) / (x[j + 1L] - x[j])
  ans
}

.empty_anchors <- function() data.frame(
  t_center = as.POSIXct(numeric(), origin = "1970-01-01", tz = "UTC"),
  offset = numeric(), corr = numeric(), n = integer())

# -----------------------------------------------------------------------------
# 0a. Normalize the averaging-period label to a common instant (default: the
#     period MIDPOINT). Run this on BOTH streams before anchoring, so a fixed
#     start-vs-end labelling difference is not mistaken for clock drift.
#     from: "start" | "mid" | "end";  interval_min: averaging interval (min).
# -----------------------------------------------------------------------------
normalize_period_label <- function(time, from = c("end", "start", "mid"),
                                   interval_min = 30, to = c("mid", "start", "end")) {
  from <- match.arg(from); to <- match.arg(to)
  half <- interval_min * 60 / 2
  full <- interval_min * 60
  to_mid   <- switch(from, start =  half, mid = 0, end = -half)  # shift label -> midpoint
  mid_to   <- switch(to,   start = -half, mid = 0, end =  half)  # midpoint -> target label
  time + to_mid + mid_to
}

# -----------------------------------------------------------------------------
# 0b. Modelled potential (extraterrestrial) radiation on a horizontal surface,
#     for solar-time anchoring. datetime_utc: POSIXct in UTC. lon EAST-positive.
#     Returns W m-2 top-of-atmosphere; its TIMING is what the anchoring uses.
# -----------------------------------------------------------------------------
potential_radiation <- function(datetime_utc, lat, lon, S0 = 1361) {
  if (length(lat) != 1L || length(lon) != 1L ||
      !is.finite(lat) || !is.finite(lon) || abs(lat) > 90 || abs(lon) > 180)
    stop("Solar anchoring requires finite lat (-90..90) and lon (-180..180).")
  lt  <- as.POSIXlt(datetime_utc, tz = "UTC")
  doy <- lt$yday + 1L
  hr  <- lt$hour + lt$min / 60 + lt$sec / 3600
  g   <- 2 * pi / 365 * (doy - 1 + (hr - 12) / 24)
  decl <- 0.006918 - 0.399912 * cos(g) + 0.070257 * sin(g) -
          0.006758 * cos(2 * g) + 0.000907 * sin(2 * g) -
          0.002697 * cos(3 * g) + 0.001480 * sin(3 * g)
  eqt  <- 229.18 * (0.000075 + 0.001868 * cos(g) - 0.032077 * sin(g) -
          0.014615 * cos(2 * g) - 0.040849 * sin(2 * g))
  E0   <- 1.00011 + 0.034221 * cos(g) + 0.001280 * sin(g) +
          0.000719 * cos(2 * g) + 0.000077 * sin(2 * g)
  tst  <- hr * 60 + eqt + 4 * lon                 # true solar time (min), UTC
  ha   <- (tst / 4 - 180) * pi / 180
  latr <- lat * pi / 180
  cosz <- sin(latr) * sin(decl) + cos(latr) * cos(decl) * cos(ha)
  pmax(S0 * E0 * cosz, 0)
}

# -----------------------------------------------------------------------------
# 1a. Sub-bin peak refinement of a correlation curve (parabolic interpolation).
# -----------------------------------------------------------------------------
.refine_peak <- function(shifts, corr) {
  j <- which.max(corr)
  if (j == 1L || j == length(corr)) return(shifts[j])
  cm1 <- corr[j - 1]; c0 <- corr[j]; cp1 <- corr[j + 1]
  denom <- (cm1 - 2 * c0 + cp1)
  if (!is.finite(denom) || denom == 0) return(shifts[j])
  delta <- max(min(0.5 * (cm1 - cp1) / denom, 1), -1)
  shifts[j] + delta * (shifts[2] - shifts[1])
}

# -----------------------------------------------------------------------------
# 1b. Offset in ONE window: shift the target series in time and maximise
#     correlation with the reference series on a fine common grid, so lags
#     finer than the native averaging interval are resolvable.
#     Returns offset (sec, to ADD to target time) and peak correlation.
# -----------------------------------------------------------------------------
estimate_offset_window <- function(ref_t, ref_v, tgt_t, tgt_v,
                                   max_lag_sec = 1800, fine_dt_sec = 30,
                                   min_points = 20, max_gap_sec = 3600) {
  if (length(ref_t) != length(ref_v) || length(tgt_t) != length(tgt_v))
    stop("Each signal must have one value per timestamp.")
  ref_t <- as.numeric(ref_t); tgt_t <- as.numeric(tgt_t)
  kr <- is.finite(ref_t); kt <- is.finite(tgt_t)
  ref_v <- ref_v[kr]; ref_t <- ref_t[kr]
  tgt_v <- tgt_v[kt]; tgt_t <- tgt_t[kt]
  if (sum(is.finite(ref_v)) < min_points || sum(is.finite(tgt_v)) < min_points)
    return(list(offset = NA_real_, corr = NA_real_, n = 0L))

  lo <- max(min(ref_t), min(tgt_t)) + max_lag_sec
  hi <- min(max(ref_t), max(tgt_t)) - max_lag_sec
  if (hi - lo < min_points * fine_dt_sec)
    return(list(offset = NA_real_, corr = NA_real_, n = 0L))

  tg    <- seq(lo, hi, by = fine_dt_sec)
  ref_g <- .gap_approx(ref_t, ref_v, tg, max_gap_sec)
  if (!is.finite(sd(ref_g, na.rm = TRUE)) || sd(ref_g, na.rm = TRUE) == 0)
    return(list(offset = NA_real_, corr = NA_real_, n = length(tg)))

  shifts <- seq(-max_lag_sec, max_lag_sec, by = fine_dt_sec)
  corr <- vapply(shifts, function(s) {
    tgt_g <- .gap_approx(tgt_t + s, tgt_v, tg, max_gap_sec)
    paired <- is.finite(ref_g) & is.finite(tgt_g)
    if (sum(paired) < min_points) return(NA_real_)
    suppressWarnings(cor(ref_g[paired], tgt_g[paired]))
  }, numeric(1))
  if (all(!is.finite(corr))) return(list(offset = NA_real_, corr = NA_real_, n = length(tg)))
  list(offset = .refine_peak(shifts, corr), corr = max(corr, na.rm = TRUE), n = length(tg))
}

# -----------------------------------------------------------------------------
# 1c. Slide the window across the record -> a series of anchor offsets.
#     Keep windows MULTI-DAY so the drift signal doesn't alias with the
#     reference variable's own diurnal cycle. ref_* = trusted time base;
#     tgt_* = stream being corrected.
# -----------------------------------------------------------------------------
estimate_offsets <- function(ref_time, ref_val, tgt_time, tgt_val,
                             window_days = 7, step_days = 2,
                             max_lag_sec = 1800, fine_dt_sec = 30,
                             min_corr = 0.7, max_gap_sec = 3600) {
  ref_time <- as.numeric(ref_time); tgt_time <- as.numeric(tgt_time)
  if (length(ref_time) != length(ref_val) || length(tgt_time) != length(tgt_val))
    stop("Each signal must have one value per timestamp.")
  if (any(!is.finite(c(window_days, step_days, max_lag_sec, fine_dt_sec))) ||
      window_days <= 0 || step_days <= 0 || max_lag_sec < 0 || fine_dt_sec <= 0)
    stop("Window, step and grid sizes must be positive; max_lag_sec must be nonnegative.")
  kr <- is.finite(ref_time); kt <- is.finite(tgt_time)
  ref_val <- ref_val[kr]; ref_time <- ref_time[kr]
  tgt_val <- tgt_val[kt]; tgt_time <- tgt_time[kt]
  if (!length(ref_time) || !length(tgt_time)) return(.empty_anchors())
  t0 <- max(min(ref_time), min(tgt_time)); t1 <- min(max(ref_time), max(tgt_time))
  win <- window_days * 86400; step <- step_days * 86400
  if (t1 - t0 < win) return(.empty_anchors())
  starts <- seq(t0, t1 - win, by = step)
  rows <- lapply(starts, function(ws) {
    we <- ws + win
    ri <- ref_time >= ws & ref_time < we
    ti <- tgt_time >= ws & tgt_time < we
    r  <- estimate_offset_window(ref_time[ri], ref_val[ri], tgt_time[ti], tgt_val[ti],
                                 max_lag_sec, fine_dt_sec, max_gap_sec = max_gap_sec)
    data.frame(t_center = ws + win / 2, offset = r$offset, corr = r$corr, n = r$n)
  })
  out <- do.call(rbind, rows)
  out <- out[is.finite(out$offset) & is.finite(out$corr) & out$corr >= min_corr, ]
  out$t_center <- as.POSIXct(out$t_center, origin = "1970-01-01", tz = "UTC")
  rownames(out) <- NULL
  out
}

# -----------------------------------------------------------------------------
# 1d. Solar-primary anchors with shared-signal gap-fill. Anchors the corrected
#     stream to ABSOLUTE solar time via its SW_in; in windows where the solar
#     anchor is weak/missing (night-heavy, overcast), it falls back to the
#     shared Ts-Ta signal. NOTE the shared fallback assumes the OTHER stream is
#     ~on true time; solar windows are preferred wherever available.
# -----------------------------------------------------------------------------
solar_blend_anchors <- function(tgt_time, tgt_sw, tgt_shared,
                                other_shared_time, other_shared_val,
                                lat, lon, ...) {
  rpot   <- potential_radiation(tgt_time, lat, lon)
  a_sol  <- estimate_offsets(tgt_time, rpot, tgt_time, tgt_sw, ...)
  a_shar <- estimate_offsets(other_shared_time, other_shared_val, tgt_time, tgt_shared, ...)
  a_sol$src <- rep("solar", nrow(a_sol)); a_shar$src <- rep("shared", nrow(a_shar))
  # prefer solar windows; add shared only where no solar anchor exists nearby
  keep_shar <- vapply(as.numeric(a_shar$t_center), function(tc)
    !length(a_sol$t_center) || min(abs(as.numeric(a_sol$t_center) - tc)) > 12 * 3600,
    logical(1))
  out <- rbind(a_sol, a_shar[keep_shar, , drop = FALSE])
  out[order(out$t_center), ]
}

# -----------------------------------------------------------------------------
# 2. Fit rate(T) = a*(T-T0)^2 + b to the anchor offsets by integrating over the
#    corrected logger's PANEL temperature (PTemp). a,b in ppm (1 ppm = 1e-6 s/s);
#    T0 degC; s0 = offset (sec) at record start. b absorbs baseline drift.
#    Anchors weighted by corr^2.
# -----------------------------------------------------------------------------
fit_temperature_drift <- function(anchors, met_time, met_ptemp,
                                   T0_start = 25, a_start = -0.03, b_start = 0,
                                   weight_by_corr = TRUE, max_ptemp_gap_sec = 7200) {
  mt <- as.numeric(met_time)
  if (length(mt) != length(met_ptemp) || length(mt) < 2L || any(!is.finite(mt)))
    stop("Panel temperature requires at least two finite timestamps and matching values.")
  ord <- order(mt); mt <- mt[ord]; Tp <- met_ptemp[ord]
  if (anyDuplicated(mt)) stop("Duplicate panel-temperature timestamps must be resolved.")
  valid <- is.finite(Tp)
  Tp <- .gap_approx(mt[valid], Tp[valid], mt, max_ptemp_gap_sec)
  if (any(!is.finite(Tp)) || any(diff(mt) > max_ptemp_gap_sec))
    stop("Panel-temperature coverage has an edge gap or exceeds max_ptemp_gap_sec; supply temperature data or split the record.")
  dt <- c(0, diff(mt))
  ac <- as.numeric(anchors$t_center)
  w <- if (weight_by_corr) anchors$corr^2 else rep(1, nrow(anchors))
  if (length(ac) < 5L || length(unique(ac)) < 5L ||
      any(!is.finite(ac)) || any(!is.finite(anchors$offset)) ||
      any(!is.finite(w) | w <= 0))
    stop("At least five distinct, finite anchors with positive weights are required.")
  if (any(ac < mt[1] | ac > tail(mt, 1)))
    stop("Anchor times must lie within the panel-temperature record.")
  # Test identifiability of the integrated quadratic before nonlinear fitting.
  basis <- cbind(1, sapply(list(Tp^2, Tp, rep(1, length(Tp))), function(z)
    approx(mt, cumsum(z * dt), xout = ac)$y))
  scales <- sqrt(colSums(basis^2))
  if (any(scales == 0) || qr(sweep(basis, 2, scales, "/"))$rank < 4L)
    stop("Temperature history cannot identify all four drift parameters.")
  cumint_fun <- function(a, b, T0) cumsum((a * (Tp - T0)^2 + b) * 1e-6 * dt)
  obj <- function(p) {
    ci   <- cumint_fun(p[1], p[2], p[3])
    pred <- p[4] + approx(mt, ci, xout = ac, rule = 2)$y
    sum(w * (anchors$offset - pred)^2)
  }
  init <- c(a_start, b_start, T0_start, median(anchors$offset, na.rm = TRUE))
  # The integrated polynomial supplies a stable starting point for BFGS.
  scaled_basis <- sweep(basis, 2, scales, "/")
  coef <- lm.wfit(scaled_basis, anchors$offset, w)$coefficients / scales
  aa <- coef[2] * 1e6; linear <- coef[3] * 1e6; constant <- coef[4] * 1e6
  if (is.finite(aa) && abs(aa) > 1e-8) {
    turn <- -linear / (2 * aa)
    init <- c(aa, constant - aa * turn^2, turn, coef[1])
  } else if (is.finite(linear) && abs(linear) < 1e-6) {
    init <- c(0, constant, T0_start, coef[1])
  }
  fit  <- optim(init, obj, method = "BFGS",
                control = list(parscale = c(0.01, 1, 5, 60), maxit = 500))
  if (fit$convergence != 0L || any(!is.finite(fit$par)) || !is.finite(fit$value))
    stop("Temperature drift optimization failed to converge; no correction was applied.")
  a <- fit$par[1]; b <- fit$par[2]; T0 <- fit$par[3]; s0 <- fit$par[4]
  pred <- s0 + approx(mt, cumint_fun(a, b, T0), xout = ac, rule = 2)$y
  list(par = c(a = a, b = b, T0 = T0, s0 = s0), convergence = fit$convergence,
       integration_time = mt, integration_offset = s0 + cumint_fun(a, b, T0),
       anchors = transform(anchors, predicted = pred, resid = anchors$offset - pred),
       rate_ppm = function(T) a * (T - T0)^2 + b)
}

# -----------------------------------------------------------------------------
# 3. Reconstruct continuous offset(t) from the stored fitted integral over the
#    full PTemp record (fills gaps in anchor coverage). Adds `offset_sec` and
#    the corrected timestamp `time_corr` to the corrected stream.
# -----------------------------------------------------------------------------
reconstruct_offset <- function(df, fit, time_col = "TIMESTAMP", ptemp_col = "PTemp") {
  mt <- as.numeric(df[[time_col]])
  if (is.null(fit$integration_time) || is.null(fit$integration_offset))
    stop("Fit lacks an integration curve; refit with fit_temperature_drift().")
  if (any(!is.finite(mt)) || any(mt < min(fit$integration_time) | mt > max(fit$integration_time)))
    stop("Reconstruction times must be finite and within the fitted temperature record.")
  # Reuse the exact fitted integral, including its origin and temperature gap policy.
  df$offset_sec <- approx(fit$integration_time, fit$integration_offset, xout = mt)$y
  df$time_corr  <- df[[time_col]] + df$offset_sec
  df
}

# -----------------------------------------------------------------------------
# 4. Reindex the drift-corrected stream onto the reference timestamps (linear
#    interpolation). Appends "<var>_corr" columns to the reference frame.
# -----------------------------------------------------------------------------
align_streams <- function(ref_df, corr_df, ref_time = "TIMESTAMP",
                          corr_time = "time_corr", vars, max_gap_sec = 3600) {
  ft <- as.numeric(ref_df[[ref_time]]); mt <- as.numeric(corr_df[[corr_time]])
  if (any(!vars %in% names(corr_df))) stop("Unknown output variable.")
  for (v in vars) ref_df[[paste0(v, "_corr")]] <-
    .gap_approx(mt, corr_df[[v]], ft, max_gap_sec)
  ref_df
}

# -----------------------------------------------------------------------------
# 5. Diagnostics (ggplot2). The residual-vs-PTemp panel is the key check the
#    temperature dependence was captured: it should be flat and structureless.
# -----------------------------------------------------------------------------
plot_drift_diagnostics <- function(fit, corr_df = NULL, ptemp_col = "PTemp", outdir = ".") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required for diagnostics.")
  library(ggplot2)
  a <- fit$anchors; pr <- fit$par
  src_aes <- if (!is.null(a$src)) aes(y = offset, size = corr, colour = src) else aes(y = offset, size = corr)
  p_offset <- ggplot(a, aes(t_center)) +
    geom_point(src_aes, alpha = 0.6) +
    geom_line(aes(y = predicted), colour = "#c53030", linewidth = 0.8) +
    labs(x = NULL, y = "Offset (s, added to corrected time)",
         title = "Reconstructed clock offset over time",
         subtitle = "points = anchors; line = temperature-driven fit") +
    theme_minimal(base_size = 12)

  Tseq <- seq(-10, 55, by = 0.5)
  p_rate <- ggplot(data.frame(T = Tseq, rate = fit$rate_ppm(Tseq)), aes(T, rate)) +
    geom_line(colour = "#805ad5", linewidth = 0.9) + geom_hline(yintercept = 0, linetype = 3) +
    labs(x = "Panel temperature (°C)", y = "Drift rate (ppm)",
         title = "Fitted temperature dependence of drift rate",
         subtitle = sprintf("a=%.4f ppm/°C²,  T0=%.1f °C,  b=%.3f ppm", pr["a"], pr["T0"], pr["b"])) +
    theme_minimal(base_size = 12)

  p_resid <- NULL
  if (!is.null(corr_df)) {
    Tp_at <- approx(as.numeric(corr_df$time_corr - corr_df$offset_sec), corr_df[[ptemp_col]],
                    xout = as.numeric(a$t_center), rule = 2)$y
    p_resid <- ggplot(data.frame(PTemp = Tp_at, resid = a$resid), aes(PTemp, resid)) +
      geom_hline(yintercept = 0, linetype = 3) + geom_point(alpha = 0.6, colour = "#2f855a") +
      geom_smooth(method = "loess", se = FALSE, colour = "#c53030", linewidth = 0.7) +
      labs(x = "Panel temperature (°C)", y = "Residual offset (s)",
           title = "Residual vs. panel temperature (the diagnostic that matters)",
           subtitle = "flat & structureless = temperature dependence captured") +
      theme_minimal(base_size = 12)
  }
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(outdir, "drift_offset_timeline.png"), p_offset, width = 9, height = 4, dpi = 150)
  ggsave(file.path(outdir, "drift_rate_vs_temp.png"),    p_rate,   width = 6, height = 4, dpi = 150)
  if (!is.null(p_resid))
    ggsave(file.path(outdir, "drift_residual_vs_temp.png"), p_resid, width = 6, height = 4, dpi = 150)
  invisible(list(offset = p_offset, rate = p_rate, resid = p_resid))
}

# -----------------------------------------------------------------------------
# 6. Column mapping. Point the pipeline at the columns of YOUR dataframe.
#    A mapping is a named list with keys: time, shared, sw, ptemp.
#      time   = timestamp column (POSIXct)
#      shared = shared signal for cross-correlation (sonic temp / air temp)
#      sw     = incoming shortwave, for solar anchoring   (NA if unused)
#      ptemp  = logger panel temperature; REQUIRED on the CORRECTED stream
#    Defaults follow AmeriFlux-ish names; override any key you need.
# -----------------------------------------------------------------------------
flux_columns <- function(time = "TIMESTAMP", shared = "Ts", sw = "SW_IN", ptemp = NA)
  list(time = time, shared = shared, sw = sw, ptemp = ptemp)
met_columns  <- function(time = "TIMESTAMP", shared = "Ta", sw = "SW_IN", ptemp = "PTemp")
  list(time = time, shared = shared, sw = sw, ptemp = ptemp)

# Auto-detect likely columns from a dataframe you provide (case-insensitive,
# common EC/AmeriFlux aliases). Returns a mapping list to review/edit, not to
# trust blindly -- always eyeball it against names(df).
suggest_columns <- function(df, panel = TRUE) {
  nm <- names(df); low <- tolower(nm)
  pick <- function(cands) { i <- which(low %in% tolower(cands))[1]; if (is.na(i)) NA else nm[i] }
  list(
    time   = pick(c("TIMESTAMP","TIMESTAMP_START","TIMESTAMP_END","datetime","date_time","time")),
    shared = pick(c("Ts","T_SONIC","TSONIC","sonic_temperature","Ta","TA","T_air","air_temperature","AirTC")),
    sw     = pick(c("SW_IN","Rg","SWIN","shortwave_in","incoming_shortwave","Rg_in","SW_IN_1_1_1")),
    ptemp  = if (panel) pick(c("PTemp","PTemp_C","panel_temp","T_panel","PTemp_Avg","LoggerT")) else NA)
}

# validate a mapping against the dataframe actually supplied
.resolve_cols <- function(df, cols, which, need_ptemp = FALSE) {
  req <- c("time", "shared")
  if (need_ptemp) req <- c(req, "ptemp")
  if (any(!req %in% names(cols))) stop("Missing required column mapping: ", which)
  for (k in names(cols)) {
    v <- cols[[k]]
    if (is.null(v) || (length(v) == 1 && is.na(v))) {
      if (k %in% req || (k == "ptemp" && need_ptemp))
        stop(sprintf("%s_cols$%s is required but not set.", which, k))
      next
    }
    if (!v %in% names(df))
      stop(sprintf("Column '%s' (mapped as %s_cols$%s) is not in the %s dataframe.\n  Available: %s",
                   v, which, k, which, paste(names(df), collapse = ", ")))
  }
  cols
}

# =============================================================================
# ORCHESTRATOR. Provide your dataframes and map their columns.
# -----------------------------------------------------------------------------
# reference:
#   "flux"   trust flux clock; correct the other stream to it (shared signal)
#   "met"    trust met clock;  correct the other stream to it (shared signal)
#   "solar"  pin the corrected stream to solar time via its SW_in vs potential
#   "shared" relative alignment between the two streams (no absolute pin)
#   "custom" YOU supply ref_time & ref_val (e.g. a GPS-synced signal)
#   "solar_blend" solar-primary with shared-signal gap-fill
# correct: which stream to drift-correct, "met" (default) or "flux".
# Timestamps assumed UTC. Period labels declared per stream (start/mid/end).
# Columns come from flux_cols / met_cols (see flux_columns()/met_columns()).
# =============================================================================
align_clock_drift <- function(
    flux, met,
    reference = "solar_blend", correct = "met",
    # --- map YOUR dataframe's columns here ---
    flux_cols = flux_columns(),
    met_cols  = met_columns(),
    # --- period-label convention per stream ---
    flux_label = "end", met_label = "start", interval_min = 30,
    # --- solar geometry (needed for solar / solar_blend) ---
    lat = NA, lon = NA,
    # --- custom reference (reference = "custom") ---
    ref_time = NULL, ref_val = NULL,
    # --- anchoring / fit knobs ---
    window_days = 7, step_days = 2, max_lag_sec = 1800, min_corr = 0.7,
    out_vars = NULL, diagnostics = TRUE, outdir = "drift_diag",
    max_gap_sec = 3600, max_ptemp_gap_sec = 7200) {

  reference <- match.arg(reference, c("flux", "met", "solar", "shared", "custom", "solar_blend"))
  correct <- match.arg(correct, c("met", "flux"))
  if (reference == correct) stop("`reference` and `correct` cannot be the same stream.")
  if (reference == "custom" && (is.null(ref_time) || is.null(ref_val)))
    stop("Custom anchoring requires ref_time and ref_val on midpoint timestamps.")

  # which stream is corrected? its mapping must include ptemp
  need_ptemp_flux <- correct == "flux"
  fc <- .resolve_cols(flux, flux_cols, "flux", need_ptemp = need_ptemp_flux)
  mc <- .resolve_cols(met,  met_cols,  "met",  need_ptemp = !need_ptemp_flux)

  # 0. normalize period labels to a common midpoint instant
  flux[[fc$time]] <- normalize_period_label(flux[[fc$time]], flux_label, interval_min, "mid")
  met[[mc$time]]  <- normalize_period_label(met[[mc$time]],  met_label,  interval_min, "mid")

  if (correct == "met") { tgt <- met;  oth <- flux; tc <- mc; oc <- fc }
  else                  { tgt <- flux; oth <- met;  tc <- fc; oc <- mc }

  if (reference %in% c("solar", "solar_blend") &&
      (is.null(tc$sw) || length(tc$sw) != 1L || is.na(tc$sw)))
    stop("Solar anchoring requires the corrected stream's sw column mapping.")

  # 1. build anchors
  anchors <- switch(reference,
    flux = , met = {
      if ((reference == "flux" && correct == "flux") ||
          (reference == "met"  && correct == "met"))
        stop("`reference` and `correct` cannot be the same stream.")
      estimate_offsets(oth[[oc$time]], oth[[oc$shared]],
                       tgt[[tc$time]], tgt[[tc$shared]],
                       window_days, step_days, max_lag_sec, min_corr = min_corr, max_gap_sec = max_gap_sec)
    },
    shared = estimate_offsets(oth[[oc$time]], oth[[oc$shared]],
                              tgt[[tc$time]], tgt[[tc$shared]],
                              window_days, step_days, max_lag_sec, min_corr = min_corr, max_gap_sec = max_gap_sec),
    solar  = estimate_offsets(tgt[[tc$time]], potential_radiation(tgt[[tc$time]], lat, lon),
                              tgt[[tc$time]], tgt[[tc$sw]],
                              window_days, step_days, max_lag_sec, min_corr = min_corr, max_gap_sec = max_gap_sec),
    custom = estimate_offsets(ref_time, ref_val, tgt[[tc$time]], tgt[[tc$shared]],
                              window_days, step_days, max_lag_sec, min_corr = min_corr, max_gap_sec = max_gap_sec),
    solar_blend = solar_blend_anchors(
                    tgt[[tc$time]], tgt[[tc$sw]], tgt[[tc$shared]],
                    oth[[oc$time]], oth[[oc$shared]], lat, lon,
                    window_days = window_days, step_days = step_days,
                    max_lag_sec = max_lag_sec, min_corr = min_corr, max_gap_sec = max_gap_sec),
    stop("unknown reference: ", reference))

  if (nrow(anchors) < 5) stop("Too few anchors (", nrow(anchors),
                              "). At least five are required; check record duration, window settings and signal quality.")

  # 2-3. fit temperature-driven drift on the CORRECTED logger's PTemp; reconstruct
  fit     <- fit_temperature_drift(anchors, tgt[[tc$time]], tgt[[tc$ptemp]],
                                     max_ptemp_gap_sec = max_ptemp_gap_sec)
  tgt_cor <- reconstruct_offset(tgt, fit, tc$time, tc$ptemp)

  # 4. reindex corrected stream onto the reference (the stream NOT corrected)
  if (is.null(out_vars))
    out_vars <- setdiff(names(tgt_cor)[vapply(tgt_cor, is.numeric, logical(1))], c(tc$time, "offset_sec", "time_corr"))
  merged <- align_streams(oth, tgt_cor, ref_time = oc$time, corr_time = "time_corr", vars = out_vars,
                          max_gap_sec = max_gap_sec)

  # 5. diagnostics + post-correction verification anchors
  plots <- if (diagnostics) plot_drift_diagnostics(fit, tgt_cor, tc$ptemp, outdir) else NULL
  solar_ref <- reference %in% c("solar", "solar_blend")
  check <- tryCatch(estimate_offsets(
      if (solar_ref) tgt_cor[["time_corr"]] else if (reference == "custom") ref_time else oth[[oc$time]],
      if (solar_ref) potential_radiation(tgt_cor[["time_corr"]], lat, lon)
        else if (reference == "custom") ref_val else oth[[oc$shared]],
      tgt_cor[["time_corr"]],
      if (solar_ref) tgt_cor[[tc$sw]] else tgt_cor[[tc$shared]],
      window_days, step_days, max_lag_sec, min_corr = min_corr, max_gap_sec = max_gap_sec),
    error = function(e) { warning("Post-correction verification failed: ", conditionMessage(e)); NULL })

  list(anchors = anchors, fit = fit, corrected = tgt_cor, merged = merged,
       residual_offsets = if (!is.null(check)) check$offset else NULL, plots = plots)
}

# =============================================================================
# EXAMPLE USAGE  (guarded: sourcing this file only DEFINES the functions).
# =============================================================================
if (FALSE) {

  flux <- read.csv("flux_30min.csv");  met <- read.csv("biomet_30min.csv")
  flux$TIMESTAMP <- ymd_hms(flux$TIMESTAMP, tz = "UTC")   # timestamps are UTC
  met$TIMESTAMP  <- ymd_hms(met$TIMESTAMP,  tz = "UTC")

  # Option 1: auto-detect column names from your dataframes, then review/edit
  suggest_columns(flux, panel = FALSE)   # -> list(time=, shared=, sw=, ptemp=)
  suggest_columns(met)

  res <- align_clock_drift(
    flux, met,
    reference = "solar_blend",           # "flux"/"met"/"solar"/"shared"/"custom"
    correct   = "met",                   # drift-correct the biomet stream
    # Option 2: map columns explicitly to whatever YOUR dataframe uses
    flux_cols = flux_columns(time = "TIMESTAMP", shared = "Ts", sw = "SW_IN"),
    met_cols  = met_columns (time = "TIMESTAMP", shared = "Ta", sw = "SW_IN", ptemp = "PTemp"),
    flux_label = "end", met_label = "start", interval_min = 30,
    lat = 25.4, lon = -80.6,             # site coords for solar anchoring
    window_days = 7, step_days = 2, min_corr = 0.7)

  print(res$fit$par); print(res$fit$convergence)
  summary(res$residual_offsets)          # post-correction residuals (sec) ~ 0 if good
  head(res$merged)                       # reference stream + "<var>_corr" columns

  # reference = "custom": supply your own trusted signal
  # res <- align_clock_drift(flux, met, reference = "custom",
  #                          ref_time = gps$TIMESTAMP, ref_val = gps$signal,
  #                          met_cols = met_columns(ptemp = "PTemp"))
}
