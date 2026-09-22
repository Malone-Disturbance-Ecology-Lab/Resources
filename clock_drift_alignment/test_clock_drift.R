args <- commandArgs(trailingOnly = TRUE)
source(if (length(args)) args[1] else 'align_clock_drift_temperature.R')
expect_error <- function(expr, pattern) {
  msg <- tryCatch({ force(expr); NA_character_ }, error = conditionMessage)
  stopifnot(!is.na(msg), grepl(pattern, msg))
}
t <- as.POSIXct('2025-01-01', tz='UTC') + (0:1919)*1800
k <- seq_along(t)
v <- sin(k*2*pi/48) + .2*sin(k*2*pi/17)
p <- 22 + 8*sin(k/310) + 3*cos(k/71)
# Empty solar/shared anchor sets and short records.
z <- solar_blend_anchors(t, rep(NA_real_,length(t)), v, t,v,25,-80)
stopifnot(nrow(z)>0, all(z$src=='shared'))
z <- solar_blend_anchors(t, potential_radiation(t,25,-80), rep(NA_real_,length(t)), t,v,25,-80)
stopifnot(nrow(z)>0, all(z$src=='solar'))
z <- solar_blend_anchors(t, rep(NA_real_,length(t)), rep(NA_real_,length(t)), t,v,25,-80)
stopifnot(nrow(z)==0, nrow(estimate_offsets(t[1:48],v[1:48],t[1:48],v[1:48]))==0)
# Explicit NAs and absent rows are not silently filled.
x <- align_streams(data.frame(TIMESTAMP=0:4),data.frame(time_corr=0:4,x=c(0,NA,NA,NA,4)),vars='x')
stopifnot(all(is.na(x$x_corr[2:4])))
stopifnot(is.na(.gap_approx(c(0,10000),c(0,1),5000)))
stopifnot(identical(.gap_approx(0:2,c(1,NA,3),0:2),c(1,NA,3)))
# Known drift with one internal temperature gap; fitted and applied curves agree.
p[100] <- NA_real_
filled <- approx(as.numeric(t)[is.finite(p)],p[is.finite(p)],xout=as.numeric(t))$y
truth <- 120+cumsum((.02*(filled-25)^2+3)*1e-6*c(0,diff(as.numeric(t))))
ai <- seq(100,length(t)-100,by=100)
a <- data.frame(t_center=t[ai],offset=truth[ai],corr=1)
f <- fit_temperature_drift(a,t,p)
r <- reconstruct_offset(data.frame(TIMESTAMP=t,PTemp=p),f)
stopifnot(max(abs(f$anchors$predicted-r$offset_sec[ai]))<1e-8,
          max(abs(f$anchors$resid))<.1)
expect_error(fit_temperature_drift(a[1:3,],t,p),'five')
expect_error(fit_temperature_drift(a,t,rep(25,length(t))),'identify')
pbad <- p; pbad[1] <- NA
expect_error(fit_temperature_drift(a,t,pbad),'coverage')
pbad <- p; pbad[100:110] <- NA
expect_error(fit_temperature_drift(a,t,pbad),'coverage')
# Custom reference must drive verification, even when the other logger is shifted.
flux <- data.frame(ft=t+600, Ts=v)
met <- data.frame(mt=t+300, Ta=v, PTemp=filled, note='metadata')
fc <- flux_columns(time='ft', sw=NA)
mc <- met_columns(time='mt', sw=NA)
expect_error(align_clock_drift(flux,met,reference='met',flux_cols=fc,met_cols=mc),'same stream')
z <- align_clock_drift(flux,met,reference='custom',ref_time=t,ref_val=v,
  flux_cols=fc,met_cols=mc,flux_label='mid',met_label='mid',diagnostics=FALSE)
stopifnot(length(z$residual_offsets)>0, max(abs(z$residual_offsets))<2,
          abs(median(z$corrected$offset_sec)+300)<2, !'note_corr' %in% names(z$merged))
# Both relative correction directions with mapped columns.
z <- align_clock_drift(flux,met,reference='flux',flux_cols=fc,met_cols=mc,
  flux_label='mid',met_label='mid',diagnostics=FALSE)
stopifnot(max(abs(z$residual_offsets))<2)
flux$panel <- filled
z <- align_clock_drift(flux,met,reference='met',correct='flux',
  flux_cols=flux_columns(time='ft',sw=NA,ptemp='panel'),met_cols=mc,
  flux_label='mid',met_label='mid',diagnostics=FALSE)
stopifnot(max(abs(z$residual_offsets))<2)
# Failed optimization is rejected before reconstruction.
original_optim <- optim
optim <- function(...) list(convergence=1L,par=c(0,0,25,0),value=0)
expect_error(fit_temperature_drift(a,t,p),'converge')
optim <- original_optim
cat('All clock-drift regression checks passed.\n')
