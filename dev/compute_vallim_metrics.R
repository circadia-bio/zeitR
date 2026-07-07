# dev/compute_vallim_metrics.R
# R port of compute_sleep_metrics() and compute_cpd_metrics()
# from pipeline_functions_fix27.py
# Run: source("dev/compute_vallim_metrics.R"); print_vallim_metrics(result)

# ── Brazilian national holidays (from pipeline_functions_fix27.py) ────────────
.FERIADOS_BR <- as.Date(c(
  "2018-09-07", "2018-10-12", "2018-11-02", "2018-11-15",
  "2018-11-20", "2018-12-25",
  "2019-01-01", "2019-01-25", "2019-03-04", "2019-03-05",
  "2019-04-19", "2019-04-21", "2019-05-01"
))

.is_free_day <- function(date) {
  d <- as.Date(date)
  weekdays(d, abbreviate = FALSE) %in% c("Saturday", "Sunday") | d %in% .FERIADOS_BR
}

.is_free_day_eve <- function(bts, tz = "America/Sao_Paulo") {
  .is_free_day(as.Date(as.POSIXct(bts, tz = tz)) + 1L)
}

# Plain decimal hours 0-24 (mirrors _series_to_h)
.to_h_plain <- function(x, tz = "America/Sao_Paulo") {
  x <- as.POSIXct(x, tz = tz)
  as.numeric(format(x, "%H", tz = tz)) +
  as.numeric(format(x, "%M", tz = tz)) / 60 +
  as.numeric(format(x, "%S", tz = tz)) / 3600
}

# Bed-time hours: if h < 12 add 24 (mirrors compute_sleep_metrics _to_h)
.to_h_bed <- function(x, tz = "America/Sao_Paulo") {
  h <- .to_h_plain(x, tz)
  ifelse(h >= 12, h, h + 24)
}

# Circular mean: values >= 12 -> subtract 24 (mirrors _mean_circular_h)
.mean_circ_h <- function(h) {
  s <- ifelse(h >= 12, h - 24, h)
  mean(s, na.rm = TRUE) %% 24
}

# Mid-sleep formula (mirrors _midsleep)
.midsleep <- function(onset_h, offset_h) {
  offset_adj <- ifelse(offset_h < onset_h, offset_h + 24, offset_h)
  so_adj     <- ifelse(onset_h >= 12, onset_h - 24, onset_h)
  (so_adj + (offset_adj - onset_h) / 2) %% 24
}

# HH:MM formatter (mirrors _hms)
.hms <- function(h) {
  h  <- h %% 24
  hh <- as.integer(h)
  mm <- as.integer(round((h %% 1) * 60))
  if (mm == 60L) { hh <- hh + 1L; mm <- 0L }
  sprintf("%02d:%02d", hh %% 24L, mm)
}

# ── compute_sleep_metrics (mirrors Cell 9) ────────────────────────────────────
compute_sleep_metrics_r <- function(nights,
                                     min_tib_h = 5.0,
                                     tz = "America/Sao_Paulo") {
  nd <- nights[!nights$is_nap & nights$tbt / 60 >= min_tib_h, ]
  if (nrow(nd) == 0) return(list())

  nd$is_free_day  <- .is_free_day(as.Date(nd$get_up_time, tz = tz))
  nd$bed_time_h   <- .to_h_bed(nd$bed_time, tz)
  nd$getup_time_h <- .to_h_plain(nd$get_up_time, tz)
  nd$onset_time_h <- nd$bed_time_h + nd$sol / 60
  nd$mid_sleep_h  <- (nd$bed_time_h + nd$tbt / 60 / 2) %% 24

  .metrics <- function(sub, dtype) {
    pf <- function(v) paste0(v, "_", dtype)
    if (nrow(sub) == 0) {
      r <- setNames(rep(NA_real_, 8),
                    pf(c("sol","tst","tbt","waso","eff","bed_time","getup_time","mid_sleep")))
      return(c(setNames(0L, pf("n")), r))
    }
    setNames(
      c(nrow(sub), mean(sub$sol), mean(sub$tst), mean(sub$tbt),
        mean(sub$waso), mean(sub$eff), mean(sub$bed_time_h),
        mean(sub$getup_time_h), mean(sub$mid_sleep_h)),
      c(pf("n"), pf(c("sol","tst","tbt","waso","eff","bed_time","getup_time","mid_sleep")))
    )
  }

  as.list(c(
    .metrics(nd,                       "overall"),
    .metrics(nd[!nd$is_free_day, ],    "weekday"),
    .metrics(nd[ nd$is_free_day, ],    "weekend")
  ))
}

# ── compute_cpd_metrics (mirrors compute_cpd_metrics + nights_to_cpd_df) ─────
compute_cpd_r <- function(nights,
                           min_tib_h     = 3.0,
                           min_tib_eve_h = 3.0,
                           tz = "America/Sao_Paulo") {
  nd <- nights[!nights$is_nap & nights$tbt / 60 >= min_tib_h, ]
  if (nrow(nd) == 0) stop("No nocturnal episodes after TIB filter.")

  onset_h  <- .to_h_plain(nd$bed_time,    tz) + nd$sol / 60
  offset_h <- .to_h_plain(nd$get_up_time, tz) - nd$soi / 60

  df <- data.frame(
    sleep_onset     = onset_h,
    sleep_offset    = offset_h,
    time_in_bed     = nd$tbt / 60,
    mid_sleep       = .midsleep(onset_h, offset_h),
    is_free_day     = .is_free_day(as.Date(nd$get_up_time, tz = tz)),
    is_free_day_eve = .is_free_day_eve(nd$bed_time, tz)
  )
  df$is_free_day_eve[df$time_in_bed < min_tib_eve_h] <- FALSE

  we  <- df$is_free_day
  wd  <- !df$is_free_day
  eve <- df$is_free_day_eve

  if (!any(we))  stop("No free days found.")
  if (!any(wd))  stop("No work days found.")
  if (!any(eve)) stop("No free-day-eve nights found.")

  msf_h <- mean(df$mid_sleep[we])
  msw_h <- mean(df$mid_sleep[wd])

  sd_f    <- mean((df$sleep_offset[we] - df$sleep_onset[we]) %% 24)
  sd_w    <- mean((df$sleep_offset[wd] - df$sleep_onset[wd]) %% 24)
  sd_week <- (5 * sd_w + 2 * sd_f) / 7
  so_f    <- .mean_circ_h(df$sleep_onset[eve])

  msfsc_h <- if (sd_f <= sd_w) (so_f + sd_f / 2) %% 24 else (so_f + sd_week / 2) %% 24
  msfsc_s <- msfsc_h * 3600

  ms_s <- df$mid_sleep * 3600
  x_i  <- msfsc_s - ms_s
  half  <- 12 * 3600; full <- 24 * 3600
  x_i  <- ifelse(x_i < -half, x_i + full, ifelse(x_i > half, x_i - full, x_i))
  y_i  <- c(0, -diff(ms_s))
  cpd_s <- mean(sqrt(x_i^2 + y_i^2))

  sjl_h  <- abs(msf_h - msw_h)
  sjla_h <- msf_h - msw_h

  list(
    n_nights_cpd = nrow(df), n_free_days = sum(we), n_workdays = sum(wd),
    msw_h = msw_h,   msw_hms   = .hms(msw_h),
    msf_h = msf_h,   msf_hms   = .hms(msf_h),
    msfsc_h = msfsc_h, msfsc_hms = .hms(msfsc_h),
    sjl_h = sjl_h,   sjl_min   = sjl_h * 60,
    sjla_h = sjla_h, sjla_min  = sjla_h * 60,
    cpd_s = cpd_s,   cpd_min   = cpd_s / 60, cpd_h = cpd_s / 3600
  )
}

# ── Pretty-print both tables ──────────────────────────────────────────────────
print_vallim_metrics <- function(result, tz = "America/Sao_Paulo") {
  sm  <- compute_sleep_metrics_r(result$nights, tz = tz)
  cpd <- tryCatch(compute_cpd_r(result$nights, tz = tz), error = function(e) {
    message("[CPD] ", conditionMessage(e)); NULL
  })

  cat("=== Sleep metrics (overall / weekday / weekend) ===\n")
  sm_df <- data.frame(
    metric = names(sm),
    value  = round(unlist(sm), 4),
    row.names = NULL
  )
  print(sm_df, row.names = FALSE)

  if (!is.null(cpd)) {
    cat("\n=== CPD / SJL metrics ===\n")
    nums <- cpd[sapply(cpd, is.numeric)]
    strs <- cpd[sapply(cpd, is.character)]
    cat(sprintf("  n_nights_cpd = %d   n_free_days = %d   n_workdays = %d\n",
                cpd$n_nights_cpd, cpd$n_free_days, cpd$n_workdays))
    cat(sprintf("  MSW  = %s  (%.4f h)\n", cpd$msw_hms,   cpd$msw_h))
    cat(sprintf("  MSF  = %s  (%.4f h)\n", cpd$msf_hms,   cpd$msf_h))
    cat(sprintf("  MSFsc= %s  (%.4f h)\n", cpd$msfsc_hms, cpd$msfsc_h))
    cat(sprintf("  SJL  = %.4f h  (%.2f min)\n", cpd$sjl_h,  cpd$sjl_min))
    cat(sprintf("  SJLa = %.4f h  (%.2f min)\n", cpd$sjla_h, cpd$sjla_min))
    cat(sprintf("  CPD  = %.2f s  (%.4f min)\n", cpd$cpd_s,  cpd$cpd_min))
  }
}
