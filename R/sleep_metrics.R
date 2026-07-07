# ── Sleep summary and CPD/SJL metrics ─────────────────────────────────────────
# R port of compute_sleep_metrics() and compute_cpd_metrics() from
# pipeline_functions_fix27.py (Julia Ribeiro da Silva Vallim, 2024).
# Column names and derived metrics mirror the Python output format exactly.

# ── Internal helpers ──────────────────────────────────────────────────────────

# Mirrors _series_to_h: plain decimal hours 0-24
.to_h_plain <- function(x, tz) {
  x <- as.POSIXct(x, tz = tz)
  as.numeric(format(x, "%H", tz = tz)) +
    as.numeric(format(x, "%M", tz = tz)) / 60 +
    as.numeric(format(x, "%S", tz = tz)) / 3600
}

# Mirrors _mean_circular_h: values >= 12 -> subtract 24 before averaging
.mean_circ_h <- function(h) {
  s <- ifelse(h >= 12, h - 24, h)
  mean(s, na.rm = TRUE) %% 24
}

# Circular SD in hours (analogous to circ_sd_h but without the 2*pi scaling)
# Kept for reference; sd() used directly in compute_sleep_metrics.
# .sd_circ_h <- function(h) { s <- ifelse(h >= 12, h - 24, h); sd(s, na.rm = TRUE) }

# Mirrors _midsleep
.midsleep <- function(onset_h, offset_h) {
  offset_adj <- ifelse(offset_h < onset_h, offset_h + 24, offset_h)
  so_adj     <- ifelse(onset_h >= 12, onset_h - 24, onset_h)
  (so_adj + (offset_adj - onset_h) / 2) %% 24
}

# Mirrors _hms
.hms_fmt <- function(h) {
  h  <- h %% 24
  hh <- as.integer(h)
  mm <- as.integer(round((h %% 1) * 60))
  if (mm == 60L) { hh <- hh + 1L; mm <- 0L }
  sprintf("%02d:%02d", hh %% 24L, mm)
}

# is_free_day applied to get-up date
.is_free_day <- function(dates, holidays) {
  d <- as.Date(dates)
  weekdays(d, abbreviate = FALSE) %in% c("Saturday", "Sunday") |
    (!is.null(holidays) & d %in% as.Date(holidays))
}

# is_free_day_eve: the bed-time date + 1 day is a free day
.is_free_day_eve <- function(bts, tz, holidays) {
  d_next <- as.Date(as.POSIXct(bts, tz = tz), tz = tz) + 1L
  .is_free_day(d_next, holidays)
}


# ── compute_sleep_metrics ─────────────────────────────────────────────────────

#' Compute sleep metrics split by day type
#'
#' Calculates mean sleep metrics separately for all nights, weekday nights,
#' and weekend/holiday nights. Output column names and derived metrics mirror
#' `compute_sleep_metrics()` from Julia Ribeiro da Silva Vallim's
#' `pipeline_functions_fix27.py`.
#'
#' @details
#' A night is assigned to the **weekend** group when the get-up date falls on a
#' Saturday, Sunday, or any date supplied in `holidays`. Weekday nights are all
#' remaining nights.
#'
#' `sleep_onset_h` uses a circular mean (values >= 12 h are shifted by −24 h
#' before averaging, then wrapped to [0, 24)). `sleep_offset_h` uses a plain
#' arithmetic mean.
#'
#' `fps_h` (free period sleep) equals `fpr_tib_h − (latencia_min + inertia_min)
#' / 60` — TBT net of sleep onset latency and sleep inertia.
#'
#' `dp_midsleep_min` and `dp_tst_min` are standard deviations of per-night
#' mid-sleep (minutes) and TST (minutes), respectively.
#'
#' @param nights A `tibble` of nightly sleep statistics as returned by
#'   [run_pipeline_native()] or [run_pipeline()]. Must contain at minimum the
#'   columns `is_nap`, `bed_time`, `get_up_time`, `tbt`, `tst`, `sol`, `soi`,
#'   `waso`, `eff`.
#' @param min_tib_h `numeric(1)`. Minimum total in-bed time (hours) for a night
#'   to be included. Default is `5.0` (matching the Python reference).
#' @param tz `character(1)`. Time zone for extracting clock hours from
#'   timestamps. Default is `"UTC"`.
#' @param holidays A `Date` vector of public holidays to treat as free days in
#'   addition to Saturdays and Sundays. Default is `NULL` (weekends only).
#'
#' @return A named list with metrics for three groups (`overall`, `wd` =
#'   weekday, `fd` = free day):
#'   \describe{
#'     \item{`n_overall`, `n_wd`, `n_fd`}{Night counts.}
#'     \item{`sleep_onset_h`, `sleep_offset_h`}{Circular mean onset and
#'       arithmetic mean offset in decimal hours.}
#'     \item{`fpr_tib_h`}{Mean TBT in hours.}
#'     \item{`fps_h`}{Mean free period sleep (TBT − SOL − SOI) in hours.}
#'     \item{`tst_h`}{Mean TST in hours.}
#'     \item{`latencia_min`, `inertia_min`}{Mean SOL and SOI in minutes.}
#'     \item{`waso_min`}{Mean WASO in minutes.}
#'     \item{`sleep_eff_pct`}{Mean sleep efficiency in percent (0–100).}
#'     \item{`tst_24h_h`}{Same as `tst_h` (24-h TST for main sleep only).}
#'     \item{`dp_midsleep_min`, `dp_tst_min`}{SD of mid-sleep and TST in
#'       minutes.}
#'   }
#'   Workday and free-day metrics carry the suffix `_wd` and `_fd`.
#'
#' @seealso [compute_cpd_metrics()], [run_pipeline_native()]
#'
#' @importFrom stats sd setNames
#' @export
#'
#' @examples
#' \dontrun{
#' result <- run_pipeline_native("recordings/P001.txt",
#'                               tz = "America/Sao_Paulo")
#'
#' sm <- compute_sleep_metrics(result$nights, tz = "America/Sao_Paulo")
#' sm$tst_h            # mean TST in hours
#' sm$sleep_onset_h    # mean sleep onset (circular, decimal hours)
#' sm$dp_midsleep_min  # within-person SD of mid-sleep in minutes
#' }
compute_sleep_metrics <- function(nights,
                                   min_tib_h = 5.0,
                                   tz        = "UTC",
                                   holidays  = NULL) {
  nd <- nights[!nights$is_nap & nights$tbt / 60 >= min_tib_h, ]
  if (nrow(nd) == 0L) {
    zeitr_warn("No nights remain after TIB filter (min_tib_h = {min_tib_h} h).")
    return(list())
  }

  nd$is_free_day <- .is_free_day(as.Date(nd$get_up_time, tz = tz), holidays)
  onset_h  <- .to_h_plain(nd$bed_time,    tz) + nd$sol / 60
  offset_h <- .to_h_plain(nd$get_up_time, tz) - nd$soi / 60
  mid_h    <- .midsleep(onset_h, offset_h)

.group_metrics <- function(sub, onset_sub, offset_sub, mid_sub, sfx) {
    pf <- function(nm) if (sfx == "") nm else paste0(nm, "_", sfx)
    if (nrow(sub) == 0L) {
      nms <- pf(c("sleep_onset_h", "sleep_offset_h", "fpr_tib_h", "fps_h",
                  "tst_h", "latencia_min", "inertia_min", "waso_min",
                  "sleep_eff_pct", "tst_24h_h", "dp_midsleep_min", "dp_tst_min"))
      return(c(setNames(0L, paste0("n", if (sfx == "") "_overall" else paste0("_", sfx))),
               setNames(rep(NA_real_, length(nms)), nms)))
    }
    tbt_h       <- mean(sub$tbt) / 60
    lat_min     <- mean(sub$sol)
    inertia_min <- mean(sub$soi)

    vals <- c(
      nrow(sub),
      .mean_circ_h(onset_sub),         # sleep_onset_h
      mean(offset_sub, na.rm = TRUE),  # sleep_offset_h
      tbt_h,                           # fpr_tib_h
      tbt_h - (lat_min + inertia_min) / 60,  # fps_h
      mean(sub$tst) / 60,              # tst_h
      lat_min,                         # latencia_min
      inertia_min,                     # inertia_min
      mean(sub$waso),                  # waso_min
      mean(sub$eff) * 100,             # sleep_eff_pct
      mean(sub$tst) / 60,              # tst_24h_h (same as tst_h for main-only)
      stats::sd(mid_sub, na.rm = TRUE) * 60,        # dp_midsleep_min
      stats::sd(sub$tst, na.rm = TRUE)        # dp_tst_min
    )
    n_sfx <- if (sfx == "") "n_overall" else paste0("n_", sfx)
    stats::setNames(vals, c(n_sfx, pf(c(
      "sleep_onset_h", "sleep_offset_h", "fpr_tib_h", "fps_h",
      "tst_h", "latencia_min", "inertia_min", "waso_min",
      "sleep_eff_pct", "tst_24h_h", "dp_midsleep_min", "dp_tst_min"
    ))))
  }

  wd <- nd$is_free_day
  as.list(c(
    .group_metrics(nd,          onset_h,       offset_h,       mid_h,       ""),
    .group_metrics(nd[!wd, ],   onset_h[!wd],  offset_h[!wd],  mid_h[!wd],  "wd"),
    .group_metrics(nd[ wd, ],   onset_h[ wd],  offset_h[ wd],  mid_h[ wd],  "fd")
  ))
}


# ── compute_cpd_metrics ───────────────────────────────────────────────────────

#' Compute CPD, MSF, MSW, MSFsc, SJL, and SJLa
#'
#' Calculates chronobiological phenotyping metrics from classified nightly sleep
#' data. Mirrors `compute_cpd_metrics()` and `nights_to_cpd_df()` from Julia
#' Ribeiro da Silva Vallim's `pipeline_functions_fix27.py`.
#'
#' @details
#' Mid-sleep is computed per night as:
#' \deqn{\text{MS} = \left(\text{SO} + \frac{\text{offset} - \text{onset}}{2}\right) \bmod 24}
#' where onset = bts + SOL and offset = gts − SOI (both in decimal hours).
#'
#' **MSW** and **MSF** use plain arithmetic means of per-night mid-sleep values.
#' **MSFsc** adjusts MSF by the free-day-eve sleep onset and the weighted
#' weekly mean sleep duration when free-day duration exceeds weekday duration.
#' **CPD** is the RMS distance of each night's mid-sleep from MSFsc in the
#' time × sequence plane.
#'
#' @param nights A `tibble` of nightly sleep statistics as returned by
#'   [run_pipeline_native()] or [run_pipeline()].
#' @param min_tib_h `numeric(1)`. Minimum TBT (hours) for inclusion. Default
#'   is `3.0`.
#' @param min_tib_eve_h `numeric(1)`. Minimum TBT (hours) for a night to
#'   qualify as a free-day-eve night. Default is `3.0`.
#' @param tz `character(1)`. Time zone for extracting clock hours. Default is
#'   `"UTC"`.
#' @param holidays A `Date` vector of public holidays to treat as free days.
#'   Default is `NULL` (weekends only).
#'
#' @return A named list with `n_nights_cpd`, `n_free_days`, `n_workdays`,
#'   `msw_h`, `msw_hms`, `msf_h`, `msf_hms`, `msfsc_h`, `msfsc_hms`,
#'   `sjl_h`, `sjl_min`, `sjla_h`, `sjla_min`, `cpd_s`, `cpd_min`, `cpd_h`.
#'
#' @seealso [compute_sleep_metrics()], [run_pipeline_native()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- run_pipeline_native("recordings/P001.txt",
#'                               tz = "America/Sao_Paulo")
#'
#' cpd <- compute_cpd_metrics(result$nights, tz = "America/Sao_Paulo")
#' cpd$sjl_min    # social jet lag in minutes
#' cpd$msf_hms    # mid-sleep on free days as HH:MM
#' cpd$cpd_min    # CPD in minutes
#' }
compute_cpd_metrics <- function(nights,
                                 min_tib_h     = 3.0,
                                 min_tib_eve_h = 3.0,
                                 tz            = "UTC",
                                 holidays      = NULL) {
  nd <- nights[!nights$is_nap & nights$tbt / 60 >= min_tib_h, ]
  if (nrow(nd) == 0L)
    zeitr_abort("No nocturnal episodes remain after TIB filter (min_tib_h = {min_tib_h} h).")

  onset_h  <- .to_h_plain(nd$bed_time,    tz) + nd$sol / 60
  offset_h <- .to_h_plain(nd$get_up_time, tz) - nd$soi / 60

  df <- data.frame(
    sleep_onset     = onset_h,
    sleep_offset    = offset_h,
    time_in_bed     = nd$tbt / 60,
    mid_sleep       = .midsleep(onset_h, offset_h),
    is_free_day     = .is_free_day(as.Date(nd$get_up_time, tz = tz), holidays),
    is_free_day_eve = .is_free_day_eve(nd$bed_time, tz, holidays)
  )
  df$is_free_day_eve[df$time_in_bed < min_tib_eve_h] <- FALSE

  we  <- df$is_free_day
  wd  <- !df$is_free_day
  eve <- df$is_free_day_eve

  if (!any(we))  zeitr_abort("No free days found. Supply public holidays via {.arg holidays}.")
  if (!any(wd))  zeitr_abort("No workdays found.")
  if (!any(eve)) zeitr_abort("No free-day-eve nights found.")

  msf_h <- mean(df$mid_sleep[we])
  msw_h <- mean(df$mid_sleep[wd])

  sd_f    <- mean((df$sleep_offset[we] - df$sleep_onset[we]) %% 24)
  sd_w    <- mean((df$sleep_offset[wd] - df$sleep_onset[wd]) %% 24)
  sd_week <- (5 * sd_w + 2 * sd_f) / 7
  so_f    <- .mean_circ_h(df$sleep_onset[eve])

  msfsc_h <- if (sd_f <= sd_w) (so_f + sd_f / 2) %% 24 else (so_f + sd_week / 2) %% 24
  msfsc_s <- msfsc_h * 3600

  ms_s  <- df$mid_sleep * 3600
  x_i   <- msfsc_s - ms_s
  half  <- 12 * 3600; full <- 24 * 3600
  x_i   <- ifelse(x_i < -half, x_i + full, ifelse(x_i > half, x_i - full, x_i))
  y_i   <- c(0, -diff(ms_s))
  cpd_s <- mean(sqrt(x_i^2 + y_i^2))

  sjl_h  <- abs(msf_h - msw_h)
  sjla_h <- msf_h - msw_h

  list(
    n_nights_cpd = nrow(df),
    n_free_days  = sum(we),
    n_workdays   = sum(wd),
    msw_h   = msw_h,   msw_hms   = .hms_fmt(msw_h),
    msf_h   = msf_h,   msf_hms   = .hms_fmt(msf_h),
    msfsc_h = msfsc_h, msfsc_hms = .hms_fmt(msfsc_h),
    sjl_h   = sjl_h,   sjl_min   = sjl_h  * 60,
    sjla_h  = sjla_h,  sjla_min  = sjla_h * 60,
    cpd_s   = cpd_s,   cpd_min   = cpd_s  / 60, cpd_h = cpd_s / 3600
  )
}
