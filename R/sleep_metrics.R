# ── Sleep summary and CPD/SJL metrics ─────────────────────────────────────────
# R port of compute_sleep_metrics() and compute_cpd_metrics() from
# pipeline_functions_fix27.py (Julia Ribeiro da Silva Vallim, 2024).
# Both functions mirror their Python counterparts exactly, including:
#   - Bed-time convention (h < 12 → h + 24, keeping overnight episodes continuous)
#   - Mid-sleep formula (_midsleep)
#   - Circular mean for sleep onset (_mean_circular_h)
#   - is_free_day applied to the get-up date (gts)
#   - is_free_day_eve applied to the bed-time date + 1 day (bts)

# ── Internal helpers ──────────────────────────────────────────────────────────

# Mirrors _series_to_h: plain decimal hours 0-24
.to_h_plain <- function(x, tz) {
  x <- as.POSIXct(x, tz = tz)
  as.numeric(format(x, "%H", tz = tz)) +
  as.numeric(format(x, "%M", tz = tz)) / 60 +
  as.numeric(format(x, "%S", tz = tz)) / 3600
}

# Mirrors compute_sleep_metrics _to_h: h < 12 -> h + 24
.to_h_bed <- function(x, tz) {
  h <- .to_h_plain(x, tz)
  ifelse(h >= 12, h, h + 24)
}

# Mirrors _mean_circular_h: values >= 12 -> subtract 24 before averaging
.mean_circ_h <- function(h) {
  s <- ifelse(h >= 12, h - 24, h)
  mean(s, na.rm = TRUE) %% 24
}

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
#' Calculates mean sleep metrics (SOL, TST, TBT, WASO, efficiency, sleep onset,
#' get-up time, and mid-sleep) separately for all nights, weekday nights, and
#' weekend/holiday nights. Mirrors `compute_sleep_metrics()` from Julia Ribeiro
#' da Silva Vallim's `pipeline_functions_fix27.py`.
#'
#' @details
#' Bed-time is expressed using an overnight convention — hours before noon are
#' shifted by +24 h so that, for example, 01:30 becomes 25.5 h. This keeps
#' consecutive overnight sleep times on a continuous scale for arithmetic
#' averaging. Get-up time is expressed in plain decimal hours (0–24).
#'
#' A night is assigned to the **weekend** group when the get-up date falls on a
#' Saturday, Sunday, or any date supplied in `holidays`. Weekday nights are all
#' remaining nights.
#'
#' @param nights A `tibble` of nightly sleep statistics as returned by
#'   [run_pipeline_native()] or [run_pipeline()]. Must contain at minimum the
#'   columns `is_nap`, `bed_time`, `get_up_time`, `tbt`, `tst`, `sol`, `waso`,
#'   `eff`.
#' @param min_tib_h `numeric(1)`. Minimum total in-bed time (hours) for a night
#'   to be included. Default is `5.0` (matching the Python reference).
#' @param tz `character(1)`. Time zone for extracting clock hours from
#'   timestamps. Default is `"UTC"`.
#' @param holidays A `Date` vector of public holidays to treat as free days in
#'   addition to Saturdays and Sundays. Default is `NULL` (weekends only).
#'
#' @return A named list with 3 × 9 elements:
#'   \describe{
#'     \item{`n_overall`, `n_weekday`, `n_weekend`}{Night counts.}
#'     \item{`sol_*`, `tst_*`, `tbt_*`, `waso_*`}{Mean SOL, TST, TBT, WASO
#'       in minutes.}
#'     \item{`eff_*`}{Mean sleep efficiency (0–1).}
#'     \item{`bed_time_*`}{Mean bed-time in overnight decimal hours (may exceed
#'       24 for post-midnight onsets).}
#'     \item{`getup_time_*`}{Mean get-up time in plain decimal hours (0–24).}
#'     \item{`mid_sleep_*`}{Mean mid-sleep in decimal hours (0–24).}
#'   }
#'
#' @seealso [compute_cpd_metrics()], [run_pipeline_native()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- run_pipeline_native("recordings/P001.txt",
#'                               tz = "America/Sao_Paulo")
#'
#' sm <- compute_sleep_metrics(result$nights, tz = "America/Sao_Paulo")
#' sm$tbt_overall / 60   # mean TBT in hours
#' sm$mid_sleep_overall  # mean mid-sleep (decimal hours, 0-24)
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

  nd$is_free_day  <- .is_free_day(as.Date(nd$get_up_time, tz = tz), holidays)
  nd$bed_time_h   <- .to_h_bed(nd$bed_time, tz)
  nd$getup_time_h <- .to_h_plain(nd$get_up_time, tz)
  nd$mid_sleep_h  <- (nd$bed_time_h + nd$tbt / 60 / 2) %% 24

  .group_metrics <- function(sub, dtype) {
    pf <- function(v) paste0(v, "_", dtype)
    if (nrow(sub) == 0L) {
      r <- setNames(rep(NA_real_, 8L),
                    pf(c("sol", "tst", "tbt", "waso", "eff",
                         "bed_time", "getup_time", "mid_sleep")))
      return(c(setNames(0L, paste0("n_", dtype)), r))
    }
    setNames(
      c(nrow(sub),
        mean(sub$sol),          mean(sub$tst),
        mean(sub$tbt),          mean(sub$waso),
        mean(sub$eff),          mean(sub$bed_time_h),
        mean(sub$getup_time_h), mean(sub$mid_sleep_h)),
      c(paste0("n_", dtype),
        pf(c("sol", "tst", "tbt", "waso", "eff",
             "bed_time", "getup_time", "mid_sleep")))
    )
  }

  as.list(c(
    .group_metrics(nd,                    "overall"),
    .group_metrics(nd[!nd$is_free_day, ], "weekday"),
    .group_metrics(nd[ nd$is_free_day, ], "weekend")
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
#' Overnight wrapping is handled before averaging.
#'
#' **MSW** (mid-sleep on workdays) and **MSF** (mid-sleep on free days) use
#' plain arithmetic means of the per-night mid-sleep values. **SJL** is the
#' absolute difference between MSF and MSW.
#'
#' **MSFsc** (MSF corrected for sleep debt) adjusts MSF by the free-day sleep
#' onset and the weighted weekly mean sleep duration when free-day duration
#' exceeds weekday duration.
#'
#' **CPD** (Circadian Profile Deviation) is the RMS distance of each night's
#' mid-sleep from MSFsc in the time × sequence plane.
#'
#' @param nights A `tibble` of nightly sleep statistics as returned by
#'   [run_pipeline_native()] or [run_pipeline()]. Must contain at minimum the
#'   columns `is_nap`, `bed_time`, `get_up_time`, `tbt`, `sol`, `soi`.
#' @param min_tib_h `numeric(1)`. Minimum total in-bed time (hours) for
#'   inclusion. Default is `3.0`.
#' @param min_tib_eve_h `numeric(1)`. Minimum TBT (hours) for a night to
#'   qualify as a free-day-eve night (used for MSFsc). Default is `3.0`.
#' @param tz `character(1)`. Time zone for extracting clock hours. Default is
#'   `"UTC"`.
#' @param holidays A `Date` vector of public holidays to treat as free days.
#'   Default is `NULL` (weekends only).
#'
#' @return A named list:
#'   \describe{
#'     \item{`n_nights_cpd`, `n_free_days`, `n_workdays`}{Night counts.}
#'     \item{`msw_h`, `msw_hms`}{Mid-sleep on workdays (decimal hours; HH:MM).}
#'     \item{`msf_h`, `msf_hms`}{Mid-sleep on free days (decimal hours; HH:MM).}
#'     \item{`msfsc_h`, `msfsc_hms`}{MSF corrected for sleep debt (decimal
#'       hours; HH:MM).}
#'     \item{`sjl_h`, `sjl_min`}{Social jet lag (hours; minutes).}
#'     \item{`sjla_h`, `sjla_min`}{Signed SJL — MSF minus MSW (hours; minutes).}
#'     \item{`cpd_s`, `cpd_min`, `cpd_h`}{CPD (seconds; minutes; hours).}
#'   }
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
