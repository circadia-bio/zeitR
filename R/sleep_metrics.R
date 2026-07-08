# ── Sleep summary and CPD/SJL metrics ─────────────────────────────────────────
# R port of compute_sleep_metrics() and compute_cpd_metrics() from
# pipeline_functions_fix27.py (Julia Ribeiro da Silva Vallim, 2024).
# Column names and derived metrics mirror the Python output format exactly.
#
# Parity fixes applied in this file:
#   fix-locale  : .is_free_day() uses format(d, "%u") (ISO 8601 weekday number)
#                 instead of weekdays() to avoid locale-dependent day names
#                 (e.g. "sabado" on pt_BR vs "Saturday" on en_US).
#   fix-circular: MSF and MSW use .mean_circ_h() (circular mean) instead of
#                 plain mean(), matching calculate_msf/msw in the fix29 notebook.
#   fix-truncate: compute_cpd_metrics() drops episodes starting after noon on
#                 the last recording day (truncated by end of file), matching
#                 the nights_to_df() filter added in the fix29 notebook.

# ── Internal helpers ──────────────────────────────────────────────────────────

# Mirrors _series_to_h: plain decimal hours 0-24
.to_h_plain <- function(x, tz) {
  x <- as.POSIXct(x, tz = tz)
  as.numeric(format(x, "%H", tz = tz)) +
    as.numeric(format(x, "%M", tz = tz)) / 60 +
    as.numeric(format(x, "%S", tz = tz)) / 3600
}

# Mirrors _mean_circular_h (fix27/fix29): values >= 12h shifted by -24 before
# averaging, then wrapped to [0, 24).  Handles midnight wrap without
# trigonometry; equivalent to fix29's sin/cos approach for typical sleep times.
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

# Parse free_days to a vector of ISO 8601 weekday integers (Mon=1 ... Sun=7).
# Accepts English day names (case-insensitive) or integers 1-7.
.parse_free_days <- function(free_days) {
  if (is.null(free_days) || length(free_days) == 0L)
    return(c(6L, 7L))  # default: Saturday + Sunday

  day_map <- c(
    monday = 1L, tuesday = 2L, wednesday = 3L, thursday = 4L,
    friday = 5L, saturday = 6L, sunday = 7L
  )

  if (is.numeric(free_days) || is.integer(free_days)) {
    wdays <- as.integer(free_days)
    bad   <- wdays[wdays < 1L | wdays > 7L]
    if (length(bad) > 0L)
      zeitr_abort("Day numbers must be between 1 (Monday) and 7 (Sunday); got: {bad}.")
    return(wdays)
  }

  wdays <- day_map[tolower(trimws(as.character(free_days)))]
  bad   <- free_days[is.na(wdays)]
  if (length(bad) > 0L)
    zeitr_abort(c(
      "Unrecognised day name(s): {.val {bad}}.",
      "i" = 'Use English day names (e.g. {.val Saturday}) or ISO integers 1-7 (Mon=1 ... Sun=7).'
    ))
  unname(wdays)
}

# is_free_day applied to get-up date.
# Uses ISO 8601 weekday number via format(d, "%u") instead of weekdays() to
# avoid locale-dependent day names: weekdays() returns "Saturday" on en_US but
# "sabado" on pt_BR, silently marking no day as a free day on non-English
# systems.  ISO 8601: Mon=1, Tue=2, ..., Sat=6, Sun=7.
#
# free_wdays: pre-parsed integer vector from .parse_free_days().
#
# holidays may contain:
#   - Date objects or "YYYY-MM-DD" strings  -> exact year-specific match
#   - "DD-MM" strings                       -> recurring annual match
# Mixing the two forms in the same vector is supported.
.is_free_day <- function(dates, holidays, free_wdays) {
  d           <- as.Date(dates)
  wd          <- as.integer(format(d, "%u"))  # ISO 8601: Mon=1 ... Sat=6, Sun=7
  is_free_wday <- wd %in% free_wdays

  if (is.null(holidays)) return(is_free_wday)

  h_chr        <- as.character(holidays)
  is_recurring <- grepl("^\\d{2}-\\d{2}$", h_chr)  # matches "DD-MM" only

  is_in_holiday <- rep(FALSE, length(d))

  if (any(!is_recurring)) {
    specific <- suppressWarnings(as.Date(h_chr[!is_recurring]))
    specific <- specific[!is.na(specific)]
    if (length(specific) > 0L)
      is_in_holiday <- is_in_holiday | (d %in% specific)
  }

  if (any(is_recurring)) {
    d_mmdd        <- format(d, "%d-%m")
    is_in_holiday <- is_in_holiday | (d_mmdd %in% h_chr[is_recurring])
  }

  is_free_wday | is_in_holiday
}

# is_free_day_eve: the bed-time date + 1 day is a free day
.is_free_day_eve <- function(bts, tz, holidays, free_wdays) {
  d_next <- as.Date(as.POSIXct(bts, tz = tz), tz = tz) + 1L
  .is_free_day(d_next, holidays, free_wdays)
}

# Emit a warning when no holidays are supplied.
# Suppressed globally with options(zeitR.no_holidays_warn = FALSE).
.warn_no_holidays <- function(fn) {
  if (!isFALSE(getOption("zeitR.no_holidays_warn", TRUE))) {
    zeitr_warn(c(
      paste0(fn, ": no {.arg holidays} supplied."),
      "i" = "Only the days in {.arg free_days} are treated as free days.",
      "i" = "If you ran {.fn run_pipeline_native} with holidays, pass {.code holidays = result$holidays}.",
      "i" = "Suppress with {.code options(zeitR.no_holidays_warn = FALSE)}."
    ))
  }
}


# ── compute_sleep_metrics ─────────────────────────────────────────────────────

#' Compute sleep metrics split by day type
#'
#' Calculates mean sleep metrics separately for all nights, weekday nights,
#' and weekend/holiday nights. Output column names and derived metrics mirror
#' `compute_sleep_metrics()` from Julia Ribeiro da Silva Vallim's
#' `pipeline_functions_fix27.py`.
#'
#' Dispatches on the class of `x`:
#' * **`data.frame` / `tibble`** -- treated as a nights table (same behaviour
#'   as the original function signature).
#' * **`zeitr_result`** -- extracts `x$nights` and inherits `x$holidays` and
#'   `x$free_days` automatically.
#'
#' @details
#' A night is assigned to the **free-day** group when the get-up date falls on
#' one of the days in `free_days` or matches an entry in `holidays`. Day-of-week
#' is determined via the ISO 8601 weekday number (`format(date, "%u")`) rather
#' than `weekdays()`, which is locale-dependent. All remaining nights are
#' **workday** nights.
#'
#' `sleep_onset_h` uses a circular mean (values >= 12 h are shifted by -24 h
#' before averaging, then wrapped to [0, 24)). `sleep_offset_h` uses a plain
#' arithmetic mean.
#'
#' `fps_h` (free period sleep) equals `fpr_tib_h - (latencia_min + inertia_min)
#' / 60` -- TBT net of sleep onset latency and sleep inertia.
#'
#' `dp_midsleep_min` and `dp_tst_min` are standard deviations of per-night
#' mid-sleep (minutes) and TST (minutes), respectively.
#'
#' @param x A `zeitr_result` object **or** a `tibble` of nightly sleep
#'   statistics as returned by [run_pipeline_native()] or [run_pipeline()].
#'   Must contain at minimum the columns `is_nap`, `bed_time`, `get_up_time`,
#'   `tbt`, `tst`, `sol`, `soi`, `waso`, `eff`.
#' @param min_tib_h `numeric(1)`. Minimum total in-bed time (hours) for a night
#'   to be included. Default is `5.0` (matching the Python reference).
#' @param tz `character(1)`. Time zone for extracting clock hours from
#'   timestamps. Default is `"UTC"`.
#' @param holidays Holidays to treat as free days in addition to the days in
#'   `free_days`. Accepts three forms, which can be mixed in the same vector:
#'   * `Date` objects or `"YYYY-MM-DD"` strings for year-specific dates (e.g.
#'     `as.Date("2019-03-04")` for one Carnival day).
#'   * `"DD-MM"` strings for dates that recur every year (e.g. `"25-12"` for
#'     Christmas, `"01-01"` for New Year).
#'   Default is `NULL`. When `x` is a `zeitr_result`, defaults to
#'   `x$holidays` automatically. A warning is emitted when `NULL`; suppress
#'   with `options(zeitR.no_holidays_warn = FALSE)`.
#' @param free_days A character vector of day names (`"Monday"` through
#'   `"Sunday"`, case-insensitive) or ISO integers (1 = Monday ... 7 = Sunday)
#'   identifying which days of the week are unconditionally treated as free
#'   days. Default is `c("Saturday", "Sunday")`. When `x` is a `zeitr_result`,
#'   defaults to `x$free_days` automatically.
#'
#' @return A named list with metrics for three groups (`overall`, `wd` =
#'   workday, `fd` = free day):
#'   \describe{
#'     \item{`n_overall`, `n_wd`, `n_fd`}{Night counts.}
#'     \item{`sleep_onset_h`, `sleep_offset_h`}{Circular mean onset and
#'       arithmetic mean offset in decimal hours.}
#'     \item{`fpr_tib_h`}{Mean TBT in hours.}
#'     \item{`fps_h`}{Mean free period sleep (TBT - SOL - SOI) in hours.}
#'     \item{`tst_h`}{Mean TST in hours.}
#'     \item{`latencia_min`, `inertia_min`}{Mean SOL and SOI in minutes.}
#'     \item{`waso_min`}{Mean WASO in minutes.}
#'     \item{`sleep_eff_pct`}{Mean sleep efficiency in percent (0-100).}
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
#' # From a zeitr_result: holidays and free_days forwarded automatically
#' result <- run_pipeline_native("recordings/P001.txt",
#'                               tz        = "America/Sao_Paulo",
#'                               holidays  = my_holidays,
#'                               free_days = c("Saturday", "Sunday"))
#' sm <- compute_sleep_metrics(result, tz = "America/Sao_Paulo")
#'
#' # Non-standard schedule: Friday + Saturday as free days
#' sm <- compute_sleep_metrics(result$nights,
#'                             tz        = "America/Sao_Paulo",
#'                             holidays  = my_holidays,
#'                             free_days = c("Friday", "Saturday"))
#'
#' sm$tst_h            # mean TST in hours
#' sm$sleep_onset_h    # mean sleep onset (circular, decimal hours)
#' sm$dp_midsleep_min  # within-person SD of mid-sleep in minutes
#' }
compute_sleep_metrics <- function(x, ...) UseMethod("compute_sleep_metrics")


#' @rdname compute_sleep_metrics
#' @export
compute_sleep_metrics.zeitr_result <- function(x,
                                                min_tib_h = 5.0,
                                                tz        = "UTC",
                                                holidays  = x$holidays,
                                                free_days = x$free_days,
                                                ...) {
  if (is.null(free_days)) free_days <- c("Saturday", "Sunday")
  compute_sleep_metrics.default(x$nights,
                                 min_tib_h = min_tib_h,
                                 tz        = tz,
                                 holidays  = holidays,
                                 free_days = free_days,
                                 ...)
}


#' @rdname compute_sleep_metrics
#' @export
compute_sleep_metrics.default <- function(x,
                                           min_tib_h = 5.0,
                                           tz        = "UTC",
                                           holidays  = NULL,
                                           free_days = c("Saturday", "Sunday"),
                                           ...) {
  if (is.null(holidays)) .warn_no_holidays("compute_sleep_metrics()")
  free_wdays <- .parse_free_days(free_days)

  nd <- x[!x$is_nap & x$tbt / 60 >= min_tib_h, ]
  if (nrow(nd) == 0L) {
    zeitr_warn("No nights remain after TIB filter (min_tib_h = {min_tib_h} h).")
    return(list())
  }

  nd$is_free_day <- .is_free_day(as.Date(nd$get_up_time, tz = tz), holidays, free_wdays)
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
      stats::sd(mid_sub, na.rm = TRUE) * 60,   # dp_midsleep_min
      stats::sd(sub$tst, na.rm = TRUE)         # dp_tst_min
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
#' Ribeiro da Silva Vallim's `pipeline_functions_fix27.py`, with two updates
#' to match the fix29 notebook:
#' * **MSF and MSW** are computed with the circular mean to correctly handle
#'   mid-sleep times that wrap midnight.
#' * **Truncated episodes** starting after noon on the last recording day are
#'   excluded before metric computation.
#'
#' Dispatches on the class of `x`:
#' * **`data.frame` / `tibble`** -- treated as a nights table (same behaviour
#'   as the original function signature).
#' * **`zeitr_result`** -- extracts `x$nights` and inherits `x$holidays` and
#'   `x$free_days` automatically.
#'
#' @details
#' Mid-sleep is computed per night as:
#' \deqn{\text{MS} = \left(\text{SO} + \frac{\text{offset} - \text{onset}}{2}\right) \bmod 24}
#' where onset = bts + SOL and offset = gts - SOI (both in decimal hours).
#'
#' **MSW** and **MSF** use the circular mean of per-night mid-sleep values
#' (values >= 12 h shifted by -24 before averaging, then wrapped to [0, 24)),
#' matching `calculate_msf()` / `calculate_msw()` from the fix29 notebook.
#' **MSFsc** adjusts MSF by the free-day-eve sleep onset and the weighted
#' weekly mean sleep duration when free-day duration exceeds weekday duration.
#' **CPD** is the RMS distance of each night's mid-sleep from MSFsc in the
#' time x sequence plane.
#'
#' @param x A `zeitr_result` object **or** a `tibble` of nightly sleep
#'   statistics as returned by [run_pipeline_native()] or [run_pipeline()].
#' @param min_tib_h `numeric(1)`. Minimum TBT (hours) for inclusion. Default
#'   is `3.0`.
#' @param min_tib_eve_h `numeric(1)`. Minimum TBT (hours) for a night to
#'   qualify as a free-day-eve night. Default is `3.0`.
#' @param tz `character(1)`. Time zone for extracting clock hours. Default is
#'   `"UTC"`.
#' @param holidays Holidays to treat as free days in addition to the days in
#'   `free_days`. Accepts three forms, which can be mixed in the same vector:
#'   * `Date` objects or `"YYYY-MM-DD"` strings for year-specific dates (e.g.
#'     `as.Date("2019-03-04")` for one Carnival day).
#'   * `"DD-MM"` strings for dates that recur every year (e.g. `"25-12"` for
#'     Christmas, `"01-01"` for New Year).
#'   Default is `NULL`. When `x` is a `zeitr_result`, defaults to
#'   `x$holidays` automatically. A warning is emitted when `NULL`; suppress
#'   with `options(zeitR.no_holidays_warn = FALSE)`.
#' @param free_days A character vector of day names (`"Monday"` through
#'   `"Sunday"`, case-insensitive) or ISO integers (1 = Monday ... 7 = Sunday)
#'   identifying which days of the week are unconditionally treated as free
#'   days. Default is `c("Saturday", "Sunday")`. When `x` is a `zeitr_result`,
#'   defaults to `x$free_days` automatically.
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
#' # From a zeitr_result: holidays and free_days forwarded automatically
#' result <- run_pipeline_native("recordings/P001.txt",
#'                               tz        = "America/Sao_Paulo",
#'                               holidays  = my_holidays,
#'                               free_days = c("Saturday", "Sunday"))
#' cpd <- compute_cpd_metrics(result, tz = "America/Sao_Paulo")
#'
#' # Non-standard schedule: Friday + Saturday as free days
#' cpd <- compute_cpd_metrics(result$nights,
#'                            tz        = "America/Sao_Paulo",
#'                            holidays  = my_holidays,
#'                            free_days = c("Friday", "Saturday"))
#'
#' cpd$sjl_min    # social jet lag in minutes
#' cpd$msf_hms    # mid-sleep on free days as HH:MM
#' cpd$cpd_min    # CPD in minutes
#' }
compute_cpd_metrics <- function(x, ...) UseMethod("compute_cpd_metrics")


#' @rdname compute_cpd_metrics
#' @export
compute_cpd_metrics.zeitr_result <- function(x,
                                              min_tib_h     = 3.0,
                                              min_tib_eve_h = 3.0,
                                              tz            = "UTC",
                                              holidays      = x$holidays,
                                              free_days     = x$free_days,
                                              ...) {
  if (is.null(free_days)) free_days <- c("Saturday", "Sunday")
  compute_cpd_metrics.default(x$nights,
                               min_tib_h     = min_tib_h,
                               min_tib_eve_h = min_tib_eve_h,
                               tz            = tz,
                               holidays      = holidays,
                               free_days     = free_days,
                               ...)
}


#' @rdname compute_cpd_metrics
#' @export
compute_cpd_metrics.default <- function(x,
                                         min_tib_h     = 3.0,
                                         min_tib_eve_h = 3.0,
                                         tz            = "UTC",
                                         holidays      = NULL,
                                         free_days     = c("Saturday", "Sunday"),
                                         ...) {
  if (is.null(holidays)) .warn_no_holidays("compute_cpd_metrics()")
  free_wdays <- .parse_free_days(free_days)

  nd <- x[!x$is_nap & x$tbt / 60 >= min_tib_h, ]
  if (nrow(nd) == 0L)
    zeitr_abort("No nocturnal episodes remain after TIB filter (min_tib_h = {min_tib_h} h).")

  # fix29: exclude episodes starting after noon on the last recording day --
  # these are truncated by end-of-file, not genuine complete nights.
  last_day  <- as.Date(max(as.POSIXct(nd$get_up_time, tz = tz)), tz = tz)
  truncated <- as.POSIXct(nd$bed_time, tz = tz) >=
    as.POSIXct(last_day, tz = tz) + 12 * 3600
  if (any(truncated)) {
    n_trunc <- sum(truncated)
    zeitr_warn(
      "Excluded {n_trunc} truncated episode(s) starting after noon on the last recording day ({last_day})."
    )
    nd <- nd[!truncated, ]
    if (nrow(nd) == 0L)
      zeitr_abort("No episodes remain after truncated-episode filter.")
  }

  onset_h  <- .to_h_plain(nd$bed_time,    tz) + nd$sol / 60
  offset_h <- .to_h_plain(nd$get_up_time, tz) - nd$soi / 60

  df <- data.frame(
    sleep_onset     = onset_h,
    sleep_offset    = offset_h,
    time_in_bed     = nd$tbt / 60,
    mid_sleep       = .midsleep(onset_h, offset_h),
    is_free_day     = .is_free_day(as.Date(nd$get_up_time, tz = tz), holidays, free_wdays),
    is_free_day_eve = .is_free_day_eve(nd$bed_time, tz, holidays, free_wdays)
  )
  df$is_free_day_eve[df$time_in_bed < min_tib_eve_h] <- FALSE

  we  <- df$is_free_day
  wd  <- !df$is_free_day
  eve <- df$is_free_day_eve

  if (!any(we))  zeitr_abort("No free days found. Supply public holidays via {.arg holidays}.")
  if (!any(wd))  zeitr_abort("No workdays found.")
  if (!any(eve)) zeitr_abort("No free-day-eve nights found.")

  # MSF and MSW: circular mean to handle midnight wrap correctly.
  # Mirrors calculate_msf() / calculate_msw() from the fix29 notebook.
  msf_h <- .mean_circ_h(df$mid_sleep[we])
  msw_h <- .mean_circ_h(df$mid_sleep[wd])

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
