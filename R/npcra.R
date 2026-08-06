#' Non-parametric circadian rhythm analysis (NPCRA)
#'
#' Computes the standard non-parametric circadian rhythm analysis variables
#' from an actigraphy recording, following Gonçalves et al. (2014) and Van
#' Someren et al. (1999). All variables are derived from the 24-hour average
#' activity profile built from **hourly means** (p = 24).
#'
#' The following variables are computed:
#'
#' \describe{
#'   \item{`IS`}{**Interdaily stability** — consistency of the 24 h rest-activity
#'     pattern across days (range 0--1; higher = more stable).}
#'   \item{`IV`}{**Intradaily variability** — fragmentation of the
#'     rest-activity rhythm (>= 0; higher = more fragmented).}
#'   \item{`RA`}{**Relative amplitude** — contrast between the most active
#'     10 h window (M10) and least active 5 h window (L5) (range 0--1).}
#'   \item{`L5`}{Mean activity during the least active 5 consecutive hours
#'     (from the 24-h mean profile).}
#'   \item{`L5_onset`}{Elapsed time ("H:MM:SS", not wrapped at 24 h) from the
#'     first recorded day's midnight to the end of the least-active window,
#'     located on a 10-min-resampled series -- mirrors the Python
#'     reference's `_lmx()`/`_td_format()` exactly. NOT the same clock as
#'     `L5` above, which stays on the coarser hourly profile.}
#'   \item{`M10`}{Mean activity during the most active 10 consecutive hours
#'     (from the 24-h mean profile).}
#'   \item{`M10_onset`}{As `L5_onset`, for the most-active window.}
#' }
#'
#' @param x A `zeitr_recording` as returned by [read_actigraphy()], or a
#'   data frame / tibble with at least `datetime` and `activity` columns.
#'   If a `state` column is present, off-wrist epochs (`state == 4`) are
#'   excluded before computing all NPCRA variables.
#' @param epoch_s `numeric(1)`. Epoch duration in seconds. If `NULL`
#'   (default), estimated automatically from the median inter-epoch interval.
#' @param L5_hours `numeric(1)`. Width of the least-active window in hours.
#'   Default is `5`.
#' @param M10_hours `numeric(1)`. Width of the most-active window in hours.
#'   Default is `10`.
#' @param window_days `numeric(1)` or `NULL`. If supplied, the recording is
#'   split into non-overlapping windows of this length (in days) and NPCRA
#'   variables are computed for each window. A `window_start` column is added
#'   to the output. Partial final windows (shorter than `window_days`) are
#'   included but flagged via a lower `n_days` value. Default `NULL` computes
#'   a single estimate over the full recording.
#' @param trim_to_d1 `logical(1)`. If `TRUE` (default), the recording is
#'   trimmed to start at 00:00 of D+1 -- the first full calendar day after
#'   recording onset -- before any NPCRA variable is computed, matching the
#'   Python reference pipeline's convention (it always starts its NPCRA
#'   window at D+1 00:00 rather than spanning the raw, typically fractional,
#'   recording length). Set to `FALSE` for the full untrimmed recording (the
#'   pre-`trim_to_d1` behaviour). If trimming would leave fewer than 2 epochs,
#'   a warning is emitted and the untrimmed recording is used instead.
#'   Off-wrist exclusion (`state == 4`) still applies either way; this does
#'   not replicate the Python pipeline's separate 30-min-threshold rule for
#'   the M10/L5 windows specifically -- only the D+1 window start.
#'
#' @return A tibble with columns `participant_id`, `window_start` (if
#'   `window_days` is set), `IS`, `IV`, `RA`, `L5`, `L5_onset`, `M10`,
#'   `M10_onset`, `n_days`, `n_epochs`.
#'
#' @references
#' Gonçalves, B. S. B., Adamowicz, T., Louzada, F. M., Moreno, C. R., &
#' Araujo, J. F. (2014). A fresh look at the use of nonparametric analysis in
#' actimetry. *Sleep Medicine Reviews*, 20, 84--91.
#' \doi{10.1016/j.smrv.2014.06.002}
#'
#' Van Someren, E. J. W., Swaab, D. F., Colenda, C. C., Cohen, W.,
#' McCall, W. V., & Rosenquist, P. B. (1999). Bright light therapy:
#' Improved sensitivity to its effects on rest-activity rhythms in
#' Alzheimer patients by application of nonparametric methods.
#' *Chronobiology International*, 16(4), 505--518.
#' \doi{10.3109/07420529908998724}
#'
#' @export
#'
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' rec   <- read_actigraphy("recordings/P001.txt")
#'
#' # Single estimate over the full recording
#' compute_npcra(rec)
#'
#' # Per-fortnight estimates
#' compute_npcra(rec, window_days = 14)
#' }
compute_npcra <- function(x, epoch_s = NULL, L5_hours = 5, M10_hours = 10,
                          window_days = NULL, trim_to_d1 = TRUE) {

  # ── Extract epochs tibble and participant_id ─────────────────────────────────
  if (inherits(x, "zeitr_recording")) {
    epochs         <- x$epochs
    participant_id <- x$metadata$participant_id %||% NA_character_
  } else if (is.data.frame(x)) {
    epochs         <- x
    participant_id <- NA_character_
  } else {
    zeitr_abort("{.arg x} must be a {.cls zeitr_recording} or a data frame.")
  }

  required <- c("datetime", "activity")
  missing  <- setdiff(required, names(epochs))
  if (length(missing) > 0L) {
    zeitr_abort("Missing required column(s): {.val {missing}}")
  }

  datetimes <- as.POSIXct(epochs$datetime)
  activity  <- as.double(epochs$activity)

  # ── Exclude off-wrist epochs if state column is present ──────────────────────
  if ("state" %in% names(epochs)) {
    keep      <- is.na(epochs$state) | epochs$state != 4L
    datetimes <- datetimes[keep]
    activity  <- activity[keep]
  }

  # ── Trim to D+1 00:00 (Python pipeline convention) ────────────────────
  # Only attempted when there are already >= 2 epochs -- an input that's
  # already too small shouldn't get a "trimming left too few" warning of
  # its own before hitting the real n < 2 abort below.
  if (isTRUE(trim_to_d1) && length(datetimes) >= 2L) {
    tz_d1     <- attr(datetimes, "tzone") %||% "UTC"
    first_day <- as.Date(min(datetimes), tz = tz_d1)
    d1_start  <- as.POSIXct(paste0(format(first_day + 1L), " 00:00:00"), tz = tz_d1)
    keep_d1   <- datetimes >= d1_start
    if (sum(keep_d1) < 2L) {
      zeitr_warn(
        "{.arg trim_to_d1} = TRUE leaves fewer than 2 epochs after trimming to D+1 00:00; using the untrimmed recording instead."
      )
    } else {
      datetimes <- datetimes[keep_d1]
      activity  <- activity[keep_d1]
    }
  }

  n <- length(activity)
  if (n < 2L) zeitr_abort("Need at least 2 epochs to compute NPCRA.")

  # ── Epoch duration ───────────────────────────────────────────────────────────
  if (is.null(epoch_s)) {
    diffs   <- as.numeric(diff(datetimes), units = "secs")
    epoch_s <- stats::median(diffs[diffs > 0], na.rm = TRUE)
  }

  epochs_per_day <- 24 * 3600 / epoch_s

  # ── Windowed vs full-recording mode ─────────────────────────────────────────
  if (!is.null(window_days)) {
    epochs_per_window <- as.integer(round(window_days * epochs_per_day))
    starts <- seq(1L, n, by = epochs_per_window)
    rows <- lapply(starts, function(s) {
      e   <- min(s + epochs_per_window - 1L, n)
      row <- .npcra_core(datetimes[s:e], activity[s:e],
                         epoch_s, L5_hours, M10_hours, participant_id)
      tibble::add_column(row, window_start = as.Date(datetimes[s]), .after = "participant_id")
    })
    return(do.call(rbind, rows))
  }

  .npcra_core(datetimes, activity, epoch_s, L5_hours, M10_hours, participant_id)
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' Core NPCRA computation (per Gonçalves et al. 2014)
#'
#' All variables are computed on the hourly-clustered series (p = 24),
#' matching the formulas in Equations (1) and (2) of the reference paper.
#'
#' @noRd
.npcra_core <- function(datetimes, activity, epoch_s, L5_hours, M10_hours,
                        participant_id) {

  tz         <- attr(datetimes, "tzone") %||% "UTC"
  n_raw      <- length(activity)
  epochs_per_day <- 24 * 3600 / epoch_s

  # ── Cluster epochs into hourly means ────────────────────────────────────────
  # date × hour-of-day gives each epoch a unique slot; mean within slot.
  local_date  <- as.Date(format(datetimes, "%Y-%m-%d", tz = tz))
  hour_of_day <- as.integer(format(datetimes, "%H", tz = tz))

  slot_key    <- paste(local_date, sprintf("%02d", hour_of_day))
  slot_means  <- tapply(activity, slot_key, mean, na.rm = TRUE)

  # Recover hour_of_day for each slot (last two chars of key)
  slot_names  <- names(slot_means)
  slot_hour   <- as.integer(substr(slot_names, nchar(slot_names) - 1L, nchar(slot_names)))
  X           <- as.double(slot_means)   # hourly series, length N
  N           <- length(X)

  if (N < 2L) zeitr_abort("Fewer than 2 hourly slots after clustering.")

  # ── 24-h mean profile (p = 24) ──────────────────────────────────────────────
  p          <- 24L
  Xm         <- mean(X, na.rm = TRUE)
  Xh         <- vapply(0L:(p - 1L), function(h) {
    vals <- X[slot_hour == h]
    if (length(vals) > 0L) mean(vals, na.rm = TRUE) else NA_real_
  }, numeric(1L))

  # ── IS (Equation 2, Gonçalves 2014) ─────────────────────────────────────────
  # IS = (N/p) * sum_h(Xh - Xm)^2 / sum_i(Xi - Xm)^2
  IS_num <- (N / p) * sum((Xh - Xm)^2, na.rm = TRUE)
  IS_den <- sum((X  - Xm)^2, na.rm = TRUE)
  IS     <- if (IS_den > 0) IS_num / IS_den else NA_real_

  # ── IV (Equation 1, Gonçalves 2014) ─────────────────────────────────────────
  # IV = N * sum_i(Xi - Xi-1)^2 / ((N-1) * sum_i(Xi - Xm)^2)
  IV_num <- N * sum(diff(X)^2, na.rm = TRUE)
  IV_den <- (N - 1L) * IS_den
  IV     <- if (IV_den > 0) IV_num / IV_den else NA_real_

  # ── L5 and M10 from the 24-h mean profile ───────────────────────────────────
  L5_result  <- .rolling_window_profile(Xh, L5_hours,  find_min = TRUE)
  M10_result <- .rolling_window_profile(Xh, M10_hours, find_min = FALSE)

  L5        <- L5_result$value
  M10       <- M10_result$value

  # Onset resolution: the *value* above stays on the p = 24 hourly profile
  # (Goncalves et al. 2014), but the onset clock time is located on a
  # 10-min-resampled series -- a faithful port of the Python reference's
  # `_lmx()` -- so it isn't locked to whole hours the way the hourly profile
  # would force it to be.
  l5_lmx    <- .lmx_window(datetimes, activity, L5_hours,  find_min = TRUE)
  m10_lmx   <- .lmx_window(datetimes, activity, M10_hours, find_min = FALSE)
  day0      <- as.POSIXct(format(as.Date(datetimes[1L], tz = tz), "%Y-%m-%d 00:00:00"), tz = tz)
  L5_onset  <- .format_elapsed_hms(l5_lmx$onset  - day0)
  M10_onset <- .format_elapsed_hms(m10_lmx$onset - day0)

  # ── RA ───────────────────────────────────────────────────────────────────────
  RA <- if (!is.na(M10) && !is.na(L5) && (M10 + L5) > 0) {
    (M10 - L5) / (M10 + L5)
  } else {
    NA_real_
  }

  tibble::tibble(
    participant_id = participant_id,
    IS             = round(IS,  4),
    IV             = round(IV,  4),
    RA             = round(RA,  4),
    L5             = round(L5,  4),
    L5_onset       = L5_onset,
    M10            = round(M10, 4),
    M10_onset      = M10_onset,
    n_days         = round(n_raw / epochs_per_day, 2),
    n_epochs       = n_raw
  )
}

#' Find least/most active window from the 24-h mean profile
#'
#' @param profile numeric(24) — hourly mean activity profile (hours 0--23)
#' @param window_hours integer window width in hours
#' @param find_min logical; TRUE for L5, FALSE for M10
#' @noRd
.rolling_window_profile <- function(profile, window_hours, find_min) {
  p    <- length(profile)   # 24
  w    <- as.integer(round(window_hours))
  if (w >= p) return(list(value = mean(profile, na.rm = TRUE), onset_hour = 0L))

  # Circular sliding window over the 24-h profile
  wrapped      <- c(profile, profile)
  window_means <- vapply(seq(0L, p - 1L), function(start) {
    mean(wrapped[(start + 1L):(start + w)], na.rm = TRUE)
  }, numeric(1L))

  onset <- if (find_min) which.min(window_means) - 1L else which.max(window_means) - 1L
  list(value = window_means[onset + 1L], onset_hour = as.integer(onset))
}

#' Locate the least/most active window at 10-minute resolution
#'
#' Faithful port of the Python reference's `_lmx()`: resamples `activity`
#' into 10-minute bins (summed), then slides a trailing window of
#' `period_hours` width across the *entire* recording (as received --
#' `compute_npcra()`'s own D+1 trim, if enabled, has already been applied
#' upstream), picking the window whose sum is smallest (L5) or largest
#' (M10). This runs independently of the hourly-profile (p = 24) value
#' calculation in `.rolling_window_profile()` -- it exists solely to locate
#' the onset at the same minute resolution `pyActigraphy`'s `_lmx()` uses,
#' rather than being locked to whole hours.
#'
#' Mirrors `_lmx()`'s `idx` exactly: the returned `onset` is the timestamp
#' at the *end* of the winning window (pandas' `rolling().sum()` labels
#' each value at the window's right edge), not the window's start. This
#' has NOT yet been validated against a `python_output`-derived fixture --
#' do that before trusting `L5_onset`/`M10_onset` in production. In
#' particular, `_td_format()` in the Python source formats elapsed time
#' from day-zero midnight and does NOT wrap at 24 h, so if the global
#' extremum window falls on a later day, the elapsed hours here can
#' legitimately exceed 24 -- if the resulting values look implausibly
#' large, check whether the actual production code (possibly patched
#' inline in `vs_condor_py_pipeline_fix29_jrsv.ipynb` rather than in
#' `pipeline_functions_fix27.py`) restricts the search window further.
#'
#' @param datetimes POSIXct vector, native epoch resolution.
#' @param activity numeric vector, same length as `datetimes`.
#' @param period_hours numeric(1). Window width in hours (5 for L5, 10 for M10).
#' @param find_min logical(1). TRUE for L5 (least active), FALSE for M10.
#' @return list(onset = POSIXct, value = numeric).
#' @noRd
.lmx_window <- function(datetimes, activity, period_hours, find_min) {
  bin_min <- 10L
  bin_s   <- bin_min * 60

  t0      <- datetimes[1L]
  bin_idx <- as.integer(floor(as.numeric(difftime(datetimes, t0, units = "secs")) / bin_s))

  binned   <- tapply(activity, bin_idx, sum, na.rm = TRUE)
  bin_seq  <- as.integer(names(binned))
  full_idx <- seq(min(bin_seq), max(bin_seq))
  vals     <- rep(0.0, length(full_idx))
  vals[match(bin_seq, full_idx)] <- as.double(binned)

  window_bins <- as.integer(round(period_hours * 60 / bin_min))
  n <- length(vals)
  if (window_bins > n) {
    return(list(onset = datetimes[length(datetimes)], value = NA_real_))
  }

  ends <- window_bins:n
  roll <- vapply(ends, function(e) sum(vals[(e - window_bins + 1L):e]), numeric(1L))
  best <- if (find_min) which.min(roll) else which.max(roll)
  end_bin <- ends[best]

  # end of that bin, matching pandas' right-labelled rolling().sum()
  onset_time <- t0 + (full_idx[end_bin] + 1L) * bin_s

  list(onset = onset_time, value = roll[best])
}

#' Format an elapsed time (difftime) as "H:MM:SS", NOT wrapped at 24 h
#'
#' Port of the Python reference's `_td_format()`.
#' @noRd
.format_elapsed_hms <- function(elapsed) {
  total <- as.numeric(elapsed, units = "secs")
  total <- as.integer(round(abs(total)))
  h <- total %/% 3600L
  r <- total %% 3600L
  m <- r %/% 60L
  s <- r %% 60L
  sprintf("%d:%02d:%02d", h, m, s)
}

#' Convert an epoch-of-day index to "HH:MM" string
#' @noRd
.epochs_to_hhmm <- function(epoch_of_day, epoch_s) {
  total_seconds <- epoch_of_day * epoch_s
  h  <- as.integer(total_seconds %/% 3600) %% 24L
  m  <- as.integer((total_seconds %% 3600) %/% 60)
  sprintf("%02d:%02d", h, m)
}
