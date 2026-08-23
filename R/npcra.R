#' Non-parametric circadian rhythm analysis (NPCRA)
#'
#' Computes the standard non-parametric circadian rhythm analysis variables
#' from an actigraphy recording. `IS` and `IV` are computed from the
#' hourly-mean activity profile (p = 24), matching pyActigraphy's actual
#' `_interdaily_stability()`/`_intradaily_variability()` (not the
#' population-variance formula in Gonçalves et al. 2014 or Van Someren et
#' al. 1999's own text -- the real implementation uses sample variance,
#' ddof = 1; see Details). `L5`/`M10` and their onsets
#' follow a different convention -- see below -- matching the notebook this
#' package's Vallim-pipeline comparisons were validated against.
#'
#' The following variables are computed:
#'
#' \describe{
#'   \item{`IS`}{**Interdaily stability** — consistency of the 24 h rest-activity
#'     pattern across days (range 0--1; higher = more stable). From the
#'     hourly-mean profile (p = 24).}
#'   \item{`IV`}{**Intradaily variability** — fragmentation of the
#'     rest-activity rhythm (>= 0; higher = more fragmented). From the
#'     hourly-mean profile (p = 24).}
#'   \item{`RA`}{**Relative amplitude** — contrast between the most active
#'     10 h window (M10) and least active 5 h window (L5) (range 0--1).}
#'   \item{`L5`}{Mean activity during the least active 5-hour window, found
#'     by a rolling mean over a 10-min-resampled series, searched globally
#'     across the whole recording (not the p = 24 hourly profile used for
#'     IS/IV).}
#'   \item{`L5_onset`}{Wall-clock time ("HH:MM") of the *end* of the
#'     least-active window -- the time of day, wrapped at 24 h regardless of
#'     which calendar day the window actually falls on.}
#'   \item{`M10`}{As `L5`, for the most active 10-hour window.}
#'   \item{`M10_onset`}{As `L5_onset`, for the most-active window.}
#' }
#'
#' @details
#' `IS`/`IV` build a 1h-resampled series `X` first: missing hourly bins get
#' a real zero (matching Python's `s_1h = s.resample('1h').mean().fillna(0)`),
#' not silent omission. `X` is grouped by hour-of-day into the p = 24 hourly
#' profile `Xh`. Both variables then use **sample variance** (divide by
#' n - 1, matching pandas' `.var()` default) rather than the population
#' variance (divide by n) that the classic Witting/Van Someren/Gonçalves
#' formulas describe on paper:
#' \deqn{IS = \frac{\sum_h(\bar{X}_h-\bar{X})^2/(p-1)}{\sum_i(X_i-\bar{X})^2/(N-1)}}
#' \deqn{IV = \frac{\sum_i(X_i-X_{i-1})^2/(N-1)}{\sum_i(X_i-\bar{X})^2/(N-1)}}
#' with `N` the number of hourly bins in the (zero-filled) recording and `p`
#' the number of hour-of-day groups present (24 for any recording spanning a
#' full day). The two formulas share the same denominator, matching
#' pyActigraphy's `d_1h = data.var()` being computed once and reused for both.
#'
#' @param x A `zeitr_recording` as returned by [read_actigraphy()], or a
#'   data frame / tibble with at least `datetime` and `activity` columns.
#'   If a `state` column is present, off-wrist epochs (`state == 4`) are
#'   used as-is by default -- see `exclude_offwrist`.
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
#'   Off-wrist exclusion (`state == 4`, if `exclude_offwrist = TRUE`) still
#'   applies either way; this does not replicate the Python pipeline's
#'   separate 30-min-threshold off-wrist-run rule for the M10/L5 windows
#'   specifically (see `exclude_offwrist`) -- only the D+1 window start.
#' @param exclude_offwrist `logical(1)`. If `TRUE`, off-wrist epochs
#'   (`state == 4`, when a `state` column is present) are deleted before
#'   computing any NPCRA variable. Default `FALSE` matches the actual
#'   production Python pipeline (Cell 16 of
#'   `vs_condor_py_pipeline_fix30_jrsv.ipynb`): `_nonparam_metrics()` is
#'   always called with `mask_series=None` for every variable (IS, IV, M10,
#'   L5, RA) -- off-wrist periods' raw device readings are used as-is, not
#'   deleted. Set `TRUE` for the more conservative (but non-Python-matching)
#'   behaviour of excluding them. This is a blunter tool than Python's own
#'   off-wrist handling for M10/L5 specifically (short runs zeroed and kept,
#'   long runs excluded via `NA` + `min_periods`), which this package does
#'   not replicate -- `TRUE` here simply deletes every off-wrist epoch
#'   outright, changing the time index rather than leaving gaps.
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
                          window_days = NULL, trim_to_d1 = TRUE,
                          exclude_offwrist = FALSE) {

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

  # ── Exclude off-wrist epochs, only if explicitly requested ─────────────
  # Default is exclude_offwrist = FALSE -- matches the actual production
  # Python (Cell 16 of vs_condor_py_pipeline_fix30_jrsv.ipynb): every NPCRA
  # variable (IS, IV, M10, L5, RA) comes from `_nonparam_metrics()`, always
  # called with `mask_series=None` ("M10/L5: sinal bruto sem mascara --
  # replica Condor ActStudio"). Off-wrist periods' raw device readings are
  # used as-is, not deleted. Excluding them by default (as this function
  # used to, unconditionally) was a reasonable general principle but not
  # what the validated pipeline actually does -- and for M10/L5
  # specifically, deleting epochs (creating index gaps) rather than
  # leaving the raw reading in place can relocate which window "wins" the
  # global min/max search entirely. Real-cohort validation against Julia's
  # Python reference showed catastrophic (near-zero or negative)
  # correlation on M10/L5 onset before this fix -- IS/IV were comparatively
  # unaffected, since summary statistics average over off-wrist stretches
  # rather than depending on a single winner-take-all window search.
  if (isTRUE(exclude_offwrist) && "state" %in% names(epochs)) {
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

#' Core NPCRA computation, matching pyActigraphy's actual IS/IV implementation
#'
#' All variables are computed on the zero-filled, hourly-resampled series
#' (p = 24). IS/IV use sample variance (ddof = 1) -- see [compute_npcra()]'s
#' Details for the exact formulas and why this differs from the classic
#' Witting/Van Someren/Gonçalves population-variance formula on paper.
#'
#' @noRd
.npcra_core <- function(datetimes, activity, epoch_s, L5_hours, M10_hours,
                        participant_id) {

  tz         <- attr(datetimes, "tzone") %||% "UTC"
  n_raw      <- length(activity)
  epochs_per_day <- 24 * 3600 / epoch_s

  # ── Build the 1h-resampled series ────────────────────────────────────────
  # Matches Python's `s_1h = s_isiv.resample('1h').mean().fillna(0)` (Cell 16):
  # missing hourly bins get a real zero, not silent omission. This is the
  # SAME series (and the same convention) `.lmx_window()` already applies
  # for M10/L5 at 10-min resolution -- IS/IV previously used tapply(), which
  # just dropped a missing hour from N/p entirely instead of zero-filling it.
  #
  # bin_times is built from sort(unique(hour_start)) directly -- NOT by
  # reconstructing timestamps from tapply()'s factor-level names via
  # as.POSIXct(names(...)). That round-trip is a real trap:
  # format.POSIXct() silently drops the time-of-day for exact midnight
  # ("2024-01-02" instead of "2024-01-02 00:00:00"), and tapply()'s default
  # character-conversion of a POSIXct grouping variable hits this for every
  # midnight bin. as.POSIXct() parsing that resulting MIXED vector of
  # date-only and full-datetime strings back in one vectorised call
  # mis-parsed multiple distinct days down to the same timestamp (silently,
  # no error) -- corrupting N and every IS/IV value derived from it.
  # Confirmed by direct execution: this bug alone took a 6-day fixture's
  # theoretical IS from 1.0362 down to 0.8782.
  hour_start  <- as.POSIXct(format(datetimes, "%Y-%m-%d %H:00:00", tz = tz), tz = tz)
  bin_times   <- sort(unique(hour_start))
  bin_means   <- vapply(bin_times, function(ht) mean(activity[hour_start == ht], na.rm = TRUE), numeric(1L))

  full_bins  <- seq(min(bin_times), max(bin_times), by = "hour")
  X          <- rep(0.0, length(full_bins))   # missing bins -> 0, matches .fillna(0)
  X[match(bin_times, full_bins)] <- as.double(bin_means)
  N          <- length(X)

  if (N < 2L) zeitr_abort("Fewer than 2 hourly slots after resampling.")

  slot_hour  <- as.integer(format(full_bins, "%H", tz = tz))

  # ── 24-h mean profile (p = 24) ──────────────────────────────────────────────
  Xm         <- mean(X, na.rm = TRUE)
  Xh         <- vapply(0L:23L, function(h) {
    vals <- X[slot_hour == h]
    if (length(vals) > 0L) mean(vals, na.rm = TRUE) else NA_real_
  }, numeric(1L))
  p          <- sum(!is.na(Xh))   # number of hour-of-day groups actually present
                                   # (== 24 for any recording spanning a full day)

  # ── IS / IV: sample variance (ddof = 1), matching pyActigraphy's actual ────
  # _interdaily_stability()/_intradaily_variability() (Cell 16), which call
  # pandas' .var() -- ddof = 1 (divide by n-1) by default. NOT the
  # population-variance (divide by n) formula their own docstrings and
  # Gonçalves et al. (2014) describe; the real production code uses ddof=1.
  # d_1h = Var(X) is the shared denominator for both IS and IV (matches
  # Python's `d_1h = data.var()`, computed once and reused for both).
  d_1h   <- sum((X  - Xm)^2, na.rm = TRUE) / (N - 1L)
  d_24h  <- sum((Xh - Xm)^2, na.rm = TRUE) / (p - 1L)
  c_1h   <- sum(diff(X)^2,   na.rm = TRUE) / (N - 1L)

  IS     <- if (d_1h > 0) d_24h / d_1h else NA_real_
  IV     <- if (d_1h > 0) c_1h  / d_1h else NA_real_

  # ── L5 and M10: rolling-mean search on a 10-min-resampled series ────────
  # Ported from the notebook's `_lmx_ow()`/`_nonparam_metrics()` (Cell 16 of
  # vs_condor_py_pipeline_fix29_jrsv.ipynb), which supersedes the shared
  # pipeline_functions.py's `_lmx()`/`compute_pyactigraphy_metrics()` and is
  # what actually produced the report's numbers. This differs from the
  # earlier Goncalves-profile approach in two ways: the value itself is a
  # rolling MEAN over 10-min bins (not the p = 24 hourly-profile mean), and
  # the search runs globally over the whole (D+1-trimmed) recording rather
  # than a single 24-h profile.
  l5_lmx    <- .lmx_window(datetimes, activity, L5_hours,  find_min = TRUE)
  m10_lmx   <- .lmx_window(datetimes, activity, M10_hours, find_min = FALSE)

  L5        <- l5_lmx$value
  M10       <- m10_lmx$value
  L5_onset  <- .format_time_of_day(l5_lmx$onset)
  M10_onset <- .format_time_of_day(m10_lmx$onset)

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
#' NOTE: no longer called from `.npcra_core()` -- `L5`/`M10`/onset are now
#' computed by `.lmx_window()` (rolling mean on a 10-min-resampled series,
#' matching the notebook's `_lmx_ow()`), which better matches the report's
#' actual production numbers. Left in place in case anything else in the
#' package or its tests calls it directly; safe to remove if not.
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

#' Locate the least/most active window via a rolling mean at 10-minute
#' resolution
#'
#' Ports the notebook's `_lmx_ow()` (Cell 16 of
#' `vs_condor_py_pipeline_fix29_jrsv.ipynb`), no-mask branch specifically --
#' `_nonparam_metrics()` always calls it with `mask_series=None` for M10/L5
#' ("replica Condor ActStudio"), so the off-wrist-run handling in `_lmx_ow()`
#' is intentionally not ported here. That branch is:
#'   r = series.resample('10min').mean().fillna(0)
#'   rolling = r.rolling(n_bins, min_periods=n_bins).mean()
#'   position = rolling.idxmin() (or idxmax()); value = rolling at that position
#' i.e. a rolling MEAN over 10-min bins (not a sum), searched globally over
#' the whole recording as received (`compute_npcra()`'s D+1 trim, if
#' enabled, has already been applied upstream).
#'
#' `min_periods = n_bins` means the first `n_bins - 1` positions (not yet a
#' full window) are excluded from the search -- mirrored here by starting
#' `ends` at `window_bins`, matching pandas' `rolling().dropna()`.
#'
#' @param datetimes POSIXct vector, native epoch resolution.
#' @param activity numeric vector, same length as `datetimes`.
#' @param period_hours numeric(1). Window width in hours (5 for L5, 10 for M10).
#' @param find_min logical(1). TRUE for L5 (least active), FALSE for M10.
#' @return list(onset = POSIXct, value = numeric). `onset` is the timestamp
#'   at the *end* of the winning window (pandas' `rolling().mean()` labels
#'   each value at the window's right edge) -- format with
#'   `.format_time_of_day()`, not as elapsed time.
#' @noRd
.lmx_window <- function(datetimes, activity, period_hours, find_min) {
  bin_min <- 10L
  bin_s   <- bin_min * 60

  t0      <- datetimes[1L]
  bin_idx <- as.integer(floor(as.numeric(difftime(datetimes, t0, units = "secs")) / bin_s))

  bin_mean <- tapply(activity, bin_idx, mean, na.rm = TRUE)
  bin_seq  <- as.integer(names(bin_mean))
  full_idx <- seq(min(bin_seq), max(bin_seq))
  means    <- rep(0.0, length(full_idx))   # missing bins -> 0, matches .fillna(0)
  means[match(bin_seq, full_idx)] <- as.double(bin_mean)

  window_bins <- as.integer(round(period_hours * 60 / bin_min))
  n <- length(means)
  if (window_bins > n) {
    return(list(onset = datetimes[length(datetimes)], value = NA_real_))
  }

  ends <- window_bins:n
  roll <- vapply(ends, function(e) mean(means[(e - window_bins + 1L):e]), numeric(1L))
  best <- if (find_min) which.min(roll) else which.max(roll)
  end_bin <- ends[best]

  # end of that bin, matching pandas' right-labelled rolling().mean()
  onset_time <- t0 + (full_idx[end_bin] + 1L) * bin_s

  list(onset = onset_time, value = roll[best])
}

#' Format a POSIXct timestamp as wall-clock "HH:MM" (time-of-day)
#'
#' Ports the notebook's onset formatting: `ts - ts.normalize()` in the
#' Python source gives elapsed time since midnight of `ts`'s OWN calendar
#' day (whichever day the winning window happens to land on) -- not the
#' recording's first day. That is what makes the Python onset wrap at 24 h
#' regardless of which day wins the global search; formatting `ts` directly
#' as a clock time here has the identical effect.
#' @noRd
.format_time_of_day <- function(ts) {
  lt <- as.POSIXlt(ts)
  sprintf("%02d:%02d", lt$hour, lt$min)
}

#' Convert an epoch-of-day index to "HH:MM" string
#' @noRd
.epochs_to_hhmm <- function(epoch_of_day, epoch_s) {
  total_seconds <- epoch_of_day * epoch_s
  h  <- as.integer(total_seconds %/% 3600) %% 24L
  m  <- as.integer((total_seconds %% 3600) %/% 60)
  sprintf("%02d:%02d", h, m)
}
