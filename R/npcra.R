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
#'   \item{`L5`}{Mean activity during the least active 5 consecutive hours.}
#'   \item{`L5_onset`}{Clock time of the L5 window onset (hh:mm).}
#'   \item{`M10`}{Mean activity during the most active 10 consecutive hours.}
#'   \item{`M10_onset`}{Clock time of the M10 window onset (hh:mm).}
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
                          window_days = NULL) {

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
  if (!is.null(epochs$state)) {
    keep      <- is.na(epochs$state) | epochs$state != 4L
    datetimes <- datetimes[keep]
    activity  <- activity[keep]
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
  L5_onset  <- sprintf("%02d:00", L5_result$onset_hour)
  M10       <- M10_result$value
  M10_onset <- sprintf("%02d:00", M10_result$onset_hour)

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

#' Convert an epoch-of-day index to "HH:MM" string
#' @noRd
.epochs_to_hhmm <- function(epoch_of_day, epoch_s) {
  total_seconds <- epoch_of_day * epoch_s
  h  <- as.integer(total_seconds %/% 3600) %% 24L
  m  <- as.integer((total_seconds %% 3600) %/% 60)
  sprintf("%02d:%02d", h, m)
}
