#' Non-parametric circadian rhythm analysis (NPCRA)
#'
#' Computes the standard non-parametric circadian rhythm analysis variables
#' from an actigraphy recording. All variables are derived from the 24-hour
#' activity profile following Van Someren et al. (1999) and Marler et al.
#' (2006).
#'
#' The following variables are computed:
#'
#' \describe{
#'   \item{`IS`}{**Interdaily stability** — consistency of the 24 h activity
#'     pattern across days (range 0--1; higher = more stable).}
#'   \item{`IV`}{**Intradaily variability** — fragmentation of the
#'     rest-activity rhythm (>= 0; higher = more fragmented).}
#'   \item{`RA`}{**Relative amplitude** — contrast between the most active
#'     10 h window (M10) and least active 5 h window (L5) (range 0--1).}
#'   \item{`L5`}{Mean activity during the least active 5 h window.}
#'   \item{`L5_onset`}{Clock time of the L5 window midpoint (hh:mm).}
#'   \item{`M10`}{Mean activity during the most active 10 h window.}
#'   \item{`M10_onset`}{Clock time of the M10 window midpoint (hh:mm).}
#' }
#'
#' @param x A `zeitr_recording` as returned by [read_actigraphy()], or a
#'   data frame / tibble with at least `datetime` and `activity` columns.
#' @param epoch_s `numeric(1)`. Epoch duration in seconds. If `NULL`
#'   (default), estimated automatically from the median inter-epoch interval.
#' @param L5_hours `numeric(1)`. Width of the least-active window in hours.
#'   Default is `5`.
#' @param M10_hours `numeric(1)`. Width of the most-active window in hours.
#'   Default is `10`.
#'
#' @return A one-row tibble with columns `participant_id`, `IS`, `IV`, `RA`,
#'   `L5`, `L5_onset`, `M10`, `M10_onset`, `n_days`, `n_epochs`.
#'
#' @references
#' Van Someren, E. J. W., Lijzenga, C., Mirmiran, M., & Swaab, D. F. (1997).
#' Long-term fitness training improves the circadian rest-activity rhythm in
#' healthy elderly males. *Journal of Biological Rhythms*, 12(2), 146--156.
#' \doi{10.1177/074873049701200206}
#'
#' Marler, M. R., Gehrman, P., Martin, J. L., & Ancoli-Israel, S. (2006).
#' The sigmoidally transformed cosine curve: a mathematical model for
#' circadian rhythms with symmetric non-sinusoidal shapes.
#' *Statistics in Medicine*, 25(22), 3893--3904.
#' \doi{10.1002/sim.2466}
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
#' @importFrom lubridate floor_date
#'
#' @examples
#' \dontrun{
#' rec   <- read_actigraphy("recordings/P001.txt")
#' npcra <- compute_npcra(rec)
#' npcra
#' }
compute_npcra <- function(x, epoch_s = NULL, L5_hours = 5, M10_hours = 10) {

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
  n         <- length(activity)

  if (n < 2L) zeitr_abort("Need at least 2 epochs to compute NPCRA.")

  # ── Epoch duration ───────────────────────────────────────────────────────────
  if (is.null(epoch_s)) {
    diffs   <- as.numeric(diff(datetimes), units = "secs")
    epoch_s <- stats::median(diffs[diffs > 0], na.rm = TRUE)
  }
  epochs_per_hour <- 3600 / epoch_s
  epochs_per_day  <- 24 * epochs_per_hour

  # ── Mean activity ────────────────────────────────────────────────────────────
  grand_mean <- mean(activity, na.rm = TRUE)
  if (grand_mean == 0) {
    zeitr_warn("All activity values are zero; NPCRA variables will be NA/NaN.")
  }

  # ── IS — Interdaily stability ────────────────────────────────────────────────
  # IS = (n / p) * sum_h(xh_bar - x_bar)^2 / sum_i(xi - x_bar)^2
  # where p = epochs per day, xh_bar = mean activity at time-of-day slot h
  #
  # floor_date must use the recording timezone so that slot 0 = local midnight.
  epoch_of_day <- as.integer(
    as.numeric(
      difftime(datetimes,
               lubridate::floor_date(datetimes, "day", tz = attr(datetimes, "tzone")),
               units = "secs")
    ) / epoch_s
  ) %% as.integer(round(epochs_per_day))

  p          <- as.integer(round(epochs_per_day))
  hourly_avg <- vapply(seq(0L, p - 1L), function(h) {
    idx <- epoch_of_day == h
    if (any(idx)) mean(activity[idx], na.rm = TRUE) else NA_real_
  }, numeric(1))

  numerator_IS   <- n / p * sum((hourly_avg - grand_mean)^2, na.rm = TRUE)
  denominator_IS <- sum((activity - grand_mean)^2, na.rm = TRUE)
  IS <- if (denominator_IS > 0) numerator_IS / denominator_IS else NA_real_

  # ── IV — Intradaily variability ──────────────────────────────────────────────
  # IV = n * sum_i(xi - x_{i-1})^2 / ((n-1) * sum_i(xi - x_bar)^2)
  diffs_activity <- diff(activity)
  numerator_IV   <- n * sum(diffs_activity^2, na.rm = TRUE)
  denominator_IV <- (n - 1L) * sum((activity - grand_mean)^2, na.rm = TRUE)
  IV <- if (denominator_IV > 0) numerator_IV / denominator_IV else NA_real_

  # ── L5 and M10 ───────────────────────────────────────────────────────────────
  L5_result  <- .rolling_window_mean(activity, epoch_of_day, p,
                                     window_hours = L5_hours,  epoch_s = epoch_s,
                                     find_min = TRUE)
  M10_result <- .rolling_window_mean(activity, epoch_of_day, p,
                                     window_hours = M10_hours, epoch_s = epoch_s,
                                     find_min = FALSE)

  L5       <- L5_result$value
  L5_onset <- .epochs_to_hhmm(L5_result$onset_epoch, epoch_s)
  M10      <- M10_result$value
  M10_onset <- .epochs_to_hhmm(M10_result$onset_epoch, epoch_s)

  # ── RA — Relative amplitude ──────────────────────────────────────────────────
  RA <- if (!is.na(M10) && !is.na(L5) && (M10 + L5) > 0) {
    (M10 - L5) / (M10 + L5)
  } else {
    NA_real_
  }

  n_days <- n / epochs_per_day

  tibble::tibble(
    participant_id = participant_id,
    IS             = round(IS,  4),
    IV             = round(IV,  4),
    RA             = round(RA,  4),
    L5             = round(L5,  4),
    L5_onset       = L5_onset,
    M10            = round(M10, 4),
    M10_onset      = M10_onset,
    n_days         = round(n_days, 2),
    n_epochs       = n
  )
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' Compute rolling window mean across all days and find min or max window
#'
#' The activity series is folded into a 24-hour profile (average activity at
#' each epoch-of-day position), then a sliding window of `window_hours` is
#' applied to find the least active (L5) or most active (M10) window.
#'
#' @param activity numeric vector
#' @param epoch_of_day integer vector of time-of-day slot for each epoch (0-based)
#' @param epochs_per_day integer; number of epochs in one 24-hour period
#' @param window_hours window width in hours
#' @param epoch_s epoch duration in seconds
#' @param find_min logical; TRUE for L5, FALSE for M10
#' @noRd
.rolling_window_mean <- function(activity, epoch_of_day, epochs_per_day,
                                 window_hours, epoch_s, find_min) {
  epochs_per_hour <- 3600 / epoch_s
  window_size     <- as.integer(round(window_hours * epochs_per_hour))

  if (window_size >= epochs_per_day) {
    return(list(value = mean(activity, na.rm = TRUE), onset_epoch = 0L))
  }

  # Build 24-hour average profile keyed by actual time-of-day slot
  profile <- vapply(seq(0L, epochs_per_day - 1L), function(h) {
    idx <- epoch_of_day == h
    if (any(idx)) mean(activity[idx], na.rm = TRUE) else 0
  }, numeric(1))

  # Wrap profile for circular sliding window
  wrapped <- c(profile, profile)

  window_means <- vapply(seq(0L, epochs_per_day - 1L), function(start) {
    mean(wrapped[(start + 1L):(start + window_size)])
  }, numeric(1))

  if (find_min) {
    onset <- which.min(window_means) - 1L
    value <- window_means[onset + 1L]
  } else {
    onset <- which.max(window_means) - 1L
    value <- window_means[onset + 1L]
  }

  list(value = value, onset_epoch = onset)
}

#' Convert an epoch-of-day index to "HH:MM" string
#' @noRd
.epochs_to_hhmm <- function(epoch_of_day, epoch_s) {
  total_seconds <- epoch_of_day * epoch_s
  h  <- as.integer(total_seconds %/% 3600) %% 24L
  m  <- as.integer((total_seconds %% 3600) %/% 60)
  sprintf("%02d:%02d", h, m)
}
