#' Detect main sleep periods using the Crespo algorithm
#'
#' Identifies the main sleep period(s) in an actigraphy recording using the
#' algorithm described in Crespo et al. (2012). The method applies an adaptive
#' median filter to the activity signal, mitigates spuriously long zero runs,
#' and thresholds the result at a quantile of the filtered signal. Morphological
#' closing and opening operations are then used to smooth the binary
#' sleep/wake estimate.
#'
#' @param x A tibble as returned by [detect_offwrist_bimodal()] (or
#'   [prepare_actigraphy()] if off-wrist detection is skipped), containing
#'   columns `datetime`, `activity`, and `state`.
#' @param epoch_h `numeric(1)`. Number of epochs per hour. If `NULL`
#'   (default), estimated automatically from the median inter-epoch interval
#'   in `datetime`.
#' @param median_filter_h `numeric(1)`. Length of the preprocessing median
#'   filter window in hours. Default is `8`.
#' @param pad_h `numeric(1)`. Padding length in hours added before the
#'   adaptive median filter. Default is `1`.
#' @param sleep_quantile `numeric(1)`. Quantile of the filtered activity used
#'   as the sleep/wake threshold. Default is `1/3`.
#' @param morph_size `integer(1)`. Size of the structuring element used in
#'   morphological closing/opening. Default is `61` epochs.
#' @param consec_zeros_thr `integer(1)`. Runs of zeros longer than this
#'   threshold are treated as invalid (zero mitigation). Default is `15`.
#' @param awake_zeros_thr `integer(1)`. Threshold for consecutive zeros within
#'   wake periods. Default is `2`.
#' @param sleep_zeros_thr `integer(1)`. Threshold for consecutive zeros within
#'   sleep periods. Default is `30`.
#' @param zero_mitigation_q `numeric(1)`. Quantile of activity used to
#'   determine the mitigation level for invalid zero runs. Default is `0.33`.
#' @param min_short_window_thr `numeric(1)`. Minimum value of the adaptive
#'   median threshold; if the fitted quantile falls below this, the threshold
#'   is clamped here. Default is `1.0`.
#'
#' @return The input tibble `x` with `state` and `sleep` columns updated.
#'   Sleep epochs have `state == 1` and `sleep == 1`; off-wrist epochs
#'   (`state == 4`) are preserved and excluded from the sleep column.
#'
#' @references
#' Crespo, C., Aboy, M., Fernández, J. R., & Mojón, A. (2012). Automatic
#' identification of activity-rest periods based on actigraphy.
#' *Journal of Medical and Biological Engineering*, 32(4), 249–256.
#' \doi{10.5405/jmbe.1033}
#'
#' @seealso [detect_naps_crespo()] for secondary sleep period (nap) detection.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rec  <- read_acttrust("recordings/P001.txt")
#' prep <- prepare_actigraphy(rec)
#' prep <- detect_offwrist_bimodal(prep)
#' prep <- detect_sleep_crespo(prep)
#' }
detect_sleep_crespo <- function(
    x,
    epoch_h            = NULL,
    median_filter_h    = 8,
    pad_h              = 1,
    sleep_quantile     = 1 / 3,
    morph_size         = 61L,
    consec_zeros_thr   = 15L,
    awake_zeros_thr    = 2L,
    sleep_zeros_thr    = 30L,
    zero_mitigation_q  = 0.33,
    min_short_window_thr = 1.0
) {
  .check_crespo_input(x)

  activity  <- as.double(x$activity)
  epoch_h   <- epoch_h %||% .estimate_epoch_h(x$datetime)

  out <- .crespo_msp(
    activity             = activity,
    epoch_h              = epoch_h,
    median_filter_h      = median_filter_h,
    pad_h                = pad_h,
    sleep_quantile       = sleep_quantile,
    morph_size           = morph_size,
    consec_zeros_thr     = consec_zeros_thr,
    awake_zeros_thr      = awake_zeros_thr,
    sleep_zeros_thr      = sleep_zeros_thr,
    zero_mitigation_q    = zero_mitigation_q,
    min_short_window_thr = min_short_window_thr
  )

  # out: 0 = sleep, 1 = wake  → invert for state (1 = sleep, 0 = wake)
  state_new        <- 1L - out
  x$state          <- ifelse(x$state == 4L, 4L, state_new)
  x$sleep          <- ifelse(x$state == 4L, 0L, state_new)

  x
}

#' Detect secondary sleep periods (naps) using the Crespo algorithm
#'
#' Identifies naps — secondary sleep periods — in an actigraphy recording
#' using the nap variant of the Crespo algorithm. This function should be
#' run *after* [detect_sleep_crespo()].
#'
#' The nap algorithm combines two criteria: a low rolling median activity
#' threshold and a high zero-activity proportion around each epoch. Epochs
#' satisfying either (or both, if `use_and = TRUE`) criteria are scored as nap
#' sleep.
#'
#' @inheritParams detect_sleep_crespo
#' @param nap_median_thr `numeric(1)`. Epochs with a rolling median activity
#'   below this value may be scored as nap sleep. Default is `2.0`.
#' @param nap_zero_prop_thr `numeric(1)`. Epochs with a rolling zero-activity
#'   proportion above this threshold may be scored as nap sleep. Default is
#'   `0.5`.
#' @param nap_zero_prop_hws `integer(1)`. Half-window size (epochs) for the
#'   rolling zero-proportion filter. Default is `5L`.
#' @param use_and `logical(1)`. If `TRUE`, both the median activity AND
#'   zero-proportion criteria must be met. If `FALSE` (default), either
#'   criterion is sufficient.
#'
#' @return The input tibble `x` with `state` and `sleep` columns updated.
#'   Nap epochs have `state == 7` and are merged into `sleep` as value `1`.
#'
#' @references
#' Crespo, C., Aboy, M., Fernández, J. R., & Mojón, A. (2012). Automatic
#' identification of activity-rest periods based on actigraphy.
#' *Journal of Medical and Biological Engineering*, 32(4), 249–256.
#' \doi{10.5405/jmbe.1033}
#'
#' @seealso [detect_sleep_crespo()] for main sleep period detection.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rec  <- read_acttrust("recordings/P001.txt")
#' prep <- prepare_actigraphy(rec)
#' prep <- detect_offwrist_bimodal(prep)
#' prep <- detect_sleep_crespo(prep)
#' prep <- detect_naps_crespo(prep)
#' }
detect_naps_crespo <- function(
    x,
    epoch_h           = NULL,
    median_filter_h   = 8,
    pad_h             = 1,
    nap_median_thr    = 2.0,
    nap_zero_prop_thr = 0.5,
    nap_zero_prop_hws = 5L,
    use_and           = FALSE
) {
  .check_crespo_input(x)

  activity <- as.double(x$activity)
  epoch_h  <- epoch_h %||% .estimate_epoch_h(x$datetime)

  # Adaptive median filter (same as MSP but no zero mitigation)
  pad_size         <- as.integer(round(epoch_h * pad_h))
  max_hws          <- as.integer(round(epoch_h * median_filter_h / 2))
  act_med_filtered <- .adaptive_median_filter(activity, pad_size, max_hws)

  # Zero-proportion filter
  zp_filtered <- zero_prop_filter(activity, nap_zero_prop_hws, pad_value = 1)

  # Threshold criteria
  low_median   <- act_med_filtered < nap_median_thr
  high_zero_p  <- zp_filtered      > nap_zero_prop_thr

  if (use_and) {
    is_nap <- low_median & high_zero_p
  } else {
    is_nap <- low_median | high_zero_p
  }

  # out: 0 = nap, 1 = wake  → state 7 for nap epochs
  state_nap <- ifelse(is_nap, 7L, x$state)
  x$state   <- ifelse(x$state == 4L, 4L, state_nap)
  x$sleep   <- ifelse(x$state == 7L, 1L, x$sleep)

  x
}

# ── Internal Crespo helpers ───────────────────────────────────────────────────

#' Core Crespo MSP algorithm
#' @noRd
.crespo_msp <- function(
    activity,
    epoch_h,
    median_filter_h,
    pad_h,
    sleep_quantile,
    morph_size,
    consec_zeros_thr,
    awake_zeros_thr,
    sleep_zeros_thr,
    zero_mitigation_q,
    min_short_window_thr
) {
  n            <- length(activity)
  pad_size     <- as.integer(round(epoch_h * pad_h))
  max_hws      <- as.integer(round(epoch_h * median_filter_h / 2))
  max_activity <- max(activity, na.rm = TRUE)

  # ── Zero mitigation ────────────────────────────────────────────────────────
  mitigated    <- activity
  zmit_level   <- as.double(stats::quantile(activity, zero_mitigation_q,
                                            names = FALSE))

  run_len <- 0L
  for (i in seq_len(n)) {
    if (activity[i] == 0) {
      run_len <- run_len + 1L
    } else {
      if (run_len > consec_zeros_thr) {
        mitigated[(i - run_len):(i - 1L)] <-
          mitigated[(i - run_len):(i - 1L)] + zmit_level
      }
      run_len <- 0L
    }
  }

  # ── Pad mitigated signal for initial coarse median filter ─────────────────
  pad_mf  <- as.integer(round(epoch_h * median_filter_h))
  padded_mit <- c(rep(max_activity, pad_mf), mitigated, rep(max_activity, pad_mf))
  coarse_med <- median_filter(padded_mit, as.integer(pad_mf / 2))
  coarse_med <- coarse_med[(pad_mf + 1L):(pad_mf + n)]

  # Initial threshold + morphological filter for invalid zero identification
  thr0          <- as.double(stats::quantile(coarse_med, sleep_quantile,
                                             names = FALSE))
  initial_det   <- as.integer(coarse_med > thr0)
  morph_det     <- .morphological_open_close(initial_det, morph_size)

  # Mark invalid zeros based on morphological result
  padded_act <- c(rep(max_activity, pad_size), activity, rep(max_activity, pad_size))

  run_len     <- 0L
  awake_run   <- 0L
  sleep_run   <- 0L
  invalid_idx <- integer(0)

  for (i in seq_len(n)) {
    in_sleep <- morph_det[i] == 0L

    if (!in_sleep) {  # wake region
      if (sleep_run > sleep_zeros_thr) {
        invalid_idx <- c(invalid_idx, (i - sleep_run):(i - 1L) + pad_size)
      }
      sleep_run <- 0L
      if (activity[i] == 0) {
        awake_run <- awake_run + 1L
      } else {
        if (awake_run > awake_zeros_thr) {
          invalid_idx <- c(invalid_idx, (i - awake_run):(i - 1L) + pad_size)
        }
        awake_run <- 0L
      }
    } else {  # sleep region
      if (awake_run > awake_zeros_thr) {
        invalid_idx <- c(invalid_idx, (i - awake_run):(i - 1L) + pad_size)
      }
      awake_run <- 0L
      if (activity[i] == 0) {
        sleep_run <- sleep_run + 1L
      } else {
        if (sleep_run > sleep_zeros_thr) {
          invalid_idx <- c(invalid_idx, (i - sleep_run):(i - 1L) + pad_size)
        }
        sleep_run <- 0L
      }
    }
  }

  if (length(invalid_idx) > 0L) {
    padded_act[unique(invalid_idx)] <- NA_real_
  }

  # ── Adaptive median filter ─────────────────────────────────────────────────
  act_med_filtered <- .adaptive_median_filter(padded_act, pad_size, max_hws,
                                              padded = TRUE, n_orig = n)

  # ── Final threshold ────────────────────────────────────────────────────────
  thr1 <- as.double(stats::quantile(act_med_filtered, sleep_quantile,
                                    names = FALSE))
  if (thr1 < min_short_window_thr) thr1 <- min_short_window_thr

  final_det <- as.integer(act_med_filtered > thr1)  # 1 = wake, 0 = sleep

  final_det
}

#' Adaptive median filter with variable window size
#' @noRd
.adaptive_median_filter <- function(
    x,
    pad_size,
    max_hws,
    padded  = FALSE,
    n_orig  = NULL
) {
  if (padded) {
    n      <- n_orig %||% (length(x) - 2L * pad_size)
    center_start <- pad_size + 1L
  } else {
    n      <- length(x)
    padded_x <- c(rep(x[1], pad_size), x, rep(x[length(x)], pad_size))
    x        <- padded_x
    center_start <- pad_size + 1L
  }

  out     <- numeric(n)
  hws     <- pad_size

  for (i in seq_len(n)) {
    center <- center_start + i - 1L
    lo     <- max(1L, center - hws)
    hi     <- min(length(x), center + hws)
    val    <- stats::median(x[lo:hi], na.rm = TRUE)

    if (is.na(val)) {
      val <- if (i > 1L) out[i - 1L] else 0
    }

    out[i] <- val

    # Variable window: grow toward max_hws then shrink near the end
    if (i < (n - max_hws + pad_size)) {
      if (hws < max_hws) hws <- hws + 1L
    } else {
      if (hws > pad_size) hws <- hws - 1L
    }
  }

  out
}

#' Binary morphological close then open (1-D, flat structuring element)
#' @param x integer vector (0/1)
#' @param size integer structuring element size (must be odd)
#' @noRd
.morphological_open_close <- function(x, size) {
  hws <- as.integer((size - 1L) / 2L)

  # Dilation then erosion (closing)
  dilated <- rolling_apply(x, hws, max, pad_value = 0L)
  closed  <- rolling_apply(dilated, hws, min, pad_value = 1L)

  # Erosion then dilation (opening)
  eroded  <- rolling_apply(closed, hws, min, pad_value = 1L)
  opened  <- rolling_apply(eroded, hws, max, pad_value = 0L)

  as.integer(round(opened))
}

# ── Shared input check ────────────────────────────────────────────────────────

#' @noRd
.check_crespo_input <- function(x) {
  required <- c("activity", "datetime", "state")
  missing  <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    zeitr_abort(
      "{.arg x} is missing required column(s): {.val {missing}}.
       Did you run {.fn prepare_actigraphy} first?"
    )
  }
}

#' Estimate epoch_h (epochs per hour) from datetime column
#' @noRd
.estimate_epoch_h <- function(datetimes) {
  if (length(datetimes) < 2L) zeitr_abort("Need at least 2 epochs to estimate epoch duration.")
  deltas_s <- as.numeric(diff(as.POSIXct(datetimes)), units = "secs")
  epoch_s  <- stats::median(deltas_s[deltas_s > 0], na.rm = TRUE)
  3600 / epoch_s
}
