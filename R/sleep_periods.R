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

  # Critical: nap detection runs ONLY on currently-awake epochs (state == 0)
  # matching the Python pipeline's nap_bool mask
  is_wake <- x$state == 0L

  if (sum(is_wake) == 0L) return(x)

  wake_activity <- activity
  wake_activity[!is_wake] <- NA_real_

  # Adaptive median filter on wake epochs only
  pad_size         <- as.integer(round(epoch_h * pad_h))
  max_hws          <- as.integer(round(epoch_h * median_filter_h / 2))
  act_med_filtered <- .adaptive_median_filter(wake_activity, pad_size, max_hws)

  # Zero-proportion filter
  zp_filtered <- zero_prop_filter(activity, nap_zero_prop_hws, pad_value = 1)

  # Threshold criteria — only within wake epochs
  low_median   <- act_med_filtered < nap_median_thr
  high_zero_p  <- zp_filtered      > nap_zero_prop_thr

  if (use_and) {
    is_nap <- is_wake & low_median & high_zero_p
  } else {
    is_nap <- is_wake & (low_median & high_zero_p)
  }

  # state 7 for nap epochs, preserve existing non-wake states
  x$state <- ifelse(is_nap, 7L, x$state)
  x$sleep <- ifelse(x$state == 7L, 1L, x$sleep)

  x
}

# ── Internal Crespo helpers ───────────────────────────────────────────────────

#' Core Crespo MSP algorithm — direct translation of Python's detect_msp
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
  # Python:
  #   median_filter_window_size = int(epoch_hour * median_filter_window_hourly_length) + 1
  #   median_filter_half_window_size = int((median_filter_window_size - 1) / 2)
  #   pad_size = int(epoch_hour * adaptive_median_filter_padding_hourly_length)
  data_length     <- length(activity)
  mf_window_size  <- as.integer(epoch_h * median_filter_h) + 1L
  mf_hws          <- as.integer((mf_window_size - 1L) / 2L)
  pad_size        <- as.integer(epoch_h * pad_h)
  maximum_activity <- max(activity, na.rm = TRUE)

  # ── padded_activity (for adaptive filter) ─────────────────────────────────
  # Python: padded_activity = insert(activity, 0, pad) + append(pad)
  # 0-indexed array of length data_length + 2*pad_size
  padded_activity <- c(rep(maximum_activity, pad_size), activity,
                       rep(maximum_activity, pad_size))

  # ── Zero mitigation (pass 1) ───────────────────────────────────────────────
  # Python loops i in range(data_length), 0-indexed
  mitigated_zeros_activity <- activity
  zero_mitigation_activity_level <- as.double(
    stats::quantile(activity, zero_mitigation_q, names = FALSE))

  zero_sequence_length <- 0L
  for (i in seq_len(data_length)) {           # i is 1-indexed in R
    i0 <- i - 1L                               # 0-indexed equivalent
    if (activity[i] == 0) {
      zero_sequence_length <- zero_sequence_length + 1L
    } else {
      if (zero_sequence_length > consec_zeros_thr) {
        # Python: mitigated[i-run:i] += level   (0-indexed, exclusive end)
        # R equivalent: [(i-run) : (i-1)] (1-indexed)
        start_r <- i - zero_sequence_length
        end_r   <- i - 1L
        mitigated_zeros_activity[start_r:end_r] <-
          mitigated_zeros_activity[start_r:end_r] + zero_mitigation_activity_level
      }
      zero_sequence_length <- 0L
    }
  }

  # ── Coarse median filter ───────────────────────────────────────────────────
  # Python: median_filter(padded_mit, mf_hws, padding='padded', center=True)
  # padding='padded' means input is already padded; no additional padding added.
  # rolling_window(padded_mit, 2*mf_hws+1) produces (len-2*mf_hws) windows
  # filt = median(rolled)[mf_hws : data_length+mf_hws]  (0-indexed)
  # = R 1-indexed: [mf_hws+1 : data_length+mf_hws]
  coarse_pad <- as.integer(epoch_h * median_filter_h)
  padded_mit <- c(rep(maximum_activity, coarse_pad),
                  mitigated_zeros_activity,
                  rep(maximum_activity, coarse_pad))
  # Direct sliding window median — matches Python rolling_window exactly
  padded_mitigated_zeros_activity <- vapply(seq_len(data_length), function(i) {
    # Python median_filter(..., padding='padded') returns, for output position
    # k (0-indexed), the CENTRED median over original position k. With
    # coarse_pad = 2*mf_hws padding on each side, the centred window for R
    # position i (1-indexed) is padded_mit[(i + mf_hws):(i + 3*mf_hws)] —
    # exactly 2*mf_hws + 1 elements. The previous (i+1):(i+2*mf_hws+1) window
    # was shifted left by mf_hws - 1 (~4 h), corrupting morph_det and hence
    # which zeros were marked invalid.
    stats::median(padded_mit[(i + mf_hws):(i + 3L * mf_hws)])
  }, numeric(1))

  # ── Coarse threshold + morphological filter ────────────────────────────────
  sleep_median_activity_threshold <- as.double(
    stats::quantile(padded_mitigated_zeros_activity, sleep_quantile, names = FALSE))
  initial_sleep_detection <- as.integer(
    padded_mitigated_zeros_activity > sleep_median_activity_threshold)
  structuring_element <- rep(1L, morph_size)
  morphological_filtered_initial_detection <-
    .morphological_open_close(initial_sleep_detection, morph_size)

  # ── Mark invalid zeros (pass 2) ───────────────────────────────────────────
  # Python loops i in range(data_length), 0-indexed
  # invalid_zero_indexes are 0-indexed signal positions
  # then shifted by +pad_size to match padded_activity (also 0-indexed)
  # padded_activity[invalid_zero_indexes] = NaN
  invalid_zero_indexes <- integer(0)
  awake_zero_sequence_length <- 0L
  sleep_zero_sequence_length <- 0L

  for (i in seq_len(data_length)) {
    i0 <- i - 1L  # 0-indexed
    if (morphological_filtered_initial_detection[i] == 1L) {  # awake
      if (sleep_zero_sequence_length > sleep_zeros_thr) {
        # Python: range(i-run, i) 0-indexed
        run_start_0 <- i0 - sleep_zero_sequence_length
        run_end_0   <- i0 - 1L
        invalid_zero_indexes <- c(invalid_zero_indexes, run_start_0:run_end_0)
      }
      sleep_zero_sequence_length <- 0L
      if (activity[i] == 0) {
        awake_zero_sequence_length <- awake_zero_sequence_length + 1L
      } else {
        if (awake_zero_sequence_length > awake_zeros_thr) {
          run_start_0 <- i0 - awake_zero_sequence_length
          run_end_0   <- i0 - 1L
          invalid_zero_indexes <- c(invalid_zero_indexes, run_start_0:run_end_0)
        }
        awake_zero_sequence_length <- 0L
      }
    } else {  # sleep
      if (awake_zero_sequence_length > awake_zeros_thr) {
        run_start_0 <- i0 - awake_zero_sequence_length
        run_end_0   <- i0 - 1L
        invalid_zero_indexes <- c(invalid_zero_indexes, run_start_0:run_end_0)
      }
      awake_zero_sequence_length <- 0L
      if (activity[i] == 0) {
        sleep_zero_sequence_length <- sleep_zero_sequence_length + 1L
      } else {
        if (sleep_zero_sequence_length > sleep_zeros_thr) {
          run_start_0 <- i0 - sleep_zero_sequence_length
          run_end_0   <- i0 - 1L
          invalid_zero_indexes <- c(invalid_zero_indexes, run_start_0:run_end_0)
        }
        sleep_zero_sequence_length <- 0L
      }
    }
  }

  # Python: invalid_zero_indexes += pad_size  (still 0-indexed padded positions)
  # Python: padded_activity[invalid_zero_indexes] = NaN
  # R: convert 0-indexed to 1-indexed by adding 1
  if (length(invalid_zero_indexes) > 0L) {
    padded_idx_r <- unique(invalid_zero_indexes) + pad_size + 1L
    padded_activity[padded_idx_r] <- NA_real_
  }

  # ── Adaptive median filter ─────────────────────────────────────────────────
  # Python: loops i in range(data_length)
  #   center = i + pad_size  (0-indexed padded position)
  #   window = padded_activity[center-hws : center+hws+1]  (exclusive end)
  #   hws starts at minimum (pad_size), grows to maximum (mf_hws), then shrinks
  #   condition to grow: i < data_length - mf_hws + pad_size - 1
  adaptive_median_filtered_activity <- numeric(data_length)
  half_window_size <- pad_size  # starts at minimum

  for (i in seq_len(data_length)) {
    i0     <- i - 1L                    # 0-indexed
    center <- i0 + pad_size             # 0-indexed center in padded_activity
    # Python slice [center-hws : center+hws+1] is 0-indexed exclusive end
    # R equivalent: [(center-hws)+1 : center+hws] = [center-hws+1 : center+hws]
    lo_r <- center - half_window_size + 1L   # 1-indexed
    hi_r <- center + half_window_size + 1L   # 1-indexed (exclusive end +1, then +1 for R)
    lo_r <- max(1L, lo_r)
    hi_r <- min(length(padded_activity), hi_r)
    val  <- stats::median(padded_activity[lo_r:hi_r], na.rm = TRUE)

    if (is.nan(val) || is.na(val)) {
      val <- if (i0 > 0L) adaptive_median_filtered_activity[i - 1L] else 0
    }
    adaptive_median_filtered_activity[i] <- val

    # Python: if i < data_length - mf_hws + pad_size - 1: grow
    if (i0 < (data_length - mf_hws + pad_size - 1L)) {
      if (half_window_size < mf_hws) half_window_size <- half_window_size + 1L
    } else {
      if (half_window_size > pad_size) half_window_size <- half_window_size - 1L
    }
  }

  # ── Final threshold ────────────────────────────────────────────────────────
  thr1 <- as.double(stats::quantile(adaptive_median_filtered_activity,
                                    sleep_quantile, names = FALSE))
  if (thr1 < min_short_window_thr) thr1 <- min_short_window_thr

  # Python: improved_sleep_detection = where(filtered > threshold, 1, 0)
  # 1 = wake, 0 = sleep
  as.integer(adaptive_median_filtered_activity > thr1)
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
  closed  <- rolling_apply(dilated, hws, min, pad_value = 0L)

  # Erosion then dilation (opening)
  eroded  <- rolling_apply(closed, hws, min, pad_value = 0L)
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
