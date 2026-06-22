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
#'   columns `datetime`, `activity`, and `state`. The detector runs on the
#'   on-wrist subset (`state != 4`) only, mirroring the Python `cspd_wrapper`.
#' @param epoch_h `numeric(1)`. Number of epochs per hour. If `NULL`
#'   (default), derived from the epoch duration (the mode of the on-wrist
#'   inter-epoch interval), as `3600 / duration`.
#' @param median_filter_h `numeric(1)`. Length of the preprocessing median
#'   filter window in hours. Default is `8`.
#' @param pad_h `numeric(1)`. Padding length in hours added before the
#'   adaptive median filter. Default is `1`.
#' @param sleep_quantile `numeric(1)`. Quantile of the filtered activity used
#'   as the MSP sleep/wake threshold. Default is `0.365` (the ActTrust CSPD
#'   value used by `cspd_wrapper`; the standalone Crespo algorithm uses `1/3`).
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
#' @param refine `logical(1)`. If `TRUE` (default), the MSP detection is
#'   refined into final sleep periods by the CSPD bed-time / get-up-time
#'   refiners (`.cspd_refine_periods`), reproducing the Python `refined_output`.
#'   If `FALSE`, the raw MSP detection is used directly.
#' @param condition `integer(1)`. Initial condition flag (default `0`). The MSP
#'   stage bumps it to `2` when its activity-median threshold clamps to
#'   `min_short_window_thr`; the refiner uses that effective condition. Affects
#'   only `refine = TRUE`.
#'
#' @return The input tibble `x` with `state` and `sleep` columns updated.
#'   Sleep epochs have `state == 1` and `sleep == 1`; off-wrist epochs
#'   (`state == 4`) are preserved and excluded from the sleep column. With
#'   `refine = TRUE` the sleep epochs delimit the refined sleep PERIODS (the
#'   Python `refined_output`); per-epoch wake/sleep within them is scored later
#'   by [compute_waso()].
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
    sleep_quantile     = 0.365,
    morph_size         = 61L,
    consec_zeros_thr   = 15L,
    awake_zeros_thr    = 2L,
    sleep_zeros_thr    = 30L,
    zero_mitigation_q  = 0.33,
    min_short_window_thr = 1.0,
    refine             = TRUE,
    condition          = 0L
) {
  .check_crespo_input(x)

  state_in <- as.integer(x$state)
  onwrist  <- state_in != 4L
  if (!any(onwrist)) {
    x$sleep <- rep(0L, nrow(x))
    return(x)
  }

  ow_activity <- as.double(x$activity[onwrist])
  ow_datetime <- x$datetime[onwrist]

  # Epoch duration = mode of the on-wrist inter-epoch gaps; epoch_h is derived
  # from it (matches CSPD.model, which sets epoch_hour = 3600 / duration).
  secs     <- as.numeric(as.POSIXct(ow_datetime))
  dd       <- diff(secs)
  duration <- if (length(dd) > 0L)
    as.numeric(names(sort(table(dd), decreasing = TRUE))[1]) else 60
  epoch_h  <- epoch_h %||% (3600 / duration)

  # CSPD.model runs detect_msp on the SCALED activity (activity / duration).
  # The detection is scale-invariant except for the absolute
  # min_short_window_thr clamp, so the scaling is required to reproduce the
  # CSPD-internal MSP that feeds the refiner.
  scaled_activity <- ow_activity / duration

  msp_res <- .crespo_msp(
    activity             = scaled_activity,
    epoch_h              = epoch_h,
    median_filter_h      = median_filter_h,
    pad_h                = pad_h,
    sleep_quantile       = sleep_quantile,
    morph_size           = morph_size,
    consec_zeros_thr     = consec_zeros_thr,
    awake_zeros_thr      = awake_zeros_thr,
    sleep_zeros_thr      = sleep_zeros_thr,
    zero_mitigation_q    = zero_mitigation_q,
    min_short_window_thr = min_short_window_thr,
    condition            = condition
  )
  msp <- msp_res$detection
  # detect_msp bumps condition to 2 when its activity-median threshold clamps;
  # the refiner uses that effective condition, not the initial one.
  eff_condition <- msp_res$condition

  # Refine the MSP into final sleep PERIODS with the CSPD bed/get-up refiners.
  if (isTRUE(refine)) {
    detection <- .cspd_refine_periods(
      activity        = ow_activity,
      datetime_stamps = ow_datetime,
      msp_detection   = msp,
      condition       = eff_condition
    )$refined_output
  } else {
    detection <- msp
  }

  # detection: 1 = wake, 0 = sleep  ->  state 1 = sleep, 0 = wake; 4 = off-wrist.
  new_state          <- state_in
  new_state[onwrist] <- ifelse(detection == 0, 1L, 0L)
  x$state            <- new_state
  x$sleep            <- ifelse(new_state == 4L, 0L, as.integer(new_state == 1L))

  x
}

#' Detect secondary sleep periods (naps) using the Crespo nap algorithm
#'
#' Faithful port of the Python `nap_wrapper`: runs the full CSPD model in nap
#' mode (`detect_naps = TRUE`) on the currently-awake epochs (`state == 0`) and
#' merges the detected naps into the sleep state. Must be run *after*
#' [detect_sleep_crespo()].
#'
#' Nap detection uses a nap-mode MSP (a high zero-proportion combined with a
#' low adaptive-median activity, `.crespo_nap_msp()`) followed by the same
#' bed-time / get-up-time refiners as the main sleep detection, with the nap
#' parameter set (`.cspd_nap_params()`) and nap-specific minimum-length
#' post-processing. Detected naps are written as `state == 1` (merged into
#' "sleep"), matching `nap_wrapper`, which assigns
#' `state[wake] = 1 - refined_output` (i.e. naps are *not* a distinct state).
#'
#' @param x A tibble as returned by [detect_sleep_crespo()], containing columns
#'   `datetime`, `activity`, and `state`. Nap detection runs on the wake
#'   (`state == 0`) subset only, mirroring the Python `nap_bool` mask.
#' @param epoch_h `numeric(1)`. Number of epochs per hour. If `NULL` (default),
#'   derived from the wake-subsequence epoch duration as `3600 / duration`.
#' @param params CSPD nap configuration list (default `.cspd_nap_params()`),
#'   the port of `nap_wrapper`'s parameter set.
#'
#' @return The input tibble `x` with `state` and `sleep` columns updated. Nap
#'   epochs become `state == 1` and `sleep == 1`; off-wrist (`state == 4`) and
#'   existing main-sleep epochs are preserved.
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
detect_naps_crespo <- function(x, epoch_h = NULL, params = .cspd_nap_params()) {
  .check_crespo_input(x)
  p <- params

  state_in <- as.integer(x$state)
  is_wake  <- state_in == 0L            # Python nap_bool = (state == 0)
  if (!any(is_wake)) return(x)

  wake_activity <- as.double(x$activity[is_wake])
  wake_datetime <- x$datetime[is_wake]

  # Epoch duration (mode of the wake-subsequence diffs) -> epoch_h + scaling.
  # CSPD.model scales activity by 1/duration before detect_msp; the large gaps
  # at removed sleep periods don't change the mode (still the device epoch).
  secs     <- as.numeric(as.POSIXct(wake_datetime))
  dd       <- diff(secs)
  duration <- if (length(dd) > 0L)
    as.numeric(names(sort(table(dd), decreasing = TRUE))[1]) else 60
  epoch_h  <- epoch_h %||% (3600 / duration)

  scaled_activity <- wake_activity / duration

  # Nap-mode MSP (detect_msp else-branch) on the scaled wake activity.
  nap_msp <- .crespo_nap_msp(
    activity             = scaled_activity,
    epoch_h              = epoch_h,
    median_filter_h      = p$nap_median_filter_h,
    pad_h                = p$nap_pad_h,
    nap_median_thr       = p$nap_median_activity_threshold,
    nap_zero_prop_thr    = p$nap_zero_proportion_threshold,
    nap_zero_prop_window = p$nap_zero_proportion_filter_window_size,
    use_and              = p$compute_output_naps_with_logical_and
  )

  # Full CSPD refiner with the nap param set (detect_naps post-processing).
  refined <- .cspd_refine_periods(
    activity        = wake_activity,
    datetime_stamps = wake_datetime,
    msp_detection   = nap_msp,
    condition       = 0L,
    params          = p
  )$refined_output                      # 1 = wake, 0 = sleep (nap)

  # nap_wrapper: state[nap_bool] = 1 - refined_output  -> naps become state 1.
  new_state          <- state_in
  new_state[is_wake] <- ifelse(refined == 0, 1L, 0L)
  x$state            <- new_state
  x$sleep            <- ifelse(new_state == 4L, 0L, as.integer(new_state == 1L))

  x
}

#' Nap-mode MSP — direct translation of detect_msp's `detect_naps = TRUE` branch.
#'
#' Computes the adaptive median filter on the (max-padded, un-marked) activity
#' exactly as the main MSP, then applies the nap thresholds: a high
#' zero-proportion (`zero_prop_filter`) combined (AND/OR) with a low
#' adaptive-median activity. Returns the 1 = wake / 0 = sleep nap detection
#' (Python `improved_sleep_detection`). No zero-mitigation / invalid-zero
#' marking (main-sleep only) and no last morphological filter (off by default).
#' @noRd
.crespo_nap_msp <- function(activity, epoch_h, median_filter_h, pad_h,
                            nap_median_thr, nap_zero_prop_thr,
                            nap_zero_prop_window, use_and = TRUE) {
  data_length      <- length(activity)
  mf_window_size   <- as.integer(epoch_h * median_filter_h) + 1L
  mf_hws           <- as.integer((mf_window_size - 1L) / 2L)
  pad_size         <- as.integer(epoch_h * pad_h)
  maximum_activity <- max(activity, na.rm = TRUE)

  padded_activity <- c(rep(maximum_activity, pad_size), activity,
                       rep(maximum_activity, pad_size))

  # Adaptive median filter (variable window), identical to the main MSP path.
  adaptive_median_filtered_activity <- numeric(data_length)
  half_window_size <- pad_size
  for (i in seq_len(data_length)) {
    i0     <- i - 1L
    center <- i0 + pad_size
    lo_r   <- max(1L, center - half_window_size + 1L)
    hi_r   <- min(length(padded_activity), center + half_window_size + 1L)
    val    <- stats::median(padded_activity[lo_r:hi_r], na.rm = TRUE)
    if (is.nan(val) || is.na(val)) {
      val <- if (i0 > 0L) adaptive_median_filtered_activity[i - 1L] else 0
    }
    adaptive_median_filtered_activity[i] <- val
    if (i0 < (data_length - mf_hws + pad_size - 1L)) {
      if (half_window_size < mf_hws) half_window_size <- half_window_size + 1L
    } else {
      if (half_window_size > pad_size) half_window_size <- half_window_size - 1L
    }
  }

  # Nap thresholds (detect_msp else-branch). zero_prop_filter pads with 1s.
  azp   <- zero_prop_filter(activity, nap_zero_prop_window, pad_value = 1)
  t_azp <- azp > nap_zero_prop_thr
  t_ma  <- adaptive_median_filtered_activity < nap_median_thr

  # Python: improved = where(and/or(t_azp, t_ma), 0, 1)  (0 = sleep, 1 = wake)
  combined <- if (isTRUE(use_and)) (t_azp & t_ma) else (t_azp | t_ma)
  as.integer(!combined)
}

# ── Internal Crespo helpers ───────────────────────────────────────────────────

#' Core Crespo MSP algorithm — direct translation of Python's detect_msp.
#'
#' Returns `list(detection, condition)`: `detection` is the 1 = wake / 0 = sleep
#' MSP vector; `condition` is the (possibly bumped) condition flag. Python's
#' detect_msp sets `condition = 2` exactly when the final activity-median
#' threshold clamps to `min_short_window_thr`; the CSPD refiner uses it.
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
    min_short_window_thr,
    condition            = 0L
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
  # The clamp (threshold < min_short_window_thr) is exactly where Python's
  # detect_msp sets self.condition = 2 (used downstream by the CSPD refiner).
  thr1 <- as.double(stats::quantile(adaptive_median_filtered_activity,
                                    sleep_quantile, names = FALSE))
  condition_out <- condition
  if (thr1 < min_short_window_thr) {
    thr1          <- min_short_window_thr
    condition_out <- 2L
  }

  # Python: improved_sleep_detection = where(filtered > threshold, 1, 0)
  # 1 = wake, 0 = sleep
  list(
    detection = as.integer(adaptive_median_filtered_activity > thr1),
    condition = condition_out
  )
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
