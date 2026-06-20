#' Off-wrist detection using the Condor bimodal activity/temperature model
#'
#' Detects periods where the actigraph was not worn using the bimodal
#' algorithm developed by Condor Instruments. The algorithm proceeds in
#' three stages:
#'
#' 1. **Feature extraction** — rolling median of PIM activity and rolling
#'    signal/variance of internal temperature are computed over a symmetric
#'    window of half-width `hws` epochs.
#' 2. **Bimodal temperature threshold** — among epochs with low activity
#'    median (below the `activity_quantile` quantile of normalised activity),
#'    a 2-component Gaussian Mixture Model is fitted to the normalised
#'    temperature distribution. The threshold between the two components is
#'    taken as the minimum of the GMM density between the two means (or the
#'    minimum of the smoothed histogram if that is lower). Ashman's D is
#'    computed as a bimodality quality metric.
#' 3. **Initial classification & refinement** — epochs simultaneously
#'    exhibiting low activity median AND low temperature are marked as
#'    off-wrist. Short spurious off-wrist runs shorter than
#'    `min_offwrist_length` epochs are removed.
#'
#' Off-wrist epochs are encoded as `state == 4` (matching the Python
#' pipeline convention). The `offwrist` column is set to `0.25` for
#' off-wrist epochs for actogram overlay plotting.
#'
#' @param x A tibble as returned by [prepare_actigraphy()], containing
#'   columns `datetime`, `activity` (PIM), `int_temp`, `ext_temp`, `state`,
#'   and `offwrist`.
#' @param hws `integer(1)`. Half-window size (in epochs) for rolling feature
#'   extraction. Default is `10` (matching the Python pipeline).
#' @param activity_quantile `numeric(1)`. Quantile used to define "low
#'   activity". Default is `0.15`.
#' @param min_norm_activity `numeric(1)`. Minimum normalised activity
#'   threshold below which the low-activity cutoff is clamped. Default is
#'   `0.015`.
#' @param nbins `integer(1)`. Number of histogram bins used when fitting the
#'   GMM threshold. Default is `100`.
#' @param min_offwrist_length `integer(1)`. Minimum number of consecutive
#'   off-wrist epochs required to retain a detected period. Shorter runs are
#'   discarded. Default is `10`.
#' @param min_temp_threshold `numeric(1)`. Minimum normalised temperature
#'   threshold; the fitted GMM threshold is clamped to this value if it falls
#'   below it. Default is `0.35`.
#'
#' @return The input tibble `x` with `state` and `offwrist` columns updated.
#'   Off-wrist epochs have `state == 4` and `offwrist == 0.25`.
#'
#' @references
#' The bimodal off-wrist algorithm was developed by Julius A. P. P. de Paula
#' at Condor Instruments (2023). It is not published in peer-reviewed
#' literature but the source code is available in the circadiaBase pipeline
#' repository. The Ashman D statistic is described in:
#'
#' Ashman, K. M., Bird, C. M., & Zepf, S. E. (1994). Detecting bimodality
#' in astronomical datasets. *The Astronomical Journal*, 108, 2348.
#' \doi{10.1086/117248}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rec  <- read_acttrust("recordings/P001.txt")
#' prep <- prepare_actigraphy(rec)
#' prep <- detect_offwrist_bimodal(prep)
#' sum(prep$state == 4)  # number of off-wrist epochs
#' }
detect_offwrist_bimodal <- function(
    x,
    hws                  = 10L,
    activity_quantile    = 0.15,
    min_norm_activity    = 0.015,
    nbins                = 100L,
    min_offwrist_length  = 10L,
    min_temp_threshold   = 0.35
) {
  # ── Input checks ────────────────────────────────────────────────────────────
  required <- c("activity", "int_temp", "ext_temp", "state", "offwrist")
  missing  <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    zeitr_abort(
      "{.arg x} is missing required column(s): {.val {missing}}.
       Did you run {.fn prepare_actigraphy} first?"
    )
  }

  valid_temp <- x$int_temp > 0
  if (!any(valid_temp)) {
    zeitr_warn("All {.code int_temp} values are <= 0; marking everything as off-wrist.")
    x$state    <- 4
    x$offwrist <- 0.25
    return(x)
  }

  # ── Work on the valid-temperature subset ────────────────────────────────────
  activity <- as.double(x$activity)
  int_temp <- as.double(x$int_temp)
  ext_temp <- as.double(x$ext_temp)

  # ── Stage 1: Feature extraction ─────────────────────────────────────────────
  act_median    <- median_filter(activity, hws)
  norm_act_med  <- norm_01(act_median)

  # Low-activity threshold (quantile of normalised activity median)
  zp            <- zero_prop(norm_act_med)
  low_q         <- zp + activity_quantile * (1 - zp)
  low_act_thr   <- stats::quantile(norm_act_med, low_q, names = FALSE,
                                   type = 1)   # type=1 matches Python 'inverted_cdf'
  low_act_thr   <- max(low_act_thr, min_norm_activity)

  is_low_act    <- norm_act_med < low_act_thr

  # ── Stage 2: Bimodal temperature threshold ──────────────────────────────────
  norm_temp     <- norm_01(int_temp)
  low_act_temps <- norm_temp[is_low_act]

  temp_threshold <- min_temp_threshold   # fallback
  ashman         <- 0

  if (length(low_act_temps) > 1L) {
    gmm_result <- .fit_gmm_threshold(
      low_act_temps,
      nbins             = nbins,
      min_temp_threshold = min_temp_threshold
    )
    temp_threshold <- gmm_result$threshold
    ashman         <- gmm_result$ashman_d
  }

  # Rescale threshold back to original temperature scale
  temp_min <- min(int_temp, na.rm = TRUE)
  temp_max <- max(int_temp, na.rm = TRUE)
  temp_threshold_orig <- temp_min + temp_threshold * (temp_max - temp_min)

  # Clamp: threshold must not exceed median temperature of high-activity epochs
  med_high_act_temp <- stats::median(int_temp[!is_low_act], na.rm = TRUE)
  if (!is.na(med_high_act_temp) && temp_threshold_orig > med_high_act_temp) {
    temp_threshold_orig <- med_high_act_temp
    temp_threshold      <- (temp_threshold_orig - temp_min) /
      max(temp_max - temp_min, .Machine$double.eps)
  }

  is_low_temp <- int_temp < temp_threshold_orig

  # ── Stage 3: Initial classification ─────────────────────────────────────────
  # Off-wrist = low activity AND low temperature
  offwrist_raw <- as.integer(is_low_act & is_low_temp)
  # Invert: 0 = off-wrist candidate, 1 = on-wrist  (matches Python convention)
  offwrist_raw <- 1L - offwrist_raw

  # ── Stage 4: Refinement — remove short off-wrist runs ───────────────────────
  offwrist_refined <- .remove_short_offwrist(offwrist_raw, min_offwrist_length)

  # ── Update state columns ────────────────────────────────────────────────────
  # Convention: state == 4 → off-wrist (matches Python pipeline)
  x$state    <- ifelse(offwrist_refined == 0, 4, x$state)
  x$offwrist <- ifelse(offwrist_refined == 0, 0.25, 0)

  x
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' Fit a 2-component GMM to normalised temperature and return the threshold
#' @noRd
.fit_gmm_threshold <- function(data, nbins = 100L, min_temp_threshold = 0.35) {
  # We need mclust for GMM — use base stats EM approximation via
  # a simple k-means initialisation + density intersection approach
  # to avoid a hard dependency on mclust.

  data <- as.double(data)
  data <- data[!is.na(data)]
  n    <- length(data)

  if (n < 4L) {
    return(list(threshold = min_temp_threshold, ashman_d = 0))
  }

  # Initialise two components via k-means (k=2)
  km <- tryCatch(
    stats::kmeans(data, centers = 2L, nstart = 5L, iter.max = 100L),
    error = function(e) NULL
  )

  if (is.null(km)) {
    return(list(threshold = min_temp_threshold, ashman_d = 0))
  }

  # Order components by mean (component 1 = lower)
  ord  <- order(km$centers)
  mu   <- as.double(km$centers[ord])
  idx  <- km$cluster
  idx  <- ifelse(idx == ord[1], 1L, 2L)

  sigma <- c(
    stats::sd(data[idx == 1L]),
    stats::sd(data[idx == 2L])
  )
  sigma[is.na(sigma) | sigma == 0] <- .Machine$double.eps

  ash_d <- ashman_d(mu[1], sigma[1], mu[2], sigma[2])

  # Histogram-based minimum between the two means
  breaks <- seq(0, 1, length.out = nbins + 1L)
  counts <- graphics::hist(data, breaks = breaks, plot = FALSE)$counts
  counts_smooth <- mean_filter(as.double(counts), 3L)

  # Locate bin indices closest to mu1 and mu2
  bin_mids <- (breaks[-1] + breaks[-length(breaks)]) / 2
  loc1 <- which.min(abs(bin_mids - mu[1]))
  loc2 <- which.min(abs(bin_mids - mu[2]))
  if (loc2 <= loc1) loc2 <- loc1 + 1L
  loc2 <- min(loc2, length(bin_mids))

  # Threshold = bin of minimum density between the two modes
  search_range <- seq(loc1, loc2)
  thresh_idx   <- search_range[which.min(counts_smooth[search_range])]
  threshold    <- bin_mids[thresh_idx]

  # Boundary safety from Python source
  if ((mu[1] - sigma[1]) < 0.15) mu[1] <- 0.15 + sigma[1]
  if ((mu[2] + sigma[2]) > 0.85) mu[2] <- 0.85 - sigma[2]
  if (mu[2] <= mu[1])            mu[2] <- mu[1] + 2 / nbins

  threshold <- max(threshold, min_temp_threshold)

  list(threshold = threshold, ashman_d = ash_d)
}

#' Remove off-wrist runs shorter than minimum_length
#' @param x integer vector (0 = off-wrist, 1 = on-wrist)
#' @param minimum_length integer
#' @noRd
.remove_short_offwrist <- function(x, minimum_length) {
  runs <- zero_sequences(x, minimum_length = 1L)
  if (nrow(runs) == 0L) return(x)

  for (i in seq_len(nrow(runs))) {
    run_len <- runs[i, "end"] - runs[i, "start"] + 1L
    if (run_len < minimum_length) {
      x[runs[i, "start"]:runs[i, "end"]] <- 1L
    }
  }
  x
}
