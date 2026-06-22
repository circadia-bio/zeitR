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

  # ── Match Python: mark temperature <= 0 as off-wrist immediately ─────────
  # Python: data = df[df["TEMPERATURE"] > 0]
  # Epochs with temp <= 0 are automatically off-wrist; algorithm runs on rest
  x$state[!valid_temp]    <- 4L
  x$offwrist[!valid_temp] <- 0.25

  # ── Work only on valid-temperature subset (matching Python exactly) ────────
  n_full   <- nrow(x)
  idx_valid <- which(valid_temp)

  activity <- as.double(x$activity[valid_temp])
  int_temp <- as.double(x$int_temp[valid_temp])
  ext_temp <- as.double(x$ext_temp[valid_temp])

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
    # Ashman's D drives is_highly_separable (D > 3 disables the description
    # report filter). Python fits the GMM with sklearn GaussianMixture(tol=1e-2);
    # that loose tolerance stops EM early, leaving wider component variances and
    # D just under 3. mclust converges fully, over-separates the components, and
    # inflates D past 3. Recompute D with a faithful sklearn-style EM so
    # is_highly_separable matches Python. The threshold is kept from the
    # (verified) mclust fit.
    em <- .gmm_em_sklearn_1d(low_act_temps)
    ashman <- if (!is.null(em)) ashman_d(em$mu1, em$sigma1, em$mu2, em$sigma2)
              else gmm_result$ashman_d
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

  # ── Stage 4: Full three-stage border refinement ────────────────────────────
  # Compute inputs needed by the refiner
  temp_var    <- var_filter(int_temp, hws)
  norm_tv     <- norm_01(temp_var)
  temp_deriv  <- diff5(int_temp)
  td_var      <- var_filter(temp_deriv, hws)
  dif_temp    <- int_temp - ext_temp

  epoch_h     <- .estimate_epoch_h(x$datetime)

  offwrist_refined <- .bimodal_refine_acttrust(
    initial_offwrist      = offwrist_raw,
    activity              = activity,
    activity_median       = act_median,
    temperature           = int_temp,
    norm_temp_variance    = norm_tv,
    temp_derivative       = temp_deriv,
    temp_derivative_variance = td_var,
    temperature_threshold = temp_threshold_orig,
    ashman                = ashman,
    activity_median_low   = as.integer(is_low_act),
    is_low_temp           = is_low_temp,
    filter_hws            = hws,
    dif_temp              = dif_temp,
    epoch_hour            = epoch_h,
    datetime_stamps       = x$datetime[valid_temp],
    do_near_all_off_detection = TRUE
  )

  # ── Update state columns ────────────────────────────────────────────────────
  # Map refined result back to full-length vector
  # Python: out = np.zeros(len(df)); out[valid_mask] = final_offwrist
  # state == 4 -> off-wrist (0 in Python convention = off-wrist)
  x$state[idx_valid]    <- ifelse(offwrist_refined == 0L, 4L, x$state[idx_valid])
  x$offwrist[idx_valid] <- ifelse(offwrist_refined == 0L, 0.25, 0)

  x
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' Fit a 2-component GMM to normalised temperature and return the threshold
#'
#' Matches Python's bimodal_thresh() using GaussianMixture with k-means++
#' initialisation. Requires the mclust package.
#' @noRd
.fit_gmm_threshold <- function(data, nbins = 100L, min_temp_threshold = 0.35) {
  data <- as.double(data)
  data <- data[!is.na(data)]
  n    <- length(data)

  if (n < 4L) {
    return(list(threshold = min_temp_threshold, ashman_d = 0))
  }

  # ── Histogram (matching Python: bins=nbins, range=(0,1)) ──────────────────
  breaks       <- seq(0, 1, length.out = nbins + 1L)
  bin_mids     <- (breaks[-1] + breaks[-length(breaks)]) / 2
  counts       <- as.double(graphics::hist(data, breaks = breaks,
                                           plot = FALSE)$counts)
  counts_filt  <- mean_filter(counts, 3L)

  # ── GMM via mclust (2 components, matching sklearn GaussianMixture) ───────
  if (!requireNamespace("mclust", quietly = TRUE)) {
    zeitr_warn(
      "{.pkg mclust} is not installed; falling back to k-means GMM approximation.
       Install it with {.code install.packages('mclust')} for better accuracy."
    )
    return(.fit_gmm_threshold_kmeans(data, nbins, min_temp_threshold,
                                     counts, counts_filt, breaks, bin_mids))
  }

  gm <- tryCatch(
    mclust::Mclust(data, G = 2, modelNames = "V", verbose = FALSE),
    error = function(e) NULL
  )

  if (is.null(gm)) {
    return(.fit_gmm_threshold_kmeans(data, nbins, min_temp_threshold,
                                     counts, counts_filt, breaks, bin_mids))
  }

  # Order components by mean (component 1 = lower)
  ord    <- order(gm$parameters$mean)
  mu     <- gm$parameters$mean[ord]
  sigma  <- sqrt(gm$parameters$variance$sigmasq[ord])
  weight <- gm$parameters$pro[ord]

  sigma[is.na(sigma) | sigma == 0] <- .Machine$double.eps

  ash_d <- ashman_d(mu[1], sigma[1], mu[2], sigma[2])

  # ── GMM density curve (matching Python dist_gm) ───────────────────────────
  x       <- seq(0, 1, length.out = nbins)
  max_cnt <- max(counts)
  a1      <- weight[1] * max_cnt
  a2      <- weight[2] * max_cnt

  gmm_density <- a1 * exp(-(x - mu[1])^2 / (2 * sigma[1]^2)) +
                 a2 * exp(-(x - mu[2])^2 / (2 * sigma[2]^2))

  # ── Boundary corrections (matching Python) ────────────────────────────────
  mu1 <- mu[1]; mu2 <- mu[2]
  if ((mu1 - sigma[1]) < 0.15) mu1 <- 0.15 + sigma[1]
  if ((mu2 + sigma[2]) > 0.85) mu2 <- 0.85 - sigma[2]
  if (mu2 <= mu1) mu2 <- mu1 + 2 / nbins

  # ── Search range ──────────────────────────────────────────────────────────
  loc1 <- which.min(abs(bin_mids - mu1))
  loc2 <- which.min(abs(bin_mids - mu2))
  if (loc2 <= loc1) loc2 <- loc1 + 1L
  loc2 <- min(loc2, length(bin_mids))

  search <- seq(loc1, loc2)

  # ── Two threshold candidates (matching Python thresh_id_0 / thresh_id_1) ──
  thresh_id_0 <- search[which.min(gmm_density[search])]
  thresh_id_1 <- search[which.min(counts_filt[search])]

  # Pick the one with lower counts_filt (matching Python logic)
  if (counts_filt[thresh_id_0] < counts_filt[thresh_id_1]) {
    thresh_id <- thresh_id_0
  } else {
    thresh_id <- thresh_id_1
  }

  threshold <- bin_mids[thresh_id]
  threshold <- max(threshold, min_temp_threshold)

  list(threshold = threshold, ashman_d = ash_d)
}

#' K-means fallback for GMM threshold (used when mclust is unavailable)
#' @noRd
.fit_gmm_threshold_kmeans <- function(data, nbins, min_temp_threshold,
                                       counts, counts_filt, breaks, bin_mids) {
  km <- tryCatch(
    stats::kmeans(data, centers = 2L, nstart = 5L, iter.max = 100L),
    error = function(e) NULL
  )

  if (is.null(km)) {
    return(list(threshold = min_temp_threshold, ashman_d = 0))
  }

  ord   <- order(km$centers)
  mu    <- as.double(km$centers[ord])
  idx   <- ifelse(km$cluster == ord[1], 1L, 2L)
  sigma <- c(stats::sd(data[idx == 1L]), stats::sd(data[idx == 2L]))
  sigma[is.na(sigma) | sigma == 0] <- .Machine$double.eps

  ash_d <- ashman_d(mu[1], sigma[1], mu[2], sigma[2])

  loc1 <- which.min(abs(bin_mids - mu[1]))
  loc2 <- which.min(abs(bin_mids - mu[2]))
  if (loc2 <= loc1) loc2 <- loc1 + 1L

  search    <- seq(loc1, loc2)
  thresh_id <- search[which.min(counts_filt[search])]
  threshold <- max(bin_mids[thresh_id], min_temp_threshold)

  list(threshold = threshold, ashman_d = ash_d)
}

#' Faithful 1-D two-component Gaussian-mixture EM, matching sklearn
#' GaussianMixture(n_components = 2, covariance_type = "full",
#' reg_covar = 1e-6, tol = 1e-2, max_iter = 100, n_init = 5).
#'
#' Returns components ordered by mean: list(mu1, sigma1, mu2, sigma2, w1, w2),
#' or NULL if the fit is not possible. Used only for Ashman's D; the loose
#' tol = 1e-2 (matching Python) means EM stops early with wider variances.
#' @noRd
.gmm_em_sklearn_1d <- function(x, n_init = 5L, max_iter = 100L,
                               tol = 1e-2, reg_covar = 1e-6) {
  x <- as.double(x[!is.na(x)])
  n <- length(x)
  if (n < 4L) return(NULL)
  LOG2PI <- log(2 * pi)
  eps10  <- 10 * .Machine$double.eps

  estep <- function(mu, v, w) {
    lp1 <- log(w[1]) - 0.5 * (LOG2PI + log(v[1]) + (x - mu[1])^2 / v[1])
    lp2 <- log(w[2]) - 0.5 * (LOG2PI + log(v[2]) + (x - mu[2])^2 / v[2])
    m   <- pmax(lp1, lp2)
    lse <- m + log(exp(lp1 - m) + exp(lp2 - m))
    list(r1 = exp(lp1 - lse), r2 = exp(lp2 - lse), lb = mean(lse))
  }
  mstep <- function(r1, r2) {
    nk <- c(sum(r1), sum(r2)) + eps10
    mu <- c(sum(r1 * x), sum(r2 * x)) / nk
    v  <- c(sum(r1 * (x - mu[1])^2), sum(r2 * (x - mu[2])^2)) / nk + reg_covar
    list(mu = mu, v = v, w = nk / n)
  }

  # Deterministic seeds approximating k-means++ initial means; the kmeans
  # restarts add spread, and best-lower-bound selection (below) mirrors n_init.
  seeds <- list(
    as.double(stats::quantile(x, c(0.25, 0.75))),
    as.double(stats::quantile(x, c(0.10, 0.90))),
    as.double(stats::quantile(x, c(0.40, 0.60)))
  )
  for (j in seq_len(max(0L, n_init - length(seeds)))) {
    km <- tryCatch(stats::kmeans(x, centers = 2L, nstart = 1L, iter.max = 50L),
                   error = function(e) NULL)
    seeds[[length(seeds) + 1L]] <-
      if (!is.null(km)) as.double(sort(km$centers))
      else as.double(stats::quantile(x, c(0.30, 0.70)))
  }

  best <- NULL; best_lb <- -Inf
  for (s in seeds) {
    mu0 <- sort(s)
    if (!all(is.finite(mu0)) || mu0[1] == mu0[2]) next
    cl  <- ifelse(abs(x - mu0[1]) <= abs(x - mu0[2]), 1L, 2L)
    par <- mstep(as.double(cl == 1L), as.double(cl == 2L))
    prev <- -Inf
    for (it in seq_len(max_iter)) {
      es  <- estep(par$mu, par$v, par$w)        # lower bound at current params
      par <- mstep(es$r1, es$r2)                # sklearn updates every iter
      if (is.finite(es$lb) && abs(es$lb - prev) < tol) break
      prev <- es$lb
    }
    lb <- estep(par$mu, par$v, par$w)$lb
    if (is.finite(lb) && lb > best_lb) { best_lb <- lb; best <- par }
  }
  if (is.null(best)) return(NULL)

  ord <- order(best$mu)
  list(mu1 = best$mu[ord[1]], sigma1 = sqrt(best$v[ord[1]]),
       mu2 = best$mu[ord[2]], sigma2 = sqrt(best$v[ord[2]]),
       w1 = best$w[ord[1]],   w2 = best$w[ord[2]])
}


