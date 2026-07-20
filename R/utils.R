# ── Internal utilities ────────────────────────────────────────────────────────
# Shared helper functions used across zeitR.
# Ported from condor_pipeline/algorithms/vendor/condor/functions.py
# (Julius A. P. P. de Paula, Condor Instruments).
# Not exported.

# ── Messages ──────────────────────────────────────────────────────────────────

#' @noRd
zeitr_abort <- function(msg, ..., .envir = parent.frame()) {
  cli::cli_abort(msg, ..., .envir = .envir)
}

#' @noRd
zeitr_warn <- function(msg, ..., .envir = parent.frame()) {
  cli::cli_warn(msg, ..., .envir = .envir)
}

#' @noRd
zeitr_inform <- function(msg, ..., .envir = parent.frame()) {
  cli::cli_inform(msg, ..., .envir = .envir)
}

# ── Scaling ───────────────────────────────────────────────────────────────────

#' Min-max scale a numeric vector to `[0, 1]`
#' @param x numeric vector
#' @noRd
norm_01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  span <- rng[2] - rng[1]
  if (span == 0) return(x - rng[1])
  (x - rng[1]) / span
}

# ── Zero proportion ───────────────────────────────────────────────────────────

#' Proportion of exact zeros in a numeric vector
#' @param x numeric vector
#' @noRd
zero_prop <- function(x) {
  if (length(x) == 0) return(0)
  sum(x == 0, na.rm = TRUE) / length(x)
}

#' Rolling zero-proportion filter
#' @param x numeric vector
#' @param hws integer half-window size
#' @param pad_value value used for boundary padding (default 1)
#' @noRd
zero_prop_filter <- function(x, hws, pad_value = 1) {
  rolling_zero_prop_cpp(as.double(x), as.integer(hws), as.double(pad_value))
}

# ── Rolling filters ───────────────────────────────────────────────────────────

#' Apply a rolling function with boundary padding
#'
#' Used as the R reference implementation in Rcpp parity tests
#' (test-crespo-cpp-parity.R) for rolling_max_cpp()/rolling_min_cpp(); not
#' called anywhere in the package's own pipeline code, which calls the Rcpp
#' rolling filters directly.
#'
#' @param x numeric vector
#' @param hws integer half-window size
#' @param FUN function to apply over each window (receives a numeric vector)
#' @param pad_value scalar pad value; if NULL uses border replication
#' @noRd
rolling_apply <- function(x, hws, FUN, pad_value = NULL) {
  n <- length(x)
  if (is.null(pad_value)) {
    pad_start <- rep(x[1],     hws)
    pad_end   <- rep(x[n],     hws)
  } else {
    pad_start <- rep(pad_value, hws)
    pad_end   <- rep(pad_value, hws)
  }
  padded <- c(pad_start, x, pad_end)
  vapply(seq_len(n), function(i) {
    FUN(padded[i:(i + 2L * hws)])
  }, numeric(1))
}

#' Rolling median filter with border replication padding
#' @param x numeric vector
#' @param hws integer half-window size
#' @noRd
median_filter <- function(x, hws) {
  rolling_median_cpp(as.double(x), as.integer(hws))
}

#' Rolling mean filter with border replication padding
#' @param x numeric vector
#' @param hws integer half-window size
#' @noRd
mean_filter <- function(x, hws) {
  rolling_mean_cpp(as.double(x), as.integer(hws))
}

#' Rolling variance filter with border replication padding
#' @param x numeric vector
#' @param hws integer half-window size
#' @noRd
var_filter <- function(x, hws) {
  rolling_var_cpp(as.double(x), as.integer(hws))
}

# ── Five-point derivative ─────────────────────────────────────────────────────

#' Five-point stencil derivative estimate
#'
#' Approximates the first derivative using the five-point central difference
#' formula, with one-sided approximations at the four boundary points.
#'
#' @param x numeric vector (length >= 5)
#' @param delta epoch duration in seconds (default 1)
#' @noRd
diff5 <- function(x, delta = 1) {
  diff5_cpp(as.double(x), as.double(delta))
}

# ── Ashman's D ────────────────────────────────────────────────────────────────

#' Ashman's D statistic for bimodality
#'
#' Measures separation between two Gaussian components.
#' Values > 2 indicate clearly separated modes.
#'
#' @param mu1,mu2 component means
#' @param sigma1,sigma2 component standard deviations
#' @noRd
ashman_d <- function(mu1, sigma1, mu2, sigma2) {
  denom <- sigma1^2 + sigma2^2
  if (denom <= 0) return(0)
  sqrt(2 / denom) * abs(mu1 - mu2)
}

# ── State labels ─────────────────────────────────────────────────────────────

#' Convert integer epoch states to a labelled factor
#'
#' Converts the integer `state` column produced by the zeitR pipeline into a
#' human-readable factor. Useful for display, plotting, and export — the
#' internal `state` column always stays integer to preserve Python reference
#' parity.
#'
#' | Integer | Label       |
#' |---------|-------------|
#' | `0`     | `"wake"`     |
#' | `1`     | `"sleep"`    |
#' | `4`     | `"off-wrist"`|
#' | `7`     | `"nap"`      |
#'
#' Any value not in the table above is silently converted to `NA`.
#'
#' @param x integer (or numeric) vector of epoch states, as found in
#'   `result$data$state`.
#'
#' @return An ordered factor with levels
#'   `c("wake", "sleep", "nap", "off-wrist")`, the same length as `x`.
#'
#' @export
#'
#' @examples
#' label_states(c(0L, 1L, 0L, 4L, 1L, 7L))
#' # [1] wake  sleep wake  off-wrist sleep nap
#' # Levels: wake < sleep < nap < off-wrist
#'
#' \dontrun{
#' result <- run_pipeline("recordings/P001.txt")
#' result$data$state_label <- label_states(result$data$state)
#' }
label_states <- function(x) {
  factor(
    c("0" = "wake", "1" = "sleep", "4" = "off-wrist", "7" = "nap")[as.character(as.integer(x))],
    levels  = c("wake", "sleep", "nap", "off-wrist"),
    ordered = TRUE
  )
}

# ── NULL coalescing operator ───────────────────────────────────────────────────

#' @noRd
`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nchar(a) > 0) a else b
