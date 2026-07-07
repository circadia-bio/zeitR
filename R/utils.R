# ── Internal utilities ────────────────────────────────────────────────────────
# Shared helper functions used across zeitR.
# Ported from condor_pipeline/algorithms/vendor/condor/functions.py
# (Julius A. P. P. de Paula, Condor Instruments).
# Not exported.

# ── Messages ──────────────────────────────────────────────────────────────────

#' @noRd
zeitr_abort <- function(msg, ...) cli::cli_abort(msg, ...)

#' @noRd
zeitr_warn <- function(msg, ...) cli::cli_warn(msg, ...)

#' @noRd
zeitr_inform <- function(msg, ...) cli::cli_inform(msg, ...)

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

#' Rolling quantile filter with border replication padding
#' @param x numeric vector
#' @param hws integer half-window size
#' @param q quantile probability (default 0.6)
#' @noRd
quantile_filter <- function(x, hws, q = 0.6) {
  rolling_quantile_cpp(as.double(x), as.integer(hws), as.double(q))
}

#' Rolling maximum filter
#' @param x numeric vector
#' @param hws integer half-window size
#' @param pad_value boundary pad value (default 0, for morphological ops)
#' @noRd
max_filter <- function(x, hws, pad_value = 0) {
  rolling_max_cpp(as.double(x), as.integer(hws), replicate = FALSE,
                  pad_value = as.double(pad_value))
}

#' Rolling minimum filter
#' @param x numeric vector
#' @param hws integer half-window size
#' @param pad_value boundary pad value (default 0, for morphological ops)
#' @noRd
min_filter <- function(x, hws, pad_value = 0) {
  rolling_min_cpp(as.double(x), as.integer(hws), replicate = FALSE,
                  pad_value = as.double(pad_value))
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

# ── Zero-sequence detection ───────────────────────────────────────────────────

#' Find contiguous runs of zeros in a binary vector
#'
#' Returns a two-column integer matrix with columns `start` and `end`
#' (1-indexed, inclusive) for each run of zeros with length >=
#' `minimum_length`.
#'
#' @param x numeric vector (treated as binary: 0 vs non-zero)
#' @param minimum_length minimum run length to return (default 1)
#' @noRd
zero_sequences <- function(x, minimum_length = 1L) {
  n <- length(x)
  is_zero <- as.integer(x == 0)

  rle_result <- rle(is_zero)
  ends   <- cumsum(rle_result$lengths)
  starts <- ends - rle_result$lengths + 1L

  zero_runs <- which(rle_result$values == 1L &
                       rle_result$lengths >= minimum_length)

  if (length(zero_runs) == 0L) {
    return(matrix(integer(0), ncol = 2L,
                  dimnames = list(NULL, c("start", "end"))))
  }

  cbind(start = starts[zero_runs], end = ends[zero_runs])
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
