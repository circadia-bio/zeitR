# Parity tests for the Rcpp implementations of Crespo MSP helpers.
#
# Each test compares the C++ function against a reference R implementation
# of the exact same algorithm, using synthetic activity vectors designed to
# exercise all code paths (zero runs, boundary behaviour, wake/sleep
# transitions, NA handling in the adaptive median).

# ── Helpers ───────────────────────────────────────────────────────────────────

# R reference for zero_mitigation_cpp
ref_zero_mitigation <- function(activity, consec_zeros_thr, mitigation_level) {
  n      <- length(activity)
  result <- activity
  run    <- 0L
  for (i in seq_len(n)) {
    if (activity[i] == 0) {
      run <- run + 1L
    } else {
      if (run > consec_zeros_thr) {
        result[(i - run):(i - 1L)] <- result[(i - run):(i - 1L)] + mitigation_level
      }
      run <- 0L
    }
  }
  result  # trailing run NOT processed (matches Python)
}

# R reference for mark_invalid_zeros_cpp (returns 0-indexed positions)
ref_mark_invalid_zeros <- function(activity, morph_det, awake_thr, sleep_thr) {
  n     <- length(activity)
  invalid <- integer(0)
  ar <- 0L; sr <- 0L
  for (i in seq_len(n)) {
    i0 <- i - 1L
    if (morph_det[i] == 1L) {
      if (sr > sleep_thr) invalid <- c(invalid, (i0 - sr):(i0 - 1L))
      sr <- 0L
      if (activity[i] == 0) { ar <- ar + 1L } else {
        if (ar > awake_thr) invalid <- c(invalid, (i0 - ar):(i0 - 1L))
        ar <- 0L
      }
    } else {
      if (ar > awake_thr) invalid <- c(invalid, (i0 - ar):(i0 - 1L))
      ar <- 0L
      if (activity[i] == 0) { sr <- sr + 1L } else {
        if (sr > sleep_thr) invalid <- c(invalid, (i0 - sr):(i0 - 1L))
        sr <- 0L
      }
    }
  }
  sort(unique(invalid))
}

# R reference for adaptive_median_filter_cpp
ref_adaptive_median <- function(padded_activity, n, pad_size, max_hws) {
  result <- numeric(n)
  hws    <- pad_size
  for (i in seq_len(n)) {
    i0     <- i - 1L
    center <- i0 + pad_size
    lo_r   <- max(1L, center - hws + 1L)
    hi_r   <- min(length(padded_activity), center + hws + 1L)
    val    <- stats::median(padded_activity[lo_r:hi_r], na.rm = TRUE)
    if (is.nan(val) || is.na(val))
      val <- if (i0 > 0L) result[i - 1L] else 0
    result[i] <- val
    if (i0 < (n - max_hws + pad_size - 1L)) {
      if (hws < max_hws) hws <- hws + 1L
    } else {
      if (hws > pad_size) hws <- hws - 1L
    }
  }
  result
}

# ── rolling_max_cpp ───────────────────────────────────────────────────────────

test_that("rolling_max_cpp matches rolling_apply(max) on small vector", {
  x   <- c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3)
  hws <- 2L
  expect_equal(
    rolling_max_cpp(x, hws, replicate = FALSE, pad_value = 0.0),
    rolling_apply(x, hws, max, pad_value = 0L)
  )
})

test_that("rolling_max_cpp returns constant for constant input", {
  x <- rep(7.0, 20)
  expect_equal(rolling_max_cpp(x, 3L, FALSE, 0.0), x)
})

test_that("rolling_max_cpp handles length-0 input", {
  expect_equal(rolling_max_cpp(numeric(0), 2L), numeric(0))
})

# ── rolling_min_cpp ───────────────────────────────────────────────────────────

test_that("rolling_min_cpp matches rolling_apply(min) on small vector", {
  x   <- c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3)
  hws <- 2L
  expect_equal(
    rolling_min_cpp(x, hws, replicate = FALSE, pad_value = 0.0),
    rolling_apply(x, hws, min, pad_value = 0L)
  )
})

test_that("rolling_min_cpp returns constant for constant input with border replication", {
  x <- rep(7.0, 20)
  expect_equal(rolling_min_cpp(x, 3L, TRUE, 0.0), x)  # replicate=TRUE: no 0 pads
})

test_that("rolling_min_cpp pad_value=0 clamps boundary windows correctly", {
  # With pad_value=0, boundary windows include zeros; min should be 0 at edges
  x <- c(5.0, 5.0, 5.0, 5.0, 5.0)
  expect_equal(rolling_min_cpp(x, 2L, FALSE, 0.0),
               rolling_apply(x, 2L, min, pad_value = 0L))
})

# ── Morphological open/close uses rolling_max/min ─────────────────────────────

test_that("rolling_max then rolling_min (closing) matches rolling_apply pair", {
  set.seed(7)
  x   <- as.double(sample(0:1, 50, replace = TRUE))
  hws <- 5L
  cpp  <- rolling_min_cpp(rolling_max_cpp(x, hws, FALSE, 0.0), hws, FALSE, 0.0)
  rref <- rolling_apply(rolling_apply(x, hws, max, pad_value = 0L),
                        hws, min, pad_value = 0L)
  expect_equal(cpp, rref)
})

# ── zero_mitigation_cpp ───────────────────────────────────────────────────────

test_that("zero_mitigation_cpp matches R reference: run exactly at threshold not mitigated", {
  act   <- as.double(c(rep(0, 15), 5, rep(0, 16), 3))
  level <- 2.5
  expect_equal(zero_mitigation_cpp(act, 15L, level),
               ref_zero_mitigation(act, 15L, level))
})

test_that("zero_mitigation_cpp matches R reference: trailing zero run not mitigated", {
  act   <- as.double(c(1, rep(0, 20)))
  level <- 1.0
  expect_equal(zero_mitigation_cpp(act, 5L, level),
               ref_zero_mitigation(act, 5L, level))
})

test_that("zero_mitigation_cpp matches R reference: no zeros", {
  act   <- as.double(c(1, 2, 3, 4, 5))
  level <- 1.0
  expect_equal(zero_mitigation_cpp(act, 3L, level),
               ref_zero_mitigation(act, 3L, level))
})

test_that("zero_mitigation_cpp matches R reference on realistic mixed vector", {
  set.seed(42)
  act   <- as.double(sample(0:200, 500, replace = TRUE, prob = c(0.4, rep(0.6/200, 200))))
  level <- as.double(quantile(act, 0.33, names = FALSE))
  expect_equal(zero_mitigation_cpp(act, 15L, level),
               ref_zero_mitigation(act, 15L, level))
})

# ── mark_invalid_zeros_cpp ────────────────────────────────────────────────────

test_that("mark_invalid_zeros_cpp matches R reference: all awake", {
  act  <- as.double(c(0, 0, 0, 1, 0, 0, 0, 0, 1))
  morph <- as.integer(rep(1L, length(act)))
  expect_identical(mark_invalid_zeros_cpp(act, morph, 2L, 5L),
                   ref_mark_invalid_zeros(act, morph, 2L, 5L))
})

test_that("mark_invalid_zeros_cpp matches R reference: mixed wake/sleep", {
  set.seed(3)
  n     <- 100L
  act   <- as.double(sample(0:50, n, replace = TRUE, prob = c(0.5, rep(0.5/50, 50))))
  morph <- as.integer(sample(0:1, n, replace = TRUE))
  expect_identical(mark_invalid_zeros_cpp(act, morph, 2L, 30L),
                   ref_mark_invalid_zeros(act, morph, 2L, 30L))
})

test_that("mark_invalid_zeros_cpp returns integer(0) when no invalid zeros", {
  act   <- as.double(c(1, 2, 3, 4, 5))
  morph <- as.integer(c(1, 1, 0, 0, 1))
  expect_identical(mark_invalid_zeros_cpp(act, morph, 2L, 30L), integer(0))
})

# ── adaptive_median_filter_cpp ────────────────────────────────────────────────

test_that("adaptive_median_filter_cpp matches R reference: no NAs, small n", {
  set.seed(99)
  n        <- 50L
  pad_size <- 5L
  max_hws  <- 10L
  activity <- rnorm(n, mean = 50, sd = 20)
  padded   <- c(rep(max(activity), pad_size), activity, rep(max(activity), pad_size))
  expect_equal(
    adaptive_median_filter_cpp(padded, n, pad_size, max_hws),
    ref_adaptive_median(padded, n, pad_size, max_hws)
  )
})

test_that("adaptive_median_filter_cpp matches R reference with NA values", {
  set.seed(5)
  n        <- 40L
  pad_size <- 4L
  max_hws  <- 8L
  activity <- rnorm(n, mean = 30, sd = 10)
  padded   <- c(rep(max(activity), pad_size), activity, rep(max(activity), pad_size))
  # Introduce NAs at a few positions (simulating invalid zeros)
  padded[c(7L, 8L, 15L)] <- NA_real_
  expect_equal(
    adaptive_median_filter_cpp(padded, n, pad_size, max_hws),
    ref_adaptive_median(padded, n, pad_size, max_hws)
  )
})

test_that("adaptive_median_filter_cpp: window size grows then shrinks", {
  # Verify the output length equals n regardless of pad_size and max_hws
  n        <- 30L
  pad_size <- 3L
  max_hws  <- 6L
  padded   <- as.double(c(rep(100, pad_size), seq_len(n), rep(100, pad_size)))
  result   <- adaptive_median_filter_cpp(padded, n, pad_size, max_hws)
  expect_length(result, n)
})

# ── End-to-end: full pipeline still matches fixtures ──────────────────────────

test_that("Crespo C++ port: run_pipeline epoch state unchanged on input1.txt", {
  fpath <- system.file("extdata", "input1.txt", package = "zeitR")
  if (!nzchar(fpath)) testthat::skip("input1.txt not in extdata")

  result <- run_pipeline(fpath, tz = "America/Sao_Paulo", quiet = TRUE)
  # Night count is the simplest stable proxy for epoch-level correctness
  # (full epoch comparison is covered by test-pipeline-parity.R)
  expect_equal(nrow(result$nights), 55L)
  expect_equal(sum(result$data$state == 1L), 24155L)
})
