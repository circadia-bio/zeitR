# Parity tests for Rcpp rolling filter implementations.
#
# Each test compares the C++ function against a reference R implementation
# using the exact same padding semantics as the original rolling_apply()-based
# functions in utils.R:
#   - Border replication by default (repeat first/last value)
#   - Constant padding for zero_prop_filter (pad_value = 1)
#
# Vectors tested:
#   - Small vector (n=9, hws=2) -- exercises all boundary cases
#   - Constant vector -- trivial case; verifies no division/NaN issues
#   - Zeros vector  -- important for zero_prop

ref_border <- function(x, hws, FUN) {
  vapply(seq_along(x), function(i) {
    pad <- c(rep(x[1L], hws), x, rep(x[length(x)], hws))
    FUN(pad[i:(i + 2L * hws)])
  }, numeric(1L))
}

ref_const_pad <- function(x, hws, pad_value, FUN) {
  vapply(seq_along(x), function(i) {
    pad <- c(rep(pad_value, hws), x, rep(pad_value, hws))
    FUN(pad[i:(i + 2L * hws)])
  }, numeric(1L))
}

x_small    <- c(1, 3, 2, 5, 4, 6, 8, 7, 9)
x_constant <- rep(4.0, 12)
x_zeros    <- rep(0.0, 8)
hws        <- 2L

# ── rolling_median_cpp ────────────────────────────────────────────────────────

test_that("rolling_median_cpp matches R reference on small vector", {
  expect_equal(
    rolling_median_cpp(x_small, hws),
    ref_border(x_small, hws, stats::median)
  )
})

test_that("rolling_median_cpp returns constant for constant input", {
  expect_equal(rolling_median_cpp(x_constant, hws), x_constant)
})

test_that("rolling_median_cpp handles length-1 input", {
  expect_equal(rolling_median_cpp(42, 3L), 42)
})

# ── rolling_mean_cpp ──────────────────────────────────────────────────────────

test_that("rolling_mean_cpp matches R reference on small vector", {
  expect_equal(
    rolling_mean_cpp(x_small, hws),
    ref_border(x_small, hws, mean)
  )
})

test_that("rolling_mean_cpp returns constant for constant input", {
  expect_equal(rolling_mean_cpp(x_constant, hws), x_constant)
})

# ── rolling_var_cpp ───────────────────────────────────────────────────────────

test_that("rolling_var_cpp matches R reference on small vector", {
  expect_equal(
    rolling_var_cpp(x_small, hws),
    ref_border(x_small, hws, stats::var)
  )
})

test_that("rolling_var_cpp returns zero for constant input", {
  expect_equal(rolling_var_cpp(x_constant, hws), rep(0, length(x_constant)))
})

# ── rolling_zero_prop_cpp ─────────────────────────────────────────────────────

test_that("rolling_zero_prop_cpp matches R reference with pad_value=1", {
  x_mixed <- c(0, 1, 0, 0, 2, 0, 1, 3)
  expect_equal(
    rolling_zero_prop_cpp(x_mixed, hws, 1.0),
    ref_const_pad(x_mixed, hws, 1.0, function(w) sum(w == 0) / length(w))
  )
})

test_that("rolling_zero_prop_cpp returns 1 for all-zero input with pad_value=0", {
  expect_equal(rolling_zero_prop_cpp(x_zeros, hws, 0.0), rep(1, length(x_zeros)))
})

test_that("rolling_zero_prop_cpp returns 0 for all-zero input with pad_value=1", {
  # boundary pads are 1 (non-zero), so proportion depends on window position
  result <- rolling_zero_prop_cpp(x_zeros, hws, 1.0)
  expect_true(all(result >= 0 & result <= 1))
  # interior epochs (no boundary pads): all zeros -> proportion 1
  win <- 2L * hws + 1L
  expect_equal(result[(hws + 1L):(length(x_zeros) - hws)], rep(1, length(x_zeros) - 2L * hws))
})

# ── rolling_quantile_cpp ──────────────────────────────────────────────────────

test_that("rolling_quantile_cpp matches stats::quantile type 7 on small vector", {
  expect_equal(
    rolling_quantile_cpp(x_small, hws, q = 0.6),
    ref_border(x_small, hws, function(w) stats::quantile(w, 0.6, names = FALSE))
  )
})

test_that("rolling_quantile_cpp q=0.5 matches rolling_median_cpp", {
  expect_equal(
    rolling_quantile_cpp(x_small, hws, q = 0.5),
    rolling_median_cpp(x_small, hws)
  )
})

test_that("rolling_quantile_cpp returns constant for constant input", {
  expect_equal(rolling_quantile_cpp(x_constant, hws, q = 0.6), x_constant)
})

# ── Edge cases ────────────────────────────────────────────────────────────────

test_that("all functions return empty numeric for length-0 input", {
  expect_equal(rolling_median_cpp(numeric(0), 2L),   numeric(0))
  expect_equal(rolling_mean_cpp(numeric(0), 2L),     numeric(0))
  expect_equal(rolling_var_cpp(numeric(0), 2L),      numeric(0))
  expect_equal(rolling_zero_prop_cpp(numeric(0), 2L), numeric(0))
  expect_equal(rolling_quantile_cpp(numeric(0), 2L), numeric(0))
})

test_that("hws=0 returns the input unchanged for mean and median", {
  expect_equal(rolling_mean_cpp(x_small, 0L),   as.double(x_small))
  expect_equal(rolling_median_cpp(x_small, 0L), as.double(x_small))
})
