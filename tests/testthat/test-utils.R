# Tests for utils.R internal helpers not already exercised (or not fully
# branch-covered) by other test files. median_filter()/mean_filter()/
# var_filter()/diff5()/zero_prop_filter() are thin single-expression Rcpp
# wrappers already called throughout the offwrist/sleep-period pipeline
# (covered indirectly via test-offwrist-parity.R etc.); this file focuses on
# the branches that aren't guaranteed to be hit by ordinary pipeline data.

test_that("norm_01() scales to [0, 1] on a normal range", {
  x <- c(0, 5, 10)
  expect_equal(norm_01(x), c(0, 0.5, 1))
})

test_that("norm_01() handles a constant vector without dividing by zero", {
  # span == 0 branch: every value is identical, so the result should be a
  # vector of zeros (x - rng[1]), not NaN from a 0/0 division.
  x <- rep(7, 5)
  expect_equal(norm_01(x), rep(0, 5))
  expect_true(all(is.finite(norm_01(x))))
})

test_that("norm_01() drops NAs from the range calculation", {
  x <- c(0, 5, 10, NA)
  expect_equal(norm_01(x), c(0, 0.5, 1, NA))
})

test_that("zero_prop() computes the proportion of exact zeros", {
  expect_equal(zero_prop(c(0, 0, 1, 2)), 0.5)
  expect_equal(zero_prop(c(1, 2, 3)), 0)
  expect_equal(zero_prop(c(0, 0, 0)), 1)
})

test_that("zero_prop() returns 0 for an empty vector rather than NaN", {
  expect_equal(zero_prop(numeric(0)), 0)
})

test_that("zero_prop() counts NAs toward the denominator but never the numerator", {
  # zero_prop() is sum(x == 0, na.rm = TRUE) / length(x) -- na.rm applies only
  # to the numerator's sum(). NAs are NOT dropped from the denominator, so an
  # NA counts as "not exactly zero" for the proportion, rather than being
  # excluded from consideration entirely.
  expect_equal(zero_prop(c(0, 0, NA)), 2 / 3)
  expect_equal(zero_prop(c(NA, NA)), 0)
})

test_that("median_filter(), mean_filter(), var_filter() run and return the right length", {
  x <- c(1, 3, 2, 5, 4, 6, 8, 7, 9)
  expect_equal(length(median_filter(x, 2L)), length(x))
  expect_equal(length(mean_filter(x, 2L)), length(x))
  expect_equal(length(var_filter(x, 2L)), length(x))
  expect_equal(mean_filter(rep(4, 10), 2L), rep(4, 10))
})

test_that("zero_prop_filter() runs and returns the right length", {
  x <- c(0, 0, 1, 0, 2, 0, 0, 3)
  out <- zero_prop_filter(x, 2L)
  expect_equal(length(out), length(x))
  expect_true(all(out >= 0 & out <= 1))
})

test_that("diff5() returns a derivative estimate the same length as the input", {
  x <- c(1, 2, 4, 7, 11, 16, 22)
  out <- diff5(x)
  expect_equal(length(out), length(x))
  expect_true(all(is.finite(out)))

  # A linear ramp has a constant derivative everywhere the five-point
  # stencil has enough support.
  ramp <- as.double(1:20)
  expect_equal(diff5(ramp)[3:18], rep(1, 16))
})

test_that("ashman_d() measures separation between two Gaussian components", {
  # Well-separated components -> D clearly above 2 (the paper's bimodality cutoff)
  expect_gt(ashman_d(0, 1, 10, 1), 2)
  # Identical components -> D = 0
  expect_equal(ashman_d(5, 1, 5, 1), 0)
})

test_that("ashman_d() returns 0 rather than dividing by zero when both sigmas are 0", {
  expect_equal(ashman_d(0, 0, 10, 0), 0)
})

test_that("label_states() maps the four known codes and NA-fills anything else", {
  x <- label_states(c(0L, 1L, 4L, 7L, 99L))
  expect_equal(as.character(x), c("wake", "sleep", "off-wrist", "nap", NA))
  expect_true(is.ordered(x))
  expect_equal(levels(x), c("wake", "sleep", "nap", "off-wrist"))
})

test_that("%||% returns the left side when it is non-NULL, non-NA, non-empty-string", {
  expect_equal("a" %||% "b", "a")
  expect_equal(1L %||% 2L, 1L)
})

test_that("%||% falls back to the right side for NULL, NA, or an empty string", {
  expect_equal(NULL %||% "b", "b")
  expect_equal(NA_character_ %||% "b", "b")
  expect_equal("" %||% "b", "b")
})

test_that("zeitr_abort()/zeitr_warn()/zeitr_inform() signal the expected condition types", {
  expect_error(zeitr_abort("boom"), "boom")
  expect_warning(zeitr_warn("careful"), "careful")
  expect_message(zeitr_inform("fyi"), "fyi")
})
