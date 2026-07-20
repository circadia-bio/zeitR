# Tests for circ_utils.R: circ_mean_h() and circ_sd_h() (circular statistics
# for clock-time variables, Fix 20 / Fix 24 in the JRSV pipeline).

# ---- Helper -------------------------------------------------------------
# 0 and 24 represent the same clock time (midnight), but are different as
# raw doubles. At the exact antipodal point (e.g. two points symmetric
# around midnight), mean(sin(theta)) involves cancelling two
# nearly-exactly-opposite values; atan2() can land on either side of the
# 0/24 boundary depending on which way the cancellation rounds, which is a
# genuine floating-point property of atan2 at a discontinuity, not a package
# bug. expect_equal_hour() treats 0 and 24 (and anything in between within
# tolerance) as equal.
expect_equal_hour <- function(actual, expected, tolerance = 1e-6) {
  diff <- abs(actual - expected) %% 24
  diff <- min(diff, 24 - diff)
  expect_lt(diff, tolerance)
}

test_that("circ_mean_h() handles midnight wrap correctly", {
  expect_equal_hour(circ_mean_h(c(23.5, 0.5)), 0)
  expect_equal_hour(circ_mean_h(c(23.0, 1.0)), 0)
  expect_equal_hour(circ_mean_h(c(22.0, 2.0)), 0)
})

test_that("circ_mean_h() matches plain mean when times don't straddle midnight", {
  expect_equal(circ_mean_h(c(7.0, 9.0)), 8)
  expect_equal(circ_mean_h(c(10.0, 12.0, 14.0)), 12)
})

test_that("circ_mean_h() drops NAs silently", {
  expect_equal(circ_mean_h(c(7.0, NA, 9.0)), 8)
  expect_equal(circ_mean_h(c(NA, NA)), NA_real_)
})

test_that("circ_mean_h() returns NA_real_ for an empty vector after NA removal", {
  expect_equal(circ_mean_h(numeric(0)), NA_real_)
})

test_that("circ_mean_h() is invariant to rotation of the whole input", {
  x <- c(1.0, 2.0, 3.0)
  shifted <- (x + 12) %% 24
  # Rotating every input by 12h should rotate the mean by 12h too (mod 24)
  expect_equal_hour((circ_mean_h(x) + 12) %% 24, circ_mean_h(shifted))
})

test_that("circ_sd_h() returns NA_real_ for fewer than 2 values", {
  expect_equal(circ_sd_h(numeric(0)), NA_real_)
  expect_equal(circ_sd_h(c(5.0)), NA_real_)
  expect_equal(circ_sd_h(c(NA, 5.0)), NA_real_)
})

test_that("circ_sd_h() is small for tightly clustered times, even across midnight", {
  # For two points 1h apart, the circular SD is close to (but slightly above,
  # by construction of the formula) half the separation, i.e. ~0.5h -- not
  # "near zero" in an absolute sense, but small relative to the 0-6h range
  # the formula can produce for widely dispersed times.
  expect_lt(circ_sd_h(c(23.5, 0.5)), 0.6)
  expect_lt(circ_sd_h(c(23.9, 0.0, 0.1)), 0.6)
})

test_that("circ_sd_h() is large for times 12h apart (maximally dispersed)", {
  # Two points exactly opposite on the clock -> resultant length ~ 0 -> SD is large
  sd_opposite <- circ_sd_h(c(6.0, 18.0))
  sd_close    <- circ_sd_h(c(6.0, 7.0))
  expect_gt(sd_opposite, sd_close)
})

test_that("circ_sd_h() drops NAs silently", {
  expect_equal(circ_sd_h(c(6.0, 7.0, NA)), circ_sd_h(c(6.0, 7.0)))
})

test_that("circ_sd_h() is non-negative", {
  expect_gte(circ_sd_h(c(1.0, 5.0, 23.0, 12.0)), 0)
})
