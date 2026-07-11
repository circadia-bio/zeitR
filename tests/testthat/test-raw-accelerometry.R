# Tests for raw_accelerometry.R: compute_activity_counts() and its internal
# helpers (.zero_crossing_indicator(), .epoch_rowsum()).
#
# compute_activity_counts() itself requires mrpheus (a cross-package Suggests
# dependency, not on CRAN) for the actual filtering step, so those tests are
# skipped via skip_if_not_installed("mrpheus") when it isn't available -- the
# internal-helper tests below don't need it and always run.

# ---- .zero_crossing_indicator() ---------------------------------------------

test_that(".zero_crossing_indicator() counts every sign change in a clean alternating signal", {
  out <- .zero_crossing_indicator(c(-1, 1, -1, 1), threshold = 0.1)
  expect_equal(out, c(TRUE, TRUE, TRUE))
})

test_that(".zero_crossing_indicator() absorbs small dead-band excursions without counting them", {
  # 0.02 stays within +/- 0.05 -- the signal never actually clears the
  # dead-band on the positive side, so this should NOT register as a
  # crossing away from and back to the negative zone.
  v   <- c(-1, -1, 0.02, -1, -1)
  out <- .zero_crossing_indicator(v, threshold = 0.05)
  expect_true(all(!out))
})

test_that(".zero_crossing_indicator() returns all-FALSE when fewer than 2 samples clear the dead-band", {
  v   <- c(0, 0.01, -0.01, 0)
  out <- .zero_crossing_indicator(v, threshold = 0.05)
  expect_equal(out, rep(FALSE, length(v) - 1))
})

test_that(".zero_crossing_indicator() treats crossings before the first defined zone as unknown", {
  # First two samples are inside the dead-band (undefined zone); the
  # crossing indicator has no reference to compare them against, so those
  # positions must not be counted, even though the vector overall does
  # contain real crossings later on.
  v   <- c(0.01, 0.01, -1, 1, -1)
  out <- .zero_crossing_indicator(v, threshold = 0.05)
  expect_length(out, 4)
  expect_equal(out, c(FALSE, FALSE, TRUE, TRUE))
})

test_that(".zero_crossing_indicator() always returns length(v) - 1, and length 0 for degenerate input", {
  expect_length(.zero_crossing_indicator(c(1, 2, 3, 4, 5), threshold = 0.1), 4)
  expect_length(.zero_crossing_indicator(numeric(0), threshold = 0.1), 0)
  expect_length(.zero_crossing_indicator(c(1), threshold = 0.1), 0)
})

# ---- .epoch_rowsum() ---------------------------------------------------------

test_that(".epoch_rowsum() sums values within each group", {
  values   <- c(1, 2, 3, 4)
  epoch_id <- c(1L, 1L, 2L, 2L)
  expect_equal(.epoch_rowsum(values, epoch_id, n_epochs = 2L), c(3, 7))
})

test_that(".epoch_rowsum() zero-fills an epoch absent from epoch_id", {
  values   <- c(1, 2, 3, 4)
  epoch_id <- c(1L, 1L, 3L, 3L)  # epoch 2 never appears
  expect_equal(.epoch_rowsum(values, epoch_id, n_epochs = 3L), c(3, 0, 7))
})

# ---- compute_activity_counts() -----------------------------------------------

test_that("compute_activity_counts() validates x/y/z length and sampling_rate/epoch_sec", {
  skip_if_not_installed("mrpheus")

  expect_error(
    compute_activity_counts(1:10, 1:9, 1:10, sampling_rate = 25),
    "same length"
  )
  expect_error(
    compute_activity_counts(rep(0, 100), rep(0, 100), rep(0, 100), sampling_rate = -1),
    "positive number"
  )
  expect_error(
    compute_activity_counts(rep(0, 100), rep(0, 100), rep(0, 100), sampling_rate = 25, epoch_sec = 0),
    "positive number"
  )
})

test_that("compute_activity_counts() errors when epoch_sec x sampling_rate rounds below 2 samples", {
  skip_if_not_installed("mrpheus")

  expect_error(
    compute_activity_counts(rep(0, 100), rep(0, 100), rep(0, 100),
                             sampling_rate = 1, epoch_sec = 1),
    "at least 2"
  )
})

test_that("compute_activity_counts() errors when there aren't enough samples for one epoch", {
  skip_if_not_installed("mrpheus")

  expect_error(
    compute_activity_counts(rep(0, 10), rep(0, 10), rep(0, 10),
                             sampling_rate = 25, epoch_sec = 60),
    "Not enough samples"
  )
})

test_that("compute_activity_counts() warns and truncates trailing samples that don't fill a full epoch", {
  skip_if_not_installed("mrpheus")

  sr        <- 25
  epoch_sec <- 60
  spe       <- sr * epoch_sec       # 1500 samples/epoch
  n         <- spe * 2L + 37L       # 2 full epochs + 37 leftover samples

  set.seed(2)
  x <- rnorm(n, sd = 0.1)
  y <- rnorm(n, sd = 0.1)
  z <- rnorm(n, sd = 0.1) + 1

  expect_warning(
    out <- compute_activity_counts(x, y, z, sampling_rate = sr, epoch_sec = epoch_sec),
    "Dropping the last 37"
  )
  expect_equal(nrow(out), 2L)
})

test_that("compute_activity_counts() returns all-zero metrics for a flat (zero-motion) signal", {
  skip_if_not_installed("mrpheus")

  n    <- 25 * 60 * 3  # 3 minutes at 25 Hz
  zero <- rep(0, n)

  out <- compute_activity_counts(zero, zero, zero, sampling_rate = 25, epoch_sec = 60)

  expect_equal(nrow(out), 3L)
  expect_true(all(out$PIM < 1e-6))
  expect_true(all(out$TAT == 0))
  expect_true(all(out$ZCM == 0))
})

test_that("compute_activity_counts() ZCM approximates 2 x frequency x epoch_sec for a clean sinusoid", {
  skip_if_not_installed("mrpheus")

  sr        <- 25
  f         <- 1    # Hz, well inside the default 0.5-2.7 Hz ActTrust passband
  epoch_sec <- 60
  n_epochs  <- 5
  n         <- sr * epoch_sec * n_epochs
  t         <- seq(0, (n - 1) / sr, by = 1 / sr)

  x <- sin(2 * pi * f * t)  # all activity on the x-axis
  y <- rep(0, n)
  z <- rep(0, n)

  out <- compute_activity_counts(x, y, z, sampling_rate = sr, epoch_sec = epoch_sec)

  expected_zcm <- 2 * f * epoch_sec  # 120 crossings/epoch, two per cycle

  # Interior epochs only: filtfilt's edge padding affects the start/end of
  # the whole recording, not the epoch boundaries in the middle, but this
  # leaves a margin against that regardless.
  interior <- out$ZCM[2:(n_epochs - 1)]
  expect_true(all(abs(interior - expected_zcm) < 15))
})

test_that("compute_activity_counts() PIM and TAT increase with signal amplitude", {
  skip_if_not_installed("mrpheus")

  sr        <- 25
  epoch_sec <- 60
  spe       <- sr * epoch_sec
  n         <- spe * 4L
  t         <- seq(0, (n - 1) / sr, by = 1 / sr)

  set.seed(42)
  quiet  <- rnorm(n, sd = 0.001)
  active <- sin(2 * pi * 1 * t) * 0.5

  # First 2 epochs quiet, last 2 epochs clearly active.
  x <- c(quiet[seq_len(spe * 2L)], active[(spe * 2L + 1L):n])
  y <- rep(0, n)
  z <- rep(0, n)

  out <- compute_activity_counts(x, y, z, sampling_rate = sr, epoch_sec = epoch_sec)

  expect_gt(mean(out$PIM[3:4]), mean(out$PIM[1:2]))
  expect_gt(mean(out$TAT[3:4]), mean(out$TAT[1:2]))
})

test_that("metrics argument restricts which columns are returned", {
  skip_if_not_installed("mrpheus")

  n <- 25 * 60 * 2
  set.seed(1)
  x <- rnorm(n, sd = 0.1)
  y <- rnorm(n, sd = 0.1)
  z <- rnorm(n, sd = 0.1) + 1

  out_all <- compute_activity_counts(x, y, z, sampling_rate = 25, metrics = c("PIM", "TAT", "ZCM"))
  expect_equal(names(out_all), c("epoch", "PIM", "TAT", "ZCM"))

  out_pim <- compute_activity_counts(x, y, z, sampling_rate = 25, metrics = "PIM")
  expect_equal(names(out_pim), c("epoch", "PIM"))

  out_zcm_tat <- compute_activity_counts(x, y, z, sampling_rate = 25, metrics = c("ZCM", "TAT"))
  expect_equal(names(out_zcm_tat), c("epoch", "TAT", "ZCM"))
})

test_that("compute_activity_counts() rejects an unknown metrics value", {
  skip_if_not_installed("mrpheus")

  n <- 25 * 60 * 2
  x <- rep(0, n); y <- rep(0, n); z <- rep(0, n)

  expect_error(
    compute_activity_counts(x, y, z, sampling_rate = 25, metrics = "BOGUS"),
    "should be one of"
  )
})
