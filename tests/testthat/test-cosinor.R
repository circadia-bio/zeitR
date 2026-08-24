# Tests for cosinor.R: compute_cosinor().
#
# Ported from pyActigraphy's actual Cosinor class as used in the reference
# notebook (Cell 16's _fit_cosinor(), period always locked). The closed-form
# OLS solution used here was verified by direct execution to match
# pyActigraphy's real lmfit (Levenberg-Marquardt) fit exactly when the
# period is locked -- see this session's transcript for the verification
# script if it needs to be re-run.

test_that(".cosinor_fit() matches pyActigraphy's real lmfit-based fit exactly on a deterministic signal", {
  # 3 days, 1-min epochs, a clean cosine peaking at 15:00 (900 min), no
  # noise -- exact recovery expected. Verified against a direct run of
  # pyActigraphy's actual Cosinor.fit() (period locked) + the notebook's
  # own sign-flip/HH:MM conversion on this exact series:
  #   acrophase_time='15:00', acrophase_time_neg=-9.0, MESOR=30.0,
  #   amplitude=15.0, period_min=1440.
  n     <- 3L * 1440L
  t     <- 0:(n - 1L)
  omega <- 2 * pi / 1440
  phi_true <- -omega * 900
  y <- 30 + 15 * cos(omega * t + phi_true)

  fit <- .cosinor_fit(y, period_min = 1440)

  expect_equal(fit$acrophase_time, "15:00")
  expect_equal(fit$acrophase_time_neg, -9.0, tolerance = 1e-6)
  expect_equal(fit$MESOR, 30.0, tolerance = 1e-6)
  expect_equal(fit$amplitude, 15.0, tolerance = 1e-6)
  expect_equal(fit$period_min, 1440)
})

test_that(".cosinor_fit() recovers a peak before noon correctly (acrophase_time_neg stays positive)", {
  # Peak at 06:00 (360 min) -- acrophase_time_neg should stay positive
  # (< 12:00, no -24h shift applied).
  n     <- 3L * 1440L
  t     <- 0:(n - 1L)
  omega <- 2 * pi / 1440
  phi_true <- -omega * 360
  y <- 40 + 10 * cos(omega * t + phi_true)

  fit <- .cosinor_fit(y, period_min = 1440)

  expect_equal(fit$acrophase_time, "06:00")
  expect_equal(fit$acrophase_time_neg, 6.0, tolerance = 1e-6)
})

test_that("compute_cosinor() is wired correctly end to end on a zeitr-shaped fixture", {
  n  <- 3L * 1440L
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  dt <- t0 + 60L * seq.int(0L, n - 1L)
  t     <- 0:(n - 1L)
  omega <- 2 * pi / 1440
  phi_true <- -omega * 900
  activity <- 30 + 15 * cos(omega * t + phi_true)

  d      <- tibble::tibble(datetime = dt, activity = activity)
  result <- compute_cosinor(d)

  expect_equal(result$acrophase_time, "15:00")
  expect_equal(result$acrophase_time_neg, -9.0, tolerance = 1e-6)
  expect_equal(result$MESOR, 30.0, tolerance = 1e-6)
  expect_equal(result$amplitude, 15.0, tolerance = 1e-6)
  expect_equal(result$period_min, 1440)
  expect_equal(result$n_epochs, n)
})

test_that("compute_cosinor() works on an arbitrary column via `col`", {
  n  <- 3L * 1440L
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  dt <- t0 + 60L * seq.int(0L, n - 1L)
  t     <- 0:(n - 1L)
  omega <- 2 * pi / 1440
  phi_true <- -omega * 900
  light <- 30 + 15 * cos(omega * t + phi_true)

  d      <- tibble::tibble(datetime = dt, light = light)
  result <- compute_cosinor(d, col = "light")

  expect_equal(result$acrophase_time, "15:00")
})

test_that("compute_cosinor() accepts a zeitr_result and uses its subject_id", {
  n  <- 3L * 1440L
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  dt <- t0 + 60L * seq.int(0L, n - 1L)
  activity <- rep(50, n)

  d <- tibble::tibble(datetime = dt, activity = activity)
  result_obj <- structure(list(data = d, subject_id = "P042"), class = "zeitr_result")

  result <- compute_cosinor(result_obj)
  expect_equal(result$participant_id, "P042")
})

test_that("compute_cosinor() warns and returns NA for fewer than 3 valid epochs", {
  d <- tibble::tibble(
    datetime = as.POSIXct(c("2024-01-01 00:00:00", "2024-01-01 00:01:00"), tz = "UTC"),
    activity = c(1, 2)
  )
  expect_warning(result <- compute_cosinor(d), "Fewer than 3")
  expect_true(is.na(result$acrophase_time))
  expect_true(is.na(result$MESOR))
})

test_that("compute_cosinor() errors on missing required columns", {
  bad <- tibble::tibble(datetime = Sys.time())
  expect_error(compute_cosinor(bad), "Missing required column")
})

test_that("compute_cosinor() errors when x is neither a recognised object nor a data frame", {
  expect_error(compute_cosinor(list(a = 1)), "must be a")
  expect_error(compute_cosinor("not a recording"), "must be a")
})
