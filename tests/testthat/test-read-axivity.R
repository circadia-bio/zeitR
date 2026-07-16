# Tests for read_axivity.R: read_axivity() (bridges axR::axivity_read_cwa()'s
# raw per-sample output into zeitR's standard epoch-level tibble shape).
#
# Requires both axR (the raw .cwa reader) and mrpheus (compute_activity_counts()'s
# filtering backend) -- skipped via skip_if_not_installed() when either is
# unavailable. axivity_read_cwa() itself is mocked throughout: these tests
# check read_axivity()'s own bridging logic (epoch grouping, sample-rate
# detection, tz handling, metadata assembly), not axR's binary file parsing.

# ---- Shared fixture ----------------------------------------------------------
# 2 epochs x 60s at 25 Hz = 1500 samples/epoch, 3000 samples total (exact,
# no trailing truncation). light/temperature_c are constant *within* each
# epoch but differ *between* epochs, so per-epoch averaging is exactly
# checkable.

make_cwa_fixture <- function(sample_rate = 25, n_epochs = 2L, epoch_sec = 60) {
  spe <- sample_rate * epoch_sec
  n   <- spe * n_epochs

  t0  <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  dt  <- t0 + (seq_len(n) - 1L) / sample_rate

  set.seed(1)
  x <- rnorm(n, sd = 0.05)
  y <- rnorm(n, sd = 0.05)
  z <- rnorm(n, sd = 0.05) + 1   # gravity on z

  epoch_id <- rep(seq_len(n_epochs), each = spe)
  light    <- epoch_id * 100          # epoch 1 -> 100, epoch 2 -> 200, ...
  temp     <- 25 + epoch_id * 1.5     # epoch 1 -> 26.5, epoch 2 -> 28, ...

  out <- data.frame(
    timestamp     = dt,
    x = x, y = y, z = z,
    light         = light,
    temperature_c = temp,
    battery_pct   = 90,
    sample_rate   = sample_rate,
    stringsAsFactors = FALSE
  )
  attr(out, "device_id")  <- "AX3-12345"
  attr(out, "session_id") <- "SESSION1"
  attr(out, "metadata")   <- "Study X protocol"

  if (requireNamespace("tibble", quietly = TRUE)) out <- tibble::as_tibble(out)
  out
}

# Creates a dummy (empty-content) file so read_axivity()'s file.exists()
# check passes; axivity_read_cwa() itself is mocked, so the file's actual
# content is never read.
make_dummy_cwa_path <- function() {
  tf <- tempfile(fileext = ".cwa")
  file.create(tf)
  withr::defer(unlink(tf), envir = parent.frame())
  tf
}

test_that("read_axivity() produces the expected column shape and types", {
  skip_if_not_installed("axR")
  skip_if_not_installed("mrpheus")
  skip_if_not_installed("withr")

  fixture <- make_cwa_fixture()
  testthat::local_mocked_bindings(
    axivity_read_cwa = function(path) fixture,
    .package = "axR"
  )

  path   <- make_dummy_cwa_path()
  result <- read_axivity(path)

  expect_s3_class(result, "zeitr_axivity")
  expect_equal(
    names(result),
    c("datetime", "activity", "int_temp", "ext_temp", "ZCMn", "light",
      "state", "offwrist", "sleep", "TAT")
  )
  expect_equal(nrow(result), 2L)
  expect_true(all(is.na(result$ext_temp)))
  expect_true(all(result$state == 0))
  expect_true(all(result$offwrist == 0))
  expect_true(all(result$sleep == 0))
})

test_that("read_axivity() averages light/int_temp per epoch correctly", {
  skip_if_not_installed("axR")
  skip_if_not_installed("mrpheus")
  skip_if_not_installed("withr")

  fixture <- make_cwa_fixture()
  testthat::local_mocked_bindings(
    axivity_read_cwa = function(path) fixture,
    .package = "axR"
  )

  path   <- make_dummy_cwa_path()
  result <- read_axivity(path)

  expect_equal(result$light,    c(100, 200))
  expect_equal(result$int_temp, c(26.5, 28))
})

test_that("read_axivity() sets datetime to epoch start times, epoch_sec apart", {
  skip_if_not_installed("axR")
  skip_if_not_installed("mrpheus")
  skip_if_not_installed("withr")

  fixture <- make_cwa_fixture(epoch_sec = 60)
  testthat::local_mocked_bindings(
    axivity_read_cwa = function(path) fixture,
    .package = "axR"
  )

  path   <- make_dummy_cwa_path()
  result <- read_axivity(path, epoch_sec = 60)

  expect_equal(as.numeric(diff(result$datetime), units = "secs"), 60)
  expect_equal(result$datetime[1], as.POSIXct("2024-01-01 00:00:00", tz = "UTC"))
})

test_that("read_axivity() re-labels the time zone without shifting the clock reading", {
  skip_if_not_installed("axR")
  skip_if_not_installed("mrpheus")
  skip_if_not_installed("withr")

  fixture <- make_cwa_fixture()
  testthat::local_mocked_bindings(
    axivity_read_cwa = function(path) fixture,
    .package = "axR"
  )

  path   <- make_dummy_cwa_path()
  result <- read_axivity(path, tz = "America/Sao_Paulo")

  expect_equal(attr(result$datetime, "tzone"), "America/Sao_Paulo")
  # Same wall-clock reading, just re-labelled -- hour/minute unchanged.
  expect_equal(format(result$datetime[1], "%H:%M:%S", tz = "America/Sao_Paulo"), "00:00:00")
})

test_that("read_axivity() attaches metadata from the cwa attributes and call arguments", {
  skip_if_not_installed("axR")
  skip_if_not_installed("mrpheus")
  skip_if_not_installed("withr")

  fixture <- make_cwa_fixture(sample_rate = 25)
  testthat::local_mocked_bindings(
    axivity_read_cwa = function(path) fixture,
    .package = "axR"
  )

  path   <- make_dummy_cwa_path()
  result <- read_axivity(path, filter_low = 0.25, filter_high = 2.5)
  meta   <- attr(result, "metadata")

  expect_equal(meta$device_id, "AX3-12345")
  expect_equal(meta$session_id, "SESSION1")
  expect_equal(meta$sample_rate, 25)
  expect_equal(meta$epoch_sec, 60)
  expect_equal(meta$filter_low, 0.25)
  expect_equal(meta$filter_high, 2.5)
  expect_equal(meta$cwa_metadata, "Study X protocol")
  expect_equal(meta$source_file, normalizePath(path, mustWork = FALSE))
})

test_that("read_axivity() uses the dominant sample_rate and warns about outliers", {
  skip_if_not_installed("axR")
  skip_if_not_installed("mrpheus")
  skip_if_not_installed("withr")

  fixture <- make_cwa_fixture(sample_rate = 25)
  # Corrupt a handful of rows with a different reported rate.
  fixture$sample_rate[1:5] <- 50

  testthat::local_mocked_bindings(
    axivity_read_cwa = function(path) fixture,
    .package = "axR"
  )

  path <- make_dummy_cwa_path()
  expect_warning(result <- read_axivity(path), "sample_rate other than")
  expect_equal(attr(result, "metadata")$sample_rate, 25)
})

test_that("read_axivity() activity/ZCMn/TAT match a direct compute_activity_counts() call", {
  skip_if_not_installed("axR")
  skip_if_not_installed("mrpheus")
  skip_if_not_installed("withr")

  fixture <- make_cwa_fixture(sample_rate = 25)
  testthat::local_mocked_bindings(
    axivity_read_cwa = function(path) fixture,
    .package = "axR"
  )

  path      <- make_dummy_cwa_path()
  result    <- read_axivity(path)
  reference <- compute_activity_counts(
    fixture$x, fixture$y, fixture$z,
    sampling_rate = 25, epoch_sec = 60,
    filter_low = 0.25, filter_high = 2.5,
    zcm_threshold = 0.01, tat_threshold = 0.05
  )

  expect_equal(result$activity, reference$PIM)
  expect_equal(result$ZCMn,     reference$ZCM)
  expect_equal(result$TAT,      reference$TAT)
})

test_that("read_axivity() errors when the file doesn't exist", {
  skip_if_not_installed("axR")

  expect_error(read_axivity("/nonexistent/path/does-not-exist.cwa"), "File not found")
})

test_that(".check_axR_pkg() errors with an installation hint when axR isn't installed", {
  skip_if(requireNamespace("axR", quietly = TRUE),
          "axR is installed locally; can't exercise the not-installed path here")
  expect_error(.check_axR_pkg(), "circadia-bio.r-universe.dev")
})
