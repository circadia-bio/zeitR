# End-to-end regression test: full pipeline vs Python reference
#
# Runs run_pipeline() on the bundled ActTrust fixture and asserts bit-exact
# agreement with `python_output.csv` (epoch-level state) and
# `python_nights.csv` (per-night statistics). Any future change to the
# pipeline that shifts even one epoch will cause this test to fail.
#
# Fixture files required in inst/extdata:
#   input1.txt         - validation ActTrust recording
#   python_output.csv  - Python per-epoch state (columns: datetime, state, ...)
#   python_nights.csv  - Python per-night stats (columns: bt, gt, nap, tbt,
#                        waso, sol, soi, tst, nw, eff, ...)
#
# Notes on column names:
#   r_nights uses: night, is_nap, bed_time, get_up_time, tbt, tst, waso,
#                  sol, soi, nw, eff
#   python_nights uses: bts, gts, bt, gt, nap, tbt, waso, sol, soi, tst, nw, eff
#   bt/gt in python_nights are 0-indexed epoch integers; bed_time/get_up_time
#   in r_nights are POSIXct timestamps — validated against datetime column.

pipeline_fixture <- function(file) {
  path <- system.file("extdata", file, package = "zeitR")
  if (!nzchar(path)) testthat::skip(paste0(file, " not available in inst/extdata"))
  path
}

# Run the pipeline with quiet = TRUE to suppress the known timestamp-gap
# warning that input1.txt always produces (5 minor gaps in the raw file).
run_pipeline_quiet <- function() {
  run_pipeline(pipeline_fixture("input1.txt"), tz = "UTC", quiet = TRUE)
}

# ---------------------------------------------------------------------------
# Epoch-level state parity
# ---------------------------------------------------------------------------

test_that("full pipeline reproduces Python epoch-level state exactly", {
  skip_if_not_installed("mclust")

  result  <- run_pipeline_quiet()
  py_out  <- read.csv(pipeline_fixture("python_output.csv"))

  r_state  <- as.integer(result$data$state)
  py_state <- as.integer(py_out$state)

  # Length guard — catches silent truncation or expansion.
  expect_equal(length(r_state), length(py_state),
               label          = "epoch count",
               expected.label = "Python epoch count")

  # Epoch-by-epoch identity.
  n_mismatch <- sum(r_state != py_state)
  expect_equal(
    n_mismatch, 0L,
    label          = sprintf("%d epoch mismatch(es)", n_mismatch),
    expected.label = "0 mismatches"
  )
})

# ---------------------------------------------------------------------------
# Per-layer epoch counts (offwrist / sleep / nap)
# ---------------------------------------------------------------------------

test_that("offwrist epoch count matches Python reference", {
  skip_if_not_installed("mclust")

  result <- run_pipeline_quiet()
  py_out <- read.csv(pipeline_fixture("python_output.csv"))

  expect_equal(
    sum(result$data$state == 4L),
    sum(py_out$state == 4L),
    label          = "R off-wrist epoch count",
    expected.label = "Python off-wrist epoch count"
  )
})

test_that("sleep epoch count matches Python reference", {
  skip_if_not_installed("mclust")

  result <- run_pipeline_quiet()
  py_out <- read.csv(pipeline_fixture("python_output.csv"))

  # States 1 (sleep) and 7 (nap-sleep) are both scored as sleep.
  expect_equal(
    sum(result$data$state %in% c(1L, 7L)),
    sum(py_out$state %in% c(1L, 7L)),
    label          = "R sleep epoch count",
    expected.label = "Python sleep epoch count"
  )
})

test_that("nap epoch count matches Python reference", {
  skip_if_not_installed("mclust")

  result <- run_pipeline_quiet()
  py_out <- read.csv(pipeline_fixture("python_output.csv"))

  expect_equal(
    sum(result$data$state == 7L),
    sum(py_out$state == 7L),
    label          = "R nap epoch count",
    expected.label = "Python nap epoch count"
  )
})

# ---------------------------------------------------------------------------
# Per-night statistics parity
# ---------------------------------------------------------------------------

test_that("nightly statistics match Python reference", {
  skip_if_not_installed("mclust")

  result    <- run_pipeline_quiet()
  py_nights <- read.csv(pipeline_fixture("python_nights.csv"))

  r_nights <- result$nights

  # Night count.
  expect_equal(nrow(r_nights), nrow(py_nights),
               label          = "R night count",
               expected.label = "Python night count")

  # nap flag: r_nights$is_nap vs python_nights$nap.
  expect_equal(r_nights$is_nap, as.logical(py_nights$nap),
               label          = "nap flags",
               expected.label = "Python nap flags")

  # Boundary timestamps: python_nights stores POSIXct in bts/gts columns;
  # r_nights stores them as bed_time/get_up_time. Compare after coercion to UTC.
  py_bed  <- as.POSIXct(py_nights$bts, tz = "UTC")
  py_getup <- as.POSIXct(py_nights$gts, tz = "UTC")

  expect_equal(as.numeric(r_nights$bed_time),  as.numeric(py_bed),
               tolerance      = 1,              # within 1 second
               label          = "bed_time",
               expected.label = "Python bts")

  expect_equal(as.numeric(r_nights$get_up_time), as.numeric(py_getup),
               tolerance      = 1,
               label          = "get_up_time",
               expected.label = "Python gts")

  # Core sleep metrics.
  expect_equal(r_nights$tbt,  py_nights$tbt,  tolerance = 1e-9, label = "TBT")
  expect_equal(r_nights$waso, py_nights$waso, tolerance = 1e-9, label = "WASO")
  expect_equal(r_nights$sol,  py_nights$sol,  tolerance = 1e-9, label = "SOL")
  expect_equal(r_nights$soi,  py_nights$soi,  tolerance = 1e-9, label = "SOI")
  expect_equal(r_nights$tst,  py_nights$tst,  tolerance = 1e-9, label = "TST")
  expect_equal(r_nights$nw,   py_nights$nw,   tolerance = 1e-9, label = "NW")
  expect_equal(r_nights$eff,  py_nights$eff,  tolerance = 1e-9, label = "efficiency")
})
