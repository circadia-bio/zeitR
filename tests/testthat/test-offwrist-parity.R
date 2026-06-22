# Regression tests for the bimodal off-wrist detector against the Condor
# Python reference, using the fixtures bundled in inst/extdata:
#   input1.txt            - validation ActTrust recording
#   python_output.csv     - Python per-epoch state (state == 4 is off-wrist)
#   python_vp_periods.csv - Python valley-peak offwrist periods [start, end)
#
# These lock in epoch-level parity with the reference implementation. The
# valley-peak and Ashman's-D checks read intermediate results captured into the
# package-internal `.zeitr_debug` environment, which is populated only while the
# `zeitR.debug_offwrist` option is TRUE.

# Resolve a bundled fixture, skipping the test if it is unavailable.
offwrist_fixture <- function(file) {
  path <- system.file("extdata", file, package = "zeitR")
  if (!nzchar(path)) testthat::skip(paste0(file, " not available in inst/extdata"))
  path
}

# Run the full detector with intermediate capture enabled, returning the
# captured `.zeitr_debug` environment. The option only needs to be set while
# the pipeline runs; the captured environment persists afterwards.
run_with_debug <- function() {
  old <- options(zeitR.debug_offwrist = TRUE)
  on.exit(options(old), add = TRUE)
  rec <- read_acttrust(offwrist_fixture("input1.txt"), tz = "UTC")
  detect_offwrist_bimodal(prepare_actigraphy(rec))
  getFromNamespace(".zeitr_debug", "zeitR")
}

test_that("bimodal off-wrist detection reproduces the Python reference exactly", {
  skip_if_not_installed("mclust")

  rec  <- read_acttrust(offwrist_fixture("input1.txt"), tz = "UTC")
  prep <- detect_offwrist_bimodal(prepare_actigraphy(rec))

  r_ow  <- as.integer(prep$state == 4)
  py_ow <- as.integer(read.csv(offwrist_fixture("python_output.csv"))$state == 4)

  expect_equal(length(r_ow), length(py_ow))
  # Exact epoch-level agreement: 0 false positives and 0 false negatives.
  expect_identical(r_ow, py_ow)
})

test_that("valley-peak offwrist periods match the Python reference set", {
  skip_if_not_installed("mclust")

  dbg <- run_with_debug()
  vp  <- dbg$vp_df
  expect_false(is.null(vp))

  got <- data.frame(start = as.integer(vp$start), end = as.integer(vp$end))
  got <- got[order(got$start), ]
  row.names(got) <- NULL

  ref <- read.csv(offwrist_fixture("python_vp_periods.csv"))
  ref <- data.frame(start = as.integer(ref$start), end = as.integer(ref$end))
  ref <- ref[order(ref$start), ]
  row.names(ref) <- NULL

  expect_equal(got, ref)
})

test_that("Ashman's D leaves is_highly_separable FALSE so the report filter runs", {
  skip_if_not_installed("mclust")

  dbg <- run_with_debug()

  # mclust's full-convergence fit regressed D to 3.27 (> 3), which disabled the
  # description-report filter; the sklearn-style EM keeps D just under 3.
  expect_false(is.null(dbg$ashman))
  expect_lt(dbg$ashman, 3)
  expect_false(isTRUE(dbg$is_highly_separable))

  # With the filter active it must remove at least one period (the spurious
  # valley-peak blocks in the 60k/65k regions).
  expect_true(nrow(dbg$post_desc) < nrow(dbg$pre_desc))
})
