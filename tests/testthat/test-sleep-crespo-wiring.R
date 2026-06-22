# Fixture-based regression test ported from dev/validate_sleep_crespo_wiring.R.
#
# Validates the detect_sleep_crespo() end-to-end wiring against the Python
# cspd_wrapper output: the on-wrist state after detect_sleep_crespo(refine=TRUE)
# must equal 1 - cspd_refined_output.csv (1 = sleep period, 0 = wake), and
# off-wrist epochs must be preserved as state 4.
#
# This test sits between test-detect_sleep_crespo.R (algorithm-level unit
# tests) and test-pipeline-parity.R (full chain): it isolates the
# detect_sleep_crespo layer without depending on offwrist detection, nap
# detection, or WASO scoring.

wiring_fixture <- function(file) {
  path <- system.file("extdata", file, package = "zeitR")
  if (!nzchar(path)) testthat::skip(paste0(file, " not available in inst/extdata"))
  path
}

test_that("detect_sleep_crespo(refine=TRUE) on-wrist state matches 1 - cspd_refined_output", {
  py          <- read.csv(wiring_fixture("python_output.csv"),
                          stringsAsFactors = FALSE)
  py$datetime <- as.POSIXct(py$datetime, tz = "UTC")

  # Reconstruct the exact pre-sleep-detection state: off-wrist preserved (4),
  # all on-wrist epochs reset to 0 (wake). This mirrors the Python cspd_wrapper
  # entry point and is the only new step relative to the validated refiner.
  x <- data.frame(
    datetime = py$datetime,
    activity = as.double(py$activity),
    state    = ifelse(py$state == 4L, 4L, 0L)
  )
  onwrist <- x$state != 4L

  out <- detect_sleep_crespo(x, refine = TRUE)

  # On-wrist state: detect_sleep_crespo encodes sleep periods as state 1,
  # which should equal 1 - cspd_refined_output (0 = sleep period in Python).
  refined_output <- as.integer(scan(wiring_fixture("cspd_refined_output.csv"),
                                    quiet = TRUE))
  expected_state <- ifelse(refined_output == 0L, 1L, 0L)

  r_ow_state <- as.integer(out$state[onwrist])

  expect_equal(length(r_ow_state), length(expected_state))
  n_mismatch <- sum(r_ow_state != expected_state)
  expect_equal(n_mismatch, 0L,
               label          = sprintf("%d on-wrist state mismatch(es)", n_mismatch),
               expected.label = "0 mismatches")

  # Off-wrist epochs untouched.
  expect_true(all(out$state[!onwrist] == 4L))

  # sleep column tracks state == 1.
  expect_identical(as.integer(out$sleep), as.integer(out$state == 1L))
})
