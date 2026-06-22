# Fixture-based regression tests ported from dev/validate_cspd_stage4.R.
#
# Validates the three stages of the CSPD refiner pipeline against Python
# intermediates exported by dev/export_cspd_intermediates.py:
#
#   cspd_stage1_in.csv    - MSP detection entering the refiner (1=wake, 0=sleep)
#   cspd_stage1_out.csv   - after .peak_valley_length_filter
#   cspd_stage2_out.csv   - after .sleep_gap_separation
#   cspd_refined_output.csv  - full .cspd_refine_periods output (1=wake, 0=sleep)
#   cspd_refined_sleep.csv   - per-night bedtime/getuptime indices
#
# Input is the on-wrist subset of the validation recording (state != 4 epochs
# from python_output.csv), mirroring dev/export_cspd_intermediates.py exactly.

refine_fixture <- function(file) {
  path <- system.file("extdata", file, package = "zeitR")
  if (!nzchar(path)) testthat::skip(paste0(file, " not available in inst/extdata"))
  path
}

# Build the on-wrist input once, shared across all tests in this file.
# Returns list(activity, datetime, msp_detection).
build_refine_input <- function() {
  py          <- read.csv(refine_fixture("python_output.csv"),
                          stringsAsFactors = FALSE)
  py$datetime <- as.POSIXct(py$datetime, tz = "UTC")
  onwrist     <- py$state != 4L
  list(
    activity      = as.numeric(py$activity[onwrist]),
    datetime      = py$datetime[onwrist],
    msp_detection = as.integer(scan(refine_fixture("cspd_stage1_in.csv"),
                                    quiet = TRUE))
  )
}

# ---------------------------------------------------------------------------
# Stage 1: peak-valley length filter
# ---------------------------------------------------------------------------

test_that(".peak_valley_length_filter matches Python stage-1 output", {
  pvlf    <- getFromNamespace(".peak_valley_length_filter", "zeitR")
  inp     <- build_refine_input()
  ref     <- as.integer(scan(refine_fixture("cspd_stage1_out.csv"), quiet = TRUE))
  p       <- getFromNamespace(".cspd_acttrust_params", "zeitR")()

  # epoch duration -> peak_valley_minimum_length (mirrors .cspd_refine_periods)
  dd       <- getFromNamespace(".datetime_diff", "zeitR")(inp$datetime)
  duration <- as.numeric(names(sort(table(dd), decreasing = TRUE))[1])
  pvml     <- as.integer(round(p$peak_valley_minimum_length * 60 / duration))

  out <- pvlf(inp$msp_detection, pvml)

  expect_equal(length(out), length(ref))
  n_mismatch <- sum(out != ref)
  expect_equal(n_mismatch, 0L,
               label          = sprintf("%d stage-1 mismatch(es)", n_mismatch),
               expected.label = "0 mismatches")
})

# ---------------------------------------------------------------------------
# Stage 2: sleep-gap separation
# ---------------------------------------------------------------------------

test_that(".sleep_gap_separation matches Python stage-2 output", {
  pvlf    <- getFromNamespace(".peak_valley_length_filter", "zeitR")
  sgs     <- getFromNamespace(".sleep_gap_separation",     "zeitR")
  dtd     <- getFromNamespace(".datetime_diff",            "zeitR")
  inp     <- build_refine_input()
  ref     <- as.integer(scan(refine_fixture("cspd_stage2_out.csv"), quiet = TRUE))
  p       <- getFromNamespace(".cspd_acttrust_params", "zeitR")()

  dd       <- dtd(inp$datetime)
  duration <- as.numeric(names(sort(table(dd), decreasing = TRUE))[1])
  pvml     <- as.integer(round(p$peak_valley_minimum_length * 60 / duration))
  sleep_gap <- p$sleep_maximum_allowed_datetime_gap * 60

  stage1 <- pvlf(inp$msp_detection, pvml)
  out    <- sgs(stage1, dd, sleep_gap)

  expect_equal(length(out), length(ref))
  n_mismatch <- sum(out != ref)
  expect_equal(n_mismatch, 0L,
               label          = sprintf("%d stage-2 mismatch(es)", n_mismatch),
               expected.label = "0 mismatches")
})

# ---------------------------------------------------------------------------
# Full refiner: .cspd_refine_periods
# ---------------------------------------------------------------------------

test_that(".cspd_refine_periods matches Python refined_output exactly", {
  refine  <- getFromNamespace(".cspd_refine_periods", "zeitR")
  inp     <- build_refine_input()
  ref     <- as.integer(scan(refine_fixture("cspd_refined_output.csv"), quiet = TRUE))

  res <- refine(
    activity        = inp$activity,
    datetime_stamps = inp$datetime,
    msp_detection   = inp$msp_detection,
    condition       = 0L
  )

  expect_equal(length(res$refined_output), length(ref))
  n_mismatch <- sum(as.integer(res$refined_output) != ref)
  expect_equal(n_mismatch, 0L,
               label          = sprintf("%d refined_output mismatch(es)", n_mismatch),
               expected.label = "0 mismatches")
})

test_that(".cspd_refine_periods bedtime/getuptime indices match Python refined_sleep", {
  refine   <- getFromNamespace(".cspd_refine_periods", "zeitR")
  inp      <- build_refine_input()
  rs_path  <- system.file("extdata", "cspd_refined_sleep.csv", package = "zeitR")
  if (!nzchar(rs_path)) testthat::skip("cspd_refined_sleep.csv not available")
  rs_ref   <- read.csv(rs_path)

  res <- refine(
    activity        = inp$activity,
    datetime_stamps = inp$datetime,
    msp_detection   = inp$msp_detection,
    condition       = 0L
  )

  expect_equal(nrow(res$refined_sleep_df), nrow(rs_ref))
  expect_equal(res$refined_sleep_df$bedtime_index,  as.integer(rs_ref$bedtime_index))
  expect_equal(res$refined_sleep_df$getuptime_index, as.integer(rs_ref$getuptime_index))
})
