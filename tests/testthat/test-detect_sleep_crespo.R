# Regression and unit tests for the Crespo main-sleep-period detector
# (`detect_sleep_crespo()` / internal `.crespo_msp()`), using fixtures in
# inst/extdata:
#   input1.txt                     - validation ActTrust recording
#   python_output.csv              - Python per-epoch output; `activity` is the
#                                    exact input the Crespo reference was run on.
#   python_crespo_intermediate.csv - Python `CrespoAlgorithm().detect_msp` run on
#                                    that activity with DEFAULT parameters
#                                    (quantile 1/3, epoch_hour 60). Columns:
#                                    `act_med` (adaptive median filter) and
#                                    `final` (1 = wake, 0 = sleep).
#
# The exact-parity test locks in the two MSP fixes:
#   * the centred coarse median-filter window (was shifted left by mf_hws - 1),
#   * the scipy-style (border_value = 0) morphological closing/opening.
# Together these took `final` from thousands of disagreements to 0.

crespo_fixture <- function(file) {
  path <- system.file("extdata", file, package = "zeitR")
  if (!nzchar(path)) testthat::skip(paste0(file, " not available in inst/extdata"))
  path
}

# ── Exact parity against the Python detect_msp reference ──────────────────────

test_that(".crespo_msp reproduces the Python detect_msp reference exactly", {
  # Run on the SAME activity vector the reference was generated from, with the
  # same default parameters (quantile 1/3, epoch_h 60). This is hermetic: it
  # does not depend on the reader or off-wrist detection, so a failure points
  # squarely at the MSP algorithm.
  ref      <- read.csv(crespo_fixture("python_crespo_intermediate.csv"))
  activity <- as.double(read.csv(crespo_fixture("python_output.csv"))$activity)
  expect_equal(length(activity), nrow(ref))

  out <- getFromNamespace(".crespo_msp", "zeitR")(
    activity = activity, epoch_h = 60, median_filter_h = 8, pad_h = 1,
    sleep_quantile = 1 / 3, morph_size = 61L,
    consec_zeros_thr = 15L, awake_zeros_thr = 2L, sleep_zeros_thr = 30L,
    zero_mitigation_q = 0.33, min_short_window_thr = 1.0
  )

  # .crespo_msp returns 1 = wake, 0 = sleep, matching the reference `final`.
  expect_identical(out, as.integer(ref$final))
  expect_equal(sum(out == 0L), 25435L)   # sleep epochs
})

test_that("read_acttrust + prepare feed .crespo_msp the reference activity", {
  # Guards the assumption above: the entry-point activity must match the vector
  # the reference was built on, byte-for-byte, or the parity target is moot.
  rec  <- read_acttrust(crespo_fixture("input1.txt"), tz = "UTC")
  prep <- prepare_actigraphy(rec)
  py   <- read.csv(crespo_fixture("python_output.csv"))

  expect_equal(as.double(prep$activity), as.double(py$activity))
})

# ── detect_sleep_crespo wiring ────────────────────────────────────────────────

test_that("detect_sleep_crespo preserves off-wrist epochs and excludes them from sleep", {
  n  <- 300L
  dt <- as.POSIXct("2021-05-27 00:00:00", tz = "UTC") + seq_len(n) * 30 - 30
  x  <- data.frame(
    datetime = dt,
    activity = as.double(rep(c(rep(0, 30), rep(200, 30)), length.out = n)),
    state    = 0L
  )
  ow <- 100:150
  x$state[ow] <- 4L   # mark an off-wrist block

  # refine = FALSE exercises the on-wrist subset + MSP->state wiring without the
  # full CSPD refiner (which needs realistic-length data; validated separately
  # against the inst/extdata fixtures via dev/validate_sleep_crespo_wiring.R).
  out <- detect_sleep_crespo(x, epoch_h = 2, refine = FALSE)

  # Off-wrist epochs stay state 4 and are never scored as sleep.
  expect_true(all(out$state[ow] == 4L))
  expect_true(all(out$sleep[ow] == 0L))
  # Everywhere else state is wake/sleep (0/1) and sleep tracks state == 1.
  expect_true(all(out$state[-ow] %in% c(0L, 1L)))
  expect_identical(as.integer(out$sleep[-ow]), as.integer(out$state[-ow] == 1L))
})

test_that("detect_sleep_crespo errors on missing required columns", {
  bad <- data.frame(activity = 1:10, datetime = Sys.time() + 1:10)  # no `state`
  expect_error(detect_sleep_crespo(bad))
})

# ── Internal helpers ──────────────────────────────────────────────────────────

test_that(".crespo_msp returns a well-formed binary detection vector", {
  msp <- getFromNamespace(".crespo_msp", "zeitR")
  act <- as.double(rep(c(rep(0, 30), rep(500, 30)), length.out = 400))

  out <- msp(
    activity = act, epoch_h = 2, median_filter_h = 8, pad_h = 1,
    sleep_quantile = 1 / 3, morph_size = 61L,
    consec_zeros_thr = 15L, awake_zeros_thr = 2L, sleep_zeros_thr = 30L,
    zero_mitigation_q = 0.33, min_short_window_thr = 1.0
  )

  expect_length(out, length(act))
  expect_true(is.integer(out))
  expect_true(all(out %in% c(0L, 1L)))
})

test_that(".morphological_open_close uses scipy-style (border_value = 0) padding", {
  morph <- getFromNamespace(".morphological_open_close", "zeitR")

  # With erosion padded by 0 (the fix), borders are NOT protected. The previous
  # convention (erosion padded by 1) returned all 1s for this input.
  x <- c(0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 0L, 1L, 0L)
  expect_identical(
    morph(x, 3L),
    c(0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L)
  )

  # Degenerate inputs: all-zeros is idempotent under closing + opening, but
  # all-ones is NOT — scipy-style erosion (border_value = 0) eats hws epochs in
  # from each border during the closing step. This is the validated reference
  # behaviour (it is what makes .crespo_msp match the Python detect_msp).
  expect_identical(morph(rep(0L, 20), 5L), rep(0L, 20))
  expect_identical(morph(rep(1L, 20), 5L), c(0L, 0L, rep(1L, 16), 0L, 0L))
})

test_that(".estimate_epoch_h recovers the epochs-per-hour from timestamps", {
  est <- getFromNamespace(".estimate_epoch_h", "zeitR")

  t60 <- as.POSIXct("2021-05-27 00:00:00", tz = "UTC") + (0:99) * 60   # 60 s
  expect_equal(est(t60), 60)

  t30 <- as.POSIXct("2021-05-27 00:00:00", tz = "UTC") + (0:99) * 30   # 30 s
  expect_equal(est(t30), 120)
})
