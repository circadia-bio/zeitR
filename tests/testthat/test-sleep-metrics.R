# Tests for sleep_metrics.R: compute_sleep_metrics() and compute_cpd_metrics().
#
# Regression coverage for three related fixes:
#   1. sleep_offset_h (and sleep_onset_h) must use circ_mean_h() (the
#      package's own exported circular mean, circ_utils.R), not a plain
#      mean(). A previous internal duplicate (.mean_circ_h(), since removed)
#      drifted into a "shift values >= 12h by -24, then plain-average"
#      shortcut that is NOT mathematically equivalent to the true
#      atan2-based circular mean -- confirmed to diverge even with zero
#      midnight wraparound.
#   2. dp_midsleep_min must use circ_sd_h() (circular SD, R-bar formula),
#      not a plain sd().
#   3. compute_cpd_metrics()'s sjla_h must be clamped to [-12, 12] before
#      abs() is taken for sjl_h.
#
# circ_mean_h()/circ_sd_h() themselves are already thoroughly unit-tested in
# test-circ-utils.R; these tests check WIRING -- that compute_sleep_metrics()/
# compute_cpd_metrics() actually call them -- not the underlying math.
#
# The sleep_offset_h values in the first test are the exact 7 real values
# from a real cohort recording, verified against a real Python execution
# earlier (circular mean 7.665164h vs a plain/shift-based 7.689325h). Only
# the derived decimal-hour values are reused here (not any raw participant
# data); they're baked into synthetic get_up_time timestamps on an arbitrary
# reference date, not the real recording dates.

test_that("compute_sleep_metrics(): sleep_offset_h uses the circular mean, not plain/shift-based", {
  offset_hours <- c(7.37027778, 10.03694444, 9.43694444, 6.77027778,
                    6.35361111, 6.90361111, 6.95361111)

  ref_date <- as.Date("2024-06-01")
  get_up   <- as.POSIXct(ref_date, tz = "UTC") + offset_hours * 3600

  nd <- data.frame(
    is_nap      = rep(FALSE, 7L),
    bed_time    = as.POSIXct(ref_date - 1L, tz = "UTC") + 22 * 3600,
    get_up_time = get_up,
    tbt         = rep(480.0, 7L),   # 8h, comfortably over min_tib_h default
    tst         = rep(440.0, 7L),
    sol         = rep(0.0,   7L),   # sol/soi = 0 -> offset_h == get_up clock hour
    soi         = rep(0.0,   7L),
    waso        = rep(20.0,  7L),
    eff         = rep(0.9,   7L)
  )

  result <- suppressWarnings(compute_sleep_metrics(nd, tz = "UTC"))

  # Run offset_hours through the SAME .to_h_plain() round-trip the function
  # itself uses (whole-second format() extraction) before computing the
  # expected circular mean -- comparing against circ_mean_h(offset_hours)
  # directly on the original decimal literals introduces sub-second
  # rounding noise that isn't present in the function's own pathway (these
  # literals don't all correspond to exact whole seconds).
  offset_h_as_fn_sees <- .to_h_plain(get_up, tz = "UTC")
  expected_circular    <- circ_mean_h(offset_h_as_fn_sees)
  plain_mean           <- mean(offset_h_as_fn_sees)

  # Sanity: confirm this fixture genuinely distinguishes the two methods --
  # if a future change to the fixture accidentally made them coincide, the
  # test below would pass for the wrong reason.
  expect_gt(abs(expected_circular - plain_mean), 1e-4)

  expect_equal(result$sleep_offset_h, expected_circular, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(result$sleep_offset_h, plain_mean, tolerance = 1e-4)))
})

test_that("compute_sleep_metrics(): dp_midsleep_min uses the circular SD, not a plain SD", {
  offset_hours <- c(7.37027778, 10.03694444, 9.43694444, 6.77027778,
                    6.35361111, 6.90361111, 6.95361111)
  onset_hours  <- rep(22.0, 7L)   # bed_time clock hour, sol = 0

  ref_date <- as.Date("2024-06-01")
  get_up   <- as.POSIXct(ref_date, tz = "UTC") + offset_hours * 3600
  bed      <- as.POSIXct(ref_date - 1L, tz = "UTC") + onset_hours * 3600

  nd <- data.frame(
    is_nap      = rep(FALSE, 7L),
    bed_time    = bed,
    get_up_time = get_up,
    tbt         = rep(480.0, 7L),
    tst         = rep(440.0, 7L),
    sol         = rep(0.0,   7L),
    soi         = rep(0.0,   7L),
    waso        = rep(20.0,  7L),
    eff         = rep(0.9,   7L)
  )

  result <- suppressWarnings(compute_sleep_metrics(nd, tz = "UTC"))

  # Same rationale as the sleep_offset_h test above: run the raw hour
  # literals through .to_h_plain() before computing expectations, matching
  # the function's own whole-second round-trip exactly.
  onset_h_as_fn_sees  <- .to_h_plain(bed,    tz = "UTC")
  offset_h_as_fn_sees <- .to_h_plain(get_up, tz = "UTC")
  mid_h               <- .midsleep(onset_h_as_fn_sees, offset_h_as_fn_sees)
  expected_circular   <- circ_sd_h(mid_h) * 60
  plain_sd_minutes    <- stats::sd(mid_h) * 60

  expect_gt(abs(expected_circular - plain_sd_minutes), 1e-4)
  expect_equal(result$dp_midsleep_min, expected_circular, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(result$dp_midsleep_min, plain_sd_minutes, tolerance = 1e-4)))
})

test_that("compute_cpd_metrics(): sjla_h is clamped to [-12, 12] before abs()", {
  # Free-day nights with mid-sleep clustered late (~23.5h); workday nights
  # with mid-sleep clustered mid-morning (~10h) -- chosen so the RAW
  # msf_h - msw_h difference comfortably exceeds 12h, forcing the clamp to
  # engage. All dates/times are synthetic, not derived from any real
  # recording.
  nd <- data.frame(
    is_nap      = rep(FALSE, 4L),
    bed_time    = as.POSIXct(c("2024-01-05 22:00:00",   # Fri -> Sat (free, eve)
                               "2024-01-06 22:30:00",   # Sat -> Sun (free)
                               "2024-01-08 08:00:00",   # Mon (workday)
                               "2024-01-09 08:30:00"),  # Tue (workday)
                             tz = "UTC"),
    get_up_time = as.POSIXct(c("2024-01-06 01:00:00",
                               "2024-01-07 00:30:00",
                               "2024-01-08 12:00:00",
                               "2024-01-09 12:30:00"),
                             tz = "UTC"),
    tbt  = rep(200.0, 4L),   # comfortably over min_tib_h/min_tib_eve_h (3h = 180 min)
    tst  = rep(180.0, 4L),
    sol  = rep(0.0,   4L),
    soi  = rep(0.0,   4L),
    waso = rep(10.0,  4L),
    eff  = rep(0.9,   4L)
  )

  result <- suppressWarnings(compute_cpd_metrics(nd, tz = "UTC"))

  # Recompute expected msf_h/msw_h the same way the function does (via
  # .midsleep() + circ_mean_h()) -- this checks the clamp logic
  # independently, not just trusting the function's own msf_h/msw_h.
  onset_h  <- c(22, 22.5, 8, 8.5)
  offset_h <- c(1, 0.5, 12, 12.5)
  mid_h    <- .midsleep(onset_h, offset_h)

  is_free <- c(TRUE, TRUE, FALSE, FALSE)
  msf_h   <- circ_mean_h(mid_h[is_free])
  msw_h   <- circ_mean_h(mid_h[!is_free])

  raw_diff <- msf_h - msw_h
  expect_gt(abs(raw_diff), 12,
            label = "fixture must produce a raw msf_h - msw_h exceeding 12h to exercise the clamp")

  expected_sjla <- raw_diff
  if (expected_sjla > 12)  expected_sjla <- expected_sjla - 24
  if (expected_sjla < -12) expected_sjla <- expected_sjla + 24
  expected_sjl <- abs(expected_sjla)

  expect_equal(result$msf_h, msf_h, tolerance = 1e-8)
  expect_equal(result$msw_h, msw_h, tolerance = 1e-8)
  expect_equal(result$sjla_h, expected_sjla, tolerance = 1e-8)
  expect_equal(result$sjl_h,  expected_sjl,  tolerance = 1e-8)

  # Regression guard: the un-clamped value would have been the raw
  # difference itself -- confirm the function's output is NOT that.
  expect_false(isTRUE(all.equal(result$sjla_h, raw_diff, tolerance = 1e-4)))
})
