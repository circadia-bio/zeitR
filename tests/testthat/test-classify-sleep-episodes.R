# Tests for sleep_classify.R: classify_sleep_episodes().
#
# Regression coverage for removing the classification-stage Fix 25 filter
# (exclude episodes starting at/after noon on the recording's last calendar
# day). It used to run as classify_sleep_episodes()'s very first step;
# removed after confirming the actual production Python (Cell 3 of the
# fix29 notebook -- the real classification logic, not the superseded
# shared pipeline_functions.py) has no equivalent step at the
# classification stage at all. Fix 25 only exists there later, specific to
# the CPD calculation, which compute_cpd_metrics() already mirrors
# correctly on the R side (see test-sleep-metrics.R and sleep_metrics.R's
# module header).
#
# .recover_fragmented_episodes() (Fix 26c) keeps its own independent
# last-day-noon guard -- test-fix26c.R's second test already exercises that
# directly and confirms it alone is sufficient to prevent the "8 main
# nights on a 7-day recording" bug (Fix 29f) the classification-stage
# filter was apparently added in response to.

test_that("classify_sleep_episodes() no longer excludes a complete episode starting late on the recording's last day", {
  # Single-day recording, 1-min epochs, 2020-01-08 00:00 -- 23:59. One
  # legitimate main sleep episode: 18:00-23:50 (5h50, well inside the
  # default nocturnal window and above min_main_tib_h), completing 9
  # minutes before the recording's own last epoch -- clearly not
  # truncated, just late in the day.
  #
  # Under the old classification-stage Fix 25, this episode's bts (18:00)
  # falls on the same calendar date as the recording's last epoch (23:59)
  # and is >= noon, so it was excluded before ever being classified --
  # even though the data proves it's complete.
  dt <- as.POSIXct("2020-01-08 00:00:00", tz = "UTC") + 60L * seq.int(0L, 1439L)

  data <- data.frame(
    datetime = dt,
    state    = rep(0L, 1440L),
    activity = rep(0.0, 1440L)
  )

  episodes <- data.frame(
    bts  = as.POSIXct("2020-01-08 18:00:00", tz = "UTC"),
    gts  = as.POSIXct("2020-01-08 23:50:00", tz = "UTC"),
    tbt  = 350.0, tst = 330.0, sol = 5.0, soi = 5.0,
    waso = 10.0,  nw  = 2L,    eff = 0.94,
    nap  = FALSE
  )

  result <- classify_sleep_episodes(episodes, data)

  expect_equal(nrow(result), 1L)
  expect_equal(result$sleep_type, "main")
  expect_false(result$is_nap)
})

test_that("classify_sleep_episodes() still classifies a normal overnight episode as main (unaffected by the Fix 25 removal)", {
  # Sanity check: a completely ordinary overnight episode, nowhere near the
  # recording's last-day boundary, should be unaffected by removing Fix 25
  # from classify_sleep_episodes().
  dt <- as.POSIXct("2020-01-01 00:00:00", tz = "UTC") + 60L * seq.int(0L, 4L * 1440L - 1L)

  data <- data.frame(
    datetime = dt,
    state    = rep(0L, length(dt)),
    activity = rep(0.0, length(dt))
  )

  episodes <- data.frame(
    bts  = as.POSIXct("2020-01-02 23:00:00", tz = "UTC"),
    gts  = as.POSIXct("2020-01-03 07:00:00", tz = "UTC"),
    tbt  = 480.0, tst = 460.0, sol = 5.0, soi = 5.0,
    waso = 10.0,  nw  = 2L,    eff = 0.96,
    nap  = FALSE
  )

  result <- classify_sleep_episodes(episodes, data)

  expect_equal(nrow(result), 1L)
  expect_equal(result$sleep_type, "main")
})
