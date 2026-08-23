# Tests for sleep_classify.R: classify_sleep_episodes().
#
# Regression coverage for Fix 29f (exclude truncated episodes at recording
# end, applied during classification itself).
#
# History: this filter used to be removed from classify_sleep_episodes()
# entirely, based on the fix29 notebook's Cell 3, which genuinely has no
# classification-stage truncation filter. That removal was wrong: fix30
# (the actual current production notebook) re-adds exactly this filter as
# "Fix 29f", with its own comment explaining the real bug it fixes
# (ID_0138: an 8th spurious "night" out of a 7-day recording, a same-day
# 19:23->23:56 artifact miscounted as a valid main night). Consistent with
# this, removing the filter earlier did not improve the cohort-wide n_main
# match rate against the real Python reference at all (stayed ~77-78%).
#
# The "last day" reference is max(gts) among the episodes themselves
# (Python: pd.to_datetime(nights_data['gts']).max().normalize()), NOT the
# raw epoch data's last timestamp -- these usually coincide but can
# diverge, and the tests below are constructed to make that distinction
# visible.

test_that("Fix 29f excludes a spurious same-day artifact night (ID_0138-style)", {
  # Two nights: a normal one crossing midnight (day1 23:00 -> day2 07:00,
  # TBT 8h, clearly a real main night), and a same-day artifact fully
  # contained within day2 (19:23 -> 23:56, TBT ~4.5h) -- matching the
  # exact shape of the confirmed real bug (ID_0138). The artifact's bts is
  # on the SAME calendar date as max(gts) (day2, from either night) and at
  # >= noon, so it should be excluded before classification ever runs.
  dt <- as.POSIXct("2020-01-01 00:00:00", tz = "UTC") + 60L * seq.int(0L, 2L * 1440L - 1L)
  data <- data.frame(
    datetime = dt,
    state    = rep(0L, length(dt)),
    activity = rep(0.0, length(dt))
  )

  episodes <- data.frame(
    bts  = as.POSIXct(c("2020-01-01 23:00:00", "2020-01-02 19:23:00"), tz = "UTC"),
    gts  = as.POSIXct(c("2020-01-02 07:00:00", "2020-01-02 23:56:00"), tz = "UTC"),
    tbt  = c(480.0, 273.0), tst = c(460.0, 260.0),
    sol  = c(5.0, 5.0),     soi = c(5.0, 5.0),
    waso = c(10.0, 5.0),    nw  = c(2L, 1L),
    eff  = c(0.96, 0.95),   nap = c(FALSE, FALSE)
  )

  result <- classify_sleep_episodes(episodes, data)

  # Only the normal night survives; the artifact is gone entirely (not
  # main, not secondary -- excluded before classification).
  expect_equal(nrow(result), 1L)
  expect_equal(result$bts, as.POSIXct("2020-01-01 23:00:00", tz = "UTC"))
  expect_equal(result$sleep_type, "main")
})

test_that("Fix 29f does not exclude a normal single night crossing midnight", {
  # A single ordinary night: bts on day1, gts on day2. max(gts) = day2, so
  # the noon threshold is day2 12:00 -- bts (day1 23:00) is comfortably
  # before that, regardless of which reference (raw epoch data or episode
  # gts) is used. Sanity check that Fix 29f doesn't over-trigger on the
  # single most common case.
  dt <- as.POSIXct("2020-01-01 00:00:00", tz = "UTC") + 60L * seq.int(0L, 2L * 1440L - 1L)
  data <- data.frame(
    datetime = dt,
    state    = rep(0L, length(dt)),
    activity = rep(0.0, length(dt))
  )

  episodes <- data.frame(
    bts  = as.POSIXct("2020-01-01 23:00:00", tz = "UTC"),
    gts  = as.POSIXct("2020-01-02 07:00:00", tz = "UTC"),
    tbt  = 480.0, tst = 460.0, sol = 5.0, soi = 5.0,
    waso = 10.0,  nw  = 2L,    eff = 0.96, nap = FALSE
  )

  result <- classify_sleep_episodes(episodes, data)

  expect_equal(nrow(result), 1L)
  expect_equal(result$sleep_type, "main")
})

test_that("Fix 29f's 'last day' reference is max(gts) among episodes, not the raw epoch data's last timestamp", {
  # Constructs a case where the two reference points would give DIFFERENT
  # answers, to confirm the episode-gts-based reference is what's actually
  # used (matching fix30's pd.to_datetime(nights_data['gts']).max()), not
  # the raw recording's last epoch.
  #
  # Recording spans through day3 12:00 (raw data continues well past the
  # last classified night's gts) -- so "last day" by raw-epoch-max would
  # be day3, but max(gts) among episodes is day2. A candidate episode
  # starting day2 20:00 (after noon on day2, the episode-gts-based last
  # day) should be EXCLUDED under the correct (episode-gts) reference,
  # even though day2 12:00 is comfortably BEFORE day3 (the raw-epoch-max
  # reference), which would NOT have excluded it.
  dt <- as.POSIXct("2020-01-01 00:00:00", tz = "UTC") +
    60L * seq.int(0L, as.integer(difftime(
      as.POSIXct("2020-01-03 12:00:00", tz = "UTC"),
      as.POSIXct("2020-01-01 00:00:00", tz = "UTC"),
      units = "mins")))
  data <- data.frame(
    datetime = dt,
    state    = rep(0L, length(dt)),
    activity = rep(0.0, length(dt))
  )

  episodes <- data.frame(
    bts  = as.POSIXct(c("2020-01-01 23:00:00", "2020-01-02 20:00:00"), tz = "UTC"),
    gts  = as.POSIXct(c("2020-01-02 07:00:00", "2020-01-02 23:30:00"), tz = "UTC"),
    tbt  = c(480.0, 210.0), tst = c(460.0, 200.0),
    sol  = c(5.0, 5.0),     soi = c(5.0, 5.0),
    waso = c(10.0, 5.0),    nw  = c(2L, 1L),
    eff  = c(0.96, 0.95),   nap = c(FALSE, FALSE)
  )

  result <- classify_sleep_episodes(episodes, data)

  # The second episode (day2 20:00) is excluded: max(gts) = day2 (from
  # either episode), and day2 20:00 >= day2 12:00 noon.
  expect_equal(nrow(result), 1L)
  expect_equal(result$bts, as.POSIXct("2020-01-01 23:00:00", tz = "UTC"))
})
