# Tests for sleep_classify.R: extract_sleep_episodes().
#
# Regression coverage for the get_up_time epoch-shift bug: gts was reading
# stamps_ow[gt0] (the *last sleep* epoch), one epoch early. The fix reads
# stamps_ow[gt0 + 1L] (the *first wake* epoch), matching Python's
# stamps[gt] exactly (condor_pipeline/detection/waso.py: gts.append(stamps[gt])).
#
# No existing test caught this: test-vallim-parity.R explicitly excludes
# timestamp/boundary checks (documented in its own header as "What is NOT
# tested"), and there was no other direct test of extract_sleep_episodes().
#
# Sleep blocks here are 150 min -- comfortably over .nights_df()'s default
# sleep_thresh = 120 min, so they register as a genuine "main" period
# rather than falling through to (and failing) the state == 7 nap check.
#
# POSIXct comparisons use as.numeric() (underlying instant) rather than
# comparing the POSIXct objects directly, sidestepping tzone-attribute
# bookkeeping noise that isn't the point of this test.

test_that("extract_sleep_episodes() sets gts to the first wake epoch, not the last sleep epoch", {
  t0 <- as.POSIXct("2020-01-01 00:00:00", tz = "UTC")
  n  <- 500L                                    # 500 x 1-min epochs
  dt <- t0 + 60L * seq.int(0L, n - 1L)

  # state == 1 for minutes 100-249 inclusive (a clean, isolated 150-min
  # sleep block -- over the 120-min sleep_thresh -- flanked by long wake
  # stretches on both sides, well clear of wake_thresh = 60L).
  state <- rep(0L, n)
  state[101:250] <- 1L                          # R 1-indexed 101:250
                                                 # = 0-indexed minutes 100:249

  # ZCMn = 0 throughout the sleep block -> Cole-Kripke scores the whole
  # window as sleep (sol = soi = 0), keeping the test focused purely on the
  # epoch-indexing fix rather than CK-scoring nuances.
  zcm <- rep(50.0, n)
  zcm[101:250] <- 0.0

  data <- data.frame(datetime = dt, ZCMn = zcm, state = state)

  result <- extract_sleep_episodes(data, wake_thresh = 60L)

  expect_equal(nrow(result), 1L)

  # Last sleep epoch (0-indexed minute 249) -- the OLD, buggy value.
  last_sleep_epoch <- t0 + 60L * 249L
  # First wake epoch (0-indexed minute 250) -- the CORRECT value, matching
  # Python's stamps[gt].
  first_wake_epoch <- t0 + 60L * 250L

  expect_equal(as.numeric(result$gts), as.numeric(first_wake_epoch))
  expect_false(isTRUE(all.equal(as.numeric(result$gts), as.numeric(last_sleep_epoch))))

  # bts is unaffected by this fix either way -- first sleep epoch (0-indexed
  # minute 100), consistent before and after.
  expect_equal(as.numeric(result$bts), as.numeric(t0 + 60L * 100L))
})

test_that("extract_sleep_episodes() gts epoch fix holds across multiple episodes", {
  # Same idea, but with two separate sleep blocks to confirm the fix applies
  # uniformly, not just to a single isolated episode.
  t0 <- as.POSIXct("2020-01-01 00:00:00", tz = "UTC")
  n  <- 1000L
  dt <- t0 + 60L * seq.int(0L, n - 1L)

  state <- rep(0L, n)
  state[101:250] <- 1L   # sleep block 1: minutes 100-249 (150 min)
  state[501:650] <- 1L   # sleep block 2: minutes 500-649 (150 min)

  zcm <- rep(50.0, n)
  zcm[101:250] <- 0.0
  zcm[501:650] <- 0.0

  data <- data.frame(datetime = dt, ZCMn = zcm, state = state)

  result <- extract_sleep_episodes(data, wake_thresh = 60L)

  expect_equal(nrow(result), 2L)
  expect_equal(as.numeric(result$gts[1]), as.numeric(t0 + 60L * 250L))   # first wake after block 1
  expect_equal(as.numeric(result$gts[2]), as.numeric(t0 + 60L * 650L))   # first wake after block 2
})
