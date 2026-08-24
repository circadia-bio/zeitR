
# Regression tests for two border-refiner internal-helper bugs found via a
# real off-wrist divergence on participant ID_0003 (see
# dev/latencia_investigation_handoff.md for the original investigation).
# Both functions are pure and fully synthetic-testable in isolation, per
# house style (synthetic fixtures preferred over real-data fixtures).

test_that(".compute_forbidden_zone_v2 excludes everything but the +/-1h edges of a sleep block, not a middle quartile", {
  # Real-bug pattern (ID_0003): a genuine off-wrist episode sat near the
  # START edge of a long, mostly-spurious estimated-sleep block. The old R
  # implementation used an invented "middle 25%-75% quantile" heuristic that
  # appears nowhere in the Python reference; it marked the block's *middle*
  # as forbidden and left the edges open, which is backwards relative to
  # what the reference implementation (analyze_sleep_borders()) actually
  # does: forbidden_start = min(sleep_start + possible_window_hours*epoch_hour, n-1)
  #       forbidden_end   = max(0, sleep_end - possible_window_hours*epoch_hour)
  #       forbidden_zone[forbidden_start:forbidden_end] = 1  (0-based, exclusive end)
  # i.e. the *edges* (first/last hour) are exempt and everything else in
  # between is forbidden -- the near-total opposite of the old formula.

  n          <- 600L
  epoch_hour <- 60L  # 1-minute epochs -> 60 epochs/hour, matching ActTrust data

  estimated_sleep <- rep(1L, n)
  # One sleep block, Python 0-based [100, 500) -> R indices 101:500
  estimated_sleep[101:500] <- 0L

  fz <- zeitR:::.compute_forbidden_zone_v2(estimated_sleep, epoch_hour, n)

  # Hand-computed expected forbidden region (possible_window_hours = 1, the
  # Python default):
  #   forbidden_start = min(100 + 1*60, 599) = 160
  #   forbidden_end   = max(0, 500 - 1*60)   = 440
  #   forbidden_zone[160:440] = 1  (Python 0-based, exclusive end)
  #   -> R indices (160+1):440 = 161:440
  expect_true(all(fz[161:440] == 1L))
  expect_true(all(fz[1:160] == 0L))
  expect_true(all(fz[441:n] == 0L))

  # A short off-wrist candidate near the block's START edge (Python epochs
  # 110-139, i.e. within the first-hour margin) must NOT be forbidden --
  # this is exactly the ID_0003 pattern that the old quartile heuristic got
  # backwards (it would have marked this region forbidden and the true
  # middle-of-block region as allowed).
  near_start_edge <- fz[111:140]  # Python [110,140) -> R 111:140
  expect_true(all(near_start_edge == 0L))

  # A short off-wrist candidate near the block's END edge (Python epochs
  # 460-489, within the last-hour margin) must also not be forbidden.
  near_end_edge <- fz[461:490]  # Python [460,490) -> R 461:490
  expect_true(all(near_end_edge == 0L))

  # A candidate deep in the middle of the block must be forbidden.
  deep_middle <- fz[251:280]  # Python [250,280)
  expect_true(all(deep_middle == 1L))
})

test_that(".compute_forbidden_zone_v2 leaves short sleep blocks entirely unforbidden", {
  # When a sleep block is shorter than 2*possible_window_hours*epoch_hour,
  # forbidden_end <= forbidden_start in Python, so the slice is empty and
  # nothing in that block is forbidden.
  n          <- 200L
  epoch_hour <- 60L
  estimated_sleep <- rep(1L, n)
  estimated_sleep[51:150] <- 0L  # Python [50,150), length 100 < 2*60

  fz <- zeitR:::.compute_forbidden_zone_v2(estimated_sleep, epoch_hour, n)
  expect_true(all(fz == 0L))
})

test_that(".compute_forbidden_zone_v2 returns all-zero when there is no estimated sleep", {
  n <- 50L
  fz <- zeitR:::.compute_forbidden_zone_v2(rep(1L, n), epoch_hour = 60L, n = n)
  expect_true(all(fz == 0L))
  expect_length(fz, n)
})

test_that(".check_valid_border_mod gates on the fixed low_activity_threshold constant, not a dynamic one", {
  # Python's check_valid_border_features() ("mod" configuration, the default
  # and the one actually used) is:
  #   activity_median_low[i] OR
  #     (temperature[i] < temperature_threshold AND activity_median[i] < 2 * self.low_activity_threshold)
  # where self.low_activity_threshold is a FIXED constructor default (500),
  # never overridden by the ActTrust wrapper -- a completely different
  # quantity from the data-driven self.activity_threshold used elsewhere.
  # An earlier R port passed the dynamic quantity in here by mistake.

  act_med_low <- c(0L, 0L)
  temperature <- c(20, 20)
  temperature_threshold <- 25
  act_median  <- c(900, 1100)  # straddles 2*500 = 1000

  # activity_median = 900 < 1000 -> valid via the temp+activity branch
  expect_true(zeitR:::.check_valid_border_mod(
    act_med_low, temperature, temperature_threshold, act_median,
    low_activity_threshold = 500, ri = 1L
  ))

  # activity_median = 1100 >= 1000 -> not valid (act_med_low is 0, so no
  # other path can succeed)
  expect_false(zeitR:::.check_valid_border_mod(
    act_med_low, temperature, temperature_threshold, act_median,
    low_activity_threshold = 500, ri = 2L
  ))

  # A larger threshold (mimicking the old, incorrect data-driven call) would
  # flip epoch 2 to valid -- confirming the fixed constant, not some other
  # dynamic value, is what actually drives the gate.
  expect_true(zeitR:::.check_valid_border_mod(
    act_med_low, temperature, temperature_threshold, act_median,
    low_activity_threshold = 600, ri = 2L
  ))
})

test_that(".check_valid_border_mod: activity_median_low always short-circuits to valid", {
  act_med_low <- c(1L, 0L)
  temperature <- c(100, 100)  # far above any plausible threshold
  temperature_threshold <- 25
  act_median  <- c(100000, 100000)  # far above 2*low_activity_threshold

  expect_true(zeitR:::.check_valid_border_mod(
    act_med_low, temperature, temperature_threshold, act_median,
    low_activity_threshold = 500, ri = 1L
  ))
  expect_false(zeitR:::.check_valid_border_mod(
    act_med_low, temperature, temperature_threshold, act_median,
    low_activity_threshold = 500, ri = 2L
  ))
})
