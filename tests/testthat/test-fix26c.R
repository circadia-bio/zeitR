# test-fix26c.R
# Synthetic unit test for Fix 26c: verifies that .recover_fragmented_episodes()
# uses Cole-Kripke epoch scoring on ZCMn to identify sleep runs, NOT the
# period-level CSPD state column.
#
# Background: before Fix 26c, R incorrectly used the period-level CSPD `state`
# to find sleep runs inside the recovery window. Because CSPD marks every epoch
# inside a detected sleep PERIOD as state=1, this collapsed an entire 19+ h
# CSPD period into one continuous sleep run, making recovery produce implausibly
# long episodes instead of the correct bounded fragment.
#
# Scenario
# --------
# Recording: 2020-01-01 00:00 – 2020-01-02 11:59 (36 h, 1-min epochs, UTC)
# CSPD state: state=1 (sleep period) from 14:00 Jan 1 to 09:00 Jan 2 (19 h).
# ZCMn:       two sleep runs separated by a 30-min warm/dark wake gap:
#               run 1  22:00 Jan 1 – 01:00 Jan 2  (3 h, ZCMn = 0)
#               gap    01:00 – 01:30 Jan 2        (30 min, ZCMn = 500,
#                                                  int_temp = 31, light = 5)
#               run 2  01:30 – 07:00 Jan 2        (5.5 h, ZCMn = 0)
#
# Expected outcome
# ----------------
# Fix 26c merges run1 + run2 via CK scoring -> TBT ~ 9 h (< 14 * 60 min).
# If the bug were present (state used instead of CK), TBT would be ~ 19 h
# (> 14 * 60 min), which the test would catch.

test_that("Fix 26c uses Cole-Kripke epoch scoring, not CSPD period-level state", {
  # ── Synthetic recording ───────────────────────────────────────────────────
  t0 <- as.POSIXct("2020-01-01 00:00:00", tz = "UTC")
  n  <- 36L * 60L                              # 2160 1-min epochs
  dt <- t0 + 60L * seq.int(0L, n - 1L)

  # R 1-indexed position of decimal hour h past midnight on 2020-01-01
  idx <- function(h) as.integer(h * 60) + 1L

  zcm   <- rep(100.0, n)   # default: active  -> CK wake
  state <- rep(0L,    n)   # default: CSPD awake
  temp  <- rep(25.0,  n)   # default: cool (off wrist)
  lux   <- rep(500.0, n)   # default: bright (daytime)

  # Sleep run 1: 22:00 Jan 1 -- 01:00 Jan 2 (3 h, ZCMn = 0 -> CK sleep)
  zcm[idx(22):(idx(25) - 1L)] <- 0.0

  # Wake gap: 01:00 -- 01:30 Jan 2 (30 min, ZCMn = 500 -> CK wake)
  # Warm + dark: meets temp/light merge condition for Fix 26c
  zcm [idx(25):(idx(25.5) - 1L)] <- 500.0
  temp[idx(25):(idx(25.5) - 1L)] <- 31.0
  lux [idx(25):(idx(25.5) - 1L)] <- 5.0

  # Sleep run 2: 01:30 -- 07:00 Jan 2 (5.5 h, ZCMn = 0 -> CK sleep)
  zcm[idx(25.5):(idx(31) - 1L)] <- 0.0

  # Bloated CSPD sleep period: 14:00 Jan 1 -- 09:00 Jan 2 (19 h, state = 1).
  # This simulates what the CSPD refiner sees before Fix 26c classification.
  # If .recover_fragmented_episodes() used `state` instead of CK scoring,
  # it would identify this entire 19-h block as one sleep run and inject
  # a TBT > 14 * 60 min -- caught by expect_lt() below.
  state[idx(14):(idx(33) - 1L)] <- 1L

  data <- data.frame(
    datetime = dt,
    ZCMn     = zcm,
    state    = state,
    int_temp = temp,
    light    = lux,
    stringsAsFactors = FALSE
  )

  # One existing main episode on 2020-01-02 (covered); 2020-01-01 is uncovered.
  episodes <- data.frame(
    bts        = as.POSIXct("2020-01-02 23:00:00", tz = "UTC"),
    gts        = as.POSIXct("2020-01-03 07:00:00", tz = "UTC"),
    tbt        = 480.0, tst = 460.0, sol = 5.0, soi = 5.0,
    waso       = 10.0,  nw  = 2L,    eff = 0.96,
    nap        = FALSE,
    sleep_date = as.Date("2020-01-02"),
    sleep_type = "main",
    stringsAsFactors = FALSE
  )

  result <- .recover_fragmented_episodes(
    data         = data,
    episodes     = episodes,
    epoch_min    = 1.0,
    temp_thresh  = 28.0,
    light_thresh = 10.0,
    min_tib_h    = 3.0,
    noc_start    = 18.0,
    noc_end      = 6.0,
    tz           = "UTC"
  )

  recovered <- result[as.character(result$sleep_date) == "2020-01-01", ,
                      drop = FALSE]

  # Exactly one episode should be recovered for the uncovered date
  expect_equal(nrow(recovered), 1L)

  # TBT must come from the merged CK runs (~9 h = ~540 min), not from the
  # bloated CSPD state window (~19 h = ~1140 min).
  # Regression check: if state were used, TBT would exceed 14 * 60.
  expect_lt(recovered$tbt, 14 * 60,
            label = "TBT from CK-merged runs, not 19-h CSPD state block")

  # Sanity lower bound: merged TBT must span both sleep runs (> 8 h = 480 min)
  expect_gt(recovered$tbt, 8 * 60,
            label = "TBT spans both merged CK sleep runs")
})


# ── Fix 25 / Fix 26c interaction ─────────────────────────────────────────────
# Regression test for the "8 main sleep episodes on a 7-day recording" bug
# reported by Julia (JRSV Fix 29f): Fix 25 correctly excludes an episode
# starting at/after noon on the recording's last calendar day (truncated by
# the end of the file, not a real wake-up) at the classify_sleep_episodes()
# level. But .recover_fragmented_episodes() previously had no knowledge of
# that boundary: once Fix 25 removed the episode, its date became
# "uncovered", and the recovery scan below reconstructed the *same* episode
# from the *same* raw epochs -- silently undoing Fix 25.
#
# Scenario (mirrors the notebook's ID_0138 case): a recording whose last
# calendar day has a short evening sleep-like run (19:23-23:56, TBT ~4.5 h)
# right where the file ends -- device removed, not a real wake-up. No other
# candidate episode covers that date, so it reaches the recovery scan.

test_that(".recover_fragmented_episodes() does not recover an episode that would be truncated on the recording's last day", {
  t0 <- as.POSIXct("2020-01-07 00:00:00", tz = "UTC")
  n  <- 24L * 60L                              # 1440 1-min epochs (one day)
  dt <- t0 + 60L * seq.int(0L, n - 1L)

  idx <- function(h) as.integer(h * 60) + 1L    # 1-indexed position of hour h

  zcm   <- rep(100.0, n)   # default: active -> CK wake
  state <- rep(0L,    n)
  temp  <- rep(25.0,  n)
  lux   <- rep(500.0, n)

  # Sleep-like run: 19:23 -- 23:56 (device removed at file end, not a real
  # wake-up). ZCMn = 0 -> CK sleep for the whole run.
  zcm[idx(19 + 23 / 60):n] <- 0.0

  data <- data.frame(
    datetime = dt,
    ZCMn     = zcm,
    state    = state,
    int_temp = temp,
    light    = lux,
    stringsAsFactors = FALSE
  )

  # No episodes cover 2020-01-07 -- as would happen if Fix 25 had just
  # excluded the only candidate episode for that date.
  episodes <- data.frame(
    bts        = as.POSIXct(character(0L)),
    gts        = as.POSIXct(character(0L)),
    tbt        = double(0L), tst = double(0L), sol = double(0L),
    soi        = double(0L), waso = double(0L), nw = integer(0L),
    eff        = double(0L), nap = logical(0L),
    sleep_date = as.Date(character(0L)),
    sleep_type = character(0L),
    stringsAsFactors = FALSE
  )

  result <- .recover_fragmented_episodes(
    data         = data,
    episodes     = episodes,
    epoch_min    = 1.0,
    temp_thresh  = 28.0,
    light_thresh = 10.0,
    min_tib_h    = 3.0,
    noc_start    = 18.0,
    noc_end      = 6.0,
    tz           = "UTC"
  )

  # No episode should be recovered for 2020-01-07: the only candidate starts
  # at 19:23, after noon on the recording's last (and only) calendar day.
  recovered <- result[as.character(result$sleep_date) == "2020-01-07", ,
                      drop = FALSE]
  expect_equal(nrow(recovered), 0L,
               label = "no recovery for an episode that would be truncated on the last day")
})
