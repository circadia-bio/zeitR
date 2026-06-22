# Tests for compute_waso and its faithful nights_df port.

test_that(".nights_df pairs and length-filters sleep periods", {
  ndf <- getFromNamespace(".nights_df", "zeitR")

  # Two 200-epoch sleep periods separated by exactly 60 wake epochs.
  # The merge rule is `gap < wake_thresh`, so a 60-epoch gap does NOT merge.
  states <- c(rep(0L, 60), rep(1L, 200), rep(0L, 60), rep(1L, 200), rep(0L, 60))
  nd <- ndf(states, wake_thresh = 60L)

  expect_equal(nrow(nd), 2L)
  expect_equal(nd$bt, c(60L, 320L))    # 0-indexed first sleep epoch
  expect_equal(nd$gt, c(260L, 520L))   # 0-indexed first wake after (exclusive)
  expect_equal(nd$nap, c(FALSE, FALSE))
})

test_that(".nights_df merges sleep periods closer than wake_thresh", {
  ndf <- getFromNamespace(".nights_df", "zeitR")

  # 59 wake epochs between periods -> 59 < 60 -> merge into one night.
  states <- c(rep(0L, 60), rep(1L, 200), rep(0L, 59), rep(1L, 200), rep(0L, 60))
  nd <- ndf(states, wake_thresh = 60L)

  expect_equal(nrow(nd), 1L)
  expect_equal(nd$bt, 60L)
  expect_equal(nd$gt, 519L)
})

test_that(".nights_df drops short non-nap periods but rescues all-nap runs", {
  ndf <- getFromNamespace(".nights_df", "zeitR")

  # A 30-epoch all-7 nap (>= nap_thresh = 20) and a 200-epoch sleep period,
  # plus a 30-epoch ordinary-sleep run that is too short to keep.
  states <- c(rep(0L, 60), rep(7L, 30), rep(0L, 60), rep(1L, 30),
              rep(0L, 60), rep(1L, 200), rep(0L, 60))
  nd <- ndf(states, wake_thresh = 60L)

  # nap kept (all 7, len 30 >= 20); short sleep run dropped (< 120); long kept.
  expect_equal(nrow(nd), 2L)
  expect_equal(nd$nap, c(TRUE, FALSE))
  expect_equal(nd$bt, c(60L, 240L))
  expect_equal(nd$gt, c(90L, 440L))
})

test_that(".nights_df handles records starting or ending asleep", {
  ndf <- getFromNamespace(".nights_df", "zeitR")

  # Starts asleep: a leading rising edge is synthesised at index 0.
  starts <- c(rep(1L, 200), rep(0L, 100))
  nd_s <- ndf(starts, wake_thresh = 60L)
  expect_equal(nd_s$bt, 0L)
  expect_equal(nd_s$gt, 200L)

  # Ends asleep: a trailing falling edge is synthesised at index n - 1
  # (so the very last epoch is the exclusive end, matching the reference).
  ends <- c(rep(0L, 100), rep(1L, 200))
  nd_e <- ndf(ends, wake_thresh = 60L)
  expect_equal(nd_e$bt, 100L)
  expect_equal(nd_e$gt, 299L)
})

test_that("compute_waso forces within-night WASO-wake epochs to wake", {
  # One 300-epoch sleep period (all on-wrist), ZCM zero everywhere except a
  # single large spike deep inside the night. The Cole-Kripke weighted sum
  # pushes the spike's NEIGHBOURS (not the spike epoch itself) over threshold,
  # producing a deterministic block of WASO-wake epochs inside the night.
  n      <- 420L
  state  <- c(rep(0L, 60), rep(1L, 300), rep(0L, 60))   # night = [60, 360)
  zcm    <- numeric(n)
  zcm[201L] <- 1e6                                       # on-wrist 0-index 200

  x <- tibble::tibble(
    datetime = as.POSIXct("2021-01-01 00:00:00", tz = "UTC") + (seq_len(n) - 1L) * 60,
    ZCMn     = zcm,
    state    = state
  )

  res <- compute_waso(x, wake_thresh = 60L)
  s   <- as.integer(res$data$state)

  # Spike epoch itself stays sleep; its neighbours become wake (the fix: the
  # old code kept the original sleep label on these epochs).
  expect_equal(s[201L], 1L)            # spike epoch -> sleep
  expect_equal(s[196L], 0L)            # within-night neighbour -> WAKE
  expect_equal(s[206L], 0L)            # within-night neighbour -> WAKE
  expect_equal(s[101L], 1L)            # quiet within-night epoch -> sleep
  expect_equal(s[11L],  0L)            # outside any night -> wake

  # Number of wake epochs around the spike: 8 (after-weights) + 9 (before).
  expect_equal(nrow(res$nights), 1L)
  expect_equal(res$nights$tbt,  300L)
  expect_equal(res$nights$waso, 17)
  expect_equal(res$nights$nw,   2L)
  expect_equal(res$nights$sol,  0L)
  expect_equal(res$nights$soi,  0L)
  expect_equal(res$nights$tst,  283)
})

test_that("compute_waso preserves off-wrist epochs as state 4", {
  n     <- 480L
  # off-wrist block in the middle; sleep period on each on-wrist side.
  state <- c(rep(0L, 30), rep(1L, 150), rep(4L, 120), rep(1L, 150), rep(0L, 30))
  zcm   <- numeric(n)

  x <- tibble::tibble(
    datetime = as.POSIXct("2021-01-01 00:00:00", tz = "UTC") + (seq_len(n) - 1L) * 60,
    ZCMn     = zcm,
    state    = state
  )

  res <- compute_waso(x, wake_thresh = 60L)
  s   <- as.integer(res$data$state)

  expect_true(all(s[state == 4L] == 4L))   # off-wrist untouched
  expect_true(all(s[state != 4L] %in% c(0L, 1L)))
  # sleep column: 1 where state in {1,7}, else 0
  expect_equal(as.integer(res$data$sleep), as.integer(s == 1L | s == 7L))
})
