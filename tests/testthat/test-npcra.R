# Tests for npcra.R: compute_npcra() (non-parametric circadian rhythm
# analysis: IS, IV, RA, L5, M10).

# ---- Shared fixture ---------------------------------------------------------
# A perfectly repeating 24h profile across 7 days at 1-minute epochs:
#   00:00-10:00  activity = 0    (10h "rest" block)
#   10:00-24:00  activity = 100  (14h "active" block)
# Because the pattern is byte-identical every day, this gives exact,
# hand-checkable values: IS = 1 (perfect between-day consistency), L5 = 0
# (the least-active 5h window sits entirely inside the rest block), M10 = 100
# (the most-active 10h window sits entirely inside the active block), and
# RA = (100 - 0) / (100 + 0) = 1.

make_npcra_fixture <- function(n_days = 7L) {
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  n  <- n_days * 24L * 60L
  dt <- t0 + 60L * seq.int(0L, n - 1L)

  hour     <- as.integer(format(dt, "%H", tz = "UTC"))
  activity <- ifelse(hour < 10L, 0, 100)

  tibble::tibble(datetime = dt, activity = activity)
}

test_that("compute_npcra() recovers exact IS/L5/M10/RA for a perfectly repeating profile", {
  d      <- make_npcra_fixture()
  result <- compute_npcra(d)

  expect_equal(result$IS, 1)
  expect_equal(result$L5, 0)
  expect_equal(result$M10, 100)
  expect_equal(result$RA, 1)
  expect_true(is.finite(result$IV))
  expect_gte(result$IV, 0)
})

test_that("compute_npcra() L5/M10 onset times are exact for a perfectly repeating profile", {
  d      <- make_npcra_fixture()
  result <- compute_npcra(d)

  # Onset is the END of the winning window (matches the production Python's
  # right-labelled pandas rolling().mean() convention -- see .lmx_window()),
  # NOT the start. Ties across the 6 identical trimmed days are broken by
  # which.min()/which.max() picking the FIRST occurrence, so the earliest
  # fully-qualifying window on day 2 (the first day after the default D+1
  # trim) wins deterministically.
  #
  # L5: earliest 5h window fully inside the 00:00-10:00 rest block is
  # [00:00, 05:00) -> onset (window end) = 05:00.
  expect_equal(result$L5_onset, "05:00")

  # M10: earliest 10h window fully inside the 10:00-24:00 active block is
  # [10:00, 20:00) -> onset (window end) = 20:00.
  expect_equal(result$M10_onset, "20:00")
})

test_that("compute_npcra() excludes off-wrist epochs when a state column is present", {
  d <- make_npcra_fixture()
  # Mark the first day entirely off-wrist with an implausible spike; if it
  # were NOT excluded this would corrupt IS/L5/M10 away from the exact values
  # checked above.
  d$state <- 0L
  d$state[1:1440] <- 4L
  d$activity[1:1440] <- 99999

  result <- compute_npcra(d)

  expect_equal(result$IS, 1)
  expect_equal(result$L5, 0)
  expect_equal(result$M10, 100)
})

test_that("compute_npcra() accepts a zeitr_recording and uses its participant_id", {
  d   <- make_npcra_fixture()
  rec <- structure(
    list(epochs = d, metadata = list(participant_id = "P042")),
    class = "zeitr_recording"
  )

  result <- compute_npcra(rec)
  expect_equal(result$participant_id, "P042")
})

test_that("compute_npcra() participant_id is NA for a bare data frame", {
  d      <- make_npcra_fixture()
  result <- compute_npcra(d)
  expect_true(is.na(result$participant_id))
})

test_that("compute_npcra() errors on missing required columns", {
  bad <- tibble::tibble(datetime = Sys.time())
  expect_error(compute_npcra(bad), "Missing required column")
})

test_that("compute_npcra() errors when x is neither a zeitr_recording nor a data frame", {
  expect_error(compute_npcra(list(a = 1)), "must be a")
  expect_error(compute_npcra("not a recording"), "must be a")
})

test_that("compute_npcra() errors with fewer than 2 epochs", {
  d <- make_npcra_fixture()[1, ]
  expect_error(compute_npcra(d), "at least 2 epochs")
})

test_that("compute_npcra(window_days = ...) splits into the expected number of windows", {
  d      <- make_npcra_fixture(n_days = 6L)
  result <- compute_npcra(d, window_days = 2)

  expect_true("window_start" %in% names(result))
  expect_equal(nrow(result), 3L)
  # Each 2-day window should still show the perfect-repeat values
  expect_equal(result$IS, rep(1, 3L))
})

test_that("compute_npcra() trims to D+1 00:00 by default", {
  d <- make_npcra_fixture(n_days = 7L)   # starts exactly at 2024-01-01 00:00

  result_default <- compute_npcra(d)                       # trim_to_d1 = TRUE
  result_notrim   <- compute_npcra(d, trim_to_d1 = FALSE)

  # Trimming removes exactly one full day (1440 epochs at 1-min resolution)
  expect_equal(result_notrim$n_epochs - result_default$n_epochs, 1440L)
  expect_equal(result_default$n_days, 6)
  expect_equal(result_notrim$n_days, 7)

  # The fixture repeats identically every day, so trimming one day off
  # should not change IS/L5/M10/RA at all.
  expect_equal(result_default$IS,  result_notrim$IS)
  expect_equal(result_default$L5,  result_notrim$L5)
  expect_equal(result_default$M10, result_notrim$M10)
  expect_equal(result_default$RA,  result_notrim$RA)
})

test_that("compute_npcra() falls back to the untrimmed recording with a warning when trim_to_d1 leaves <2 epochs", {
  # A recording that stays entirely within a single calendar day (08:00 --
  # 23:59, never crossing midnight) has nothing left after trimming to
  # D+1 00:00.
  t0 <- as.POSIXct("2024-01-01 08:00:00", tz = "UTC")
  d  <- tibble::tibble(
    datetime = t0 + 60L * 0:959,    # 2024-01-01 08:00 -- 2024-01-01 23:59
    activity = rep(c(0, 100), each = 480L)
  )

  expect_warning(
    result <- compute_npcra(d, trim_to_d1 = TRUE),
    "fewer than 2 epochs"
  )
  result_notrim <- compute_npcra(d, trim_to_d1 = FALSE)
  expect_equal(result$n_epochs, result_notrim$n_epochs)
})

test_that("compute_npcra() respects a manually supplied epoch_s", {
  d <- make_npcra_fixture()
  # Supplying the correct epoch length explicitly should give the same
  # result as letting it auto-estimate from the datetime column.
  result_auto <- compute_npcra(d)
  result_manual <- compute_npcra(d, epoch_s = 60)
  expect_equal(result_auto$IS, result_manual$IS)
  expect_equal(result_auto$n_epochs, result_manual$n_epochs)
})
