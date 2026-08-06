# Tests for npcra.R: compute_npcra() (non-parametric circadian rhythm
# analysis: IS, IV, RA, L5, M10).

# ---- Shared fixture ---------------------------------------------------------
# A perfectly repeating 24h profile across 7 days at 1-minute epochs:
#   00:00-10:00  activity = 0    (10h "rest" block)
#   10:00-24:00  activity = 100  (14h "active" block)
# Because the pattern is byte-identical every day, this gives exact,
# hand-checkable values for L5/M10/RA -- and for IS, an exact value derived
# from the ddof = 1 formula (see .theoretical_is_perfect() below; it is NOT
# exactly 1, even for a perfectly repeating pattern). L5 = 0 (the
# least-active 5h window sits entirely inside the rest block), M10 = 100
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

# Theoretical IS for a perfectly repeating 24h profile over `n_days` full
# days, under the ddof = 1 (sample variance) formula compute_npcra() now
# uses to match pyActigraphy's actual _interdaily_stability(). Derived
# algebraically: since X is just the 24-value profile Xh tiled n_days
# times, Sum((X-Xm)^2) = n_days * Sum((Xh-Xm)^2) exactly, so the
# Sum((Xh-Xm)^2) factor cancels in the IS ratio, leaving:
#   IS = [1/(p-1)] / [n_days/(N-1)],  N = 24*n_days, p = 24
# Notably this does NOT converge to 1 as n_days -> Inf -- it converges to
# 24/23 ~= 1.0435. That's a real property of the ddof = 1 formula (the
# p-1 vs N-1 denominators never fully cancel, even in the limit), not an
# artifact of this fixture -- confirmed against real pyActigraphy source
# (_interdaily_stability() calls pandas' .var(), ddof = 1 by default).
.theoretical_is_perfect <- function(n_days) {
  N <- 24 * n_days
  p <- 24
  (1 / (p - 1)) / (n_days / (N - 1))
}

test_that("compute_npcra() recovers exact IS/L5/M10/RA for a perfectly repeating profile", {
  d      <- make_npcra_fixture()
  result <- compute_npcra(d)

  # d has n_days = 7; default trim_to_d1 = TRUE removes one day -> 6 days
  # remain. IS is no longer exactly 1 under ddof = 1 -- see
  # .theoretical_is_perfect()'s comment above.
  expect_equal(result$IS, round(.theoretical_is_perfect(6L), 4))
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

  # Off-wrist exclusion removes day 1 entirely (days 2-7 remain); the
  # default D+1 trim then removes day 2 too, on top of that -> 5 days left.
  expect_equal(result$IS, round(.theoretical_is_perfect(5L), 4))
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
  # 7-day fixture; default trim_to_d1 removes 1 day -> 6 remain, dividing
  # evenly into three clean 2-day windows (avoids a ragged final window).
  d      <- make_npcra_fixture(n_days = 7L)
  result <- compute_npcra(d, window_days = 2)

  expect_true("window_start" %in% names(result))
  expect_equal(nrow(result), 3L)
  # Each 2-day window should still show the perfect-repeat IS value.
  expect_equal(result$IS, rep(round(.theoretical_is_perfect(2L), 4), 3L))
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
  # should not change L5/M10/RA at all. IS DOES change slightly -- it's a
  # function of n_days under the ddof = 1 formula (6 vs 7 days here), not a
  # fixed constant -- so check both against the theoretical value rather
  # than against each other.
  expect_equal(result_default$IS, round(.theoretical_is_perfect(6L), 4))
  expect_equal(result_notrim$IS,  round(.theoretical_is_perfect(7L), 4))
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

test_that("compute_npcra() zero-fills a genuinely missing hourly bin, matching an explicit zero", {
  # Regression coverage for the missing-hour fix: IS/IV used to build the
  # hourly series via tapply(), which just DROPS a bin with no data at all
  # (shrinking N/p) instead of inserting a real zero -- matching Python's
  # `s_1h = s.resample('1h').mean().fillna(0)` (Cell 16). If a whole hour's
  # epochs are off-wrist-excluded (genuinely missing from the data reaching
  # .npcra_core()), the result should be IDENTICAL to that same hour having
  # an explicit activity = 0 with no off-wrist exclusion at all -- both
  # should zero-fill to the same X.
  d_missing <- make_npcra_fixture(n_days = 3L)
  hour      <- as.integer(format(d_missing$datetime, "%H", tz = "UTC"))
  day       <- as.Date(d_missing$datetime, tz = "UTC")
  target    <- day == as.Date("2024-01-02") & hour == 15L   # an "active" hour

  d_missing$state <- 0L
  d_missing$state[target] <- 4L   # off-wrist -> genuinely absent from binning

  d_explicit_zero <- make_npcra_fixture(n_days = 3L)
  d_explicit_zero$activity[target] <- 0.0   # present, just recorded as zero

  result_missing <- compute_npcra(d_missing,       trim_to_d1 = FALSE)
  result_zero     <- compute_npcra(d_explicit_zero, trim_to_d1 = FALSE)

  expect_equal(result_missing$IS, result_zero$IS)
  expect_equal(result_missing$IV, result_zero$IV)
})
