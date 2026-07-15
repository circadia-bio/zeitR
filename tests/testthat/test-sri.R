# Tests for sri.R: compute_sri() (Sleep Regularity Index, Fix 30 -- STATE +
# off-wrist based, ported from the Python reference pipeline's SRI_vallim).

# ---- Shared fixture ---------------------------------------------------------
# A perfectly regular 7-day sleep/wake pattern at 1-minute epochs:
#   00:00-08:00  state = 1 (sleep, 8 h)
#   08:00-24:00  state = 0 (wake, 16 h)
# Identical every day -> every 24h-apart pair matches -> SRI = 100 exactly.

make_sri_fixture <- function(n_days = 7L) {
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  n  <- n_days * 24L * 60L
  dt <- t0 + 60L * seq.int(0L, n - 1L)

  hour  <- as.integer(format(dt, "%H", tz = "UTC"))
  state <- ifelse(hour < 8L, 1L, 0L)

  tibble::tibble(datetime = dt, state = state)
}

test_that("compute_sri() returns exactly 100 for a perfectly regular pattern", {
  d      <- make_sri_fixture()
  result <- compute_sri(d)

  expect_equal(result$sri, 100)
  expect_gt(result$n_pairs, 0L)
})

test_that("compute_sri() returns close to -100 for a perfectly inverted alternating pattern", {
  # Sleep/wake flips every day: day 1 sleeps 00:00-08:00, day 2 sleeps
  # 12:00-20:00, day 3 back to 00:00-08:00, etc. Every 24h-apart pair should
  # mismatch.
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  n_days <- 6L
  n  <- n_days * 24L * 60L
  dt <- t0 + 60L * seq.int(0L, n - 1L)

  day   <- as.integer(format(dt, "%d", tz = "UTC")) - 1L
  hour  <- as.integer(format(dt, "%H", tz = "UTC"))
  even_day_sleep <- hour < 8L
  odd_day_sleep  <- hour >= 12L & hour < 20L
  state <- ifelse(day %% 2L == 0L, as.integer(even_day_sleep), as.integer(odd_day_sleep))

  d      <- tibble::tibble(datetime = dt, state = state)
  result <- compute_sri(d)

  expect_equal(result$sri, -100)
})

test_that("compute_sri() excludes off-wrist gaps longer than max_gap_min", {
  d <- make_sri_fixture()
  # Knock out day 4 entirely as off-wrist (> 30 min) -- those epochs should
  # be excluded from the comparison, not counted as mismatches, so SRI
  # should remain exactly 100 (all remaining pairs still match perfectly).
  day4 <- as.Date(d$datetime, tz = "UTC") == as.Date("2024-01-04")
  d$state[day4] <- 4L

  result <- compute_sri(d)
  expect_equal(result$sri, 100)
  full_result <- compute_sri(make_sri_fixture())
  expect_lt(result$n_pairs, full_result$n_pairs)
})

test_that("compute_sri() interpolates off-wrist gaps of max_gap_min or less", {
  d <- make_sri_fixture()
  # A 10-minute off-wrist gap in the middle of a sleep block should be
  # interpolated (ffill) back to sleep, not excluded -- SRI stays 100.
  gap_idx <- which(format(d$datetime, "%H:%M") == "03:00")[1]
  d$state[gap_idx:(gap_idx + 9L)] <- 4L

  result <- compute_sri(d)
  expect_equal(result$sri, 100)
})

test_that("compute_sri() accepts a zeitr_result and uses its subject_id", {
  d      <- make_sri_fixture()
  result_obj <- structure(
    list(data = d, subject_id = "P042"),
    class = "zeitr_result"
  )

  result <- compute_sri(result_obj)
  expect_equal(result$participant_id, "P042")
  expect_equal(result$sri, 100)
})

test_that("compute_sri() participant_id is NA for a bare data frame", {
  d      <- make_sri_fixture()
  result <- compute_sri(d)
  expect_true(is.na(result$participant_id))
})

test_that("compute_sri() warns and returns NA for a recording shorter than 24h", {
  d <- make_sri_fixture()[1:(12L * 60L), ]   # 12 h only

  expect_warning(result <- compute_sri(d), "less than 24")
  expect_true(is.na(result$sri))
  expect_equal(result$n_pairs, 0L)
})

test_that("compute_sri() errors on missing required columns", {
  bad <- tibble::tibble(datetime = Sys.time())
  expect_error(compute_sri(bad), "Missing required column")
})

test_that("compute_sri() errors when x is neither a recognised object nor a data frame", {
  expect_error(compute_sri(list(a = 1)), "must be a")
  expect_error(compute_sri("not a recording"), "must be a")
})

test_that("compute_sri() errors with fewer than 2 epochs", {
  d <- make_sri_fixture()[1, ]
  expect_error(compute_sri(d), "at least 2 epochs")
})

test_that(".interpolate_short_gaps() forward-fills short gaps and leaves long gaps as NA", {
  x <- c(1, 1, NA, NA, 1, 0, NA, NA, NA, NA, 0)
  out <- .interpolate_short_gaps(x, max_gap_epochs = 2L)

  expect_equal(out[3:4], c(1, 1))       # short gap ffilled from x[2]
  expect_true(all(is.na(out[7:10])))    # long gap (4 > 2) left untouched
})

test_that(".interpolate_short_gaps() back-fills a short gap at the start of the vector", {
  x <- c(NA, NA, 1, 0, 0)
  out <- .interpolate_short_gaps(x, max_gap_epochs = 2L)
  expect_equal(out[1:2], c(1, 1))
})

test_that(".interpolate_short_gaps() leaves a gap spanning the whole vector as NA", {
  x <- c(NA, NA, NA)
  out <- .interpolate_short_gaps(x, max_gap_epochs = 5L)
  expect_true(all(is.na(out)))
})
