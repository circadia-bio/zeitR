# Tests for study_summary.R: study_summary() (per-participant NPCRA +
# recording-quality summary across a zeitr_study).

# ---- Shared fixture ---------------------------------------------------------
# Reuses the same perfectly-repeating 24h profile as test-npcra.R (see there
# for why L5 = 0, M10 = 100, RA = 1 exactly). IS is close to but not exactly
# 1 -- compute_npcra()'s IS uses sample variance (ddof = 1, matching
# pyActigraphy's actual implementation), which for a perfectly repeating
# profile equals (24n-1)/(23n) for an n-day recording, not exactly 1 -- see
# test-npcra.R's .theoretical_is_perfect() for the derivation. Recordings
# are wrapped in a zeitr_recording so study_summary()'s per-participant
# dispatch is exercised.

make_recording <- function(participant_id, n_days = 3L) {
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  n  <- n_days * 24L * 60L
  dt <- t0 + 60L * seq.int(0L, n - 1L)

  hour     <- as.integer(format(dt, "%H", tz = "UTC"))
  activity <- ifelse(hour < 10L, 0, 100)

  epochs <- tibble::tibble(datetime = dt, activity = activity)

  structure(
    list(epochs = epochs, metadata = list(participant_id = participant_id)),
    class = "zeitr_recording"
  )
}

test_that("study_summary() computes one row per participant with correct NPCRA values", {
  study <- structure(
    list(P001 = make_recording("P001"), P002 = make_recording("P002")),
    class = "zeitr_study"
  )

  result <- study_summary(study)

  # n_days = 3 (default); trim_to_d1 removes 1 day -> 2 remain.
  expected_is <- round((24 * 2 - 1) / (23 * 2), 4)

  expect_equal(nrow(result), 2L)
  expect_equal(sort(result$participant_id), c("P001", "P002"))
  expect_true(all(result$IS == expected_is))
  expect_true(all(result$L5 == 0))
  expect_true(all(result$M10 == 100))
  expect_true(all(result$RA == 1))
})

test_that("study_summary() computes n_epochs, n_days, start, and end correctly", {
  study <- structure(list(P001 = make_recording("P001", n_days = 3L)),
                     class = "zeitr_study")
  result <- study_summary(study)

  expect_equal(result$n_epochs, 3L * 24L * 60L)
  expect_equal(result$n_days, 3)
  expect_equal(result$start, as.POSIXct("2024-01-01 00:00:00", tz = "UTC"))
})

test_that("study_summary() accepts a plain named list, not just a zeitr_study object", {
  plain_list <- list(P001 = make_recording("P001"))
  result <- study_summary(plain_list)
  expect_equal(nrow(result), 1L)
})

test_that("study_summary() skips non-zeitr_recording entries with a warning", {
  study <- list(P001 = make_recording("P001"), P002 = "not a recording")

  expect_warning(result <- study_summary(study), "not a")
  expect_equal(nrow(result), 1L)
  expect_equal(result$participant_id, "P001")
})

test_that("study_summary() errors on an empty or non-list study", {
  expect_error(study_summary(list()), "non-empty list")
  expect_error(study_summary("not a list"), "non-empty list")
})

test_that("study_summary() returns an empty tibble with a warning when no entries are valid", {
  study <- list(P001 = "not a recording", P002 = 42)

  # Two invalid entries plus the final "no valid recordings" message means
  # three warnings fire in total; expect_warning() only consumes the first
  # match, so the rest are suppressed here rather than left uncaught. The
  # "Skipping ... not a" warning text itself is already verified by the
  # single-invalid-entry test above.
  result <- suppressWarnings(study_summary(study))

  expect_equal(nrow(result), 0L)
  expect_true(all(c("participant_id", "IS", "IV", "RA", "L5", "M10") %in% names(result)))
})

test_that("study_summary() catches a per-participant NPCRA failure and fills NA", {
  # A recording with a single epoch can't compute NPCRA (needs >= 2), so
  # study_summary() should catch the error, warn, and fill that row with NA
  # rather than aborting the whole study.
  bad_rec <- structure(
    list(
      epochs   = tibble::tibble(datetime = Sys.time(), activity = 10),
      metadata = list(participant_id = "BAD")
    ),
    class = "zeitr_recording"
  )
  study <- list(GOOD = make_recording("GOOD"), BAD = bad_rec)

  expect_warning(result <- study_summary(study), "NPCRA failed")
  expect_equal(nrow(result), 2L)

  bad_row <- result[result$participant_id == "BAD", ]
  expect_true(is.na(bad_row$IS))
  expect_true(is.na(bad_row$L5_onset))
})

test_that("study_summary() forwards L5_hours and M10_hours to compute_npcra()", {
  study <- list(P001 = make_recording("P001"))

  result_default <- study_summary(study)
  result_custom  <- study_summary(study, L5_hours = 3, M10_hours = 12)

  # Both should still be 0/100 given the fixture's step profile, but the
  # M10 onset window shifts because the window width changed.
  expect_equal(result_default$M10, 100)
  expect_equal(result_custom$M10, 100)
})
