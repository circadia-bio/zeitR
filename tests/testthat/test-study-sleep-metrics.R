# Tests for study_sleep_metrics.R: batch wrapper stacking
# compute_sleep_metrics() and compute_cpd_metrics() across a study into one
# sync-ready tibble (the sleep-timing/chronotype analogue of study_summary()).

# ---- Shared fixture ---------------------------------------------------------
# Synthetic zeitr_result objects built directly from a nights tibble (no real
# pipeline run needed -- compute_sleep_metrics()/compute_cpd_metrics() only
# need is_nap, bed_time, get_up_time, tbt, tst, sol, soi, waso, eff). Bed
# 23:00, get-up 07:00 (7h TBT, well above both min_tib_h defaults). Weekday
# membership is computed via ISO 8601 weekday number, matching the package's
# own locale-independent convention (see .parse_free_days()/.is_free_day()).

make_result <- function(subject_id, dates, free_days = c("Saturday", "Sunday"),
                        holidays = "01-01") {
  n   <- length(dates)
  bts <- as.POSIXct(paste0(dates, " 23:00:00"), tz = "UTC")
  gts <- bts + 7L * 3600L

  nights <- tibble::tibble(
    night       = seq_len(n),
    is_nap      = FALSE,
    sleep_type  = "main",
    bed_time    = bts,
    get_up_time = gts,
    tbt         = 420,
    tst         = 360,
    waso        = 30,
    sol         = 15,
    soi         = 10,
    nw          = 5L,
    eff         = 360 / 420
  )

  structure(
    list(nights = nights, subject_id = subject_id,
        holidays = holidays, free_days = free_days),
    class = "zeitr_result"
  )
}

# Two full Mon-Sun weeks: includes workdays, free days (Sat/Sun), and
# free-day-eve nights -- everything compute_cpd_metrics() needs to succeed.
two_week_dates <- seq(as.Date("2024-01-08"), by = "day", length.out = 14L)

# Weekdays only, and restricted to Monday-Thursday bed-dates specifically:
# since bed_time is 23:00 and get_up_time is bed_time + 7h, a Friday bed-date
# would roll over into a SATURDAY get-up (23:00 + 7h crosses midnight), which
# would still count as a free day via is_free_day()'s get-up-date check. Only
# Mon-Thu bed-dates guarantee every get-up date is also Tue-Fri, i.e. zero
# free-day nights at all. compute_sleep_metrics() still succeeds on this
# fixture (just NA-fills the _fd subgroup), but compute_cpd_metrics() aborts
# with "No free days found" -- a targeted per-metric partial failure.
wd_num              <- as.integer(format(two_week_dates, "%u"))
weekdays_only_dates <- two_week_dates[wd_num <= 4L]

test_that("study_sleep_metrics() computes one row per participant with both metric sets", {
  results <- list(
    P001 = make_result("P001", two_week_dates),
    P002 = make_result("P002", two_week_dates)
  )

  out <- study_sleep_metrics(results)

  expect_equal(nrow(out), 2L)
  expect_equal(sort(out$participant_id), c("P001", "P002"))

  # compute_sleep_metrics() columns present and sane
  expect_true(all(c("n_overall", "n_wd", "n_fd", "tst_h", "sleep_onset_h") %in% names(out)))
  expect_equal(out$n_overall, c(14, 14))
  expect_equal(out$n_fd, c(4, 4))    # Sat/Sun get-ups in a 2-week Mon-Sun fixture
  expect_equal(out$n_wd, c(10, 10))

  # compute_cpd_metrics() columns present and sane
  expect_true(all(c("msw_h", "msf_h", "sjl_h", "cpd_s") %in% names(out)))
  expect_true(all(is.finite(out$msw_h)))
  expect_true(all(is.finite(out$cpd_s)))
})

test_that("study_sleep_metrics() forwards a participant's holidays into the free/workday split", {
  # get-up date 2024-01-11 is a Thursday (normally a workday). Declaring it a
  # recurring holiday via "11-01" should move that one night into the
  # free-day group: 4 weekend nights (default fixture) + 1 holiday = 5.
  results <- list(P001 = make_result("P001", two_week_dates, holidays = "11-01"))

  out <- study_sleep_metrics(results)

  expect_equal(out$n_fd, 5)
  expect_equal(out$n_wd, 9)
})

test_that("study_sleep_metrics() falls back per-participant to result$holidays/free_days", {
  # P001 uses a Friday+Saturday weekend, P002 uses the default Saturday+Sunday.
  # Both are set on the result objects, not passed to study_sleep_metrics(),
  # so each should use its own schedule.
  results <- list(
    P001 = make_result("P001", two_week_dates, free_days = c("Friday", "Saturday")),
    P002 = make_result("P002", two_week_dates, free_days = c("Saturday", "Sunday"))
  )

  out <- study_sleep_metrics(results)

  # Fri+Sat get-ups (Jan 12/13, 19/20) -> 4 free nights; Sat+Sun likewise -> 4.
  expect_equal(out$n_fd[out$participant_id == "P001"], 4)
  expect_equal(out$n_fd[out$participant_id == "P002"], 4)
})

test_that("study_sleep_metrics() free_days argument overrides every participant's own default", {
  results <- list(
    P001 = make_result("P001", two_week_dates, free_days = c("Friday", "Saturday")),
    P002 = make_result("P002", two_week_dates, free_days = c("Saturday", "Sunday"))
  )

  out <- study_sleep_metrics(results, free_days = c("Saturday", "Sunday"))

  # Both forced onto Sat+Sun regardless of their own stored free_days.
  expect_equal(out$n_fd, c(4, 4))
})

test_that("study_sleep_metrics() NA-fills only the failing metric set for a participant", {
  # GOOD has a normal two-week fixture; NOFREE has weekdays only, so
  # compute_cpd_metrics() aborts (no free days) but compute_sleep_metrics()
  # still succeeds (n_fd = 0, _fd columns NA, everything else populated).
  results <- list(
    GOOD   = make_result("GOOD",   two_week_dates),
    NOFREE = make_result("NOFREE", weekdays_only_dates)
  )

  expect_warning(
    out <- study_sleep_metrics(results),
    "compute_cpd_metrics\\(\\) failed"
  )

  expect_equal(nrow(out), 2L)

  good   <- out[out$participant_id == "GOOD", ]
  nofree <- out[out$participant_id == "NOFREE", ]

  # GOOD: both metric sets populated
  expect_true(is.finite(good$msw_h))
  expect_true(is.finite(good$n_overall))

  # NOFREE: compute_sleep_metrics() succeeded (n_overall/n_wd populated,
  # n_fd = 0 since there were no free days at all)...
  expect_equal(nofree$n_overall, 8)
  expect_equal(nofree$n_fd, 0)
  expect_true(is.na(nofree$tst_h_fd))

  # ...but compute_cpd_metrics() failed entirely -> NA, not an aborted batch.
  expect_true(is.na(nofree$msw_h))
  expect_true(is.na(nofree$cpd_s))
})

test_that("study_sleep_metrics() skips non-zeitr_result entries with a warning", {
  results <- list(P001 = make_result("P001", two_week_dates), P002 = "not a result")

  expect_warning(out <- study_sleep_metrics(results), "not a")
  expect_equal(nrow(out), 1L)
  expect_equal(out$participant_id, "P001")
})

test_that("study_sleep_metrics() errors on an empty or non-list results argument", {
  expect_error(study_sleep_metrics(list()), "non-empty list")
  expect_error(study_sleep_metrics("not a list"), "non-empty list")
})

test_that("study_sleep_metrics() returns an empty, correctly-shaped tibble when nothing is valid", {
  results <- list(P001 = "not a result", P002 = 42)

  out <- suppressWarnings(study_sleep_metrics(results))

  expect_equal(nrow(out), 0L)
  expect_true(all(c("participant_id", "n_overall", "msw_h", "cpd_s") %in% names(out)))
})

test_that("study_sleep_metrics() falls back to the list name when subject_id is missing", {
  bare <- make_result("P001", two_week_dates)
  bare$subject_id <- NULL
  results <- list(FALLBACK_NAME = bare)

  out <- study_sleep_metrics(results)
  expect_equal(out$participant_id, "FALLBACK_NAME")
})
