# Tests for free_days / holidays classification infrastructure (v0.1.3)
#
# Covers:
#   1. .parse_free_days()  -- input validation and day-to-integer mapping
#   2. .is_free_day()      -- locale-independent weekday detection and all
#                             three holiday input forms (Date, YYYY-MM-DD, DD-MM)
#   3. No-holidays warning -- zeitR.no_holidays_warn option
#   4. S3 dispatch         -- zeitr_result methods forward free_days / holidays
#
# Calendar anchors used throughout (all verified, no locale dependency):
#   2024-01-08  Monday   (ISO 1)
#   2024-01-09  Tuesday  (ISO 2)
#   2024-01-10  Wednesday(ISO 3)
#   2024-01-11  Thursday (ISO 4)
#   2024-01-12  Friday   (ISO 5)
#   2024-01-13  Saturday (ISO 6)
#   2024-01-14  Sunday   (ISO 7)

# ---- Shared fixture ---------------------------------------------------------
# Minimal nights tibble: two Mon-Sun weeks, bed at 23:00, get-up at 06:00.
# get_up dates: 2024-01-09 (Tue) through 2024-01-22 (Mon).
# Default free-day schedule (Sat + Sun) -> n_fd = 4, n_wd = 10.

make_nights <- function(start_date = "2024-01-08", n_weeks = 2L) {
  dates <- seq(as.Date(start_date), by = "day", length.out = 7L * n_weeks)
  n     <- length(dates)
  bts   <- as.POSIXct(paste0(dates, " 23:00:00"), tz = "UTC")
  gts   <- bts + 7L * 3600L
  tibble::tibble(
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
}

# ---- 1. .parse_free_days() --------------------------------------------------

test_that(".parse_free_days: NULL returns Saturday + Sunday (6, 7)", {
  expect_equal(zeitR:::.parse_free_days(NULL), c(6L, 7L))
})

test_that(".parse_free_days: empty vector returns Saturday + Sunday", {
  expect_equal(zeitR:::.parse_free_days(character(0)), c(6L, 7L))
})

test_that(".parse_free_days: English day names -> correct ISO integers", {
  expect_equal(zeitR:::.parse_free_days(c("Saturday", "Sunday")), c(6L, 7L))
  expect_equal(zeitR:::.parse_free_days(c("Friday",   "Saturday")), c(5L, 6L))
  expect_equal(zeitR:::.parse_free_days("Monday"), 1L)
  expect_equal(zeitR:::.parse_free_days("Wednesday"), 3L)
})

test_that(".parse_free_days: names are case-insensitive", {
  expect_equal(zeitR:::.parse_free_days(c("SATURDAY", "sunday")), c(6L, 7L))
  expect_equal(zeitR:::.parse_free_days("FRIDAY"), 5L)
})

test_that(".parse_free_days: ISO integers pass through unchanged", {
  expect_equal(zeitR:::.parse_free_days(c(6L, 7L)), c(6L, 7L))
  expect_equal(zeitR:::.parse_free_days(c(5L, 6L)), c(5L, 6L))
  expect_equal(zeitR:::.parse_free_days(1L), 1L)
})

test_that(".parse_free_days: integer out of range errors", {
  expect_error(zeitR:::.parse_free_days(0L))
  expect_error(zeitR:::.parse_free_days(8L))
  expect_error(zeitR:::.parse_free_days(c(6L, 9L)))
})

test_that(".parse_free_days: unrecognised day name errors", {
  expect_error(zeitR:::.parse_free_days("Samstag"))
  expect_error(zeitR:::.parse_free_days("sabado"))
})

# ---- 2. .is_free_day() -- locale-independent weekday detection --------------

test_that(".is_free_day: Saturday and Sunday are free on default schedule", {
  fw  <- zeitR:::.parse_free_days(c("Saturday", "Sunday"))
  expect_false(zeitR:::.is_free_day(as.Date("2024-01-08"), NULL, fw)) # Monday
  expect_false(zeitR:::.is_free_day(as.Date("2024-01-12"), NULL, fw)) # Friday
  expect_true( zeitR:::.is_free_day(as.Date("2024-01-13"), NULL, fw)) # Saturday
  expect_true( zeitR:::.is_free_day(as.Date("2024-01-14"), NULL, fw)) # Sunday
})

test_that(".is_free_day: vectorised over a full week", {
  week <- seq(as.Date("2024-01-08"), by = "day", length.out = 7L)
  fw   <- zeitR:::.parse_free_days(c("Saturday", "Sunday"))
  # Mon Tue Wed Thu Fri Sat Sun
  expect_equal(
    zeitR:::.is_free_day(week, NULL, fw),
    c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)
  )
})

test_that(".is_free_day: custom schedule (Friday + Saturday)", {
  fw <- zeitR:::.parse_free_days(c("Friday", "Saturday"))
  expect_false(zeitR:::.is_free_day(as.Date("2024-01-14"), NULL, fw)) # Sunday
  expect_false(zeitR:::.is_free_day(as.Date("2024-01-08"), NULL, fw)) # Monday
  expect_true( zeitR:::.is_free_day(as.Date("2024-01-12"), NULL, fw)) # Friday
  expect_true( zeitR:::.is_free_day(as.Date("2024-01-13"), NULL, fw)) # Saturday
})

test_that(".is_free_day: custom schedule vectorised over week", {
  week <- seq(as.Date("2024-01-08"), by = "day", length.out = 7L)
  fw   <- zeitR:::.parse_free_days(c("Friday", "Saturday"))
  # Mon Tue Wed Thu Fri Sat Sun
  expect_equal(
    zeitR:::.is_free_day(week, NULL, fw),
    c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE)
  )
})

# ---- 3. Holiday input forms -------------------------------------------------

test_that(".is_free_day: Date object holiday overrides weekday", {
  fw  <- zeitR:::.parse_free_days(c("Saturday", "Sunday"))
  mon <- as.Date("2024-01-08") # Monday -- normally a workday
  expect_true(zeitR:::.is_free_day(mon, as.Date("2024-01-08"), fw))
})

test_that(".is_free_day: YYYY-MM-DD string holiday overrides weekday", {
  fw  <- zeitR:::.parse_free_days(c("Saturday", "Sunday"))
  mon <- as.Date("2024-01-08")
  expect_true(zeitR:::.is_free_day(mon, "2024-01-08", fw))
})

test_that(".is_free_day: YYYY-MM-DD holiday does not match a different year", {
  fw       <- zeitR:::.parse_free_days(c("Saturday", "Sunday"))
  mon_2024 <- as.Date("2024-01-08") # Monday
  mon_2025 <- as.Date("2025-01-13") # Monday
  expect_true( zeitR:::.is_free_day(mon_2024, "2024-01-08", fw))
  expect_false(zeitR:::.is_free_day(mon_2025, "2024-01-08", fw))
})

test_that(".is_free_day: DD-MM recurring holiday matches any year", {
  fw       <- zeitR:::.parse_free_days(c("Saturday", "Sunday"))
  # January 8th is a Monday in 2024 -- always a workday without a holiday override
  mon_2024 <- as.Date("2024-01-08")
  # January 13th is a Monday in 2025
  mon_2025 <- as.Date("2025-01-13")
  expect_true(zeitR:::.is_free_day(mon_2024, "08-01", fw))
  expect_true(zeitR:::.is_free_day(mon_2025, "13-01", fw))
})

test_that(".is_free_day: DD-MM holiday does not match a different day", {
  fw  <- zeitR:::.parse_free_days(c("Saturday", "Sunday"))
  mon <- as.Date("2024-01-08") # Monday January 8
  expect_false(zeitR:::.is_free_day(mon, "09-01", fw)) # January 9th
})

test_that(".is_free_day: mixed holiday forms can be combined in one vector", {
  fw    <- zeitR:::.parse_free_days(c("Saturday", "Sunday"))
  dates <- as.Date(c("2024-01-08", "2024-01-09", "2024-01-10"))
  # Monday  = YYYY-MM-DD string; Tuesday = YYYY-MM-DD string; Wednesday = DD-MM string.
  # NOTE: do NOT mix Date objects and character strings in the same c() call --
  # lubridate's c.Date method coerces everything to Date and silently NAs the
  # "DD-MM" form. Use an all-character vector; YYYY-MM-DD strings work identically
  # to Date objects once as.character() is applied inside .is_free_day().
  holidays <- c("2024-01-08", "2024-01-09", "10-01")
  expect_equal(
    zeitR:::.is_free_day(dates, holidays, fw),
    c(TRUE, TRUE, TRUE)
  )
})

# ---- 4. No-holidays warning -------------------------------------------------

test_that("compute_sleep_metrics warns when holidays = NULL", {
  nd <- make_nights()
  expect_warning(
    compute_sleep_metrics(nd, tz = "UTC", holidays = NULL),
    regexp = "holidays"
  )
})

test_that("compute_cpd_metrics warns when holidays = NULL", {
  nd <- make_nights()
  expect_warning(
    compute_cpd_metrics(nd, tz = "UTC", holidays = NULL),
    regexp = "holidays"
  )
})

test_that("no-holidays warning is suppressible with zeitR.no_holidays_warn = FALSE", {
  nd  <- make_nights()
  old <- options(zeitR.no_holidays_warn = FALSE)
  on.exit(options(old), add = TRUE)
  expect_no_warning(
    compute_sleep_metrics(nd, tz = "UTC", holidays = NULL)
  )
})

# ---- 5. compute_sleep_metrics: free_days affects n_fd / n_wd counts ---------

test_that("compute_sleep_metrics: default free_days gives 4 free nights in 2-week fixture", {
  nd  <- make_nights()
  old <- options(zeitR.no_holidays_warn = FALSE)
  on.exit(options(old), add = TRUE)
  sm <- compute_sleep_metrics(nd, tz = "UTC", holidays = NULL)
  # get_up dates land on Sat (Jan 13, 20) and Sun (Jan 14, 21) -> 4 free nights
  expect_equal(sm$n_fd, 4)
  expect_equal(sm$n_wd, 10)
  expect_equal(sm$n_overall, 14)
})

test_that("compute_sleep_metrics: Friday + Saturday schedule shifts counts", {
  nd  <- make_nights()
  old <- options(zeitR.no_holidays_warn = FALSE)
  on.exit(options(old), add = TRUE)
  sm <- compute_sleep_metrics(nd, tz = "UTC", holidays = NULL,
                              free_days = c("Friday", "Saturday"))
  # get_up on Fri (Jan 12, 19) and Sat (Jan 13, 20) -> 4 free nights, 10 workdays
  expect_equal(sm$n_fd, 4)
  expect_equal(sm$n_wd, 10)
})

test_that("compute_sleep_metrics: DD-MM holiday adds a workday to free days", {
  nd  <- make_nights()
  old <- options(zeitR.no_holidays_warn = FALSE)
  on.exit(options(old), add = TRUE)
  # get_up on 2024-01-09 (Tuesday) is normally a workday; "09-01" makes it a holiday
  sm_no_hol <- compute_sleep_metrics(nd, tz = "UTC", holidays = NULL)
  sm_hol    <- compute_sleep_metrics(nd, tz = "UTC", holidays = "09-01")
  expect_equal(sm_hol$n_fd, sm_no_hol$n_fd + 1)
  expect_equal(sm_hol$n_wd, sm_no_hol$n_wd - 1)
})

# ---- 6. S3 dispatch from zeitr_result forwards free_days and holidays -------

test_that("compute_sleep_metrics.zeitr_result forwards free_days", {
  nd     <- make_nights()
  result <- structure(
    list(nights = nd, holidays = NULL, free_days = c("Friday", "Saturday")),
    class = "zeitr_result"
  )
  old <- options(zeitR.no_holidays_warn = FALSE)
  on.exit(options(old), add = TRUE)

  sm_dispatch <- compute_sleep_metrics(result, tz = "UTC")
  sm_explicit <- compute_sleep_metrics(nd, tz = "UTC",
                                       free_days = c("Friday", "Saturday"),
                                       holidays  = NULL)
  expect_equal(sm_dispatch$n_fd, sm_explicit$n_fd)
  expect_equal(sm_dispatch$n_wd, sm_explicit$n_wd)
})

test_that("compute_sleep_metrics.zeitr_result forwards holidays", {
  nd     <- make_nights()
  result <- structure(
    list(nights = nd, holidays = "09-01", free_days = c("Saturday", "Sunday")),
    class = "zeitr_result"
  )
  old <- options(zeitR.no_holidays_warn = FALSE)
  on.exit(options(old), add = TRUE)

  sm_dispatch <- compute_sleep_metrics(result, tz = "UTC")
  sm_explicit <- compute_sleep_metrics(nd, tz = "UTC",
                                       free_days = c("Saturday", "Sunday"),
                                       holidays  = "09-01")
  expect_equal(sm_dispatch$n_fd, sm_explicit$n_fd)
})

test_that("compute_cpd_metrics.zeitr_result forwards free_days", {
  nd     <- make_nights()
  result <- structure(
    list(nights = nd, holidays = NULL, free_days = c("Friday", "Saturday")),
    class = "zeitr_result"
  )
  old <- options(zeitR.no_holidays_warn = FALSE)
  on.exit(options(old), add = TRUE)

  cpd_dispatch <- compute_cpd_metrics(result, tz = "UTC")
  cpd_explicit <- compute_cpd_metrics(nd, tz = "UTC",
                                      free_days = c("Friday", "Saturday"),
                                      holidays  = NULL)
  expect_equal(cpd_dispatch$n_free_days, cpd_explicit$n_free_days)
  expect_equal(cpd_dispatch$n_workdays,  cpd_explicit$n_workdays)
})

test_that("compute_cpd_metrics.zeitr_result forwards holidays", {
  nd     <- make_nights()
  result <- structure(
    list(nights = nd, holidays = "09-01", free_days = c("Saturday", "Sunday")),
    class = "zeitr_result"
  )
  old <- options(zeitR.no_holidays_warn = FALSE)
  on.exit(options(old), add = TRUE)

  cpd_dispatch <- compute_cpd_metrics(result, tz = "UTC")
  cpd_explicit <- compute_cpd_metrics(nd, tz = "UTC",
                                      free_days = c("Saturday", "Sunday"),
                                      holidays  = "09-01")
  expect_equal(cpd_dispatch$n_free_days, cpd_explicit$n_free_days)
})
