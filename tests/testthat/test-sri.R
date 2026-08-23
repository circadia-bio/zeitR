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

test_that("compute_sri() returns exactly -100 for a perfectly inverted alternating pattern", {
  # Whole calendar days alternate entirely sleep / entirely wake, so every
  # epoch differs from the epoch exactly 24h later at every single minute
  # (not just during a shifted sleep window, which would still agree on
  # wake outside of it).
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  n_days <- 6L
  n  <- n_days * 24L * 60L
  dt <- t0 + 60L * seq.int(0L, n - 1L)

  day   <- as.integer(format(dt, "%d", tz = "UTC")) - 1L
  state <- ifelse(day %% 2L == 0L, 1L, 0L)   # even days all-sleep, odd all-wake

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

# ---- algo = "sadeh" -----------------------------------------------------
# Regression coverage ported directly from pyActigraphy's real source
# (pyActigraphy/sleep/scoring_base.py's _sadeh(), and pyActigraphy/sleep/
# scoring/sri.py's sri()/prob_stability()). Expected values below were
# computed by running the ACTUAL pandas-based algorithm on these exact
# fixtures (not hand-derived), to avoid arithmetic mistakes on a
# moderately complex multi-term formula -- see this session's transcript
# for the verification script if it needs to be re-run.

test_that(".sadeh_score() matches pyActigraphy's real _sadeh() exactly on a synthetic fixture", {
  # A single spike (60, within the 50-100 NAT band) surrounded by zeros.
  # Verified against a direct run of pyActigraphy's actual _sadeh() on
  # this exact vector.
  activity <- c(0,0,0,0,0,10,60,10,0,0,0,0,0,0,0,0,0,0,0,0)
  expected <- c(0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0)

  expect_equal(.sadeh_score(activity), as.integer(expected))
})

test_that(".rolling_centered_mean()/.rolling_centered_count_between()/.rolling_trailing_sd() match pandas' intermediate columns exactly", {
  # Same fixture and verification method as the .sadeh_score() test above --
  # these are pandas' own mean_W5/NAT/sd_Last6 columns for this input.
  activity <- c(0,0,0,0,0,10,60,10,0,0,0,0,0,0,0,0,0,0,0,0)

  mean_W5 <- .rolling_centered_mean(activity, halfwin = 5L)
  expect_true(all(is.na(mean_W5[1:5])))
  expect_true(all(is.na(mean_W5[16:20])))
  expect_equal(mean_W5[6], 80 / 11, tolerance = 1e-8)    # window = activity[1:11]
  expect_equal(mean_W5[12], 70 / 11, tolerance = 1e-8)   # window = activity[7:17]

  NAT <- .rolling_centered_count_between(activity, halfwin = 5L, lo = 50, hi = 100)
  expect_true(all(is.na(NAT[1:5])))
  expect_equal(NAT[6], 1)     # window = activity[1:11], one value (60) in (50,100)
  expect_equal(NAT[13], 0)    # window = activity[8:18], no values in (50,100)

  sd_Last6 <- .rolling_trailing_sd(activity, win = 6L)
  expect_true(all(is.na(sd_Last6[1:5])))
  expect_equal(sd_Last6[6], stats::sd(activity[1:6]), tolerance = 1e-8)
  expect_equal(sd_Last6[7], stats::sd(activity[2:7]), tolerance = 1e-8)
})

test_that(".sri_pyactigraphy() matches pyActigraphy's real sri()/prob_stability() exactly", {
  # 3 days, hourly epochs, a perfectly repeating sleep(1)/wake(0) pattern
  # (hours 0-7 sleep, 8-23 wake) with a SINGLE mismatch injected at day 2,
  # hour 3. Verified against a direct run of pyActigraphy's actual sri()
  # on this exact series: SRI = 91.66666666666669.
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  dt <- t0 + 3600 * 0:71
  vals <- rep(c(rep(1, 8), rep(0, 16)), 3)
  vals[24 + 3 + 1] <- 0   # day 2 (0-indexed), hour 3 -> flip 1 to 0 (R 1-indexed)

  result <- .sri_pyactigraphy(dt, vals, tz = "UTC", threshold = NULL)
  expect_equal(result, 91.66666666666669, tolerance = 1e-8)
})

test_that("compute_sri(algo = 'sadeh') is wired correctly end to end", {
  # Reuses the perfectly-regular activity-free SRI fixture's TIMING, but
  # needs an activity column instead of state -- construct a simple
  # repeating activity pattern (active by day, quiet by night) over
  # several days, and confirm the function runs, returns a finite SRI in
  # [-100, 100], and that its n_pairs is NA (pyActigraphy's two-step
  # average has no single pooled pair count -- see compute_sri()'s docs).
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  n  <- 4L * 24L * 60L
  dt <- t0 + 60L * seq.int(0L, n - 1L)
  hour <- as.integer(format(dt, "%H", tz = "UTC"))
  activity <- ifelse(hour < 8L, 0, 50)   # quiet at night, active by day

  d      <- tibble::tibble(datetime = dt, activity = activity)
  result <- compute_sri(d, algo = "sadeh")

  expect_true(is.finite(result$sri))
  expect_gte(result$sri, -100)
  expect_lte(result$sri, 100)
  expect_true(is.na(result$n_pairs))
})

test_that("compute_sri()'s default algo is 'vallim', unchanged from before the algo argument existed", {
  d <- make_sri_fixture()
  expect_equal(compute_sri(d), compute_sri(d, algo = "vallim"))
})

test_that("compute_sri(algo = 'sadeh') errors on missing activity column", {
  d <- make_sri_fixture()   # has state, not activity
  expect_error(compute_sri(d, algo = "sadeh"), "Missing required column")
})

# ---- algo = "ck" ----------------------------------------------------------
# Regression coverage ported directly from pyActigraphy's real source
# (pyActigraphy/sleep/scoring_base.py's _cole_kripke()/CK(), and
# pyActigraphy/sleep/scoring/utils.py's consecutive_values()/
# rescore_if_preceded()/rescore_if_surrounded()/rescore()). Expected
# values below were computed by running the ACTUAL algorithm on these
# exact fixtures, not hand-derived.

test_that(".consecutive_values() matches pyActigraphy's real consecutive_values() exactly, both targets", {
  x <- c(0,0,1,1,1,0,0,1,1,1,1,1,0,1,0,0)

  runs1 <- .consecutive_values(x, target = 1, min_length = 2)
  expect_equal(unname(runs1[, "start"]), c(3L, 8L))
  expect_equal(unname(runs1[, "end"]),   c(5L, 12L))

  runs0 <- .consecutive_values(x, target = 0, min_length = 2)
  expect_equal(unname(runs0[, "start"]), c(1L, 6L, 15L))
  expect_equal(unname(runs0[, "end"]),   c(2L, 7L, 16L))
})

test_that(".consecutive_values() returns 0 rows when nothing qualifies", {
  x <- c(1, 0, 1, 0, 1)
  runs <- .consecutive_values(x, target = 1, min_length = 2)
  expect_equal(nrow(runs), 0L)
})

test_that(".consecutive_values() treats NA as not-a-match, matching numpy's np.equal(NaN, target) == False", {
  # Required for Roenneberg's seed-finding, where the categorized series
  # can genuinely contain NA (trend undefined at the edges) -- R's own
  # `NA == target` would otherwise propagate NA and corrupt the whole
  # computation, unlike numpy's np.equal(NaN, x) which is simply False.
  x <- c(1, 1, NA, 1, 1, 1, 0, 0)
  # The NA breaks what would otherwise be one run of five 1's into two
  # shorter runs (length 2 and length 3) -- neither should be silently
  # merged across the NA, and NA itself should never count as a match.
  runs <- .consecutive_values(x, target = 1, min_length = 2)
  expect_equal(unname(runs[, "start"]), c(1L, 4L))
  expect_equal(unname(runs[, "end"]),   c(2L, 6L))
})

test_that(".rescore_if_preceded() matches pyActigraphy exactly", {
  # 4 wake epochs then 1 isolated sleep epoch -> rescored to wake.
  x1 <- c(0,0,0,0,1,0,0,0,0)
  expect_equal(.rescore_if_preceded(x1, n_periods = 1L, n_previous = 4L),
               c(1,1,1,1,0,1,1,1,1))

  # Only 3 preceding wake epochs (not enough) -> NOT rescored.
  x2 <- c(0,0,0,1,0,0,0,0)
  expect_equal(.rescore_if_preceded(x2, n_periods = 1L, n_previous = 4L),
               c(1,1,1,1,1,1,1,1))
})

test_that(".rescore_if_surrounded() matches pyActigraphy exactly", {
  # wake(4) sleep(2) wake(4): gap (2) <= n_periods (3) -> gap rescored to wake.
  x3 <- c(0,0,0,0,1,1,0,0,0,0)
  expect_equal(.rescore_if_surrounded(x3, n_periods = 3L, n_surround = 4L),
               c(1,1,1,1,0,0,1,1,1,1))

  # wake(4) sleep(4) wake(4): gap (4) > n_periods (3) -> NOT rescored.
  x4 <- c(0,0,0,0,1,1,1,1,0,0,0,0)
  expect_equal(.rescore_if_surrounded(x4, n_periods = 3L, n_surround = 4L),
               rep(1, 12))
})

test_that(".ck_native_score() matches pyActigraphy's real CK()/_cole_kripke()/rescore() exactly on a synthetic fixture", {
  # Mostly quiet activity with a few noisy bursts, long enough to exercise
  # both the raw scoring and Webster's rescoring rules. Verified against a
  # direct run of pyActigraphy's actual _cole_kripke() + rescore() on this
  # exact vector.
  activity <- rep(0, 200)
  activity[21:25]   <- c(80, 5, 90, 3, 70)   # short noisy burst
  activity[61:65]   <- 100                    # a real wake period
  activity[101:103] <- c(60, 2, 55)           # another short blip
  activity[151:170] <- 100                    # long wake stretch

  scored <- .ck_native_score(activity)

  # Only these 5 positions (1-indexed) should differ from the raw
  # (pre-rescoring) scoring, per the real Python run.
  changed <- c(5L, 172L, 173L, 174L, 175L)
  raw     <- .ck_native_score(activity, rescoring = FALSE)

  expect_equal(which(scored != raw), changed)
  expect_true(all(scored[changed] == 0L))   # all rescored TO wake, matching Webster's rules
})

test_that("compute_sri(algo = 'ck') is wired correctly end to end", {
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  n  <- 4L * 24L * 60L
  dt <- t0 + 60L * seq.int(0L, n - 1L)
  hour <- as.integer(format(dt, "%H", tz = "UTC"))
  activity <- ifelse(hour < 8L, 0, 50)

  d      <- tibble::tibble(datetime = dt, activity = activity)
  result <- compute_sri(d, algo = "ck")

  expect_true(is.finite(result$sri))
  expect_gte(result$sri, -100)
  expect_lte(result$sri, 100)
  expect_true(is.na(result$n_pairs))
})

test_that("compute_sri(algo = 'ck') errors on missing activity column", {
  d <- make_sri_fixture()
  expect_error(compute_sri(d, algo = "ck"), "Missing required column")
})

# ---- algo = "scripps" -------------------------------------------------
# Regression coverage ported directly from pyActigraphy's real source
# (pyActigraphy/sleep/scoring_base.py's _scripps()/Scripps()). Expected
# values computed by running the ACTUAL algorithm on this exact fixture.

test_that(".scripps_score() matches pyActigraphy's real _scripps() exactly on a synthetic fixture", {
  activity <- rep(0, 60)
  activity[21:25] <- c(80, 5, 90, 3, 70)
  activity[41:45] <- 100

  expected <- c(rep(0,10), rep(1,11), rep(0,5), rep(1,14), rep(0,20))
  expect_equal(.scripps_score(activity), as.integer(expected))
})

test_that("compute_sri(algo = 'scripps') is wired correctly end to end", {
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  n  <- 4L * 24L * 60L
  dt <- t0 + 60L * seq.int(0L, n - 1L)
  hour <- as.integer(format(dt, "%H", tz = "UTC"))
  activity <- ifelse(hour < 8L, 0, 50)

  d      <- tibble::tibble(datetime = dt, activity = activity)
  result <- compute_sri(d, algo = "scripps")

  expect_true(is.finite(result$sri))
  expect_gte(result$sri, -100)
  expect_lte(result$sri, 100)
  expect_true(is.na(result$n_pairs))
})

test_that("compute_sri(algo = 'scripps') errors on missing activity column", {
  d <- make_sri_fixture()
  expect_error(compute_sri(d, algo = "scripps"), "Missing required column")
})

# ---- algo = "roenneberg" -------------------------------------------------
# Regression coverage ported directly from pyActigraphy's real source
# (pyActigraphy/sleep/scoring/roenneberg.py's roenneberg() and its
# sub-functions, pyActigraphy/sleep/scoring/utils.py's pearsonr()/
# correlation_series()/find_first_peak_idx()). By far the most involved
# of the four scoring algorithms -- verified in stages: the trend's exact
# centered-window edge behaviour, then the full pipeline end to end on a
# deterministic (not random -- R and numpy don't share a RNG) synthetic
# fixture with one clear 7h quiet period.

test_that(".roenneberg_trend() matches pandas' exact centered-window/min_periods edge behaviour", {
  # data = 0..19, win_size = 8 (even), min_win_size = 4. Verified against a
  # direct run of pandas' actual .rolling(8, center=True, min_periods=4,
  # closed='right').mean() on this exact input -- confirms pandas puts the
  # extra element on the LEFT for even windows (4 before, self, 3 after in
  # the interior), and that edge positions use whatever's available once
  # the count reaches min_periods.
  x <- 0:19
  trend <- .roenneberg_trend(x, win_size = 8L, min_win_size = 4L)
  expected <- c(1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5,
                9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.0, 16.5, 17.0)
  expect_equal(trend, expected, tolerance = 1e-8)
})

test_that(".pearsonr() matches pyActigraphy's real pearsonr() on simple cases", {
  expect_equal(.pearsonr(c(1,2,3,4), c(1,2,3,4)), 1.0, tolerance = 1e-8)
  expect_equal(.pearsonr(c(1,2,3,4), c(4,3,2,1)), -1.0, tolerance = 1e-8)
  expect_equal(.pearsonr(c(1,2,3,4), c(2,4,6,8)), 1.0, tolerance = 1e-8)
})

test_that(".find_highest_peak_idx() matches pyActigraphy's real find_highest_peak_idx() exactly", {
  # Position 3 (1-indexed) is a LOCAL peak over the next 2 values, but
  # position 6 (value 9) is the GLOBAL highest among all qualifying
  # peaks -- find_highest_peak_idx() must return 6, not 3 (which is what
  # the WRONG find_first_peak_idx() would have returned; this exact
  # example is verified against a direct run of the real function).
  x <- c(1, 2, 5, 1, 1, 9, 3, 3, 3)
  expect_equal(.find_highest_peak_idx(x, n_succ = 2L), 6L)

  # No peak long enough -> NA.
  x2 <- c(1, 2, 3, 4, 5)
  expect_true(is.na(.find_highest_peak_idx(x2, n_succ = 3L)))
})

test_that(".roenneberg_score() matches pyActigraphy's real roenneberg() exactly on a deterministic fixture", {
  # Fully deterministic (no randomness -- R and numpy don't share a RNG,
  # so a random fixture couldn't be reproduced identically in both
  # languages): 2 days at 1-min epochs, a smooth baseline wiggle
  # (40 + 5*sin(t*0.05)) with ONE clear 7h quiet period (1 + 0.3*sin(t*0.2))
  # from minute 1380 to 1799 (0-indexed) -- day 1 23:00 to day 2 06:00.
  # Verified against a direct run of pyActigraphy's actual roenneberg()
  # on this exact array: exactly one detected bout, onset at minute 1380
  # (R position 1381), offset at minute 1800 (R position 1801), value
  # counts {sleep: 420, wake: 2460}.
  n <- 2L * 24L * 60L
  t <- 0:(n - 1L)
  activity <- 40 + 5 * sin(t * 0.05)
  activity[1381:1800] <- 1 + 0.3 * sin(t[1381:1800] * 0.2)

  scoring <- .roenneberg_score(
    activity,
    trend_period_min     = 1440L,
    min_trend_period_min = 720L,
    threshold             = 0.15,
    min_seed_period_min  = 30L,
    max_test_period_min  = 720L,
    r_consec_below_min   = 30L
  )

  expect_equal(sum(scoring == 1L), 420L)
  expect_equal(sum(scoring == 0L), 2460L)
  expect_equal(scoring[1381], 1L)
  expect_equal(scoring[1380], 0L)   # epoch just before onset
  expect_equal(scoring[1801], 0L)   # epoch just after offset (offset itself is 1)
  expect_equal(scoring[1800], 1L)   # last sleep epoch
  expect_true(all(scoring[1381:1800] == 1L))
  expect_true(all(scoring[1:1380] == 0L))
  expect_true(all(scoring[1801:n] == 0L))
})

test_that("compute_sri(algo = 'roenneberg') is wired correctly end to end", {
  n <- 3L * 24L * 60L
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  dt <- t0 + 60L * seq.int(0L, n - 1L)
  t  <- 0:(n - 1L)
  activity <- 40 + 5 * sin(t * 0.05)
  activity[1381:1800] <- 1 + 0.3 * sin(t[1381:1800] * 0.2)

  d      <- tibble::tibble(datetime = dt, activity = activity)
  result <- compute_sri(d, algo = "roenneberg")

  expect_true(is.finite(result$sri))
  expect_gte(result$sri, -100)
  expect_lte(result$sri, 100)
  expect_true(is.na(result$n_pairs))
})

test_that("compute_sri(algo = 'roenneberg') errors on missing activity column", {
  d <- make_sri_fixture()
  expect_error(compute_sri(d, algo = "roenneberg"), "Missing required column")
})
