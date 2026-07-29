# Tests for lids.R: lids_transform(), fit_lids(), detect_lids_bouts(),
# compute_lids(), study_lids_metrics() -- LIDS ultradian-rhythm pipeline
# (Winnebeck et al. 2018; Hammad et al. 2026).

# ---- Shared fixtures ---------------------------------------------------------

# A clean synthetic sloped-cosine LIDS profile: period 60 min, 1-min epochs,
# 5 h long, small noise -- should recover period/amplitude/offset closely.
make_lids_fixture <- function(period = 60, n_min = 300, amplitude = 15,
                               offset = 85, slope = -0.05, phase = 0, noise_sd = 1,
                               seed = 1) {
  set.seed(seed)
  t <- seq(0, n_min - 1)
  amplitude * cos(2 * pi * t / period + phase) + offset + slope * t +
    stats::rnorm(length(t), sd = noise_sd)
}

# A raw activity signal with a clear ~5h nighttime quiet bout (18:00-23:00)
# embedded in noisy daytime activity, at 1-min epochs over 2 days.
make_activity_fixture <- function() {
  t0 <- as.POSIXct("2024-01-01 12:00:00", tz = "UTC")
  n  <- 2L * 24L * 60L
  dt <- t0 + 60L * seq.int(0L, n - 1L)
  hour <- as.numeric(format(dt, "%H", tz = "UTC")) +
    as.numeric(format(dt, "%M", tz = "UTC")) / 60

  set.seed(42)
  is_night <- (hour >= 18 & hour < 23)
  activity <- ifelse(is_night, stats::runif(n, 0, 2), stats::runif(n, 20, 60))

  list(datetime = dt, activity = activity)
}

# ---- lids_transform() --------------------------------------------------------

test_that("lids_transform() maps zero activity to 100 and increases toward it", {
  out <- lids_transform(c(0, 0, 0, 0, 0, 0), method = "mva", win_min = 3)
  expect_true(all(out > 90))
})

test_that("lids_transform() decreases as activity increases", {
  low  <- lids_transform(rep(0, 20), method = "mva", win_min = 5)
  high <- lids_transform(rep(50, 20), method = "mva", win_min = 5)
  expect_true(mean(low) > mean(high))
})

test_that("lids_transform() interpolates internal NAs", {
  x <- c(0, 0, NA, NA, 0, 0, 0, 0, 0, 0)
  out <- lids_transform(x, method = "mva", win_min = 3)
  expect_false(anyNA(out))
})

test_that("lids_transform() errors with fewer than 2 epochs", {
  expect_error(lids_transform(1), "at least 2 epochs")
})

test_that("lids_transform() gaussian and mva methods both return finite, same-length output", {
  x <- make_lids_fixture()
  g <- lids_transform(x, method = "gaussian")
  m <- lids_transform(x, method = "mva")
  expect_length(g, length(x))
  expect_length(m, length(x))
  expect_true(all(is.finite(g)))
  expect_true(all(is.finite(m)))
})

# ---- fit_lids() ---------------------------------------------------------------

test_that("fit_lids() recovers the period of a clean synthetic sloped cosine", {
  lids <- make_lids_fixture(period = 60, amplitude = 15, offset = 85, slope = -0.05, noise_sd = 0.5)
  fit  <- fit_lids(lids, period_range = c(30, 120), period_step = 1)

  expect_equal(fit$period_min, 60, tolerance = 3)
  expect_gt(fit$pearson_r, 0.9)
  expect_lt(fit$p_value, 0.05)
  expect_equal(fit$amplitude, 15, tolerance = 2)
  expect_equal(fit$slope_per_60min, -0.05 * 60, tolerance = 3)
})

test_that("fit_lids() phase_rad is near zero when the cosine peaks at bout start", {
  lids <- make_lids_fixture(period = 60, phase = 0, noise_sd = 0.1)
  fit  <- fit_lids(lids, period_range = c(30, 120), period_step = 1)
  expect_equal(fit$phase_rad, 0, tolerance = 0.3)
})

test_that("fit_lids() errors on short or NA-containing input", {
  expect_error(fit_lids(c(1, 2, 3)), "at least 4 epochs")
  expect_error(fit_lids(c(1, 2, NA, 4, 5)), "missing values")
})

# ---- detect_lids_bouts() -----------------------------------------------------

test_that("detect_lids_bouts() finds the embedded nighttime quiet bout", {
  fx <- make_activity_fixture()
  bouts <- detect_lids_bouts(
    fx$datetime, fx$activity,
    duration_range = c(3, 6), main_window = c("17:00", "09:00")
  )

  expect_gt(nrow(bouts), 0L)
  expect_true(all(bouts$duration_h >= 3 & bouts$duration_h <= 6))
})

test_that("detect_lids_bouts() one_per_night keeps only the longest bout per night", {
  fx <- make_activity_fixture()
  bouts <- detect_lids_bouts(
    fx$datetime, fx$activity,
    duration_range = c(1, 6), main_window = c("17:00", "09:00"), one_per_night = TRUE
  )
  nights <- as.Date(bouts$bout_start) + ifelse(
    as.numeric(format(bouts$bout_start, "%H")) < 12, 0L, 1L
  )
  expect_equal(anyDuplicated(nights), 0L)
})

test_that("detect_lids_bouts() returns an empty tibble when nothing qualifies", {
  fx <- make_activity_fixture()
  bouts <- detect_lids_bouts(fx$datetime, fx$activity, duration_range = c(100, 200))
  expect_equal(nrow(bouts), 0L)
  expect_named(bouts, c("bout_id", "bout_start", "bout_end", "duration_h", "n_epochs"))
})

test_that("detect_lids_bouts() errors on mismatched lengths", {
  expect_error(detect_lids_bouts(Sys.time() + 0:5, 1:3), "same length")
})

test_that(".fuse_bouts() merges bouts separated by a short gap and preserves longer gaps", {
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  bouts <- list(
    list(start = t0, end = t0 + 3600, n_epochs = 60L),
    list(start = t0 + 3600 + 300, end = t0 + 7200, n_epochs = 60L),   # 5-min gap
    list(start = t0 + 20000, end = t0 + 23600, n_epochs = 60L)         # far away
  )
  fused <- .fuse_bouts(bouts, max_gap_s = 900)
  expect_length(fused, 2L)
  expect_equal(fused[[1]]$start, t0)
  expect_equal(fused[[1]]$end, t0 + 7200)
})

# ---- compute_lids() ----------------------------------------------------------

test_that("compute_lids() with bout_source = 'state' extracts and fits sleep-state bouts", {
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  n  <- 300L
  dt <- t0 + 60 * seq.int(0L, n - 1L)
  lids_like <- make_lids_fixture(n_min = n, period = 60, noise_sd = 0.5)
  # invert LIDS back to an activity-like signal: high LIDS -> low activity
  activity <- pmax(0, 100 / lids_like - 1)

  d <- tibble::tibble(datetime = dt, activity = activity, state = 1L)

  result <- compute_lids(d, bout_source = "state", duration_range = c(1, 24),
                          period_range = c(30, 120), period_step = 2)

  expect_equal(nrow(result), 1L)
  expect_equal(result$period_min, 60, tolerance = 5)
  expect_true(is.na(result$participant_id))
})

test_that("compute_lids() bout_source = 'auto' picks 'state' when present, 'roenneberg' otherwise", {
  fx <- make_activity_fixture()
  d_state <- tibble::tibble(datetime = fx$datetime, activity = fx$activity, state = 0L)
  d_raw   <- tibble::tibble(datetime = fx$datetime, activity = fx$activity)

  # state-only: no state == 1 runs -> empty result, no error
  r1 <- compute_lids(d_state, duration_range = c(1, 6))
  expect_equal(nrow(r1), 0L)

  # no state column -> falls through to roenneberg detection
  r2 <- compute_lids(d_raw, duration_range = c(3, 6),
                      bout_args = list(main_window = c("17:00", "09:00")))
  expect_true(is.data.frame(r2))
})

test_that("compute_lids() errors when bout_source = 'state' is forced without a state column", {
  fx <- make_activity_fixture()
  d  <- tibble::tibble(datetime = fx$datetime, activity = fx$activity)
  expect_error(compute_lids(d, bout_source = "state"), "state")
})

test_that("compute_lids() picks up participant_id from a zeitr_result", {
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  n  <- 200L
  dt <- t0 + 60 * seq.int(0L, n - 1L)
  d  <- tibble::tibble(datetime = dt, activity = rep(0, n), state = 1L)
  result_obj <- structure(list(data = d, subject_id = "P099"), class = "zeitr_result")

  out <- compute_lids(result_obj, duration_range = c(1, 24))
  expect_equal(unique(out$participant_id), "P099")
})

test_that("compute_lids() errors on missing required columns or unsupported input", {
  expect_error(compute_lids(list(a = 1)), "must be a")
  expect_error(compute_lids(tibble::tibble(datetime = Sys.time())), "Missing required column")
})

# ---- study_lids_metrics() ----------------------------------------------------

test_that("study_lids_metrics() summarises quality-filtered bouts per participant", {
  make_result <- function(pid, n_bouts = 4L) {
    all_dt  <- list(); all_act <- list(); all_state <- list()
    for (i in seq_len(n_bouts)) {
      n <- 300L
      lids_like <- make_lids_fixture(n_min = n, period = 60, noise_sd = 0.3, seed = i)
      act <- pmax(0, 100 / lids_like - 1)
      t0  <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC") + (i - 1L) * 86400
      all_dt[[i]]    <- t0 + 60 * seq.int(0L, n - 1L)
      all_act[[i]]   <- act
      all_state[[i]] <- rep(1L, n)
    }
    d <- tibble::tibble(
      datetime = do.call(c, all_dt), activity = unlist(all_act), state = unlist(all_state)
    )
    structure(list(data = d, subject_id = pid), class = "zeitr_result")
  }

  results <- list(P1 = make_result("P1"), P2 = make_result("P2"))
  summary <- study_lids_metrics(results, min_bouts = 2, duration_range = c(1, 24))

  expect_equal(nrow(summary), 2L)
  expect_true(all(summary$n_bouts_passing >= 0L))
  expect_true("period_min_median" %in% names(summary))
})

test_that("study_lids_metrics() skips non-zeitr_result entries and errors on empty input", {
  expect_error(study_lids_metrics(list()), "non-empty list")
  expect_warning(
    expect_warning(
      out <- study_lids_metrics(list(not_a_result = data.frame(x = 1))),
      "No valid results"
    ),
    "not a"
  )
  expect_equal(nrow(out), 0L)
})
