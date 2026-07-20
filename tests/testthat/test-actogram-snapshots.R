# Snapshot tests for the three actogram plot functions, using vdiffr.
#
# Entirely skipped if vdiffr or ggplot2 are not installed (both Suggests-only,
# matching the package's existing optional-dependency pattern for mclust and
# future.apply). vdiffr snapshots are SVGs stored under
# tests/testthat/_snaps/actogram-snapshots/ and are compared on subsequent
# runs; a changed plot shows as a diff for review via vdiffr::manage_cases()
# or testthat::snapshot_review(), not a silent pass.
#
# NOTE (first run only): the baseline SVGs do not exist yet. The first
# `devtools::test()` run after adding this file will report these tests as
# skipped with "New file snapshot ... has been written", not passed. Review
# the generated SVGs under tests/testthat/_snaps/actogram-snapshots/, then
# commit them alongside this file -- from then on the tests do a real
# comparison against the committed baseline.

# ---- Shared fixture ---------------------------------------------------------
# 2 days of 1-min epochs (2880 epochs), deterministic state + activity pattern
# so the rendered plots are fully reproducible:
#   00:00-06:00  wake        ZCMn ramps with a slow sine (never exactly 0)
#   06:00-14:00  sleep       ZCMn = 0
#   14:00-22:00  wake        ZCMn ramps with a slow sine
#   22:00-24:00  off-wrist   ZCMn = 0

make_actogram_fixture <- function() {
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  n  <- 2L * 24L * 60L
  dt <- t0 + 60L * seq.int(0L, n - 1L)

  hour  <- as.integer(format(dt, "%H", tz = "UTC"))
  state <- integer(n)
  state[hour >= 0L  & hour < 6L]  <- 0L  # wake
  state[hour >= 6L  & hour < 14L] <- 1L  # sleep
  state[hour >= 14L & hour < 22L] <- 0L  # wake
  state[hour >= 22L]              <- 4L  # off-wrist

  zcm <- ifelse(state == 0L, 200 + 50 * sin(seq_len(n) / 30), 0)

  tibble::tibble(datetime = dt, state = state, ZCMn = zcm)
}

test_that("plot_actogram() renders consistently", {
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")

  d <- make_actogram_fixture()
  vdiffr::expect_doppelganger(
    "actogram-single",
    plot_actogram(d, tz = "UTC", title = "Actogram")
  )
})

test_that("plot_actogram_double() renders consistently", {
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")

  d <- make_actogram_fixture()
  vdiffr::expect_doppelganger(
    "actogram-double",
    plot_actogram_double(d, tz = "UTC", title = "Double actogram")
  )
})

test_that("plot_actogram_activity() renders consistently", {
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")

  d <- make_actogram_fixture()
  vdiffr::expect_doppelganger(
    "actogram-activity",
    plot_actogram_activity(d, tz = "UTC", title = "Activity actogram")
  )
})

test_that("plot_actogram_activity() with a custom activity_cap_quantile renders consistently", {
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")

  d <- make_actogram_fixture()
  vdiffr::expect_doppelganger(
    "actogram-activity-cap95",
    plot_actogram_activity(d, tz = "UTC", title = "Activity actogram (95th pct cap)",
                           activity_cap_quantile = 0.95)
  )
})

test_that("plot_actogram_activity() with log_scale = TRUE renders consistently", {
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")

  d <- make_actogram_fixture()
  vdiffr::expect_doppelganger(
    "actogram-activity-log",
    plot_actogram_activity(d, tz = "UTC", title = "Activity actogram (log scale)",
                           log_scale = TRUE)
  )
})

# ---- Non-visual regression checks (run regardless of vdiffr) ---------------
# These don't need vdiffr: they check the error paths and don't depend on
# comparing rendered pixels, so they run whenever ggplot2 is available.

test_that("actogram functions error clearly when required columns are missing", {
  skip_if_not_installed("ggplot2")

  bad <- tibble::tibble(datetime = Sys.time())
  expect_error(plot_actogram(bad), "missing column")
  expect_error(plot_actogram_double(bad), "missing column")
  expect_error(plot_actogram_activity(bad), "missing column")
})

test_that("plot_actogram_activity() errors clearly on a missing activity_col", {
  skip_if_not_installed("ggplot2")

  d <- make_actogram_fixture()
  expect_error(
    plot_actogram_activity(d, activity_col = "not_a_real_column"),
    "not_a_real_column"
  )
})

test_that("all three actogram functions accept a zeitr_result list, not just a bare tibble", {
  skip_if_not_installed("ggplot2")

  d      <- make_actogram_fixture()
  result <- structure(list(data = d, subject_id = "P999"), class = "zeitr_result")

  expect_s3_class(plot_actogram(result, tz = "UTC"), "ggplot")
  expect_s3_class(plot_actogram_double(result, tz = "UTC"), "ggplot")
  expect_s3_class(plot_actogram_activity(result, tz = "UTC"), "ggplot")
})

test_that("log_scale = FALSE (default) is unchanged from pre-log_scale behaviour", {
  skip_if_not_installed("ggplot2")

  d <- make_actogram_fixture()
  p_default        <- plot_actogram_activity(d, tz = "UTC")
  p_explicit_false <- plot_actogram_activity(d, tz = "UTC", log_scale = FALSE)

  expect_identical(p_default$data$act_norm, p_explicit_false$data$act_norm)
})

test_that("log_scale = TRUE changes bar heights relative to the linear scale", {
  skip_if_not_installed("ggplot2")

  d <- make_actogram_fixture()
  p_linear <- plot_actogram_activity(d, tz = "UTC")
  p_log    <- plot_actogram_activity(d, tz = "UTC", log_scale = TRUE)

  # Different normalised heights overall...
  expect_false(isTRUE(all.equal(p_linear$data$act_norm, p_log$data$act_norm)))

  # ...but zero-activity epochs (sleep, off-wrist) never become non-finite --
  # log1p(0) = 0, not -Inf.
  expect_true(all(is.finite(p_log$data$act_norm)))
})
