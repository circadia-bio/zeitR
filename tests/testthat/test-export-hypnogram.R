# Tests for export.R: export_hypnogram() (zeitR pipeline result -> hypnoR
# hypnogram format).

# ---- Shared fixture ---------------------------------------------------------
# 60 epochs (1-min): wake / sleep (ZCMn = 0, "Quiet sleep") /
# sleep (ZCMn > 0, "Sleep") / off-wrist, so all four stage-mapping branches
# are exercised.

make_export_fixture <- function() {
  t0 <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  n  <- 40L
  dt <- t0 + 60L * seq.int(0L, n - 1L)

  state <- c(rep(0L, 10L), rep(1L, 10L), rep(7L, 10L), rep(4L, 10L))
  zcm   <- c(rep(50, 10L), rep(0, 10L), rep(20, 10L), rep(0, 10L))

  tibble::tibble(datetime = dt, state = state, ZCMn = zcm)
}

test_that("export_hypnogram() maps states to the correct stage labels", {
  d   <- make_export_fixture()
  hyp <- export_hypnogram(d)

  expect_s3_class(hyp$stage, "ordered")
  expect_equal(levels(hyp$stage), c("W", "Sleep", "Quiet sleep"))

  # state 0 (wake) -> "W"
  expect_true(all(hyp$stage[1:10] == "W"))
  # state 1, ZCMn == 0 -> "Quiet sleep"
  expect_true(all(hyp$stage[11:20] == "Quiet sleep"))
  # state 7 (nap), ZCMn > 0 -> "Sleep"
  expect_true(all(hyp$stage[21:30] == "Sleep"))
  # state 4 (off-wrist) -> "W"
  expect_true(all(hyp$stage[31:40] == "W"))
})

test_that("export_hypnogram() maps all sleep to \"Sleep\" when ZCMn is absent", {
  d <- make_export_fixture()
  d$ZCMn <- NULL

  hyp <- export_hypnogram(d)

  # state 1 and state 7 (both "in_sleep") -> all "Sleep", no "Quiet sleep"
  expect_true(all(hyp$stage[11:30] == "Sleep"))
})

test_that("export_hypnogram() resolves subject_id: explicit > result$subject_id > omitted", {
  d      <- make_export_fixture()
  result <- structure(list(data = d, subject_id = "FROM_RESULT"), class = "zeitr_result")

  # Explicit arg takes precedence
  hyp1 <- export_hypnogram(result, subject_id = "EXPLICIT")
  expect_equal(unique(hyp1$subject_id), "EXPLICIT")

  # Falls back to result$subject_id
  hyp2 <- export_hypnogram(result)
  expect_equal(unique(hyp2$subject_id), "FROM_RESULT")

  # Bare tibble, no subject_id anywhere -> column omitted entirely
  hyp3 <- export_hypnogram(d)
  expect_false("subject_id" %in% names(hyp3))
})

test_that("export_hypnogram() accepts a zeitr_result list with $data", {
  d      <- make_export_fixture()
  result <- structure(list(data = d), class = "zeitr_result")

  hyp <- export_hypnogram(result)
  expect_equal(nrow(hyp), nrow(d))
})

test_that("export_hypnogram() errors when a zeitr_result list has no $data element", {
  bad_result <- structure(list(subject_id = "P001"), class = "zeitr_result")
  expect_error(export_hypnogram(bad_result), "must be a")
})

test_that("export_hypnogram() errors clearly when required columns are missing", {
  bad <- tibble::tibble(datetime = Sys.time())
  expect_error(export_hypnogram(bad), "missing required column")
})

test_that("export_hypnogram(drop_offwrist = TRUE) removes off-wrist epochs and re-indexes", {
  d   <- make_export_fixture()
  hyp <- export_hypnogram(d, drop_offwrist = TRUE)

  expect_equal(nrow(hyp), 30L)
  expect_equal(hyp$epoch, seq_len(30L))
})

test_that("export_hypnogram() warns when epoch_sec differs from the observed interval", {
  # Regression test for the .envir interpolation bug fixed alongside the
  # actogram error-path bugs: this warning previously failed with
  # 'object mode_gap/epoch_sec not found' instead of firing correctly.
  d <- make_export_fixture()  # 1-minute (60s) epochs

  expect_warning(
    export_hypnogram(d, epoch_sec = 30),
    "differs from"
  )

  # No warning when epoch_sec matches the observed interval.
  # (Written manually rather than via expect_no_warning() for compatibility
  # with the package's stated testthat >= 3.0.0 minimum.)
  warnings_raised <- character(0)
  withCallingHandlers(
    export_hypnogram(d, epoch_sec = 60),
    warning = function(w) {
      warnings_raised <<- c(warnings_raised, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings_raised, 0L)
})

test_that("export_hypnogram() source column defaults to \"zeitR\" and can be overridden", {
  d <- make_export_fixture()

  hyp1 <- export_hypnogram(d)
  expect_true(all(hyp1$source == "zeitR"))

  hyp2 <- export_hypnogram(d, source = "manual")
  expect_true(all(hyp2$source == "manual"))
})
