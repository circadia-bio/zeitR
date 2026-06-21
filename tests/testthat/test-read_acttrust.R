test_that("read_acttrust parses input1.txt correctly", {
  path <- system.file("extdata", "input1.txt", package = "zeitR")
  skip_if_not(nzchar(path), "input1.txt not available")

  rec <- read_acttrust(path, tz = "UTC")

  # ── Dimensions ──────────────────────────────────────────────────────────────
  expect_equal(nrow(rec), 76196L)

  # ── Required columns present ─────────────────────────────────────────────────
  expect_true(all(c("datetime", "activity", "int_temp",
                    "ext_temp", "ZCMn", "state", "offwrist", "sleep") %in%
                    names(rec)))

  # ── Datetime parsing ─────────────────────────────────────────────────────────
  expect_s3_class(rec$datetime, "POSIXct")
  expect_equal(format(rec$datetime[1]),     "2021-05-27 11:10:15")
  expect_equal(format(rec$datetime[nrow(rec)]), "2021-07-19 10:56:25")

  # ── Epoch interval ───────────────────────────────────────────────────────────
  interval_s <- as.numeric(difftime(rec$datetime[2], rec$datetime[1],
                                    units = "secs"))
  expect_equal(interval_s, 60)

  # ── Activity range ───────────────────────────────────────────────────────────
  expect_equal(min(rec$activity), 0)
  expect_equal(max(rec$activity), 103688)

  # ── First three rows match Python reference ──────────────────────────────────
  expect_equal(rec$activity[1:3], c(4856, 4483, 425))
  expect_equal(round(rec$int_temp[1:3], 2), c(24.24, 24.36, 24.42))
  expect_equal(round(rec$ext_temp[1:3], 2), c(23.94, 24.06, 23.94))

  # ── State columns initialised to zero ────────────────────────────────────────
  expect_true(all(rec$state    == 0))
  expect_true(all(rec$offwrist == 0))
  expect_true(all(rec$sleep    == 0))

  # ── No NAs in ZCMn ───────────────────────────────────────────────────────────
  expect_equal(sum(is.na(rec$ZCMn)), 0L)
})
