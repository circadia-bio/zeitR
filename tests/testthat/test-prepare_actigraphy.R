test_that("prepare_actigraphy matches Python reference output", {
  path <- system.file("extdata", "input1.txt", package = "zeitR")
  skip_if_not(nzchar(path), "input1.txt not available")

  rec  <- read_acttrust(path, tz = "UTC")
  prep <- prepare_actigraphy(rec)

  # ── Temperature clamping ─────────────────────────────────────────────────────
  expect_equal(range(prep$int_temp, na.rm = TRUE), c(0, 42))
  expect_equal(range(prep$ext_temp, na.rm = TRUE), c(15.31, 40.56))

  # ── Scaled temperature columns ───────────────────────────────────────────────
  expect_equal(range(prep$int_temp_, na.rm = TRUE), c(0, 1))
  expect_equal(round(range(prep$ext_temp_, na.rm = TRUE), 7),
               c(0.3645238, 0.9657143))

  # ── Zero int_temp epochs (off-wrist candidates) ──────────────────────────────
  expect_equal(sum(prep$int_temp == 0), 188L)

  # ── State columns initialised to zero ────────────────────────────────────────
  expect_true(all(prep$state    == 0))
  expect_true(all(prep$offwrist == 0))
  expect_true(all(prep$sleep    == 0))

  # ── First 3 rows match Python reference ──────────────────────────────────────
  expect_equal(round(prep$int_temp_[1:3], 6), c(0.577143, 0.580000, 0.581429))
  expect_equal(round(prep$ext_temp_[1:3], 6), c(0.570000, 0.572857, 0.570000))
})
