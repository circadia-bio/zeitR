test_that("detect_offwrist_bimodal matches Python reference output exactly", {
  path    <- system.file("extdata", "input1.txt",        package = "zeitR")
  py_path <- system.file("extdata", "python_output.csv", package = "zeitR")
  skip_if_not(nzchar(path),    "input1.txt not available")
  skip_if_not(nzchar(py_path), "python_output.csv not available")

  rec  <- read_acttrust(path, tz = "UTC")
  prep <- prepare_actigraphy(rec)
  prep <- detect_offwrist_bimodal(prep)

  r_ow  <- as.integer(prep$state == 4)
  py_ow <- as.integer(read.csv(py_path)$state == 4)

  # ── Epoch counts ────────────────────────────────────────────────────────────
  expect_equal(sum(r_ow),  sum(py_ow))

  # ── Exact epoch-level agreement ─────────────────────────────────────────────
  expect_equal(r_ow, py_ow)
})
