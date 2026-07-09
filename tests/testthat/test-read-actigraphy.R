# Tests for read_actigraphy.R: read_actigraphy(), read_actigraphy_dir(), and
# the print.zeitr_recording / print.zeitr_study S3 methods.
#
# Reuses inst/extdata/input1.txt, the same real (committed, non-privacy-
# restricted) ActTrust fixture used by test-read_acttrust.R.

test_that("read_actigraphy() returns a well-formed zeitr_recording", {
  path <- system.file("extdata", "input1.txt", package = "zeitR")
  skip_if_not(nzchar(path), "input1.txt not available")

  rec <- read_actigraphy(path, tz = "UTC")

  expect_s3_class(rec, "zeitr_recording")
  expect_true(all(c("epochs", "metadata") %in% names(rec)))
  expect_equal(nrow(rec$epochs), 76196L)
  expect_true(all(c("datetime", "activity", "int_temp", "ext_temp",
                    "ZCMn", "state", "offwrist", "sleep") %in% names(rec$epochs)))

  # participant_id derived from the filename stem
  expect_equal(rec$metadata$participant_id, "input1")

  # the raw tibble's own zeitr_recording class/metadata attr should not leak
  # into the embedded $epochs tibble
  expect_false(inherits(rec$epochs, "zeitr_recording"))
  expect_null(attr(rec$epochs, "metadata"))
})

test_that("read_actigraphy() errors clearly on an unsupported device", {
  path <- system.file("extdata", "input1.txt", package = "zeitR")
  skip_if_not(nzchar(path), "input1.txt not available")

  expect_error(
    read_actigraphy(path, device = "some_other_device"),
    "Unsupported device"
  )
})

test_that("read_actigraphy() device argument is case- and whitespace-insensitive", {
  path <- system.file("extdata", "input1.txt", package = "zeitR")
  skip_if_not(nzchar(path), "input1.txt not available")

  rec <- read_actigraphy(path, device = "  ActTrust  ", tz = "UTC")
  expect_s3_class(rec, "zeitr_recording")
})

test_that("read_actigraphy_dir() reads every matching file into a zeitr_study", {
  path <- system.file("extdata", "input1.txt", package = "zeitR")
  skip_if_not(nzchar(path), "input1.txt not available")

  tmp <- tempfile("zeitr-test-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  file.copy(path, file.path(tmp, "P001.txt"))
  file.copy(path, file.path(tmp, "P002.txt"))

  study <- read_actigraphy_dir(tmp, tz = "UTC")

  expect_s3_class(study, "zeitr_study")
  expect_equal(sort(names(study)), c("P001", "P002"))
  expect_true(all(vapply(study, inherits, logical(1), "zeitr_recording")))
})

test_that("read_actigraphy_dir() returns an empty zeitr_study with a warning when nothing matches", {
  tmp <- tempfile("zeitr-test-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  expect_warning(study <- read_actigraphy_dir(tmp), "No files matching")
  expect_s3_class(study, "zeitr_study")
  expect_length(study, 0L)
})

test_that("read_actigraphy_dir() skips unparseable files with a warning and keeps the rest", {
  path <- system.file("extdata", "input1.txt", package = "zeitR")
  skip_if_not(nzchar(path), "input1.txt not available")

  tmp <- tempfile("zeitr-test-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  file.copy(path, file.path(tmp, "GOOD.txt"))
  writeLines("this is not a valid ActTrust file", file.path(tmp, "BAD.txt"))

  expect_warning(study <- read_actigraphy_dir(tmp, tz = "UTC"), "Skipping")
  expect_equal(names(study), "GOOD")
})

test_that("print.zeitr_recording() and print.zeitr_study() run without error", {
  path <- system.file("extdata", "input1.txt", package = "zeitR")
  skip_if_not(nzchar(path), "input1.txt not available")

  # print.zeitr_recording()/print.zeitr_study() build their output with
  # cli::cli_h1()/cli_bullets(), which signal R conditions rather than
  # writing to stdout via cat() -- so this needs expect_message(), not
  # expect_output().
  rec <- read_actigraphy(path, tz = "UTC")
  expect_message(print(rec), "zeitr_recording")

  study <- structure(list(input1 = rec), class = "zeitr_study")
  expect_message(print(study), "zeitr_study")
})
