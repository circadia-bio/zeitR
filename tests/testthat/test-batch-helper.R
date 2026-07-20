# Tests for .run_pipeline_over_files(), the shared implementation behind
# run_pipeline_batch() and run_pipeline_native_batch() (parallel option, v0.1.4).
#
# .run_pipeline_over_files() does not care what run_fn actually does with a
# file -- it only needs a path in and a result (or an error) out -- so these
# tests use a small fake run_fn instead of real ActTrust fixtures.

# ---- Shared fixture ---------------------------------------------------------
# Dummy paths; content is irrelevant since .fake_run_fn() never reads them.
# A path containing "bad" simulates a file that fails to process.

.fake_run_fn <- function(path, ...) {
  if (grepl("bad", basename(path), fixed = TRUE)) {
    stop("simulated pipeline failure")
  }
  list(file = path)
}

test_that("sequential (parallel = FALSE) processes all files and names by stem", {
  files <- c("P001.txt", "P002.txt")

  result <- .run_pipeline_over_files(files, .fake_run_fn, parallel = FALSE)

  expect_type(result, "list")
  expect_named(result, c("P001", "P002"))
  expect_equal(result$P001$file, "P001.txt")
  expect_equal(result$P002$file, "P002.txt")
})

test_that("a failing file is skipped with a warning, not aborting the batch", {
  files <- c("P001.txt", "bad_subject.txt", "P002.txt")

  expect_warning(
    result <- .run_pipeline_over_files(files, .fake_run_fn, parallel = FALSE),
    "Failed to process"
  )

  expect_named(result, c("P001", "P002"))
  expect_length(result, 2L)
})

test_that("an all-failing batch returns an empty named list, no error", {
  files <- c("bad_one.txt", "bad_two.txt")

  # Two failing files means two warnings; expect_warning() only matches the
  # first occurrence, so wrap in suppressWarnings() here. The warning
  # behaviour itself (message content, one warning per failure) is already
  # covered by the partial-failure test above.
  result <- suppressWarnings(
    .run_pipeline_over_files(files, .fake_run_fn, parallel = FALSE)
  )

  expect_type(result, "list")
  expect_length(result, 0L)
})

test_that("parallel = TRUE without future.apply falls back to sequential with a warning", {
  skip_if(
    requireNamespace("future.apply", quietly = TRUE),
    "future.apply is installed; fallback path is not exercised here"
  )

  files <- c("P001.txt", "P002.txt")

  expect_warning(
    result <- .run_pipeline_over_files(files, .fake_run_fn, parallel = TRUE),
    "future.apply"
  )

  expect_named(result, c("P001", "P002"))
})

test_that("parallel = TRUE with future.apply installed dispatches via future_lapply", {
  skip_if_not_installed("future.apply")
  skip_if_not_installed("future")

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential())

  files <- c("P001.txt", "bad_subject.txt", "P002.txt")

  # One failing file among two good ones -> exactly one warning.
  expect_warning(
    result <- .run_pipeline_over_files(files, .fake_run_fn, parallel = TRUE),
    "Failed to process"
  )

  expect_named(result, c("P001", "P002"))
  expect_length(result, 2L)
})

test_that("run_pipeline_batch() and run_pipeline_native_batch() forward parallel = FALSE by default", {
  # Regression guard: confirms the exported wrappers still default to
  # sequential processing (backwards compatible with pre-0.1.4 behaviour).
  expect_false(formals(run_pipeline_batch)$parallel)
  expect_false(formals(run_pipeline_native_batch)$parallel)
})
