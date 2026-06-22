# Fixture-based regression tests ported from dev/validate_waso.R.
#
# Validates compute_waso() in isolation against the Python reference by
# reconstructing the exact pre-WASO state and comparing per-night statistics
# to python_nights.csv. This isolates the WASO / Cole-Kripke / boundary layer
# from the full pipeline.
#
# Strategy (mirrors validate_waso.R):
#   - off-wrist epochs: preserved from python_output.csv (state == 4)
#   - on-wrist epochs:  1 - cspd_refined_output.csv (1 = sleep period, 0 = wake)
#   - nap detection is NOT applied here, so nap nights are absent from R's
#     .nights_df; the test uses a matched-night join (by bt index) to compare
#     only the nights that R can see.

waso_fixture <- function(file) {
  path <- system.file("extdata", file, package = "zeitR")
  if (!nzchar(path)) testthat::skip(paste0(file, " not available in inst/extdata"))
  path
}

# Build the pre-WASO tibble shared across tests.
build_waso_input <- function() {
  py          <- read.csv(waso_fixture("python_output.csv"),
                          stringsAsFactors = FALSE)
  py$datetime <- as.POSIXct(py$datetime, tz = "UTC")
  n           <- nrow(py)
  offwrist    <- py$state == 4L

  refined     <- as.numeric(scan(waso_fixture("cspd_refined_output.csv"),
                                 quiet = TRUE))   # 1=wake, 0=sleep (on-wrist)

  out_pre            <- integer(n)
  out_pre[offwrist]  <- 4L
  out_pre[!offwrist] <- ifelse(refined == 0L, 1L, 0L)

  tibble::tibble(
    datetime = py$datetime,
    ZCMn     = as.double(py$ZCMn),
    state    = out_pre
  )
}

# ---------------------------------------------------------------------------
# .nights_df boundaries
# ---------------------------------------------------------------------------

test_that(".nights_df bt/gt boundaries match Python reference (non-nap nights)", {
  ndf     <- getFromNamespace(".nights_df", "zeitR")
  inp     <- build_waso_input()
  pn      <- read.csv(waso_fixture("python_nights.csv"),
                      stringsAsFactors = FALSE)
  onwrist <- inp$state != 4L

  nd <- ndf(as.integer(inp$state[onwrist]), wake_thresh = 60L)

  # Join on bt. Some matched nights will have different gt because Python's
  # python_nights.csv was produced after nap detection, which can extend a
  # night's gt boundary when a nap is appended immediately after. We assert
  # only that the majority of matched nights agree on gt, and that at least
  # one night matches exactly — the full pipeline parity test covers the
  # post-nap boundaries end-to-end.
  matched <- merge(
    data.frame(bt = nd$bt, gt_r = nd$gt),
    data.frame(bt = as.integer(pn$bt), gt_py = as.integer(pn$gt)),
    by = "bt"
  )

  expect_gt(nrow(matched), 0L)
  gt_agree <- matched$gt_r == matched$gt_py
  expect_gt(mean(gt_agree), 0.9,
            label = sprintf("%d / %d matched nights agree on gt",
                            sum(gt_agree), nrow(matched)))
})

# ---------------------------------------------------------------------------
# Per-night statistics on matched nights
# ---------------------------------------------------------------------------

test_that("compute_waso per-night stats match Python reference on matched nights", {
  inp  <- build_waso_input()
  pn   <- read.csv(waso_fixture("python_nights.csv"),
                   stringsAsFactors = FALSE)
  ndf  <- getFromNamespace(".nights_df", "zeitR")
  onwrist <- inp$state != 4L
  nd   <- ndf(as.integer(inp$state[onwrist]), wake_thresh = 60L)

  res  <- compute_waso(inp, wake_thresh = 60L)
  rn   <- res$nights
  rn$bt <- nd$bt
  rn$gt <- nd$gt

  # Join on both bt and gt: nights where Python's boundary shifted due to nap
  # detection will not match gt and are correctly excluded here. The full
  # pipeline parity test covers those nights end-to-end.
  m <- merge(
    rn[, c("bt", "gt", "tbt", "waso", "sol", "soi", "tst", "nw", "eff")],
    pn[, c("bt", "gt", "tbt", "waso", "sol", "soi", "tst", "nw", "eff")],
    by = c("bt", "gt"), suffixes = c("_r", "_py")
  )

  expect_gt(nrow(m), 40L,
            label = sprintf("%d boundary-matched nights", nrow(m)))

  expect_equal(m$tbt_r,  m$tbt_py,  tolerance = 1e-9, label = "TBT")
  expect_equal(m$waso_r, m$waso_py, tolerance = 1e-9, label = "WASO")
  expect_equal(m$sol_r,  m$sol_py,  tolerance = 1e-9, label = "SOL")
  expect_equal(m$soi_r,  m$soi_py,  tolerance = 1e-9, label = "SOI")
  expect_equal(m$tst_r,  m$tst_py,  tolerance = 1e-9, label = "TST")
  expect_equal(m$nw_r,   m$nw_py,   tolerance = 1e-9, label = "NW")
  expect_equal(m$eff_r,  m$eff_py,  tolerance = 1e-9, label = "efficiency")
})

# ---------------------------------------------------------------------------
# Epoch-level state (matched nights only)
# ---------------------------------------------------------------------------

test_that("compute_waso epoch-level state agrees with Python on non-nap nights", {
  inp         <- build_waso_input()
  pn          <- read.csv(waso_fixture("python_nights.csv"),
                          stringsAsFactors = FALSE)
  py          <- read.csv(waso_fixture("python_output.csv"),
                          stringsAsFactors = FALSE)
  ndf         <- getFromNamespace(".nights_df", "zeitR")
  onwrist     <- inp$state != 4L

  nd      <- ndf(as.integer(inp$state[onwrist]), wake_thresh = 60L)
  res     <- compute_waso(inp, wake_thresh = 60L)
  r_state  <- as.integer(res$data$state)
  py_state <- as.integer(py$state)

  # Only test nights where both bt and gt match Python exactly (i.e. nap
  # detection did not shift the boundary). The 4 residual mismatches in the
  # original validate_waso.R all fall inside the boundary-shifted nights.
  matched_bt_gt <- merge(
    data.frame(bt = nd$bt, gt = nd$gt),
    data.frame(bt = as.integer(pn$bt), gt = as.integer(pn$gt)),
    by = c("bt", "gt")
  )

  if (nrow(matched_bt_gt) == 0L) testthat::skip("no boundary-matched nights")

  total_mismatch <- 0L
  total_epochs   <- 0L

  for (i in seq_len(nrow(matched_bt_gt))) {
    bt <- matched_bt_gt$bt[i]
    gt <- matched_bt_gt$gt[i]
    full_positions <- which(onwrist)[seq(bt + 1L, gt)]
    total_mismatch <- total_mismatch + sum(r_state[full_positions] != py_state[full_positions])
    total_epochs   <- total_epochs   + length(full_positions)
  }

  expect_equal(total_mismatch, 0L,
               label          = sprintf("%d within-night epoch mismatch(es) across %d epochs",
                                        total_mismatch, total_epochs),
               expected.label = "0 mismatches")
})
