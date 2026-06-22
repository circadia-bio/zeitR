# dev/validate_waso.R
#
# Validate compute_waso (per-epoch Cole-Kripke WASO scoring) against the Python
# reference. Run from the package root:
#
#   source("dev/validate_waso.R")
#
# Strategy: reconstruct the EXACT pre-WASO state that the Python detect_waso
# received, feed it through compute_waso, and compare the epoch-level state to
# python_output$state (post-WASO ground truth). This isolates compute_waso from
# the upstream chain (which is separately validated bit-exact).
#
#   pre-WASO state:  off-wrist (state == 4) preserved from python_output;
#                    on-wrist  = 1 - cspd_refined_output (1 = sleep period,
#                    0 = wake). No naps in this recording (python_nights nap is
#                    False throughout), so no state == 7 to reconstruct.
#
# Checks: (1) .nights_df boundaries vs python_nights bt/gt/nap;
#         (2) per-night stats vs python_nights;
#         (3) full epoch-level state vs python_output$state.

suppressMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    stop("install devtools, or adapt this script to source the R/ files directly")
  }
})

ext <- "inst/extdata"
rd1 <- function(f) as.numeric(scan(file.path(ext, f), quiet = TRUE))

py          <- read.csv(file.path(ext, "python_output.csv"), stringsAsFactors = FALSE)
py$datetime <- as.POSIXct(py$datetime, tz = "UTC")
n           <- nrow(py)
state_final <- as.integer(py$state)

cat(sprintf("rows: %d   unique final states: {%s}\n",
            n, paste(sort(unique(state_final)), collapse = ", ")))

offwrist <- state_final == 4L
onwrist  <- !offwrist
cat(sprintf("on-wrist: %d   off-wrist: %d\n", sum(onwrist), sum(offwrist)))

# ── Reconstruct the pre-WASO state ───────────────────────────────────────────
refined <- rd1("cspd_refined_output.csv")          # 1 = wake, 0 = sleep
stopifnot(length(refined) == sum(onwrist))
out_pre            <- integer(n)
out_pre[offwrist]  <- 4L
out_pre[onwrist]   <- ifelse(refined == 0, 1L, 0L) # 1 = sleep period, 0 = wake

x <- tibble::tibble(
  datetime = py$datetime,
  ZCMn     = as.double(py$ZCMn),
  state    = out_pre
)

# ── nights_df + compute_waso, then a MATCHED-night comparison ────────────────
# NOTE: this reconstruction has no naps (detect_naps not run), so R will be
# missing every nap-detected night that the Python pipeline adds as state 1.
# The matched nights isolate the WASO / Cole-Kripke / boundary layer.
ndf <- get(".nights_df", asNamespace("zeitR"))
nd  <- ndf(out_pre[onwrist], wake_thresh = 60L)
pn  <- read.csv(file.path(ext, "python_nights.csv"), stringsAsFactors = FALSE)

res     <- compute_waso(x, wake_thresh = 60L)
state_r <- as.integer(res$data$state)

cat(sprintf("\nnights: R = %d   Python = %d\n", nrow(nd), nrow(pn)))

# res$nights rows correspond 1:1 (same order) with nd rows.
rstats <- data.frame(
  bt  = nd$bt,           gt  = nd$gt,
  tbt = res$nights$tbt,  waso = res$nights$waso, sol = res$nights$sol,
  soi = res$nights$soi,  tst  = res$nights$tst,  nw  = res$nights$nw,
  eff = res$nights$eff
)
m <- merge(rstats, pn, by = "bt", suffixes = c("_r", "_py"))
cat(sprintf("matched nights (by bt): %d  (R-only %d, Python-only %d)\n",
            nrow(m), sum(!(rstats$bt %in% pn$bt)), sum(!(pn$bt %in% rstats$bt))))
if (nrow(m) > 0L) {
  cat(sprintf("  matched-night agreement -> gt:%s tbt:%s waso:%s sol:%s soi:%s tst:%s nw:%s  eff_maxdiff:%.2e\n",
              all(m$gt_r == m$gt_py), all(m$tbt_r == m$tbt_py),
              all(m$waso_r == m$waso_py), all(m$sol_r == m$sol_py),
              all(m$soi_r == m$soi_py), all(m$tst_r == m$tst_py),
              all(m$nw_r == m$nw_py), max(abs(m$eff_r - m$eff_py))))
}

py_only <- pn[!(pn$bt %in% rstats$bt), c("bt", "gt", "tbt", "nap")]
if (nrow(py_only) > 0L) {
  cat("Python-only nights (missing from R = nap-detected periods):\n")
  print(py_only)
}
r_only <- rstats[!(rstats$bt %in% pn$bt), c("bt", "gt", "tbt")]
if (nrow(r_only) > 0L) {
  cat("R-only nights:\n"); print(r_only)
}

# ── full epoch-level state vs python_output$state ────────────────────────────
mism <- sum(state_r != state_final)
cat(sprintf("\ncompute_waso state vs python_output$state: %d / %d mismatches  (agreement %.5f)\n",
            mism, n, mean(state_r == state_final)))
if (mism > 0L) {
  cat("mismatch breakdown (rows = python state, cols = R state):\n")
  print(table(python = state_final[state_r != state_final],
              R      = state_r[state_r != state_final]))
}

cat("\nNote: residual is expected until detect_naps_crespo faithfully ports\n")
cat("nap_wrapper (full CSPD nap-mode, naps encoded as state 1).\n")
cat("\ndone.\n")
