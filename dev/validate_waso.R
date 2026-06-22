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

# ── 1. .nights_df boundaries vs python_nights ────────────────────────────────
ndf <- get(".nights_df", asNamespace("zeitR"))
nd  <- ndf(out_pre[onwrist], wake_thresh = 60L)
pn  <- read.csv(file.path(ext, "python_nights.csv"), stringsAsFactors = FALSE)
pn_nap <- tolower(as.character(pn$nap)) == "true"

cat(sprintf("\nnights: R = %d   Python = %d\n", nrow(nd), nrow(pn)))
if (nrow(nd) == nrow(pn)) {
  cat(sprintf("  bt  match: %s\n", all(nd$bt  == pn$bt)))
  cat(sprintf("  gt  match: %s\n", all(nd$gt  == pn$gt)))
  cat(sprintf("  nap match: %s\n", all(nd$nap == pn_nap)))
} else {
  cat("  (night count differs - inspect nd vs pn)\n")
  print(utils::head(nd)); print(utils::head(pn[, c("bt", "gt", "nap")]))
}

# ── 2. compute_waso state vs python_output$state ─────────────────────────────
res     <- compute_waso(x, wake_thresh = 60L)
state_r <- as.integer(res$data$state)

mism <- sum(state_r != state_final)
cat(sprintf("\ncompute_waso state vs python_output$state: %d / %d mismatches  (agreement %.5f)\n",
            mism, n, mean(state_r == state_final)))
if (mism > 0L) {
  cat("mismatch breakdown (rows = python state, cols = R state):\n")
  print(table(python = state_final[state_r != state_final],
              R      = state_r[state_r != state_final]))
}

# ── 3. per-night stats vs python_nights ──────────────────────────────────────
ns <- res$nights
if (nrow(ns) == nrow(pn)) {
  cat(sprintf("\nper-night stat agreement (n = %d):\n", nrow(ns)))
  cat(sprintf("  tbt:  %s\n", all(ns$tbt  == pn$tbt)))
  cat(sprintf("  waso: %s\n", all(ns$waso == pn$waso)))
  cat(sprintf("  sol:  %s\n", all(ns$sol  == pn$sol)))
  cat(sprintf("  soi:  %s\n", all(ns$soi  == pn$soi)))
  cat(sprintf("  tst:  %s\n", all(ns$tst  == pn$tst)))
  cat(sprintf("  nw:   %s\n", all(ns$nw   == pn$nw)))
  cat(sprintf("  eff:  max abs diff = %.3e\n", max(abs(ns$eff - pn$eff))))
}

cat("\ndone.\n")
