# dev/validate_naps.R
#
# Validate the faithful nap_wrapper port (detect_naps_crespo) end-to-end against
# the Python reference. Run from the package root:
#
#   source("dev/validate_naps.R")
#
# Strategy: reconstruct the EXACT pre-NAP state the Python detect_naps received
# (off-wrist = 4; on-wrist = 1 - cspd_refined_output, i.e. the detect_sleep
# _periods output with NO naps), run detect_naps_crespo (which adds nap-detected
# sleep periods as state 1, matching nap_wrapper), then compute_waso, and compare
# the epoch-level state to python_output$state (post-WASO ground truth).
#
# detect_sleep_crespo is separately validated bit-exact (1 - refined matches its
# on-wrist state 0/70654), so out_pre is its output and this isolates the nap
# port + WASO layer.
#
# Target: 0 / 76196 mismatches, and 55 nights matching python_nights.

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

offwrist <- state_final == 4L
onwrist  <- !offwrist
cat(sprintf("rows: %d   on-wrist: %d   off-wrist: %d\n",
            n, sum(onwrist), sum(offwrist)))

# ── Reconstruct the pre-NAP state (detect_sleep_crespo output, no naps) ──────
refined <- rd1("cspd_refined_output.csv")          # 1 = wake, 0 = sleep
stopifnot(length(refined) == sum(onwrist))
out_pre           <- integer(n)
out_pre[offwrist] <- 4L
out_pre[onwrist]  <- ifelse(refined == 0, 1L, 0L)  # 1 = sleep period, 0 = wake

cat(sprintf("pre-nap state: sleep(1)=%d  wake(0)=%d  offwrist(4)=%d\n",
            sum(out_pre == 1L), sum(out_pre == 0L), sum(out_pre == 4L)))

x <- tibble::tibble(
  datetime = py$datetime,
  activity = as.double(py$activity),
  ZCMn     = as.double(py$ZCMn),
  state    = out_pre
)

# ── 1. detect_naps_crespo: add nap-detected periods as state 1 ───────────────
xn        <- detect_naps_crespo(x)
state_nap <- as.integer(xn$state)
n_added   <- sum(state_nap == 1L) - sum(out_pre == 1L)
cat(sprintf("\ndetect_naps_crespo: sleep(1)=%d  wake(0)=%d  (nap epochs added: %d)\n",
            sum(state_nap == 1L), sum(state_nap == 0L), n_added))
cat(sprintf("unique post-nap states: {%s}\n",
            paste(sort(unique(state_nap)), collapse = ", ")))

# ── 2. nights after naps vs python_nights ────────────────────────────────────
ndf <- get(".nights_df", asNamespace("zeitR"))
nd  <- ndf(state_nap[onwrist], wake_thresh = 60L)
pn  <- read.csv(file.path(ext, "python_nights.csv"), stringsAsFactors = FALSE)
cat(sprintf("\nnights: R = %d   Python = %d\n", nrow(nd), nrow(pn)))

bt_match  <- intersect(nd$bt, pn$bt)
py_only   <- pn[!(pn$bt %in% nd$bt), c("bt", "gt", "tbt", "nap")]
r_only    <- nd[!(nd$bt %in% pn$bt), c("bt", "gt")]
cat(sprintf("  matched bt: %d   Python-only: %d   R-only: %d\n",
            length(bt_match), nrow(py_only), nrow(r_only)))
if (nrow(py_only) > 0L) { cat("Python-only nights:\n"); print(py_only) }
if (nrow(r_only)  > 0L) { cat("R-only nights:\n");     print(r_only)  }
if (nrow(nd) == nrow(pn) && all(nd$bt == pn$bt)) {
  cat(sprintf("  gt match: %s\n", all(nd$gt == pn$gt)))
}

# ── 3. full chain: detect_naps_crespo -> compute_waso vs python_output$state ─
res     <- compute_waso(xn, wake_thresh = 60L)
state_r <- as.integer(res$data$state)

mism <- sum(state_r != state_final)
cat(sprintf("\nfull chain (naps + WASO) vs python_output$state: %d / %d mismatches  (agreement %.5f)\n",
            mism, n, mean(state_r == state_final)))
if (mism > 0L) {
  cat("mismatch breakdown (rows = python state, cols = R state):\n")
  print(table(python = state_final[state_r != state_final],
              R      = state_r[state_r != state_final]))
  # Where do the mismatches fall (on-wrist index ranges)?
  ow_idx  <- cumsum(onwrist)                 # on-wrist index at each position
  bad     <- which(state_r != state_final)
  bad_ow  <- ow_idx[bad][onwrist[bad]]
  if (length(bad_ow) > 0L) {
    cat(sprintf("  mismatch on-wrist index range: [%d, %d]   (%d on-wrist, %d off-wrist)\n",
                min(bad_ow), max(bad_ow), sum(onwrist[bad]), sum(!onwrist[bad])))
  }
}

cat("\ndone.\n")
