# dev/validate_cspd_stage4.R
#
# Validate the stage-4 CSPD refiner wiring (.cspd_refine_periods) against the
# Python fixtures in inst/extdata/. Run locally from the package root:
#
#   source("dev/validate_cspd_stage4.R")
#
# Reconstructs the on-wrist input exactly as dev/export_cspd_intermediates.py:
# python_output.csv, on-wrist = state != 4, raw activity + datetime. The full
# refiner output is compared to cspd_refined_output.csv (1 = wake, 0 = sleep
# period, post boolean-length-filter + border heuristic), and the per-night
# bedtime/getuptime indices to cspd_refined_sleep.csv. Intermediate stages are
# checked first so any divergence can be localised.

suppressMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    stop("install devtools, or adapt this script to source the R/ files directly")
  }
})

ext <- "inst/extdata"
rd1 <- function(f) as.numeric(scan(file.path(ext, f), quiet = TRUE))   # single-column %d fixtures

# ── reconstruct the on-wrist input (mirrors the exporter) ────────────────────
out <- read.csv(file.path(ext, "python_output.csv"), stringsAsFactors = FALSE)
out$datetime <- as.POSIXct(out$datetime, tz = "UTC")
onwrist  <- out$state != 4
activity <- as.numeric(out$activity[onwrist])
datetime <- out$datetime[onwrist]
cat(sprintf("rows: %d   on-wrist (state != 4): %d\n", nrow(out), sum(onwrist)))

msp_detection <- rd1("cspd_stage1_in.csv")
stopifnot(length(msp_detection) == sum(onwrist))

# internal accessors
ns    <- asNamespace("zeitR")
refine <- get(".cspd_refine_periods",        ns)
pvlf   <- get(".peak_valley_length_filter",  ns)
sgs    <- get(".sleep_gap_separation",       ns)
ddiff  <- get(".datetime_diff",              ns)

# ── staged intermediates (duration = 60 -> peak_valley_minimum_length = 11) ──
cat("\n── intermediate stages ──\n")
s1 <- pvlf(msp_detection, 11L)
cat(sprintf("stage1 (peak-valley length filter) mismatches: %d\n",
            sum(s1 != rd1("cspd_stage1_out.csv"))))

s2 <- sgs(s1, ddiff(datetime), 3600)
cat(sprintf("stage2 (sleep-gap separation)      mismatches: %d\n",
            sum(s2 != rd1("cspd_stage2_out.csv"))))

# ── full refiner ─────────────────────────────────────────────────────────────
cat("\n── full refiner (.cspd_refine_periods) ──\n")
res    <- refine(activity, datetime, msp_detection, condition = 0L)
ro_ref <- rd1("cspd_refined_output.csv")

cat(sprintf("refined_output length: R = %d   py = %d\n",
            length(res$refined_output), length(ro_ref)))
if (length(res$refined_output) == length(ro_ref)) {
  mm <- which(res$refined_output != ro_ref)
  cat(sprintf("refined_output mismatches: %d / %d  (%.4f%%)\n",
              length(mm), length(ro_ref), 100 * length(mm) / length(ro_ref)))
  cat(sprintf("epoch-level agreement: %.5f\n", mean(res$refined_output == ro_ref)))
  if (length(mm) > 0) cat("  first mismatch epochs (0-based):", head(mm - 1L, 25), "\n")
}

# ── transitions vs cspd_pre_transitions.csv ──────────────────────────────────
tr_ref <- read.csv(file.path(ext, "cspd_pre_transitions.csv"))
cat(sprintf("\ntransitions: R = %d   py = %d\n", nrow(res$transitions), nrow(tr_ref)))
if (nrow(res$transitions) == nrow(tr_ref)) {
  cat(sprintf("  index mismatches: %d   direction-sign mismatches: %d\n",
              sum(res$transitions$index != tr_ref$index),
              sum(sign(res$transitions$direction) != sign(tr_ref$direction))))
}

# ── refined_sleep_df vs cspd_refined_sleep.csv ───────────────────────────────
rs_path <- file.path(ext, "cspd_refined_sleep.csv")
if (file.exists(rs_path)) {
  rs_ref <- read.csv(rs_path)
  cat(sprintf("\nrefined_sleep pairs: R = %d   py = %d\n",
              nrow(res$refined_sleep_df), nrow(rs_ref)))
  k <- min(nrow(res$refined_sleep_df), nrow(rs_ref))
  if (k > 0) {
    cat(sprintf("  bedtime_index mismatches:  %d / %d\n",
                sum(res$refined_sleep_df$bedtime_index[1:k]  != rs_ref$bedtime_index[1:k]), k))
    cat(sprintf("  getuptime_index mismatches: %d / %d\n",
                sum(res$refined_sleep_df$getuptime_index[1:k] != rs_ref$getuptime_index[1:k]), k))
  }
}

cat("\ndone.\n")
