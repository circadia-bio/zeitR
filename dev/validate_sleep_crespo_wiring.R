# dev/validate_sleep_crespo_wiring.R
#
# Validate the full detect_sleep_crespo wiring (on-wrist MSP at q = 0.365 ->
# CSPD refiner -> sleep periods) against the Python fixtures. Run from the
# package root:
#
#   source("dev/validate_sleep_crespo_wiring.R")
#
# Two checks:
#   1. the on-wrist MSP (scaled activity, q = 0.365) reproduces the
#      CSPD-internal detect_msp (cspd_stage1_in.csv) — this is the only NEW
#      step relative to the already-validated stage-4 refiner.
#   2. detect_sleep_crespo(refine = TRUE) writes on-wrist state equal to
#      1 - cspd_refined_output.csv (the refined sleep periods), and preserves
#      off-wrist epochs.
#
# NOTE: this checks the refined sleep PERIODS, not python_output$state (which
# is post-WASO Cole-Kripke and is the next milestone).

suppressMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    stop("install devtools, or adapt this script to source the R/ files directly")
  }
})

ext <- "inst/extdata"
rd1 <- function(f) as.numeric(scan(file.path(ext, f), quiet = TRUE))

# ── reconstruct the detector input (off-wrist marked, on-wrist = state 0) ────
py <- read.csv(file.path(ext, "python_output.csv"), stringsAsFactors = FALSE)
py$datetime <- as.POSIXct(py$datetime, tz = "UTC")
x <- data.frame(
  datetime = py$datetime,
  activity = as.double(py$activity),
  state    = ifelse(py$state == 4L, 4L, 0L)
)
onwrist <- x$state != 4L
cat(sprintf("rows: %d   on-wrist (state != 4): %d\n", nrow(x), sum(onwrist)))

ns         <- asNamespace("zeitR")
crespo_msp <- get(".crespo_msp", ns)

# ── 1. on-wrist MSP (scaled, q = 0.365) vs cspd_stage1_in ────────────────────
ow_activity <- x$activity[onwrist]
secs        <- as.numeric(x$datetime[onwrist])
duration    <- as.numeric(names(sort(table(diff(secs)), decreasing = TRUE))[1])
epoch_h     <- 3600 / duration
cat(sprintf("duration = %g s   epoch_h = %g\n", duration, epoch_h))

msp <- crespo_msp(
  activity = ow_activity / duration, epoch_h = epoch_h,
  median_filter_h = 8, pad_h = 1, sleep_quantile = 0.365, morph_size = 61L,
  consec_zeros_thr = 15L, awake_zeros_thr = 2L, sleep_zeros_thr = 30L,
  zero_mitigation_q = 0.33, min_short_window_thr = 1.0
)
stage1_in <- rd1("cspd_stage1_in.csv")
cat(sprintf("MSP vs cspd_stage1_in mismatches: %d / %d\n",
            sum(msp != stage1_in), length(stage1_in)))

# ── 2. detect_sleep_crespo(refine = TRUE) vs refined sleep periods ───────────
out <- detect_sleep_crespo(x, refine = TRUE)

refined_output <- rd1("cspd_refined_output.csv")          # 1 = wake, 0 = sleep
expected_state <- ifelse(refined_output == 0, 1L, 0L)     # 1 = sleep, 0 = wake
ow_state       <- as.integer(out$state[onwrist])

cat(sprintf("on-wrist state vs refined periods mismatches: %d / %d  (agreement %.5f)\n",
            sum(ow_state != expected_state), length(expected_state),
            mean(ow_state == expected_state)))
cat(sprintf("off-wrist epochs preserved as state 4: %s\n",
            all(out$state[!onwrist] == 4L)))
cat(sprintf("sleep column tracks state == 1: %s\n",
            identical(as.integer(out$sleep), as.integer(out$state == 1L))))

cat("\ndone.\n")
