# dev/desc_trace.R  — trace the description-report filter on the two FP regions.
# Run: source("dev/desc_trace.R")  then paste output.

devtools::load_all(quiet = TRUE)
options(zeitR.debug_offwrist = TRUE)
rec  <- read_acttrust(system.file("extdata","input1.txt",package="zeitR"), tz="UTC")
prep <- detect_offwrist_bimodal(prepare_actigraphy(rec))

d <- if (exists(".zeitR_dbg2", envir = .GlobalEnv)) get(".zeitR_dbg2", envir = .GlobalEnv) else NULL
if (is.null(d)) { cat("No .zeitR_dbg2 captured — description filter block did not run (is_bimodal FALSE or no periods).\n") } else {
  cat(sprintf("ashman = %.4f   is_highly_separable = %s   is_bimodal = %s\n",
              d$ashman, d$is_highly_separable, d$is_bimodal))
  cat(sprintf("activity_thr = %.3f   sleep_lat = %.3f   long_offwrist_length = %d   offwrist_min_temp_dif_med = %.3f\n",
              d$activity_thr, d$sleep_lat, d$long_offwrist_length, d$offwrist_min_temp_dif_med))

  pre  <- d$pre_desc
  post <- d$post_desc
  rep  <- d$report
  cat(sprintf("\npre-filter periods: %d   post-filter periods: %d\n", nrow(pre), nrow(post)))

  # attach start/end to report for readability
  rep2 <- cbind(start = pre$start, end = pre$end, rep)

  show_region <- function(lo, hi, label) {
    cat(sprintf("\n=== %s  (periods overlapping [%d,%d]) ===\n", label, lo, hi))
    idx <- which(pre$start < hi & pre$end > lo)
    if (length(idx) == 0) { cat("  (none)\n"); return(invisible()) }
    for (i in idx) {
      r <- rep2[i, ]
      kept_post <- any(post$start == pre$start[i] & post$end == pre$end[i])
      cat(sprintf("  [%d,%d] len=%d vp=%s | zero_act=%.3f low_act=%.3f | hi_act_bef=%.3f hi_act_aft=%.3f | bord_act=%.3f bord_temp=%.3f | low_temp=%.3f tdm=%.3f | KEPT=%s\n",
                  pre$start[i], pre$end[i], pre$length[i],
                  isTRUE(r$valley_peak),
                  r$activity_zero_prop, r$low_act_prop,
                  r$high_act_before, r$high_act_after,
                  r$border_act_conc, r$border_temp_conc,
                  r$low_temp_prop, r$temp_dif_median, kept_post))
    }
  }
  show_region(60000, 60600, "60k FP region")
  show_region(65000, 65600, "65k FP region")
}
options(zeitR.debug_offwrist = FALSE)
