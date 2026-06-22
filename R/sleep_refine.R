# ── CSPD sleep-period refinement ──────────────────────────────────────────────
# Faithful R port of the Condor CSPD refinement stack
# (condor_pipeline/algorithms/vendor/condor/cspd_without_prints.py and the
#  cspd_bt_refine_/cspd_gt_refine_ refiners).
#
# Pipeline position: detect_sleep_crespo() produces the MSP detection; this
# stack refines it into the final sleep PERIODS (the Python `refined_output`).
# Per-epoch wake/sleep within those periods is then scored downstream by
# compute_waso() (Cole-Kripke), which is what python_output.csv$state reflects.
#
# Indexing follows the Python source: region `start`/`end` are 0-indexed with
# EXCLUSIVE ends. Convert to R vector positions with (start + 1L):end. This
# matches the convention already used in offwrist_refiner.R / .rle_periods().
#
# Stages (in CSPD.model order):
#   1. peak-valley length filter        -> .peak_valley_length_filter()   [DONE]
#   2. sleep-gap separation             -> .sleep_gap_separation()        [DONE]
#   3. bed-time / get-up-time refiners  -> .bt_refine() / .gt_refine()    [DONE]
#   4. minimum-length filter + wiring   -> .boolean_length_filter() and
#                                          .cspd_refine_periods()         [DONE]

# ── Peak/valley region helpers (port of cspd_functions_without_prints.py) ─────

#' Compute features of a peak/valley region (port of compute_features).
#'
#' `start`/`end` are 0-indexed with an exclusive end (Python convention); the
#' region is `activity[(start + 1L):end]` in R.
#' @noRd
.compute_features <- function(activity, threshold, start, end) {
  region <- activity[(start + 1L):end]
  len    <- end - start
  list(
    length                     = len,
    mean                       = mean(region),
    median                     = stats::median(region),
    zero_proportion            = zero_prop(region),
    above_threshold_proportion = sum(region > threshold) / len
  )
}

# ── Stage 3: get-up-time refiner (port of CSPD_GetUpTime_Refiner) ────────────
# The get-up refiner mirrors the bed-time one but reversed: it walks the
# peaks/valleys from the awake-most (last) region inward. It reuses the generic
# helpers .bt_datetime_gap_check / .bt_compute_refinement_window_median /
# .bt_compute_metric (the constructor stores metric_method/metric_parameter so
# those work unchanged) and adds get-up-specific methods. Indices follow the
# Python convention (0-indexed, exclusive ends).

#' Construct a get-up-time refiner object (mirrors CSPD_GetUpTime_Refiner.__init__).
#' @noRd
.cspd_getuptime_refiner <- function(activity,
                                    datetime_stamps,
                                    short_window_activity_median,
                                    minimum_short_window_activity_median_threshold,
                                    short_window_activity_median_minimum_high_epochs,
                                    half_window_around_border,
                                    activity_median_analysis_window,
                                    maximum_allowed_gap,
                                    quantile_threshold,
                                    median_filter_half_window_size,
                                    getuptime_metric_method,
                                    getuptime_metric_parameter,
                                    condition,
                                    getuptime_first_candidate_score,
                                    getuptime_best_median_difference_candidate_score,
                                    getuptime_best_crossing_distance_candidate_score,
                                    getuptime_best_epochs_above_metric_after_score,
                                    getuptime_thresholded_candidate_score_amplitude,
                                    getuptime_thresholded_candidate_score_minimum,
                                    getuptime_minimum_peak_length,
                                    getuptime_minimum_valley_length,
                                    after_candidate_window,
                                    getuptime_maximum_epochs_above_metric_after_candidate,
                                    zero_proportion_threshold,
                                    do_remove_after_long_tall_peak,
                                    getuptime_high_probability_awake_peak_length,
                                    getuptime_high_probability_sleep_valley_length,
                                    do_remove_before_long_valley,
                                    update_peaks_and_valleys,
                                    getuptime_score_first_candidate) {
  self <- new.env(parent = emptyenv())
  self$activity                       <- as.numeric(activity)
  self$datetime_stamps                <- datetime_stamps
  self$secs                           <- as.numeric(datetime_stamps)
  self$short_window_activity_median   <- as.numeric(short_window_activity_median)
  self$minimum_short_window_activity_median_threshold <- minimum_short_window_activity_median_threshold
  self$short_window_activity_median_minimum_high_epochs <- short_window_activity_median_minimum_high_epochs
  self$half_window_around_border      <- half_window_around_border
  self$activity_median_analysis_window <- activity_median_analysis_window
  self$maximum_allowed_gap            <- maximum_allowed_gap
  self$quantile_threshold             <- quantile_threshold
  self$median_filter_half_window_size <- median_filter_half_window_size
  # stored generically so the shared .bt_compute_metric works unchanged
  self$metric_method                  <- getuptime_metric_method
  self$metric_parameter               <- getuptime_metric_parameter
  self$condition                      <- condition
  self$getuptime_first_candidate_score                <- getuptime_first_candidate_score
  self$getuptime_best_median_difference_candidate_score <- getuptime_best_median_difference_candidate_score
  self$getuptime_best_crossing_distance_candidate_score <- getuptime_best_crossing_distance_candidate_score
  self$getuptime_best_epochs_above_metric_after_score   <- getuptime_best_epochs_above_metric_after_score
  self$getuptime_thresholded_candidate_score_amplitude  <- getuptime_thresholded_candidate_score_amplitude
  self$getuptime_thresholded_candidate_score_minimum    <- getuptime_thresholded_candidate_score_minimum
  self$getuptime_minimum_peak_length   <- getuptime_minimum_peak_length
  self$getuptime_minimum_valley_length <- getuptime_minimum_valley_length
  self$after_candidate_window          <- after_candidate_window
  self$getuptime_maximum_epochs_above_metric_after_candidate <- getuptime_maximum_epochs_above_metric_after_candidate
  self$zero_proportion_threshold       <- zero_proportion_threshold
  self$do_remove_after_long_tall_peak  <- do_remove_after_long_tall_peak
  self$getuptime_high_probability_awake_peak_length    <- getuptime_high_probability_awake_peak_length
  self$getuptime_high_probability_sleep_valley_length  <- getuptime_high_probability_sleep_valley_length
  self$do_remove_before_long_valley    <- do_remove_before_long_valley
  self$update_peaks_and_valleys        <- update_peaks_and_valleys
  self$getuptime_score_first_candidate <- getuptime_score_first_candidate
  self$data_length                     <- length(self$activity)
  self
}

#' compute_zero_proportion_around_start: zero-proportion in a window centred on
#' the window start. `start` is 0-indexed.
#' @noRd
.gt_compute_zero_proportion_around_start <- function(self, start) {
  hwab <- self$half_window_around_border
  N    <- self$data_length
  wend   <- if ((start + hwab) < N) start + hwab + 1L else N
  wstart <- if ((start - hwab) > 0L) start - hwab else 0L
  zero_prop(self$activity[(wstart + 1L):wend])
}

#' compute_refinement_window_end: stretch the window end forward until a
#' sustained high level of preceding short-window-median activity is reached
#' (subject to valid gaps and the next transition). Uses self$next_transition.
#' @noRd
.gt_compute_refinement_window_end <- function(self, initial_candidate) {
  rwe  <- initial_candidate
  she  <- self$short_window_activity_median_minimum_high_epochs
  qt   <- self$quantile_threshold
  swam <- self$short_window_activity_median

  valid_gap <- .bt_datetime_gap_check(self, rwe)
  before <- .np_slice(swam, rwe - she, rwe)
  high   <- 1 - .below_prop(before, qt)
  while ((rwe + 1L < self$data_length) && (rwe + 1L < self$next_transition) && valid_gap && (high < 1)) {
    rwe <- rwe + 1L
    valid_gap <- .bt_datetime_gap_check(self, rwe)
    before <- .np_slice(swam, rwe - she, rwe)
    high   <- 1 - .below_prop(before, qt)
  }
  rwe
}

#' compute_refinement_window_start: pull the start back while a low-activity /
#' high-zero-proportion region persists (subject to valid gaps and the previous
#' transition). Uses self$previous_transition.
#' @noRd
.gt_compute_refinement_window_start <- function(self, initial_candidate, refinement_window_end) {
  rws  <- initial_candidate
  if (rws < 0L) rws <- 0L
  amaw <- self$activity_median_analysis_window
  zpt  <- self$zero_proportion_threshold

  valid_gap <- .bt_datetime_gap_check(self, rws)
  zp  <- .gt_compute_zero_proportion_around_start(self, rws)
  med <- .bt_compute_refinement_window_median(self, rws, refinement_window_end)
  metric <- .bt_compute_metric(self, med)
  above  <- 1 - .below_prop(.np_slice(med, 0L, amaw), metric)

  while ((rws > 0L) && (rws - 1L > self$previous_transition) && valid_gap &&
         ((zp < zpt) || (above > 0))) {
    rws <- rws - 1L
    if (rws > 0L) valid_gap <- .bt_datetime_gap_check(self, rws)
    med <- .bt_compute_refinement_window_median(self, rws, refinement_window_end)
    metric <- .bt_compute_metric(self, med)
    above  <- 1 - .below_prop(.np_slice(med, 0L, amaw), metric)
    zp  <- .gt_compute_zero_proportion_around_start(self, rws)
  }
  rws
}

#' bridge_gap_validation: when the start search stops on an invalid gap, jump
#' before it and search for a fresh window start. Returns list(start, end).
#' @noRd
.gt_bridge_gap_validation <- function(self, refinement_window_start) {
  rwe <- refinement_window_start - 1L
  ic  <- refinement_window_start - 2L
  rws <- .gt_compute_refinement_window_start(self, ic, rwe)
  list(start = rws, end = rwe)
}

#' remove_after_long_tall_peak: drop everything after the first long AND tall
#' peak (length >= high-prob awake length and mean > 10 * metric), when metric > 5.
#' @noRd
.gt_remove_after_long_tall_peak <- function(self, pv) {
  metric <- self$metric
  invalid <- if (metric > 5) {
    which(pv$class == "p" &
          pv$length >= self$getuptime_high_probability_awake_peak_length &
          pv$mean > 10 * metric)
  } else integer(0)
  if (isTRUE(self$do_remove_after_long_tall_peak) && length(invalid) > 0) {
    pv <- pv[seq_len(invalid[1]), , drop = FALSE]    # keep 0..first invalid inclusive
  }
  rownames(pv) <- NULL
  pv
}

#' remove_before_long_valley: drop everything before the last "long" valley
#' (length >= high-prob sleep length), provided it is before the penultimate
#' region (last invalid index < count - 2).
#' @noRd
.gt_remove_before_long_valley <- function(self, pv) {
  cnt <- nrow(pv)
  thr <- self$getuptime_high_probability_sleep_valley_length
  invalid <- which(pv$class == "v" & !(pv$length < thr))
  if (isTRUE(self$do_remove_before_long_valley) && length(invalid) > 0) {
    last <- invalid[length(invalid)]
    if ((last - 1L) < cnt - 2L) pv <- pv[last:cnt, , drop = FALSE]   # keep last..end
  }
  rownames(pv) <- NULL
  pv
}

#' filter_peaks_and_valleys: reversed feature-based merging (mirror of the
#' bed-time loop, walking inward from the awake-most region). region indices are
#' 0-indexed via last_index - t; R rows are that + 1. Includes the upstream
#' "probable error" class-reassignment in the 3-region valley-peak-valley case,
#' replicated faithfully.
#' @noRd
.gt_filter_peaks_and_valleys <- function(self, pv) {
  gmvl   <- self$getuptime_minimum_valley_length
  gmpl   <- self$getuptime_minimum_peak_length
  rwa    <- self$refinement_window_activity
  metric <- self$metric

  cnt  <- nrow(pv)
  rownames(pv) <- NULL
  last <- cnt - 1L                                  # 0-indexed last_index

  if (cnt > 2) {
    if (pv$class[last + 1L] == "v") {
      if (pv$length[last + 1L] < gmvl) {
        pv <- .remove_peak_valley(pv, last, rwa, metric); cnt <- nrow(pv); last <- cnt - 1L
      } else if (cnt == 3) {
        if (pv$length[(last - 2L) + 1L] > 60 ||
            4 * pv$length[(last - 1L) + 1L] < pv$length[(last - 2L) + 1L]) {
          pv <- .remove_peak_valley(pv, last - 1L, rwa, metric); cnt <- nrow(pv); last <- cnt - 1L
          pv$class[last + 1L] <- "p"                # upstream "probable error" quirk
        }
      } else {
        if (pv$length[(last - 1L) + 1L] >= 20 &&
            3 * pv$median[last + 1L] < pv$median[(last - 1L) + 1L]) {
          pv <- .remove_peak_valley(pv, last, rwa, metric); cnt <- nrow(pv); last <- cnt - 1L
        }
      }
    }

    if (cnt > 3) {
      t <- 1L
      while (t < cnt && (cnt > 2 && cnt != 3)) {
        remove <- FALSE
        idx <- last - t                            # 0-indexed region
        ir  <- idx + 1L                            # R row
        if (pv$class[ir] == "p") {
          if (pv$length[ir] < gmpl) {
            if (t > 1L) {
              if (t < cnt - 2L) {
                remove <- TRUE
              } else if (pv$above_threshold_proportion[ir] < 0.8 &&
                         (pv$length[ir + 1L] >= 70 ||
                          (pv$length[ir + 1L] > gmvl &&
                           (pv$zero_proportion[ir + 1L] > 0.6 ||
                            (pv$above_threshold_proportion[ir + 1L] < 0.25 &&
                             pv$zero_proportion[ir + 1L] > 0.4))))) {
                remove <- TRUE
              }
            } else {
              if (pv$length[ir + 1L] >= 30 && pv$zero_proportion[ir + 1L] > 0.6) remove <- TRUE
            }
          }
        } else {                                   # valley
          if (t < cnt - 1L) {
            if (pv$length[ir] < gmvl) {
              if (t > 1L) {
                remove <- TRUE
              } else if (pv$above_fixed_threshold_proportion[ir] > 0.1 &&
                         pv$zero_proportion[ir - 1L] < 0.1) {
                remove <- TRUE
              }
            } else {
              remove_points <- 0
              if (pv$above_threshold_proportion[ir] > 0.1) remove_points <- remove_points + 1
              if (pv$zero_proportion[ir] < 0.25) remove_points <- remove_points + 1
              if (pv$mean[ir] >= 0.66 * metric) remove_points <- remove_points + 0.5
              if (pv$length[ir - 1L] > pv$length[ir] && pv$mean[ir - 1L] > 3 * pv$mean[ir]) {
                remove_points <- remove_points + 1
              }
              if (pv$length[ir + 1L] > 2 * pv$length[ir] && pv$mean[ir + 1L] > 3 * pv$mean[ir]) {
                remove_points <- remove_points + 1
              }
              if (pv$length[ir] / length(rwa) >= 0.3) remove_points <- remove_points - 1
              if (remove_points > 2) remove <- TRUE
            }
          }
        }

        if (remove) {
          pv <- .remove_peak_valley(pv, last - t, rwa, metric); cnt <- nrow(pv); last <- cnt - 1L
        } else {
          t <- t + 1L
        }
      }
    }
  }
  rownames(pv) <- NULL
  pv
}

#' identify_getuptime_candidates: peak starts (for non-first regions) are
#' candidate get-up times; if the last region is a valley, its end is also a
#' candidate. Mutates self$refinement_window_levels. Returns window-relative,
#' 0-indexed candidate positions.
#' @noRd
.gt_identify_getuptime_candidates <- function(self, pv) {
  cnt <- nrow(pv)
  rws <- self$refinement_window_start
  cands <- integer(0)
  for (t in seq_len(cnt)) {
    t0    <- t - 1L
    start <- as.integer(rws + pv$start[t]); end <- as.integer(rws + pv$end[t])
    if (end > start) self$refinement_window_levels[(start + 1L):end] <- pv$mean[t]
    if (pv$class[t] == "p" && t0 > 0L) cands <- c(cands, as.integer(pv$start[t]))
  }
  if (pv$class[cnt] == "v") {
    new <- as.integer(pv$end[cnt])
    if (rws + new == self$data_length) new <- new - 1L
    cands <- c(cands, new)
  }
  cands
}

#' choose_best_getuptime_candidate: score and pick the refined get-up time
#' (port of choose_best_getuptime_candidate). Mirrors the bed-time scorer but:
#' rewards the FIRST candidate, sorts by the activity-median-difference PEAK
#' (rise) descending, uses the FIRST upward crossing, counts epochs BELOW the
#' metric, applies bonuses by candidate value, and has no second-best override.
#' @noRd
.gt_choose_best_getuptime_candidate <- function(self, candidates, up) {
  n <- length(candidates)
  if (n == 0L) return(as.integer(self$initial_transition_candidate))

  rws      <- self$refinement_window_start
  smoothed <- self$refinement_window_activity_median_difference_smoothed

  amd <- vapply(candidates, function(c0) .get_peak(smoothed, c0, valley = FALSE), numeric(1))
  sc <- data.frame(
    candidate   = as.integer(candidates),
    amd         = as.numeric(amd),
    epochs      = rep(self$after_candidate_window, n),
    thresholded = rep(FALSE, n),
    gap_after   = rep(TRUE, n),
    score       = rep(0, n),
    stringsAsFactors = FALSE
  )

  if (isTRUE(self$getuptime_score_first_candidate)) {
    sc$score[1] <- sc$score[1] + self$getuptime_first_candidate_score
  }

  sc <- sc[order(-sc$amd, -sc$candidate), , drop = FALSE]       # amd desc, candidate desc
  sc$score[1] <- sc$score[1] + self$getuptime_best_median_difference_candidate_score

  if (length(up) > 0) {
    cu   <- up[1]
    bcdc <- which.min(abs(as.integer(candidates) - cu))         # nearest the first up-crossing
    tgt  <- as.integer(candidates[bcdc])
    sc$score[sc$candidate == tgt] <- sc$score[sc$candidate == tgt] +
      self$getuptime_best_crossing_distance_candidate_score
  }

  acw    <- self$after_candidate_window
  maxep  <- self$getuptime_maximum_epochs_above_metric_after_candidate
  metric <- self$metric
  for (ci in seq_len(n)) {
    cand <- rws + as.integer(candidates[ci])
    end  <- cand + acw
    if (end > self$data_length) end <- self$data_length
    valid <- TRUE
    g <- cand
    while (valid && (g + 1L < end)) {
      valid <- .bt_datetime_gap_check(self, g, direction = "forward")
      g <- g + 1L
    }
    if (valid) {
      tgt <- as.integer(candidates[ci])
      ea  <- if (end > cand) sum(self$activity[(cand + 1L):end] <= metric) else 0L
      sc$gap_after[sc$candidate == tgt] <- FALSE
      sc$epochs[sc$candidate == tgt] <- ea
      if (ea <= maxep) sc$thresholded[sc$candidate == tgt] <- TRUE
    }
  }

  if (sum(!sc$gap_after) > 0) {
    if (sum(sc$thresholded) == 0) {
      sc <- sc[order(sc$gap_after, sc$epochs), , drop = FALSE]
      sc$score[1] <- sc$score[1] + self$getuptime_best_epochs_above_metric_after_score
    } else {
      factor <- if (maxep > 0) {
        -self$getuptime_thresholded_candidate_score_amplitude / maxep
      } else {
        -self$getuptime_thresholded_candidate_score_amplitude / 1e-3
      }
      thr <- which(sc$thresholded)
      sc$score[thr] <- sc$score[thr] +
        factor * (sc$epochs[thr] - maxep) + self$getuptime_thresholded_candidate_score_minimum
    }
  }

  sc <- sc[order(-sc$score, -sc$candidate), , drop = FALSE]
  as.integer(rws + as.integer(sc$candidate[1]))
}

#' refine: full get-up-time refinement for one transition (port of
#' CSPD_GetUpTime_Refiner.refine). Returns a list mirroring the Python tuple:
#' refined_getuptime, refinement_window_start, refinement_window_end,
#' refinement_window_activity_median,
#' refinement_window_activity_median_difference_smoothed,
#' refinement_window_levels, metric. Validated bit-exact (refined index, window
#' start, window end) against the real refiner across all 52 transitions.
#' @noRd
.gt_refine <- function(self, refinement_window_levels, initial_transition_candidate,
                       previous_transition, next_transition) {
  self$refinement_window_levels     <- refinement_window_levels
  self$initial_transition_candidate <- initial_transition_candidate
  self$previous_transition          <- previous_transition
  self$next_transition              <- next_transition

  rwe <- .gt_compute_refinement_window_end(self, initial_transition_candidate + 1L)
  rws <- .gt_compute_refinement_window_start(self, initial_transition_candidate - 1L, rwe)

  if (!.bt_datetime_gap_check(self, rws)) {
    br  <- .gt_bridge_gap_validation(self, rws)
    rws <- br$start; rwe <- br$end
  }

  rwdd   <- .datetime_diff(self$secs[(rws + 1L):(rwe + 1L)])
  rwgaps <- as.integer(rwdd >= self$maximum_allowed_gap)
  if (sum(rwgaps) > 0) {
    fg <- which(rwgaps == 1L)[1] - 1L
    if (fg > length(rwdd) - fg) rwe <- rws + fg - 1L else rws <- rws + fg + 1L
  }
  if (rws >= rwe) {
    rwe <- rws + 60L
    if (rwe >= self$data_length) rwe <- self$data_length - 1L
  }

  self$refinement_window_start <- rws
  self$refinement_window_end   <- rwe
  self$refinement_window_activity_median <- .bt_compute_refinement_window_median(self, rws, rwe)
  self$metric <- .bt_compute_metric(self, self$refinement_window_activity_median)
  self$refinement_window_activity_median_difference <- diff(self$refinement_window_activity_median)
  self$refinement_window_activity_median_difference_smoothed <-
    .cspd_median_filter(self$refinement_window_activity_median_difference, self$median_filter_half_window_size)
  self$refinement_window_activity <- self$activity[(rws + 1L):(rwe + 1L)]

  pv  <- .identify_peaks_and_valleys(signal = self$refinement_window_activity_median,
                                     activity = self$refinement_window_activity,
                                     threshold = self$metric)
  rwa <- self$refinement_window_activity
  pv$above_fixed_threshold_proportion <- vapply(seq_len(nrow(pv)), function(i) {
    s <- pv$start[i]; e <- pv$end[i]
    sum(rwa[(s + 1L):e] > 10) / (e - s)
  }, numeric(1))

  pv <- .gt_remove_after_long_tall_peak(self, pv)
  pv <- .gt_remove_before_long_valley(self, pv)
  self$peaks_and_valleys <- .gt_filter_peaks_and_valleys(self, pv)
  cnt <- nrow(self$peaks_and_valleys)

  new_rws <- self$refinement_window_start + self$peaks_and_valleys$start[1]
  self$refinement_window_end <- self$refinement_window_start + self$peaks_and_valleys$end[cnt]
  if (self$refinement_window_end >= self$data_length) self$refinement_window_end <- self$data_length - 1L
  self$refinement_window_start <- new_rws

  self$refinement_window_activity_median <-
    .bt_compute_refinement_window_median(self, self$refinement_window_start, self$refinement_window_end)
  self$metric <- .bt_compute_metric(self, self$refinement_window_activity_median)
  self$refinement_window_activity_median_difference <- diff(self$refinement_window_activity_median)
  self$refinement_window_activity_median_difference_smoothed <-
    .cspd_median_filter(self$refinement_window_activity_median_difference, self$median_filter_half_window_size)
  self$refinement_window_activity <-
    self$activity[(self$refinement_window_start + 1L):(self$refinement_window_end + 1L)]

  if (isTRUE(self$update_peaks_and_valleys)) {
    self$peaks_and_valleys <- .identify_peaks_and_valleys(signal = self$refinement_window_activity_median,
                                                          activity = self$refinement_window_activity,
                                                          threshold = self$metric)
  } else {
    self$peaks_and_valleys$start[1] <- 0L
    self$peaks_and_valleys$end[1]   <- self$peaks_and_valleys$start[1] + self$peaks_and_valleys$length[1]
    if (cnt >= 2L) {
      for (ii in 2:cnt) {
        self$peaks_and_valleys$start[ii] <- self$peaks_and_valleys$end[ii - 1L]
        self$peaks_and_valleys$end[ii]   <- self$peaks_and_valleys$start[ii] + self$peaks_and_valleys$length[ii]
      }
    }
  }

  cands <- .gt_identify_getuptime_candidates(self, self$peaks_and_valleys)

  seg <- self$short_window_activity_median[(self$refinement_window_start + 1L):(self$refinement_window_end + 1L)]
  cr  <- diff(c(0L, as.integer(seg >= self$metric)))
  up  <- which(cr > 0) - 1L

  refined_getuptime <- .gt_choose_best_getuptime_candidate(self, cands, up)

  list(
    refined_getuptime                      = refined_getuptime,
    refinement_window_start                = self$refinement_window_start,
    refinement_window_end                  = self$refinement_window_end,
    refinement_window_activity_median      = self$refinement_window_activity_median,
    refinement_window_activity_median_difference_smoothed =
      self$refinement_window_activity_median_difference_smoothed,
    refinement_window_levels               = self$refinement_window_levels,
    metric                                 = self$metric
  )
}

#' Identify contiguous peak ("p") and valley ("v") regions
#' (port of identify_peaks_and_valleys).
#'
#' `signal` is thresholded at `threshold` (`signal > threshold` -> peak). Runs
#' of high values are peaks ("p"), runs of low values are valleys ("v").
#' Returns a data.frame with 0-indexed, exclusive-end `start`/`end`. Activity
#' features (mean/median/zero_proportion/above_threshold_proportion) are
#' computed only when `activity` is supplied — the length filter does not need
#' them, but the bed/get-up refiners do.
#' @noRd
.identify_peaks_and_valleys <- function(signal, activity = NULL, threshold = 0.5) {
  n  <- length(signal)
  th <- as.integer(signal > threshold)

  cls <- character(0); starts <- integer(0); ends <- integer(0)
  t <- 1L                                  # R 1-indexed cursor (Python t, 0-indexed)
  while (t <= n) {
    start0 <- t - 1L                       # 0-indexed start of this region
    if (th[t] == 1L) {                     # peak (sustained high)
      t <- t + 1L
      while (t <= n && th[t] == 1L) t <- t + 1L
      region_class <- "p"
    } else {                               # valley (sustained low)
      t <- t + 1L
      while (t <= n && th[t] == 0L) t <- t + 1L
      region_class <- "v"
    }
    end0   <- t - 1L                       # 0-indexed exclusive end
    cls    <- c(cls, region_class)
    starts <- c(starts, start0)
    ends   <- c(ends, end0)
  }

  df <- data.frame(
    class  = cls,
    start  = starts,
    end    = ends,
    length = ends - starts,
    stringsAsFactors = FALSE
  )

  if (!is.null(activity)) {
    feats <- lapply(seq_len(nrow(df)), function(i)
      .compute_features(activity, threshold, df$start[i], df$end[i]))
    df$mean                       <- vapply(feats, `[[`, numeric(1), "mean")
    df$median                     <- vapply(feats, `[[`, numeric(1), "median")
    df$zero_proportion            <- vapply(feats, `[[`, numeric(1), "zero_proportion")
    df$above_threshold_proportion <- vapply(feats, `[[`, numeric(1), "above_threshold_proportion")
  }

  df
}

# ── Stage 1: peak-valley length filter ────────────────────────────────────────

#' Peak-valley length filter (inline stage-1 filter from CSPD.model).
#'
#' Short peaks/valleys (`length <= min_length`) are merged into their
#' neighbours, then the signal is rebuilt as 1 = wake, 0 = sleep (valleys).
#' Faithful to the inline loop in `cspd_without_prints.py`: the first region is
#' special-cased (merged forward into the second), the main pass starts at the
#' second region, and merging updates only `start`/`end`/`length` (activity
#' features are not recomputed, unlike `remove_peak_valley`).
#'
#' @param detection integer 0/1 vector (1 = wake, 0 = sleep) from `.crespo_msp`.
#' @param min_length integer; regions this short or shorter are removed.
#' @return integer 0/1 vector (1 = wake, 0 = sleep) of `length(detection)`.
#' @noRd
.peak_valley_length_filter <- function(detection, min_length) {
  n  <- length(detection)
  pv <- .identify_peaks_and_valleys(detection)   # features not needed here
  cls <- pv$class; start <- pv$start; end <- pv$end; len <- pv$length
  cnt <- length(cls)

  # Special-case the first region: merge it forward into the second, keeping the
  # second region's class (Python drops index 0, reindexes, then loops from 1).
  if (cnt > 1L && len[1] <= min_length) {
    s <- start[1]; e <- end[2]
    start[2] <- s; len[2] <- e - s
    keep <- -1L
    cls <- cls[keep]; start <- start[keep]; end <- end[keep]; len <- len[keep]
    cnt <- length(cls)
  }

  ri <- 2L                                        # R index == Python region_index 1
  while (ri <= cnt) {
    if (len[ri] <= min_length) {
      s <- start[ri - 1L]
      if (ri < cnt) {                             # not the last region: merge ri-1, ri, ri+1
        e    <- end[ri + 1L]
        drop <- c(ri, ri + 1L)
      } else {                                    # last region: merge ri-1, ri
        e    <- end[ri]
        drop <- ri
      }
      end[ri - 1L] <- e
      len[ri - 1L] <- e - s
      keep <- setdiff(seq_len(cnt), drop)
      cls <- cls[keep]; start <- start[keep]; end <- end[keep]; len <- len[keep]
      cnt <- length(cls)
      # region_index unchanged: re-check the same position (now a different region)
    } else {
      ri <- ri + 1L
    }
  }

  out <- rep(1L, n)
  for (i in which(cls == "v")) {
    out[(start[i] + 1L):end[i]] <- 0L
  }
  out
}

# ── Stage 2: sleep-gap separation ─────────────────────────────────────────────

#' Per-epoch datetime gap in seconds (port of datetime_diff).
#'
#' Returns a numeric vector the same length as `stamps`; element 1 is 0 and
#' element i (i > 1) is the seconds between `stamps[i]` and `stamps[i - 1]`.
#' @param stamps POSIXct (or numeric seconds) vector of epoch timestamps.
#' @noRd
.datetime_diff <- function(stamps) {
  s <- as.numeric(stamps)            # POSIXct -> seconds since epoch
  c(0, diff(s))
}

#' numpy-style slice assignment: x[start:stop] <- value with 0-indexed,
#' exclusive-end semantics, including negative-index wrap and out-of-range
#' clamping (matches NumPy exactly, so boundary gaps behave as in Python).
#' @noRd
.np_slice_assign <- function(x, start, stop, value) {
  n <- length(x)
  if (start < 0L) start <- start + n
  if (stop  < 0L) stop  <- stop  + n
  start <- max(0L, min(start, n))
  stop  <- max(0L, min(stop,  n))
  if (start < stop) x[(start + 1L):stop] <- value   # 0-indexed [start, stop) -> R
  x
}

#' Sleep-gap separation (stage 2 from CSPD.model).
#'
#' A "big" datetime gap inside a sleep period (gap > `max_gap_seconds` at an
#' epoch currently scored as sleep) probably joins two distinct sleep periods
#' across removed/off-wrist time, so a +/-10 epoch window around it is forced to
#' wake. Faithful to the Python loop `final_sleep_detection[l-10:l+11] = 1`,
#' including NumPy's slice semantics at the array boundaries.
#'
#' @param detection integer 0/1 vector (1 = wake, 0 = sleep) after stage 1.
#' @param datetime_diff numeric per-epoch gap in seconds (see `.datetime_diff`).
#' @param max_gap_seconds numeric threshold (CSPD default 60 * 60 = 3600).
#' @return integer 0/1 vector (1 = wake, 0 = sleep) of `length(detection)`.
#' @noRd
.sleep_gap_separation <- function(detection, datetime_diff, max_gap_seconds) {
  out <- detection
  # Python: for l in range(n) if dd[l] > thr and det[l] == 0  (0-indexed l)
  gaps0 <- which(datetime_diff > max_gap_seconds & out == 0) - 1L
  for (l in gaps0) {
    out <- .np_slice_assign(out, l - 10L, l + 11L, 1L)
  }
  out
}

# ── Stage 3: shared refiner helpers ───────────────────────────────────────────
# Deterministic building blocks the bed-time / get-up-time refiners depend on.
# Ported from functions.py and cspd_functions_without_prints.py. Indices follow
# the Python convention (0-indexed, exclusive ends); convert at access.

#' Proportion of elements <= threshold (port of below_prop).
#' @noRd
.below_prop <- function(a, thresh) {
  if (length(a) == 0L) return(0)
  sum(a <= thresh) / length(a)
}

#' Local extremum value in a +/- peak_hws window around `center` (port of
#' get_peak). `center` is 0-indexed. Returns the max (or min, if `valley`)
#' value of the clamped window; argmax/argmin take the first index, matching
#' NumPy and R's which.max/which.min.
#' @noRd
.get_peak <- function(signal, center, valley = FALSE, peak_hws = 10L) {
  n  <- length(signal)
  lo <- center - peak_hws            # 0-indexed inclusive start
  hi <- center + peak_hws + 1L       # 0-indexed exclusive end
  if (lo >= 0L) {
    if (hi <= n) {
      ps <- signal[(lo + 1L):hi]
    } else {
      ps <- signal[(lo + 1L):n]
    }
  } else {
    ps <- signal[1L:hi]              # signal[0:center+hws+1]
  }
  if (valley) ps[which.min(ps)] else ps[which.max(ps)]
}

#' Median filter matching functions.py median_filter (center = TRUE).
#'
#' Distinct from utils.R's `median_filter`: this is the CSPD-refiner variant
#' with selectable padding. `padding = NULL` replicates the first/last value
#' for 2*hws on each side; `padding = "padded"` assumes the signal is already
#' padded and trims 4*hws from the output length; a length-2 numeric vector
#' pads with constant front/back values.
#' @noRd
.cspd_median_filter <- function(signal, hws, padding = NULL) {
  s <- as.numeric(signal)
  n <- length(s)
  if (n <= hws) return(s)

  med_win <- function(v, count) {
    vapply(seq_len(count), function(i) {
      stats::median(v[(hws + i):(hws + i + 2L * hws)])
    }, numeric(1))
  }

  if (is.null(padding)) {
    filt <- c(rep(s[1], 2L * hws), s, rep(s[n], 2L * hws))
    med_win(filt, n)
  } else if (is.character(padding) && padding == "padded") {
    med_win(s, n - 4L * hws)
  } else {
    filt <- c(rep(padding[1], 2L * hws), s, rep(padding[2], 2L * hws))
    med_win(filt, n)
  }
}

#' Merge a peak/valley into its neighbours, recomputing features
#' (port of remove_peak_valley). `index_to_remove` is 0-indexed.
#'
#' First region merges forward into the second; the last merges back into the
#' penultimate; an interior region merges its two neighbours around it into one
#' (keeping the preceding region's class). Unlike the inline stage-1 filter,
#' this recomputes mean/median/zero_proportion/above_threshold_proportion via
#' `.compute_features`.
#' @noRd
.remove_peak_valley <- function(peaks_and_valleys, index_to_remove, activity, threshold) {
  pv  <- peaks_and_valleys
  cnt <- nrow(pv)
  i0  <- index_to_remove             # 0-indexed

  set_features <- function(row, start, end) {
    f <- .compute_features(activity, threshold, start, end)
    pv$length[row]                     <<- as.integer(f$length)
    pv$mean[row]                       <<- f$mean
    pv$median[row]                     <<- f$median
    pv$zero_proportion[row]            <<- f$zero_proportion
    pv$above_threshold_proportion[row] <<- f$above_threshold_proportion
  }

  if (i0 == 0L) {
    start <- as.integer(pv$start[1]); end <- as.integer(pv$end[2])
    pv$start[2] <- start
    set_features(2L, start, end)
    pv <- pv[-1L, , drop = FALSE]
  } else {
    r <- i0 + 1L                       # R row of index_to_remove
    start <- as.integer(pv$start[r - 1L])
    if (i0 < cnt - 1L) {
      end  <- as.integer(pv$end[r + 1L])
      drop <- c(r, r + 1L)
    } else {
      end  <- as.integer(pv$end[r])
      drop <- r
    }
    pv$end[r - 1L] <- end
    set_features(r - 1L, start, end)
    pv <- pv[-drop, , drop = FALSE]
  }

  rownames(pv) <- NULL
  pv
}

#' numpy-style slice read: x[start:stop] with 0-indexed, exclusive-end
#' semantics, negative-index wrap and out-of-range clamping. Read counterpart
#' of `.np_slice_assign`.
#' @noRd
.np_slice <- function(x, start, stop) {
  n <- length(x)
  if (start < 0L) start <- start + n
  if (stop  < 0L) stop  <- stop  + n
  start <- max(0L, min(start, n))
  stop  <- max(0L, min(stop,  n))
  if (start < stop) x[(start + 1L):stop] else x[integer(0)]
}

# ── Stage 3: bed-time refiner ────────────────────────────────────────────────
# Port of CSPD_BedTime_Refiner (cspd_bt_refine_without_prints.py). The Python
# class is mirrored by an environment-based object: `.cspd_bedtime_refiner()`
# builds `self`; methods are `.bt_*(self, ...)` functions that read config and
# the activity / datetime / short-window-median arrays, and (in later stages)
# mutate refinement state on `self`. All indices follow the Python convention
# (0-indexed, exclusive ends); convert at vector access with (i + 1L).
#
# Unit A (this section): refinement-window discovery + metric.
#   datetime_gap_check, compute_refinement_window_median,
#   compute_zero_proportion_around_end, compute_metric,
#   compute_initial_refinement_window_start, compute_refinement_window_end,
#   bridge_gap_validation, compute_improved_refinement_window_start.

#' Construct a bed-time refiner object (mirrors CSPD_BedTime_Refiner.__init__).
#' @noRd
.cspd_bedtime_refiner <- function(activity,
                                  datetime_stamps,
                                  short_window_activity_median,
                                  minimum_short_window_activity_median_threshold,
                                  short_window_activity_median_minimum_high_epochs,
                                  half_window_around_border,
                                  activity_median_analysis_window,
                                  maximum_allowed_gap,
                                  quantile_threshold,
                                  median_filter_half_window_size,
                                  metric_method,
                                  metric_parameter,
                                  condition,
                                  bedtime_score_last_candidate,
                                  bedtime_last_candidate_score,
                                  bedtime_best_median_difference_candidate_score,
                                  bedtime_best_crossing_distance_candidate_score,
                                  bedtime_best_epochs_above_metric_after_score,
                                  bedtime_thresholded_candidate_score_amplitude,
                                  bedtime_thresholded_candidate_score_minimum,
                                  bedtime_minimum_peak_length,
                                  bedtime_minimum_valley_length,
                                  after_candidate_window,
                                  bedtime_maximum_epochs_above_metric_after_candidate,
                                  zero_proportion_threshold,
                                  do_remove_after_long_valley,
                                  do_remove_before_long_peak,
                                  do_remove_before_tall_peak,
                                  update_peaks_and_valleys,
                                  do_bedtime_candidates_crossings_filter,
                                  consider_second_best_candidate,
                                  bedtime_high_probability_awake_peak_length,
                                  bedtime_high_probability_sleep_valley_length) {
  self <- new.env(parent = emptyenv())
  self$activity                       <- as.numeric(activity)
  self$datetime_stamps                <- datetime_stamps
  self$secs                           <- as.numeric(datetime_stamps)   # seconds for gap arithmetic
  self$short_window_activity_median   <- as.numeric(short_window_activity_median)
  self$minimum_short_window_activity_median_threshold <- minimum_short_window_activity_median_threshold
  self$short_window_activity_median_minimum_high_epochs <- short_window_activity_median_minimum_high_epochs
  self$half_window_around_border      <- half_window_around_border
  self$activity_median_analysis_window <- activity_median_analysis_window
  self$maximum_allowed_gap            <- maximum_allowed_gap
  self$quantile_threshold             <- quantile_threshold
  self$median_filter_half_window_size <- median_filter_half_window_size
  self$metric_method                  <- metric_method
  self$metric_parameter               <- metric_parameter
  self$condition                      <- condition
  self$bedtime_score_last_candidate   <- bedtime_score_last_candidate
  self$bedtime_last_candidate_score   <- bedtime_last_candidate_score
  self$bedtime_best_median_difference_candidate_score  <- bedtime_best_median_difference_candidate_score
  self$bedtime_best_crossing_distance_candidate_score  <- bedtime_best_crossing_distance_candidate_score
  self$bedtime_best_epochs_above_metric_after_score    <- bedtime_best_epochs_above_metric_after_score
  self$bedtime_thresholded_candidate_score_amplitude   <- bedtime_thresholded_candidate_score_amplitude
  self$bedtime_thresholded_candidate_score_minimum     <- bedtime_thresholded_candidate_score_minimum
  self$bedtime_minimum_peak_length    <- bedtime_minimum_peak_length
  self$bedtime_minimum_valley_length  <- bedtime_minimum_valley_length
  self$after_candidate_window         <- after_candidate_window
  self$bedtime_maximum_epochs_above_metric_after_candidate <- bedtime_maximum_epochs_above_metric_after_candidate
  self$zero_proportion_threshold      <- zero_proportion_threshold
  self$do_remove_after_long_valley    <- do_remove_after_long_valley
  self$do_remove_before_long_peak     <- do_remove_before_long_peak
  self$do_remove_before_tall_peak     <- do_remove_before_tall_peak
  self$update_peaks_and_valleys       <- update_peaks_and_valleys
  self$do_bedtime_candidates_crossings_filter <- do_bedtime_candidates_crossings_filter
  self$consider_second_best_candidate <- consider_second_best_candidate
  self$bedtime_high_probability_awake_peak_length      <- bedtime_high_probability_awake_peak_length
  self$bedtime_high_probability_sleep_valley_length    <- bedtime_high_probability_sleep_valley_length
  self$data_length                    <- length(self$activity)
  self
}

#' datetime_gap_check: is the gap to the previous/next epoch valid?
#' `location` is 0-indexed. Returns logical, or list(valid, gap) if return_gap.
#' @noRd
.bt_datetime_gap_check <- function(self, location, direction = "backward", return_gap = FALSE) {
  if (direction == "backward") {
    gap <- self$secs[location + 1L] - self$secs[location]          # secs[loc] - secs[loc-1] (0-indexed)
  } else {
    gap <- self$secs[location + 2L] - self$secs[location + 1L]     # secs[loc+1] - secs[loc]
  }
  valid <- !(gap > self$maximum_allowed_gap)
  if (return_gap) list(valid = valid, gap = gap) else valid
}

#' compute_refinement_window_median: median-filtered activity in a window,
#' pre-padding the raw activity then using the "padded" median filter.
#' @noRd
.bt_compute_refinement_window_median <- function(self, start, end) {
  hws <- self$median_filter_half_window_size
  N   <- self$data_length
  act <- self$activity
  if (start - 2L * hws >= 0L) {
    if (end + 2L * hws + 1L <= N) {
      rwa <- act[(start - 2L * hws + 1L):(end + 2L * hws + 1L)]
    } else {
      rwa <- act[(start - hws + 1L):N]
      mx  <- max(act[(start + 1L):(end + 1L)])
      need <- end + 1L + 4L * hws - start
      if (length(rwa) < need) rwa <- c(rwa, rep(mx, need - length(rwa)))
    }
  } else {
    rwa <- act[1L:(end + 2L * hws + 1L)]
    mx  <- max(act[(start + 1L):(end + 1L)])
    need <- end + 1L + 4L * hws - start
    if (length(rwa) < need) rwa <- c(rep(mx, need - length(rwa)), rwa)
  }
  .cspd_median_filter(rwa, hws, padding = "padded")
}

#' compute_zero_proportion_around_end: zero-proportion in a window centred on
#' the window end. `end` is 0-indexed.
#' @noRd
.bt_compute_zero_proportion_around_end <- function(self, end) {
  hwab <- self$half_window_around_border
  N    <- self$data_length
  wend   <- if ((end + hwab) < N) end + hwab + 1L else N
  wstart <- if ((end - hwab) > 0L) end - hwab else 0L
  zero_prop(self$activity[(wstart + 1L):wend])
}

#' compute_metric: high/low activity separation threshold within a window.
#' @noRd
.bt_compute_metric <- function(self, refinement_window_activity_median) {
  m   <- refinement_window_activity_median
  pos <- m[m > 0]
  if (length(pos) == 0L) return(0)
  if (self$metric_method == 1) {
    self$metric_parameter * mean(pos)
  } else if (self$metric_method == 2) {
    stats::quantile(pos, self$metric_parameter, names = FALSE, type = 7)
  } else {
    0
  }
}

#' compute_initial_refinement_window_start: stretch the window start back to
#' the first epoch with a high level of preceding short-window-median activity
#' (subject to valid gaps and the previous transition). Uses self$previous_transition.
#' @noRd
.bt_compute_initial_refinement_window_start <- function(self, initial_candidate) {
  rws <- initial_candidate
  if (rws < 0L) rws <- 0L
  swam <- self$short_window_activity_median
  she  <- self$short_window_activity_median_minimum_high_epochs
  qt   <- self$quantile_threshold

  valid_gap <- .bt_datetime_gap_check(self, rws)
  before <- if (rws >= she) swam[(rws - she + 1L):rws] else swam[seq_len(rws)]
  high_prop <- 1 - .below_prop(before, qt)

  while ((rws > 0L) && (rws - 1L > self$previous_transition) && valid_gap && (high_prop < 1)) {
    rws <- rws - 1L
    if (rws > 0L) valid_gap <- .bt_datetime_gap_check(self, rws)
    before <- if (rws >= she) swam[(rws - she + 1L):rws] else swam[seq_len(rws)]
    high_prop <- 1 - .below_prop(before, qt)
  }
  rws
}

#' compute_refinement_window_end: stretch the window end forward until a
#' sustained low-activity / high-zero-proportion region is reached (subject to
#' valid gaps and the next transition). Uses self$next_transition.
#' @noRd
.bt_compute_refinement_window_end <- function(self, initial_candidate, refinement_window_start) {
  rwe  <- initial_candidate
  amaw <- self$activity_median_analysis_window
  zpt  <- self$zero_proportion_threshold

  valid_gap <- .bt_datetime_gap_check(self, rwe)
  zp  <- .bt_compute_zero_proportion_around_end(self, rwe)
  med <- .bt_compute_refinement_window_median(self, refinement_window_start, rwe)
  metric <- .bt_compute_metric(self, med)
  ending <- .np_slice(med, length(med) - amaw, length(med))
  above  <- 1 - .below_prop(ending, metric)

  while ((rwe + 1L < self$data_length) && (rwe + 1L < self$next_transition) && valid_gap &&
         ((zp < zpt) || (above > 0))) {
    rwe <- rwe + 1L
    zp  <- .bt_compute_zero_proportion_around_end(self, rwe)
    valid_gap <- .bt_datetime_gap_check(self, rwe)
    med <- .bt_compute_refinement_window_median(self, refinement_window_start, rwe)
    metric <- .bt_compute_metric(self, med)
    ending <- .np_slice(med, length(med) - amaw, length(med))
    above  <- 1 - .below_prop(ending, metric)
  }
  rwe
}

#' bridge_gap_validation: when the end search stops on an invalid gap, jump past
#' it and search for a fresh window end. Returns list(start, end).
#' @noRd
.bt_bridge_gap_validation <- function(self, refinement_window_end) {
  rws <- refinement_window_end + 1L
  ic  <- rws + 1L
  rwe <- .bt_compute_refinement_window_end(self, ic, rws)
  list(start = rws, end = rwe)
}

#' compute_improved_refinement_window_start: with a full window known, pull the
#' start back while the leading median stays below the metric. Uses
#' self$previous_transition.
#' @noRd
.bt_compute_improved_refinement_window_start <- function(self, refinement_window_start, refinement_window_end) {
  rws  <- refinement_window_start
  amaw <- self$activity_median_analysis_window

  valid_gap <- .bt_datetime_gap_check(self, rws)
  med <- .bt_compute_refinement_window_median(self, rws, refinement_window_end)
  metric <- .bt_compute_metric(self, med)
  ending <- .np_slice(med, 0L, amaw)
  below  <- .below_prop(ending, metric)

  while ((rws > 0L) && (rws - 1L > self$previous_transition) && valid_gap && (below > 0)) {
    rws <- rws - 1L
    if (rws > 0L) valid_gap <- .bt_datetime_gap_check(self, rws)
    med <- .bt_compute_refinement_window_median(self, rws, refinement_window_end)
    metric <- .bt_compute_metric(self, med)
    ending <- .np_slice(med, 0L, amaw)
    below  <- .below_prop(ending, metric)
  }
  rws
}

# ── Stage 3: bed-time refiner — Unit B (peak/valley filtering + candidates) ─────
# Operate on a peaks_and_valleys data.frame (from .identify_peaks_and_valleys)
# and refinement state on `self` (metric, refinement_window_start/end,
# refinement_window_activity, refinement_window_levels). Region indices in the
# table follow the Python 0-indexed convention.

#' remove_after_long_valley: drop everything after the first "long" valley
#' (a valley at least bedtime_high_probability_sleep_valley_length long).
#' @noRd
.bt_remove_after_long_valley <- function(self, pv) {
  thr <- self$bedtime_high_probability_sleep_valley_length
  if (isTRUE(self$do_remove_after_long_valley)) {
    invalid <- which(pv$class == "v" & !(pv$length < thr))   # 1-indexed rows
    if (length(invalid) > 0) pv <- pv[seq_len(invalid[1]), , drop = FALSE]
  }
  rownames(pv) <- NULL
  pv
}

#' remove_before_long_peak: drop everything before the last "long" peak
#' (a peak at least bedtime_high_probability_awake_peak_length long), unless it
#' is already the final region.
#' @noRd
.bt_remove_before_long_peak <- function(self, pv) {
  cnt <- nrow(pv)
  thr <- self$bedtime_high_probability_awake_peak_length
  if (isTRUE(self$do_remove_before_long_peak)) {
    invalid <- which(pv$class == "p" & !(pv$length < thr))
    if (length(invalid) > 0) {
      last <- invalid[length(invalid)]
      if (last < cnt) pv <- pv[last:cnt, , drop = FALSE]
    }
  }
  rownames(pv) <- NULL
  pv
}

#' remove_before_tall_peak: drop everything before the last "tall" peak
#' (mean > 10 * metric), only when metric > 5.
#' @noRd
.bt_remove_before_tall_peak <- function(self, pv) {
  cnt <- nrow(pv)
  metric <- self$metric
  if (isTRUE(self$do_remove_before_tall_peak)) {
    invalid <- if (metric > 5) which(pv$class == "p" & !(pv$mean <= 10 * metric)) else integer(0)
    if (length(invalid) > 0) {
      last <- invalid[length(invalid)]
      pv <- pv[last:cnt, , drop = FALSE]
    }
  }
  rownames(pv) <- NULL
  pv
}

#' filter_peaks_and_valleys: sequential feature-based merging of short / weak
#' peaks and valleys (the main bed-time cleaning loop). region_index is
#' 0-indexed (Python); R rows are region_index + 1.
#' @noRd
.bt_filter_peaks_and_valleys <- function(self, pv) {
  bmvl   <- self$bedtime_minimum_valley_length
  bmpl   <- self$bedtime_minimum_peak_length
  cond   <- self$condition
  rwa    <- self$refinement_window_activity
  metric <- self$metric

  cnt <- nrow(pv)
  rownames(pv) <- NULL
  if (!(cnt > 2 && cnt != 3)) return(pv)

  # First region, if a short valley, merges into its successor.
  if (pv$class[1] == "v" && pv$length[1] < bmvl) {
    pv  <- .remove_peak_valley(pv, 0L, rwa, metric)
    cnt <- nrow(pv)
  }

  ri <- 1L                                   # 0-indexed region_index
  while (ri < cnt && (cnt > 2 && cnt != 3)) {
    remove <- FALSE
    r <- ri + 1L                             # R row index
    if (pv$class[r] == "p") {
      if (ri < cnt - 2L) {
        if (pv$length[r] < bmpl) {
          if (pv$length[r + 1L] < bmvl) {
            # short peak right before a short valley is spared
          } else {
            remove <- TRUE
          }
        } else {
          if (pv$mean[r] <= 1.33 * metric) {
            remove <- TRUE
          } else if (cond == 2 && ri > 0 && pv$length[r - 1L] >= 15 * pv$length[r]) {
            remove <- TRUE
          }
        }
      } else {
        if (pv$length[r - 1L] >= 30 &&
            (pv$zero_proportion[r - 1L] > 2/3 ||
             pv$above_threshold_proportion[r - 1L] < 0.1)) {
          remove <- TRUE
        }
      }
    } else {                                 # valley
      if (ri < cnt - 1L) {
        if (pv$length[r] < bmvl) {
          if (ri > 1L) {
            remove <- TRUE
          } else {
            if (pv$length[r + 1L] > bmpl && pv$zero_proportion[r + 1L] < 1/3) {
              remove <- TRUE
            }
          }
        } else {
          remove_points <- 0
          if (pv$above_threshold_proportion[r] >= 0.33) remove_points <- remove_points + 1
          if (pv$zero_proportion[r] < 0.45) {
            if (pv$zero_proportion[r] > 0.1) remove_points <- remove_points + 1
            else                             remove_points <- remove_points + 1.5
          }
          if (pv$mean[r] >= 0.66 * metric) remove_points <- remove_points + 0.5
          if ((pv$length[r] / length(rwa) >= 0.3) ||
              (ri > 0 && pv$length[r] >= 1.5 * pv$length[r - 1L])) {
            remove_points <- remove_points - 1
          }
          if (remove_points > 1.5) remove <- TRUE
        }
      }
    }

    if (remove) {
      pv  <- .remove_peak_valley(pv, ri, rwa, metric)
      cnt <- nrow(pv)
    } else {
      ri <- ri + 1L
    }
  }
  pv
}

#' identify_bedtime_candidates: valley starts are candidate bed-times, with a
#' relative-depth test for interior valleys. Mutates self$refinement_window_levels.
#' Returns an integer vector of window-relative candidate starts.
#' @noRd
.bt_identify_bedtime_candidates <- function(self, pv) {
  cnt <- nrow(pv)
  rws <- self$refinement_window_start
  cands <- integer(0)
  for (ri in seq_len(cnt)) {
    r0    <- ri - 1L                         # 0-indexed region_index
    start <- as.integer(rws + pv$start[ri])
    end   <- as.integer(rws + pv$end[ri])
    if (end > start) self$refinement_window_levels[(start + 1L):end] <- pv$mean[ri]
    if (pv$class[ri] == "v") {
      if (r0 > 1L && r0 < cnt - 1L) {
        if (pv$mean[ri] < 0.5 * pv$mean[ri - 1L]) cands <- c(cands, as.integer(pv$start[ri]))
      } else {
        cands <- c(cands, as.integer(pv$start[ri]))
      }
    }
  }
  if (length(cands) == 0L) cands <- as.integer(pv$start[pv$class == "v"])
  cands
}

#' bedtime_candidates_crossings_filter: drop candidates before the last upward
#' crossing of the short-window median past the metric. Returns
#' list(candidates, down) where `down` are the downward-crossing positions
#' (window-relative, 0-indexed) used later in scoring.
#' @noRd
.bt_bedtime_candidates_crossings_filter <- function(self, candidates) {
  swam   <- self$short_window_activity_median
  rws    <- self$refinement_window_start
  rwe    <- self$refinement_window_end
  metric <- self$metric

  seg     <- swam[(rws + 1L):(rwe + 1L)]                  # swam[rws:rwe+1]
  cross01 <- as.integer(seg >= metric)
  d       <- diff(c(0L, cross01))                         # diff(concat([0], cross))
  up      <- which(d > 0) - 1L                            # 0-indexed positions
  down    <- which(d < 0) - 1L

  cands <- candidates
  if (length(up) > 0) {
    last_up <- up[length(up)]
    cands   <- cands[cands >= last_up]
  }
  list(candidates = cands, down = down)
}

# ── Stage 3: bed-time refiner — Unit C (scoring + refine orchestrator) ─────────

#' choose_best_bedtime_candidate: score each candidate and pick the refined
#' bed-time (port of CSPD_BedTime_Refiner.choose_best_bedtime_candidate).
#'
#' Scoring, in order: an optional bonus for the last candidate; a bonus for the
#' candidate nearest the last downward metric crossing; after sorting by the
#' activity-median-difference valley depth (ascending), a bonus for the sharpest
#' drop; then, for condition 0/2, a forward gap / epochs-above-metric analysis
#' that either rewards the fewest-epochs candidate or applies a graded
#' "thresholded" bonus. The final pick is the highest score (candidate index
#' breaks ties), with an optional second-best override.
#'
#' pandas `sort_values` tie-breaks are reproduced with R's stable `order()`;
#' this was verified bit-exact against the real refiner across all transitions,
#' including amd-tie cases. `candidates` and `down` are window-relative,
#' 0-indexed (from .bt_identify_bedtime_candidates / crossings filter).
#' @noRd
.bt_choose_best_bedtime_candidate <- function(self, candidates, down) {
  n <- length(candidates)
  if (n == 0L) return(as.integer(self$initial_transition_candidate))

  rws      <- self$refinement_window_start
  smoothed <- self$refinement_window_activity_median_difference_smoothed

  amd <- vapply(candidates, function(c0) .get_peak(smoothed, c0, valley = TRUE), numeric(1))
  sc <- data.frame(
    candidate   = as.integer(candidates),
    amd         = as.numeric(amd),
    epochs      = rep(self$after_candidate_window, n),
    thresholded = rep(FALSE, n),
    gap_after   = rep(TRUE, n),
    score       = rep(0, n),
    stringsAsFactors = FALSE
  )

  if (isTRUE(self$bedtime_score_last_candidate)) {
    sc$score[n] <- sc$score[n] + self$bedtime_last_candidate_score
  }
  if (length(down) > 0) {
    last_down <- down[length(down)]
    bcdc <- which.min(abs(as.integer(candidates) - last_down))   # numpy argmin = first min
    sc$score[bcdc] <- sc$score[bcdc] + self$bedtime_best_crossing_distance_candidate_score
  }

  sc <- sc[order(sc$amd), , drop = FALSE]                        # stable ascending
  sc$score[1] <- sc$score[1] + self$bedtime_best_median_difference_candidate_score

  if (self$condition %in% c(0, 2)) {
    acw    <- self$after_candidate_window
    maxep  <- self$bedtime_maximum_epochs_above_metric_after_candidate
    metric <- self$metric
    for (i in seq_len(n)) {
      cand <- rws + sc$candidate[i]
      end  <- cand + acw
      if (end > self$data_length) end <- self$data_length
      valid <- TRUE
      g <- cand
      while (valid && (g + 1L < end)) {
        valid <- .bt_datetime_gap_check(self, g, direction = "forward")
        g <- g + 1L
      }
      if (valid) {
        sc$gap_after[i] <- FALSE
        ea <- sum(self$activity[(cand + 1L):end] >= metric)
        sc$epochs[i] <- ea
        if (ea <= maxep) sc$thresholded[i] <- TRUE
      }
    }

    if (sum(!sc$gap_after) > 0) {
      if (sum(sc$thresholded) == 0) {
        sc <- sc[order(sc$gap_after, sc$epochs), , drop = FALSE]  # FALSE first, epochs asc
        sc$score[1] <- sc$score[1] + self$bedtime_best_epochs_above_metric_after_score
      } else {
        factor <- if (maxep > 0) {
          -self$bedtime_thresholded_candidate_score_amplitude / maxep
        } else {
          -self$bedtime_thresholded_candidate_score_amplitude / 1e-3
        }
        thr <- which(sc$thresholded)
        sc$score[thr] <- sc$score[thr] +
          factor * (sc$epochs[thr] - maxep) + self$bedtime_thresholded_candidate_score_minimum
      }
    }

    sc <- sc[order(-sc$score, -sc$candidate), , drop = FALSE]     # score desc, candidate desc

    if (isTRUE(self$consider_second_best_candidate) && n > 1L) {
      if ((sc$candidate[1] - sc$candidate[2] > 60) &&
          (sc$thresholded[2] && !sc$thresholded[1])) {
        refined <- rws + as.integer(sc$candidate[2])
      } else {
        refined <- rws + as.integer(sc$candidate[1])
      }
    } else {
      refined <- rws + as.integer(sc$candidate[1])
    }
  } else {
    refined <- rws + as.integer(candidates[n])
  }
  as.integer(refined)
}

#' refine: full bed-time refinement for one transition (port of
#' CSPD_BedTime_Refiner.refine). Discovers the refinement window (Unit A),
#' cleans the peaks/valleys (Unit B), then scores candidates (Unit C).
#'
#' Returns a list mirroring the Python tuple: refined_bedtime,
#' refinement_window_start, refinement_window_end, refinement_window_activity_median,
#' refinement_window_activity_median_difference_smoothed, refinement_window_levels,
#' metric. All indices are 0-indexed. Validated bit-exact (refined index, window
#' start, window end) against the real refiner across all 52 transitions.
#' @noRd
.bt_refine <- function(self, refinement_window_levels, initial_transition_candidate,
                       previous_transition, next_transition) {
  self$refinement_window_levels     <- refinement_window_levels
  self$initial_transition_candidate <- initial_transition_candidate
  self$previous_transition          <- previous_transition
  self$next_transition              <- next_transition

  rws <- .bt_compute_initial_refinement_window_start(self, initial_transition_candidate - 1L)
  rwe <- .bt_compute_refinement_window_end(self, initial_transition_candidate + 1L, rws)

  if (!.bt_datetime_gap_check(self, rwe)) {
    br  <- .bt_bridge_gap_validation(self, rwe)
    rws <- br$start; rwe <- br$end
  } else {
    rws <- .bt_compute_improved_refinement_window_start(self, rws, rwe)
  }

  # In-window datetime-gap trim: keep the larger side of the first big gap.
  rwdd   <- .datetime_diff(self$secs[(rws + 1L):(rwe + 1L)])
  rwgaps <- as.integer(rwdd >= self$maximum_allowed_gap)
  if (sum(rwgaps) > 0) {
    fg <- which(rwgaps == 1L)[1] - 1L                  # 0-indexed first gap location
    if (fg > length(rwdd) - fg) {
      rwe <- rws + fg - 1L
    } else {
      rws <- rws + fg + 1L
    }
  }
  if (rws >= rwe) {
    rwe <- rws + 60L
    if (rwe >= self$data_length) rwe <- self$data_length - 1L
  }

  self$refinement_window_start <- rws
  self$refinement_window_end   <- rwe
  self$refinement_window_activity_median <- .bt_compute_refinement_window_median(self, rws, rwe)
  self$metric <- .bt_compute_metric(self, self$refinement_window_activity_median)
  self$refinement_window_activity_median_difference <- diff(self$refinement_window_activity_median)
  self$refinement_window_activity_median_difference_smoothed <-
    .cspd_median_filter(self$refinement_window_activity_median_difference, self$median_filter_half_window_size)
  self$refinement_window_activity <- self$activity[(rws + 1L):(rwe + 1L)]

  pv  <- .identify_peaks_and_valleys(signal = self$refinement_window_activity_median,
                                     activity = self$refinement_window_activity,
                                     threshold = self$metric)
  cnt <- nrow(pv)
  for (ri in seq_len(cnt)) {
    s <- as.integer(rws + pv$start[ri]); e <- as.integer(rws + pv$end[ri])
    if (e > s) self$refinement_window_levels[(s + 1L):e] <- pv$mean[ri]
  }

  pv <- .bt_remove_after_long_valley(self, pv)
  pv <- .bt_remove_before_long_peak(self, pv)
  pv <- .bt_remove_before_tall_peak(self, pv)
  self$peaks_and_valleys <- .bt_filter_peaks_and_valleys(self, pv)
  cnt <- nrow(self$peaks_and_valleys)

  new_rws <- self$refinement_window_start + self$peaks_and_valleys$start[1]
  self$refinement_window_end <- self$refinement_window_start + self$peaks_and_valleys$end[cnt]
  if (self$refinement_window_end >= self$data_length) self$refinement_window_end <- self$data_length - 1L
  self$refinement_window_start <- new_rws

  self$refinement_window_activity_median <-
    .bt_compute_refinement_window_median(self, self$refinement_window_start, self$refinement_window_end)
  self$metric <- .bt_compute_metric(self, self$refinement_window_activity_median)
  self$refinement_window_activity_median_difference <- diff(self$refinement_window_activity_median)
  self$refinement_window_activity_median_difference_smoothed <-
    .cspd_median_filter(self$refinement_window_activity_median_difference, self$median_filter_half_window_size)
  self$refinement_window_activity <-
    self$activity[(self$refinement_window_start + 1L):(self$refinement_window_end + 1L)]

  if (isTRUE(self$update_peaks_and_valleys)) {
    self$peaks_and_valleys <- .identify_peaks_and_valleys(signal = self$refinement_window_activity_median,
                                                          activity = self$refinement_window_activity,
                                                          threshold = self$metric)
  } else {
    self$peaks_and_valleys$start[1] <- 0L
    self$peaks_and_valleys$end[1]   <- self$peaks_and_valleys$start[1] + self$peaks_and_valleys$length[1]
    if (cnt >= 2L) {
      for (ii in 2:cnt) {
        self$peaks_and_valleys$start[ii] <- self$peaks_and_valleys$end[ii - 1L]
        self$peaks_and_valleys$end[ii]   <- self$peaks_and_valleys$start[ii] + self$peaks_and_valleys$length[ii]
      }
    }
  }

  cands <- .bt_identify_bedtime_candidates(self, self$peaks_and_valleys)
  cf    <- .bt_bedtime_candidates_crossings_filter(self, cands)
  if (isTRUE(self$do_bedtime_candidates_crossings_filter)) cands <- cf$candidates

  refined_bedtime <- .bt_choose_best_bedtime_candidate(self, cands, cf$down)

  list(
    refined_bedtime                        = refined_bedtime,
    refinement_window_start                = self$refinement_window_start,
    refinement_window_end                  = self$refinement_window_end,
    refinement_window_activity_median      = self$refinement_window_activity_median,
    refinement_window_activity_median_difference_smoothed =
      self$refinement_window_activity_median_difference_smoothed,
    refinement_window_levels               = self$refinement_window_levels,
    metric                                 = self$metric
  )
}

# ── Stage 4: minimum-length filter + CSPD.model wiring ───────────────────────
# Ports the tail of CSPD.model: the boolean length filter (with its faithful
# second-pass quirk) and the transition-refinement driver that threads the
# bed-time / get-up-time refiners over every MSP transition to build the final
# `refined_output` (1 = wake, 0 = sleep PERIOD). Indices follow the Python
# convention (0-indexed, exclusive ends); convert at vector access.

#' boolean_length_filter (port of cspd_functions_without_prints.boolean_length_filter).
#'
#' Removes runs of `class_to_filter` ("v" = sleep, "p" = wake) shorter than
#' `minimum_length` by merging them into their neighbours, then rebuilds the
#' 0/1 signal. The activity basis is all-zeros (as in Python: a boolean signal
#' carries no activity features), so only `length`/`class` drive the merges.
#'
#' NOTE: the upstream second filtering pass contains a bug — an `it = 0` that is
#' never used, and a `while (t < count)` that reuses the *stale* `t` left over
#' from the signal-rebuild loop instead of resetting it to 0. In almost all
#' cases that stale `t` is already >= the new region count, so the second pass
#' is a no-op. This is replicated faithfully (the rebuild that follows it is a
#' fresh `for` loop, so it always runs).
#' @noRd
.boolean_length_filter <- function(minimum_length, signal, class_to_filter = "v") {
  n     <- length(signal)
  act0  <- numeric(n)                # zero activity: only length/class are used
  thr   <- 0.5
  pv    <- .identify_peaks_and_valleys(signal = signal, activity = act0, threshold = thr)
  count <- nrow(pv)
  filtered_signal <- rep(1, n)

  if (count > 1L) {
    # First pass: merge short class_to_filter regions into their neighbours.
    t <- 0L
    while (t < count) {
      remove <- (pv$length[t + 1L] <= minimum_length) &&
                (pv$class[t + 1L] == class_to_filter)
      if (remove) {
        pv <- .remove_peak_valley(pv, t, act0, thr); count <- nrow(pv)
      } else {
        t <- t + 1L
      }
    }

    # Rebuild the signal from the surviving valleys (t ends == count).
    count <- nrow(pv)
    t <- 0L
    while (t < count) {
      if (pv$class[t + 1L] == "v")
        filtered_signal[(pv$start[t + 1L] + 1L):pv$end[t + 1L]] <- 0
      t <- t + 1L
    }

    # Second pass — re-identify and (faithfully) reuse the stale `t`.
    pv    <- .identify_peaks_and_valleys(signal = filtered_signal, activity = act0, threshold = thr)
    count <- nrow(pv)
    if (count > 1L) {
      # `it <- 0L` in the source is intentionally unused; `t` is NOT reset here.
      while (t < count) {
        remove <- (pv$length[t + 1L] <= minimum_length) &&
                  (pv$class[t + 1L] == class_to_filter)
        if (remove) {
          pv <- .remove_peak_valley(pv, t, act0, thr); count <- nrow(pv)
        } else {
          t <- t + 1L
        }
      }
      filtered_signal <- rep(1, n)
      count <- nrow(pv)
      for (t in seq_len(count)) {
        if (pv$class[t] == "v")
          filtered_signal[(pv$start[t] + 1L):pv$end[t]] <- 0
      }
    }
  } else if (count == 1L) {
    filtered_signal <- as.numeric(signal)
  }

  filtered_signal
}

#' ActTrust CSPD configuration (the cspd_wrapper param_set + CSPD defaults).
#'
#' This is the production parameter set `cspd_wrapper` uses to generate the
#' validation fixtures. Values are pre-conversion; the minutes->epochs and score
#' unpacking happen inside `.cspd_refine_periods`, mirroring CSPD.model.
#' @noRd
.cspd_acttrust_params <- function() {
  list(
    # bedtime scores  = [last, best_median, best_crossing, best_epochs, thr_amp, thr_min]
    bedtime_scores   = c(0.4081634489795918, 0.551020306122449, 0.7142852857142857,
                         0.6122446734693877, 0.48979593877551014, 0.44897969387755093),
    # getuptime scores = [first, best_median, best_crossing, best_epochs, thr_amp, thr_min]
    getuptime_scores = c(0.26530659183673466, 0.48979593877551014, 0.795917775510204,
                         0.24489846938775509, 0.999999, 0.24489846938775509),
    length_thresholds    = c(8, 20, 17, 16),   # bt peak, bt valley, gt peak, gt valley
    candidate_thresholds = c(0.4328571428571429, 0.37244897959183676, 0.5619047619047619),
    short_window_activity_median_threshold_quantile = 0.43163269387755104,

    peak_valley_minimum_length                       = 11,
    median_filter_short_window                       = 41,
    after_candidate_window                           = 47,
    half_window_around_border                        = 28,
    median_filter_half_window_size                   = 7,
    short_window_activity_median_minimum_high_epochs = 9,
    activity_median_analysis_window                  = 9,
    minimum_short_window_activity_median_threshold   = 1.0,
    sleep_minimum_length                             = 120,
    nap_minimum_length                               = 20,
    refinement_maximum_allowed_datetime_gap          = 10,   # minutes; *60 -> seconds
    sleep_maximum_allowed_datetime_gap               = 60,   # minutes; *60 -> seconds
    do_peak_valley_length_filter                     = TRUE,
    detect_naps                                      = FALSE,

    bedtime_metric_method                            = 1L,
    bedtime_metric_parameter                         = 0.42755102040816323,
    bedtime_do_remove_before_long_peak               = TRUE,
    bedtime_do_remove_before_tall_peak               = FALSE,
    bedtime_do_remove_after_long_valley              = FALSE,
    bedtime_update_peaks_and_valleys                 = FALSE,
    bedtime_do_bedtime_candidates_crossings_filter   = FALSE,
    bedtime_consider_second_best_candidate           = TRUE,
    bedtime_score_last_candidate                     = TRUE,
    bedtime_high_probability_awake_peak_length       = 45,
    bedtime_high_probability_sleep_valley_length     = 45,

    getuptime_metric_method                          = 1L,
    getuptime_metric_parameter                       = 0.4071428571428571,
    getuptime_do_remove_after_long_tall_peak         = TRUE,
    getuptime_do_remove_before_long_valley           = FALSE,
    getuptime_update_peaks_and_valleys               = FALSE,
    getuptime_score_first_candidate                  = TRUE,
    getuptime_high_probability_awake_peak_length     = 45,
    getuptime_high_probability_sleep_valley_length   = 45
  )
}

#' Nap-mode CSPD configuration (the nap_wrapper param_set + CSPD defaults).
#'
#' Mirrors `nap_wrapper`'s `CSPD(detect_naps=True, ...)` construction. Values
#' are pre-conversion (minutes->epochs and score unpacking happen inside
#' `.cspd_refine_periods`). The trailing `nap_*` fields configure the nap-mode
#' MSP (`.crespo_nap_msp`) rather than the refiner.
#' @noRd
.cspd_nap_params <- function() {
  list(
    bedtime_scores   = c(0.21052689, 0.84210458, 0.10526395,
                         0.21052689, 0.21052689, 0.57894721),
    getuptime_scores = c(0.63157868, 0.63157868, 0.57894721,
                         0.68421016, 0.89473605, 0.68421016),
    length_thresholds    = c(12, 16, 7, 12),
    candidate_thresholds = c(0.67210526, 0.43842105, 0.58421053),
    short_window_activity_median_threshold_quantile = 0.45,

    peak_valley_minimum_length                       = 4,
    median_filter_short_window                       = 20,
    after_candidate_window                           = 15,
    half_window_around_border                        = 60,   # CSPD default
    median_filter_half_window_size                   = 4,    # CSPD default
    short_window_activity_median_minimum_high_epochs = 3,    # CSPD default
    activity_median_analysis_window                  = 20,
    minimum_short_window_activity_median_threshold   = 1.0,
    sleep_minimum_length                             = 120,  # unused (nap path)
    nap_minimum_length                               = 20,
    refinement_maximum_allowed_datetime_gap          = 10,
    sleep_maximum_allowed_datetime_gap               = 60,
    do_peak_valley_length_filter                     = TRUE,
    detect_naps                                      = TRUE,

    bedtime_metric_method                            = 1L,
    bedtime_metric_parameter                         = 0.50263158,
    bedtime_do_remove_before_long_peak               = TRUE,
    bedtime_do_remove_before_tall_peak               = TRUE,
    bedtime_do_remove_after_long_valley              = FALSE,
    bedtime_update_peaks_and_valleys                 = FALSE,
    bedtime_do_bedtime_candidates_crossings_filter   = FALSE,
    bedtime_consider_second_best_candidate           = FALSE,
    bedtime_score_last_candidate                     = FALSE,
    bedtime_high_probability_awake_peak_length       = 20,
    bedtime_high_probability_sleep_valley_length     = 20,

    getuptime_metric_method                          = 1L,
    getuptime_metric_parameter                       = 0.41578947,
    getuptime_do_remove_after_long_tall_peak         = TRUE,
    getuptime_do_remove_before_long_valley           = FALSE,
    getuptime_update_peaks_and_valleys               = FALSE,
    getuptime_score_first_candidate                  = TRUE,
    getuptime_high_probability_awake_peak_length     = 20,
    getuptime_high_probability_sleep_valley_length   = 20,

    # ── Nap-mode MSP params (consumed by .crespo_nap_msp, not the refiner) ──
    nap_median_filter_h                    = 0.97894737,
    nap_pad_h                              = 0.05263247,
    nap_median_activity_threshold          = 3.15789563,
    nap_zero_proportion_threshold          = 0.22105263,
    nap_zero_proportion_filter_window_size = 5L,
    compute_output_naps_with_logical_and   = TRUE
  )
}

#' Peak/valley zero-proportion filter (port of peak_valley_zero_proportion_filter).
#'
#' Faithful port of cspd_functions.peak_valley_zero_proportion_filter: removes
#' regions of `class_to_filter` whose zero-proportion is `<= minimum_zero_
#' proportion`. As in the reference, features are computed against a zeros
#' activity array, so every region's zero-proportion is 1.0 and (for the nap
#' threshold 1/3) nothing is removed — the filter is an identity here, ported
#' for faithfulness. Unlike `.boolean_length_filter`, the second pass resets
#' the loop index (the reference does not carry the stale `t`).
#' @noRd
.peak_valley_zero_proportion_filter <- function(minimum_zero_proportion, signal,
                                                class_to_filter = "v") {
  n    <- length(signal)
  act0 <- numeric(n)        # zeros (matches reference: activity = np.zeros)
  thr  <- 0.5
  filtered_signal <- rep(1, n)

  pv    <- .identify_peaks_and_valleys(signal = signal, activity = act0, threshold = thr)
  count <- nrow(pv)

  if (count > 1L) {
    t <- 0L
    while (t < count) {
      remove <- (pv$zero_proportion[t + 1L] <= minimum_zero_proportion) &&
                (pv$class[t + 1L] == class_to_filter)
      if (remove) {
        pv <- .remove_peak_valley(pv, t, act0, thr); count <- nrow(pv)
      } else {
        t <- t + 1L
      }
    }

    count <- nrow(pv); t <- 0L
    while (t < count) {
      if (pv$class[t + 1L] == "v")
        filtered_signal[(pv$start[t + 1L] + 1L):pv$end[t + 1L]] <- 0
      t <- t + 1L
    }

    # Second pass — re-identify; the reference resets t = 0 here.
    pv    <- .identify_peaks_and_valleys(signal = filtered_signal, activity = act0, threshold = thr)
    count <- nrow(pv)
    if (count > 1L) {
      t <- 0L
      while (t < count) {
        remove <- (pv$zero_proportion[t + 1L] <= minimum_zero_proportion) &&
                  (pv$class[t + 1L] == class_to_filter)
        if (remove) {
          pv <- .remove_peak_valley(pv, t, act0, thr); count <- nrow(pv)
        } else {
          t <- t + 1L
        }
      }
      filtered_signal <- rep(1, n)
      count <- nrow(pv)
      for (t in seq_len(count)) {
        if (pv$class[t] == "v")
          filtered_signal[(pv$start[t] + 1L):pv$end[t]] <- 0
      }
    }
  } else if (count == 1L) {
    filtered_signal <- as.numeric(signal)
  }

  filtered_signal
}

#' Refine an MSP detection into final sleep periods (port of the CSPD.model tail).
#'
#' Given the (on-wrist) activity, timestamps and the MSP sleep detection from
#' `detect_sleep_crespo`, this reproduces CSPD.model from the parameter
#' derivation onward: scale activity, derive the epoch-converted parameters,
#' apply the stage-1 peak-valley length filter and stage-2 sleep-gap
#' separation, build the transition list, refine every transition with the
#' bed-time / get-up-time refiners (threading the shared
#' `refinement_window_levels`), apply the stage-4 minimum-length filter, and
#' snap the recording borders with the Crespo heuristic.
#'
#' @param activity numeric raw activity (PIM), on-wrist subset.
#' @param datetime_stamps POSIXct timestamps, same length as `activity`.
#' @param msp_detection integer 0/1 vector (1 = wake, 0 = sleep) — the MSP
#'   detection entering the refiner (Python `final_sleep_detection`).
#' @param condition integer condition flag from MSP (default 0).
#' @param params CSPD configuration list (default `.cspd_acttrust_params()`).
#' @return list(refined_output, refined_sleep_df, refined_bedtimes,
#'   refined_getuptimes, transitions). `refined_output` is 1 = wake, 0 = sleep
#'   PERIOD (post boolean-length-filter and post border heuristic).
#' @noRd
.cspd_refine_periods <- function(activity, datetime_stamps, msp_detection,
                                 condition = 0L, params = .cspd_acttrust_params()) {
  p           <- params
  data_length <- length(activity)
  activity    <- as.numeric(activity)

  # ── Epoch duration (mode of the datetime differences) + scaled activity ──
  dd_all   <- .datetime_diff(datetime_stamps)
  tb       <- sort(table(dd_all), decreasing = TRUE)
  duration <- as.numeric(names(tb)[1])
  scaled_activity <- (1 / duration) * activity

  # ── Score / length / candidate unpacking ──
  b  <- p$bedtime_scores
  g  <- p$getuptime_scores
  l  <- as.integer(round(p$length_thresholds))
  cc <- p$candidate_thresholds

  bedtime_minimum_peak_length     <- l[1]
  bedtime_minimum_valley_length   <- l[2]
  getuptime_minimum_peak_length   <- l[3]
  getuptime_minimum_valley_length <- l[4]

  bedtime_max_epochs        <- cc[1]
  getuptime_max_epochs      <- cc[2]
  zero_proportion_threshold <- cc[3]

  swam_quantile <- p$short_window_activity_median_threshold_quantile
  if (condition == 2L) {
    zero_proportion_threshold <- zero_proportion_threshold * 0.666
    swam_quantile             <- 1.333 * swam_quantile
  }

  refinement_gap <- p$refinement_maximum_allowed_datetime_gap * 60
  sleep_gap      <- p$sleep_maximum_allowed_datetime_gap * 60

  # ── Minutes -> epochs (round) and truncating int() conversions ──
  conv_round <- function(x) as.integer(round(x * 60 / duration))
  conv_int   <- function(x) as.integer(x * 60 / duration)

  bedtime_minimum_peak_length     <- conv_round(bedtime_minimum_peak_length)
  bedtime_minimum_valley_length   <- conv_round(bedtime_minimum_valley_length)
  getuptime_minimum_peak_length   <- conv_round(getuptime_minimum_peak_length)
  getuptime_minimum_valley_length <- conv_round(getuptime_minimum_valley_length)
  after_candidate_window          <- conv_round(p$after_candidate_window)
  # NB: the *original* after_candidate_window (minutes) scales the max-epochs
  bedtime_max_epochs   <- as.integer(round(p$after_candidate_window * bedtime_max_epochs))
  getuptime_max_epochs <- as.integer(round(p$after_candidate_window * getuptime_max_epochs))
  half_window_around_border       <- conv_round(p$half_window_around_border)
  activity_median_analysis_window <- conv_round(p$activity_median_analysis_window)
  median_filter_half_window_size  <- conv_round(p$median_filter_half_window_size)
  peak_valley_minimum_length      <- conv_round(p$peak_valley_minimum_length)
  median_filter_short_window      <- conv_int(p$median_filter_short_window)
  swamm_high_epochs               <- conv_int(p$short_window_activity_median_minimum_high_epochs)

  # ── Auxiliary short-window activity median + threshold (functions.py median_filter) ──
  swam <- .cspd_median_filter(scaled_activity, median_filter_short_window)   # padding = NULL
  swam_threshold <- as.numeric(stats::quantile(swam, swam_quantile, names = FALSE, type = 7))

  # ── Stage 1: peak-valley length filter ──
  detection <- as.numeric(msp_detection)
  if (isTRUE(p$do_peak_valley_length_filter)) {
    detection <- .peak_valley_length_filter(detection, peak_valley_minimum_length)
  }

  # ── Stage 2: sleep-gap separation ──
  detection <- .sleep_gap_separation(detection, .datetime_diff(datetime_stamps), sleep_gap)

  # ── Build transitions: diff -> nonzero, prepend a synthetic bedtime if the
  #    recording starts asleep. Reproduces the source faithfully, including a
  #    possible duplicate index-0 entry. ──
  borders <- diff(detection)
  idx0    <- which(borders != 0) - 1L                 # 0-indexed positions in the diff
  if (detection[1] == 0) {
    borders[1] <- -1
    idx0       <- c(0L, idx0)
  }
  trans_i   <- idx0
  trans_dir <- borders[idx0 + 1L]
  num       <- length(trans_i)
  refined_i <- trans_i                                # refined_transitions[*][0], mutable

  refinement_window_levels <- numeric(data_length)
  refined_bedtimes   <- integer(0)
  refined_getuptimes <- integer(0)
  refined_output     <- numeric(data_length)          # zeros (= all sleep)

  if (num > 0L) {
    bt <- .cspd_bedtime_refiner(
      activity = scaled_activity, datetime_stamps = datetime_stamps,
      short_window_activity_median = swam,
      minimum_short_window_activity_median_threshold = p$minimum_short_window_activity_median_threshold,
      short_window_activity_median_minimum_high_epochs = swamm_high_epochs,
      half_window_around_border = half_window_around_border,
      activity_median_analysis_window = activity_median_analysis_window,
      maximum_allowed_gap = refinement_gap,
      quantile_threshold = swam_threshold,
      median_filter_half_window_size = median_filter_half_window_size,
      metric_method = p$bedtime_metric_method,
      metric_parameter = p$bedtime_metric_parameter,
      condition = condition,
      bedtime_score_last_candidate = p$bedtime_score_last_candidate,
      bedtime_last_candidate_score = b[1],
      bedtime_best_median_difference_candidate_score = b[2],
      bedtime_best_crossing_distance_candidate_score = b[3],
      bedtime_best_epochs_above_metric_after_score = b[4],
      bedtime_thresholded_candidate_score_amplitude = b[5],
      bedtime_thresholded_candidate_score_minimum = b[6],
      bedtime_minimum_peak_length = bedtime_minimum_peak_length,
      bedtime_minimum_valley_length = bedtime_minimum_valley_length,
      after_candidate_window = after_candidate_window,
      bedtime_maximum_epochs_above_metric_after_candidate = bedtime_max_epochs,
      zero_proportion_threshold = zero_proportion_threshold,
      do_remove_after_long_valley = p$bedtime_do_remove_after_long_valley,
      do_remove_before_long_peak = p$bedtime_do_remove_before_long_peak,
      do_remove_before_tall_peak = p$bedtime_do_remove_before_tall_peak,
      update_peaks_and_valleys = p$bedtime_update_peaks_and_valleys,
      do_bedtime_candidates_crossings_filter = p$bedtime_do_bedtime_candidates_crossings_filter,
      consider_second_best_candidate = p$bedtime_consider_second_best_candidate,
      bedtime_high_probability_awake_peak_length = p$bedtime_high_probability_awake_peak_length,
      bedtime_high_probability_sleep_valley_length = p$bedtime_high_probability_sleep_valley_length
    )

    gt <- .cspd_getuptime_refiner(
      activity = scaled_activity, datetime_stamps = datetime_stamps,
      short_window_activity_median = swam,
      minimum_short_window_activity_median_threshold = p$minimum_short_window_activity_median_threshold,
      short_window_activity_median_minimum_high_epochs = swamm_high_epochs,
      half_window_around_border = half_window_around_border,
      activity_median_analysis_window = activity_median_analysis_window,
      maximum_allowed_gap = refinement_gap,
      quantile_threshold = swam_threshold,
      median_filter_half_window_size = median_filter_half_window_size,
      getuptime_metric_method = p$getuptime_metric_method,
      getuptime_metric_parameter = p$getuptime_metric_parameter,
      condition = condition,
      getuptime_first_candidate_score = g[1],
      getuptime_best_median_difference_candidate_score = g[2],
      getuptime_best_crossing_distance_candidate_score = g[3],
      getuptime_best_epochs_above_metric_after_score = g[4],
      getuptime_thresholded_candidate_score_amplitude = g[5],
      getuptime_thresholded_candidate_score_minimum = g[6],
      getuptime_minimum_peak_length = getuptime_minimum_peak_length,
      getuptime_minimum_valley_length = getuptime_minimum_valley_length,
      after_candidate_window = after_candidate_window,
      getuptime_maximum_epochs_above_metric_after_candidate = getuptime_max_epochs,
      zero_proportion_threshold = zero_proportion_threshold,
      do_remove_after_long_tall_peak = p$getuptime_do_remove_after_long_tall_peak,
      getuptime_high_probability_awake_peak_length = p$getuptime_high_probability_awake_peak_length,
      getuptime_high_probability_sleep_valley_length = p$getuptime_high_probability_sleep_valley_length,
      do_remove_before_long_valley = p$getuptime_do_remove_before_long_valley,
      update_peaks_and_valleys = p$getuptime_update_peaks_and_valleys,
      getuptime_score_first_candidate = p$getuptime_score_first_candidate
    )

    for (k in seq_len(num)) {
      itc    <- trans_i[k]
      prev_t <- if (k >= 2L) refined_i[k - 1L] else 0L
      next_t <- if (k < num) refined_i[k + 1L] else (data_length + 1L)

      if (trans_dir[k] < 0) {                         # bed-time refinement
        res <- .bt_refine(bt, refinement_window_levels, itc, prev_t, next_t)
        refinement_window_levels <- res$refinement_window_levels
        rbt <- res$refined_bedtime
        refined_bedtimes <- c(refined_bedtimes, rbt)
        refined_i[k] <- rbt
        if (prev_t < rbt) refined_output[(prev_t + 1L):rbt] <- 1
      } else {                                        # get-up-time refinement
        res <- .gt_refine(gt, refinement_window_levels, itc, prev_t, next_t)
        refinement_window_levels <- res$refinement_window_levels
        rgt <- res$refined_getuptime
        refined_getuptimes <- c(refined_getuptimes, rgt)
        refined_i[k] <- rgt
      }
    }

    if (trans_dir[num] > 0) {
      last <- refined_i[num]
      if (last < data_length) refined_output[(last + 1L):data_length] <- 1
    }
  } else {
    refined_output <- rep(1, data_length)
  }

  # ── Stage 4: minimum-length filter ──
  # CSPD.model converts nap_minimum_length (minutes) to epochs here, then
  # branches on detect_naps for the post-processing.
  nap_minimum_length <- as.integer(round(p$nap_minimum_length * 60 / duration))
  if (isTRUE(p$detect_naps)) {
    # Python: boolean_length_filter(nap_min, "v") -> peak_valley_zero_proportion
    # _filter(1/3) -> boolean_length_filter(0.5*nap_min, "p").
    fsd <- .boolean_length_filter(nap_minimum_length, refined_output)
    fsd <- .peak_valley_zero_proportion_filter(1 / 3, fsd)
    fsd <- .boolean_length_filter(as.integer(0.5 * nap_minimum_length), fsd,
                                  class_to_filter = "p")
    final_sleep_detection <- fsd
  } else {
    final_sleep_detection <- .boolean_length_filter(p$sleep_minimum_length, refined_output)
  }
  refined_output <- final_sleep_detection

  # ── Per-night transitions (pre-heuristic) -> refined_sleep_df ──
  rb_borders <- diff(refined_output)
  ridx0      <- which(rb_borders != 0) - 1L
  npairs     <- length(ridx0) %/% 2L
  if (npairs > 0L) {
    dts <- as.POSIXct(datetime_stamps)
    bi  <- ridx0[seq(1L, by = 2L, length.out = npairs)]   # even (0-indexed) -> bedtime
    gi  <- ridx0[seq(2L, by = 2L, length.out = npairs)]   # odd  -> getuptime
    refined_sleep_df <- data.frame(
      bedtime         = dts[bi + 1L],
      getuptime       = dts[gi + 1L],
      bedtime_index   = as.integer(bi),
      getuptime_index = as.integer(gi)
    )
    refined_sleep_df$hour_length <-
      as.numeric(difftime(refined_sleep_df$getuptime, refined_sleep_df$bedtime, units = "hours"))
  } else {
    refined_sleep_df <- data.frame()
  }

  # ── Crespo heuristic: force the recording borders awake ──
  final_sleep_detection[1] <- 1
  final_sleep_detection[data_length] <- 1
  refined_output[1] <- 1
  refined_output[data_length] <- 1

  list(
    refined_output     = refined_output,
    refined_sleep_df   = refined_sleep_df,
    refined_bedtimes   = refined_bedtimes,
    refined_getuptimes = refined_getuptimes,
    transitions        = data.frame(index = as.integer(trans_i),
                                    direction = as.numeric(trans_dir))
  )
}
