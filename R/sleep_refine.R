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
#   2. sleep-gap separation             -> (pending)
#   3. bed-time / get-up-time refiners  -> (pending)
#   4. minimum-length filter            -> (pending; boolean_length_filter)

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
