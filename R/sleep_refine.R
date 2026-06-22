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
