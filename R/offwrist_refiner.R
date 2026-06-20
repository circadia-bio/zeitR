# ── Bimodal off-wrist refiner ─────────────────────────────────────────────────
# R port of BimodalOffwristRefiner from
# condor_pipeline/algorithms/vendor/condor/bimodal_offwrist_refine_without_prints.py
# Author of original Python: Julius A. P. P. de Paula (Condor Instruments, 2023)
#
# Only the ActTrust configuration is implemented here — i.e. the subset of
# filters that are active when called from bimodal_offwrist_wrapper_acttrust.
# Capacitive-sensor branches (ActLumus only) are omitted.
#
# Parameter values fixed to the ActTrust wrapper defaults:
#   do_onwrist_length_filter              = TRUE
#   do_onwrist_zero_activity_proportion_filter = FALSE
#   search_border_configuration           = "mod"
#   do_low_activity_filter                = TRUE
#   do_offwrist_length_filter             = TRUE
#   do_sleep_filter                       = TRUE
#   do_description_report_based_filter    = TRUE
#   do_valley_peak_algorithm              = TRUE
#   do_temperature_difference_filter      = TRUE
#   offwrists_need_decreasing_temperature = FALSE
#   do_forbidden_zone                     = TRUE
#   do_valley_peak_filter                 = FALSE
#   skip_description_filters              = FALSE
# ─────────────────────────────────────────────────────────────────────────────

# ── Public entry point ────────────────────────────────────────────────────────

#' Refine an initial off-wrist detection using the Condor three-stage algorithm
#'
#' Implements the `BimodalOffwristRefiner` from the circadiaBase pipeline,
#' configured for Condor ActTrust devices. The refiner applies three sequential
#' stages to improve an initial off-wrist estimate:
#'
#' 1. **First stage** — filters out short on-wrist periods, merging the
#'    surrounding off-wrist intervals.
#' 2. **Second stage** — refines the start and end of each off-wrist period
#'    by searching for temperature-variance peaks near the initial border
#'    candidates.
#' 3. **Third stage** — applies a battery of feature-based filters (low
#'    activity, off-wrist length, temperature difference median, sleep overlap,
#'    valley-peak algorithm, forbidden zones, description-report filters).
#'
#' @param initial_offwrist `integer` vector. Initial off-wrist estimate:
#'   `0` = off-wrist, `1` = on-wrist.
#' @param activity `numeric` vector. PIM activity counts.
#' @param activity_median `numeric` vector. Rolling median of activity.
#' @param temperature `numeric` vector. Internal (on-body) temperature.
#' @param norm_temp_variance `numeric` vector. Normalised rolling variance
#'   of temperature.
#' @param temp_derivative `numeric` vector. Five-point stencil temperature
#'   derivative.
#' @param temp_derivative_variance `numeric` vector. Rolling variance of the
#'   temperature derivative.
#' @param temperature_threshold `numeric(1)`. Temperature threshold separating
#'   on-wrist (above) from off-wrist (below) epochs.
#' @param ashman `numeric(1)`. Ashman's D statistic for the temperature
#'   distribution bimodality.
#' @param activity_median_low `integer` vector. `1` = low activity median
#'   epoch, `0` otherwise.
#' @param is_low_temp `logical` vector. `TRUE` = epoch below temperature
#'   threshold.
#' @param filter_hws `integer(1)`. Half-window size used in feature
#'   extraction filters. Default is `10`.
#' @param dif_temp `numeric` vector. Internal minus external temperature
#'   difference (required for the temperature difference median filter).
#' @param epoch_hour `numeric(1)`. Number of epochs per hour.
#' @param do_near_all_off_detection `logical(1)`. If `TRUE`, apply near-all-
#'   off detection when the signal is unimodal. Default is `TRUE`.
#'
#' @return `integer` vector, same length as `initial_offwrist`:
#'   `0` = off-wrist (refined), `1` = on-wrist.
#'
#' @references
#' Algorithm by Julius A. P. P. de Paula, Condor Instruments (2023).
#' Unpublished; source available in the circadiaBase pipeline repository.
#'
#' @noRd
.bimodal_refine_acttrust <- function(
    initial_offwrist,
    activity,
    activity_median,
    temperature,
    norm_temp_variance,
    temp_derivative,
    temp_derivative_variance,
    temperature_threshold,
    ashman,
    activity_median_low,
    is_low_temp,
    filter_hws         = 10L,
    dif_temp,
    epoch_hour         = 60L,
    do_near_all_off_detection = TRUE
) {
  n <- length(initial_offwrist)

  # ── ActTrust wrapper parameter values ─────────────────────────────────────
  minimum_onwrist_length               <- 20L
  minimum_offwrist_length              <- 10L
  minimum_preceding_onwrist_length     <- 40L
  activity_threshold_quantile          <- 0.24   # from wrapper
  minimum_low_activity_proportion      <- 1.0    # bimodal_minimum_low_activity_proportion
  max_low_temp_var_prop_border         <- 0.3
  long_offwrist_length                 <- 4L * as.integer(epoch_hour)
  short_offwrist_length                <- 20L
  temp_var_thr_quantile                <- 0.9
  sleep_act_hws                        <- 120L
  sleep_all_q                          <- 0.4
  sleep_pos_q                          <- 0.05
  max_offwrist_sleep_prop              <- 0.4
  ashman_d_min                         <- 1.5
  ashman_d_max                         <- 2.6
  bimodal_max_offwrist_prop            <- 0.0    # from wrapper
  bimodal_min_low_act_prop             <- 1.0    # from wrapper
  report_zero_act_prop_min             <- 0.35
  border_conc_min                      <- 0.5
  report_act_around_min                <- 0.1
  report_low_act_prop_min              <- 0.5
  offwrist_max_temp_dif_median         <- 0.8    # from wrapper
  offwrist_min_temp_dif_median         <- 0.65   # from wrapper
  valley_quantile                      <- 0.99   # from wrapper
  peak_quantile                        <- 0.99   # from wrapper
  valley_peak_low_temp_min             <- 0.5
  short_vp_decrease_ratio_min         <- 1.25
  decrease_ratio_min                   <- 1.0    # from wrapper
  sleep_low_temp_prop_max              <- 0.5
  half_filter_hws                      <- as.integer(filter_hws / 2)

  # ── Trivial case ───────────────────────────────────────────────────────────
  offwrist_periods <- .rle_periods(initial_offwrist == 0L)
  if (nrow(offwrist_periods) == 0L) {
    return(rep(1L, n))
  }

  # ── Activity threshold ─────────────────────────────────────────────────────
  act_zero_prop  <- zero_prop(activity)
  act_thr_q      <- act_zero_prop + (1 - act_zero_prop) * activity_threshold_quantile
  activity_thr   <- as.double(stats::quantile(activity, act_thr_q,
                                               names = FALSE, type = 1))
  if (activity_thr < 200) activity_thr <- 200

  pos_act_q      <- act_zero_prop + (1 - act_zero_prop) * 0.5
  positive_act_median <- as.double(stats::quantile(activity, pos_act_q,
                                                    names = FALSE, type = 1))

  temp_var_thr   <- as.double(stats::quantile(norm_temp_variance,
                                               temp_var_thr_quantile,
                                               names = FALSE))

  # ── Stage 1: filter short on-wrist periods ─────────────────────────────────
  onwrist_periods <- .rle_periods(initial_offwrist == 1L)
  onwrist_periods <- onwrist_periods[
    (onwrist_periods$end - onwrist_periods$start) > minimum_onwrist_length, ]
  onwrist_periods <- onwrist_periods[
    order(onwrist_periods$start), ]
  row.names(onwrist_periods) <- NULL

  # Rebuild off-wrist from surviving on-wrist periods
  stage1_onwrist <- integer(n)
  for (i in seq_len(nrow(onwrist_periods))) {
    stage1_onwrist[onwrist_periods$start[i]:onwrist_periods$end[i]] <- 1L
  }
  offwrist_periods <- .rle_periods(stage1_onwrist == 0L)

  if (nrow(offwrist_periods) == 0L) {
    return(rep(1L, n))
  }

  # ── Stage 2: border refinement ─────────────────────────────────────────────
  refined_offwrist_periods <- .refine_borders(
    offwrist_periods      = offwrist_periods,
    norm_temp_variance    = norm_temp_variance,
    activity_median       = activity_median,
    activity_median_low   = activity_median_low,
    temperature           = temperature,
    temperature_threshold = temperature_threshold,
    activity_thr          = activity_thr,
    minimum_low_act_prop  = minimum_low_activity_proportion,
    min_preceding_onwrist = minimum_preceding_onwrist_length,
    filter_hws            = filter_hws,
    half_filter_hws       = half_filter_hws,
    max_low_tv_border     = max_low_temp_var_prop_border,
    temp_var_thr          = temp_var_thr,
    activity              = activity,
    n                     = n
  )

  if (nrow(refined_offwrist_periods) == 0L) {
    return(rep(1L, n))
  }

  # ── Stage 3 initial filters ────────────────────────────────────────────────

  # Compute features per period
  refined_offwrist_periods$length <-
    refined_offwrist_periods$end - refined_offwrist_periods$start

  refined_offwrist_periods$low_temp_prop <- vapply(
    seq_len(nrow(refined_offwrist_periods)), function(i) {
      below_prop(temperature[refined_offwrist_periods$start[i]:
                               refined_offwrist_periods$end[i]],
                 temperature_threshold)
    }, numeric(1))

  refined_offwrist_periods$low_act_prop <- vapply(
    seq_len(nrow(refined_offwrist_periods)), function(i) {
      below_prop(activity[refined_offwrist_periods$start[i]:
                            refined_offwrist_periods$end[i]],
                 activity_thr)
    }, numeric(1))

  refined_offwrist_periods$zero_act_prop <- vapply(
    seq_len(nrow(refined_offwrist_periods)), function(i) {
      zero_prop(activity[refined_offwrist_periods$start[i]:
                           refined_offwrist_periods$end[i]])
    }, numeric(1))

  # Low-activity filter
  refined_offwrist_periods <- refined_offwrist_periods[
    refined_offwrist_periods$low_act_prop >= minimum_low_activity_proportion, ]
  row.names(refined_offwrist_periods) <- NULL

  # Off-wrist length filter
  refined_offwrist_periods <- refined_offwrist_periods[
    refined_offwrist_periods$length >= minimum_offwrist_length, ]
  row.names(refined_offwrist_periods) <- NULL

  if (nrow(refined_offwrist_periods) == 0L) {
    return(rep(1L, n))
  }

  # Temperature difference median filter (ActTrust specific)
  refined_offwrist_periods$temp_dif_median <- vapply(
    seq_len(nrow(refined_offwrist_periods)), function(i) {
      stats::median(dif_temp[refined_offwrist_periods$start[i]:
                               refined_offwrist_periods$end[i]],
                    na.rm = TRUE)
    }, numeric(1))

  refined_offwrist_periods <- refined_offwrist_periods[
    refined_offwrist_periods$temp_dif_median < offwrist_max_temp_dif_median, ]
  row.names(refined_offwrist_periods) <- NULL

  if (nrow(refined_offwrist_periods) == 0L) {
    return(rep(1L, n))
  }

  # ── Mark long off-wrist periods (>= 4 h) ──────────────────────────────────
  # Long periods are exempt from sleep filter
  is_long <- refined_offwrist_periods$length >= long_offwrist_length
  short_only_mask <- rep(TRUE, n)
  for (i in which(is_long)) {
    short_only_mask[refined_offwrist_periods$start[i]:
                      refined_offwrist_periods$end[i]] <- FALSE
  }
  valid_activity    <- activity[short_only_mask]
  valid_temperature <- temperature[short_only_mask]

  # ── Sleep estimation for sleep filter ─────────────────────────────────────
  estimated_sleep <- .estimate_sleep(
    activity          = activity,
    valid_activity    = valid_activity,
    n                 = n,
    short_only_mask   = short_only_mask,
    sleep_act_hws     = sleep_act_hws,
    sleep_all_q       = sleep_all_q,
    sleep_pos_q       = sleep_pos_q,
    temperature       = temperature,
    temperature_threshold = temperature_threshold,
    sleep_low_temp_prop_max = sleep_low_temp_prop_max,
    epoch_hour        = epoch_hour
  )

  # ── Sleep filter ───────────────────────────────────────────────────────────
  refined_offwrist_periods$sleep_prop <- vapply(
    seq_len(nrow(refined_offwrist_periods)), function(i) {
      seg <- estimated_sleep[refined_offwrist_periods$start[i]:
                               refined_offwrist_periods$end[i]]
      mean(seg == 1L)
    }, numeric(1))

  refined_offwrist_periods <- refined_offwrist_periods[
    refined_offwrist_periods$sleep_prop <= max_offwrist_sleep_prop |
      is_long[seq_len(nrow(refined_offwrist_periods))], ]
  row.names(refined_offwrist_periods) <- NULL

  if (nrow(refined_offwrist_periods) == 0L) {
    return(rep(1L, n))
  }

  # ── Forbidden zones ────────────────────────────────────────────────────────
  forbidden_zone <- .compute_forbidden_zone(estimated_sleep, epoch_hour, n)

  refined_offwrist_periods$in_forbidden <- vapply(
    seq_len(nrow(refined_offwrist_periods)), function(i) {
      any(forbidden_zone[refined_offwrist_periods$start[i]:
                           refined_offwrist_periods$end[i]] == 1L)
    }, logical(1))

  refined_offwrist_periods <- refined_offwrist_periods[
    !refined_offwrist_periods$in_forbidden, ]
  row.names(refined_offwrist_periods) <- NULL

  # ── Bimodality check ───────────────────────────────────────────────────────
  is_bimodal <- .check_bimodality(
    initial_offwrist        = initial_offwrist,
    activity_median         = activity_median,
    activity_thr            = activity_thr,
    ashman                  = ashman,
    ashman_d_min            = ashman_d_min,
    ashman_d_max            = ashman_d_max,
    bimodal_max_offwrist_p  = bimodal_max_offwrist_prop,
    bimodal_min_low_act_p   = bimodal_min_low_act_prop,
    n                       = n
  )

  # ── Valley-peak algorithm ─────────────────────────────────────────────────
  if (is_bimodal) {
    vp_periods <- .valley_peak_algorithm(
      temp_derivative       = temp_derivative,
      temperature           = temperature,
      temperature_threshold = temperature_threshold,
      estimated_sleep       = estimated_sleep,
      valley_quantile       = valley_quantile,
      peak_quantile         = peak_quantile,
      minimum_offwrist_length = minimum_offwrist_length,
      short_offwrist_length   = short_offwrist_length,
      short_vp_decrease_ratio_min = short_vp_decrease_ratio_min,
      valley_peak_low_temp_min   = valley_peak_low_temp_min,
      epoch_hour            = epoch_hour,
      n                     = n,
      next_possible_length  = TRUE,
      short_criteria        = TRUE,
      medium_criteria       = TRUE
    )

    # Merge: long periods from sleep-filtered + short from both + vp periods
    long_df  <- refined_offwrist_periods[
      refined_offwrist_periods$length >= long_offwrist_length, ]
    short_df <- refined_offwrist_periods[
      refined_offwrist_periods$length <  long_offwrist_length, ]

    if (nrow(vp_periods) > 0) {
      vp_periods$length       <- vp_periods$end - vp_periods$start
      vp_periods$valley_peak  <- TRUE
      short_df$valley_peak    <- FALSE
      long_df$valley_peak     <- FALSE

      combined <- rbind(
        long_df[, c("start", "end", "length", "valley_peak")],
        short_df[, c("start", "end", "length", "valley_peak")],
        vp_periods[, c("start", "end", "length", "valley_peak")]
      )
    } else {
      short_df$valley_peak   <- FALSE
      long_df$valley_peak    <- FALSE
      combined <- rbind(
        long_df[, c("start", "end", "length", "valley_peak")],
        short_df[, c("start", "end", "length", "valley_peak")]
      )
    }

    combined <- combined[order(combined$start), ]
    row.names(combined) <- NULL
    refined_offwrist_periods <- combined

    # Recompute features after merge
    refined_offwrist_periods$low_act_prop <- vapply(
      seq_len(nrow(refined_offwrist_periods)), function(i) {
        below_prop(activity[refined_offwrist_periods$start[i]:
                              refined_offwrist_periods$end[i]], activity_thr)
      }, numeric(1))
    refined_offwrist_periods$zero_act_prop <- vapply(
      seq_len(nrow(refined_offwrist_periods)), function(i) {
        zero_prop(activity[refined_offwrist_periods$start[i]:
                             refined_offwrist_periods$end[i]])
      }, numeric(1))
    refined_offwrist_periods$low_temp_prop <- vapply(
      seq_len(nrow(refined_offwrist_periods)), function(i) {
        below_prop(temperature[refined_offwrist_periods$start[i]:
                                 refined_offwrist_periods$end[i]],
                   temperature_threshold)
      }, numeric(1))
    refined_offwrist_periods$length <-
      refined_offwrist_periods$end - refined_offwrist_periods$start

  } else {
    # Unimodal path
    if (do_near_all_off_detection) {
      refined_offwrist_periods <- .near_all_off_detection(
        activity        = activity,
        activity_thr    = activity_thr,
        positive_act_median = positive_act_median,
        n               = n,
        epoch_hour      = epoch_hour
      )
    } else {
      refined_offwrist_periods <- data.frame(
        start = integer(), end = integer(), length = integer(),
        valley_peak = logical()
      )
    }

    if (nrow(refined_offwrist_periods) == 0L) {
      return(rep(1L, n))
    }
  }

  # ── Description-report-based filter ───────────────────────────────────────
  if (nrow(refined_offwrist_periods) > 0L) {
    report <- .describe_offwrist_periods(
      offwrist_periods      = refined_offwrist_periods,
      activity              = activity,
      temperature           = temperature,
      temperature_variance  = norm_temp_variance,
      temperature_threshold = temperature_threshold,
      activity_thr          = activity_thr,
      dif_temp              = dif_temp,
      n                     = n
    )

    refined_offwrist_periods <- .description_report_filter(
      periods               = refined_offwrist_periods,
      report                = report,
      report_zero_act_min   = report_zero_act_prop_min,
      border_conc_min       = border_conc_min,
      act_around_min        = report_act_around_min,
      low_act_prop_min      = report_low_act_prop_min,
      temp_dif_min          = offwrist_min_temp_dif_median,
      temp_dif_max          = offwrist_max_temp_dif_median
    )
    row.names(refined_offwrist_periods) <- NULL
  }

  # ── Surrounded on-wrist filter ─────────────────────────────────────────────
  refined_offwrist_periods <- .surrounded_onwrist_filter(
    periods = refined_offwrist_periods,
    n       = n
  )

  # ── Final border snap ─────────────────────────────────────────────────────
  if (nrow(refined_offwrist_periods) > 0L) {
    refined_offwrist_periods$length <-
      refined_offwrist_periods$end - refined_offwrist_periods$start

    if (refined_offwrist_periods$start[1] <= minimum_offwrist_length) {
      refined_offwrist_periods$start[1] <- 0L
    }
    last <- nrow(refined_offwrist_periods)
    if ((n - refined_offwrist_periods$end[last]) <= minimum_offwrist_length) {
      refined_offwrist_periods$end[last] <- n
    }
  }

  # ── Assemble final output ─────────────────────────────────────────────────
  refined <- rep(1L, n)
  for (i in seq_len(nrow(refined_offwrist_periods))) {
    s <- refined_offwrist_periods$start[i]
    e <- min(refined_offwrist_periods$end[i], n)
    if (s < e) refined[s:e] <- 0L
  }

  refined
}

# ── Internal stage helpers ────────────────────────────────────────────────────

#' Convert a binary vector to a data frame of contiguous runs of 1s
#' @noRd
.rle_periods <- function(x) {
  r      <- rle(as.integer(x))
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  keep   <- r$values == 1L
  if (!any(keep)) return(data.frame(start = integer(), end = integer()))
  data.frame(start = starts[keep], end = ends[keep])
}

#' Proportion of values in x below threshold
#' @noRd
below_prop <- function(x, thr) {
  if (length(x) == 0L) return(0)
  sum(x < thr, na.rm = TRUE) / length(x)
}

#' Refine borders of each off-wrist period using temperature variance peaks
#' @noRd
.refine_borders <- function(
    offwrist_periods, norm_temp_variance, activity_median, activity_median_low,
    temperature, temperature_threshold, activity_thr, minimum_low_act_prop,
    min_preceding_onwrist, filter_hws, half_filter_hws, max_low_tv_border,
    temp_var_thr, activity, n
) {
  refined <- list()
  periods <- as.matrix(offwrist_periods[, c("start", "end")])
  n_off   <- nrow(periods)
  idx     <- 1L
  prev_end <- 0L

  while (idx <= n_off) {
    start <- periods[idx, 1L]
    end   <- periods[idx, 2L]

    # ── Refine start ────────────────────────────────────────────────────────
    if (start < prev_end) start <- prev_end + 2L
    start_ok <- FALSE

    if (start >= filter_hws) {
      if (norm_temp_variance[start] >= temp_var_thr) {
        lo <- max(1L, start - filter_hws)
        ltvp <- below_prop(norm_temp_variance[lo:start], temp_var_thr)
        if (ltvp <= max_low_tv_border) {
          # peak found at start
          new_start <- lo + which.max(norm_temp_variance[lo:start]) - 1L
          new_start <- max(new_start, prev_end + 1L)
          refined[[length(refined) + 1L]] <- c(new_start, 0L)
          start_ok <- TRUE
        } else {
          # look backwards for peak base
          res <- .find_peak_base_start(start, prev_end, norm_temp_variance,
                                       activity_median, activity_median_low,
                                       temperature, temperature_threshold,
                                       activity_thr, temp_var_thr, filter_hws)
          if (!is.null(res)) {
            refined[[length(refined) + 1L]] <- c(res, 0L)
            start_ok <- TRUE
          } else {
            idx <- idx + 1L
            next
          }
        }
      } else {
        res <- .find_peak_base_start(start, prev_end, norm_temp_variance,
                                     activity_median, activity_median_low,
                                     temperature, temperature_threshold,
                                     activity_thr, temp_var_thr, filter_hws)
        if (!is.null(res)) {
          refined[[length(refined) + 1L]] <- c(res, 0L)
          start_ok <- TRUE
        } else {
          idx <- idx + 1L
          next
        }
      }
    } else {
      new_start <- which.max(norm_temp_variance[1:max(1L, start)]) - 1L
      refined[[length(refined) + 1L]] <- c(new_start, 0L)
      start_ok <- TRUE
    }

    if (!start_ok) { idx <- idx + 1L; next }

    # ── Refine end ──────────────────────────────────────────────────────────
    end_ok   <- FALSE
    end_done <- FALSE

    if (norm_temp_variance[min(end, n)] >= temp_var_thr) {
      hi <- min(n, end + filter_hws)
      ltvp <- below_prop(norm_temp_variance[end:hi], temp_var_thr)
      if (ltvp <= max_low_tv_border) {
        # peak found at end
        new_end <- end - 1L + which.max(norm_temp_variance[end:hi])
        refined[[length(refined)]][2L] <- new_end
        prev_end <- new_end
        idx <- idx + 1L
        end_done <- TRUE
      }
    }

    if (!end_done) {
      # Try to find peak base forward
      if (idx < n_off) {
        next_start <- periods[idx + 1L, 1L]
        gap <- next_start - end
        # Check if too close
        too_close <- FALSE
        if (gap <= min_preceding_onwrist) {
          act_seg <- activity[end:next_start]
          if (below_prop(act_seg, activity_thr) > minimum_low_act_prop) {
            too_close <- TRUE
          }
        }
        if (too_close) {
          # merge: remove current end entry, continue to next period as end
          periods <- periods[-idx, , drop = FALSE]
          n_off   <- n_off - 1L
        } else {
          # search forward for peak
          res <- .find_peak_base_end(end, next_start, norm_temp_variance,
                                     activity_median, activity_median_low,
                                     temperature, temperature_threshold,
                                     activity_thr, temp_var_thr, filter_hws)
          if (!is.null(res)) {
            refined[[length(refined)]][2L] <- res
            prev_end <- res
          } else {
            # Remove this period — no valid end found
            refined[[length(refined)]] <- NULL
          }
          idx <- idx + 1L
        }
      } else {
        # Last period — use max variance in remaining signal
        hi <- min(n, end + filter_hws)
        new_end <- end - 1L + which.max(norm_temp_variance[end:hi])
        refined[[length(refined)]][2L] <- new_end
        prev_end <- new_end
        idx <- idx + 1L
      }
    }
  }

  # Convert to data frame, drop any with end == 0
  if (length(refined) == 0L) {
    return(data.frame(start = integer(), end = integer()))
  }
  m <- do.call(rbind, refined)
  m <- m[m[, 2L] > 0L, , drop = FALSE]
  if (nrow(m) == 0L) return(data.frame(start = integer(), end = integer()))
  data.frame(start = m[, 1L], end = m[, 2L])
}

#' Search backwards for a temperature-variance peak base
#' @noRd
.find_peak_base_start <- function(start, prev_end, ntv, act_median,
                                   act_med_low, temperature,
                                   temperature_threshold, activity_thr,
                                   temp_var_thr, filter_hws) {
  i <- start - 1L
  last_valid <- start
  peak_base  <- start

  while (i > prev_end) {
    valid <- .check_valid_border_mod(act_med_low, temperature,
                                     temperature_threshold, act_median,
                                     activity_thr, i)
    if (valid) {
      if (ntv[i] >= temp_var_thr) {
        peak_base <- i
        break
      }
      i <- i - 1L
    } else {
      last_valid <- i
      break
    }
  }

  if (last_valid == start) {
    if (peak_base < start) {
      # Walk to actual peak
      peak <- peak_base - 1L
      while (peak > 1L && ntv[peak] >= ntv[peak + 1L]) peak <- peak - 1L
      return(peak)
    } else {
      return(NULL)  # filter out
    }
  } else {
    return(last_valid)
  }
}

#' Search forwards for a temperature-variance peak base
#' @noRd
.find_peak_base_end <- function(end, next_start, ntv, act_median,
                                 act_med_low, temperature,
                                 temperature_threshold, activity_thr,
                                 temp_var_thr, filter_hws) {
  i <- end + 1L
  last_valid <- end
  peak_base  <- end

  while (i < next_start) {
    valid <- .check_valid_border_mod(act_med_low, temperature,
                                     temperature_threshold, act_median,
                                     activity_thr, i)
    if (valid) {
      if (ntv[i] >= temp_var_thr) {
        peak_base <- i
        break
      }
      i <- i + 1L
    } else {
      last_valid <- i
      break
    }
  }

  if (last_valid == end) {
    if (peak_base > end) {
      peak <- peak_base + 1L
      while (peak < length(ntv) && ntv[peak] >= ntv[peak - 1L]) peak <- peak + 1L
      return(peak)
    } else {
      return(NULL)
    }
  } else {
    return(last_valid)
  }
}

#' "mod" border validity check (activity median low OR low temp + low act)
#' @noRd
.check_valid_border_mod <- function(act_med_low, temperature,
                                     temperature_threshold, act_median,
                                     activity_thr, i) {
  act_med_low[i] == 1L ||
    (temperature[i] < temperature_threshold &&
       act_median[i] < 2 * activity_thr)
}

#' Estimate sleep periods for the sleep filter
#' @noRd
.estimate_sleep <- function(activity, valid_activity, n, short_only_mask,
                             sleep_act_hws, sleep_all_q, sleep_pos_q,
                             temperature, temperature_threshold,
                             sleep_low_temp_prop_max, epoch_hour) {
  # Compute threshold from valid (non-long-offwrist) activity
  act_zp     <- zero_prop(valid_activity)
  all_q_idx  <- act_zp + (1 - act_zp) * sleep_all_q
  pos_q_idx  <- act_zp + (1 - act_zp) * sleep_pos_q
  thr_all    <- as.double(stats::quantile(valid_activity, all_q_idx,
                                          names = FALSE))
  thr_pos    <- as.double(stats::quantile(valid_activity, pos_q_idx,
                                          names = FALSE))
  sleep_thr  <- max(thr_all, thr_pos)

  # Rolling median filter
  act_med_sleep <- median_filter(activity, sleep_act_hws)

  # Threshold
  is_sleep <- as.integer(act_med_sleep < sleep_thr)

  # Filter short sleep periods (< 1 hour)
  sleep_periods <- .rle_periods(is_sleep == 1L)
  min_sleep_len <- as.integer(epoch_hour)
  if (nrow(sleep_periods) > 0L) {
    sleep_periods <- sleep_periods[
      (sleep_periods$end - sleep_periods$start) >= min_sleep_len, ]
    is_sleep <- integer(n)
    for (i in seq_len(nrow(sleep_periods))) {
      is_sleep[sleep_periods$start[i]:sleep_periods$end[i]] <- 1L
    }
  }

  # Remove sleep periods with too high low-temperature proportion
  sleep_periods <- .rle_periods(is_sleep == 1L)
  if (nrow(sleep_periods) > 0L) {
    valid_sleep <- vapply(seq_len(nrow(sleep_periods)), function(i) {
      below_prop(temperature[sleep_periods$start[i]:sleep_periods$end[i]],
                 temperature_threshold) <= sleep_low_temp_prop_max
    }, logical(1))
    sleep_periods <- sleep_periods[valid_sleep, ]
    is_sleep <- integer(n)
    for (i in seq_len(nrow(sleep_periods))) {
      is_sleep[sleep_periods$start[i]:sleep_periods$end[i]] <- 1L
    }
  }

  is_sleep
}

#' Compute forbidden zones at the centre of sleep periods
#' @noRd
.compute_forbidden_zone <- function(estimated_sleep, epoch_hour, n) {
  forbidden <- integer(n)
  sleep_periods <- .rle_periods(estimated_sleep == 1L)
  if (nrow(sleep_periods) == 0L) return(forbidden)

  # Forbidden zone = central 50% of each sleep period
  for (i in seq_len(nrow(sleep_periods))) {
    s   <- sleep_periods$start[i]
    e   <- sleep_periods$end[i]
    len <- e - s
    if (len > 0L) {
      q1 <- s + as.integer(len * 0.25)
      q3 <- s + as.integer(len * 0.75)
      forbidden[q1:q3] <- 1L
    }
  }
  forbidden
}

#' Check if the temperature distribution is bimodal
#' @noRd
.check_bimodality <- function(initial_offwrist, activity_median, activity_thr,
                               ashman, ashman_d_min, ashman_d_max,
                               bimodal_max_offwrist_p, bimodal_min_low_act_p, n) {
  offwrist_prop  <- mean(initial_offwrist == 0L)
  low_act_prop   <- mean(activity_median < activity_thr)

  if (offwrist_prop > bimodal_max_offwrist_p) return(FALSE)
  if (low_act_prop  > bimodal_min_low_act_p)  return(FALSE)

  ashman >= ashman_d_min
}

#' Valley-peak off-wrist detection from temperature derivative
#' @noRd
.valley_peak_algorithm <- function(
    temp_derivative, temperature, temperature_threshold, estimated_sleep,
    valley_quantile, peak_quantile, minimum_offwrist_length,
    short_offwrist_length, short_vp_decrease_ratio_min,
    valley_peak_low_temp_min, epoch_hour, n,
    next_possible_length, short_criteria, medium_criteria
) {
  # Identify daily valleys (sharp drops) and peaks (sharp rises)
  valley_thr <- as.double(stats::quantile(temp_derivative, 1 - valley_quantile,
                                          names = FALSE))
  peak_thr   <- as.double(stats::quantile(temp_derivative, peak_quantile,
                                          names = FALSE))

  valley_idx <- which(temp_derivative <= valley_thr)
  peak_idx   <- which(temp_derivative >= peak_thr)

  if (length(valley_idx) == 0L || length(peak_idx) == 0L) {
    return(data.frame(start = integer(), end = integer(), valley_peak = logical()))
  }

  candidates <- list()

  for (v in valley_idx) {
    # Find next peak after this valley
    following_peaks <- peak_idx[peak_idx > v]
    if (length(following_peaks) == 0L) next

    p <- following_peaks[1L]
    length_vp <- p - v

    if (length_vp < minimum_offwrist_length) {
      # Try next valley if next_possible_length is TRUE
      if (next_possible_length) {
        next_valleys <- valley_idx[valley_idx > v & valley_idx < p]
        if (length(next_valleys) > 0L) {
          v2 <- next_valleys[length(next_valleys)]
          if ((p - v2) >= minimum_offwrist_length) {
            v <- v2
            length_vp <- p - v
          }
        }
      }
      if (length_vp < minimum_offwrist_length) next
    }

    seg_temp <- temperature[v:p]
    low_temp_p <- below_prop(seg_temp, temperature_threshold)

    if (low_temp_p < valley_peak_low_temp_min) next

    # Check overlap with estimated sleep
    sleep_overlap <- mean(estimated_sleep[v:p] == 1L)
    if (sleep_overlap > 0.4) next

    if (length_vp < short_offwrist_length && short_criteria) {
      # Short off-wrist: require temperature decrease ratio
      decrease_ratio <- .compute_decrease_ratio(temp_derivative, v, p)
      if (decrease_ratio < short_vp_decrease_ratio_min) next
    }

    candidates[[length(candidates) + 1L]] <- c(v, p)
  }

  if (length(candidates) == 0L) {
    return(data.frame(start = integer(), end = integer(), valley_peak = logical()))
  }

  m <- do.call(rbind, candidates)
  data.frame(start = m[, 1L], end = m[, 2L], valley_peak = TRUE)
}

#' Compute temperature derivative decrease ratio for a valley-peak segment
#' @noRd
.compute_decrease_ratio <- function(temp_derivative, v, p) {
  seg <- temp_derivative[v:p]
  neg <- seg[seg < 0]
  pos <- seg[seg > 0]
  if (length(pos) == 0L || sum(pos) == 0) return(0)
  abs(sum(neg)) / sum(pos)
}

#' Near-all-off detection for unimodal recordings
#' @noRd
.near_all_off_detection <- function(activity, activity_thr,
                                     positive_act_median, n, epoch_hour) {
  # If activity is very low overall, the whole recording may be off-wrist
  act_zp <- zero_prop(activity)
  if (act_zp > 0.9 || positive_act_median < activity_thr) {
    return(data.frame(start = 0L, end = n, length = n, valley_peak = FALSE))
  }
  data.frame(start = integer(), end = integer(),
             length = integer(), valley_peak = logical())
}

#' Compute description report for off-wrist periods
#' @noRd
.describe_offwrist_periods <- function(offwrist_periods, activity, temperature,
                                        temperature_variance, temperature_threshold,
                                        activity_thr, dif_temp, n,
                                        segments = 7L, window = 60L) {
  np <- nrow(offwrist_periods)
  if (np == 0L) return(data.frame())

  act_zp   <- zero_prop(activity)
  act_tq   <- act_zp + (1 - act_zp) * 0.05
  act_thr2 <- as.double(stats::quantile(activity, act_tq, names = FALSE, type = 1))

  dif_temp_var <- var_filter(dif_temp, 3L)

  report <- data.frame(
    activity_zero_prop         = numeric(np),
    low_act_prop               = numeric(np),
    high_act_before            = numeric(np),
    high_act_after             = numeric(np),
    start_act_weight           = numeric(np),
    end_act_weight             = numeric(np),
    border_act_conc            = numeric(np),
    start_temp_weight          = numeric(np),
    end_temp_weight            = numeric(np),
    border_temp_conc           = numeric(np),
    low_temp_prop              = numeric(np),
    temp_dif_median            = numeric(np),
    temp_dif_variance          = numeric(np)
  )

  for (i in seq_len(np)) {
    s <- offwrist_periods$start[i]
    e <- offwrist_periods$end[i]
    e <- min(e, n)

    seg_act      <- activity[s:e]
    seg_temp     <- temperature[s:e]
    seg_tv       <- temperature_variance[s:e]
    seg_dif      <- dif_temp[s:e]
    seg_dif_var  <- dif_temp_var[s:e]

    act_segs <- .segmentation(seg_act, segments)
    tv_segs  <- .segmentation(seg_tv,  segments)

    before_s <- max(1L, s - window)
    after_e  <- min(n, e + window)

    report$activity_zero_prop[i]  <- zero_prop(seg_act)
    report$low_act_prop[i]        <- below_prop(seg_act, activity_thr)
    report$high_act_before[i]     <- 1 - below_prop(activity[before_s:(s - 1L)],
                                                      act_thr2)
    report$high_act_after[i]      <- 1 - below_prop(activity[(e + 1L):after_e],
                                                      act_thr2)
    report$start_act_weight[i]    <- act_segs[1L]
    report$end_act_weight[i]      <- act_segs[segments]
    report$start_temp_weight[i]   <- tv_segs[1L]
    report$end_temp_weight[i]     <- tv_segs[segments]
    report$low_temp_prop[i]       <- below_prop(seg_temp, temperature_threshold)
    report$temp_dif_median[i]     <- stats::median(seg_dif, na.rm = TRUE)
    report$temp_dif_variance[i]   <- stats::median(seg_dif_var, na.rm = TRUE)
  }

  report$border_act_conc  <- report$start_act_weight + report$end_act_weight
  report$border_temp_conc <- report$start_temp_weight + report$end_temp_weight

  report
}

#' Segment a vector and return the proportion of total in each segment
#' @noRd
.segmentation <- function(x, n_segs = 7L) {
  total <- sum(x, na.rm = TRUE)
  if (total == 0) return(rep(0, n_segs))
  len   <- length(x)
  bins  <- round(seq(0, len, length.out = n_segs + 1L))
  vapply(seq_len(n_segs), function(i) {
    sum(x[(bins[i] + 1L):bins[i + 1L]], na.rm = TRUE) / total
  }, numeric(1))
}

#' Description-report-based filter
#' @noRd
.description_report_filter <- function(periods, report,
                                        report_zero_act_min,
                                        border_conc_min,
                                        act_around_min,
                                        low_act_prop_min,
                                        temp_dif_min,
                                        temp_dif_max) {
  if (nrow(periods) == 0L || nrow(report) == 0L) return(periods)

  keep <- rep(TRUE, nrow(periods))

  for (i in seq_len(nrow(periods))) {
    r <- report[i, ]

    # Zero activity proportion filter
    if (r$activity_zero_prop < report_zero_act_min) {
      keep[i] <- FALSE; next
    }

    # Border concentration filter
    if (r$border_act_conc < border_conc_min ||
        r$border_temp_conc < border_conc_min) {
      keep[i] <- FALSE; next
    }

    # Activity around borders filter
    act_around <- (r$high_act_before + r$high_act_after) / 2
    if (act_around < act_around_min) {
      keep[i] <- FALSE; next
    }

    # Low activity proportion filter
    if (r$low_act_prop < low_act_prop_min) {
      keep[i] <- FALSE; next
    }

    # Temperature difference median filter (ActTrust)
    if (r$temp_dif_median > temp_dif_max ||
        r$temp_dif_median < temp_dif_min) {
      keep[i] <- FALSE; next
    }
  }

  periods[keep, , drop = FALSE]
}

#' Remove on-wrist periods surrounded by off-wrist on both sides
#' @noRd
.surrounded_onwrist_filter <- function(periods, n) {
  if (nrow(periods) < 2L) return(periods)

  # Rebuild off-wrist binary
  offwrist_bin <- integer(n)
  for (i in seq_len(nrow(periods))) {
    s <- periods$start[i]; e <- min(periods$end[i], n)
    if (s < e) offwrist_bin[s:e] <- 1L
  }

  onwrist_periods <- .rle_periods(offwrist_bin == 0L)
  if (nrow(onwrist_periods) == 0L) return(periods)

  # An on-wrist period is "surrounded" if both neighbouring off-wrist periods
  # exist and it is shorter than the minimum meaningful on-wrist duration (20)
  remove_onwrist <- vapply(seq_len(nrow(onwrist_periods)), function(i) {
    len <- onwrist_periods$end[i] - onwrist_periods$start[i]
    s   <- onwrist_periods$start[i]
    e   <- onwrist_periods$end[i]
    preceded  <- any(periods$end   <= s)
    succeeded <- any(periods$start >= e)
    len < 20L && preceded && succeeded
  }, logical(1))

  if (!any(remove_onwrist)) return(periods)

  # Expand off-wrist periods to swallow surrounded on-wrist
  for (j in which(remove_onwrist)) {
    s <- onwrist_periods$start[j]
    e <- onwrist_periods$end[j]
    # Find preceding and succeeding off-wrist periods
    prec <- which(periods$end == s)
    succ <- which(periods$start == e)
    if (length(prec) == 1L && length(succ) == 1L) {
      periods$end[prec] <- periods$end[succ]
      periods <- periods[-succ, , drop = FALSE]
      row.names(periods) <- NULL
    }
  }

  periods
}
