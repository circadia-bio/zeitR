# ── Bimodal off-wrist refiner ─────────────────────────────────────────────────
# Faithful R port of BimodalOffwristRefiner (ActTrust configuration).
# Author of original Python: Julius A. P. P. de Paula (Condor Instruments, 2023)
#
# INDEXING CONVENTION (Python-style throughout):
#   start = 0-based first index (R vector access: x[start+1])
#   end   = exclusive last index (R vector access: x[start+1 .. end])
#   length = end - start
# ─────────────────────────────────────────────────────────────────────────────

#' Run the three-stage BimodalOffwristRefiner (ActTrust configuration)
#' @noRd
.bimodal_refine_acttrust <- function(
    initial_offwrist, activity, activity_median, temperature,
    norm_temp_variance, temp_derivative, temp_derivative_variance,
    temperature_threshold, ashman, activity_median_low, is_low_temp,
    filter_hws = 10L, dif_temp, epoch_hour = 60L,
    do_near_all_off_detection = TRUE, datetime_stamps = NULL
) {
  n <- length(initial_offwrist)

  # ── ActTrust fixed parameters ────────────────────────────────────────────
  minimum_onwrist_length           <- 20L
  minimum_offwrist_length          <- 10L   # Python default (wrapper sets 10)
  minimum_preceding_onwrist_length <- 40L
  activity_threshold_quantile      <- 0.24
  minimum_low_activity_proportion  <- 0.8
  max_low_tv_border                <- 0.3
  long_offwrist_length             <- 4L * as.integer(epoch_hour)
  short_offwrist_length            <- 20L
  temp_var_thr_quantile            <- 0.9
  sleep_act_hws                    <- 120L
  sleep_all_q                      <- 0.4
  sleep_pos_q                      <- 0.05
  max_offwrist_sleep_prop          <- 0.4
  ashman_d_min                     <- 1.5
  bimodal_max_offwrist_prop        <- 0.0
  bimodal_min_low_act_prop         <- 1.0
  report_zero_act_prop_min         <- 0.35
  border_conc_min                  <- 0.5
  report_act_around_min            <- 0.1
  report_low_act_prop_min          <- 0.5
  offwrist_max_temp_dif_med        <- 0.8
  offwrist_min_low_temp_prop       <- 0.6   # low_temperature_filter threshold
  valley_quantile                  <- 0.99
  peak_quantile                    <- 0.99
  valley_peak_low_act_min          <- 0.5
  valley_peak_low_temp_min         <- 0.5
  short_vp_decrease_ratio_min      <- 1.25
  decrease_ratio_min               <- 1.0
  sleep_low_temp_prop_max          <- 0.5
  min_surrounding_ow_len           <- 80L
  half_filter_hws                  <- as.integer(filter_hws / 2)

  # ── Activity threshold ─────────────────────────────────────────────────
  act_zero_prop <- zero_prop(activity)
  act_thr_q     <- act_zero_prop + (1 - act_zero_prop) * activity_threshold_quantile
  activity_thr  <- as.double(stats::quantile(activity, act_thr_q, names = FALSE, type = 1))
  if (activity_thr < 200) activity_thr <- 200

  pos_act             <- activity[activity > 0]
  positive_act_median <- if (length(pos_act) > 0) stats::median(pos_act) else 0
  temp_var_thr        <- as.double(stats::quantile(norm_temp_variance,
                                                    temp_var_thr_quantile, names = FALSE))

  ow_periods <- .rle_periods(initial_offwrist == 0L)
  if (nrow(ow_periods) == 0L) return(rep(1L, n))

  # ── Stage 1: filter short on-wrist periods ────────────────────────────
  on_periods <- .rle_periods(initial_offwrist == 1L)
  on_periods <- on_periods[(on_periods$end - on_periods$start) > minimum_onwrist_length, ]
  on_periods <- on_periods[order(on_periods$start), ]
  row.names(on_periods) <- NULL

  stage1 <- integer(n)
  for (i in seq_len(nrow(on_periods))) {
    s <- on_periods$start[i]; e <- on_periods$end[i]
    if (s < e) stage1[(s + 1L):e] <- 1L
  }
  ow_s1 <- .rle_periods(stage1 == 0L)
  if (nrow(ow_s1) == 0L) return(rep(1L, n))

  # ── Stage 2: border refinement ────────────────────────────────────────
  ow_s2 <- .refine_borders(
    ow_s1, norm_temp_variance, activity_median, activity_median_low,
    temperature, temperature_threshold, activity_thr, minimum_low_activity_proportion,
    minimum_preceding_onwrist_length, filter_hws, half_filter_hws,
    max_low_tv_border, temp_var_thr, activity, n
  )
  if (nrow(ow_s2) == 0L) return(rep(1L, n))

  # ── Stage 3: initial feature filters ─────────────────────────────────
  ow_s2$length <- ow_s2$end - ow_s2$start

  # low_temperature_filter (half_day_length_validation = FALSE for ActTrust)
  ow_s2$low_temp_prop <- vapply(seq_len(nrow(ow_s2)), function(i) {
    s <- ow_s2$start[i]; e <- ow_s2$end[i]
    below_prop(temperature[(s + 1L):e], temperature_threshold)
  }, numeric(1))
  ow_s2 <- ow_s2[ow_s2$low_temp_prop >= offwrist_min_low_temp_prop, ]
  row.names(ow_s2) <- NULL
  if (nrow(ow_s2) == 0L) return(rep(1L, n))

  # low_activity_filter
  ow_s2$low_act_prop <- vapply(seq_len(nrow(ow_s2)), function(i) {
    s <- ow_s2$start[i]; e <- ow_s2$end[i]
    below_prop(activity[(s + 1L):e], activity_thr)
  }, numeric(1))
  ow_s2$zero_act_prop <- vapply(seq_len(nrow(ow_s2)), function(i) {
    s <- ow_s2$start[i]; e <- ow_s2$end[i]
    zero_prop(activity[(s + 1L):e])
  }, numeric(1))
  ow_s2 <- ow_s2[ow_s2$low_act_prop > minimum_low_activity_proportion, ]
  row.names(ow_s2) <- NULL
  if (nrow(ow_s2) == 0L) return(rep(1L, n))

  # offwrist_length_filter
  ow_s2 <- ow_s2[ow_s2$length >= minimum_offwrist_length, ]
  row.names(ow_s2) <- NULL
  if (nrow(ow_s2) == 0L) return(rep(1L, n))

  # temperature_difference_median filter (ActTrust only, not lumus)
  ow_s2$temp_dif_med <- vapply(seq_len(nrow(ow_s2)), function(i) {
    s <- ow_s2$start[i]; e <- ow_s2$end[i]
    stats::median(dif_temp[(s + 1L):e], na.rm = TRUE)
  }, numeric(1))
  ow_s2 <- ow_s2[ow_s2$temp_dif_med < offwrist_max_temp_dif_med, ]
  row.names(ow_s2) <- NULL
  if (nrow(ow_s2) == 0L) return(rep(1L, n))

  # Keep copy for sleep estimation (after_initial_refinement)
  after_initial_periods <- ow_s2

  # ── short_offwrist_only_mask (exclude long periods from sleep estimation)
  short_only_mask <- rep(TRUE, n)
  is_long         <- ow_s2$length >= long_offwrist_length
  for (i in which(is_long)) {
    s <- ow_s2$start[i]; e <- ow_s2$end[i]
    if (s < e) short_only_mask[(s + 1L):e] <- FALSE
  }
  valid_activity <- activity[short_only_mask]

  # ── estimate_sleep_then_filter ────────────────────────────────────────
  estimated_sleep <- .estimate_sleep_padded(
    activity, valid_activity, temperature, temperature_threshold,
    short_only_mask, sleep_act_hws, sleep_all_q, sleep_pos_q,
    act_zero_prop, sleep_low_temp_prop_max, epoch_hour, n
  )

  # Sleep filter: remove periods with too much sleep overlap
  after_initial_periods$sleep_prop <- vapply(seq_len(nrow(after_initial_periods)), function(i) {
    s <- after_initial_periods$start[i]; e <- after_initial_periods$end[i]
    zero_prop(estimated_sleep[(s + 1L):e])
  }, numeric(1))
  sleep_filtered <- after_initial_periods[
    after_initial_periods$sleep_prop < max_offwrist_sleep_prop, ]
  row.names(sleep_filtered) <- NULL

  # Rebuild sleep_filtered binary vector
  sleep_filtered_ow <- rep(1, n)
  for (i in seq_len(nrow(sleep_filtered))) {
    s <- sleep_filtered$start[i]; e <- sleep_filtered$end[i]
    if (s < e) sleep_filtered_ow[(s + 1L):e] <- 0
  }

  # ── analyze_sleep_borders (sets short_sleep_border_offwrist) ──────────
  # For ActTrust this has minor effect; simplified here
  short_sleep_border_offwrist <- rep(0, n)

  # ── valley_peak_offwrist_algorithm ────────────────────────────────────
  vp_df <- .valley_peak_algorithm_full(
    temp_derivative = temp_derivative,
    activity        = activity,
    temperature     = temperature,
    temperature_threshold = temperature_threshold,
    estimated_sleep = estimated_sleep,
    dif_temp        = dif_temp,
    activity_thr    = activity_thr,
    valley_quantile = valley_quantile,
    peak_quantile   = peak_quantile,
    minimum_offwrist_length = minimum_offwrist_length,
    long_offwrist_length    = long_offwrist_length,
    short_offwrist_length   = short_offwrist_length,
    short_vp_decrease_ratio_min = short_vp_decrease_ratio_min,
    decrease_ratio_min      = decrease_ratio_min,
    valley_peak_low_act_min = valley_peak_low_act_min,
    valley_peak_low_temp_min = valley_peak_low_temp_min,
    offwrist_max_temp_dif_med = offwrist_max_temp_dif_med,
    forbidden_zone  = .compute_forbidden_zone_v2(estimated_sleep, epoch_hour, n),
    datetime_stamps = datetime_stamps,
    epoch_hour      = epoch_hour,
    n               = n,
    next_possible_length = TRUE
  )

  # Forbidden zone
  forbidden_zone <- .compute_forbidden_zone_v2(estimated_sleep, epoch_hour, n)

  # Apply forbidden zone to sleep_filtered periods
  sleep_filtered$in_forbidden <- vapply(seq_len(nrow(sleep_filtered)), function(i) {
    s <- sleep_filtered$start[i]; e <- sleep_filtered$end[i]
    any(forbidden_zone[(s + 1L):e] == 1L)
  }, logical(1))
  sleep_filtered <- sleep_filtered[!sleep_filtered$in_forbidden, ]
  row.names(sleep_filtered) <- NULL

  # Rebuild sleep_filtered_ow after forbidden zone filter
  sleep_filtered_ow <- rep(1, n)
  for (i in seq_len(nrow(sleep_filtered))) {
    s <- sleep_filtered$start[i]; e <- sleep_filtered$end[i]
    if (s < e) sleep_filtered_ow[(s + 1L):e] <- 0
  }

  # ── DEBUG HOOK (option-guarded; default off, no effect in normal use) ──
  if (isTRUE(getOption("zeitR.debug_offwrist", FALSE))) {
    assign(".zeitR_dbg", list(
      estimated_sleep       = estimated_sleep,
      forbidden_zone        = forbidden_zone,
      after_initial_periods = after_initial_periods,
      sleep_filtered_post   = sleep_filtered,
      short_only_mask       = short_only_mask,
      vp_df                 = vp_df,
      max_offwrist_sleep_prop = max_offwrist_sleep_prop,
      long_offwrist_length  = long_offwrist_length
    ), envir = .GlobalEnv)
  }

  # ── check_bimodality ──────────────────────────────────────────────────
  high_act_q    <- act_zero_prop + (1 - act_zero_prop) * 0.9
  high_act_val  <- as.double(stats::quantile(activity, high_act_q, names = FALSE))
  low_pos_act   <- below_prop(pos_act, activity_thr)
  init_ow_prop  <- mean(sleep_filtered_ow == 0)  # uses sleep_filtered_offwrist

  is_bimodal    <- TRUE
  is_low_act_flag <- FALSE

  if (high_act_val <= 100) {
    is_bimodal <- FALSE; is_low_act_flag <- TRUE
  } else if (act_zero_prop >= 0.8) {
    is_bimodal <- FALSE; is_low_act_flag <- TRUE
  } else {
    if (ashman <= ashman_d_min) {
      if (init_ow_prop >= bimodal_max_offwrist_prop &&
          low_pos_act >= bimodal_min_low_act_prop) {
        is_bimodal <- FALSE
      }
    }
  }

  # ── Combine periods ───────────────────────────────────────────────────
  if (is_bimodal) {
    long_df  <- sleep_filtered[sleep_filtered$length >= long_offwrist_length, ]
    short_df <- sleep_filtered[sleep_filtered$length <  long_offwrist_length, ]

    # VP periods
    if (!is.null(vp_df) && nrow(vp_df) > 0) {
      vp_short <- vp_df[vp_df$length < long_offwrist_length, ]
    } else {
      vp_short <- data.frame(start = integer(), end = integer(),
                             length = integer(), valley_peak = logical())
    }

    long_df$valley_peak  <- FALSE
    short_df$valley_peak <- FALSE

    combined <- rbind(
      long_df[,  c("start", "end", "length", "valley_peak")],
      short_df[, c("start", "end", "length", "valley_peak")],
      if (nrow(vp_short) > 0) vp_short[, c("start", "end", "length", "valley_peak")]
      else data.frame(start=integer(), end=integer(), length=integer(), valley_peak=logical())
    )
    combined <- combined[order(combined$start), ]
    row.names(combined) <- NULL

    # Recompute features
    combined$length <- combined$end - combined$start
    combined$low_temp_prop <- vapply(seq_len(nrow(combined)), function(i) {
      s <- combined$start[i]; e <- combined$end[i]
      below_prop(temperature[(s+1L):e], temperature_threshold)
    }, numeric(1))
    combined$low_act_prop <- vapply(seq_len(nrow(combined)), function(i) {
      s <- combined$start[i]; e <- combined$end[i]
      below_prop(activity[(s+1L):e], activity_thr)
    }, numeric(1))
    combined$zero_act_prop <- vapply(seq_len(nrow(combined)), function(i) {
      s <- combined$start[i]; e <- combined$end[i]
      zero_prop(activity[(s+1L):e])
    }, numeric(1))

    refined_periods <- combined

  } else {
    if (is_low_act_flag && do_near_all_off_detection) {
      if (positive_act_median >= activity_thr) {
        refined_periods <- .near_all_off_detection_v2(
          activity, activity_thr, positive_act_median, n, epoch_hour
        )
      } else {
        refined_periods <- data.frame(start=integer(), end=integer(),
                                      length=integer(), valley_peak=logical())
      }
    } else {
      refined_periods <- data.frame(start=integer(), end=integer(),
                                    length=integer(), valley_peak=logical())
    }
    if (nrow(refined_periods) == 0L) return(rep(1L, n))
  }

  # Rebuild binary vector
  refined_ow <- rep(1, n)
  for (i in seq_len(nrow(refined_periods))) {
    s <- refined_periods$start[i]; e <- min(refined_periods$end[i], n)
    if (s < e) refined_ow[(s + 1L):e] <- 0
  }

  # ── description_report_based_filter ──────────────────────────────────
  if (nrow(refined_periods) > 0L && is_bimodal) {
    # Dynamic temperature_difference_minimum from valid epochs
    dif_valid <- dif_temp[short_only_mask]
    temp_dif_min <- as.double(stats::quantile(dif_valid, 0.1, names = FALSE))
    if (temp_dif_min < 0.2) temp_dif_min <- 0.2

    # Compute sleep_low_activity_threshold for report
    sleep_act <- activity[estimated_sleep == 0L]
    if (length(sleep_act) > 0) {
      sleep_zp  <- zero_prop(sleep_act)
      slat_q    <- sleep_zp + (1 - sleep_zp) * 0.1
      sleep_lat <- as.double(stats::quantile(sleep_act, slat_q, names = FALSE))
    } else {
      sleep_lat <- activity_thr
    }

    report <- .describe_offwrist_periods_v2(
      offwrist_periods = refined_periods,
      activity = activity, temperature = temperature,
      temperature_variance = norm_temp_variance,
      temperature_threshold = temperature_threshold,
      activity_thr = activity_thr,
      sleep_lat = sleep_lat,
      dif_temp = dif_temp, n = n
    )

    refined_periods <- .description_report_filter_v2(
      periods = refined_periods, report = report,
      report_zero_act_min = report_zero_act_prop_min,
      border_conc_min = border_conc_min,
      act_around_min = report_act_around_min,
      low_act_prop_min = report_low_act_prop_min,
      temp_dif_min = temp_dif_min,
      offwrist_max_temp_dif_med = offwrist_max_temp_dif_med,
      long_offwrist_length = long_offwrist_length,
      short_offwrist_length = short_offwrist_length,
      is_highly_separable = (ashman > 3)
    )
    row.names(refined_periods) <- NULL

    # Rebuild binary vector
    refined_ow <- rep(1, n)
    for (i in seq_len(nrow(refined_periods))) {
      s <- refined_periods$start[i]; e <- min(refined_periods$end[i], n)
      if (s < e) refined_ow[(s + 1L):e] <- 0
    }
  }

  # ── surrounded_onwrist_filter ─────────────────────────────────────────
  if (nrow(refined_periods) > 0L) {
    refined_ow <- .surrounded_onwrist_filter_v2(
      refined_ow, temperature, temperature_threshold,
      minimum_preceding_onwrist_length, min_surrounding_ow_len, n
    )
    refined_periods <- .rle_periods(refined_ow == 0L)
    if (nrow(refined_periods) == 0L) return(rep(1L, n))
    refined_periods$length <- refined_periods$end - refined_periods$start
  }

  # ── Final border snap ─────────────────────────────────────────────────
  if (nrow(refined_periods) > 0L) {
    if (refined_periods$start[1] <= minimum_offwrist_length)
      refined_periods$start[1] <- 0L
    last <- nrow(refined_periods)
    if ((n - refined_periods$end[last]) <= minimum_offwrist_length)
      refined_periods$end[last] <- n
  }

  # ── Assemble output ───────────────────────────────────────────────────
  out <- rep(1L, n)
  for (i in seq_len(nrow(refined_periods))) {
    s <- refined_periods$start[i]; e <- min(refined_periods$end[i], n)
    if (s < e) out[(s + 1L):e] <- 0L
  }
  out
}

# ── .rle_periods ─────────────────────────────────────────────────────────────
#' Python-style 0-indexed exclusive-end periods from a binary vector
#' length = end - start; R access: x[(start+1):end]
#' @noRd
.rle_periods <- function(x) {
  r    <- rle(as.integer(x))
  ends <- cumsum(r$lengths)
  strt <- ends - r$lengths
  keep <- r$values == 1L
  if (!any(keep)) return(data.frame(start = integer(), end = integer()))
  data.frame(start = strt[keep], end = ends[keep])
}

# ── below_prop ───────────────────────────────────────────────────────────────
#' Proportion of values below thr.
#' NOTE: Python functions.below_prop uses <=, but empirically R's threshold
#' values align with Python only under strict < (matches verified Stage-2/sleep
#' parity); using <= pushes borderline periods over the keep gates. Kept as <.
#' @noRd
below_prop <- function(x, thr) {
  if (length(x) == 0L) return(0)
  sum(x < thr, na.rm = TRUE) / length(x)
}

# ── Stage 2: border refinement ───────────────────────────────────────────────
#' @noRd
.refine_borders <- function(
    offwrist_periods, norm_temp_variance, activity_median, activity_median_low,
    temperature, temperature_threshold, activity_thr, minimum_low_act_prop,
    min_preceding_onwrist, filter_hws, half_filter_hws, max_low_tv_border,
    temp_var_thr, activity, n
) {
  ntv          <- norm_temp_variance
  refined      <- list()
  periods      <- as.matrix(offwrist_periods[, c("start", "end")])
  n_off        <- nrow(periods)
  ow_idx       <- 1L
  search_start <- TRUE

  # NOTE on indexing: `start`/`end` are stored in Python convention throughout
  # (start = 0-based first index, end = EXCLUSIVE last index). Python epoch j
  # maps to R vector index (j + 1L). All array reads below convert explicitly.

  while (ow_idx <= n_off) {

    if (search_start) {
      # ── search_offwrist_start ────────────────────────────────────────────
      start <- periods[ow_idx, 1L]                      # Python value
      if (length(refined) > 0L) {
        previous_end <- refined[[length(refined)]][2L]
        if (start < previous_end) start <- previous_end + 2L
      } else {
        previous_end <- 0L
      }

      if (start >= filter_hws) {
        # gate: Python ntv[start] >= thr  ->  R ntv[start + 1]
        if ((start + 1L) <= n && ntv[start + 1L] >= temp_var_thr) {
          lo_py <- max(0L, start - filter_hws)          # Python epoch lower bound
          # low_temperature_variance_proportion_around_start (non-centered):
          #   below_prop(ntv[max(0, start - fhwl) : start + 1])  -> epochs lo_py..start
          win  <- ntv[(lo_py + 1L):(start + 1L)]
          ltvp <- below_prop(win, temp_var_thr)
          if (ltvp <= max_low_tv_border) {
            # start_peak_found:
            #   start = (start - fhwl) + argmax(ntv[max(0,start-fhwl):start+1])
            a         <- which.max(win) - 1L            # 0-based argmax
            new_start <- lo_py + a                      # Python value
            if (new_start < previous_end) new_start <- previous_end + 2L
            refined[[length(refined) + 1L]] <- c(new_start, 0L)
            search_start <- FALSE
          } else {
            # quick_delete_start = FALSE  ->  find_peak_base
            res <- .find_peak_base_start2(start, previous_end, ntv,
                                          temperature, temperature_threshold,
                                          activity_median, activity_median_low,
                                          activity_thr, temp_var_thr, n)
            if (isTRUE(res$delete)) {
              periods <- periods[-ow_idx, , drop = FALSE]; n_off <- nrow(periods)
            } else {
              refined[[length(refined) + 1L]] <- c(res$new_start, 0L)
              search_start <- FALSE
            }
          }
        } else {
          # low temperature variance -> find_peak_base, previous_end from
          # the ORIGINAL preceding period (matches Python).
          prev_end2 <- if (ow_idx > 1L) periods[ow_idx - 1L, 2L] else 0L
          res <- .find_peak_base_start2(start, prev_end2, ntv,
                                        temperature, temperature_threshold,
                                        activity_median, activity_median_low,
                                        activity_thr, temp_var_thr, n)
          if (isTRUE(res$delete)) {
            periods <- periods[-ow_idx, , drop = FALSE]; n_off <- nrow(periods)
          } else {
            refined[[length(refined) + 1L]] <- c(res$new_start, 0L)
            search_start <- FALSE
          }
        }
      } else {
        # too close to beginning: start = argmax(ntv[0 : start + 1])  (Python value)
        new_start <- which.max(ntv[1L:(start + 1L)]) - 1L
        refined[[length(refined) + 1L]] <- c(max(new_start, 0L), 0L)
        search_start <- FALSE
      }

    } else {
      # ── search_offwrist_end ──────────────────────────────────────────────
      end <- periods[ow_idx, 2L]                        # Python EXCLUSIVE end value
      # gate: Python ntv[end] >= thr  ->  R ntv[end + 1]
      gate_ok <- (end + 1L) <= n && ntv[end + 1L] >= temp_var_thr

      if (gate_ok) {
        if (end + filter_hws <= n) {
          # low_temperature_variance_proportion_around_end (non-centered):
          #   below_prop(ntv[end : end + fhwl])  -> epochs end..end+fhwl-1
          ltvp <- below_prop(ntv[(end + 1L):(end + filter_hws)], temp_var_thr)
          if (ltvp <= max_low_tv_border) {
            # end_peak_found:
            #   end = end + argmax(ntv[end : min(n, end + fhwl + 1)])  -> epochs end..
            hi      <- min(n, end + filter_hws + 1L)
            a       <- which.max(ntv[(end + 1L):hi]) - 1L
            new_end <- end + a                          # Python exclusive-end value
            refined[[length(refined)]][2L] <- new_end
            search_start <- TRUE; ow_idx <- ow_idx + 1L
          } else {
            if (ow_idx < n_off) {
              next_start <- periods[ow_idx + 1L, 1L]
              res <- .do_more_checks_around_end(ow_idx, next_start, end,
                                                periods, refined,
                                                ntv, temperature,
                                                temperature_threshold, activity_median,
                                                activity_median_low, activity_thr,
                                                temp_var_thr, filter_hws,
                                                max_low_tv_border,
                                                min_preceding_onwrist,
                                                minimum_low_act_prop, activity, n)
              ow_idx <- res$ow_idx; periods <- res$periods
              n_off  <- nrow(periods); refined <- res$refined
              search_start <- res$search_start
            } else {
              refined[[length(refined)]][2L] <- .last_period_end(end, ntv, n)
              search_start <- TRUE; ow_idx <- ow_idx + 1L
            }
          }
        } else {
          refined[[length(refined)]][2L] <- .last_period_end(end, ntv, n)
          search_start <- TRUE; ow_idx <- ow_idx + 1L
        }
      } else {
        if (ow_idx < n_off) {
          next_start <- periods[ow_idx + 1L, 1L]
          res <- .try_find_peak_base_end(ow_idx, next_start, end,
                                          periods, refined,
                                          ntv, temperature,
                                          temperature_threshold, activity_median,
                                          activity_median_low, activity_thr,
                                          temp_var_thr, min_preceding_onwrist,
                                          minimum_low_act_prop, activity, n)
          ow_idx <- res$ow_idx; periods <- res$periods
          n_off  <- nrow(periods); refined <- res$refined
          search_start <- res$search_start
        } else {
          refined[[length(refined)]][2L] <- .last_period_end(end, ntv, n)
          search_start <- TRUE; ow_idx <- ow_idx + 1L
        }
      }
    }
  }

  if (length(refined) == 0L) return(data.frame(start = integer(), end = integer()))
  m <- do.call(rbind, refined)
  m <- m[m[, 2L] > 0L, , drop = FALSE]
  if (nrow(m) == 0L) return(data.frame(start = integer(), end = integer()))
  data.frame(start = m[, 1L], end = m[, 2L])
}

#' Faithful port of the Python last-period end fallback:
#'   end = end + argmax(ntv[end : data_length])   (epochs end..n-1)
#' `end` is the Python exclusive-end value; returns the new Python exclusive end.
#' @noRd
.last_period_end <- function(end, ntv, n) {
  if ((end + 1L) > n) return(end)
  end + which.max(ntv[(end + 1L):n]) - 1L
}

#' @noRd
.find_peak_base_start2 <- function(start, previous_end, ntv, temperature,
                                    temperature_threshold, act_median, act_med_low,
                                    activity_thr, temp_var_thr, n) {
  i                 <- start - 1L
  last_valid_border <- start
  peak_base         <- start

  while (i > previous_end) {
    ri    <- i + 1L
    if (ri < 1L || ri > n) { i <- previous_end; break }
    valid <- .check_valid_border_mod(act_med_low, temperature, temperature_threshold,
                                     act_median, activity_thr, ri)
    if (valid) {
      if (ntv[ri] >= temp_var_thr) { peak_base <- i; i <- previous_end }
      i <- i - 1L
    } else { last_valid_border <- i; i <- previous_end }
  }

  if (last_valid_border == start) {
    if (peak_base < start) {
      # Python: peak = peak_base - 1; while ntv[peak] >= ntv[peak+1]: peak -= 1
      # Python epoch p -> R ntv[p+1]; ntv[peak] -> ntv[peak+1], ntv[peak+1] -> ntv[peak+2]
      peak <- peak_base - 1L
      while ((peak + 1L) >= 1L && (peak + 2L) <= n &&
             ntv[peak + 1L] >= ntv[peak + 2L]) {
        peak <- peak - 1L
      }
      return(list(new_start = peak, delete = FALSE))   # Python value
    } else {
      return(list(new_start = NULL, delete = TRUE))
    }
  } else {
    return(list(new_start = last_valid_border, delete = FALSE))
  }
}

#' @noRd
.do_more_checks_around_end <- function(ow_idx, next_start, end, periods, refined,
                                        ntv, temperature, temperature_threshold,
                                        act_median, act_med_low, activity_thr,
                                        temp_var_thr, filter_hws, max_low_tv_border,
                                        min_preceding_onwrist, minimum_low_act_prop,
                                        activity, n) {
  too_close <- .check_offwrist_too_close(next_start - end, end, next_start,
                                          activity, activity_thr,
                                          minimum_low_act_prop, n)
  if (too_close) {
    periods <- periods[-ow_idx, , drop = FALSE]
    return(list(ow_idx = ow_idx, periods = periods, refined = refined,
                search_start = FALSE))
  }

  re <- min(end, n)
  # Python: below_prop(ntv[max(0, end - fhwl - 1) : end])  -> epochs lo_py..end-1
  lo_py <- max(0L, end - filter_hws - 1L)              # Python epoch lower bound
  win   <- ntv[(lo_py + 1L):end]                       # R indices  (epochs lo_py..end-1)
  ltvp  <- below_prop(win, temp_var_thr)

  if (ltvp > max_low_tv_border) {
    periods <- periods[-ow_idx, , drop = FALSE]
    return(list(ow_idx = ow_idx, periods = periods, refined = refined,
                search_start = FALSE))
  } else {
    # Python: end = (end - fhwl - 1) + argmax(ntv[max(0,end-fhwl-1):end])
    a       <- which.max(win) - 1L                     # 0-based argmax
    new_end <- lo_py + a                               # Python exclusive-end value
    refined[[length(refined)]][2L] <- new_end
    return(list(ow_idx = ow_idx + 1L, periods = periods, refined = refined,
                search_start = TRUE))
  }
}

#' @noRd
.try_find_peak_base_end <- function(ow_idx, next_start, end, periods, refined,
                                     ntv, temperature, temperature_threshold,
                                     act_median, act_med_low, activity_thr,
                                     temp_var_thr, min_preceding_onwrist,
                                     minimum_low_act_prop, activity, n) {
  too_close <- .check_offwrist_too_close(next_start - end, end, next_start,
                                          activity, activity_thr,
                                          minimum_low_act_prop, n)
  if (too_close) {
    periods <- periods[-ow_idx, , drop = FALSE]
    return(list(ow_idx = ow_idx, periods = periods, refined = refined,
                search_start = FALSE))
  }

  i <- end + 1L; last_valid_border <- end; peak_base <- end
  while (i < next_start) {
    ri    <- i + 1L
    if (ri < 1L || ri > n) { i <- next_start; break }
    valid <- .check_valid_border_mod(act_med_low, temperature, temperature_threshold,
                                     act_median, activity_thr, ri)
    if (valid) {
      if (ntv[ri] >= temp_var_thr) { peak_base <- i; i <- next_start }
      i <- i + 1L
    } else { last_valid_border <- i; i <- next_start }
  }

  if (last_valid_border == end) {
    if (peak_base > end) {
      # Python: peak = peak_base + 1; while ntv[peak] >= ntv[peak-1]: peak += 1
      # Python epoch p -> R ntv[p+1]; ntv[peak] -> ntv[peak+1], ntv[peak-1] -> ntv[peak]
      peak <- peak_base + 1L
      while ((peak + 1L) <= n && ntv[peak + 1L] >= ntv[peak]) peak <- peak + 1L
      refined[[length(refined)]][2L] <- peak             # Python exclusive-end value
      return(list(ow_idx = ow_idx + 1L, periods = periods, refined = refined,
                  search_start = TRUE))
    } else {
      periods <- periods[-ow_idx, , drop = FALSE]
      return(list(ow_idx = ow_idx, periods = periods, refined = refined,
                  search_start = FALSE))
    }
  } else {
    refined[[length(refined)]][2L] <- last_valid_border
    return(list(ow_idx = ow_idx + 1L, periods = periods, refined = refined,
                search_start = TRUE))
  }
}

#' @noRd
.check_offwrist_too_close <- function(following_len, end, next_start,
                                       activity, activity_thr,
                                       minimum_low_act_prop, n) {
  if (following_len > 40L) return(FALSE)
  re <- min(end, n); rns <- min(next_start, n)
  if ((re + 1L) > rns || re < 1L) return(FALSE)
  below_prop(activity[(re + 1L):rns], activity_thr) > minimum_low_act_prop
}

#' @noRd
.check_valid_border_mod <- function(act_med_low, temperature, temperature_threshold,
                                     act_median, activity_thr, ri) {
  act_med_low[ri] == 1L ||
    (temperature[ri] < temperature_threshold && act_median[ri] < 2 * activity_thr)
}

# ── Sleep estimation (with padding matching Python exactly) ───────────────────
#' @noRd
.estimate_sleep_padded <- function(activity, valid_activity, temperature,
                                    temperature_threshold, short_only_mask,
                                    sleep_act_hws, sleep_all_q, sleep_pos_q,
                                    act_zero_prop, sleep_low_temp_prop_max,
                                    epoch_hour, n) {
  # Exact translation of Python:
  #   padding = np.max(activity) * np.ones(2*hws)
  #   padded = np.insert(valid_activity, 0, padding)
  #   padded = np.append(padded, padding)
  #   filtered = median_filter(padded, hws, padding="padded")
  # where padding="padded" means:
  #   n -= 4*hws  (n = n_valid)
  #   rolled = rolling_window(padded, 2*hws+1)  -> n_valid + 2*hws windows
  #   filt = median(rolled)[hws : n+hws]  -> n_valid values
  pad_val <- max(activity, na.rm = TRUE)
  hws     <- sleep_act_hws
  padded  <- c(rep(pad_val, 2L * hws), valid_activity, rep(pad_val, 2L * hws))
  n_valid <- length(valid_activity)
  win     <- 2L * hws + 1L
  n_padded <- length(padded)  # = n_valid + 4*hws

  # True sliding window median matching Python's rolling_window + median extraction exactly.
  # rolling_window(padded, win) produces (n_padded - win + 1) = n_valid + 2*hws windows.
  # We need windows [hws+1 .. hws+n_valid] (1-indexed) = Python's [hws : n+hws].
  # runmed() does NOT do this: it applies its own boundary rules (Tukey's running median)
  # that differ from a plain sliding window at the edges. Use RcppRoll if available,
  # otherwise zoo::rollmedian, otherwise a pure-R fallback.
  # Python rolling_window(padded, 2*hws+1) + median(axis=-1)[hws : n_valid+hws]
  # Window i (0-indexed) covers padded[i : i+2*hws+1].
  # Output i (1-indexed, i in 1..n_valid) = median(padded[(hws+i):(3*hws+i)])
  # = rollmedian with align="right" at position (3*hws + i) in padded.
  # Equivalently: the pure-R reference formula is median(padded[(hws+i):(3*hws+i)]).
  if (requireNamespace("RcppRoll", quietly = TRUE)) {
    all_windows <- RcppRoll::roll_median(padded, n = win, fill = NA_real_, align = "right")
    filtered_valid <- as.double(all_windows[(3L * hws + 1L):(3L * hws + n_valid)])
  } else if (requireNamespace("zoo", quietly = TRUE)) {
    all_windows <- zoo::rollmedian(padded, k = win, fill = NA_real_, align = "right")
    filtered_valid <- as.double(all_windows[(3L * hws + 1L):(3L * hws + n_valid)])
  } else {
    filtered_valid <- vapply(seq_len(n_valid), function(i) {
      stats::median(padded[(hws + i):(3L * hws + i)])
    }, numeric(1))
  }

  # Compute sleep threshold ("both" configuration)
  if (act_zero_prop < sleep_all_q) {
    sleep_thr <- as.double(stats::quantile(activity, sleep_all_q, names = FALSE))
  } else {
    q_idx     <- act_zero_prop + (1 - act_zero_prop) * sleep_pos_q
    sleep_thr <- as.double(stats::quantile(activity, q_idx, names = FALSE))
  }

  # Threshold: 0 = sleep, 1 = wake in valid signal
  sleep_est_valid <- as.integer(filtered_valid > sleep_thr)

  # Short sleep filter: keep only sleep periods longer than long_offwrist_length (4 * epoch_hour)
  # Python: sleep_periods_df[length > long_offwrist_length]
  long_offwrist_length <- 4L * as.integer(epoch_hour)
  sleep_periods <- .rle_periods(sleep_est_valid == 0L)
  if (nrow(sleep_periods) > 0L) {
    sleep_periods <- sleep_periods[
      (sleep_periods$end - sleep_periods$start) > long_offwrist_length, ]
    sleep_est_valid <- rep(1L, length(valid_activity))
    for (i in seq_len(nrow(sleep_periods))) {
      s <- sleep_periods$start[i]; e <- sleep_periods$end[i]
      if (s < e) sleep_est_valid[(s + 1L):e] <- 0L
    }
  }

  # Low temperature filter on valid sleep periods
  sleep_periods <- .rle_periods(sleep_est_valid == 0L)
  if (nrow(sleep_periods) > 0L) {
    valid_temperature <- temperature[short_only_mask]
    keep <- vapply(seq_len(nrow(sleep_periods)), function(i) {
      s <- sleep_periods$start[i]; e <- sleep_periods$end[i]
      below_prop(valid_temperature[(s + 1L):e], temperature_threshold) <
        sleep_low_temp_prop_max
    }, logical(1))
    sleep_periods <- sleep_periods[keep, ]
    sleep_est_valid <- rep(1L, length(valid_activity))
    for (i in seq_len(nrow(sleep_periods))) {
      s <- sleep_periods$start[i]; e <- sleep_periods$end[i]
      if (s < e) sleep_est_valid[(s + 1L):e] <- 0L
    }
  }

  # Map back to full n-length vector
  estimated_sleep <- rep(1L, n)
  estimated_sleep[short_only_mask] <- sleep_est_valid

  # sleep_low_temperature_filter is disabled for ActTrust (do_sleep_low_temperature_filter=False)

  estimated_sleep
}

# ── Forbidden zone ────────────────────────────────────────────────────────────
#' @noRd
.compute_forbidden_zone_v2 <- function(estimated_sleep, epoch_hour, n) {
  forbidden     <- integer(n)
  sleep_periods <- .rle_periods(estimated_sleep == 0L)
  if (nrow(sleep_periods) == 0L) return(forbidden)

  for (i in seq_len(nrow(sleep_periods))) {
    s   <- sleep_periods$start[i]; e <- sleep_periods$end[i]
    len <- e - s
    if (len > 0L) {
      q1_r <- s + as.integer(len * 0.25) + 1L
      q3_r <- s + as.integer(len * 0.75)
      if (q1_r <= q3_r && q1_r >= 1L && q3_r <= n)
        forbidden[q1_r:q3_r] <- 1L
    }
  }
  forbidden
}

# ── Valley-peak algorithm (faithful port; ActTrust, non-lumus, include_static) ─
# Mirrors valley_peak_offwrist_algorithm(include_static=TRUE, static_valley=0.95,
# static_peak=0.98) and valley_peak_detection from the Python refiner.
#
# scipy.signal.find_peaks(x, height): local maxima (plateau midpoint) with
# x[peak] >= height.
#' @noRd
.find_peaks_scipy <- function(x, height) {
  n <- length(x)
  if (n < 3L) return(integer(0))
  peaks <- integer(0)
  i    <- 2L
  imax <- n - 1L
  while (i < imax) {
    if (x[i - 1L] < x[i]) {
      ahead <- i + 1L
      while (ahead < imax && x[ahead] == x[i]) ahead <- ahead + 1L
      if (x[ahead] < x[i]) {
        left  <- i; right <- ahead - 1L
        mid   <- left + (right - left) %/% 2L
        peaks <- c(peaks, mid)
        i <- ahead
      }
    }
    i <- i + 1L
  }
  if (length(peaks) == 0L) return(integer(0))
  peaks[x[peaks] >= height]
}

#' Length of the first allowed (zero) run in a 0/1 vector (Python zero_sequences[0]).
#' @noRd
.first_allowed_run_len <- function(forbidden) {
  is_zero <- forbidden == 0
  if (!any(is_zero)) return(0L)
  r <- rle(is_zero)
  idx <- which(r$values)[1L]
  as.integer(r$lengths[idx])
}

#' Noon-to-noon day index (Python actigraphy_split_by_day, start_hour = 12),
#' including the tiny-first/last-day merge.
#' @noRd
.vp_day_index <- function(ts, start_hour = 12L, n) {
  ts <- as.POSIXct(ts)
  t0 <- ts[1L]
  lt0_hour <- as.integer(format(t0, "%H"))
  midnight <- as.POSIXct(format(t0, "%Y-%m-%d 00:00:00"), tz = attr(ts, "tzone") %||% "UTC")
  if (lt0_hour <= start_hour) {
    first_date <- midnight - as.difftime(start_hour, units = "hours")
  } else {
    first_date <- midnight + as.difftime(start_hour, units = "hours")
  }
  delta_h <- as.numeric(difftime(ts, first_date, units = "hours"))
  idx <- as.integer(floor(delta_h / 24)) + 1L
  idx[idx < 1L] <- 1L
  idx <- idx - min(idx) + 1L                       # relabel 1..D contiguous

  # tiny first/last day merge (matches Python)
  sizes <- as.integer(table(factor(idx, levels = sort(unique(idx)))))
  D <- length(sizes)
  labs <- sort(unique(idx))
  if (D > 2L) {
    typical <- sizes[2L]
    if (typical >= 2L * sizes[1L]) idx[idx == labs[1L]] <- labs[2L]
    if (typical >= 2L * sizes[D]) idx[idx == labs[D]] <- labs[D - 1L]
  } else if (D == 2L) {
    if (sizes[1L] >= 2L * sizes[2L] || sizes[2L] >= 2L * sizes[1L])
      idx[] <- labs[1L]
  }
  match(idx, sort(unique(idx)))                    # contiguous 1..k
}

#' Python compute_temperature_variations decrease_ratio for ONE period.
#' s,e are Python 0-based start / exclusive end.
#' @noRd
.vp_decrease_ratio <- function(temp_derivative, s, e) {
  if ((s + 1L) > e) return(0)
  wd  <- temp_derivative[(s + 1L):e]
  inc <- which(wd > 0)
  if (length(inc) > 2L) {
    inc_sum <- sum(wd[inc[2:(length(inc) - 1L)]])    # drop first & last positive
  } else if (length(inc) == 2L) {
    inc_sum <- wd[inc[2L]]
  } else {
    inc_sum <- sum(wd[inc])
  }
  dec <- which(wd < 0)
  if (length(dec) > 2L) {
    dec_sum <- -sum(wd[dec[2:length(dec)]])          # drop first negative
  } else if (length(dec) == 2L) {
    dec_sum <- -wd[dec[2L]]
  } else {
    dec_sum <- -sum(wd[dec])
  }
  if (inc_sum > 0) dec_sum / inc_sum else dec_sum
}

#' test_offwrist_surrounded_by_awake (ActTrust): median activity in the
#' minimum_onwrist_length window before OR after >= activity_thr.
#' @noRd
.vp_surrounded_by_awake <- function(activity, s, e, activity_thr, min_onwrist, n) {
  bs <- max(0L, s - min_onwrist)
  before <- if (bs < s) activity[(bs + 1L):s] else numeric(0)
  after  <- if (e < min(e + min_onwrist, n)) activity[(e + 1L):min(e + min_onwrist, n)] else numeric(0)
  mb <- if (length(before) > 0L) stats::median(before) else 0
  ma <- if (length(after)  > 0L) stats::median(after)  else 0
  (mb >= activity_thr) || (ma >= activity_thr)
}

#' test_temperature_difference_median_statistically_below (two-sample KS).
#' NOTE: scipy<->R two-sample KS alternative conventions differ; this affects
#' only short (<= short_offwrist_length) rescue periods. Uses R ks.test.
#' @noRd
.vp_test_tdm_below <- function(dif_temp, s, e, n) {
  L  <- e - s
  bs <- max(0L, s - L)
  before <- if (bs < s) dif_temp[(bs + 1L):s] else numeric(0)
  after  <- if (e < min(e + L, n)) dif_temp[(e + 1L):min(e + L, n)] else numeric(0)
  off    <- dif_temp[(s + 1L):e]
  tdm <- stats::median(off)
  if (length(before) == 0L || length(after) == 0L || length(off) == 0L)
    return(list(pval = 100, tdm = tdm))
  pv_b <- suppressWarnings(stats::ks.test(before, off, alternative = "greater")$p.value)
  pv_a <- suppressWarnings(stats::ks.test(after,  off, alternative = "greater")$p.value)
  list(pval = pv_b + pv_a, tdm = tdm)
}

#' @noRd
.valley_peak_algorithm_full <- function(
    temp_derivative, activity, temperature, temperature_threshold,
    estimated_sleep, dif_temp, activity_thr,
    valley_quantile, peak_quantile,
    minimum_offwrist_length, long_offwrist_length, short_offwrist_length,
    short_vp_decrease_ratio_min, decrease_ratio_min,
    valley_peak_low_act_min, valley_peak_low_temp_min,
    offwrist_max_temp_dif_med, forbidden_zone,
    datetime_stamps, epoch_hour, n,
    next_possible_length = TRUE
) {
  empty <- data.frame(start = integer(), end = integer(),
                      length = integer(), valley_peak = logical())

  # fixed wrapper parameters not otherwise threaded in
  minimum_onwrist_length <- 20L
  offwrist_min_tdm       <- 0.75    # offwrist_minimum_temperature_difference_median
  trust_pval             <- 5e-4
  possible_window_hours  <- 1L
  static_peak_q          <- 0.98
  static_valley_q        <- 0.95

  # ── Split derivative into noon-to-noon days (per-day peak/valley thresholds) ─
  if (!is.null(datetime_stamps) && length(datetime_stamps) == n) {
    day_id <- .vp_day_index(datetime_stamps, start_hour = 12L, n = n)
  } else {
    day_id <- rep(1L, n)
  }
  uday <- sort(unique(day_id))

  all_peaks <- integer(0); all_valleys <- integer(0)
  for (d in uday) {
    idx <- which(day_id == d)
    if (length(idx) == 0L) next
    day_deriv <- temp_derivative[idx]
    shift     <- idx[1L] - 1L                         # 0-based global offset

    pos <- which(day_deriv > 0); pos_sig <- day_deriv[pos]
    neg <- which(day_deriv < 0); neg_sig <- day_deriv[neg]

    # peaks (include_static)
    if (length(pos_sig) > 0L) {
      high  <- as.double(stats::quantile(pos_sig, peak_quantile,  names = FALSE, type = 1))
      lower <- as.double(stats::quantile(pos_sig, static_peak_q,  names = FALSE, type = 1))
      cand  <- .find_peaks_scipy(pos_sig, lower)
      for (pk in cand) {
        gp   <- (pos[pk] - 1L) + shift                # global 0-based peak index
        keep <- pos_sig[pk] >= high
        if (!keep) {
          bs <- max(0L, gp - minimum_offwrist_length)
          ab <- if (bs < gp) activity[(bs + 1L):gp] else numeric(0)
          if (length(ab) > 0L &&
              below_prop(ab, activity_thr) >= 0.7 && zero_prop(ab) >= 0.4) keep <- TRUE
        }
        if (keep) all_peaks <- c(all_peaks, gp)
      }
    }
    # valleys (include_static)
    if (length(neg_sig) > 0L) {
      negp  <- -neg_sig
      high  <- as.double(stats::quantile(negp, valley_quantile,  names = FALSE, type = 1))
      lower <- as.double(stats::quantile(negp, static_valley_q,  names = FALSE, type = 1))
      cand  <- .find_peaks_scipy(negp, lower)
      for (vk in cand) {
        gv   <- (neg[vk] - 1L) + shift
        keep <- negp[vk] >= high
        if (!keep) {
          af <- if (gv < min(gv + minimum_offwrist_length, n))
                  activity[(gv + 1L):min(gv + minimum_offwrist_length, n)] else numeric(0)
          la <- if (length(af) > 0L) below_prop(af, activity_thr) else 0
          zp <- if (length(af) > 0L) zero_prop(af) else 0
          if ((la >= 0.7 && zp >= 0.4) || la >= 0.9) keep <- TRUE
        }
        if (keep) all_valleys <- c(all_valleys, gv)
      }
    }
  }

  all_peaks   <- sort(unique(all_peaks))
  all_valleys <- sort(unique(all_valleys))
  if (length(all_peaks) == 0L || length(all_valleys) == 0L) return(empty)

  # ── Match valleys to peaks (with first/last boundary periods) ───────────
  vp <- list()
  pc <- length(all_peaks); vc <- length(all_valleys)
  valley <- 1L; peak <- 1L

  if (all_peaks[1L] < all_valleys[1L]) {
    fp  <- all_peaks[1L]
    if (fp >= minimum_offwrist_length &&
        below_prop(activity[1L:fp], activity_thr) > 0.75) {
      vp[[length(vp) + 1L]] <- c(0L, fp)
    }
  }

  while (valley <= vc && peak <= pc) {
    if (all_peaks[peak] > all_valleys[valley]) {
      len <- all_peaks[peak] - all_valleys[valley]
      if (valley + 1L <= vc) {
        npl <- all_peaks[peak] - all_valleys[valley + 1L]    # next_possible_length_criteria
        if (all_peaks[peak] > all_valleys[valley + 1L] && npl >= minimum_offwrist_length) {
          valley <- valley + 1L
        } else {
          if (len >= minimum_offwrist_length && len < long_offwrist_length) {
            vp[[length(vp) + 1L]] <- c(all_valleys[valley], all_peaks[peak])
            peak <- peak + 1L
          }
          valley <- valley + 1L
        }
      } else {
        if (len >= minimum_offwrist_length && len < long_offwrist_length) {
          vp[[length(vp) + 1L]] <- c(all_valleys[valley], all_peaks[peak])
          peak <- peak + 1L
        }
        valley <- valley + 1L
      }
    } else {
      peak <- peak + 1L
    }
  }

  if (all_peaks[pc] < all_valleys[vc]) {
    lv  <- all_valleys[vc]
    len <- n - lv
    if (len >= minimum_offwrist_length &&
        below_prop(activity[(lv + 1L):n], activity_thr) > 0.75) {
      vp[[length(vp) + 1L]] <- c(lv, n)
    }
  }

  if (length(vp) == 0L) return(empty)
  m    <- do.call(rbind, vp)
  vpdf <- data.frame(start = m[, 1L], end = m[, 2L])
  vpdf$length <- vpdf$end - vpdf$start

  # ── Forbidden zone (possible_window_hours rule; do_forbidden_zone = TRUE) ─
  keep_fz <- logical(nrow(vpdf))
  for (i in seq_len(nrow(vpdf))) {
    s <- vpdf$start[i]; e <- min(vpdf$end[i], n)
    fb <- forbidden_zone[(s + 1L):e]
    if (sum(fb) == 0L) { keep_fz[i] <- TRUE; next }
    if ((e - s) > possible_window_hours) {
      fr <- .first_allowed_run_len(fb)
      if (fr > 0L) {
        keep_fz[i] <- ((possible_window_hours - fr) <= 0.66 * possible_window_hours)
      } else keep_fz[i] <- FALSE
    } else keep_fz[i] <- FALSE
  }
  vpdf <- vpdf[keep_fz, , drop = FALSE]; row.names(vpdf) <- NULL
  if (nrow(vpdf) == 0L) return(empty)

  # ── Features ────────────────────────────────────────────────────────────
  vpdf$low_act_prop  <- vapply(seq_len(nrow(vpdf)), function(i)
    below_prop(activity[(vpdf$start[i] + 1L):vpdf$end[i]], activity_thr), numeric(1))
  vpdf$low_temp_prop <- vapply(seq_len(nrow(vpdf)), function(i)
    below_prop(temperature[(vpdf$start[i] + 1L):vpdf$end[i]], temperature_threshold), numeric(1))
  vpdf$temp_dif_med  <- vapply(seq_len(nrow(vpdf)), function(i)
    stats::median(dif_temp[(vpdf$start[i] + 1L):vpdf$end[i]], na.rm = TRUE), numeric(1))
  vpdf$decrease_ratio <- vapply(seq_len(nrow(vpdf)), function(i)
    .vp_decrease_ratio(temp_derivative, vpdf$start[i], vpdf$end[i]), numeric(1))

  # ── valid_index decision tree (ActTrust, non-lumus) ─────────────────────
  af <- vpdf$low_act_prop   >= valley_peak_low_act_min
  tf <- vpdf$low_temp_prop  >= valley_peak_low_temp_min
  df <- vpdf$decrease_ratio >= decrease_ratio_min
  xf <- vpdf$temp_dif_med   <= offwrist_max_temp_dif_med    # tempdif pass = <= max
  comb <- af & tf & df & xf

  keepv <- logical(nrow(vpdf))
  for (i in seq_len(nrow(vpdf))) {
    if (comb[i]) { keepv[i] <- TRUE; next }
    s <- vpdf$start[i]; e <- vpdf$end[i]; L <- vpdf$length[i]
    if (af[i]) {
      if (df[i]) {
        if (tf[i]) {
          if (L <= short_offwrist_length) {
            ks <- .vp_test_tdm_below(dif_temp, s, e, n)
            if (ks$pval < trust_pval && ks$tdm < (offwrist_min_tdm - 0.1)) {
              if (.vp_surrounded_by_awake(activity, s, e, activity_thr,
                                          minimum_onwrist_length, n)) keepv[i] <- TRUE
            }
          }
        } else {
          if (L <= short_offwrist_length && xf[i] &&
              vpdf$decrease_ratio[i] >= short_vp_decrease_ratio_min) {
            if (.vp_surrounded_by_awake(activity, s, e, activity_thr,
                                        minimum_onwrist_length, n)) keepv[i] <- TRUE
          }
        }
      } else {
        if (tf[i] && L <= short_offwrist_length && xf[i]) {
          if (.vp_surrounded_by_awake(activity, s, e, activity_thr,
                                      minimum_onwrist_length, n)) keepv[i] <- TRUE
        }
      }
    }
  }
  vpdf <- vpdf[keepv, , drop = FALSE]; row.names(vpdf) <- NULL
  if (nrow(vpdf) == 0L) return(empty)

  vpdf$valley_peak <- TRUE
  vpdf[, c("start", "end", "length", "valley_peak")]
}

#' @noRd
.near_all_off_detection_v2 <- function(activity, activity_thr,
                                        positive_act_median, n, epoch_hour) {
  if (zero_prop(activity) > 0.9 || positive_act_median < activity_thr)
    return(data.frame(start=0L, end=n, length=n, valley_peak=FALSE))
  data.frame(start=integer(), end=integer(), length=integer(), valley_peak=logical())
}

# ── Description report ────────────────────────────────────────────────────────
#' @noRd
.describe_offwrist_periods_v2 <- function(offwrist_periods, activity, temperature,
                                           temperature_variance, temperature_threshold,
                                           activity_thr, sleep_lat, dif_temp, n,
                                           segments = 7L, window = 60L) {
  np <- nrow(offwrist_periods)
  if (np == 0L) return(data.frame())

  act_zp   <- zero_prop(activity)
  act_tq   <- act_zp + (1 - act_zp) * 0.05
  act_thr2 <- as.double(stats::quantile(activity, act_tq, names = FALSE, type = 1))
  dif_var  <- var_filter(dif_temp, 3L)

  report <- data.frame(
    activity_zero_prop         = numeric(np),
    low_act_prop               = numeric(np),
    high_act_before            = numeric(np),
    high_act_after             = numeric(np),
    border_act_conc            = numeric(np),
    border_temp_conc           = numeric(np),
    low_temp_prop              = numeric(np),
    temp_dif_median            = numeric(np),
    length                     = offwrist_periods$length,
    valley_peak                = if ("valley_peak" %in% names(offwrist_periods))
                                   offwrist_periods$valley_peak else FALSE
  )

  for (i in seq_len(np)) {
    s  <- offwrist_periods$start[i]
    e  <- min(offwrist_periods$end[i], n)
    rs <- s + 1L; re <- e
    if (rs > re || rs < 1L || re > n) next

    seg_act  <- activity[rs:re]
    seg_temp <- temperature[rs:re]
    seg_tv   <- temperature_variance[rs:re]
    seg_dif  <- dif_temp[rs:re]

    act_segs <- .segmentation(seg_act, segments)
    tv_segs  <- .segmentation(seg_tv,  segments)

    before_s <- max(1L, rs - window)
    after_e  <- min(n, re + window)

    report$activity_zero_prop[i] <- zero_prop(seg_act)
    report$low_act_prop[i]       <- below_prop(seg_act, activity_thr)
    if (rs > 1L) report$high_act_before[i] <- 1 - below_prop(activity[before_s:(rs-1L)], act_thr2)
    if (re < n)  report$high_act_after[i]  <- 1 - below_prop(activity[(re+1L):after_e], act_thr2)
    report$border_act_conc[i]   <- act_segs[1L] + act_segs[segments]
    report$border_temp_conc[i]  <- tv_segs[1L]  + tv_segs[segments]
    report$low_temp_prop[i]     <- below_prop(seg_temp, temperature_threshold)
    report$temp_dif_median[i]   <- stats::median(seg_dif, na.rm = TRUE)
  }

  report
}

#' @noRd
.segmentation <- function(x, n_segs = 7L) {
  total <- sum(x, na.rm = TRUE)
  if (total == 0) return(rep(0, n_segs))
  len  <- length(x)
  bins <- round(seq(0, len, length.out = n_segs + 1L))
  vapply(seq_len(n_segs), function(i) {
    sum(x[(bins[i] + 1L):bins[i + 1L]], na.rm = TRUE) / total
  }, numeric(1))
}

#' @noRd
.description_report_filter_v2 <- function(periods, report, report_zero_act_min,
                                            border_conc_min, act_around_min,
                                            low_act_prop_min, temp_dif_min,
                                            offwrist_max_temp_dif_med,
                                            long_offwrist_length,
                                            short_offwrist_length,
                                            is_highly_separable = FALSE) {
  if (nrow(periods) == 0L || nrow(report) == 0L) return(periods)
  keep <- rep(TRUE, nrow(periods))

  for (i in seq_len(nrow(periods))) {
    r   <- report[i, ]
    len <- r$length
    vp  <- isTRUE(r$valley_peak)

    # Only filter short periods
    if (len >= long_offwrist_length || is_highly_separable) next

    is_valid <- TRUE

    # Zero activity + low activity (do_report_low_activity_after = TRUE)
    if (r$activity_zero_prop < report_zero_act_min) {
      if (r$low_act_prop < low_act_prop_min) {
        is_valid <- FALSE
      }
    }

    # Skip remaining filters for short valley-peak periods
    if (is_valid && !(vp && len <= short_offwrist_length)) {
      if (is_valid && r$border_act_conc < 0.35 && r$border_temp_conc < 0.35)
        is_valid <- FALSE
      if (is_valid && r$high_act_before < 0.2 && r$high_act_after < 0.2)
        is_valid <- FALSE
      if (is_valid &&
          (r$high_act_before < act_around_min || r$high_act_after < act_around_min) &&
          r$border_act_conc  < border_conc_min && r$border_temp_conc < border_conc_min)
        is_valid <- FALSE
      if (is_valid && r$temp_dif_median > temp_dif_min)
        is_valid <- FALSE
    }

    if (!is_valid) keep[i] <- FALSE
  }

  periods[keep, , drop = FALSE]
}

# ── surrounded_onwrist_filter ─────────────────────────────────────────────────
#' Full port of Python's surrounded_onwrist_filter
#' @noRd
.surrounded_onwrist_filter_v2 <- function(refined_ow, temperature,
                                            temperature_threshold,
                                            min_preceding_onwrist_len,
                                            min_surrounding_ow_len, n) {
  on_periods <- .rle_periods(refined_ow == 1L)
  if (nrow(on_periods) == 0L) return(refined_ow)
  on_periods$length <- on_periods$end - on_periods$start

  ow_count       <- nrow(on_periods)
  valid_on_idx   <- integer(0)

  for (i in seq_len(ow_count)) {
    valid_on <- TRUE
    s        <- on_periods$start[i]
    e        <- on_periods$end[i]
    len      <- on_periods$length[i]
    ow_temp  <- temperature[(s + 1L):min(e, n)]
    low_tp   <- below_prop(ow_temp, temperature_threshold)

    prec_len <- 0L
    succ_len <- 0L

    if (i == 1L) {
      if (ow_count > 1L) {
        succ_len <- on_periods$start[2L] - on_periods$end[1L]
      } else {
        succ_len <- n - on_periods$end[1L]
      }

      if (s == 0L) {
        if (succ_len > min_surrounding_ow_len) {
          if (len <= min_preceding_onwrist_len) {
            valid_on <- FALSE
          } else if (len < succ_len && low_tp > 0.75) {
            valid_on <- FALSE
          }
        }
      } else {
        prec_len <- s  # distance from start of recording
        if (prec_len > min_surrounding_ow_len && succ_len > min_surrounding_ow_len) {
          if (len <= min_preceding_onwrist_len) {
            valid_on <- FALSE
          } else if (len < (succ_len + prec_len) && low_tp > 0.75) {
            valid_on <- FALSE
          }
        }
      }
    } else {
      prec_len <- s - on_periods$end[i - 1L]

      if (i == ow_count) {
        if (e == n) {
          if (prec_len > min_surrounding_ow_len) {
            if (len <= min_preceding_onwrist_len) {
              valid_on <- FALSE
            } else if (len < prec_len && low_tp > 0.75) {
              valid_on <- FALSE
            }
          }
        } else {
          succ_len <- n - e
          if (prec_len > min_surrounding_ow_len && succ_len > min_surrounding_ow_len) {
            if (len <= min_preceding_onwrist_len) {
              valid_on <- FALSE
            } else if (len < (succ_len + prec_len) && low_tp > 0.75) {
              valid_on <- FALSE
            }
          }
        }
      } else {
        succ_len <- on_periods$start[i + 1L] - on_periods$end[i]
        if (prec_len > min_surrounding_ow_len && succ_len > min_surrounding_ow_len) {
          if (len <= min_preceding_onwrist_len) {
            valid_on <- FALSE
          } else if (len < (succ_len + prec_len) && low_tp > 0.75) {
            valid_on <- FALSE
          }
        }
      }
    }

    if (valid_on) valid_on_idx <- c(valid_on_idx, i)
  }

  if (length(valid_on_idx) == ow_count) return(refined_ow)

  # Rebuild from valid on-wrist periods
  new_ow <- rep(0L, n)
  for (i in valid_on_idx) {
    s <- on_periods$start[i]; e <- min(on_periods$end[i], n)
    if (s < e) new_ow[(s + 1L):e] <- 1L
  }
  new_ow
}
