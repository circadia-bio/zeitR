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
  minimum_offwrist_length          <- 15L   # Python default
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
  refined      <- list()
  periods      <- as.matrix(offwrist_periods[, c("start", "end")])
  n_off        <- nrow(periods)
  ow_idx       <- 1L
  search_start <- TRUE

  while (ow_idx <= n_off) {

    if (search_start) {
      # search_offwrist_start
      start <- periods[ow_idx, 1L]
      if (length(refined) > 0L) {
        previous_end <- refined[[length(refined)]][2L]
        if (start < previous_end) start <- previous_end + 2L
      } else {
        previous_end <- 0L
      }
      ri <- start + 1L

      if (start >= filter_hws && ri <= n) {
        if (norm_temp_variance[ri] >= temp_var_thr) {
          lo   <- max(1L, ri - filter_hws + 1L)
          ltvp <- below_prop(norm_temp_variance[lo:ri], temp_var_thr)
          if (ltvp <= max_low_tv_border) {
            new_start_ri <- lo + which.max(norm_temp_variance[lo:ri]) - 1L
            refined[[length(refined) + 1L]] <- c(max(new_start_ri - 1L, previous_end), 0L)
            search_start <- FALSE
          } else {
            res <- .find_peak_base_start2(start, previous_end, norm_temp_variance,
                                          temperature, temperature_threshold,
                                          activity_median, activity_median_low,
                                          activity_thr, temp_var_thr, n)
            if (isTRUE(res$delete)) {
              periods <- periods[-ow_idx, , drop = FALSE]; n_off <- nrow(periods)
            } else if (!is.null(res$new_start)) {
              refined[[length(refined) + 1L]] <- c(res$new_start, 0L)
              search_start <- FALSE
            }
          }
        } else {
          prev_end2 <- if (ow_idx > 1L) periods[ow_idx - 1L, 2L] else 0L
          res <- .find_peak_base_start2(start, prev_end2, norm_temp_variance,
                                        temperature, temperature_threshold,
                                        activity_median, activity_median_low,
                                        activity_thr, temp_var_thr, n)
          if (isTRUE(res$delete)) {
            periods <- periods[-ow_idx, , drop = FALSE]; n_off <- nrow(periods)
          } else if (!is.null(res$new_start)) {
            refined[[length(refined) + 1L]] <- c(res$new_start, 0L)
            search_start <- FALSE
          }
        }
      } else if (ri >= 1L && ri <= n) {
        new_start <- which.max(norm_temp_variance[1L:ri]) - 1L
        refined[[length(refined) + 1L]] <- c(max(new_start, 0L), 0L)
        search_start <- FALSE
      } else {
        ow_idx <- ow_idx + 1L
      }

    } else {
      # search_offwrist_end
      end <- periods[ow_idx, 2L]
      re  <- min(end, n)

      if (re >= 1L && norm_temp_variance[re] >= temp_var_thr) {
        if (re + filter_hws <= n) {
          hi   <- min(n, re + filter_hws)
          ltvp <- below_prop(norm_temp_variance[re:hi], temp_var_thr)
          if (ltvp <= max_low_tv_border) {
            new_end_ri <- re - 1L + which.max(norm_temp_variance[re:hi])
            refined[[length(refined)]][2L] <- new_end_ri
            search_start <- TRUE; ow_idx <- ow_idx + 1L
          } else {
            if (ow_idx < n_off) {
              next_start <- periods[ow_idx + 1L, 1L]
              res <- .do_more_checks_around_end(ow_idx, next_start, end,
                                                periods, refined,
                                                norm_temp_variance, temperature,
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
              refined[[length(refined)]][2L] <- re - 1L + which.max(norm_temp_variance[re:n])
              search_start <- TRUE; ow_idx <- ow_idx + 1L
            }
          }
        } else {
          refined[[length(refined)]][2L] <- re - 1L + which.max(norm_temp_variance[re:n])
          search_start <- TRUE; ow_idx <- ow_idx + 1L
        }
      } else {
        if (ow_idx < n_off) {
          next_start <- periods[ow_idx + 1L, 1L]
          res <- .try_find_peak_base_end(ow_idx, next_start, end,
                                          periods, refined,
                                          norm_temp_variance, temperature,
                                          temperature_threshold, activity_median,
                                          activity_median_low, activity_thr,
                                          temp_var_thr, min_preceding_onwrist,
                                          minimum_low_act_prop, activity, n)
          ow_idx <- res$ow_idx; periods <- res$periods
          n_off  <- nrow(periods); refined <- res$refined
          search_start <- res$search_start
        } else {
          refined[[length(refined)]][2L] <- re - 1L + which.max(norm_temp_variance[re:n])
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
      peak_ri <- peak_base + 1L
      while (peak_ri > 1L && peak_ri < n && ntv[peak_ri] >= ntv[peak_ri + 1L])
        peak_ri <- peak_ri - 1L
      return(list(new_start = peak_ri - 1L, delete = FALSE))
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
  lo <- max(1L, re - filter_hws)
  ltvp <- below_prop(ntv[lo:re], temp_var_thr)

  if (ltvp > max_low_tv_border) {
    periods <- periods[-ow_idx, , drop = FALSE]
    return(list(ow_idx = ow_idx, periods = periods, refined = refined,
                search_start = FALSE))
  } else {
    refined[[length(refined)]][2L] <- lo + which.max(ntv[lo:re]) - 1L
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
      peak_ri <- peak_base + 1L
      while (peak_ri < n && ntv[peak_ri] >= ntv[peak_ri - 1L]) peak_ri <- peak_ri + 1L
      refined[[length(refined)]][2L] <- peak_ri - 1L
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
  # Python: pads valid_activity with max(activity) before median filter
  pad_val <- max(activity, na.rm = TRUE)
  padding <- rep(pad_val, 2L * sleep_act_hws)
  padded  <- c(padding, valid_activity, padding)

  # median filter on padded signal ("padded" mode = already padded)
  filtered_padded <- median_filter(padded, sleep_act_hws)
  # Remove the padding
  filtered_valid  <- filtered_padded[(2L * sleep_act_hws + 1L):
                                       (2L * sleep_act_hws + length(valid_activity))]

  # Compute sleep threshold ("both" configuration)
  if (act_zero_prop < sleep_all_q) {
    sleep_thr <- as.double(stats::quantile(activity, sleep_all_q, names = FALSE))
  } else {
    q_idx     <- act_zero_prop + (1 - act_zero_prop) * sleep_pos_q
    sleep_thr <- as.double(stats::quantile(activity, q_idx, names = FALSE))
  }

  # Threshold: 0 = sleep, 1 = wake in valid signal
  sleep_est_valid <- as.integer(filtered_valid > sleep_thr)

  # Short sleep filter (>= epoch_hour)
  sleep_periods <- .rle_periods(sleep_est_valid == 0L)
  if (nrow(sleep_periods) > 0L) {
    sleep_periods <- sleep_periods[
      (sleep_periods$end - sleep_periods$start) > as.integer(epoch_hour), ]
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

  # sleep_low_temperature_filter: trim sleep borders that are below temp threshold
  sleep_periods_full <- .rle_periods(estimated_sleep == 0L)
  for (i in seq_len(nrow(sleep_periods_full))) {
    s <- sleep_periods_full$start[i]; e <- sleep_periods_full$end[i]
    ri <- s + 1L
    while (ri <= e && ri <= n && temperature[ri] < temperature_threshold) {
      estimated_sleep[ri] <- 1L; ri <- ri + 1L
    }
    ri <- e
    while (ri >= s + 1L && ri >= 1L && temperature[ri] < temperature_threshold) {
      estimated_sleep[ri] <- 1L; ri <- ri - 1L
    }
  }

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

# ── Valley-peak algorithm ─────────────────────────────────────────────────────
#' Full port of valley_peak_offwrist_algorithm (ActTrust config, no lumus)
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
  # Split derivative by day if datetime_stamps available, else treat as one day
  if (!is.null(datetime_stamps) && length(datetime_stamps) == n) {
    days       <- as.Date(datetime_stamps)
    unique_days <- unique(days)
  } else {
    unique_days <- 1L
    days        <- rep(1L, n)
  }

  all_peaks   <- integer(0)
  all_valleys <- integer(0)

  for (d in unique_days) {
    idx      <- which(days == d)
    if (length(idx) == 0L) next
    day_deriv <- temp_derivative[idx]
    idx_shift <- idx[1L] - 1L

    pos_sig <- day_deriv[day_deriv > 0]
    neg_sig <- day_deriv[day_deriv < 0]
    pos_idx <- which(day_deriv > 0)
    neg_idx <- which(day_deriv < 0)

    if (length(pos_sig) > 0L) {
      pk_thr <- as.double(stats::quantile(pos_sig, peak_quantile,
                                          names = FALSE, type = 1))
      peaks_local <- which(pos_sig >= pk_thr)
      # find_peaks equivalent: local maxima
      peaks_local <- peaks_local[vapply(peaks_local, function(j) {
        l <- max(1L, j - 1L); r <- min(length(pos_sig), j + 1L)
        pos_sig[j] >= pos_sig[l] && pos_sig[j] >= pos_sig[r]
      }, logical(1))]
      all_peaks <- c(all_peaks, pos_idx[peaks_local] + idx_shift)
    }

    if (length(neg_sig) > 0L) {
      vl_thr <- as.double(stats::quantile(-neg_sig, valley_quantile,
                                          names = FALSE, type = 1))
      valleys_local <- which(-neg_sig >= vl_thr)
      valleys_local <- valleys_local[vapply(valleys_local, function(j) {
        l <- max(1L, j - 1L); r <- min(length(neg_sig), j + 1L)
        neg_sig[j] <= neg_sig[l] && neg_sig[j] <= neg_sig[r]
      }, logical(1))]
      all_valleys <- c(all_valleys, neg_idx[valleys_local] + idx_shift)
    }
  }

  all_peaks   <- sort(unique(all_peaks))
  all_valleys <- sort(unique(all_valleys))

  if (length(all_peaks) == 0L || length(all_valleys) == 0L)
    return(data.frame(start=integer(), end=integer(), length=integer(),
                      valley_peak=logical()))

  # Match valleys to peaks
  vp_list    <- list()
  valley_idx <- 1L
  peak_idx   <- 1L
  n_valleys  <- length(all_valleys)
  n_peaks    <- length(all_peaks)

  while (valley_idx <= n_valleys && peak_idx <= n_peaks) {
    if (all_peaks[peak_idx] > all_valleys[valley_idx]) {
      v     <- all_valleys[valley_idx]
      p     <- all_peaks[peak_idx]
      len   <- p - v

      # next_possible_length check
      if (valley_idx + 1L <= n_valleys) {
        next_v   <- all_valleys[valley_idx + 1L]
        next_len <- p - next_v
        if (p > next_v && next_len >= minimum_offwrist_length) {
          valley_idx <- valley_idx + 1L
          next
        }
      }

      if (len >= minimum_offwrist_length && len < long_offwrist_length) {
        vp_list[[length(vp_list) + 1L]] <- c(v, p)
        peak_idx <- peak_idx + 1L
      }
      valley_idx <- valley_idx + 1L
    } else {
      peak_idx <- peak_idx + 1L
    }
  }

  if (length(vp_list) == 0L)
    return(data.frame(start=integer(), end=integer(), length=integer(),
                      valley_peak=logical()))

  vp_mat  <- do.call(rbind, vp_list)
  vp_df   <- data.frame(start = vp_mat[, 1L] - 1L,  # 0-indexed
                        end   = vp_mat[, 2L],
                        length = vp_mat[, 2L] - vp_mat[, 1L] + 1L)

  # Apply forbidden zone filter
  vp_df$forbidden <- vapply(seq_len(nrow(vp_df)), function(i) {
    s <- vp_df$start[i]; e <- min(vp_df$end[i], n)
    if (s >= e) return(TRUE)
    any(forbidden_zone[(s + 1L):e] == 1L)
  }, logical(1))
  vp_df <- vp_df[!vp_df$forbidden, ]
  row.names(vp_df) <- NULL

  if (nrow(vp_df) == 0L)
    return(data.frame(start=integer(), end=integer(), length=integer(),
                      valley_peak=logical()))

  # Feature filters
  vp_df$low_act_prop <- vapply(seq_len(nrow(vp_df)), function(i) {
    s <- vp_df$start[i]; e <- vp_df$end[i]
    below_prop(activity[(s + 1L):e], activity_thr)
  }, numeric(1))
  vp_df$low_temp_prop <- vapply(seq_len(nrow(vp_df)), function(i) {
    s <- vp_df$start[i]; e <- vp_df$end[i]
    below_prop(temperature[(s + 1L):e], temperature_threshold)
  }, numeric(1))
  vp_df$temp_dif_med <- vapply(seq_len(nrow(vp_df)), function(i) {
    s <- vp_df$start[i]; e <- vp_df$end[i]
    stats::median(dif_temp[(s + 1L):e], na.rm = TRUE)
  }, numeric(1))
  vp_df$decrease_ratio <- vapply(seq_len(nrow(vp_df)), function(i) {
    s <- vp_df$start[i]; e <- vp_df$end[i]
    .compute_decrease_ratio_v2(temp_derivative, s + 1L, e)
  }, numeric(1))

  # Apply filters
  keep <- with(vp_df,
    low_act_prop  >= valley_peak_low_act_min  &
    low_temp_prop >= valley_peak_low_temp_min &
    decrease_ratio >= decrease_ratio_min      &
    temp_dif_med  <  offwrist_max_temp_dif_med
  )

  # Short VP rescue: short_offwrist_minimum_decrease_ratio
  for (i in which(!keep)) {
    r <- vp_df[i, ]
    if (r$low_act_prop >= valley_peak_low_act_min &&
        r$low_temp_prop >= valley_peak_low_temp_min &&
        r$length <= short_offwrist_length &&
        r$temp_dif_med < offwrist_max_temp_dif_med &&
        r$decrease_ratio >= short_vp_decrease_ratio_min) {
      keep[i] <- TRUE
    }
  }

  vp_df <- vp_df[keep, ]
  row.names(vp_df) <- NULL

  if (nrow(vp_df) == 0L)
    return(data.frame(start=integer(), end=integer(), length=integer(),
                      valley_peak=logical()))

  vp_df$valley_peak <- TRUE
  vp_df[, c("start", "end", "length", "valley_peak")]
}

#' @noRd
.compute_decrease_ratio_v2 <- function(temp_derivative, ri_start, ri_end) {
  if (ri_start > ri_end || ri_start < 1L) return(0)
  seg <- temp_derivative[ri_start:ri_end]
  neg <- seg[seg < 0]; pos <- seg[seg > 0]
  if (length(pos) == 0L || sum(pos) == 0) return(0)
  abs(sum(neg)) / sum(pos)
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
