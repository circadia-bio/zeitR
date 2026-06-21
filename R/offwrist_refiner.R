# ── Bimodal off-wrist refiner ─────────────────────────────────────────────────
# R port of BimodalOffwristRefiner from
# condor_pipeline/algorithms/vendor/condor/bimodal_offwrist_refine_without_prints.py
# Author of original Python: Julius A. P. P. de Paula (Condor Instruments, 2023)
#
# INDEXING CONVENTION: All period data frames use Python-style indexing:
#   start = 0-based first index
#   end   = exclusive last index  (i.e. length = end - start)
# R vector indexing therefore uses (start+1):end throughout.
# ─────────────────────────────────────────────────────────────────────────────

#' Refine an initial off-wrist detection using the Condor three-stage algorithm
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

  minimum_onwrist_length           <- 20L
  minimum_offwrist_length          <- 10L
  minimum_preceding_onwrist_length <- 40L
  activity_threshold_quantile      <- 0.24
  minimum_low_activity_proportion  <- 0.8   # minimum_low_activity_proportion (NOT bimodal variant)
  max_low_temp_var_prop_border     <- 0.3
  long_offwrist_length             <- 4L * as.integer(epoch_hour)
  short_offwrist_length            <- 20L
  temp_var_thr_quantile            <- 0.9
  sleep_act_hws                    <- 120L
  sleep_all_q                      <- 0.4
  sleep_pos_q                      <- 0.05
  max_offwrist_sleep_prop          <- 0.4
  ashman_d_min                     <- 1.5
  ashman_d_max                     <- 2.6
  bimodal_max_offwrist_prop        <- 0.0
  bimodal_min_low_act_prop         <- 1.0
  report_zero_act_prop_min         <- 0.35
  border_conc_min                  <- 0.5
  report_act_around_min            <- 0.1
  report_low_act_prop_min          <- 0.5
  offwrist_max_temp_dif_median     <- 0.8
  offwrist_min_temp_dif_median     <- 0.65
  valley_quantile                  <- 0.99
  peak_quantile                    <- 0.99
  valley_peak_low_temp_min         <- 0.5
  short_vp_decrease_ratio_min      <- 1.25
  sleep_low_temp_prop_max          <- 0.5
  half_filter_hws                  <- as.integer(filter_hws / 2)

  offwrist_periods <- .rle_periods(initial_offwrist == 0L)
  if (nrow(offwrist_periods) == 0L) return(rep(1L, n))

  act_zero_prop <- zero_prop(activity)
  act_thr_q     <- act_zero_prop + (1 - act_zero_prop) * activity_threshold_quantile
  activity_thr  <- as.double(stats::quantile(activity, act_thr_q, names = FALSE, type = 1))
  if (activity_thr < 200) activity_thr <- 200

  pos_act_q           <- act_zero_prop + (1 - act_zero_prop) * 0.5
  positive_act_median <- as.double(stats::quantile(activity, pos_act_q, names = FALSE, type = 1))
  temp_var_thr        <- as.double(stats::quantile(norm_temp_variance, temp_var_thr_quantile, names = FALSE))

  # ── Stage 1: filter short on-wrist periods ───────────────────────────────
  onwrist_periods <- .rle_periods(initial_offwrist == 1L)
  onwrist_periods <- onwrist_periods[
    (onwrist_periods$end - onwrist_periods$start) > minimum_onwrist_length, ]
  onwrist_periods <- onwrist_periods[order(onwrist_periods$start), ]
  row.names(onwrist_periods) <- NULL

  stage1_onwrist <- integer(n)
  for (i in seq_len(nrow(onwrist_periods))) {
    s <- onwrist_periods$start[i]; e <- onwrist_periods$end[i]
    if (s < e) stage1_onwrist[(s + 1L):e] <- 1L
  }
  offwrist_periods <- .rle_periods(stage1_onwrist == 0L)
  if (nrow(offwrist_periods) == 0L) return(rep(1L, n))

  # ── Stage 2: border refinement ────────────────────────────────────────────
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
  if (nrow(refined_offwrist_periods) == 0L) return(rep(1L, n))

  # ── Stage 3: per-period features ─────────────────────────────────────────
  refined_offwrist_periods$length <-
    refined_offwrist_periods$end - refined_offwrist_periods$start

  refined_offwrist_periods$low_temp_prop <- vapply(
    seq_len(nrow(refined_offwrist_periods)), function(i) {
      s <- refined_offwrist_periods$start[i]; e <- refined_offwrist_periods$end[i]
      below_prop(temperature[(s + 1L):e], temperature_threshold)
    }, numeric(1))

  refined_offwrist_periods$low_act_prop <- vapply(
    seq_len(nrow(refined_offwrist_periods)), function(i) {
      s <- refined_offwrist_periods$start[i]; e <- refined_offwrist_periods$end[i]
      below_prop(activity[(s + 1L):e], activity_thr)
    }, numeric(1))

  refined_offwrist_periods$zero_act_prop <- vapply(
    seq_len(nrow(refined_offwrist_periods)), function(i) {
      s <- refined_offwrist_periods$start[i]; e <- refined_offwrist_periods$end[i]
      zero_prop(activity[(s + 1L):e])
    }, numeric(1))

  refined_offwrist_periods <- refined_offwrist_periods[
    refined_offwrist_periods$low_act_prop > minimum_low_activity_proportion, ]
  row.names(refined_offwrist_periods) <- NULL

  refined_offwrist_periods <- refined_offwrist_periods[
    refined_offwrist_periods$length >= minimum_offwrist_length, ]
  row.names(refined_offwrist_periods) <- NULL
  if (nrow(refined_offwrist_periods) == 0L) return(rep(1L, n))

  refined_offwrist_periods$temp_dif_median <- vapply(
    seq_len(nrow(refined_offwrist_periods)), function(i) {
      s <- refined_offwrist_periods$start[i]; e <- refined_offwrist_periods$end[i]
      stats::median(dif_temp[(s + 1L):e], na.rm = TRUE)
    }, numeric(1))

  refined_offwrist_periods <- refined_offwrist_periods[
    refined_offwrist_periods$temp_dif_median < offwrist_max_temp_dif_median, ]
  row.names(refined_offwrist_periods) <- NULL
  if (nrow(refined_offwrist_periods) == 0L) return(rep(1L, n))

  # ── Mark long off-wrist periods (>= 4h) ──────────────────────────────────
  is_long         <- refined_offwrist_periods$length >= long_offwrist_length
  short_only_mask <- rep(TRUE, n)
  for (i in which(is_long)) {
    s <- refined_offwrist_periods$start[i]; e <- refined_offwrist_periods$end[i]
    if (s < e) short_only_mask[(s + 1L):e] <- FALSE
  }
  valid_activity <- activity[short_only_mask]

  # ── Sleep estimation ──────────────────────────────────────────────────────
  estimated_sleep <- .estimate_sleep(
    activity = activity, valid_activity = valid_activity, n = n,
    short_only_mask = short_only_mask, sleep_act_hws = sleep_act_hws,
    sleep_all_q = sleep_all_q, sleep_pos_q = sleep_pos_q,
    temperature = temperature, temperature_threshold = temperature_threshold,
    sleep_low_temp_prop_max = sleep_low_temp_prop_max, epoch_hour = epoch_hour
  )

  # ── Sleep filter ──────────────────────────────────────────────────────────
  refined_offwrist_periods$sleep_prop <- vapply(
    seq_len(nrow(refined_offwrist_periods)), function(i) {
      s <- refined_offwrist_periods$start[i]; e <- refined_offwrist_periods$end[i]
      mean(estimated_sleep[(s + 1L):e] == 1L)
    }, numeric(1))

  refined_offwrist_periods <- refined_offwrist_periods[
    refined_offwrist_periods$sleep_prop <= max_offwrist_sleep_prop |
      is_long[seq_len(nrow(refined_offwrist_periods))], ]
  row.names(refined_offwrist_periods) <- NULL
  if (nrow(refined_offwrist_periods) == 0L) return(rep(1L, n))

  # ── Forbidden zones ───────────────────────────────────────────────────────
  forbidden_zone <- .compute_forbidden_zone(estimated_sleep, epoch_hour, n)

  refined_offwrist_periods$in_forbidden <- vapply(
    seq_len(nrow(refined_offwrist_periods)), function(i) {
      s <- refined_offwrist_periods$start[i]; e <- refined_offwrist_periods$end[i]
      any(forbidden_zone[(s + 1L):e] == 1L)
    }, logical(1))

  refined_offwrist_periods <- refined_offwrist_periods[
    !refined_offwrist_periods$in_forbidden, ]
  row.names(refined_offwrist_periods) <- NULL

  # ── Bimodality check ──────────────────────────────────────────────────────
  is_bimodal <- .check_bimodality(
    initial_offwrist = initial_offwrist, activity = activity,
    activity_thr = activity_thr, ashman = ashman,
    ashman_d_min = ashman_d_min, ashman_d_max = ashman_d_max,
    bimodal_max_offwrist_p = bimodal_max_offwrist_prop,
    bimodal_min_low_act_p = bimodal_min_low_act_prop, n = n
  )

  # ── Valley-peak algorithm ─────────────────────────────────────────────────
  if (is_bimodal) {
    vp_periods <- .valley_peak_algorithm(
      temp_derivative = temp_derivative, temperature = temperature,
      temperature_threshold = temperature_threshold,
      estimated_sleep = estimated_sleep,
      valley_quantile = valley_quantile, peak_quantile = peak_quantile,
      minimum_offwrist_length = minimum_offwrist_length,
      short_offwrist_length = short_offwrist_length,
      short_vp_decrease_ratio_min = short_vp_decrease_ratio_min,
      valley_peak_low_temp_min = valley_peak_low_temp_min,
      epoch_hour = epoch_hour, n = n,
      next_possible_length = TRUE, short_criteria = TRUE, medium_criteria = TRUE
    )

    long_df  <- refined_offwrist_periods[
      refined_offwrist_periods$length >= long_offwrist_length, ]
    short_df <- refined_offwrist_periods[
      refined_offwrist_periods$length <  long_offwrist_length, ]

    if (nrow(vp_periods) > 0) {
      vp_periods$length      <- vp_periods$end - vp_periods$start
      vp_periods$valley_peak <- TRUE
      short_df$valley_peak   <- FALSE
      long_df$valley_peak    <- FALSE
      combined <- rbind(
        long_df[,  c("start", "end", "length", "valley_peak")],
        short_df[, c("start", "end", "length", "valley_peak")],
        vp_periods[, c("start", "end", "length", "valley_peak")]
      )
    } else {
      short_df$valley_peak <- FALSE
      long_df$valley_peak  <- FALSE
      combined <- rbind(
        long_df[,  c("start", "end", "length", "valley_peak")],
        short_df[, c("start", "end", "length", "valley_peak")]
      )
    }

    combined <- combined[order(combined$start), ]
    row.names(combined) <- NULL
    refined_offwrist_periods <- combined

    refined_offwrist_periods$low_act_prop <- vapply(
      seq_len(nrow(refined_offwrist_periods)), function(i) {
        s <- refined_offwrist_periods$start[i]; e <- refined_offwrist_periods$end[i]
        below_prop(activity[(s + 1L):e], activity_thr)
      }, numeric(1))
    refined_offwrist_periods$zero_act_prop <- vapply(
      seq_len(nrow(refined_offwrist_periods)), function(i) {
        s <- refined_offwrist_periods$start[i]; e <- refined_offwrist_periods$end[i]
        zero_prop(activity[(s + 1L):e])
      }, numeric(1))
    refined_offwrist_periods$low_temp_prop <- vapply(
      seq_len(nrow(refined_offwrist_periods)), function(i) {
        s <- refined_offwrist_periods$start[i]; e <- refined_offwrist_periods$end[i]
        below_prop(temperature[(s + 1L):e], temperature_threshold)
      }, numeric(1))
    refined_offwrist_periods$length <-
      refined_offwrist_periods$end - refined_offwrist_periods$start

  } else {
    if (do_near_all_off_detection) {
      refined_offwrist_periods <- .near_all_off_detection(
        activity = activity, activity_thr = activity_thr,
        positive_act_median = positive_act_median, n = n, epoch_hour = epoch_hour
      )
    } else {
      refined_offwrist_periods <- data.frame(
        start = integer(), end = integer(), length = integer(), valley_peak = logical()
      )
    }
    if (nrow(refined_offwrist_periods) == 0L) return(rep(1L, n))
  }

  # ── Description-report filter ─────────────────────────────────────────────
  if (nrow(refined_offwrist_periods) > 0L) {
    report <- .describe_offwrist_periods(
      offwrist_periods = refined_offwrist_periods, activity = activity,
      temperature = temperature, temperature_variance = norm_temp_variance,
      temperature_threshold = temperature_threshold, activity_thr = activity_thr,
      dif_temp = dif_temp, n = n
    )
    refined_offwrist_periods <- .description_report_filter(
      periods = refined_offwrist_periods, report = report,
      report_zero_act_min = report_zero_act_prop_min,
      border_conc_min = border_conc_min, act_around_min = report_act_around_min,
      low_act_prop_min = report_low_act_prop_min,
      temp_dif_min = offwrist_min_temp_dif_median,
      temp_dif_max = offwrist_max_temp_dif_median,
      long_offwrist_length = long_offwrist_length,
      short_offwrist_length = short_offwrist_length,
      is_highly_separable = (ashman > 3)
    )
    row.names(refined_offwrist_periods) <- NULL
  }

  # ── Surrounded on-wrist filter ────────────────────────────────────────────
  refined_offwrist_periods <- .surrounded_onwrist_filter(
    periods = refined_offwrist_periods, n = n
  )

  # ── Final border snap ─────────────────────────────────────────────────────
  if (nrow(refined_offwrist_periods) > 0L) {
    refined_offwrist_periods$length <-
      refined_offwrist_periods$end - refined_offwrist_periods$start
    if (refined_offwrist_periods$start[1] <= minimum_offwrist_length)
      refined_offwrist_periods$start[1] <- 0L
    last <- nrow(refined_offwrist_periods)
    if ((n - refined_offwrist_periods$end[last]) <= minimum_offwrist_length)
      refined_offwrist_periods$end[last] <- n
  }

  # ── Assemble final output ─────────────────────────────────────────────────
  refined <- rep(1L, n)
  for (i in seq_len(nrow(refined_offwrist_periods))) {
    s <- refined_offwrist_periods$start[i]
    e <- min(refined_offwrist_periods$end[i], n)
    if (s < e) refined[(s + 1L):e] <- 0L
  }
  refined
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' Convert a binary vector to runs of 1s using Python-style 0-indexed exclusive ends
#' length = end - start; R vector access: x[(start+1):end]
#' @noRd
.rle_periods <- function(x) {
  r    <- rle(as.integer(x))
  ends <- cumsum(r$lengths)      # 1-indexed inclusive (R native) = Python exclusive end
  strt <- ends - r$lengths       # 0-indexed start
  keep <- r$values == 1L
  if (!any(keep)) return(data.frame(start = integer(), end = integer()))
  data.frame(start = strt[keep], end = ends[keep])
}

#' Proportion of values in x below threshold
#' @noRd
below_prop <- function(x, thr) {
  if (length(x) == 0L) return(0)
  sum(x < thr, na.rm = TRUE) / length(x)
}

#' Refine borders of each off-wrist period using temperature variance peaks
#' All indices are 0-based exclusive (Python convention)
#' @noRd
.refine_borders <- function(
    offwrist_periods, norm_temp_variance, activity_median, activity_median_low,
    temperature, temperature_threshold, activity_thr, minimum_low_act_prop,
    min_preceding_onwrist, filter_hws, half_filter_hws, max_low_tv_border,
    temp_var_thr, activity, n
) {
  refined  <- list()
  periods  <- as.matrix(offwrist_periods[, c("start", "end")])
  n_off    <- nrow(periods)
  idx      <- 1L
  prev_end <- 0L

  while (idx <= n_off) {
    start <- periods[idx, 1L]   # 0-indexed
    end   <- periods[idx, 2L]   # exclusive

    if (start < prev_end) start <- prev_end + 2L
    start_ok <- FALSE

    ri <- start + 1L  # R 1-indexed equivalent of start
    if (ri >= filter_hws && ri <= n) {
      if (norm_temp_variance[ri] >= temp_var_thr) {
        lo   <- max(1L, ri - filter_hws)
        ltvp <- below_prop(norm_temp_variance[lo:ri], temp_var_thr)
        if (ltvp <= max_low_tv_border) {
          new_start <- lo + which.max(norm_temp_variance[lo:ri]) - 2L
          new_start <- max(new_start, prev_end)
          refined[[length(refined) + 1L]] <- c(new_start, 0L)
          start_ok <- TRUE
        } else {
          res <- .find_peak_base_start(start, prev_end, norm_temp_variance,
                                       activity_median, activity_median_low,
                                       temperature, temperature_threshold,
                                       activity_thr, temp_var_thr, filter_hws, n)
          if (!is.null(res)) { refined[[length(refined) + 1L]] <- c(res, 0L); start_ok <- TRUE }
          else { idx <- idx + 1L; next }
        }
      } else {
        res <- .find_peak_base_start(start, prev_end, norm_temp_variance,
                                     activity_median, activity_median_low,
                                     temperature, temperature_threshold,
                                     activity_thr, temp_var_thr, filter_hws, n)
        if (!is.null(res)) { refined[[length(refined) + 1L]] <- c(res, 0L); start_ok <- TRUE }
        else { idx <- idx + 1L; next }
      }
    } else if (ri >= 1L) {
      new_start <- which.max(norm_temp_variance[1:ri]) - 2L
      new_start <- max(new_start, 0L)
      refined[[length(refined) + 1L]] <- c(new_start, 0L)
      start_ok <- TRUE
    }

    if (!start_ok) { idx <- idx + 1L; next }

    end_done <- FALSE
    re <- min(end, n)  # R index (end is exclusive so re IS the last valid R index)

    if (re >= 1L && norm_temp_variance[re] >= temp_var_thr) {
      hi   <- min(n, re + filter_hws)
      ltvp <- below_prop(norm_temp_variance[re:hi], temp_var_thr)
      if (ltvp <= max_low_tv_border) {
        new_end <- re - 1L + which.max(norm_temp_variance[re:hi])
        refined[[length(refined)]][2L] <- new_end
        prev_end <- new_end
        idx <- idx + 1L
        end_done <- TRUE
      }
    }

    if (!end_done) {
      if (idx < n_off) {
        next_start <- periods[idx + 1L, 1L]
        gap        <- next_start - end
        too_close  <- FALSE
        if (gap <= min_preceding_onwrist) {
          # activity between end and next_start: R indices (end+1):next_start
          if ((end + 1L) <= next_start && next_start <= n) {
            act_seg <- activity[(end + 1L):next_start]
            if (below_prop(act_seg, activity_thr) > minimum_low_act_prop) too_close <- TRUE
          }
        }
        if (too_close) {
          periods <- periods[-idx, , drop = FALSE]; n_off <- n_off - 1L
        } else {
          res <- .find_peak_base_end(end, next_start, norm_temp_variance,
                                     activity_median, activity_median_low,
                                     temperature, temperature_threshold,
                                     activity_thr, temp_var_thr, filter_hws, n)
          if (!is.null(res)) { refined[[length(refined)]][2L] <- res; prev_end <- res }
          else refined[[length(refined)]] <- NULL
          idx <- idx + 1L
        }
      } else {
        hi      <- min(n, re + filter_hws)
        new_end <- re - 1L + which.max(norm_temp_variance[re:hi])
        refined[[length(refined)]][2L] <- new_end
        prev_end <- new_end
        idx <- idx + 1L
      }
    }
  }

  if (length(refined) == 0L) return(data.frame(start = integer(), end = integer()))
  m <- do.call(rbind, refined)
  m <- m[m[, 2L] > 0L, , drop = FALSE]
  if (nrow(m) == 0L) return(data.frame(start = integer(), end = integer()))
  data.frame(start = m[, 1L], end = m[, 2L])
}

#' Search backwards for a temperature-variance peak base
#' start is 0-indexed; returns 0-indexed result
#' @noRd
.find_peak_base_start <- function(start, prev_end, ntv, act_median,
                                   act_med_low, temperature,
                                   temperature_threshold, activity_thr,
                                   temp_var_thr, filter_hws, n) {
  i          <- start
  last_valid <- start
  peak_base  <- start

  while (i > prev_end) {
    ri <- i + 1L
    if (ri < 1L || ri > n) break
    valid <- .check_valid_border_mod(act_med_low, temperature, temperature_threshold,
                                     act_median, activity_thr, ri)
    if (valid) {
      if (ntv[ri] >= temp_var_thr) { peak_base <- i; break }
      i <- i - 1L
    } else { last_valid <- i; break }
  }

  if (last_valid == start) {
    if (peak_base < start) {
      peak_ri <- peak_base + 1L
      while (peak_ri > 2L && peak_ri <= n && ntv[peak_ri] >= ntv[peak_ri + 1L])
        peak_ri <- peak_ri - 1L
      return(peak_ri - 1L)
    } else return(NULL)
  } else return(last_valid)
}

#' Search forwards for a temperature-variance peak base
#' end and next_start are 0-indexed exclusive; returns 0-indexed result
#' @noRd
.find_peak_base_end <- function(end, next_start, ntv, act_median,
                                 act_med_low, temperature,
                                 temperature_threshold, activity_thr,
                                 temp_var_thr, filter_hws, n) {
  i          <- end
  last_valid <- end
  peak_base  <- end

  while (i < next_start) {
    ri <- i + 1L
    if (ri < 1L || ri > n) break
    valid <- .check_valid_border_mod(act_med_low, temperature, temperature_threshold,
                                     act_median, activity_thr, ri)
    if (valid) {
      if (ntv[ri] >= temp_var_thr) { peak_base <- i; break }
      i <- i + 1L
    } else { last_valid <- i; break }
  }

  if (last_valid == end) {
    if (peak_base > end) {
      peak_ri <- peak_base + 1L
      while (peak_ri < n && ntv[peak_ri] >= ntv[peak_ri - 1L]) peak_ri <- peak_ri + 1L
      return(peak_ri - 1L)
    } else return(NULL)
  } else return(last_valid)
}

#' Border validity check (uses R 1-indexed position ri)
#' @noRd
.check_valid_border_mod <- function(act_med_low, temperature,
                                     temperature_threshold, act_median,
                                     activity_thr, ri) {
  act_med_low[ri] == 1L ||
    (temperature[ri] < temperature_threshold && act_median[ri] < 2 * activity_thr)
}

#' Estimate sleep periods for the sleep filter
#' @noRd
.estimate_sleep <- function(activity, valid_activity, n, short_only_mask,
                             sleep_act_hws, sleep_all_q, sleep_pos_q,
                             temperature, temperature_threshold,
                             sleep_low_temp_prop_max, epoch_hour) {
  act_zp    <- zero_prop(valid_activity)
  all_q_idx <- act_zp + (1 - act_zp) * sleep_all_q
  pos_q_idx <- act_zp + (1 - act_zp) * sleep_pos_q
  thr_all   <- as.double(stats::quantile(valid_activity, all_q_idx, names = FALSE))
  thr_pos   <- as.double(stats::quantile(valid_activity, pos_q_idx, names = FALSE))
  sleep_thr <- max(thr_all, thr_pos)

  act_med_sleep <- median_filter(activity, sleep_act_hws)
  is_sleep      <- as.integer(act_med_sleep < sleep_thr)

  sleep_periods <- .rle_periods(is_sleep == 1L)
  min_sleep_len <- as.integer(epoch_hour)
  if (nrow(sleep_periods) > 0L) {
    sleep_periods <- sleep_periods[
      (sleep_periods$end - sleep_periods$start) >= min_sleep_len, ]
    is_sleep <- integer(n)
    for (i in seq_len(nrow(sleep_periods))) {
      s <- sleep_periods$start[i]; e <- sleep_periods$end[i]
      if (s < e) is_sleep[(s + 1L):e] <- 1L
    }
  }

  sleep_periods <- .rle_periods(is_sleep == 1L)
  if (nrow(sleep_periods) > 0L) {
    valid_sleep <- vapply(seq_len(nrow(sleep_periods)), function(i) {
      s <- sleep_periods$start[i]; e <- sleep_periods$end[i]
      below_prop(temperature[(s + 1L):e], temperature_threshold) <= sleep_low_temp_prop_max
    }, logical(1))
    sleep_periods <- sleep_periods[valid_sleep, ]
    is_sleep <- integer(n)
    for (i in seq_len(nrow(sleep_periods))) {
      s <- sleep_periods$start[i]; e <- sleep_periods$end[i]
      if (s < e) is_sleep[(s + 1L):e] <- 1L
    }
  }
  is_sleep
}

#' Compute forbidden zones at the centre of sleep periods
#' @noRd
.compute_forbidden_zone <- function(estimated_sleep, epoch_hour, n) {
  forbidden     <- integer(n)
  sleep_periods <- .rle_periods(estimated_sleep == 1L)
  if (nrow(sleep_periods) == 0L) return(forbidden)

  for (i in seq_len(nrow(sleep_periods))) {
    s   <- sleep_periods$start[i]; e <- sleep_periods$end[i]
    len <- e - s  # Python-style length
    if (len > 0L) {
      q1_r <- s + as.integer(len * 0.25) + 1L  # R index
      q3_r <- s + as.integer(len * 0.75)        # R index (exclusive boundary)
      if (q1_r <= q3_r && q1_r >= 1L && q3_r <= n)
        forbidden[q1_r:q3_r] <- 1L
    }
  }
  forbidden
}

#' Check if the temperature distribution is bimodal
#' Matches Python's check_bimodality() exactly:
#' starts True; only False if:
#'   1. activity_zero_proportion >= 0.8
#'   2. ashman_d <= ashman_d_min AND offwrist_prop >= bimodal_max AND low_pos_act_prop >= bimodal_min
#' @noRd
.check_bimodality <- function(initial_offwrist, activity, activity_thr,
                               ashman, ashman_d_min, ashman_d_max,
                               bimodal_max_offwrist_p, bimodal_min_low_act_p, n) {
  # Condition 1: mostly zero activity -> unimodal
  if (zero_prop(activity) >= 0.8) return(FALSE)

  # Condition 2: only tested when ashman_d is low
  if (ashman <= ashman_d_min) {
    offwrist_prop    <- mean(initial_offwrist == 0L)
    pos_act          <- activity[activity > 0]
    low_pos_act_prop <- if (length(pos_act) > 0) mean(pos_act < activity_thr) else 0
    if (offwrist_prop >= bimodal_max_offwrist_p &&
        low_pos_act_prop >= bimodal_min_low_act_p) return(FALSE)
  }

  TRUE
}

#' Valley-peak off-wrist detection from temperature derivative
#' Returns periods in 0-indexed exclusive convention
#' @noRd
.valley_peak_algorithm <- function(
    temp_derivative, temperature, temperature_threshold, estimated_sleep,
    valley_quantile, peak_quantile, minimum_offwrist_length,
    short_offwrist_length, short_vp_decrease_ratio_min,
    valley_peak_low_temp_min, epoch_hour, n,
    next_possible_length, short_criteria, medium_criteria
) {
  valley_thr <- as.double(stats::quantile(temp_derivative, 1 - valley_quantile, names = FALSE))
  peak_thr   <- as.double(stats::quantile(temp_derivative, peak_quantile, names = FALSE))

  valley_idx <- which(temp_derivative <= valley_thr)  # 1-indexed R positions
  peak_idx   <- which(temp_derivative >= peak_thr)

  if (length(valley_idx) == 0L || length(peak_idx) == 0L)
    return(data.frame(start = integer(), end = integer(), valley_peak = logical()))

  candidates <- list()

  for (v_ri in valley_idx) {
    following_peaks <- peak_idx[peak_idx > v_ri]
    if (length(following_peaks) == 0L) next
    p_ri      <- following_peaks[1L]
    length_vp <- p_ri - v_ri  # matches Python: both 1-indexed, difference = Python length

    if (length_vp < minimum_offwrist_length) {
      if (next_possible_length) {
        next_valleys <- valley_idx[valley_idx > v_ri & valley_idx < p_ri]
        if (length(next_valleys) > 0L) {
          v2_ri <- next_valleys[length(next_valleys)]
          if ((p_ri - v2_ri) >= minimum_offwrist_length) {
            v_ri <- v2_ri; length_vp <- p_ri - v_ri
          }
        }
      }
      if (length_vp < minimum_offwrist_length) next
    }

    low_temp_p    <- below_prop(temperature[v_ri:p_ri], temperature_threshold)
    if (low_temp_p < valley_peak_low_temp_min) next

    sleep_overlap <- mean(estimated_sleep[v_ri:p_ri] == 1L)
    if (sleep_overlap > 0.4) next

    if (length_vp < short_offwrist_length && short_criteria) {
      decrease_ratio <- .compute_decrease_ratio(temp_derivative, v_ri, p_ri)
      if (decrease_ratio < short_vp_decrease_ratio_min) next
    }

    # Store as 0-indexed [start, end) — start = v_ri-1, end = p_ri
    candidates[[length(candidates) + 1L]] <- c(v_ri - 1L, p_ri)
  }

  if (length(candidates) == 0L)
    return(data.frame(start = integer(), end = integer(), valley_peak = logical()))

  m <- do.call(rbind, candidates)
  data.frame(start = m[, 1L], end = m[, 2L], valley_peak = TRUE)
}

#' Compute temperature derivative decrease ratio (1-indexed R positions)
#' @noRd
.compute_decrease_ratio <- function(temp_derivative, v_ri, p_ri) {
  seg <- temp_derivative[v_ri:p_ri]
  neg <- seg[seg < 0]; pos <- seg[seg > 0]
  if (length(pos) == 0L || sum(pos) == 0) return(0)
  abs(sum(neg)) / sum(pos)
}

#' Near-all-off detection for unimodal recordings
#' @noRd
.near_all_off_detection <- function(activity, activity_thr,
                                     positive_act_median, n, epoch_hour) {
  act_zp <- zero_prop(activity)
  if (act_zp > 0.9 || positive_act_median < activity_thr)
    return(data.frame(start = 0L, end = n, length = n, valley_peak = FALSE))
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
    activity_zero_prop = numeric(np), low_act_prop       = numeric(np),
    high_act_before    = numeric(np), high_act_after     = numeric(np),
    start_act_weight   = numeric(np), end_act_weight     = numeric(np),
    border_act_conc    = numeric(np), start_temp_weight  = numeric(np),
    end_temp_weight    = numeric(np), border_temp_conc   = numeric(np),
    low_temp_prop      = numeric(np), temp_dif_median    = numeric(np),
    temp_dif_variance  = numeric(np)
  )

  for (i in seq_len(np)) {
    s  <- offwrist_periods$start[i]        # 0-indexed
    e  <- min(offwrist_periods$end[i], n)  # exclusive
    rs <- s + 1L; re <- e                  # R indices

    if (rs > re || rs < 1L || re > n) next

    seg_act  <- activity[rs:re]
    seg_temp <- temperature[rs:re]
    seg_tv   <- temperature_variance[rs:re]
    seg_dif  <- dif_temp[rs:re]
    seg_dv   <- dif_temp_var[rs:re]

    act_segs <- .segmentation(seg_act, segments)
    tv_segs  <- .segmentation(seg_tv,  segments)

    before_s <- max(1L, rs - window)
    after_e  <- min(n, re + window)

    report$activity_zero_prop[i] <- zero_prop(seg_act)
    report$low_act_prop[i]       <- below_prop(seg_act, activity_thr)
    if (rs > 1L) report$high_act_before[i] <- 1 - below_prop(activity[before_s:(rs - 1L)], act_thr2)
    if (re < n)  report$high_act_after[i]  <- 1 - below_prop(activity[(re + 1L):after_e], act_thr2)
    report$start_act_weight[i]   <- act_segs[1L]
    report$end_act_weight[i]     <- act_segs[segments]
    report$start_temp_weight[i]  <- tv_segs[1L]
    report$end_temp_weight[i]    <- tv_segs[segments]
    report$low_temp_prop[i]      <- below_prop(seg_temp, temperature_threshold)
    report$temp_dif_median[i]    <- stats::median(seg_dif, na.rm = TRUE)
    report$temp_dif_variance[i]  <- stats::median(seg_dv,  na.rm = TRUE)
  }

  report$border_act_conc  <- report$start_act_weight + report$end_act_weight
  report$border_temp_conc <- report$start_temp_weight + report$end_temp_weight
  report
}

#' Segment a vector and return proportion of total in each segment
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

#' Description-report-based filter
#' Matches Python's third_stage_refinement description filter exactly.
#' Only short periods (length < long_offwrist_length) are filtered.
#' Valley-peak short periods are exempt from most filters.
#' @noRd
.description_report_filter <- function(periods, report, report_zero_act_min,
                                        border_conc_min, act_around_min,
                                        low_act_prop_min, temp_dif_min, temp_dif_max,
                                        long_offwrist_length, short_offwrist_length,
                                        is_highly_separable = FALSE) {
  if (nrow(periods) == 0L || nrow(report) == 0L) return(periods)

  keep <- rep(TRUE, nrow(periods))

  for (i in seq_len(nrow(periods))) {
    r   <- report[i, ]
    len <- periods$length[i]
    vp  <- if ("valley_peak" %in% names(periods)) periods$valley_peak[i] else FALSE

    # Only filter short periods (long ones are always kept)
    if (len >= long_offwrist_length || is_highly_separable) next

    is_valid <- TRUE

    # Zero activity proportion filter (with low_activity_after check)
    # do_report_low_activity_after = TRUE in ActTrust -> only filter if ALSO low_act_prop < min
    if (r$activity_zero_prop < report_zero_act_min) {
      if (r$low_act_prop < low_act_prop_min) {
        is_valid <- FALSE
      }
    }

    # Skip remaining filters for short valley-peak periods
    if (is_valid && !(vp && len <= short_offwrist_length)) {

      # Border concentration: both must be low
      if (is_valid && r$border_act_conc < 0.35 && r$border_temp_conc < 0.35) {
        is_valid <- FALSE
      }

      # Activity around: both must be low (< 0.2)
      if (is_valid && r$high_act_before < 0.2 && r$high_act_after < 0.2) {
        is_valid <- FALSE
      }

      # Border activity filter: (before OR after < min) AND (both concentrations < min)
      if (is_valid &&
          (r$high_act_before < act_around_min || r$high_act_after < act_around_min) &&
          r$border_act_conc  < border_conc_min &&
          r$border_temp_conc < border_conc_min) {
        is_valid <- FALSE
      }

      # Temperature difference filter: remove if dif_median > min threshold
      if (is_valid && r$temp_dif_median > temp_dif_min) {
        is_valid <- FALSE
      }
    }

    if (!is_valid) keep[i] <- FALSE
  }

  periods[keep, , drop = FALSE]
}

#' Remove on-wrist periods surrounded by off-wrist on both sides
#' @noRd
.surrounded_onwrist_filter <- function(periods, n) {
  if (nrow(periods) < 2L) return(periods)

  offwrist_bin <- integer(n)
  for (i in seq_len(nrow(periods))) {
    s <- periods$start[i]; e <- min(periods$end[i], n)
    if (s < e) offwrist_bin[(s + 1L):e] <- 1L
  }

  onwrist_periods <- .rle_periods(offwrist_bin == 0L)
  if (nrow(onwrist_periods) == 0L) return(periods)

  remove_onwrist <- vapply(seq_len(nrow(onwrist_periods)), function(i) {
    len <- onwrist_periods$end[i] - onwrist_periods$start[i]  # Python-style length
    s   <- onwrist_periods$start[i]
    e   <- onwrist_periods$end[i]
    preceded  <- any(periods$end   <= s)
    succeeded <- any(periods$start >= e)
    len < 20L && preceded && succeeded
  }, logical(1))

  if (!any(remove_onwrist)) return(periods)

  for (j in which(remove_onwrist)) {
    s    <- onwrist_periods$start[j]
    e    <- onwrist_periods$end[j]
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
