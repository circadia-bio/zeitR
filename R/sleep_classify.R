# ── Vallim pipeline: native sleep episode extraction and classification ─────────
# R port of Julia Ribeiro da Silva Vallim's (JRSV) condor_pipeline post-processing
# rules. Codename: Vallim.
# Source: pipeline_functions_fix27.py + vs_condor_py_pipeline_fix29_jrsv.ipynb
#
# Fixes ported:
#   Fix 25  — exclude truncated episodes at recording end
#   Fix 26a — adaptive nocturnal window (infer_nocturnal_window)
#   Fix 26b — sleep-date collision fix (fix_sleep_date_collision)
#   Fix 26c — fragmented episode recovery (recover_fragmented_episodes)
#   Fix 27  — execution order audit + 14 h TBT ceiling
#   Fix 29  — corrected split/exclude logic for long episodes
#
# Execution order in classify_sleep_episodes() follows Fix 27:
#   Fix 25 -> Fix 26a -> Fix 29 / Rule 2 -> Fix 26c ->
#   Rules 3-5 -> Fix 26b -> Rule 6 -> Rule 7

# ══════════════════════════════════════════════════════════════════════════════
# Exported: extract_sleep_episodes
# ══════════════════════════════════════════════════════════════════════════════

#' Extract sleep episodes from a CSPD-scored epoch table
#'
#' Converts the epoch-level `state` vector produced by [detect_sleep_crespo()]
#' into a per-episode tibble analogous to the Python SleepPipeline's
#' `results.nights`. Contiguous on-wrist sleep periods are delimited by
#' `.nights_df()`; per-epoch wake/sleep within each period is scored by
#' [score_epochs_cole_kripke()] to derive WASO, SOL, SOI, TST, EFF, and NW.
#' All time metrics are returned in **minutes**.
#'
#' @param data A tibble as returned by [detect_sleep_crespo()], containing
#'   columns `datetime`, `ZCMn`, and `state`.
#' @param wake_thresh `integer(1)`. Minimum wake-bout length (epochs) used
#'   by `.nights_df()` to separate episode boundaries. Default `60L`.
#'
#' @return A tibble with one row per detected episode and columns `bts`
#'   (POSIXct), `gts` (POSIXct), `tbt`, `tst`, `sol`, `soi`, `waso`
#'   (minutes), `nw` (integer), `eff` (0-1), `nap` (logical, always `FALSE`
#'   -- classification into main/secondary happens in
#'   [classify_sleep_episodes()]).
#'
#' @seealso [classify_sleep_episodes()], [run_pipeline_native()]
#' @export
extract_sleep_episodes <- function(data, wake_thresh = 60L) {
  required <- c("datetime", "ZCMn", "state")
  missing  <- setdiff(required, names(data))
  if (length(missing) > 0L)
    zeitr_abort("{.arg data} is missing required column(s): {.val {missing}}.")

  state     <- as.integer(data$state)
  datetimes <- as.POSIXct(data$datetime)
  zcm       <- as.double(data$ZCMn)

  onwrist   <- state != 4L
  stamps_ow <- datetimes[onwrist]
  zcm_ow    <- zcm[onwrist]
  state_ow  <- state[onwrist]

  # Epoch duration in minutes (mode of on-wrist inter-epoch gaps)
  secs_ow   <- as.numeric(stamps_ow)
  dd        <- diff(secs_ow)
  dd        <- dd[dd > 0]
  epoch_s   <- as.numeric(names(sort(table(dd), decreasing = TRUE))[1L])
  epoch_min <- epoch_s / 60.0

  nd <- .nights_df(state_ow, wake_thresh = as.integer(wake_thresh))

  if (nrow(nd) == 0L) return(.empty_episodes_tbl())

  out <- vector("list", nrow(nd))

  for (i in seq_len(nrow(nd))) {
    bt0 <- nd$bt[i]
    gt0 <- nd$gt[i]
    idx <- (bt0 + 1L):gt0
    ck  <- score_epochs_cole_kripke(zcm_ow[idx])   # 1 = wake, 0 = sleep
    m   <- length(ck)

    lat <- 0L
    while (lat < m && ck[lat + 1L] > 0L) lat <- lat + 1L
    ine <- m - 1L
    while (ine > 0L && ck[ine + 1L] > 0L) ine <- ine - 1L

    sol  <- lat
    soi  <- (m - 1L) - ine
    waso <- if ((lat + 1L) <= ine) sum(ck[(lat + 1L):ine]) else 0L
    nw   <- sum(diff(ck) > 0L)
    tbt  <- gt0 - bt0
    tst  <- tbt - waso - soi - sol
    eff  <- if (tbt > 0L) tst / tbt else NA_real_

    out[[i]] <- list(
      bts  = stamps_ow[bt0 + 1L],   # first sleep epoch
      gts  = stamps_ow[gt0],         # last sleep epoch (Julia's convention)
      tbt  = tbt  * epoch_min,
      tst  = max(0.0, tst * epoch_min),
      sol  = sol  * epoch_min,
      soi  = soi  * epoch_min,
      waso = waso * epoch_min,
      nw   = as.integer(nw),
      eff  = as.double(eff),
      nap  = nd$nap[i]
    )
  }

  tibble::tibble(
    bts  = do.call(c, lapply(out, `[[`, "bts")),
    gts  = do.call(c, lapply(out, `[[`, "gts")),
    tbt  = vapply(out, `[[`, double(1L),  "tbt"),
    tst  = vapply(out, `[[`, double(1L),  "tst"),
    sol  = vapply(out, `[[`, double(1L),  "sol"),
    soi  = vapply(out, `[[`, double(1L),  "soi"),
    waso = vapply(out, `[[`, double(1L),  "waso"),
    nw   = vapply(out, `[[`, integer(1L), "nw"),
    eff  = vapply(out, `[[`, double(1L),  "eff"),
    nap  = vapply(out, `[[`, logical(1L), "nap")
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# Exported: classify_sleep_episodes
# ══════════════════════════════════════════════════════════════════════════════

#' Classify sleep episodes as main or secondary (native pipeline)
#'
#' Applies the JRSV classification rule set to a raw episode table from
#' [extract_sleep_episodes()], assigning each episode a `sleep_type` of
#' `"main"` or `"secondary"`. The execution order follows Fix 27:
#'
#' \enumerate{
#'   \item **Fix 25** -- exclude truncated episodes at the recording end.
#'   \item **Fix 26a** -- infer adaptive nocturnal window from `int_temp`
#'     and `light`.
#'   \item **Fix 29** -- split TBT 14-16 h episodes first; exclude TBT > 16 h
#'     directly without attempting a split.
#'   \item **Fix 26c** -- recover fragmented sleep nights missed by the scorer.
#'   \item **Rules 3-5** -- classify as main or secondary using the nocturnal
#'     window and `min_main_tib_h`.
#'   \item **Fix 26b** -- resolve sleep-date collisions from the noon threshold.
#'   \item **Rule 6** -- keep the longest main episode per sleep date.
#'   \item **Rule 7** -- exclude all episodes on days with no main sleep.
#' }
#'
#' @param episodes A tibble as returned by [extract_sleep_episodes()].
#' @param data A tibble as returned by [detect_sleep_crespo()], containing
#'   at minimum `datetime`, `state`, `activity`. The ActTrust channels
#'   `int_temp` (wrist temperature) and `light` (ambient lux) must be
#'   present for Fix 26a/c; the function falls back to defaults when absent.
#' @param max_tib_h `numeric(1)`. Episodes with TBT > this value are excluded
#'   directly. Default `16`.
#' @param max_main_tib_h `numeric(1)`. Episodes with TBT in
#'   `(max_main_tib_h, max_tib_h]` are split first. Default `14`.
#' @param min_main_tib_h `numeric(1)`. Minimum TBT (hours) for a main
#'   episode. Default `4`.
#' @param nocturnal_onset_start `numeric(1)`. Default nocturnal window start
#'   (decimal hours). Overridden by Fix 26a. Default `18`.
#' @param nocturnal_onset_end `numeric(1)`. Default nocturnal window end.
#'   Overridden by Fix 26a. Default `6`.
#' @param temp_thresh `numeric(1)`. Minimum wrist temperature (degC) for
#'   candidate episodes in Fix 26a/c. Default `28`.
#' @param light_thresh_window `numeric(1)`. Maximum ambient light (lux) for
#'   the nocturnal window inference (Fix 26a). Default `5`.
#' @param light_thresh_recovery `numeric(1)`. Maximum ambient light (lux) for
#'   gap merging in fragment recovery (Fix 26c). Default `10`.
#' @param min_fragment_tib_h `numeric(1)`. Minimum TBT (hours) for each
#'   fragment produced by an episode split (Fix 29). Default `1`.
#' @param rolling_window_min `integer(1)`. Smoothing window (epochs) for the
#'   activity signal used to locate episode split points. Default `15L`.
#' @param max_split_iterations `integer(1)`. Maximum recursive splits per
#'   episode. Default `3L`.
#' @param collision_gap_h `numeric(1)`. Minimum gap (hours) between two
#'   same-date main episodes to trigger sleep-date reassignment (Fix 26b).
#'   Default `4`.
#' @param verbose `logical(1)`. Print step-by-step diagnostics. Default
#'   `FALSE`.
#'
#' @return A tibble with the same columns as `episodes` plus `sleep_type`
#'   (`"main"` or `"secondary"`) and `is_nap` (logical, `TRUE` when
#'   `sleep_type == "secondary"`, for backwards compatibility with the
#'   `zeitr_result$nights` schema).
#'
#' @seealso [extract_sleep_episodes()], [run_pipeline_native()]
#' @export
classify_sleep_episodes <- function(
    episodes,
    data,
    max_tib_h             = 16.0,
    max_main_tib_h        = 14.0,
    min_main_tib_h        = 4.0,
    nocturnal_onset_start = 18.0,
    nocturnal_onset_end   = 6.0,
    temp_thresh           = 28.0,
    light_thresh_window   = 5.0,
    light_thresh_recovery = 10.0,
    min_fragment_tib_h    = 1.0,
    rolling_window_min    = 15L,
    max_split_iterations  = 3L,
    collision_gap_h       = 4.0,
    verbose               = FALSE
) {
  if (nrow(episodes) == 0L) {
    return(tibble::add_column(episodes,
                              sleep_type = character(0L),
                              is_nap     = logical(0L)))
  }

  ep <- as.data.frame(episodes, stringsAsFactors = FALSE)
  ep$bts <- as.POSIXct(ep$bts)
  ep$gts <- as.POSIXct(ep$gts)

  dt     <- as.POSIXct(data$datetime)
  ow     <- as.integer(data$state) != 4L
  secs_v <- as.numeric(dt[ow])
  dd_v   <- diff(secs_v); dd_v <- dd_v[dd_v > 0]
  epoch_s   <- as.numeric(names(sort(table(dd_v), decreasing = TRUE))[1L])
  epoch_min <- epoch_s / 60.0
  tz        <- attr(dt[1L], "tzone") %||% "UTC"

  # ── Fix 25: exclude truncated episodes at recording end ──────────────────
  last_date      <- as.Date(max(dt), tz = tz)
  last_day_noon  <- as.POSIXct(paste0(format(last_date), " 12:00:00"), tz = tz)
  keep           <- ep$bts < last_day_noon
  n_trunc        <- sum(!keep)
  if (n_trunc > 0L) {
    if (verbose) cli::cli_inform("  [Fix 25] Removed {n_trunc} truncated episode(s) at recording end.")
    ep <- ep[keep, , drop = FALSE]
  }
  if (nrow(ep) == 0L) return(.empty_classified())

  # ── Fix 26a: infer adaptive nocturnal window ──────────────────────────────
  noc <- .infer_nocturnal_window(
    data = data, episodes = ep,
    temp_thresh   = temp_thresh,
    light_thresh  = light_thresh_window,
    default_start = nocturnal_onset_start,
    default_end   = nocturnal_onset_end,
    verbose       = verbose
  )
  noc_start <- noc[1L]; noc_end <- noc[2L]

  # ── Fix 29 (Rule 1): split TBT 14-16 h; exclude >16 h directly ──────────
  tbt_h    <- ep$tbt / 60.0
  in_split <- tbt_h > max_main_tib_h & tbt_h <= max_tib_h

  if (any(in_split)) {
    new_rows <- vector("list", nrow(ep))
    for (i in seq_len(nrow(ep))) {
      if (in_split[i]) {
        frags <- .split_recursive(
          ep[i, , drop = FALSE], data, epoch_min,
          min_frag_h = min_fragment_tib_h,
          max_tib_h  = max_main_tib_h,
          noc_start  = noc_start, noc_end = noc_end,
          rolling_min = as.integer(rolling_window_min),
          max_iter    = as.integer(max_split_iterations),
          verbose     = verbose
        )
        new_rows[[i]] <- do.call(rbind, frags)
      } else {
        new_rows[[i]] <- ep[i, , drop = FALSE]
      }
    }
    ep <- do.call(rbind, new_rows)
    rownames(ep) <- NULL

    before  <- nrow(ep)
    ep      <- ep[ep$tbt / 60.0 <= max_main_tib_h, , drop = FALSE]
    dropped <- before - nrow(ep)
    if (dropped > 0L && verbose)
      cli::cli_inform("  [Fix 29] Excluded {dropped} episode(s) still > {max_main_tib_h} h after split.")
  }

  # Rule 2: exclude TBT > max_tib_h directly (biologically implausible)
  before  <- nrow(ep)
  ep      <- ep[ep$tbt / 60.0 <= max_tib_h, , drop = FALSE]
  dropped <- before - nrow(ep)
  if (dropped > 0L && verbose)
    cli::cli_inform("  [Rule 2] Excluded {dropped} episode(s) with TBT > {max_tib_h} h.")

  if (nrow(ep) == 0L) return(.empty_classified())

  # Sleep date (noon threshold) -- needed before Fix 26c for covered-dates check
  ep$sleep_date <- .noon_date(ep$bts, tz = tz)

  # Preliminary classification for Fix 26c covered-dates calculation
  bts_h_pre  <- .decimal_h(ep$bts)
  in_noc_pre <- (bts_h_pre >= noc_start) | (bts_h_pre < noc_end)
  ep$sleep_type <- ifelse(in_noc_pre & (ep$tbt / 60.0 >= min_main_tib_h), "main", "secondary")

  # ── Fix 26c: recover fragmented episodes ─────────────────────────────────
  ep <- .recover_fragmented_episodes(
    data         = data, episodes = ep,
    epoch_min    = epoch_min,
    temp_thresh  = temp_thresh,
    light_thresh = light_thresh_recovery,
    min_tib_h    = 3.0,
    noc_start    = noc_start, noc_end = noc_end,
    tz           = tz, verbose = verbose
  )
  ep$bts <- as.POSIXct(ep$bts)
  ep$gts <- as.POSIXct(ep$gts)

  # ── Rules 3-5: final classification ──────────────────────────────────────
  bts_h  <- .decimal_h(ep$bts)
  in_noc <- (bts_h >= noc_start) | (bts_h < noc_end)
  ep$sleep_type <- "secondary"
  ep$sleep_type[in_noc & (ep$tbt / 60.0 >= min_main_tib_h)] <- "main"

  # Re-stamp sleep_date (covers recovered episodes and any type changes)
  ep$sleep_date <- .noon_date(ep$bts, tz = tz)

  # ── Fix 26b: resolve sleep-date collisions ────────────────────────────────
  ep <- .fix_sleep_date_collision(ep, gap_h = collision_gap_h, verbose = verbose)

  # ── Rule 6: keep longest main per sleep_date ─────────────────────────────
  ep <- .keep_longest_main(ep, verbose = verbose)

  # ── Rule 7: exclude days with no main sleep ───────────────────────────────
  days_with_main <- unique(ep$sleep_date[ep$sleep_type == "main"])
  excluded_days  <- setdiff(unique(ep$sleep_date), days_with_main)
  if (length(excluded_days) > 0L) {
    if (verbose) {
      for (d in excluded_days)
        cli::cli_inform("  [Rule 7] Excluded sleep-date {format(d)} -- no valid main episode.")
    }
    ep <- ep[ep$sleep_date %in% days_with_main, , drop = FALSE]
  }

  ep$is_nap     <- ep$sleep_type != "main"
  ep$sleep_date <- NULL
  rownames(ep)  <- NULL
  tibble::as_tibble(ep)
}


# ══════════════════════════════════════════════════════════════════════════════
# Internal: Fix 26a -- infer adaptive nocturnal window
# ══════════════════════════════════════════════════════════════════════════════

#' @noRd
.infer_nocturnal_window <- function(
    data, episodes,
    temp_thresh   = 28.0,
    light_thresh  = 5.0,
    min_tib_h     = 3.0,
    margin_h      = 4.0,
    outlier_h     = 6.0,
    default_start = 18.0,
    default_end   = 6.0,
    temp_col      = "int_temp",
    light_col     = "light",
    verbose       = FALSE
) {
  dt        <- as.POSIXct(data$datetime)
  has_temp  <- temp_col  %in% names(data)
  has_light <- light_col %in% names(data)

  cands <- episodes[episodes$tbt / 60.0 >= min_tib_h, , drop = FALSE]
  if (nrow(cands) == 0L) {
    if (verbose) cli::cli_inform("  [Fix 26a] No episodes >= {min_tib_h} h -> defaults ({default_start} h-{default_end} h).")
    return(c(default_start, default_end))
  }

  valid <- list()
  for (i in seq_len(nrow(cands))) {
    bts_i <- as.POSIXct(cands$bts[i])
    gts_i <- as.POSIXct(cands$gts[i])
    mask  <- dt >= bts_i & dt <= gts_i
    seg   <- data[mask, , drop = FALSE]
    if (nrow(seg) == 0L) next

    t_med <- if (has_temp)  stats::median(seg[[temp_col]],  na.rm = TRUE) else NA_real_
    l_med <- if (has_light) stats::median(seg[[light_col]], na.rm = TRUE) else NA_real_

    ok <- (!is.na(t_med) && t_med >= temp_thresh) &&
          (!is.na(l_med) && l_med <= light_thresh)
    if (ok) valid <- c(valid, list(list(onset_h = .decimal_h(bts_i))))
  }

  if (length(valid) < 2L) {
    if (verbose)
      cli::cli_inform("  [Fix 26a] {length(valid)}/{nrow(cands)} candidates passed TEMP+LIGHT filter -> defaults.")
    return(c(default_start, default_end))
  }

  onsets  <- vapply(valid, `[[`, double(1L), "onset_h")
  anchor0 <- .circ_mean_raw(onsets)

  diffs <- pmin(abs(onsets - anchor0), 24.0 - abs(onsets - anchor0))
  clean <- onsets[diffs <= outlier_h]

  if (length(clean) < 2L) {
    if (verbose) cli::cli_inform("  [Fix 26a] After outlier removal: {length(clean)} candidates -> defaults.")
    return(c(default_start, default_end))
  }

  anchor <- .circ_mean_raw(clean)
  start  <- (anchor - margin_h + 24.0) %% 24.0
  end    <- (anchor + margin_h)         %% 24.0

  if (verbose)
    cli::cli_inform("  [Fix 26a] anchor={round(anchor,1)} h -> nocturnal window {round(start,1)} h-{round(end,1)} h")

  c(start, end)
}


# ══════════════════════════════════════════════════════════════════════════════
# Internal: Fix 26b -- sleep-date collision fix
# ══════════════════════════════════════════════════════════════════════════════

#' @noRd
.fix_sleep_date_collision <- function(episodes, gap_h = 4.0, verbose = FALSE) {
  ep <- episodes[order(as.POSIXct(episodes$bts)), , drop = FALSE]
  n  <- nrow(ep)
  if (n <= 1L) return(ep)

  fixed <- 0L
  for (i in 2:n) {
    if (!identical(ep$sleep_date[i], ep$sleep_date[i - 1L])) next
    if (ep$sleep_type[i] != "main" || ep$sleep_type[i - 1L] != "main") next
    gap <- as.numeric(difftime(as.POSIXct(ep$bts[i]),
                               as.POSIXct(ep$gts[i - 1L]), units = "hours"))
    if (gap >= gap_h) {
      new_date <- as.Date(as.POSIXct(ep$bts[i]))
      if (verbose)
        cli::cli_inform("  [Fix 26b] {format(ep$sleep_date[i])} -> {format(new_date)} (gap={round(gap,1)} h)")
      ep$sleep_date[i] <- new_date
      fixed <- fixed + 1L
    }
  }
  if (fixed == 0L && verbose) cli::cli_inform("  [Fix 26b] No sleep-date collisions detected.")
  ep
}


# ══════════════════════════════════════════════════════════════════════════════
# Internal: Fix 26c -- recover fragmented episodes
# ══════════════════════════════════════════════════════════════════════════════

#' @noRd
.recover_fragmented_episodes <- function(
    data, episodes, epoch_min,
    temp_thresh  = 28.0,
    light_thresh = 10.0,
    min_tib_h    = 3.0,
    noc_start    = 18.0,
    noc_end      = 6.0,
    temp_col     = "int_temp",
    light_col    = "light",
    tz           = "UTC",
    verbose      = FALSE
) {
  dt        <- as.POSIXct(data$datetime)
  has_temp  <- temp_col  %in% names(data)
  has_light <- light_col %in% names(data)

  ep      <- episodes
  covered <- unique(ep$sleep_date[ep$sleep_type == "main"])

  first_date <- as.Date(min(dt), tz = tz)
  last_date  <- as.Date(max(dt), tz = tz)
  all_dates  <- seq(first_date, last_date, by = "day")

  new_eps <- list()

  for (i in seq_along(all_dates)) {
    sd <- all_dates[i]
    if (sd %in% covered) next

    win_start <- as.POSIXct(format(sd), tz = tz) + noc_start * 3600.0
    win_end   <- as.POSIXct(format(sd), tz = tz) + (24.0 + noc_end + 4.0) * 3600.0

    seg_mask <- dt >= win_start & dt <= win_end
    if (sum(seg_mask) < 60L) next

    seg_dt  <- dt[seg_mask]
    # Use Cole-Kripke epoch scoring on ZCMn (0=sleep, 1=wake) to mirror
    # Python's epoch-level state column. The CSPD state column is period-level
    # (all epochs inside a detected period are state=1) -- using it here would
    # treat entire 19+ h periods as one sleep run, incorrectly triggering recovery.
    seg_zcm  <- as.double(data$ZCMn)[seg_mask]
    ck_score <- score_epochs_cole_kripke(seg_zcm)
    is_sleep <- ck_score == 0L   # CK: 0=sleep, 1=wake
    if (!any(is_sleep)) next

    rle_res    <- rle(is_sleep)
    ends_r     <- cumsum(rle_res$lengths)
    starts_r   <- ends_r - rle_res$lengths + 1L
    sleep_runs <- which(rle_res$values)

    runs <- lapply(sleep_runs, function(k) {
      list(bts = seg_dt[starts_r[k]], gts = seg_dt[ends_r[k]])
    })
    if (length(runs) == 0L) next

    # Sliding-window merge
    merged <- list(runs[[1L]])
    if (length(runs) > 1L) {
      for (j in 2:length(runs)) {
        gap_mask <- dt >= merged[[length(merged)]]$gts & dt <= runs[[j]]$bts
        gap_seg  <- data[gap_mask, , drop = FALSE]

        t_ok <- has_temp && nrow(gap_seg) > 0L &&
          stats::median(gap_seg[[temp_col]], na.rm = TRUE) >= temp_thresh
        l_ok <- !has_light || nrow(gap_seg) == 0L ||
          stats::median(gap_seg[[light_col]], na.rm = TRUE) <= light_thresh

        if (t_ok && l_ok) {
          merged[[length(merged)]]$gts <- runs[[j]]$gts
        } else {
          merged <- c(merged, list(runs[[j]]))
        }
      }
    }

    # Inject first merged episode meeting min_tib_h
    for (ep_i in merged) {
      tbt_h <- as.numeric(difftime(ep_i$gts, ep_i$bts, units = "hours"))
      if (tbt_h < min_tib_h) next

      stats_i <- .episode_stats_from_state(data, ep_i$bts, ep_i$gts, epoch_min)
      if (is.null(stats_i)) next

      if (verbose)
        cli::cli_inform(
          "  [Fix 26c] Recovered: {format(ep_i$bts,'%d/%m %H:%M')} -> {format(ep_i$gts,'%d/%m %H:%M')} TBT={round(tbt_h,1)} h"
        )

      new_eps <- c(new_eps, list(data.frame(
        bts        = ep_i$bts,
        gts        = ep_i$gts,
        tbt        = stats_i$tbt,
        tst        = stats_i$tst,
        sol        = stats_i$sol,
        soi        = stats_i$soi,
        waso       = stats_i$waso,
        nw         = stats_i$nw,
        eff        = stats_i$eff,
        nap        = FALSE,
        sleep_date = sd,
        sleep_type = "main",   # overwritten by Rules 3-5
        stringsAsFactors = FALSE
      )))
      break   # one recovery per sleep-date
    }
  }

  if (length(new_eps) == 0L) {
    if (verbose) cli::cli_inform("  [Fix 26c] No fragmented episodes recovered.")
    return(ep)
  }

  new_df <- do.call(rbind, new_eps)
  result <- rbind(ep, new_df)
  result <- result[order(as.POSIXct(result$bts)), , drop = FALSE]
  rownames(result) <- NULL
  if (verbose) cli::cli_inform("  [Fix 26c] {length(new_eps)} episode(s) recovered.")
  tibble::as_tibble(result)
}


# ══════════════════════════════════════════════════════════════════════════════
# Internal: episode stats from raw STATE
# ══════════════════════════════════════════════════════════════════════════════

#' Port of Julia's compute_episode_stats() -- uses raw STATE column directly
#' @noRd
.episode_stats_from_state <- function(data, bts, gts, epoch_min) {
  dt   <- as.POSIXct(data$datetime)
  mask <- dt >= as.POSIXct(bts) & dt <= as.POSIXct(gts)
  seg  <- data[mask, , drop = FALSE]
  if (nrow(seg) == 0L) return(NULL)

  tbt    <- nrow(seg)
  is_slp <- as.integer(seg$state) == 1L
  if (!any(is_slp)) return(NULL)

  sl_idx   <- which(is_slp)
  onset_i  <- sl_idx[1L]
  offset_i <- sl_idx[length(sl_idx)]

  sol <- onset_i - 1L
  soi <- tbt - offset_i

  win  <- is_slp[onset_i:offset_i]
  waso <- sum(!win)
  tst  <- sum(win)
  eff  <- if (tbt > 0L) round(tst / tbt, 4L) else NA_real_

  if (waso > 0L) {
    wr <- rle(!win)
    nw <- as.integer(sum(wr$values & wr$lengths >= 1L))
  } else {
    nw <- 0L
  }

  list(
    tbt  = tbt  * epoch_min,
    tst  = tst  * epoch_min,
    sol  = round(sol  * epoch_min, 1L),
    soi  = round(soi  * epoch_min, 1L),
    waso = waso * epoch_min,
    nw   = nw,
    eff  = as.double(eff)
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# Internal: Fix 29 -- episode splitting
# ══════════════════════════════════════════════════════════════════════════════

#' @noRd
.get_search_pct <- function(bts, noc_start, noc_end) {
  h <- .decimal_h(as.POSIXct(bts))
  if ((h >= noc_start) || (h < noc_end)) c(0.50, 0.90) else c(0.10, 0.50)
}

#' @noRd
.split_episode_by_activity <- function(row, data, epoch_min, min_frag_h,
                                        rolling_min = 15L,
                                        search_pct  = c(0.10, 0.50),
                                        verbose     = FALSE) {
  bts_ts  <- as.POSIXct(row$bts)
  gts_ts  <- as.POSIXct(row$gts)
  total_s <- as.numeric(difftime(gts_ts, bts_ts, units = "secs"))

  s_start <- bts_ts + total_s * search_pct[1L]
  s_end   <- bts_ts + total_s * search_pct[2L]

  dt   <- as.POSIXct(data$datetime)
  mask <- dt >= s_start & dt <= s_end & as.integer(data$state) != 4L

  if (!any(mask)) return(list(row))

  act      <- as.numeric(data$activity[mask])
  hws      <- max(1L, as.integer(rolling_min / 2L))
  smoothed <- mean_filter(act, hws)
  peak_i   <- which.max(smoothed)
  anchor   <- dt[mask][peak_i]

  tbt_before <- as.numeric(difftime(anchor, bts_ts, units = "mins"))
  tbt_after  <- as.numeric(difftime(gts_ts, anchor, units = "mins"))

  if (tbt_before < min_frag_h * 60.0 || tbt_after < min_frag_h * 60.0)
    return(list(row))

  r_bef     <- row; r_bef$gts <- anchor; r_bef$tbt <- tbt_before
  r_bef$sol <- 0.0; r_bef$soi <- 0.0
  r_aft     <- row; r_aft$bts <- anchor; r_aft$tbt <- tbt_after
  r_aft$sol <- 0.0; r_aft$soi <- 0.0

  if (verbose)
    cli::cli_inform(
      "  [Fix 29] Split at {format(anchor,'%d/%m %H:%M')} -> {round(tbt_before/60,1)} h + {round(tbt_after/60,1)} h"
    )

  list(r_bef, r_aft)
}

#' @noRd
.split_recursive <- function(row, data, epoch_min, min_frag_h, max_tib_h,
                              noc_start, noc_end, rolling_min = 15L,
                              iteration = 0L, max_iter = 3L, verbose = FALSE) {
  if (row$tbt / 60.0 <= max_tib_h || iteration >= max_iter) return(list(row))
  pct   <- .get_search_pct(row$bts, noc_start, noc_end)
  frags <- .split_episode_by_activity(row, data, epoch_min, min_frag_h,
                                      rolling_min = rolling_min,
                                      search_pct  = pct,
                                      verbose     = verbose)
  if (length(frags) == 1L) return(frags)
  result <- list()
  for (f in frags)
    result <- c(result, .split_recursive(f, data, epoch_min, min_frag_h, max_tib_h,
                                          noc_start, noc_end, rolling_min,
                                          iteration + 1L, max_iter, verbose))
  result
}


# ══════════════════════════════════════════════════════════════════════════════
# Internal: Rule 6 -- keep longest main per sleep_date
# ══════════════════════════════════════════════════════════════════════════════

#' @noRd
.keep_longest_main <- function(episodes, verbose = FALSE) {
  ep   <- episodes
  dups <- names(which(table(as.character(ep$sleep_date[ep$sleep_type == "main"])) > 1L))

  for (d in dups) {
    mains    <- which(ep$sleep_type == "main" & as.character(ep$sleep_date) == d)
    h_onset  <- .decimal_h(as.POSIXct(ep$bts[mains]))
    eve_idx  <- mains[h_onset >= 12.0]
    keep_idx <- if (length(eve_idx) >= 1L)
      eve_idx[which.max(ep$tbt[eve_idx])]
    else
      mains[which.max(ep$tbt[mains])]

    demote <- setdiff(mains, keep_idx)
    ep$sleep_type[demote] <- "secondary"
    if (verbose) cli::cli_inform("  [Rule 6] Demoted {length(demote)} duplicate main(s) on {d}.")
  }
  ep
}


# ══════════════════════════════════════════════════════════════════════════════
# Internal: shared utilities
# ══════════════════════════════════════════════════════════════════════════════

#' Decimal hours from a POSIXct scalar or vector
#' @noRd
.decimal_h <- function(x) {
  lt <- as.POSIXlt(x)
  lt$hour + lt$min / 60.0 + lt$sec / 3600.0
}

#' Noon-threshold sleep date: if bts hour >= 12, use bts date; else bts - 1 day
#' @noRd
.noon_date <- function(bts, tz = "UTC") {
  bts <- as.POSIXct(bts)
  h   <- as.integer(format(bts, "%H", tz = tz))
  d   <- as.Date(bts, tz = tz)
  as.Date(ifelse(h >= 12L, as.character(d), as.character(d - 1L)))
}

#' Internal circular mean (no NA handling; used for small validated vectors)
#' @noRd
.circ_mean_raw <- function(x) {
  theta <- x * (2.0 * pi / 24.0)
  ((atan2(mean(sin(theta)), mean(cos(theta))) * 24.0 / (2.0 * pi)) + 24.0) %% 24.0
}

#' Empty episodes tibble (zero rows, correct types)
#' @noRd
.empty_episodes_tbl <- function() {
  tibble::tibble(
    bts  = as.POSIXct(character(0L)),
    gts  = as.POSIXct(character(0L)),
    tbt  = double(0L),
    tst  = double(0L),
    sol  = double(0L),
    soi  = double(0L),
    waso = double(0L),
    nw   = integer(0L),
    eff  = double(0L),
    nap  = logical(0L)
  )
}

#' Empty classified tibble (zero rows, correct types)
#' @noRd
.empty_classified <- function() {
  tibble::tibble(
    bts        = as.POSIXct(character(0L)),
    gts        = as.POSIXct(character(0L)),
    tbt        = double(0L),
    tst        = double(0L),
    sol        = double(0L),
    soi        = double(0L),
    waso       = double(0L),
    nw         = integer(0L),
    eff        = double(0L),
    nap        = logical(0L),
    sleep_type = character(0L),
    is_nap     = logical(0L)
  )
}
