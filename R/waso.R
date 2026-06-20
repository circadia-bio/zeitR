#' Compute WASO and nightly sleep statistics
#'
#' Scores each epoch within detected sleep periods as wake or sleep using
#' [score_epochs_cole_kripke()], then computes per-night statistics.
#'
#' The following metrics are computed for each detected night / nap:
#'
#' | Metric | Definition |
#' |--------|-----------|
#' | `tbt`  | Total Bed Time — epochs from bed time to get-up time |
#' | `tst`  | Total Sleep Time — `tbt − waso − sol − soi` |
#' | `waso` | Wake After Sleep Onset — wake epochs between sleep onset and final wake |
#' | `sol`  | Sleep Onset Latency — epochs from bed time to first sleep epoch |
#' | `soi`  | Sleep Offset Inertia — trailing wake epochs at end of sleep period |
#' | `nw`   | Number of awakenings — count of wake-onset transitions |
#' | `eff`  | Sleep efficiency — `tst / tbt` |
#'
#' @param x A tibble as returned by [detect_naps_crespo()] (or
#'   [detect_sleep_crespo()] if nap detection is skipped), containing columns
#'   `datetime`, `ZCMn`, and `state`.
#' @param wake_thresh `integer(1)`. Minimum duration in epochs of a wake bout
#'   required to delimit a new sleep period boundary. Default is `60`.
#' @param search_gap `logical(1)`. If `TRUE`, allows a gap search between
#'   consecutive sleep periods when identifying night boundaries. Default is
#'   `FALSE`.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{`nights`}{A tibble with one row per detected night/nap and
#'       columns `night`, `is_nap`, `bed_time`, `get_up_time`, `tbt`, `tst`,
#'       `waso`, `sol`, `soi`, `nw`, `eff`.}
#'     \item{`data`}{The input tibble `x` with `state` and `sleep` updated
#'       with epoch-level Cole-Kripke wake/sleep scores.}
#'   }
#'
#' @references
#' Cole, R. J., Kripke, D. F., Gruen, W., Mullaney, D. J., & Gillin, J. C.
#' (1992). Automatic sleep/wake identification from wrist activity.
#' *Sleep*, 15(5), 461–469. \doi{10.1093/sleep/15.5.461}
#'
#' @export
#'
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' rec    <- read_acttrust("recordings/P001.txt")
#' prep   <- prepare_actigraphy(rec)
#' prep   <- detect_offwrist_bimodal(prep)
#' prep   <- detect_sleep_crespo(prep)
#' prep   <- detect_naps_crespo(prep)
#' result <- compute_waso(prep)
#' result$nights
#' }
compute_waso <- function(x, wake_thresh = 60L, search_gap = FALSE) {
  required <- c("datetime", "ZCMn", "state")
  missing  <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    zeitr_abort(
      "{.arg x} is missing required column(s): {.val {missing}}."
    )
  }

  state    <- as.integer(x$state)
  zcm      <- as.double(x$ZCMn)
  datetimes <- as.POSIXct(x$datetime)

  # On-wrist mask (exclude state == 4)
  onwrist  <- state != 4L
  stamps   <- datetimes[onwrist]
  zcm_ow   <- zcm[onwrist]
  state_ow <- state[onwrist]
  n_ow     <- sum(onwrist)

  # ── Identify sleep period boundaries ────────────────────────────────────────
  # A sleep period is a contiguous block of state %in% c(1, 7)
  # separated by >= wake_thresh wake epochs (state == 0)
  nights_info <- .find_sleep_periods(state_ow, wake_thresh = wake_thresh)

  if (nrow(nights_info) == 0L) {
    zeitr_warn("No sleep periods detected.")
    return(list(
      nights = tibble::tibble(
        night = integer(), is_nap = logical(),
        bed_time = as.POSIXct(character()), get_up_time = as.POSIXct(character()),
        tbt = integer(), tst = double(), waso = double(),
        sol = integer(), soi = integer(), nw = integer(), eff = double()
      ),
      data = x
    ))
  }

  n_nights <- nrow(nights_info)
  state_scored <- numeric(n_ow)

  results <- vector("list", n_nights)

  for (i in seq_len(n_nights)) {
    bt    <- nights_info$bt[i]
    gt    <- nights_info$gt[i]
    is_nap <- nights_info$is_nap[i]

    ck_pred <- score_epochs_cole_kripke(zcm_ow[bt:gt])  # 1 = wake, 0 = sleep

    # Sleep onset latency — leading wake epochs
    sol <- 0L
    while (sol < length(ck_pred) && ck_pred[sol + 1L] > 0L) sol <- sol + 1L

    # Sleep offset inertia — trailing wake epochs
    soi_end <- length(ck_pred)
    while (soi_end > 1L && ck_pred[soi_end] > 0L) soi_end <- soi_end - 1L
    soi <- length(ck_pred) - soi_end

    tbt  <- gt - bt + 1L
    waso <- sum(ck_pred[(sol + 1L):max(soi_end, sol + 1L)])
    nw   <- sum(diff(ck_pred) > 0L)
    tst  <- tbt - waso - sol - soi
    eff  <- if (tbt > 0L) tst / tbt else NA_real_

    # Encode scored state back
    scored <- if (is_nap) {
      ifelse(ck_pred == 0L, 7L, 0L)
    } else {
      ifelse(ck_pred == 0L, 1L, 0L)
    }
    state_scored[bt:gt] <- scored

    results[[i]] <- list(
      night       = i,
      is_nap      = is_nap,
      bed_time    = stamps[bt],
      get_up_time = stamps[gt],
      tbt         = tbt,
      tst         = tst,
      waso        = waso,
      sol         = sol,
      soi         = soi,
      nw          = nw,
      eff         = eff
    )
  }

  # ── Write scored states back to full-length vectors ─────────────────────────
  full_state         <- state
  full_state[onwrist] <- ifelse(state_scored == 0L, state[onwrist], state_scored)
  x$state            <- full_state

  sleep_col          <- ifelse(full_state == 4L, 0L,
                        ifelse(full_state %in% c(1L, 7L), 1L, 0L))
  x$sleep            <- sleep_col

  nights_tbl <- tibble::tibble(
    night       = vapply(results, `[[`, integer(1),  "night"),
    is_nap      = vapply(results, `[[`, logical(1),  "is_nap"),
    bed_time    = do.call(c, lapply(results, `[[`, "bed_time")),
    get_up_time = do.call(c, lapply(results, `[[`, "get_up_time")),
    tbt         = vapply(results, `[[`, integer(1),  "tbt"),
    tst         = vapply(results, `[[`, double(1),   "tst"),
    waso        = vapply(results, `[[`, double(1),   "waso"),
    sol         = vapply(results, `[[`, integer(1),  "sol"),
    soi         = vapply(results, `[[`, integer(1),  "soi"),
    nw          = vapply(results, `[[`, integer(1),  "nw"),
    eff         = vapply(results, `[[`, double(1),   "eff")
  )

  list(nights = nights_tbl, data = x)
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' Identify sleep period start/end indices from state vector
#' @param state integer vector (0 = wake, 1 = sleep, 4 = offwrist, 7 = nap)
#' @param wake_thresh minimum wake-bout length to separate sleep periods
#' @noRd
.find_sleep_periods <- function(state, wake_thresh = 60L) {
  in_sleep <- as.integer(state %in% c(1L, 7L))
  n        <- length(in_sleep)

  if (sum(in_sleep) == 0L) {
    return(data.frame(bt = integer(), gt = integer(), is_nap = logical()))
  }

  # RLE to find runs of sleep
  r      <- rle(in_sleep)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L

  sleep_runs <- which(r$values == 1L)
  if (length(sleep_runs) == 0L) {
    return(data.frame(bt = integer(), gt = integer(), is_nap = logical()))
  }

  # Merge runs separated by fewer than wake_thresh wake epochs
  bt_all <- starts[sleep_runs]
  gt_all <- ends[sleep_runs]

  merged_bt <- bt_all[1L]
  merged_gt <- gt_all[1L]
  bt_out    <- integer()
  gt_out    <- integer()

  for (i in seq_along(bt_all)[-1L]) {
    gap <- bt_all[i] - merged_gt - 1L
    if (gap < wake_thresh) {
      merged_gt <- gt_all[i]
    } else {
      bt_out <- c(bt_out, merged_bt)
      gt_out <- c(gt_out, merged_gt)
      merged_bt <- bt_all[i]
      merged_gt <- gt_all[i]
    }
  }
  bt_out <- c(bt_out, merged_bt)
  gt_out <- c(gt_out, merged_gt)

  # Determine if each period is a nap (majority state == 7)
  is_nap <- vapply(seq_along(bt_out), function(i) {
    seg <- state[bt_out[i]:gt_out[i]]
    mean(seg == 7L) > 0.5
  }, logical(1))

  data.frame(bt = bt_out, gt = gt_out, is_nap = is_nap)
}
