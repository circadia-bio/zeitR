#' Compute WASO and nightly sleep statistics
#'
#' Scores each epoch within detected sleep periods as wake or sleep using
#' [score_epochs_cole_kripke()], then computes per-night statistics. This is a
#' faithful port of the Python `condor_pipeline` `detect_waso`: night
#' boundaries are built by [.nights_df()] (the `search_gap = FALSE` path of the
#' reference `nights_df`), each night's ZCM is scored with Cole-Kripke, and the
#' epoch-level `state` is rebuilt from a zero (wake) base so that within-night
#' WASO-wake epochs and all epochs outside detected nights are scored as wake.
#'
#' The following metrics are computed for each detected night / nap:
#'
#' | Metric | Definition |
#' |--------|-----------|
#' | `tbt`  | Total Bed Time — epochs from bed time to get-up time |
#' | `tst`  | Total Sleep Time — `tbt - waso - sol - soi` |
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
#' @param search_gap `logical(1)`. Reserved for API compatibility with the
#'   reference pipeline. The reference `detect_waso` always builds boundaries
#'   with `search_gap = FALSE`, so this argument currently has no effect.
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

  wake_thresh <- as.integer(wake_thresh)

  state     <- as.integer(x$state)
  zcm       <- as.double(x$ZCMn)
  datetimes <- as.POSIXct(x$datetime)

  # ── On-wrist subset (exclude state == 4) ────────────────────────────────────
  # Python: onwrist = np.where(out == 4, False, True)
  onwrist <- state != 4L
  stamps  <- datetimes[onwrist]
  zcm_ow  <- zcm[onwrist]
  in_bed  <- state[onwrist]       # 1 = sleep period, 0 = wake, 7 = nap
  n_ow    <- sum(onwrist)

  # Python: state = np.zeros(n)  — everything defaults to wake (0)
  state_base <- numeric(n_ow)

  # ── Night boundaries (faithful nights_df, search_gap = FALSE) ───────────────
  nd       <- .nights_df(in_bed, wake_thresh = wake_thresh)
  n_nights <- nrow(nd)

  results <- vector("list", n_nights)

  for (i in seq_len(n_nights)) {
    bt0    <- nd$bt[i]            # 0-indexed first sleep/nap epoch
    gt0    <- nd$gt[i]            # 0-indexed first wake epoch after (exclusive)
    is_nap <- nd$nap[i]

    # Python zcm[bt:gt] (exclusive gt) -> R (bt0 + 1):gt0, length gt0 - bt0
    idx   <- (bt0 + 1L):gt0
    cpred <- score_epochs_cole_kripke(zcm_ow[idx])   # 1 = wake, 0 = sleep
    m     <- length(cpred)

    # Sleep onset latency — leading wake epochs
    latency <- 0L
    while (latency < m && cpred[latency + 1L] > 0L) latency <- latency + 1L

    # Sleep offset inertia — walk back from the last epoch while wake
    # (inertia is the 0-indexed position of the last sleep epoch)
    inertia <- m - 1L
    while (inertia > 0L && cpred[inertia + 1L] > 0L) inertia <- inertia - 1L

    sol  <- latency
    soi  <- (m - 1L) - inertia
    # waso = sum(cpred[latency:inertia]) — Python slice, exclusive of inertia
    waso <- if ((latency + 1L) <= inertia) sum(cpred[(latency + 1L):inertia]) else 0
    nw   <- sum(diff(cpred) > 0L)
    tbt  <- gt0 - bt0
    tst  <- tbt - waso - soi - sol
    eff  <- if (tbt > 0L) tst / tbt else NA_real_

    # Python: state[bt:gt] = (7 - 7*cpred) if nap else (1 - cpred)
    state_base[idx] <- if (is_nap) 7 - 7 * cpred else 1 - cpred

    results[[i]] <- list(
      night       = i,
      is_nap      = is_nap,
      bed_time    = stamps[bt0 + 1L],
      get_up_time = stamps[gt0 + 1L],
      tbt         = as.integer(tbt),
      tst         = as.double(tst),
      waso        = as.double(waso),
      sol         = as.integer(sol),
      soi         = as.integer(soi),
      nw          = as.integer(nw),
      eff         = as.double(eff)
    )
  }

  # ── Write scored states back (Python: out[onwrist] = state) ─────────────────
  full_state          <- state            # off-wrist (4) preserved
  full_state[onwrist] <- as.integer(state_base)
  x$state             <- full_state

  # Python: sleep = where(out == 4, 0, out); sleep = where(sleep == 7, 1, sleep)
  sleep_col <- full_state
  sleep_col[sleep_col == 4L] <- 0L
  sleep_col[sleep_col == 7L] <- 1L
  x$sleep   <- as.integer(sleep_col)

  if (n_nights == 0L) {
    zeitr_warn("No sleep periods detected.")
    nights_tbl <- tibble::tibble(
      night = integer(), is_nap = logical(),
      bed_time = stamps[0L], get_up_time = stamps[0L],
      tbt = integer(), tst = double(), waso = double(),
      sol = integer(), soi = integer(), nw = integer(), eff = double()
    )
  } else {
    nights_tbl <- tibble::tibble(
      night       = vapply(results, `[[`, integer(1), "night"),
      is_nap      = vapply(results, `[[`, logical(1), "is_nap"),
      bed_time    = do.call(c, lapply(results, `[[`, "bed_time")),
      get_up_time = do.call(c, lapply(results, `[[`, "get_up_time")),
      tbt         = vapply(results, `[[`, integer(1), "tbt"),
      tst         = vapply(results, `[[`, double(1),  "tst"),
      waso        = vapply(results, `[[`, double(1),  "waso"),
      sol         = vapply(results, `[[`, integer(1), "sol"),
      soi         = vapply(results, `[[`, integer(1), "soi"),
      nw          = vapply(results, `[[`, integer(1), "nw"),
      eff         = vapply(results, `[[`, double(1),  "eff")
    )
  }

  list(nights = nights_tbl, data = x)
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' Build night/nap boundaries from an on-wrist state vector
#'
#' Faithful port of the `search_gap = FALSE` (no-gap) path of the reference
#' `condor/nights_df.py`. Returns a data.frame with **0-indexed** `bt` (first
#' sleep/nap epoch) and `gt` (first wake epoch after the period, exclusive),
#' plus a logical `nap` flag (`TRUE` iff every epoch in `[bt, gt)` equals `7`).
#'
#' Mirrors the reference exactly: transition pairing (with a leading rising
#' edge prepended when the record starts asleep and a trailing falling edge
#' appended when it ends asleep), a wake-gap merge (`next_bt - prev_gt <
#' wake_thresh`), a length filter (`gt - bt >= sleep_thresh`, with an all-`7`
#' nap rescue when `gt - bt >= nap_thresh`), and a final wake-gap merge.
#'
#' @param states integer vector (0 = wake, 1 = sleep, 7 = nap), on-wrist.
#' @param wake_thresh minimum wake-bout length (epochs) to separate periods.
#' @param sleep_thresh minimum sleep-period length to keep (epochs).
#' @param nap_thresh minimum all-nap-period length to keep (epochs).
#' @noRd
.nights_df <- function(states, wake_thresh = 60L, sleep_thresh = 120L,
                       nap_thresh = 20L) {
  states <- as.integer(states)
  n      <- length(states)
  if (n == 0L) {
    return(data.frame(bt = integer(), gt = integer(), nap = logical()))
  }

  # Forward wake-gap merge: merges consecutive [bt, gt) pairs whose wake gap
  # (next_bt - current_gt) is below wake_thresh. Equivalent to the in-place
  # while-loop merge in the reference.
  .merge_gaps <- function(bt, gt) {
    if (length(bt) <= 1L) return(list(bt = bt, gt = gt))
    out_bt <- integer(0)
    out_gt <- integer(0)
    cur_bt <- bt[1L]
    cur_gt <- gt[1L]
    for (k in seq_along(bt)[-1L]) {
      if (bt[k] - cur_gt < wake_thresh) {
        cur_gt <- gt[k]
      } else {
        out_bt <- c(out_bt, cur_bt)
        out_gt <- c(out_gt, cur_gt)
        cur_bt <- bt[k]
        cur_gt <- gt[k]
      }
    }
    list(bt = c(out_bt, cur_bt), gt = c(out_gt, cur_gt))
  }

  # edges = concat([0], diff(states)); transitions at nonzero edges.
  edges <- c(0L, diff(states))
  pos0  <- which(edges != 0L) - 1L            # 0-indexed transition positions

  bt_keep <- integer(0)
  gt_keep <- integer(0)

  if (length(pos0) > 0L) {
    dir <- edges[pos0 + 1L]                    # transition direction (sign used)

    # Record starts asleep -> prepend a rising edge at index 0
    if (dir[1L] < 0L) {
      pos0 <- c(0L, pos0)
      dir  <- c(1L, dir)
    }
    # Record ends asleep -> append a falling edge at index n - 1
    if (dir[length(dir)] > 0L) {
      pos0 <- c(pos0, n - 1L)
      dir  <- c(dir, -1L)
    }

    # Rising edges open a period (bt); falling edges close it (gt, exclusive).
    bt_raw <- pos0[dir > 0L]
    gt_raw <- pos0[dir < 0L]

    merged <- .merge_gaps(bt_raw, gt_raw)

    # Length filter, with all-nap rescue.
    for (j in seq_along(merged$bt)) {
      blen <- merged$gt[j] - merged$bt[j]
      if (blen >= sleep_thresh) {
        bt_keep <- c(bt_keep, merged$bt[j])
        gt_keep <- c(gt_keep, merged$gt[j])
      } else if (blen >= nap_thresh) {
        seg <- states[(merged$bt[j] + 1L):merged$gt[j]]
        if (all(seg == 7L)) {
          bt_keep <- c(bt_keep, merged$bt[j])
          gt_keep <- c(gt_keep, merged$gt[j])
        }
      }
    }
  } else {
    # No transitions: the reference appends a single [0, n - 1] block.
    bt_keep <- 0L
    gt_keep <- n - 1L
  }

  # Final wake-gap merge on the kept boundaries.
  final <- .merge_gaps(bt_keep, gt_keep)
  bt    <- final$bt
  gt    <- final$gt

  nap <- vapply(seq_along(bt), function(j) {
    if (gt[j] <= bt[j]) return(FALSE)
    all(states[(bt[j] + 1L):gt[j]] == 7L)
  }, logical(1))

  data.frame(bt = bt, gt = gt, nap = nap)
}
