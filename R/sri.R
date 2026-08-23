#' Sleep Regularity Index (SRI)
#'
#' Computes the Sleep Regularity Index (Phillips et al. 2017): a measure of
#' day-to-day consistency in the sleep/wake pattern, based on the probability
#' that an epoch's sleep/wake state matches the state at the same clock time
#' exactly 24 h later (or earlier), averaged across the whole recording.
#' Ranges from -100 (perfectly inverted day-to-day) to +100 (perfectly
#' regular); 0 corresponds to chance-level agreement.
#'
#' @details
#' `algo = "vallim"` (default) ports **Fix 30** of the Python reference
#' pipeline (`SRI_vallim`): sleep/wake comes directly from the epoch-level
#' `state` column already produced by zeitR's own pipelines
#' ([run_pipeline()] / [run_pipeline_native()]) -- the same classification
#' `compute_sleep_metrics()` and `compute_cpd_metrics()` already treat as
#' ground truth for this recording. Off-wrist handling mirrors Fix 30
#' exactly: off-wrist epochs (`state == 4`) are treated as missing. Gaps of
#' `max_gap_min` minutes or less are interpolated (forward-filled from the
#' last valid epoch, or back-filled from the next valid epoch when the gap
#' starts at the very beginning of the recording); longer gaps are left as
#' missing and excluded from the day-to-day comparison entirely, rather than
#' being counted as a non-match. The aggregation is a single flat average
#' over every valid 24h-apart epoch pair, pooled across the whole recording:
#' \deqn{SRI = -100 + 200 \times \frac{1}{M}\sum_{t} \Psi(t, t + 24h)}
#' where \eqn{\Psi(t, t+24h) = 1} if the sleep/wake state at epoch \eqn{t}
#' matches the state 24 h later, \eqn{0} otherwise, and \eqn{M} is the number
#' of epoch pairs where both epochs have a valid (non-missing) state. This
#' matches Julia's own `compute_sri_vallim()` replica (also a flat pooled
#' ratio, not a per-time-of-day average -- see below).
#'
#' `algo = "sadeh"` instead derives sleep/wake from the raw `activity`
#' signal via the Sadeh et al. (1994) algorithm, matching pyActigraphy's
#' actual `Sadeh()`/`SleepRegularityIndex()` (Cell 16 of
#' `vs_condor_py_pipeline_fix30_jrsv.ipynb`) exactly -- both sourced from
#' pyActigraphy's real code (`pyActigraphy/sleep/scoring_base.py`,
#' `_sadeh()`; `pyActigraphy/sleep/scoring/sri.py`, `sri()`), not
#' reconstructed from documentation. Two precise details that differ from
#' the `"vallim"` path above:
#' \itemize{
#'   \item **Aggregation is a two-step average, not a flat pooled one.**
#'     pyActigraphy's `sri()` first groups epochs by time-of-day (hour,
#'     minute, second) across all days, averages day-to-day stability
#'     *within* each time-of-day slot, and only then averages *across*
#'     slots. This is mathematically different from a flat pooled average
#'     whenever slots have unequal numbers of valid day-pairs (e.g. a
#'     partial first/last day) -- both `"vallim"` and `"sadeh"` are ported
#'     faithfully to their own respective real source, not reconciled to
#'     use the same aggregation.
#'   \item **No off-wrist handling.** Sadeh scores whatever raw activity it
#'     is given; there is no off-wrist masking/interpolation step in
#'     pyActigraphy's own code for this path. `max_gap_min` has no effect
#'     when `algo = "sadeh"`.
#' }
#' Sadeh's own edge behaviour is also reproduced exactly: `mean_W5`/`NAT`
#' need a full centered 11-epoch window (5 before, self, 5 after) and are
#' `NA` for the first/last 5 epochs; `sd_Last6` needs a full trailing
#' 6-epoch window (self + 5 before) and is `NA` for the first 5 epochs;
#' `logAct` uses the *following* epoch's activity (`shift(-1)`) and is `NA`
#' for the last epoch. Where the resulting `PS` score is `NA`, pandas'
#' `NaN > threshold` evaluates to `False` (not propagated as missing), so
#' that epoch is scored wake (`0`), not excluded -- reproduced here
#' explicitly, since R's `NA > threshold` gives `NA`, not `FALSE`.
#'
#' @param x A `zeitr_recording`/`zeitr_result`, or a data frame / tibble with
#'   at least `datetime` and `state` columns (`algo = "vallim"`) or
#'   `datetime` and `activity` columns (`algo = "sadeh"`). For `state`:
#'   `state == 1` or `state == 7` is treated as sleep, `state == 4` as
#'   off-wrist (missing), and any other value as wake -- matching the coding
#'   already used across zeitR's pipeline output ([run_pipeline()],
#'   [run_pipeline_native()], [export_hypnogram()]).
#' @param epoch_s `numeric(1)`. Epoch duration in seconds. If `NULL`
#'   (default), estimated automatically from the median inter-epoch interval.
#' @param max_gap_min `numeric(1)`. Off-wrist gaps of this many minutes or
#'   less are interpolated rather than excluded. Default `30`, matching
#'   Fix 30 and Fix 14's 30-minute threshold elsewhere in the pipeline. Only
#'   used when `algo = "vallim"`.
#' @param algo `character(1)`. Which sleep/wake source and SRI aggregation
#'   to use: `"vallim"` (default, uses the pipeline's own `state` column) or
#'   `"sadeh"` (scores raw `activity` via the Sadeh algorithm). See Details
#'   for the precise differences between the two beyond just the scoring
#'   algorithm.
#'
#' @return A tibble with columns `participant_id`, `sri`, `n_pairs` (number
#'   of valid 24h-apart epoch comparisons used; `NA` for `algo = "sadeh"`,
#'   whose aggregation isn't a single pooled pair count), and `n_epochs`.
#'   `sri` is `NA` if the recording is shorter than 24 h, or (for
#'   `algo = "vallim"`) if no valid pairs remain after off-wrist exclusion.
#'
#' @references
#' Phillips, A. J. K., Clerx, W. M., O'Brien, C. S., Sano, A., Barger, L. K.,
#' Picard, R. W., Lockley, S. W., Klerman, E. B., & Czeisler, C. A. (2017).
#' Irregular sleep/wake patterns are associated with poorer academic
#' performance and delayed circadian and sleep/wake timing. *Scientific
#' Reports*, 7, 3216. \doi{10.1038/s41598-017-03171-4}
#'
#' Sadeh, A., Sharkey, M., & Carskadon, M. A. (1994). Activity-Based
#' Sleep-Wake Identification: An Empirical Test of Methodological Issues.
#' *Sleep*, 17(3), 201-207. \doi{10.1093/sleep/17.3.201}
#'
#' @seealso [compute_npcra()], [compute_sleep_metrics()],
#'   [run_pipeline_native()]
#'
#' @export
#'
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' result <- run_pipeline_native("recordings/P001.txt", tz = "America/Sao_Paulo")
#' compute_sri(result)
#' compute_sri(result, algo = "sadeh")
#' }
compute_sri <- function(x, epoch_s = NULL, max_gap_min = 30, algo = "vallim") {
  algo <- match.arg(algo, c("vallim", "sadeh"))

  # ── Extract epochs tibble and participant_id ──────────────────────────────
  if (inherits(x, "zeitr_result")) {
    epochs         <- x$data
    participant_id <- x$subject_id %||% NA_character_
  } else if (inherits(x, "zeitr_recording")) {
    epochs         <- x$epochs
    participant_id <- x$metadata$participant_id %||% NA_character_
  } else if (is.data.frame(x)) {
    epochs         <- x
    participant_id <- NA_character_
  } else {
    zeitr_abort("{.arg x} must be a {.cls zeitr_result}, {.cls zeitr_recording}, or a data frame.")
  }

  if (algo == "vallim") {
    .compute_sri_vallim(epochs, participant_id, epoch_s, max_gap_min)
  } else {
    .compute_sri_sadeh(epochs, participant_id, epoch_s)
  }
}

# ── Internal: algo = "vallim" (Fix 30, state-column-based) ───────────────────

#' @noRd
.compute_sri_vallim <- function(epochs, participant_id, epoch_s, max_gap_min) {
  required <- c("datetime", "state")
  missing  <- setdiff(required, names(epochs))
  if (length(missing) > 0L) {
    zeitr_abort("Missing required column(s): {.val {missing}}")
  }

  datetimes <- as.POSIXct(epochs$datetime)
  state     <- as.integer(epochs$state)

  ord <- order(datetimes)
  datetimes <- datetimes[ord]
  state     <- state[ord]

  n0 <- length(state)
  if (n0 < 2L) zeitr_abort("Need at least 2 epochs to compute SRI.")

  # ── Epoch duration ─────────────────────────────────────────────────────────
  if (is.null(epoch_s)) {
    diffs   <- as.numeric(diff(datetimes), units = "secs")
    epoch_s <- stats::median(diffs[diffs > 0], na.rm = TRUE)
  }

  # ── Sleep/wake/missing coding ────────────────────────────────────────────
  # state == 1 (main sleep) or 7 (nap, run_pipeline()'s vendor output) -> sleep
  # state == 4 (off-wrist)                                              -> NA
  # anything else (0 = wake, or unrecognised)                           -> wake
  sleep_bin <- ifelse(state == 4L, NA_real_,
                      ifelse(state == 1L | state == 7L, 1, 0))

  # ── Interpolate off-wrist gaps <= max_gap_min; leave longer gaps as NA ──
  max_gap_epochs <- round(max_gap_min * 60 / epoch_s)
  sleep_bin      <- .interpolate_short_gaps(sleep_bin, max_gap_epochs)

  n <- length(sleep_bin)

  # ── 24h-apart comparison ─────────────────────────────────────────────────
  lag_epochs <- round(24 * 3600 / epoch_s)
  if (lag_epochs >= n) {
    zeitr_warn("Recording spans less than 24 h; {.fn compute_sri} returns NA.")
    return(tibble::tibble(
      participant_id = participant_id,
      sri             = NA_real_,
      n_pairs         = 0L,
      n_epochs        = n
    ))
  }

  a     <- sleep_bin[1:(n - lag_epochs)]
  b     <- sleep_bin[(1L + lag_epochs):n]
  valid <- !is.na(a) & !is.na(b)
  n_pairs <- sum(valid)

  sri <- if (n_pairs == 0L) {
    NA_real_
  } else {
    -100 + 200 * (sum(a[valid] == b[valid]) / n_pairs)
  }

  tibble::tibble(
    participant_id = participant_id,
    sri             = round(sri, 4),
    n_pairs         = n_pairs,
    n_epochs        = n
  )
}

# ── Internal: algo = "sadeh" (raw-activity-based, pyActigraphy-faithful) ─────

#' @noRd
.compute_sri_sadeh <- function(epochs, participant_id, epoch_s) {
  required <- c("datetime", "activity")
  missing  <- setdiff(required, names(epochs))
  if (length(missing) > 0L) {
    zeitr_abort("Missing required column(s): {.val {missing}}")
  }

  datetimes <- as.POSIXct(epochs$datetime)
  activity  <- as.double(epochs$activity)

  ord <- order(datetimes)
  datetimes <- datetimes[ord]
  activity  <- activity[ord]

  n <- length(activity)
  if (n < 2L) zeitr_abort("Need at least 2 epochs to compute SRI.")

  tz <- attr(datetimes, "tzone") %||% "UTC"

  if (is.null(epoch_s)) {
    diffs   <- as.numeric(diff(datetimes), units = "secs")
    epoch_s <- stats::median(diffs[diffs > 0], na.rm = TRUE)
  }

  lag_epochs <- round(24 * 3600 / epoch_s)
  if (lag_epochs >= n) {
    zeitr_warn("Recording spans less than 24 h; {.fn compute_sri} returns NA.")
    return(tibble::tibble(
      participant_id = participant_id,
      sri             = NA_real_,
      n_pairs         = NA_integer_,
      n_epochs        = n
    ))
  }

  scoring <- .sadeh_score(activity)
  sri_val <- .sri_pyactigraphy(datetimes, scoring, tz = tz, threshold = NULL)

  tibble::tibble(
    participant_id = participant_id,
    sri             = round(sri_val, 4),
    n_pairs         = NA_integer_,   # pyActigraphy's two-step average has no single pooled pair count
    n_epochs        = n
  )
}

#' Sadeh algorithm for sleep/wake scoring
#'
#' Ports pyActigraphy's `_sadeh()` exactly (`pyActigraphy/sleep/
#' scoring_base.py`). `PS = offset + weights . [mean_W5, NAT, sd_Last6,
#' logAct]`; `PS > threshold` -> wake (1), else sleep (0) -- matching
#' pyActigraphy's own sleep=0/wake=1 convention for scoring functions
#' (opposite of zeitR's own `state` convention; the caller converts).
#'
#' `mean_W5`: centered 11-epoch rolling mean (5 before, self, 5 after);
#' `NA` at both ends where the full window isn't available (pandas
#' `.rolling(11, center=True)` default `min_periods` requires the full
#' window).
#' `NAT`: count, in the same centered 11-epoch window, of epochs with
#' activity strictly between 50 and 100.
#' `sd_Last6`: TRAILING (not centered) 6-epoch rolling SD (self + 5
#' before); `NA` for the first 5 epochs.
#' `logAct`: `log(1 + activity)` of the FOLLOWING epoch (`shift(-1)` in
#' pandas); `NA` for the last epoch.
#'
#' Where `PS` is `NA` (edge epochs), pandas' `NaN > threshold` evaluates to
#' `False`, scoring that epoch as sleep (`0`) via `.astype(int)` -- NOT
#' propagated as missing. R's `NA > threshold` gives `NA`, not `FALSE`, so
#' this is handled explicitly below to match.
#'
#' @param activity numeric vector of raw activity counts, one per epoch.
#' @param offset `numeric(1)`. Default `7.601` (pyActigraphy's default).
#' @param weights `numeric(4)`. Weights for `mean_W5`, `NAT`, `sd_Last6`,
#'   `logAct` respectively. Default `c(-0.065, -1.08, -0.056, -0.703)`.
#' @param threshold `numeric(1)`. Default `0.0`.
#' @return Integer vector, same length as `activity`: `0` = sleep,
#'   `1` = wake.
#' @noRd
.sadeh_score <- function(activity, offset = 7.601,
                         weights = c(-0.065, -1.08, -0.056, -0.703),
                         threshold = 0.0) {
  mean_W5  <- .rolling_centered_mean(activity, halfwin = 5L)
  NAT      <- .rolling_centered_count_between(activity, halfwin = 5L, lo = 50, hi = 100)
  sd_Last6 <- .rolling_trailing_sd(activity, win = 6L)
  logAct   <- c(log1p(activity[-1]), NA_real_)

  PS <- offset + weights[1L] * mean_W5 + weights[2L] * NAT +
        weights[3L] * sd_Last6 + weights[4L] * logAct

  # pandas: NaN > threshold is False -> astype(int) -> 0 (sleep). R's
  # NA > threshold gives NA, not FALSE -- handled explicitly to match.
  ifelse(is.na(PS), 0L, as.integer(PS > threshold))
}

#' Centered rolling mean, NA at the edges where the full window is unavailable
#'
#' Matches pandas' `.rolling(2*halfwin+1, center=True).mean()` default
#' behaviour (requires the full window; no partial-window fallback).
#' @noRd
.rolling_centered_mean <- function(x, halfwin) {
  n   <- length(x)
  win <- 2L * halfwin + 1L
  out <- rep(NA_real_, n)
  if (win > n) return(out)
  cs  <- c(0, cumsum(x))
  idx <- (halfwin + 1L):(n - halfwin)
  out[idx] <- (cs[idx + halfwin + 1L] - cs[idx - halfwin]) / win
  out
}

#' Centered rolling count of values strictly between lo and hi
#' @noRd
.rolling_centered_count_between <- function(x, halfwin, lo, hi) {
  ind <- as.numeric(x > lo & x < hi)
  n   <- length(ind)
  win <- 2L * halfwin + 1L
  out <- rep(NA_real_, n)
  if (win > n) return(out)
  cs  <- c(0, cumsum(ind))
  idx <- (halfwin + 1L):(n - halfwin)
  out[idx] <- cs[idx + halfwin + 1L] - cs[idx - halfwin]
  out
}

#' Trailing (right-aligned) rolling SD, NA where the full window is unavailable
#'
#' Matches pandas' `.rolling(win).std()` (not centered; window covers the
#' current epoch and the win-1 preceding it). Uses R's `sd()` (ddof = 1,
#' matching pandas' default).
#' @noRd
.rolling_trailing_sd <- function(x, win) {
  n   <- length(x)
  out <- rep(NA_real_, n)
  if (win > n) return(out)
  for (i in win:n) {
    out[i] <- stats::sd(x[(i - win + 1L):i])
  }
  out
}

# ── Internal: pyActigraphy-faithful SRI aggregation (algo != "vallim") ──────

#' Probability that consecutive values in a time-of-day slot are stable
#'
#' Ports pyActigraphy's `prob_stability()` (`pyActigraphy/sleep/scoring/
#' sri.py`) exactly: `mean(1 - abs(diff(data)))`, where `data` is
#' optionally re-binarized against `threshold` first (`data > threshold`).
#' @noRd
.sri_prob_stability <- function(ts, threshold = NULL) {
  data <- if (!is.null(threshold)) as.numeric(ts > threshold) else ts
  mean(1 - abs(diff(data)), na.rm = TRUE)
}

#' pyActigraphy-faithful SRI: two-step average (per time-of-day slot, then
#' across slots), NOT a flat pooled average over all 24h-apart pairs
#'
#' Ports `sri_profile()` + `sri()` (`pyActigraphy/sleep/scoring/sri.py`)
#' exactly: groups `scoring` by (hour, minute, second) -- i.e. by clock
#' time-of-day, across all days in the recording -- computes
#' `.sri_prob_stability()` within each time-of-day slot (across
#' consecutive days at that slot, which are exactly 24h apart), then
#' averages those per-slot values and rescales to [-100, 100]. This is
#' NOT mathematically equivalent to a flat pooled average over all pairs
#' unless every slot has the same number of valid day-pairs (see
#' `compute_sri()`'s Details).
#'
#' `split()` preserves each slot's original (chronological) order, matching
#' pandas' groupby -- required since `.sri_prob_stability()`'s `diff()`
#' depends on day-to-day order within the slot.
#'
#' @param datetimes POSIXct vector, native epoch resolution.
#' @param scoring integer/numeric vector (0/1 sleep-wake scores), same
#'   length as `datetimes`.
#' @param tz Timezone string.
#' @param threshold Passed through to `.sri_prob_stability()`; `NULL`
#'   (default) for an already-binary `scoring` input.
#' @noRd
.sri_pyactigraphy <- function(datetimes, scoring, tz, threshold = NULL) {
  hms   <- format(datetimes, "%H:%M:%S", tz = tz)
  slots <- split(scoring, hms)
  slot_probs <- vapply(slots, .sri_prob_stability, numeric(1L), threshold = threshold)
  200 * mean(slot_probs, na.rm = TRUE) - 100
}

# ── Internal helpers ───────────────────────────────────────────────────────────

#' Forward-fill (falling back to backward-fill) short runs of NA
#'
#' Runs of `NA` up to `max_gap_epochs` long are filled from the last valid
#' value before the run (matching Fix 30's "ffill/bfill" off-wrist
#' interpolation); if the run starts at the very beginning of the vector
#' (no prior valid value), it is back-filled from the next valid value
#' instead. Runs longer than `max_gap_epochs` are left untouched.
#'
#' @param x numeric vector, possibly containing `NA`
#' @param max_gap_epochs `integer(1)`. Maximum gap length (in epochs) to
#'   interpolate.
#' @noRd
.interpolate_short_gaps <- function(x, max_gap_epochs) {
  n <- length(x)
  is_na <- is.na(x)
  if (!any(is_na) || max_gap_epochs <= 0L) return(x)

  r      <- rle(is_na)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L

  for (k in seq_along(r$values)) {
    if (!r$values[k]) next            # not a gap
    if (r$lengths[k] > max_gap_epochs) next   # too long to interpolate

    s <- starts[k]; e <- ends[k]
    prev_val <- if (s > 1L) x[s - 1L] else NA_real_
    next_val <- if (e < n) x[e + 1L] else NA_real_

    if (!is.na(prev_val)) {
      x[s:e] <- prev_val
    } else if (!is.na(next_val)) {
      x[s:e] <- next_val
    }
    # else: gap spans the entire vector -- nothing to fill from, left as NA
  }

  x
}
