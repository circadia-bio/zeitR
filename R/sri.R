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
#' Ported from **Fix 30** of the Python reference pipeline (`SRI_vallim`):
#' rather than deriving sleep/wake from a pyActigraphy scoring algorithm
#' (Sadeh, Cole-Kripke, Roenneberg, Scripps -- all of which showed
#' substantially worse agreement with manual reference scoring in Julia's
#' concordance analysis, ICC 0.19-0.67 vs 0.82 here), `compute_sri()` derives
#' sleep/wake directly from the epoch-level `state` column already produced
#' by zeitR's own pipelines ([run_pipeline()] / [run_pipeline_native()]) --
#' the same classification `compute_sleep_metrics()` and
#' `compute_cpd_metrics()` already treat as ground truth for this recording.
#'
#' Off-wrist handling mirrors Fix 30 exactly: off-wrist epochs
#' (`state == 4`) are treated as missing. Gaps of `max_gap_min` minutes or
#' less are interpolated (forward-filled from the last valid epoch, or
#' back-filled from the next valid epoch when the gap starts at the very
#' beginning of the recording); longer gaps are left as missing and excluded
#' from the day-to-day comparison entirely, rather than being counted as a
#' non-match.
#'
#' \deqn{SRI = -100 + 200 \times \frac{1}{M}\sum_{t} \Psi(t, t + 24h)}
#'
#' where \eqn{\Psi(t, t+24h) = 1} if the sleep/wake state at epoch \eqn{t}
#' matches the state 24 h later, \eqn{0} otherwise, and \eqn{M} is the number
#' of epoch pairs where both epochs have a valid (non-missing) state.
#'
#' @param x A `zeitr_recording`/`zeitr_result`, or a data frame / tibble with
#'   at least `datetime` and `state` columns. `state == 1` or `state == 7`
#'   is treated as sleep, `state == 4` as off-wrist (missing), and any other
#'   value as wake -- matching the coding already used across zeitR's
#'   pipeline output ([run_pipeline()], [run_pipeline_native()],
#'   [export_hypnogram()]).
#' @param epoch_s `numeric(1)`. Epoch duration in seconds. If `NULL`
#'   (default), estimated automatically from the median inter-epoch interval.
#' @param max_gap_min `numeric(1)`. Off-wrist gaps of this many minutes or
#'   less are interpolated rather than excluded. Default `30`, matching
#'   Fix 30 and Fix 14's 30-minute threshold elsewhere in the pipeline.
#'
#' @return A tibble with columns `participant_id`, `sri`, `n_pairs` (number
#'   of valid 24h-apart epoch comparisons used), and `n_epochs` (total
#'   epochs after off-wrist gap interpolation, before the 24h-pairing step).
#'   `sri` is `NA` if the recording is shorter than 24 h, or if no valid
#'   pairs remain after off-wrist exclusion.
#'
#' @references
#' Phillips, A. J. K., Clerx, W. M., O'Brien, C. S., Sano, A., Barger, L. K.,
#' Picard, R. W., Lockley, S. W., Klerman, E. B., & Czeisler, C. A. (2017).
#' Irregular sleep/wake patterns are associated with poorer academic
#' performance and delayed circadian and sleep/wake timing. *Scientific
#' Reports*, 7, 3216. \doi{10.1038/s41598-017-03171-4}
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
#' }
compute_sri <- function(x, epoch_s = NULL, max_gap_min = 30) {

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

  # ── Sleep/wake/missing coding ──────────────────────────────────────────────
  # state == 1 (main sleep) or 7 (nap, run_pipeline()'s vendor output) -> sleep
  # state == 4 (off-wrist)                                              -> NA
  # anything else (0 = wake, or unrecognised)                           -> wake
  sleep_bin <- ifelse(state == 4L, NA_real_,
                      ifelse(state == 1L | state == 7L, 1, 0))

  # ── Interpolate off-wrist gaps <= max_gap_min; leave longer gaps as NA ────
  max_gap_epochs <- round(max_gap_min * 60 / epoch_s)
  sleep_bin      <- .interpolate_short_gaps(sleep_bin, max_gap_epochs)

  n <- length(sleep_bin)

  # ── 24h-apart comparison ───────────────────────────────────────────────────
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
