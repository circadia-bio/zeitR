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
#' `algo = "ck"` derives sleep/wake via pyActigraphy's native Cole-Kripke
#' implementation (`.CK()`, default `settings = "30sec_max_non_overlap"`)
#' -- a DIFFERENT weight set from the Condor-native `ColeKripke` class
#' used elsewhere in zeitR's pipeline (`R/cole_kripke.R`); the two share an
#' algorithm family name but are otherwise unrelated. Shares `"sadeh"`'s
#' two-step SRI aggregation and lack of off-wrist handling (see above).
#' Despite the reference notebook resampling to 30-second bins before
#' calling this, that round-trip is a mathematical no-op on data that was
#' only ever 1-minute resolution (verified by direct execution -- see
#' `?.ck_native_score`), so this operates directly on native-resolution
#' `activity` with no resampling needed. Also applies Webster's (1982)
#' rescoring rules afterward, matching pyActigraphy's default. Uses the
#' *opposite* threshold polarity from Sadeh (`D < threshold` = sleep here,
#' vs `PS > threshold` = sleep for Sadeh) -- both faithful to their own
#' respective source.
#'
#' `algo = "scripps"` derives sleep/wake via pyActigraphy's native Scripps
#' Clinic implementation (`Scripps()`/`_scripps()`). Structurally identical
#' to `algo = "ck"` -- same centered rolling weighted dot product, same
#' `D < threshold` = sleep polarity, same two-step SRI aggregation and
#' lack of off-wrist handling -- just a different scale/window/threshold,
#' and no Webster rescoring step (pyActigraphy's `Scripps()` doesn't call
#' `rescore()`, unlike `CK()`).
#'
#' `algo = "roenneberg"` derives sleep/wake via pyActigraphy's native
#' Roenneberg et al. algorithm (`roenneberg()`) -- by far the most
#' involved of the four: trend extraction (a 24h centered rolling mean,
#' allowing a partial window down to 12h at the recording's edges), a
#' 15%-of-trend threshold categorization, seed-finding (candidate
#' sleep-onset runs at least 30 min long), then an iterative
#' correlation-based bout-cleaning loop (each candidate onset is tested
#' against a family of triangular "sleep bout ending at position i"
#' templates over the following 12h, accepting the HIGHEST correlation
#' peak as the bout's offset -- not simply the first qualifying one; see
#' `?.find_highest_peak_idx` for a real version mismatch this caught).
#' Shares the two-step SRI aggregation and lack of off-wrist handling
#' with the other three raw-activity algorithms above. No rescoring step
#' (rescoring is specific to `CK()`).
#'
#' @param x A `zeitr_recording`/`zeitr_result`, or a data frame / tibble with
#'   at least `datetime` and `state` columns (`algo = "vallim"`) or
#'   `datetime` and `activity` columns (`algo = "sadeh"`, `"ck"`,
#'   `"scripps"`, or `"roenneberg"`). For `state`: `state == 1` or
#'   `state == 7` is treated as sleep, `state == 4` as off-wrist (missing),
#'   and any other value as wake -- matching the coding already used
#'   across zeitR's pipeline output ([run_pipeline()],
#'   [run_pipeline_native()], [export_hypnogram()]).
#' @param epoch_s `numeric(1)`. Epoch duration in seconds. If `NULL`
#'   (default), estimated automatically from the median inter-epoch interval.
#' @param max_gap_min `numeric(1)`. Off-wrist gaps of this many minutes or
#'   less are interpolated rather than excluded. Default `30`, matching
#'   Fix 30 and Fix 14's 30-minute threshold elsewhere in the pipeline. Only
#'   used when `algo = "vallim"`.
#' @param algo `character(1)`. Which sleep/wake source and SRI aggregation
#'   to use: `"vallim"` (default, uses the pipeline's own `state` column),
#'   `"sadeh"`, `"ck"`, `"scripps"`, or `"roenneberg"` (all four score raw
#'   `activity`). See Details for the precise differences beyond just the
#'   scoring algorithm.
#'
#' @return A tibble with columns `participant_id`, `sri`, `n_pairs` (number
#'   of valid 24h-apart epoch comparisons used; `NA` for `algo != "vallim"`,
#'   whose aggregation isn't a single pooled pair count), and
#'   `n_epochs`. `sri` is `NA` if the recording is shorter than 24 h, or
#'   (for `algo = "vallim"`) if no valid pairs remain after off-wrist
#'   exclusion.
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
  algo <- match.arg(algo, c("vallim", "sadeh", "ck", "scripps", "roenneberg"))

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
  } else if (algo == "sadeh") {
    .compute_sri_sadeh(epochs, participant_id, epoch_s)
  } else if (algo == "ck") {
    .compute_sri_ck(epochs, participant_id, epoch_s)
  } else if (algo == "scripps") {
    .compute_sri_scripps(epochs, participant_id, epoch_s)
  } else {
    .compute_sri_roenneberg(epochs, participant_id, epoch_s)
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

# ── Internal: algo = "ck" (native pyActigraphy Cole-Kripke, raw-activity-based) ──

#' @noRd
.compute_sri_ck <- function(epochs, participant_id, epoch_s) {
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

  scoring <- .ck_native_score(activity)
  sri_val <- .sri_pyactigraphy(datetimes, scoring, tz = tz, threshold = NULL)

  tibble::tibble(
    participant_id = participant_id,
    sri             = round(sri_val, 4),
    n_pairs         = NA_integer_,
    n_epochs        = n
  )
}

#' pyActigraphy-native Cole-Kripke algorithm for sleep/wake scoring
#'
#' Ports pyActigraphy's `.CK()` method with its default `settings =
#' "30sec_max_non_overlap"` (`pyActigraphy/sleep/scoring_base.py`'s
#' `_cole_kripke()`, called from `ScoringMixin.CK()`). This is a
#' DIFFERENT weight set than the Condor-native `ColeKripke` class already
#' ported elsewhere in zeitR (`R/cole_kripke.R`) -- the two are unrelated
#' beyond sharing an algorithm family name.
#'
#' The reference notebook's actual call
#' (`rawATR_ck.SleepRegularityIndex(algo='CK')`) first artificially
#' resamples already-1-minute-resolution data to 30-second bins
#' (`.resample('30s').sum()`) before `.CK()`'s own
#' `"30sec_max_non_overlap"` branch resamples it back to 60-second bins via
#' `.max()`. Verified by direct execution that this round-trip is a
#' mathematical no-op on genuinely-1-minute data: every artificially
#' created 30-second bin that doesn't align with a real 1-minute timestamp
#' becomes `0` (`.resample().sum()` on an empty group), and since activity
#' counts are never negative, `max(original_value, 0)` always equals the
#' original value. So this applies the `"30sec_max_non_overlap"` weights
#' directly to the native 1-minute `activity` values, with no resampling
#' needed at all.
#'
#' `D = scale * dot([A_{-4},...,A_0,...,A_{+4}], window)`, computed over a
#' centered 9-epoch window (`window`'s last two elements are `0`, so only
#' `A_{-4}` through `A_{+2}` actually contribute -- matching the original
#' Cole-Kripke algorithm's asymmetric window, just padded to 9 elements for
#' a symmetric centered rolling window). `D < threshold` -> sleep (`1`),
#' else wake (`0`) -- matching pyActigraphy's own docstring ("D < 1 ==
#' sleep, D >= 1 == wake"). Note this is the OPPOSITE polarity from
#' `.sadeh_score()` (`PS > threshold` = sleep there).
#'
#' Edge epochs (first/last 4, where the full 9-epoch window isn't
#' available) are `NA` before thresholding; pandas' `NaN < threshold`
#' evaluates to `False`, scoring them wake (`0`) -- reproduced explicitly,
#' same reasoning as `.sadeh_score()`'s edge handling.
#'
#' If `rescoring = TRUE` (pyActigraphy's default), Webster's rescoring
#' rules are applied afterward via `.rescore()`.
#'
#' @param activity numeric vector of raw activity counts, one per epoch.
#' @param scale `numeric(1)`. Default `0.0001`.
#' @param window `numeric(9)`. Default
#'   `c(50, 30, 14, 28, 12, 8, 50, 0, 0)`.
#' @param threshold `numeric(1)`. Default `1.0`.
#' @param rescoring `logical(1)`. Apply Webster's rescoring rules. Default
#'   `TRUE`.
#' @return Integer vector, same length as `activity`: `1` = sleep,
#'   `0` = wake.
#' @noRd
.ck_native_score <- function(activity, scale = 0.0001,
                             window = c(50, 30, 14, 28, 12, 8, 50, 0, 0),
                             threshold = 1.0, rescoring = TRUE) {
  win_size <- length(window)
  halfwin  <- (win_size - 1L) %/% 2L

  D  <- .rolling_centered_dot(activity, window, halfwin, scale)
  ck <- ifelse(is.na(D), 0L, as.integer(D < threshold))

  if (isTRUE(rescoring)) {
    mask <- .rescore(ck, sleep_score = 1L)
    ck   <- as.integer(ck * mask)
  }
  ck
}

#' Centered rolling weighted dot product, NA where the full window is
#' unavailable
#'
#' For position `i`, `scale * sum(x[(i-halfwin):(i+halfwin)] * window)`,
#' where `window` is applied positionally in temporal order (`window[1]`
#' pairs with `x[i-halfwin]`, `window[length(window)]` pairs with
#' `x[i+halfwin]`). Matches pandas' `.rolling(win_size, center=True)
#' .apply(_window_convolution, ...)` default behaviour (requires the full
#' window; no partial-window fallback).
#' @noRd
.rolling_centered_dot <- function(x, window, halfwin, scale) {
  n   <- length(x)
  win <- 2L * halfwin + 1L
  out <- rep(NA_real_, n)
  if (win > n) return(out)
  for (i in (halfwin + 1L):(n - halfwin)) {
    out[i] <- scale * sum(x[(i - halfwin):(i + halfwin)] * window)
  }
  out
}

# ── Internal: algo = "scripps" (native pyActigraphy Scripps, raw-activity-based) ──

#' @noRd
.compute_sri_scripps <- function(epochs, participant_id, epoch_s) {
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

  scoring <- .scripps_score(activity)
  sri_val <- .sri_pyactigraphy(datetimes, scoring, tz = tz, threshold = NULL)

  tibble::tibble(
    participant_id = participant_id,
    sri             = round(sri_val, 4),
    n_pairs         = NA_integer_,
    n_epochs        = n
  )
}

#' pyActigraphy-native Scripps Clinic algorithm for sleep/wake scoring
#'
#' Ports pyActigraphy's `Scripps()`/`_scripps()` exactly
#' (`pyActigraphy/sleep/scoring_base.py`). Structurally identical to
#' `.ck_native_score()` (same centered rolling weighted dot product, same
#' `D < threshold` = sleep polarity) -- just a different scale/window/
#' threshold, and Scripps has no rescoring step in pyActigraphy's own code
#' (unlike `CK()`, which applies Webster's rules by default).
#'
#' `window` has 21 elements (10 before, self, 10 after) but the last 10
#' are `0`, so only `A_{-10}` through `A_{+2}` actually contribute --
#' matching the original Scripps algorithm's asymmetric window, padded to
#' 21 elements for a symmetric centered rolling window.
#'
#' Edge epochs (first/last 10, where the full 21-epoch window isn't
#' available) are `NA` before thresholding; pandas' `NaN < threshold`
#' evaluates to `False`, scoring them wake (`0`) -- same reasoning as
#' `.ck_native_score()`'s edge handling.
#'
#' @param activity numeric vector of raw activity counts, one per epoch.
#' @param scale `numeric(1)`. Default `0.204`.
#' @param window `numeric(21)`. Default matches pyActigraphy's published
#'   Scripps weights.
#' @param threshold `numeric(1)`. Default `1.0`.
#' @return Integer vector, same length as `activity`: `1` = sleep,
#'   `0` = wake.
#' @noRd
.scripps_score <- function(activity, scale = 0.204,
                           window = c(0.0064, 0.0074, 0.0112, 0.0112, 0.0118,
                                      0.0118, 0.0128, 0.0188, 0.0280, 0.0664,
                                      0.0300, 0.0112, 0.0100, 0.0000, 0.0000,
                                      0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                      0.0000),
                           threshold = 1.0) {
  win_size <- length(window)
  halfwin  <- (win_size - 1L) %/% 2L

  D <- .rolling_centered_dot(activity, window, halfwin, scale)
  ifelse(is.na(D), 0L, as.integer(D < threshold))
}

# ── Internal: algo = "roenneberg" (native pyActigraphy Roenneberg, raw-activity-based) ──

#' @noRd
.compute_sri_roenneberg <- function(epochs, participant_id, epoch_s) {
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

  # Convert pyActigraphy's default time-string parameters to epoch counts
  # at this recording's native resolution, matching
  # `int(pd.Timedelta(period) / data.index.freq)` exactly.
  epoch_min           <- epoch_s / 60
  trend_period_min     <- round(24 * 60 / epoch_min)   # '24h'
  min_trend_period_min <- round(12 * 60 / epoch_min)   # '12h'
  min_seed_period_min  <- round(30 / epoch_min)         # '30Min'
  max_test_period_min  <- round(12 * 60 / epoch_min)   # '12h'
  r_consec_below_min   <- round(30 / epoch_min)         # '30Min'

  scoring <- .roenneberg_score(
    activity,
    trend_period_min     = trend_period_min,
    min_trend_period_min = min_trend_period_min,
    threshold             = 0.15,
    min_seed_period_min  = min_seed_period_min,
    max_test_period_min  = max_test_period_min,
    r_consec_below_min   = r_consec_below_min
  )

  # Cell 16 calls this with bin_threshold=False -- in Python, `False is not
  # None` is True, triggering re-binarization via `ts > threshold`(=False,
  # i.e. 0). A no-op on Roenneberg's already-{0,1} output (`x > 0` on a
  # {0,1} vector reproduces it exactly) -- so passing NULL here (no
  # re-binarization at all) gives the identical numeric result.
  sri_val <- .sri_pyactigraphy(datetimes, scoring, tz = tz, threshold = NULL)

  tibble::tibble(
    participant_id = participant_id,
    sri             = round(sri_val, 4),
    n_pairs         = NA_integer_,
    n_epochs        = n
  )
}

#' pyActigraphy-native Roenneberg algorithm for sleep/wake scoring
#'
#' Ports pyActigraphy's `roenneberg()` exactly (`pyActigraphy/sleep/
#' scoring/roenneberg.py`) -- by far the most involved of the four SRI
#' scoring algorithms ported here: trend extraction, threshold
#' categorization, seed-finding, then iterative correlation-based bout
#' cleaning. All parameters (`trend_period_min`, etc.) are pre-converted to
#' epoch counts at the recording's native resolution by the caller,
#' matching Python's `int(pd.Timedelta(period) / data.index.freq)`
#' exactly. Operates entirely on integer epoch POSITIONS rather than
#' datetime-index lookups (`.loc`/`.iloc`/`.get_loc()` in the Python
#' source) -- equivalent given a regularly-spaced recording, which the
#' whole algorithm already assumes via `data.index.freq`.
#'
#' Sub-steps, each individually verified against a real Python execution
#' before assembly (see this session's transcript for the verification
#' scripts if they need to be re-run):
#' \itemize{
#'   \item `.roenneberg_trend()`: centered rolling mean with
#'     `min_periods < window size` -- i.e. a partial window at the edges
#'     still produces a value once at least `min_trend_period_min` epochs
#'     are available, not just a full-or-nothing window like Sadeh/CK/
#'     Scripps. For an EVEN window size `W` (guaranteed here: 24h/12h in
#'     minutes is always even for any whole-minute epoch), pandas' actual
#'     centered convention puts the extra element on the LEFT:
#'     `floor(W/2)` epochs before, self, `W - 1 - floor(W/2)` after --
#'     verified against real pandas output, not assumed from
#'     documentation.
#'   \item `.roenneberg_categorize()`: `sw = 1` (sleep) where
#'     `activity <= threshold * trend`, else `0`; `NA` wherever the trend
#'     itself is `NA` (edges shorter than `min_trend_period_min`).
#'   \item `.roenneberg_seeds()`: candidate sleep-onset positions -- start
#'     of every run of `sw == 1` at least `min_seed_period_min` long.
#'     Reuses `.consecutive_values()`, which needed the `NA`-safety fix
#'     above specifically for this caller (`sw` can genuinely contain
#'     `NA`, unlike Webster rescoring's always-clean input).
#'   \item `.clean_sleep_bout()`: for a candidate seed, tests correlation
#'     of the raw `sw` values against a family of triangular "sleep bout
#'     ending at position i" templates over the next
#'     `max_test_period_min` epochs, and returns the position of the
#'     HIGHEST correlation peak (see `.find_highest_peak_idx()`) as the
#'     accepted sleep offset -- or `NA` if no qualifying peak exists.
#' }
#'
#' Main loop: walks the seeds in order. Any seed already consumed by a
#' previously-accepted bout is skipped; anything scored `1` between the
#' last accepted bout's end and the next candidate seed is reset to wake
#' (`0`) -- these were `sw == 1` runs too short to be their own seed, or
#' left over from a seed that failed the correlation test. Everything
#' before the first seed, and everything after the last accepted bout's
#' offset, is also forced to wake. Any remaining `NA` (trend undefined at
#' the very edges) resolves to wake in the final output, since this
#' function's return value feeds directly into `.sri_pyactigraphy()`,
#' which expects a clean `{0, 1}` series.
#'
#' @param activity numeric vector of raw activity counts, one per epoch.
#' @param trend_period_min,min_trend_period_min,min_seed_period_min,
#'   max_test_period_min,r_consec_below_min Epoch counts (already
#'   converted from pyActigraphy's default time-string parameters at this
#'   recording's native resolution).
#' @param threshold `numeric(1)`. Fraction of the trend used as the
#'   sleep/wake threshold. Default `0.15`.
#' @return Integer vector, same length as `activity`: `1` = sleep,
#'   `0` = wake.
#' @noRd
.roenneberg_score <- function(activity, trend_period_min, min_trend_period_min,
                              threshold, min_seed_period_min,
                              max_test_period_min, r_consec_below_min) {
  n <- length(activity)

  trend <- .roenneberg_trend(activity, win_size = trend_period_min,
                             min_win_size = min_trend_period_min)
  sw    <- .roenneberg_categorize(activity, trend, threshold)

  seeds <- .roenneberg_seeds(sw, win_size_seed = min_seed_period_min)
  if (length(seeds) == 0L) {
    return(rep(0L, n))   # no candidate sleep bout at all -- everything wake
  }

  # Score all potential sleep epochs (1) before the first seed as wake (0).
  if (seeds[1L] > 1L) {
    idx <- 1L:(seeds[1L] - 1L)
    sw[idx][!is.na(sw[idx]) & sw[idx] == 1L] <- 0L
  }

  sot    <- list()   # list of c(onset, offset) integer position pairs
  n_succ <- r_consec_below_min + 1L

  for (seed in seeds) {
    if (length(sot) > 0L && seed < sot[[length(sot)]][2L]) next

    if (length(sot) > 0L) {
      lo <- sot[[length(sot)]][2L] + 1L
      hi <- seed - 1L
      if (hi >= lo) {
        idx <- lo:hi
        sw[idx][!is.na(sw[idx]) & sw[idx] == 1L] <- 0L
      }
    }

    sleep_onset <- seed
    # Python: uncleaned_binary_data = sw.loc[onset : onset + max_test_period]
    # (inclusive both ends), then _test_sleep_bout() immediately truncates
    # to .iloc[:win_size] (max_test_period_min epochs) -- so the net window
    # is exactly max_test_period_min epochs starting at onset, clipped to
    # the recording's end if shorter.
    window_end <- min(n, sleep_onset + max_test_period_min - 1L)
    sw_window  <- sw[sleep_onset:window_end]

    offset_rel <- .clean_sleep_bout(sw_window, win_size_test = max_test_period_min,
                                    n_succ = n_succ)

    if (!is.na(offset_rel)) {
      sleep_offset <- sleep_onset + offset_rel - 1L
      sw[sleep_onset:sleep_offset] <- 1L
      sot[[length(sot) + 1L]] <- c(sleep_onset, sleep_offset)
    }
  }

  # Score all potential sleep epochs (1) after the last accepted bout's
  # offset as wake (0).
  if (length(sot) > 0L) {
    last_offset <- sot[[length(sot)]][2L]
    if (last_offset < n) {
      idx <- (last_offset + 1L):n
      sw[idx][!is.na(sw[idx]) & sw[idx] == 1L] <- 0L
    }
  }

  sw[is.na(sw)] <- 0L
  as.integer(sw)
}

#' Centered rolling mean with `min_periods < window size` (partial windows
#' allowed at the edges)
#'
#' Ports pyActigraphy's `_extract_trend()` exactly. For an EVEN `win_size`
#' `W`, the theoretical (uncapped) window for position `i` is
#' `[i - floor(W/2), i + W - 1 - floor(W/2)]` -- confirmed against real
#' pandas `.rolling(W, center=True).mean()` output, not assumed (pandas
#' puts the extra element on the LEFT for even windows, not documented
#' explicitly anywhere obvious). At the edges, this window is clipped to
#' `[1, n]`; if the resulting count of available epochs is at least
#' `min_win_size`, the mean is computed over just those -- otherwise `NA`.
#' @noRd
.roenneberg_trend <- function(x, win_size, min_win_size) {
  n      <- length(x)
  before <- win_size %/% 2L
  after  <- win_size - 1L - before
  cs     <- c(0, cumsum(x))
  out    <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    lo  <- max(1L, i - before)
    hi  <- min(n, i + after)
    cnt <- hi - lo + 1L
    if (cnt >= min_win_size) {
      out[i] <- (cs[hi + 1L] - cs[lo]) / cnt
    }
  }
  out
}

#' Roenneberg's sleep/wake threshold categorization
#'
#' Ports pyActigraphy's `_sleep_wake_categorization()` exactly: `1`
#' (sleep) where `x <= threshold * trend`, `0` otherwise, `NA` wherever
#' `trend` itself is `NA`.
#' @noRd
.roenneberg_categorize <- function(x, trend, threshold) {
  sw <- ifelse(x <= threshold * trend, 1L, 0L)
  sw[is.na(trend)] <- NA_integer_
  sw
}

#' Roenneberg's sleep-bout seed positions
#'
#' Ports pyActigraphy's `_find_sleep_bout_seeds()` exactly: the START
#' position of every run of `sw == 1` at least `win_size_seed` epochs
#' long.
#' @noRd
.roenneberg_seeds <- function(sw, win_size_seed) {
  runs <- .consecutive_values(sw, target = 1L, min_length = win_size_seed)
  if (nrow(runs) == 0L) return(integer(0))
  unname(runs[, "start"])
}

#' Roenneberg's correlation-based sleep-bout offset test
#'
#' Ports pyActigraphy's `_test_sleep_bout()` + `_clean_sleep_bout()`
#' exactly: correlates `sw_window` against a family of `m` triangular
#' "sleep bout ending at position i" templates (`m = min(win_size_test,
#' length(sw_window))`; template row `i` is `1` for positions `1..i`, `0`
#' after), then finds the position of the HIGHEST correlation peak (see
#' `.find_highest_peak_idx()`).
#'
#' @param sw_window Binary (or `NA`-containing) vector, the candidate bout
#'   window starting at the seed position.
#' @param win_size_test `integer(1)`. Maximum template/window size
#'   (`max_test_period_min`).
#' @param n_succ `integer(1)`. Number of subsequent values a peak must
#'   exceed (`r_consec_below_min + 1`).
#' @return `integer(1)` position (1-indexed, relative to the start of
#'   `sw_window`) of the accepted sleep offset, or `NA` if no peak is
#'   found.
#' @noRd
.clean_sleep_bout <- function(sw_window, win_size_test, n_succ) {
  m <- min(win_size_test, length(sw_window))
  test_data   <- sw_window[seq_len(m)]
  test_series <- outer(seq_len(m), seq_len(m), function(i, j) as.numeric(j <= i))
  corr <- .correlation_series(test_data, test_series)
  .find_highest_peak_idx(corr, n_succ = n_succ)
}

#' Pearson correlation, clipped to `[-1, 1]`
#'
#' Ports pyActigraphy's `pearsonr()` exactly (`pyActigraphy/sleep/
#' scoring/utils.py`, numba-jitted there; plain R here).
#' @noRd
.pearsonr <- function(x, y) {
  xm <- x - mean(x)
  ym <- y - mean(y)
  normxm <- sqrt(sum(xm^2))
  normym <- sqrt(sum(ym^2))
  r <- sum((xm / normxm) * (ym / normym))
  max(min(r, 1.0), -1.0)
}

#' Correlation between `x` and each row of `Y`
#'
#' Ports pyActigraphy's `correlation_series()` exactly.
#' @noRd
.correlation_series <- function(x, Y) {
  apply(Y, 1L, function(row) .pearsonr(x, row))
}

#' `TRUE` if `x[1]` is strictly greater than every other element of `x`
#'
#' Ports pyActigraphy's `is_a_peak()` exactly, including its numpy-NaN
#' comparison semantics: `NaN > y` is `False` in numpy (never `NaN`), so a
#' window containing a `NaN` correlation value (possible from `.pearsonr()`
#' on a zero-variance window) simply fails the peak test rather than
#' propagating a missing value. R's own `NA > y` gives `NA`, not `FALSE`,
#' which would otherwise crash the `if()` in `.find_highest_peak_idx()`'s
#' calling loop -- handled explicitly here to match.
#' @noRd
.is_a_peak <- function(x) {
  cmp <- x[1L] > x[-1L]
  cmp[is.na(cmp)] <- FALSE
  all(cmp)
}

#' Position of the highest rolling-window peak
#'
#' Ports pyActigraphy's `find_highest_peak_idx()` exactly (from the
#' `artvalencio/pyActigraphy` fork -- confirmed via the Docker image's
#' `jupyter/Dockerfile`, which installs `pip install git+https://
#' github.com/artvalencio/pyActigraphy` specifically, NOT the official
#' `ghammad/pyActigraphy` repo. An earlier version of this port matched
#' the official repo's `find_first_peak_idx()` instead -- a genuinely
#' different algorithm (stops at the FIRST qualifying local peak, not the
#' best one), confirmed wrong by a real end-to-end comparison against 4
#' real participant recordings: `find_first_peak_idx()` fragmented long
#' sleep bouts into several short ones, each stopping at a locally-good-
#' enough peak, while the real reference finds one long, best-fit bout.
#'
#' Two precise differences from `find_first_peak_idx()`, both confirmed
#' against the actual fork source, not assumed from the name alone:
#' \itemize{
#'   \item The internal window size is `n_succ + 1`, not `n_succ` --
#'     `find_first_peak_idx()` uses `rolling_window(x, n_succ)` directly.
#'   \item Among all positions satisfying `.is_a_peak()`, picks the one
#'     with the maximum value -- not simply the first qualifying position.
#'     Ties are broken by re-searching the WHOLE array (not just the
#'     peak candidates) for the first occurrence of that maximum value,
#'     matching the source's own `np.where(x == np.max(x[peak_candidate_
#'     idx]))[0][0]` exactly (irrelevant in practice for continuous
#'     correlation values, but reproduced faithfully regardless).
#' }
#' `NA` if no position satisfies `.is_a_peak()` at all (including when
#' `n_succ + 1 > length(x)`, which would make the peak-candidate search
#' window search space empty).
#' @noRd
.find_highest_peak_idx <- function(x, n_succ) {
  n   <- length(x)
  win <- n_succ + 1L
  if (win > n) return(NA_integer_)

  peak_candidates <- integer(0L)
  for (i in seq_len(n - win + 1L)) {
    if (.is_a_peak(x[i:(i + win - 1L)])) peak_candidates <- c(peak_candidates, i)
  }
  if (length(peak_candidates) == 0L) return(NA_integer_)

  max_val <- max(x[peak_candidates])
  which(x == max_val)[1L]
}

#' Find runs of `x == target` of at least `min_length`, as 1-indexed
#' inclusive `[start, end]` pairs
#'
#' Ports pyActigraphy's `consecutive_values()` (`pyActigraphy/sleep/
#' scoring/utils.py`) exactly, translated from Python's 0-indexed
#' half-open `[a, b)` ranges to R's 1-indexed inclusive `[start, end]`
#' (`start = a + 1`, `end = b`) -- verified against the real Python
#' function's output on several test arrays, not derived from the source
#' alone, given how easy an off-by-one error would be to introduce here.
#'
#' `NA` in `x` is treated as "not a match", matching numpy's
#' `np.equal(NaN, target)` evaluating to `False` -- R's own `NA == target`
#' would otherwise propagate `NA` into `targets`/`absdiff` and corrupt the
#' whole computation. Harmless for Webster rescoring (`scoring` is always
#' clean 0/1 there) but required for Roenneberg's categorized series,
#' which can genuinely contain `NA` (where the trend is undefined).
#'
#' @param x integer/numeric vector, may contain `NA`.
#' @param target Value to find runs of.
#' @param min_length `integer(1)`. Minimum run length to include.
#' @return A 2-column matrix (`start`, `end`), one row per qualifying run,
#'   0 rows if none found.
#' @noRd
.consecutive_values <- function(x, target, min_length) {
  n <- length(x)
  targets <- c(0L, as.integer(!is.na(x) & x == target), 0L)
  absdiff <- abs(diff(targets))
  pos <- which(absdiff == 1L)
  if (length(pos) == 0L) return(matrix(integer(0), ncol = 2, dimnames = list(NULL, c("start", "end"))))
  ranges <- matrix(pos, ncol = 2, byrow = TRUE)
  starts <- ranges[, 1L]
  ends   <- ranges[, 2L] - 1L
  keep   <- (ends - starts + 1L) >= min_length
  cbind(start = starts[keep], end = ends[keep])
}

#' Webster's rescoring rule: rescore the start of a sleep run preceded by
#' wake
#'
#' Ports pyActigraphy's `rescore_if_preceded()` (`pyActigraphy/sleep/
#' scoring/utils.py`) exactly: for every run of `sleep_score` at least
#' `n_periods` long, if the `n_previous` epochs immediately before it are
#' all wake (`0`), rescore the run's first `n_periods` epochs to wake.
#' @noRd
.rescore_if_preceded <- function(scoring, n_periods, n_previous, sleep_score = 1L) {
  n    <- length(scoring)
  mask <- rep(1L, n)
  runs <- .consecutive_values(scoring, target = sleep_score, min_length = n_periods)
  if (nrow(runs) == 0L) return(mask)
  for (r in seq_len(nrow(runs))) {
    start <- runs[r, "start"]; end <- runs[r, "end"]
    if (start <= n_previous) next   # not enough preceding epochs
    if (end > n) next               # defensive, mirrors Python's own guard
    if (all(scoring[(start - n_previous):(start - 1L)] == 0L)) {
      mask[start:(start + n_periods - 1L)] <- 0L
    }
  }
  mask
}

#' Webster's rescoring rule: rescore a short sleep gap surrounded by long
#' wake runs on both sides
#'
#' Ports pyActigraphy's `rescore_if_surrounded()` (`pyActigraphy/sleep/
#' scoring/utils.py`) exactly: for every pair of consecutive wake runs
#' (each at least `n_surround` long), if the gap between them is at most
#' `n_periods` epochs, rescore that whole gap to wake -- regardless of its
#' actual sleep/wake composition.
#' @noRd
.rescore_if_surrounded <- function(scoring, n_periods, n_surround, sleep_score = 1L) {
  n    <- length(scoring)
  mask <- rep(1L, n)
  wake_runs <- .consecutive_values(scoring, target = abs(sleep_score - 1L), min_length = n_surround)
  if (nrow(wake_runs) < 2L) return(mask)
  for (i in seq_len(nrow(wake_runs) - 1L)) {
    end_current    <- wake_runs[i, "end"]
    start_next     <- wake_runs[i + 1L, "start"]
    sleep_duration <- (start_next - 1L) - end_current
    if (sleep_duration <= n_periods) {
      gap_start <- end_current + 1L
      gap_end   <- start_next - 1L
      if (gap_end >= gap_start) mask[gap_start:gap_end] <- 0L
    }
  }
  mask
}

#' Webster's (1982) rescoring rules, combined
#'
#' Ports pyActigraphy's `rescore()` (`pyActigraphy/sleep/scoring/utils.py`)
#' exactly: the five rules below, multiplied together elementwise (`0`
#' from any rule wins).
#' @noRd
.rescore <- function(scoring, sleep_score = 1L) {
  m1 <- .rescore_if_preceded(scoring, n_periods = 1L, n_previous = 4L,  sleep_score = sleep_score)
  m2 <- .rescore_if_preceded(scoring, n_periods = 3L, n_previous = 10L, sleep_score = sleep_score)
  m3 <- .rescore_if_preceded(scoring, n_periods = 4L, n_previous = 15L, sleep_score = sleep_score)
  m4 <- .rescore_if_surrounded(scoring, n_periods = 6L,  n_surround = 10L, sleep_score = sleep_score)
  m5 <- .rescore_if_surrounded(scoring, n_periods = 10L, n_surround = 20L, sleep_score = sleep_score)
  m1 * m2 * m3 * m4 * m5
}

#' Sadeh algorithm for sleep/wake scoring
#'
#' Ports pyActigraphy's `_sadeh()` exactly (`pyActigraphy/sleep/
#' scoring_base.py`). `PS = offset + weights . [mean_W5, NAT, sd_Last6,
#' logAct]`; `PS > threshold` -> sleep (1), else wake (0) -- matching
#' pyActigraphy's own docstring convention ("PS >= 0 == sleep, PS < 0 ==
#' wake") and zeitR's own `state` convention (`1` = sleep). Note this is
#' the OPPOSITE polarity from `.ck_native_score()` below (`D < threshold`
#' = sleep there) -- both are faithful to their own respective source,
#' just opposite conventions in the original algorithms.
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
#' `False`, scoring that epoch as wake (`0`) via `.astype(int)` -- NOT
#' propagated as missing. R's `NA > threshold` gives `NA`, not `FALSE`, so
#' this is handled explicitly below to match.
#'
#' @param activity numeric vector of raw activity counts, one per epoch.
#' @param offset `numeric(1)`. Default `7.601` (pyActigraphy's default).
#' @param weights `numeric(4)`. Weights for `mean_W5`, `NAT`, `sd_Last6`,
#'   `logAct` respectively. Default `c(-0.065, -1.08, -0.056, -0.703)`.
#' @param threshold `numeric(1)`. Default `0.0`.
#' @return Integer vector, same length as `activity`: `1` = sleep,
#'   `0` = wake.
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

  # pandas: NaN > threshold is False -> astype(int) -> 0 (wake). R's
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
