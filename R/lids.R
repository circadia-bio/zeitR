# ── LIDS (Locomotor Inactivity During Sleep) ──────────────────────────────────
# Ports the LIDS ultradian-rhythm methodology from:
#
#   Winnebeck, Fischer, Leise & Roenneberg (2018), Current Biology 28, 49-59
#     -- original method: LIDS non-linear transform, plain cosine fit, Munich
#        Rhythmicity Index (MRI) period scan, moving-average smoothing.
#        Mirrored (where the two papers overlap) by pyActigraphy's LIDS class:
#        https://ghammad.github.io/pyActigraphy/LIDS.html
#
#   Hammad, Schoch, Engelmann, Spock, Kurth & Winnebeck (2026), SLEEP
#     -- infant extension: adds a linear slope term to the cosine (sloped
#        cosine), Gaussian smoothing (sigma = 5 min), and a 30-180 min /
#        2-min period scan tuned for shorter (~60 min) ultradian cycles.
#        https://zenodo.org/records/18199381
#
# Ported from a prototype R notebook draft by Mario Leocadio-Miguel (itself
# adapted from an older MATLAB script), informed by both papers above. One
# bug in that draft is fixed here: its bout-fusing step
# (`fused.append(...)`) was truncated mid-statement in the source notebook
# and would not have run as written.
#
# Not yet validated against pyActigraphy's LIDS output or an external
# reference dataset -- see individual function docs for what's Hammad-2026
# vs. Winnebeck-2018 behaviour, and treat results accordingly until a parity
# check is run (cf. the zeitR parity-first philosophy used everywhere else
# in this package).

# ── LIDS transform ─────────────────────────────────────────────────────────

#' Apply the LIDS (Locomotor Inactivity During Sleep) transform
#'
#' Converts an activity signal into "inactivity" via the LIDS non-linear
#' transform, then smooths it -- the first step of the LIDS ultradian-rhythm
#' pipeline (Winnebeck et al. 2018; Hammad et al. 2026). LIDS is only ever
#' computed *inside* an already-identified sleep bout; extract bouts first
#' via [compute_lids()] or [detect_lids_bouts()].
#'
#' @details
#' \deqn{\text{LIDS}_i = \frac{100}{1+x_i}}
#'
#' where \eqn{x_i} is the raw activity count at epoch \eqn{i}. A LIDS value of
#' 100 means zero movement; it falls toward 0 as movement increases.
#'
#' Two smoothing methods are available, both over a nominal 30-min window:
#' * `"gaussian"` (default) -- Gaussian kernel, standard deviation
#'   `sigma_min` (default 5 min), truncated at +/- 3 sigma (Hammad et al.
#'   2026).
#' * `"mva"` -- centered moving average (Winnebeck et al. 2018;
#'   `pyActigraphy`'s default). Uses zeitR's border-replicated internal
#'   `rolling_mean_cpp()`, which replicates the edge value rather than
#'   shrinking the window near the bout boundary (pandas'
#'   `min_periods=1` behaviour) -- a minor difference confined to the first/
#'   last ~15 min of each bout.
#'
#' Any `NA` in `activity` (e.g. a brief off-wrist gap inside an otherwise
#' valid bout) is linearly interpolated first -- [fit_lids()] cannot handle
#' missing values.
#'
#' @param activity Numeric vector of activity counts for a single sleep bout,
#'   in chronological order at a constant epoch length.
#' @param epoch_min `numeric(1)`. Epoch duration in minutes. Default `1`.
#' @param method `character(1)`. `"gaussian"` (default) or `"mva"`.
#' @param win_min `numeric(1)`. Smoothing window width in minutes (full
#'   width for `"mva"`; +/- 3 sigma width for `"gaussian"`). Default `30`.
#' @param sigma_min `numeric(1)`. Gaussian kernel standard deviation in
#'   minutes, used only when `method = "gaussian"`. Default `5`
#'   (Hammad et al. 2026); ignored for `"mva"`.
#'
#' @return Numeric vector of smoothed LIDS values, same length as `activity`.
#'
#' @references
#' Winnebeck, E. C., Fischer, D., Leise, T., & Roenneberg, T. (2018). Dynamics
#' and Ultradian Structure of Human Sleep in Real Life. *Current Biology*,
#' 28(1), 49-59.e5. \doi{10.1016/j.cub.2017.11.063}
#'
#' Hammad, G., Schoch, S. F., Engelmann, M., Spock, Z., Kurth, S., & Winnebeck,
#' E. C. (2026). Charting infant sleep cycle development using actigraphy:
#' Longitudinal evidence for ultradian cycle lengthening within the first
#' year of life. *SLEEP*.
#'
#' @seealso [fit_lids()], [detect_lids_bouts()], [compute_lids()]
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' activity <- pmax(0, 20 + 15 * sin(seq(0, 6 * pi, length.out = 360)) +
#'                     rnorm(360, sd = 5))
#' lids <- lids_transform(activity, method = "gaussian")
lids_transform <- function(activity,
                            epoch_min = 1,
                            method    = c("gaussian", "mva"),
                            win_min   = 30,
                            sigma_min = 5) {
  method   <- match.arg(method)
  activity <- as.double(activity)
  n <- length(activity)
  if (n < 2L) zeitr_abort("{.arg activity} must have at least 2 epochs.")

  if (anyNA(activity)) {
    activity <- stats::approx(seq_len(n), activity, xout = seq_len(n), rule = 2)$y
  }

  lids_raw <- 100 / (1 + activity)
  win_bins <- max(1L, round(win_min / epoch_min))

  if (method == "mva") {
    rolling_mean_cpp(lids_raw, hws = win_bins %/% 2L, replicate = TRUE)
  } else {
    .gaussian_smooth(lids_raw, sigma_bins = sigma_min / epoch_min, win_bins = win_bins)
  }
}

#' Gaussian-kernel smoothing truncated at +/- 3 sigma, border-replicated
#' @noRd
.gaussian_smooth <- function(x, sigma_bins, win_bins) {
  n   <- length(x)
  hws <- max(1L, round(win_bins / 2))
  idx <- seq(-hws, hws)
  kernel <- exp(-(idx^2) / (2 * sigma_bins^2))
  kernel <- kernel / sum(kernel)

  padded <- c(rep(x[1], hws), x, rep(x[n], hws))
  out    <- stats::filter(padded, kernel, sides = 2)
  as.numeric(out[(hws + 1L):(hws + n)])
}

# ── Sloped cosine fit ──────────────────────────────────────────────────────

#' Fit a sloped cosine to a LIDS profile, scanning candidate periods
#'
#' For a fixed period \eqn{T}, the sloped-cosine model
#' \deqn{f(t) = \alpha\cos(2\pi t/T) + \beta\sin(2\pi t/T) + b + s\,t}
#' is linear in \eqn{(\alpha,\beta,b,s)}, so it is solved by ordinary least
#' squares rather than a non-linear optimiser (Hammad et al. 2026). Candidate
#' periods are scanned over `period_range` in steps of `period_step`, and for
#' each one the fit's Munich Rhythmicity Index
#' (\eqn{\text{MRI} = 2 \times \text{amplitude} \times r}, Winnebeck et al.
#' 2018) is computed. The period with the highest MRI is returned as the
#' bout's estimated ultradian cycle length -- the same selection rule used by
#' `pyActigraphy.analysis.LIDS.lids_fit()`, generalised here to include the
#' linear slope term.
#'
#' @param lids Numeric vector of (smoothed) LIDS values for one sleep bout,
#'   as returned by [lids_transform()]. Must not contain `NA`.
#' @param epoch_min `numeric(1)`. Epoch duration in minutes. Default `1`.
#' @param period_range `numeric(2)`. Candidate period bounds in minutes.
#'   Default `c(30, 180)` (Hammad et al. 2026, tuned for infant/ultradian
#'   cycles); use `c(60, 180)` with `period_step = 5` to match Winnebeck et
#'   al. (2018)'s adult/adolescent scan.
#' @param period_step `numeric(1)`. Step size in minutes for the period scan.
#'   Default `2`.
#'
#' @return A named list:
#'   \describe{
#'     \item{`period_min`}{Estimated cycle length (minutes) at peak MRI.}
#'     \item{`amplitude`}{\eqn{\sqrt{\alpha^2+\beta^2}}.}
#'     \item{`phase_rad`}{\eqn{-\mathrm{atan2}(\beta,\alpha)}, radians; `0` =
#'       LIDS peak at bout start.}
#'     \item{`offset`}{Inactivity level at bout start (\eqn{b}).}
#'     \item{`slope_per_60min`}{Linear trend, rescaled to LIDS units per hour.}
#'     \item{`pearson_r`}{Correlation between fitted and observed LIDS.}
#'     \item{`p_value`}{Two-sided p-value for `pearson_r` ([stats::cor.test()]).}
#'     \item{`mri`}{Munich Rhythmicity Index at the selected period.}
#'   }
#'
#' @references
#' Winnebeck, E. C., Fischer, D., Leise, T., & Roenneberg, T. (2018). Dynamics
#' and Ultradian Structure of Human Sleep in Real Life. *Current Biology*,
#' 28(1), 49-59.e5. \doi{10.1016/j.cub.2017.11.063}
#'
#' Hammad, G., Schoch, S. F., Engelmann, M., Spock, Z., Kurth, S., & Winnebeck,
#' E. C. (2026). Charting infant sleep cycle development using actigraphy.
#' *SLEEP*.
#'
#' @seealso [lids_transform()], [compute_lids()]
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' t <- seq(0, 300, by = 1)
#' lids <- 85 + 15 * cos(2 * pi * t / 60) - 0.05 * t + rnorm(length(t), sd = 2)
#' fit_lids(lids)
fit_lids <- function(lids,
                      epoch_min    = 1,
                      period_range = c(30, 180),
                      period_step  = 2) {
  lids <- as.double(lids)
  n <- length(lids)
  if (n < 4L) zeitr_abort("{.arg lids} must have at least 4 epochs to fit a cosine.")
  if (anyNA(lids)) {
    zeitr_abort(
      "{.arg lids} contains missing values; interpolate first (see {.fn lids_transform})."
    )
  }

  t <- (seq_len(n) - 1L) * epoch_min
  periods <- seq(period_range[1], period_range[2], by = period_step)

  best <- NULL
  for (period in periods) {
    X <- cbind(cos(2 * pi * t / period), sin(2 * pi * t / period), 1, t)
    coefs <- tryCatch(unname(qr.solve(X, lids)), error = function(e) NULL)
    if (is.null(coefs)) next

    alpha <- coefs[1]; beta <- coefs[2]; b <- coefs[3]; s <- coefs[4]
    fitted_vals <- as.vector(X %*% coefs)

    r <- suppressWarnings(stats::cor(lids, fitted_vals))
    if (is.na(r)) r <- 0
    amplitude <- sqrt(alpha^2 + beta^2)
    mri <- 2 * amplitude * r

    if (is.null(best) || mri > best$mri) {
      best <- list(
        period_min      = period,
        amplitude       = amplitude,
        phase_rad       = -atan2(beta, alpha),
        offset          = b,
        slope_per_60min = s * 60 / epoch_min,
        pearson_r       = r,
        mri             = mri,
        fitted          = fitted_vals
      )
    }
  }

  if (is.null(best)) zeitr_abort("Cosine fit failed for all candidate periods.")

  best$p_value <- suppressWarnings(tryCatch(
    stats::cor.test(lids, best$fitted)$p.value,
    error = function(e) NA_real_
  ))
  best$fitted <- NULL
  best
}

# ── Independent bout detection (Roenneberg relative-immobility) ────────────

#' Detect nighttime sleep bouts via the Roenneberg relative-immobility method
#'
#' Standalone sleep-bout detector for raw actigraphy, independent of zeitR's
#' main Crespo-based pipeline ([detect_sleep_crespo()]) or the Vallim native
#' pipeline ([run_pipeline_native()]). Intended for LIDS analysis on activity
#' recordings that haven't been run through either of those. If you already
#' have a `zeitr_result`, pass it straight to [compute_lids()] instead
#' (`bout_source = "state"`, the default there) and skip this function.
#'
#' @details
#' Ported (with one bug fixed -- see below) from a prototype R notebook
#' translation of Mario Leocadio-Miguel's method, itself adapted from the
#' relative-immobility algorithm used in Winnebeck et al. (2018) and
#' Hammad et al. (2026):
#' \enumerate{
#'   \item Compute a `ma_window_h`-centered moving average of activity (the
#'     recording's own baseline).
#'   \item Flag epochs where activity < `relative_threshold` times that
#'     moving average as *candidate sleep*.
#'   \item **Consolidate**: brief active blips (<= `bridge_min` minutes),
#'     surrounded on both sides by candidate-sleep epochs, are relabelled as
#'     sleep.
#'   \item Keep only consolidated runs lasting >= `min_bout_min`.
#'   \item **Fuse** bouts separated by a gap <= `max_gap_min`. (This is the
#'     step truncated mid-statement in the source notebook,
#'     `fused.a...` -- reimplemented in full here as the internal
#'     `.fuse_bouts()`. Unlike the notebook, which pads fused gaps with
#'     `NaN` to preserve a regular time axis, the fused bout here simply
#'     spans `start` to `end`
#'     over the original activity values -- fine for bout *timing*, but if
#'     you need the gap epochs excluded from the LIDS fit itself, mask them
#'     to `NA` in `activity` beforehand.)
#'   \item Restrict to bouts starting within `main_window` (default
#'     18:00-08:00), excluding daytime naps / recording artefacts.
#'   \item Filter by total duration (`duration_range`) and, if
#'     `one_per_night = TRUE`, keep only the longest bout per calendar night.
#' }
#'
#' The moving-average step uses zeitR's border-replicated
#' `rolling_mean_cpp()` rather than pandas' `min_periods=1` edge behaviour
#' (which shrinks the window near the recording boundary instead of
#' replicating the edge value) -- a minor difference confined to the first/
#' last `ma_window_h`/2 hours of the recording.
#'
#' @param datetime `POSIXct` vector of epoch timestamps, regularly spaced.
#' @param activity Numeric vector of activity counts, same length as
#'   `datetime`.
#' @param relative_threshold `numeric(1)`. Default `0.15`.
#' @param ma_window_h `numeric(1)`. Moving-average window in hours. Default
#'   `24`.
#' @param bridge_min `numeric(1)`. Maximum active-blip length (minutes) to
#'   bridge during consolidation. Default `5`.
#' @param min_bout_min `numeric(1)`. Minimum consolidated-run length
#'   (minutes) to keep as a candidate bout. Default `30`.
#' @param max_gap_min `numeric(1)`. Maximum gap (minutes) between candidate
#'   bouts to fuse them into one. Default `15`.
#' @param main_window `character(2)`. `c(start, end)` clock times
#'   (`"HH:MM"`) defining the window a bout must *start* within. Default
#'   `c("18:00", "08:00")` (wraps midnight).
#' @param duration_range `numeric(2)`. Final bout duration bounds, in hours.
#'   Default `c(3, 12)` (Winnebeck et al. 2018).
#' @param one_per_night `logical(1)`. If `TRUE` (default), keep only the
#'   longest surviving bout per calendar night (the "main sleep episode").
#'
#' @return A tibble with one row per detected bout: `bout_id`, `bout_start`,
#'   `bout_end`, `duration_h`, `n_epochs`.
#'
#' @references
#' Winnebeck, E. C., Fischer, D., Leise, T., & Roenneberg, T. (2018). Dynamics
#' and Ultradian Structure of Human Sleep in Real Life. *Current Biology*,
#' 28(1), 49-59.e5. \doi{10.1016/j.cub.2017.11.063}
#'
#' @seealso [compute_lids()], [lids_transform()], [fit_lids()]
#'
#' @export
#'
#' @importFrom tibble tibble
detect_lids_bouts <- function(datetime,
                               activity,
                               relative_threshold = 0.15,
                               ma_window_h        = 24,
                               bridge_min          = 5,
                               min_bout_min        = 30,
                               max_gap_min         = 15,
                               main_window         = c("18:00", "08:00"),
                               duration_range      = c(3, 12),
                               one_per_night       = TRUE) {

  datetime <- as.POSIXct(datetime)
  activity <- as.double(activity)
  n <- length(activity)
  if (n != length(datetime)) {
    zeitr_abort("{.arg datetime} and {.arg activity} must be the same length.")
  }
  if (n < 2L) zeitr_abort("Need at least 2 epochs.")

  empty <- tibble::tibble(
    bout_id = integer(), bout_start = as.POSIXct(character()),
    bout_end = as.POSIXct(character()), duration_h = double(), n_epochs = integer()
  )

  ord      <- order(datetime)
  datetime <- datetime[ord]
  activity <- activity[ord]

  diffs     <- as.numeric(diff(datetime), units = "secs")
  epoch_s   <- stats::median(diffs[diffs > 0], na.rm = TRUE)
  epoch_min <- epoch_s / 60

  # 1. Centred moving-average baseline
  ma_bins    <- max(1L, round(ma_window_h * 60 / epoch_min))
  activity0  <- ifelse(is.na(activity), 0, activity)
  moving_avg <- rolling_mean_cpp(activity0, hws = ma_bins %/% 2L, replicate = TRUE)

  # 2. Candidate sleep mask
  candidate <- activity0 < (relative_threshold * moving_avg) & !is.na(activity)

  # 3. Consolidate brief active blips
  bridge_bins <- max(0L, round(bridge_min / epoch_min))
  candidate   <- .bridge_short_gaps(candidate, bridge_bins)

  # 4. Keep runs >= min_bout_min
  min_bins  <- max(1L, round(min_bout_min / epoch_min))
  candidate <- .drop_short_runs(candidate, min_bins)

  bouts <- .runs_to_bouts(candidate, datetime)
  if (length(bouts) == 0L) return(empty)

  # 5. Fuse bouts separated by a short gap
  bouts <- .fuse_bouts(bouts, max_gap_min * 60)

  # 6. Restrict to bouts starting within the main nighttime window
  bouts <- Filter(function(b) .starts_in_window(b$start, main_window), bouts)

  # 7. Filter by duration
  dmin <- duration_range[1] * 3600
  dmax <- duration_range[2] * 3600
  bouts <- Filter(function(b) {
    dur <- as.numeric(difftime(b$end, b$start, units = "secs"))
    dur >= dmin && dur <= dmax
  }, bouts)

  # 8. One bout per calendar night
  if (isTRUE(one_per_night) && length(bouts) > 0L) {
    bouts <- .select_main_bout_per_night(bouts)
  }

  if (length(bouts) == 0L) return(empty)

  tibble::tibble(
    bout_id    = seq_along(bouts),
    bout_start = do.call(c, lapply(bouts, `[[`, "start")),
    bout_end   = do.call(c, lapply(bouts, `[[`, "end")),
    duration_h = vapply(bouts, function(b) {
      as.numeric(difftime(b$end, b$start, units = "hours"))
    }, double(1)),
    n_epochs   = vapply(bouts, function(b) b$n_epochs, integer(1))
  )
}

# ── Internal helpers for detect_lids_bouts() ────────────────────────────────

#' @noRd
.bridge_short_gaps <- function(candidate, bridge_bins) {
  if (bridge_bins <= 0L) return(candidate)
  r      <- rle(candidate)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  out    <- candidate

  for (k in seq_along(r$values)) {
    if (r$values[k]) next                      # already candidate-sleep
    if (r$lengths[k] > bridge_bins) next        # too long to bridge
    surrounded <- k > 1L && k < length(r$values) && r$values[k - 1L] && r$values[k + 1L]
    if (surrounded) out[starts[k]:ends[k]] <- TRUE
  }
  out
}

#' @noRd
.drop_short_runs <- function(candidate, min_bins) {
  r <- rle(candidate)
  r$values[r$values & r$lengths < min_bins] <- FALSE
  inverse.rle(r)
}

#' @noRd
.runs_to_bouts <- function(candidate, datetime) {
  r      <- rle(candidate)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  keep   <- which(r$values)

  lapply(keep, function(k) {
    list(start = datetime[starts[k]], end = datetime[ends[k]], n_epochs = r$lengths[k])
  })
}

#' Merge bouts separated by a gap of at most `max_gap_s` seconds
#'
#' The R reimplementation of the notebook's truncated `fused.a...` line.
#' @noRd
.fuse_bouts <- function(bouts, max_gap_s) {
  if (length(bouts) <= 1L) return(bouts)

  fused <- list(bouts[[1L]])
  for (i in 2:length(bouts)) {
    last <- fused[[length(fused)]]
    gap  <- as.numeric(difftime(bouts[[i]]$start, last$end, units = "secs"))
    if (gap <= max_gap_s) {
      fused[[length(fused)]] <- list(
        start    = last$start,
        end      = bouts[[i]]$end,
        n_epochs = last$n_epochs + bouts[[i]]$n_epochs
      )
    } else {
      fused[[length(fused) + 1L]] <- bouts[[i]]
    }
  }
  fused
}

#' @noRd
.hm_to_h <- function(hm) {
  parts <- as.numeric(strsplit(hm, ":")[[1]])
  parts[1] + parts[2] / 60
}

#' @noRd
.starts_in_window <- function(dt, window) {
  h  <- as.numeric(format(dt, "%H")) + as.numeric(format(dt, "%M")) / 60
  w1 <- .hm_to_h(window[1])
  w2 <- .hm_to_h(window[2])
  if (w1 <= w2) h >= w1 & h <= w2 else h >= w1 | h <= w2
}

#' @noRd
.select_main_bout_per_night <- function(bouts) {
  night_of <- function(b) {
    d <- as.Date(b$start)
    h <- as.numeric(format(b$start, "%H"))
    as.character(if (h < 12) d else d + 1L)
  }
  nights <- vapply(bouts, night_of, character(1))
  durs   <- vapply(bouts, function(b) {
    as.numeric(difftime(b$end, b$start, units = "secs"))
  }, double(1))

  keep <- tapply(seq_along(bouts), nights, function(idx) idx[which.max(durs[idx])])
  bouts[unlist(keep, use.names = FALSE)]
}

# ── Full pipeline: compute_lids() ───────────────────────────────────────────

#' @noRd
.lids_na_tibble <- function(participant_id) {
  tibble::tibble(
    participant_id = participant_id, bout_id = integer(),
    bout_start = as.POSIXct(character()), bout_end = as.POSIXct(character()),
    duration_h = double(), period_min = double(), amplitude = double(),
    offset = double(), slope_per_60min = double(), phase_rad = double(),
    pearson_r = double(), p_value = double(), mri = double(),
    passes_quality_filter = logical()
  )
}

#' @noRd
.state_bouts <- function(datetime, state, duration_range) {
  if (is.null(state)) {
    zeitr_abort(
      'No {.val state} column available for {.code bout_source = "state"};
       use {.code bout_source = "roenneberg"} for raw activity input.'
    )
  }
  is_sleep <- state == 1L
  r        <- rle(is_sleep)
  ends     <- cumsum(r$lengths)
  starts   <- ends - r$lengths + 1L
  keep     <- which(r$values)

  dmin <- duration_range[1] * 3600
  dmax <- duration_range[2] * 3600

  bouts <- lapply(keep, function(k) {
    list(start = datetime[starts[k]], end = datetime[ends[k]], idx = starts[k]:ends[k])
  })
  Filter(function(b) {
    dur <- as.numeric(difftime(b$end, b$start, units = "secs"))
    dur >= dmin && dur <= dmax
  }, bouts)
}

#' Compute LIDS ultradian-rhythm parameters for every sleep bout
#'
#' The main entry point for the LIDS pipeline: extracts sleep bouts, applies
#' [lids_transform()] and [fit_lids()] to each, and returns one row per bout
#' with cosine-fit parameters and a quality-filter flag. Ports the full
#' pipeline described in Winnebeck et al. (2018) and Hammad et al. (2026).
#'
#' @details
#' # Where bouts come from (`bout_source`)
#' * `"state"` -- uses the epoch-level `state` column already produced by
#'   zeitR's own pipelines ([run_pipeline()] / [run_pipeline_native()]):
#'   contiguous `state == 1` runs are treated as bouts, filtered by
#'   `duration_range`. Off-wrist (`state == 4`) epochs break a run rather
#'   than being bridged.
#' * `"roenneberg"` -- ignores any existing `state` column and runs the
#'   independent [detect_lids_bouts()] relative-immobility detector directly
#'   on the activity signal. Use this for standalone recordings that haven't
#'   been run through zeitR's Crespo/Vallim pipelines, or to reproduce
#'   Winnebeck/Hammad's own bout-detection method rather than zeitR's.
#' * `"auto"` (the default) -- `"state"` if a `state` column is present in
#'   `x`, otherwise `"roenneberg"`.
#'
#' # Quality filtering
#' Following Winnebeck et al. (2018) and Hammad et al. (2026), a bout
#' *passes* quality filtering when all of:
#' * `pearson_r >= min_r` (default `0.4` -- a soft data-quality threshold,
#'   not a hard significance test; ~75% of adult bouts cleared this bar in
#'   Winnebeck et al. 2018),
#' * `p_value <= max_p` (default `0.05`),
#' * `offset_bounds[1] < offset < offset_bounds[2]` (default `1 < offset <
#'   99`), excluding spuriously flat bouts (e.g. a lost/removed device).
#'
#' Bouts failing quality filtering are still returned (with
#' `passes_quality_filter = FALSE`) rather than dropped, so callers can
#' inspect what was excluded.
#'
#' @param x A `zeitr_result`, `zeitr_recording`, or a data frame / tibble
#'   with at least `datetime` and `activity` columns (and, for
#'   `bout_source = "state"`/`"auto"`, a `state` column).
#' @param bout_source `character(1)`. `"auto"` (default), `"state"`, or
#'   `"roenneberg"`. See Details.
#' @param activity_col `character(1)`. Name of the activity column in `x`.
#'   Default `"activity"`.
#' @param method,win_min,sigma_min Forwarded to [lids_transform()].
#' @param period_range,period_step Forwarded to [fit_lids()].
#' @param duration_range `numeric(2)`, hours. Bout duration bounds (both
#'   `bout_source` paths). Default `c(3, 12)`.
#' @param min_r,max_p,offset_bounds Quality-filter thresholds; see Details.
#' @param bout_args Named list of additional arguments forwarded to
#'   [detect_lids_bouts()] when `bout_source = "roenneberg"` (e.g.
#'   `relative_threshold`, `main_window`). Default `list()`.
#'
#' @return A tibble with one row per bout: `participant_id`, `bout_id`,
#'   `bout_start`, `bout_end`, `duration_h`, `period_min`, `amplitude`,
#'   `offset`, `slope_per_60min`, `phase_rad`, `pearson_r`, `p_value`, `mri`,
#'   `passes_quality_filter`.
#'
#' @references
#' Winnebeck, E. C., Fischer, D., Leise, T., & Roenneberg, T. (2018). Dynamics
#' and Ultradian Structure of Human Sleep in Real Life. *Current Biology*,
#' 28(1), 49-59.e5. \doi{10.1016/j.cub.2017.11.063}
#'
#' Hammad, G., Schoch, S. F., Engelmann, M., Spock, Z., Kurth, S., & Winnebeck,
#' E. C. (2026). Charting infant sleep cycle development using actigraphy:
#' Longitudinal evidence for ultradian cycle lengthening within the first
#' year of life. *SLEEP*.
#'
#' @seealso [lids_transform()], [fit_lids()], [detect_lids_bouts()],
#'   [study_lids_metrics()]
#'
#' @export
#'
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' result <- run_pipeline_native("recordings/P001.txt", tz = "America/Sao_Paulo")
#' compute_lids(result)
#' }
compute_lids <- function(x,
                          bout_source    = c("auto", "state", "roenneberg"),
                          activity_col   = "activity",
                          method         = c("gaussian", "mva"),
                          win_min        = 30,
                          sigma_min      = 5,
                          period_range   = c(30, 180),
                          period_step    = 2,
                          duration_range = c(3, 12),
                          min_r          = 0.4,
                          max_p          = 0.05,
                          offset_bounds  = c(1, 99),
                          bout_args      = list()) {

  bout_source <- match.arg(bout_source)
  method      <- match.arg(method)

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

  required <- c("datetime", activity_col)
  missing  <- setdiff(required, names(epochs))
  if (length(missing) > 0L) zeitr_abort("Missing required column(s): {.val {missing}}")

  has_state <- "state" %in% names(epochs)
  if (identical(bout_source, "auto")) bout_source <- if (has_state) "state" else "roenneberg"
  if (identical(bout_source, "state") && !has_state) {
    zeitr_abort(
      'No {.val state} column available for {.code bout_source = "state"};
       use {.code bout_source = "roenneberg"} for raw activity input.'
    )
  }

  datetime <- as.POSIXct(epochs$datetime)
  activity <- as.double(epochs[[activity_col]])
  ord      <- order(datetime)
  datetime <- datetime[ord]
  activity <- activity[ord]
  state    <- if (has_state) as.integer(epochs$state)[ord] else NULL

  diffs     <- as.numeric(diff(datetime), units = "secs")
  epoch_s   <- stats::median(diffs[diffs > 0], na.rm = TRUE)
  epoch_min <- epoch_s / 60

  empty_out <- .lids_na_tibble(participant_id)

  # ── Extract bouts ──────────────────────────────────────────────────────────
  if (bout_source == "state") {
    bouts <- .state_bouts(datetime, state, duration_range)
  } else {
    args      <- utils::modifyList(list(duration_range = duration_range), bout_args)
    bouts_tbl <- do.call(
      detect_lids_bouts,
      c(list(datetime = datetime, activity = activity), args)
    )
    if (nrow(bouts_tbl) == 0L) return(empty_out)
    bouts <- lapply(seq_len(nrow(bouts_tbl)), function(i) {
      idx <- which(datetime >= bouts_tbl$bout_start[i] & datetime <= bouts_tbl$bout_end[i])
      list(start = bouts_tbl$bout_start[i], end = bouts_tbl$bout_end[i], idx = idx)
    })
  }

  if (length(bouts) == 0L) return(empty_out)

  rows <- lapply(seq_along(bouts), function(i) {
    b     <- bouts[[i]]
    act_b <- activity[b$idx]
    if (length(act_b) < 4L || all(is.na(act_b))) return(NULL)

    lids <- lids_transform(act_b, epoch_min = epoch_min, method = method,
                            win_min = win_min, sigma_min = sigma_min)
    fit  <- tryCatch(
      fit_lids(lids, epoch_min = epoch_min, period_range = period_range,
                period_step = period_step),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)

    passes <- isTRUE(fit$pearson_r >= min_r) &&
      isTRUE(fit$p_value <= max_p) &&
      isTRUE(fit$offset > offset_bounds[1] && fit$offset < offset_bounds[2])

    tibble::tibble(
      participant_id         = participant_id,
      bout_id                = i,
      bout_start             = b$start,
      bout_end               = b$end,
      duration_h             = as.numeric(difftime(b$end, b$start, units = "hours")),
      period_min             = fit$period_min,
      amplitude              = fit$amplitude,
      offset                 = fit$offset,
      slope_per_60min        = fit$slope_per_60min,
      phase_rad              = fit$phase_rad,
      pearson_r              = fit$pearson_r,
      p_value                = fit$p_value,
      mri                    = fit$mri,
      passes_quality_filter  = passes
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(empty_out)
  do.call(rbind, rows)
}

# ── Study-level wrapper ──────────────────────────────────────────────────────

#' @noRd
.lids_summary_na_row <- function(pid, n_bouts, n_pass = 0L) {
  tibble::tibble(
    participant_id = pid, n_bouts = n_bouts, n_bouts_passing = n_pass,
    period_min_median = NA_real_, period_min_iqr = NA_real_,
    amplitude_median = NA_real_, amplitude_iqr = NA_real_,
    offset_median = NA_real_, offset_iqr = NA_real_,
    slope_per_60min_median = NA_real_, slope_per_60min_iqr = NA_real_
  )
}

#' Batch LIDS ultradian-rhythm summary across a study
#'
#' Computes [compute_lids()] for every participant in a batch of pipeline
#' results and summarises each participant's quality-filtered bouts into a
#' single row -- the LIDS counterpart to [study_sleep_metrics()] and
#' [study_summary()], making per-participant median cycle length, amplitude,
#' offset, and slope database-ready for tools like `syncR::sync()`.
#'
#' @param results A named list of `zeitr_result` objects, as returned by
#'   [run_pipeline_batch()] or [run_pipeline_native_batch()].
#' @param min_bouts `integer(1)`. Minimum number of quality-filtered bouts
#'   required for a participant's summary metrics to be reported (rather
#'   than `NA`, with `n_bouts`/`n_bouts_passing` still populated). Default
#'   `3`.
#' @param ... Forwarded to [compute_lids()] for every participant.
#'
#' @return A tibble with one row per participant: `participant_id`,
#'   `n_bouts`, `n_bouts_passing`, and (computed only over
#'   `passes_quality_filter == TRUE` bouts, and only when there are at
#'   least `min_bouts` of them) `period_min_median`, `amplitude_median`,
#'   `offset_median`, `slope_per_60min_median`, plus an `_iqr` variant of
#'   each.
#'
#' @seealso [compute_lids()], [study_sleep_metrics()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results <- run_pipeline_native_batch("recordings/", tz = "America/Sao_Paulo")
#' study_lids_metrics(results)
#' }
study_lids_metrics <- function(results, min_bouts = 3, ...) {
  if (!is.list(results) || length(results) == 0L) {
    zeitr_abort(
      "{.arg results} must be a non-empty list of {.cls zeitr_result} objects."
    )
  }

  list_names <- names(results)
  if (is.null(list_names)) list_names <- rep(NA_character_, length(results))

  rows <- lapply(seq_along(results), function(i) {
    result <- results[[i]]
    if (!inherits(result, "zeitr_result")) {
      zeitr_warn("Skipping {.val {list_names[i]}}: not a {.cls zeitr_result}.")
      return(NULL)
    }
    pid <- result$subject_id %||% list_names[i]

    bouts <- tryCatch(compute_lids(result, ...), error = function(e) {
      zeitr_warn("compute_lids() failed for {.val {pid}}: {conditionMessage(e)}")
      NULL
    })

    if (is.null(bouts) || nrow(bouts) == 0L) return(.lids_summary_na_row(pid, n_bouts = 0L))

    passing <- bouts[bouts$passes_quality_filter, ]
    n_bouts <- nrow(bouts)
    n_pass  <- nrow(passing)

    if (n_pass < min_bouts) return(.lids_summary_na_row(pid, n_bouts = n_bouts, n_pass = n_pass))

    tibble::tibble(
      participant_id         = pid,
      n_bouts                = n_bouts,
      n_bouts_passing        = n_pass,
      period_min_median      = stats::median(passing$period_min),
      period_min_iqr         = stats::IQR(passing$period_min),
      amplitude_median       = stats::median(passing$amplitude),
      amplitude_iqr          = stats::IQR(passing$amplitude),
      offset_median          = stats::median(passing$offset),
      offset_iqr             = stats::IQR(passing$offset),
      slope_per_60min_median = stats::median(passing$slope_per_60min),
      slope_per_60min_iqr    = stats::IQR(passing$slope_per_60min)
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    zeitr_warn("No valid results found.")
    return(tibble::tibble(
      participant_id = character(), n_bouts = integer(), n_bouts_passing = integer(),
      period_min_median = double(), period_min_iqr = double(),
      amplitude_median = double(), amplitude_iqr = double(),
      offset_median = double(), offset_iqr = double(),
      slope_per_60min_median = double(), slope_per_60min_iqr = double()
    ))
  }
  do.call(rbind, rows)
}
