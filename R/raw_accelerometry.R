# ── Raw accelerometry -> PIM / TAT / ZCM activity counts ─────────────────────
#
# Converts raw triaxial acceleration samples into the epoch-level activity
# metrics (PIM, TAT, ZCM) that wrist actigraphs like ActTrust and GT3X+
# normally compute onboard and export directly (see read_acttrust()).
#
# Filtering reuses mrpheus::remove_dc() and mrpheus::bandpass_filter() --
# the same zero-phase Butterworth implementation validated as part of
# mrpheus's YASA-parity PSG pipeline -- rather than reimplementing filtering
# from scratch with a new package dependency. This makes mrpheus a
# cross-package Suggests dependency of zeitR for this function specifically
# (precedented by hypnoR, which already Suggests both mrpheus and zeitR);
# the rest of zeitR remains fully independent of mrpheus.
#
# IMPORTANT: reusing a validated filter primitive does not make the whole
# function validated. The PIM/TAT/ZCM epoch-aggregation logic below (what
# happens *after* filtering) is still an original implementation with no
# reference to check against -- no raw-to-counts converter exists in
# Lucas's local ecosystem (condor_pipeline/circadiaBase_Docker's vendored
# Condor scripts only operate on already-epoched PIM/ZCM/TAT), and Condor's/
# ActiGraph's exact onboard threshold constants are proprietary and
# unpublished. Built from the general processing description in Batista et
# al. (2026, PLoS ONE, doi.org/10.1371/journal.pone.0348631) and standard
# actigraphy literature (Bouten et al. 1997; Ancoli-Israel et al. 2003).

# ── Internal helpers ──────────────────────────────────────────────────────────

#' @noRd
.check_mrpheus_pkg <- function() {
  if (!requireNamespace("mrpheus", quietly = TRUE))
    zeitr_abort(c(
      "Package {.pkg mrpheus} is required for raw-accelerometry filtering (reuses {.fn mrpheus::bandpass_filter} / {.fn mrpheus::remove_dc}).",
      "i" = 'Install from r-universe: {.code install.packages("mrpheus", repos = c("https://circadia-bio.r-universe.dev", "https://cloud.r-project.org"))}'
    ))
}

# Indicator of a zero crossing (beyond a dead-band) between each pair of
# consecutive samples. Returns a logical vector of length length(v) - 1;
# crossing[i] is TRUE if a crossing occurred between v[i] and v[i + 1].
# Samples inside the dead-band inherit the last defined zone (+1/-1), so a
# genuine crossing must clear the dead-band on both sides to be counted --
# this prevents high-frequency noise near zero from inflating the count.
#' @noRd
.zero_crossing_indicator <- function(v, threshold) {
  n    <- length(v)
  zone <- rep(NA_integer_, n)
  zone[v >  threshold] <-  1L
  zone[v < -threshold] <- -1L

  idx <- which(!is.na(zone))
  if (length(idx) < 2L) return(rep(FALSE, n - 1L))

  filled <- rep(NA_integer_, n)
  filled[idx[1]:n] <- rep(zone[idx], diff(c(idx, n + 1L)))

  d <- diff(filled)
  !is.na(d) & d != 0L
}

# Sum `values` within each epoch_id group, returned as a length-n_epochs
# vector (0 for any epoch with no contributing values, e.g. a ZCM epoch with
# zero crossings).
#' @noRd
.epoch_rowsum <- function(values, epoch_id, n_epochs) {
  s   <- rowsum(values, group = epoch_id)
  out <- numeric(n_epochs)
  out[as.integer(rownames(s))] <- s[, 1]
  out
}

# ── compute_activity_counts ───────────────────────────────────────────────────

#' Compute PIM, TAT, and ZCM activity counts from raw triaxial acceleration
#'
#' Reproduces the general wrist-actigraph processing chain -- per-axis
#' DC removal and band-pass filter, vector norm, epoch-level integration --
#' that converts raw acceleration samples into the epoch-level activity
#' metrics used throughout the rest of zeitR (e.g. [read_acttrust()]'s
#' `activity`/`ZCMn` columns) and the wider PA-intensity literature (see
#' [pa_equations()]).
#'
#' @details
#' # What this is (and isn't)
#' Filtering reuses [mrpheus::remove_dc()] and [mrpheus::bandpass_filter()]
#' -- the same zero-phase Butterworth implementation already validated as
#' part of mrpheus's YASA-parity PSG pipeline -- rather than a new,
#' unvalidated filter written from scratch. This makes **mrpheus a required
#' package for this function specifically** (a `Suggests` dependency of
#' zeitR, checked at runtime; the rest of zeitR does not need it).
#'
#' Reusing a validated filter primitive does not make the whole function
#' validated, though: the PIM/TAT/ZCM epoch-aggregation logic below (what
#' happens *after* filtering) is still an original implementation with
#' nothing to check it against -- no raw-to-counts converter exists to
#' compare it to, and Condor's/ActiGraph's exact onboard threshold constants
#' are proprietary and unpublished. Built from the general processing
#' description in Batista et al. (2026, PLoS ONE) and standard actigraphy
#' literature. Treat its output as a reasonable approximation of
#' device-computed counts, not a validated reproduction -- and prefer
#' device-computed counts ([read_acttrust()]) over this function whenever
#' they're available.
#'
#' # Processing steps
#' 1. Each of `x`, `y`, `z` has its DC offset removed
#'    ([mrpheus::remove_dc()]) -- important for the z-axis in particular,
#'    which typically carries a ~1 g gravity offset -- then is independently
#'    band-pass filtered with a zero-phase ([mrpheus::bandpass_filter()])
#'    Butterworth filter (`filter_low`-`filter_high` Hz).
#' 2. A vector norm is computed per sample: `sqrt(x^2 + y^2 + z^2)` on the
#'    *filtered* axes, matching the order described for ActTrust in Batista
#'    et al. (2026): filter first, then combine axes.
#' 3. Samples are grouped into non-overlapping epochs of
#'    `epoch_sec * sampling_rate` samples. Trailing samples that don't fill
#'    a complete final epoch are dropped, with a warning.
#' 4. Per epoch:
#'    - **PIM** (proportional integration mode) -- sum of the absolute
#'      filtered norm across the epoch.
#'    - **TAT** (time above threshold) -- seconds within the epoch where the
#'      absolute filtered norm exceeds `tat_threshold`.
#'    - **ZCM** (zero crossing mode) -- number of sign changes beyond a
#'      `zcm_threshold` dead-band, summed across the three filtered axes,
#'      within the epoch.
#'
#' `zcm_threshold` and `tat_threshold` have no published reference value --
#' they are exposed as tunable parameters rather than given a false-precision
#' default. The values shipped here are small, physically motivated starting
#' points (see Arguments), not calibrated constants.
#'
#' @param x,y,z Numeric vectors of raw triaxial acceleration samples (same
#'   units -- typically g -- and length, uniformly sampled at
#'   `sampling_rate`).
#' @param sampling_rate `numeric(1)`. Sampling rate of `x`/`y`/`z` in Hz.
#' @param epoch_sec `numeric(1)`. Epoch length in seconds. Default `60`
#'   (1-minute epochs, the ActTrust standard used elsewhere in zeitR).
#' @param filter_low,filter_high `numeric(1)`. Band-pass cutoffs in Hz,
#'   passed to [mrpheus::bandpass_filter()]. Default `0.5`/`2.7`
#'   (ActTrust-style, per Batista et al. 2026); use `0.25`/`2.5` for
#'   GT3X+-style processing instead.
#' @param zcm_threshold `numeric(1)`. Dead-band around zero (same units as
#'   `x`/`y`/`z`) a filtered axis must clear to register a zero crossing.
#'   Default `0.01`.
#' @param tat_threshold `numeric(1)`. Amplitude threshold (same units) the
#'   filtered norm must exceed to count toward TAT. Default `0.05`.
#' @param metrics `character()`. Which metrics to compute. Default
#'   `c("PIM", "TAT", "ZCM")` (all three).
#'
#' @return A tibble with one row per epoch: `epoch` (integer index, `1`-based)
#'   plus one column per requested metric.
#'
#' @seealso [pa_equations()], [estimate_ee()], [classify_pa_counts()] to go
#'   from PIM counts to METs and PA intensity bands; [read_acttrust()] for
#'   device-computed counts (preferred over this function when available);
#'   [mrpheus::bandpass_filter()], [mrpheus::remove_dc()] for the underlying
#'   filter primitives.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n  <- 25 * 60 * 5  # 5 minutes at 25 Hz
#' x  <- rnorm(n, sd = 0.05)
#' y  <- rnorm(n, sd = 0.05)
#' z  <- rnorm(n, sd = 0.05) + 1  # gravity on the z-axis
#'
#' \dontrun{
#' compute_activity_counts(x, y, z, sampling_rate = 25, epoch_sec = 60)
#' }
compute_activity_counts <- function(x, y, z,
                                     sampling_rate,
                                     epoch_sec     = 60,
                                     filter_low    = 0.5,
                                     filter_high   = 2.7,
                                     zcm_threshold = 0.01,
                                     tat_threshold = 0.05,
                                     metrics       = c("PIM", "TAT", "ZCM")) {

  .check_mrpheus_pkg()
  metrics <- match.arg(metrics, several.ok = TRUE)

  if (length(x) != length(y) || length(x) != length(z)) {
    zeitr_abort("{.arg x}, {.arg y}, and {.arg z} must be the same length.")
  }
  if (!is.numeric(sampling_rate) || length(sampling_rate) != 1L || sampling_rate <= 0) {
    zeitr_abort("{.arg sampling_rate} must be a single positive number.")
  }
  if (!is.numeric(epoch_sec) || length(epoch_sec) != 1L || epoch_sec <= 0) {
    zeitr_abort("{.arg epoch_sec} must be a single positive number.")
  }

  samples_per_epoch <- round(epoch_sec * sampling_rate)
  if (samples_per_epoch < 2L) {
    zeitr_abort(
      "epoch_sec x sampling_rate = {samples_per_epoch} sample(s) per epoch; need at least 2 to compute zero crossings meaningfully."
    )
  }

  n <- length(x)
  if (n < samples_per_epoch) {
    zeitr_abort(
      "Not enough samples ({n}) for even one {epoch_sec}s epoch at {sampling_rate} Hz ({samples_per_epoch} samples/epoch)."
    )
  }

  n_epochs <- n %/% samples_per_epoch
  n_use    <- n_epochs * samples_per_epoch
  if (n_use < n) {
    zeitr_warn("Dropping the last {n - n_use} sample(s) that don't fill a complete final epoch.")
  }

  keep <- seq_len(n_use)
  xf <- mrpheus::bandpass_filter(mrpheus::remove_dc(x[keep]), sr = sampling_rate, low_hz = filter_low, high_hz = filter_high)
  yf <- mrpheus::bandpass_filter(mrpheus::remove_dc(y[keep]), sr = sampling_rate, low_hz = filter_low, high_hz = filter_high)
  zf <- mrpheus::bandpass_filter(mrpheus::remove_dc(z[keep]), sr = sampling_rate, low_hz = filter_low, high_hz = filter_high)

  norm     <- sqrt(xf^2 + yf^2 + zf^2)
  epoch_id <- rep(seq_len(n_epochs), each = samples_per_epoch)

  out <- tibble::tibble(epoch = seq_len(n_epochs))

  if ("PIM" %in% metrics) {
    out$PIM <- .epoch_rowsum(abs(norm), epoch_id, n_epochs)
  }

  if ("TAT" %in% metrics) {
    above   <- as.numeric(abs(norm) > tat_threshold)
    out$TAT <- .epoch_rowsum(above, epoch_id, n_epochs) / sampling_rate
  }

  if ("ZCM" %in% metrics) {
    epoch_id_pairs <- epoch_id[-1]
    zx <- .epoch_rowsum(as.numeric(.zero_crossing_indicator(xf, zcm_threshold)), epoch_id_pairs, n_epochs)
    zy <- .epoch_rowsum(as.numeric(.zero_crossing_indicator(yf, zcm_threshold)), epoch_id_pairs, n_epochs)
    zz <- .epoch_rowsum(as.numeric(.zero_crossing_indicator(zf, zcm_threshold)), epoch_id_pairs, n_epochs)
    out$ZCM <- zx + zy + zz
  }

  out
}
