#' Read an Axivity AX3/AX6 .cwa file into a zeitR-standard epoch tibble
#'
#' Bridges [axR::axivity_read_cwa()]'s raw per-sample output (triaxial
#' acceleration at the device's native sampling rate) into the same
#' epoch-level, 9-column shape [read_acttrust()] produces -- so an Axivity
#' recording can flow through the rest of zeitR ([run_pipeline()],
#' [compute_npcra()], [compute_sri()], etc.) exactly like an ActTrust one.
#'
#' @details
#' # What this does (and doesn't) validate
#' The raw-to-counts conversion itself is [compute_activity_counts()] --
#' already documented there as an unvalidated approximation of onboard
#' device processing (no reference converter exists to check it against).
#' Axivity devices additionally have **no published or validated
#' filter/threshold preset** at all (unlike ActTrust and GT3X+, which at
#' least have a documented processing description to approximate). The
#' `filter_low`/`filter_high` defaults here (`0.25`/`2.5` Hz) are
#' [compute_activity_counts()]'s GT3X+-style preset, reused because it's
#' the closer starting point of the two existing options for a
#' research-grade wrist accelerometer like the AX3 -- **not** because it has
#' been checked against real Axivity output. Treat `activity`/`ZCMn` from
#' this function as a rough approximation only; validate against a
#' reference (e.g. GGIR) before relying on it for any published analysis.
#'
#' `ZCMn` is `compute_activity_counts()`'s raw `ZCM` count with no
#' additional normalisation applied -- named `ZCMn` only for column-name
#' compatibility with [read_acttrust()]'s CK-scoring input, not because a
#' normalisation step has actually been performed.
#'
#' # Sampling rate
#' `axivity_read_cwa()` reports `sample_rate` per sample (block-level, as
#' stored in the `.cwa` file). This function takes the single most common
#' value across the whole recording and uses it for the entire conversion;
#' if any samples report a different rate (a genuine rate change mid
#' recording, or a corrupt block), a warning names the discrepancy but the
#' dominant rate is still used throughout. Epoch boundaries and grouping for
#' `light`/`int_temp` averaging are derived from this same dominant rate,
#' matching [compute_activity_counts()]'s own epoch grouping exactly
#' (including which trailing samples get dropped).
#'
#' # Time zone
#' `axivity_read_cwa()` tags `timestamp` as UTC by convention (the device's
#' own real-time clock, not a true UTC source) -- exactly like
#' [read_acttrust()] does for ActTrust's `DATE/TIME` column. `tz` here
#' re-labels the same clock reading under the recording's actual local time
#' zone (via [lubridate::force_tz()]; no shift in wall-clock value), rather
#' than converting it -- set it to the time zone the device's clock was
#' actually set to for correct circadian alignment downstream.
#'
#' @param path `character(1)`. Path to a `.cwa`/AX6 file, forwarded to
#'   [axR::axivity_read_cwa()].
#' @param tz `character(1)`. Time zone the device's clock was set to.
#'   Default `"UTC"`.
#' @param epoch_sec `numeric(1)`. Epoch length in seconds, forwarded to
#'   [compute_activity_counts()]. Default `60`, matching the rest of zeitR.
#' @param filter_low,filter_high `numeric(1)`. Band-pass cutoffs in Hz,
#'   forwarded to [compute_activity_counts()]. Default `0.25`/`2.5`
#'   (GT3X+-style preset -- see Details; no validated Axivity-specific
#'   preset exists).
#' @param zcm_threshold,tat_threshold `numeric(1)`. Forwarded to
#'   [compute_activity_counts()]. Defaults `0.01`/`0.05`.
#'
#' @return A tibble with one row per epoch and the same columns as
#'   [read_acttrust()] (`datetime`, `activity`, `int_temp`, `ext_temp`,
#'   `ZCMn`, `light`, `state`, `offwrist`, `sleep`), plus one extra column
#'   not present there: `TAT` (time above threshold, seconds/epoch, from
#'   [compute_activity_counts()]). `activity` is PIM. `ext_temp` is always
#'   `NA` -- Axivity devices have a single on-body temperature sensor, no
#'   separate ambient sensor. The tibble carries a `"zeitr_axivity"` class
#'   and a `metadata` attribute (a named list with `device_id`,
#'   `session_id`, `sample_rate`, `epoch_sec`, `filter_low`, `filter_high`,
#'   `cwa_metadata` (axR's raw device metadata string), `source_file`).
#'
#' @seealso [read_actigraphy()] (`device = "axivity"`),
#'   [compute_activity_counts()] for the underlying raw-to-counts
#'   conversion and its validation caveats, [read_acttrust()] for the
#'   column shape this matches.
#'
#' @export
#'
#' @importFrom lubridate force_tz
#' @importFrom tibble as_tibble
#'
#' @examples
#' \dontrun{
#' rec <- read_axivity("recordings/P001.cwa", tz = "America/Sao_Paulo")
#' rec
#' attr(rec, "metadata")
#' }
read_axivity <- function(path, tz = "UTC",
                          epoch_sec     = 60,
                          filter_low    = 0.25,
                          filter_high   = 2.5,
                          zcm_threshold = 0.01,
                          tat_threshold = 0.05) {

  .check_axR_pkg()

  path <- as.character(path)
  if (!file.exists(path)) {
    zeitr_abort("File not found: {.path {path}}")
  }

  raw <- axR::axivity_read_cwa(path)
  if (nrow(raw) < 2L) {
    zeitr_abort("{.path {path}} contains fewer than 2 samples.")
  }

  # ── Dominant sampling rate ──────────────────────────────────────────────
  rate_tbl      <- table(raw$sample_rate)
  sampling_rate <- as.numeric(names(rate_tbl)[which.max(rate_tbl)])
  n_other_rates <- sum(raw$sample_rate != sampling_rate)
  if (n_other_rates > 0L) {
    zeitr_warn(
      "{n_other_rates} sample(s) in {.path {path}} report a sample_rate other than the dominant {sampling_rate} Hz; using {sampling_rate} Hz for the whole conversion."
    )
  }

  # ── Raw -> epoch-level PIM/TAT/ZCM ───────────────────────────────────────
  counts <- compute_activity_counts(
    raw$x, raw$y, raw$z,
    sampling_rate = sampling_rate,
    epoch_sec     = epoch_sec,
    filter_low    = filter_low,
    filter_high   = filter_high,
    zcm_threshold = zcm_threshold,
    tat_threshold = tat_threshold,
    metrics       = c("PIM", "TAT", "ZCM")
  )

  # ── Matching epoch grouping for light / temperature, and epoch datetimes ──
  # Mirrors compute_activity_counts()'s own truncation of trailing samples
  # exactly, so light/int_temp line up with the same epochs as PIM/TAT/ZCM.
  samples_per_epoch <- round(epoch_sec * sampling_rate)
  n_epochs          <- nrow(counts)
  n_use             <- n_epochs * samples_per_epoch
  epoch_id          <- rep(seq_len(n_epochs), each = samples_per_epoch)

  light    <- tapply(raw$light[seq_len(n_use)],           epoch_id, mean, na.rm = TRUE)
  int_temp <- tapply(raw$temperature_c[seq_len(n_use)],    epoch_id, mean, na.rm = TRUE)

  datetime <- raw$timestamp[1L] + epoch_sec * (seq_len(n_epochs) - 1L)
  datetime <- lubridate::force_tz(datetime, tzone = tz)

  out <- tibble::tibble(
    datetime = datetime,
    activity = counts$PIM,
    int_temp = as.numeric(int_temp),
    ext_temp = NA_real_,
    ZCMn     = counts$ZCM,
    light    = as.numeric(light),
    state    = 0,
    offwrist = 0,
    sleep    = 0,
    TAT      = counts$TAT
  )

  metadata <- list(
    device_id    = attr(raw, "device_id")  %||% NA_character_,
    session_id   = attr(raw, "session_id") %||% NA_character_,
    sample_rate  = sampling_rate,
    epoch_sec    = epoch_sec,
    filter_low   = filter_low,
    filter_high  = filter_high,
    cwa_metadata = attr(raw, "metadata")   %||% NA_character_,
    source_file  = normalizePath(path, mustWork = FALSE)
  )

  class(out) <- c("zeitr_axivity", class(out))
  attr(out, "metadata") <- metadata

  out
}

#' @noRd
.check_axR_pkg <- function() {
  if (!requireNamespace("axR", quietly = TRUE))
    zeitr_abort(c(
      "Package {.pkg axR} is required to read Axivity {.file .cwa} files ({.fn axR::axivity_read_cwa}).",
      "i" = 'Install from r-universe: {.code install.packages("axR", repos = c("https://circadia-bio.r-universe.dev", "https://cloud.r-project.org"))}'
    ))
}
