#' Summarise a zeitr_study across participants
#'
#' Computes per-participant summary statistics from a `zeitr_study` object
#' (as returned by [read_actigraphy_dir()]). For each recording, the function
#' computes NPCRA variables (IS, IV, RA, L5, M10) and basic recording quality
#' metrics.
#'
#' @param study A `zeitr_study` object as returned by [read_actigraphy_dir()],
#'   or a named list of `zeitr_recording` objects.
#' @param epoch_s `numeric(1)`. Epoch duration in seconds. If `NULL`
#'   (default), estimated separately for each recording.
#' @param L5_hours `numeric(1)`. Width of the L5 window in hours.
#'   Default is `5`.
#' @param M10_hours `numeric(1)`. Width of the M10 window in hours.
#'   Default is `10`.
#'
#' @return A tibble with one row per participant and columns:
#'   \describe{
#'     \item{`participant_id`}{Participant identifier (filename stem).}
#'     \item{`n_epochs`}{Total number of epochs in the recording.}
#'     \item{`n_days`}{Recording duration in days.}
#'     \item{`start`}{`POSIXct` — first epoch timestamp.}
#'     \item{`end`}{`POSIXct` — last epoch timestamp.}
#'     \item{`IS`}{Interdaily stability.}
#'     \item{`IV`}{Intradaily variability.}
#'     \item{`RA`}{Relative amplitude.}
#'     \item{`L5`}{Mean activity in the least active 5 h window.}
#'     \item{`L5_onset`}{Clock time of L5 midpoint.}
#'     \item{`M10`}{Mean activity in the most active 10 h window.}
#'     \item{`M10_onset`}{Clock time of M10 midpoint.}
#'   }
#'
#' @seealso [compute_npcra()] for single-recording NPCRA, [read_actigraphy_dir()]
#'   to create a `zeitr_study`.
#'
#' @export
#'
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' study <- read_actigraphy_dir("recordings/", tz = "America/Sao_Paulo")
#' study_summary(study)
#' }
study_summary <- function(study, epoch_s = NULL, L5_hours = 5, M10_hours = 10) {
  if (!is.list(study) || length(study) == 0L) {
    zeitr_abort(
      "{.arg study} must be a non-empty list of {.cls zeitr_recording} objects.
       Did you run {.fn read_actigraphy_dir}?"
    )
  }

  rows <- lapply(names(study), function(pid) {
    rec <- study[[pid]]

    if (!inherits(rec, "zeitr_recording")) {
      zeitr_warn("Skipping {.val {pid}}: not a {.cls zeitr_recording}.")
      return(NULL)
    }

    epochs    <- rec$epochs
    datetimes <- as.POSIXct(epochs$datetime)
    n         <- nrow(epochs)

    # Estimate epoch_s if not supplied
    ep_s <- epoch_s
    if (is.null(ep_s)) {
      diffs <- as.numeric(diff(datetimes), units = "secs")
      ep_s  <- stats::median(diffs[diffs > 0], na.rm = TRUE)
    }

    n_days <- n / (24 * 3600 / ep_s)

    npcra <- tryCatch(
      compute_npcra(rec, epoch_s = ep_s,
                    L5_hours = L5_hours, M10_hours = M10_hours),
      error = function(e) {
        zeitr_warn("NPCRA failed for {.val {pid}}: {conditionMessage(e)}")
        tibble::tibble(IS = NA_real_, IV = NA_real_, RA = NA_real_,
                       L5 = NA_real_, L5_onset = NA_character_,
                       M10 = NA_real_, M10_onset = NA_character_)
      }
    )

    tibble::tibble(
      participant_id = pid,
      n_epochs       = n,
      n_days         = round(n_days, 2),
      start          = datetimes[1L],
      end            = datetimes[n],
      IS             = npcra$IS,
      IV             = npcra$IV,
      RA             = npcra$RA,
      L5             = npcra$L5,
      L5_onset       = npcra$L5_onset,
      M10            = npcra$M10,
      M10_onset      = npcra$M10_onset
    )
  })

  rows <- Filter(Negate(is.null), rows)

  if (length(rows) == 0L) {
    zeitr_warn("No valid recordings found in study.")
    return(tibble::tibble(
      participant_id = character(),
      n_epochs       = integer(),
      n_days         = double(),
      start          = as.POSIXct(character()),
      end            = as.POSIXct(character()),
      IS = double(), IV = double(), RA = double(),
      L5 = double(), L5_onset = character(),
      M10 = double(), M10_onset = character()
    ))
  }

  do.call(rbind, rows)
}
