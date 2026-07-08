# ── Export to hypnoR ──────────────────────────────────────────────────────────

#' Export a zeitR pipeline result as a hypnoR-compatible hypnogram
#'
#' Converts the epoch-level `data` tibble from [run_pipeline()] or
#' [run_pipeline_native()] into the tidy hypnogram format expected by
#' `hypnoR` metric functions.
#'
#' @details
#' The coarse (3-state) stage mapping used by `zeitR` is:
#'
#' | `state` | `ZCMn`  | Stage         |
#' |---------|---------|---------------|
#' | `0`     | any     | `"W"`         |
#' | `1`, `7`| `> 0`   | `"Sleep"`     |
#' | `1`, `7`| `== 0`  | `"Quiet sleep"`|
#' | `4`     | any     | `"W"`         |
#'
#' Zero-count epochs within sleep (`ZCMn == 0`) are mapped to `"Quiet sleep"`
#' as the standard actigraphy proxy for quiet/deep sleep. When `ZCMn` is not
#' present in the data, all sleep epochs are mapped to `"Sleep"`.
#'
#' The `stage` column is an ordered factor with levels
#' `c("W", "Sleep", "Quiet sleep")`, matching `hypnoR`'s coarse resolution
#' contract. `"Quiet sleep"` is not produced by actigraphy but the level is
#' present so downstream `hypnoR` functions do not throw factor-level errors.
#'
#' `subject_id` is taken from `result$subject_id` when not explicitly
#' supplied. Both [run_pipeline()] and [run_pipeline_native()] derive this
#' from the input filename stem automatically, so in most cases no manual
#' override is needed.
#'
#' @param result A `zeitr_result` list as returned by [run_pipeline()] or
#'   [run_pipeline_native()], or a tibble with at minimum the columns
#'   `datetime` and `state`.
#' @param subject_id `character(1)` or `NULL`. Written to the `subject_id`
#'   column and always takes precedence. When `NULL` (default),
#'   `result$subject_id` is used if present; the column is omitted entirely
#'   when no ID is available from either source.
#' @param source `character(1)`. Label written to the `source` column.
#'   Default `"zeitR"`.
#' @param drop_offwrist `logical(1)`. Remove off-wrist epochs
#'   (`state == 4`) from the output and re-index epochs. Default `FALSE`.
#' @param epoch_sec `numeric(1)`. When supplied, warns if the observed epoch
#'   duration differs from this value. Default `NULL`.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{`epoch`}{Integer epoch index, 1-based.}
#'     \item{`time`}{`POSIXct` timestamp for the start of each epoch.}
#'     \item{`stage`}{Ordered factor: `c("W", "Sleep", "Quiet sleep")`.}
#'     \item{`subject_id`}{Character subject identifier (omitted when not
#'       available).}
#'     \item{`source`}{Character scorer label.}
#'   }
#'
#' @seealso [run_pipeline()], [run_pipeline_native()], [label_states()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- run_pipeline("recordings/P001.txt", tz = "America/Sao_Paulo")
#'
#' # subject_id inferred from filename automatically
#' hyp <- export_hypnogram(result)
#' hyp$subject_id  # "P001"
#'
#' # Override when filename does not match study code
#' hyp <- export_hypnogram(result, subject_id = "STUDY_001")
#'
#' # Batch: subject_id inferred for every file automatically
#' results <- run_pipeline_batch("recordings/", tz = "America/Sao_Paulo")
#' hyps    <- lapply(results, export_hypnogram)
#'
#' # Drop off-wrist epochs
#' hyp_clean <- export_hypnogram(result, drop_offwrist = TRUE)
#' }
export_hypnogram <- function(result,
                             subject_id    = NULL,
                             source        = "zeitR",
                             drop_offwrist = FALSE,
                             epoch_sec     = NULL) {
  # Accept either a zeitr_result list or a bare data tibble
  data <- if (is.list(result) && !is.data.frame(result)) {
    if (is.null(result$data))
      zeitr_abort(
        "{.arg result} must be a `zeitr_result` list with a `$data` element, ",
        "or a tibble with `datetime` and `state` columns."
      )
    result$data
  } else {
    result
  }

  required <- c("datetime", "state")
  missing  <- setdiff(required, names(data))
  if (length(missing) > 0L)
    zeitr_abort("{.arg result} is missing required column(s): {.val {missing}}.")

  # Resolve subject_id: explicit arg > result$subject_id > NULL (omit column)
  sid <- if (!is.null(subject_id)) {
    subject_id
  } else if (is.list(result) && !is.null(result$subject_id)) {
    result$subject_id
  } else {
    NULL
  }

  state <- as.integer(data$state)
  zcm   <- if ("ZCMn" %in% names(data)) as.double(data$ZCMn) else NULL

  in_sleep <- state == 1L | state == 7L

  stage_chr <- dplyr::case_when(
    in_sleep & !is.null(zcm) & zcm == 0 ~ "Quiet sleep",
    in_sleep                             ~ "Sleep",
    TRUE                                 ~ "W"
  )

  stage_fct <- factor(stage_chr,
                      levels  = c("W", "Sleep", "Quiet sleep"),
                      ordered = TRUE)

  out <- tibble::tibble(
    epoch  = seq_len(nrow(data)),
    time   = as.POSIXct(data$datetime),
    stage  = stage_fct,
    source = source
  )

  if (!is.null(sid)) out$subject_id <- sid

  if (isTRUE(drop_offwrist)) {
    out       <- out[state != 4L, , drop = FALSE]
    out$epoch <- seq_len(nrow(out))
  }

  if (!is.null(epoch_sec)) {
    gaps     <- as.numeric(diff(as.POSIXct(data$datetime[state != 4L])), units = "secs")
    gaps     <- gaps[gaps > 0]
    if (length(gaps) > 0L) {
      mode_gap <- as.numeric(names(sort(table(gaps), decreasing = TRUE))[1L])
      if (abs(mode_gap - epoch_sec) > 1)
        zeitr_warn(
          "Epoch duration in data ({round(mode_gap)} s) differs from ",
          "`epoch_sec` ({epoch_sec} s)."
        )
    }
  }

  out
}
