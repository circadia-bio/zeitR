#' Run the full actigraphy sleep analysis pipeline
#'
#' Orchestrates the complete pipeline for a single ActTrust recording:
#'
#' 1. **Read** — [read_acttrust()]
#' 2. **Consistency check** — [check_consistency()]
#' 3. **Prepare** — [prepare_actigraphy()]
#' 4. **Off-wrist detection** — [detect_offwrist_bimodal()]
#' 5. **Main sleep period detection** — [detect_sleep_crespo()]
#' 6. **Nap detection** — [detect_naps_crespo()]
#' 7. **WASO + nightly statistics** — [compute_waso()]
#'
#' @param path `character(1)`. Path to the ActTrust `.txt` file.
#' @param tz `character(1)`. Recording time zone. Passed to [read_acttrust()].
#'   Default is `"UTC"`.
#' @param wake_thresh `integer(1)`. Minimum wake-bout length (epochs) used to
#'   separate sleep periods in [compute_waso()]. Default is `60`.
#' @param gap_s `numeric(1)`. Gap threshold (seconds) for [check_consistency()].
#'   Default is `120`.
#' @param offwrist_args `list`. Additional arguments passed to
#'   [detect_offwrist_bimodal()]. Default is an empty list.
#' @param sleep_args `list`. Additional arguments passed to
#'   [detect_sleep_crespo()]. Default is an empty list.
#' @param nap_args `list`. Additional arguments passed to
#'   [detect_naps_crespo()]. Default is an empty list.
#'
#' @return A `zeitr_result` S3 object — a named list with:
#'   \describe{
#'     \item{`subject_id`}{`character` — derived from the input filename stem.}
#'     \item{`source_file`}{`character` — absolute path to the input file.}
#'     \item{`data`}{`tibble` — final epoch-level data frame with all state
#'       columns populated.}
#'     \item{`nights`}{`tibble` — per-night sleep statistics.}
#'     \item{`issues`}{`tibble` — timestamp consistency issues (0 rows if none).}
#'     \item{`metadata`}{`list` — device and subject metadata from the file header.}
#'   }
#'
#' @seealso [run_pipeline_batch()] for processing a directory of files.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- run_pipeline("recordings/P001.txt", tz = "America/Sao_Paulo")
#' result$nights
#' result$data
#' }
run_pipeline <- function(
    path,
    tz             = "UTC",
    wake_thresh    = 60L,
    gap_s          = 120,
    offwrist_args  = list(),
    sleep_args     = list(),
    nap_args       = list()
) {
  path       <- as.character(path)
  subject_id <- tools::file_path_sans_ext(basename(path))

  zeitr_inform("Reading {.path {basename(path)}} ...")

  # 1. Read
  rec <- read_acttrust(path, tz = tz)

  # 2. Consistency check
  issues <- check_consistency(rec, gap_s = gap_s)
  if (nrow(issues) > 0L) {
    zeitr_warn(
      "[{subject_id}] {nrow(issues)} timestamp issue(s) detected. \\
       Check {.code result$issues} for details."
    )
  }

  # 3. Prepare
  prep <- prepare_actigraphy(rec)

  # 4. Off-wrist detection
  prep <- do.call(detect_offwrist_bimodal, c(list(x = prep), offwrist_args))

  # 5. Main sleep periods
  prep <- do.call(detect_sleep_crespo, c(list(x = prep), sleep_args))

  # 6. Naps
  prep <- do.call(detect_naps_crespo, c(list(x = prep), nap_args))

  # 7. WASO + nightly stats
  waso_result <- compute_waso(prep, wake_thresh = wake_thresh)

  zeitr_inform(
    "[{subject_id}] Done. {nrow(waso_result$nights)} night(s) detected."
  )

  result <- structure(
    list(
      subject_id  = subject_id,
      source_file = normalizePath(path, mustWork = FALSE),
      data        = waso_result$data,
      nights      = waso_result$nights,
      issues      = issues,
      metadata    = attr(rec, "metadata")
    ),
    class = "zeitr_result"
  )

  result
}

#' Run the pipeline on all files in a directory
#'
#' Applies [run_pipeline()] to every file matching `pattern` in `folder`,
#' returning a list of `zeitr_result` objects. Files that fail are skipped
#' with a warning.
#'
#' @param folder `character(1)`. Path to a directory containing ActTrust files.
#' @param pattern `character(1)`. Glob pattern for file discovery. Default is
#'   `"*.txt"`.
#' @param ... Additional arguments forwarded to [run_pipeline()].
#'
#' @return A named list of `zeitr_result` objects, one per successfully
#'   processed file. Names are the file stem (subject IDs).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results <- run_pipeline_batch("recordings/", tz = "America/Sao_Paulo")
#' lapply(results, function(r) r$nights)
#' }
run_pipeline_batch <- function(folder, pattern = "*.txt", ...) {
  folder <- as.character(folder)
  files  <- Sys.glob(file.path(folder, pattern))

  if (length(files) == 0L) {
    zeitr_warn("No files matching {.val {pattern}} found in {.path {folder}}.")
    return(list())
  }

  results <- list()

  for (f in files) {
    subject_id <- tools::file_path_sans_ext(basename(f))
    tryCatch(
      {
        res <- run_pipeline(f, ...)
        results[[subject_id]] <- res
      },
      error = function(e) {
        zeitr_warn("Failed to process {.path {basename(f)}}: {conditionMessage(e)}")
      }
    )
  }

  results
}

# ── S3 methods ───────────────────────────────────────────────────────────────

#' @export
print.zeitr_result <- function(x, ...) {
  cli::cli_h1("zeitr_result: {x$subject_id}")
  cli::cli_bullets(c(
    "*" = "Source:  {x$source_file}",
    "*" = "Epochs:  {nrow(x$data)}",
    "*" = "Nights:  {nrow(x$nights)}",
    "*" = "Issues:  {nrow(x$issues)}"
  ))
  invisible(x)
}
