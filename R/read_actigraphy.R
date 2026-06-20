#' Read an actigraphy file into a zeitr_recording object
#'
#' A device-agnostic wrapper that reads a raw actigraphy file and returns a
#' `zeitr_recording` object with `$epochs` (a tidy tibble) and `$metadata`
#' (a named list of device and recording information).
#'
#' Currently supported devices:
#' \itemize{
#'   \item `"acttrust"` — Condor Instruments ActTrust / ActTrust2 (`.txt`)
#' }
#'
#' @param path `character(1)`. Path to the raw actigraphy file.
#' @param device `character(1)`. Device type. One of `"acttrust"` (default).
#'   Additional devices will be added in future versions.
#' @param tz `character(1)`. Recording time zone passed to the underlying
#'   reader. Default is `"UTC"`.
#' @param ... Additional arguments forwarded to the device-specific reader
#'   (e.g. `encoding` for [read_acttrust()]).
#'
#' @return A `zeitr_recording` S3 object — a named list with:
#'   \describe{
#'     \item{`$epochs`}{A tibble with one row per epoch and columns
#'       `datetime`, `activity`, `int_temp`, `ext_temp`, `ZCMn`,
#'       `state`, `offwrist`, `sleep`.}
#'     \item{`$metadata`}{A named list with `subject`, `device_id`,
#'       `device_model`, `firmware_version`, `interval_s`, `source_file`,
#'       and `participant_id` (derived from the filename stem).}
#'   }
#'
#' @seealso [read_actigraphy_dir()] to read a whole directory at once.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rec <- read_actigraphy("recordings/P001.txt")
#' rec$epochs
#' rec$metadata
#' }
read_actigraphy <- function(path, device = "acttrust", tz = "UTC", ...) {
  device <- tolower(trimws(device))

  raw <- switch(device,
    acttrust = read_acttrust(path, tz = tz, ...),
    zeitr_abort(
      "Unsupported device {.val {device}}.
       Currently supported: {.val acttrust}."
    )
  )

  participant_id <- tools::file_path_sans_ext(basename(as.character(path)))

  meta           <- attr(raw, "metadata")
  meta$participant_id <- participant_id

  # Strip the zeitr_recording class from the raw tibble before embedding
  class(raw) <- setdiff(class(raw), "zeitr_recording")
  attr(raw, "metadata") <- NULL

  rec <- structure(
    list(
      epochs   = raw,
      metadata = meta
    ),
    class = "zeitr_recording"
  )

  rec
}

#' Read all actigraphy files in a directory
#'
#' Applies [read_actigraphy()] to every file matching `pattern` in `folder`.
#' Returns a `zeitr_study` object — a named list of `zeitr_recording` objects,
#' one per file. Files that fail to parse are skipped with a warning.
#'
#' @param folder `character(1)`. Path to a directory containing actigraphy
#'   files.
#' @param device `character(1)`. Device type passed to [read_actigraphy()].
#'   Default is `"acttrust"`.
#' @param pattern `character(1)`. Glob pattern for file discovery. Default is
#'   `"*.txt"`.
#' @param tz `character(1)`. Recording time zone. Default is `"UTC"`.
#' @param ... Additional arguments forwarded to [read_actigraphy()].
#'
#' @return A `zeitr_study` S3 object — a named list of `zeitr_recording`
#'   objects. Names are participant IDs (filename stems).
#'
#' @seealso [study_summary()] to summarise a `zeitr_study`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' study <- read_actigraphy_dir("recordings/", tz = "America/Sao_Paulo")
#' study_summary(study)
#' }
read_actigraphy_dir <- function(
    folder,
    device  = "acttrust",
    pattern = "*.txt",
    tz      = "UTC",
    ...
) {
  folder <- as.character(folder)
  files  <- Sys.glob(file.path(folder, pattern))

  if (length(files) == 0L) {
    zeitr_warn(
      "No files matching {.val {pattern}} found in {.path {folder}}."
    )
    return(structure(list(), class = "zeitr_study"))
  }

  recordings <- list()

  for (f in files) {
    pid <- tools::file_path_sans_ext(basename(f))
    tryCatch(
      {
        rec <- read_actigraphy(f, device = device, tz = tz, ...)
        recordings[[pid]] <- rec
        zeitr_inform("Read {.val {pid}} ({nrow(rec$epochs)} epochs).")
      },
      error = function(e) {
        zeitr_warn(
          "Skipping {.path {basename(f)}}: {conditionMessage(e)}"
        )
      }
    )
  }

  structure(recordings, class = "zeitr_study")
}

# ── S3 print methods ─────────────────────────────────────────────────────────

#' @export
print.zeitr_recording <- function(x, ...) {
  m <- x$metadata
  cli::cli_h1("zeitr_recording: {m$participant_id %||% 'unknown'}")
  cli::cli_bullets(c(
    "*" = "Device:   {m$device_model %||% 'unknown'} (ID {m$device_id %||% '?'})",
    "*" = "Firmware: {m$firmware_version %||% 'unknown'}",
    "*" = "Interval: {m$interval_s %||% '?'} s",
    "*" = "Epochs:   {nrow(x$epochs)}",
    "*" = "From:     {format(x$epochs$datetime[1])}",
    "*" = "To:       {format(x$epochs$datetime[nrow(x$epochs)])}"
  ))
  invisible(x)
}

#' @export
print.zeitr_study <- function(x, ...) {
  cli::cli_h1("zeitr_study: {length(x)} recording(s)")
  for (nm in names(x)) {
    cli::cli_bullets(c("*" = "{nm}: {nrow(x[[nm]]$epochs)} epochs"))
  }
  invisible(x)
}
