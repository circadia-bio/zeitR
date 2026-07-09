#' Check actigraphy timestamps for consistency issues
#'
#' Scans a recording tibble for three classes of timestamp problem:
#'
#' * **Gaps** — intervals between consecutive epochs longer than `gap_s`
#'   seconds.
#' * **Backward jumps** — timestamps that go backwards in time.
#' * **Year artefacts** — timestamps in the years 1970 or 2000, which
#'   typically indicate firmware epoch-counter rollover bugs.
#'
#' @param x A tibble as returned by [read_acttrust()] or
#'   [prepare_actigraphy()], containing a `datetime` column.
#' @param gap_s `numeric(1)`. Gap threshold in seconds. Intervals longer than
#'   this are flagged. Default is `120` (2 minutes).
#' @param datetime_col `character(1)`. Name of the datetime column.
#'   Default is `"datetime"`.
#'
#' @return A tibble with one row per detected issue and columns:
#'   \describe{
#'     \item{`row`}{`integer` — row index in `x` where the issue occurs.}
#'     \item{`datetime`}{`POSIXct` — timestamp at that row.}
#'     \item{`issue`}{`character` — one of `"gap"`, `"backward_jump"`, or
#'       `"year_artefact"`.}
#'     \item{`detail`}{`character` — human-readable description.}
#'   }
#'   Returns a zero-row tibble if no issues are found.
#'
#' @export
#'
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' rec    <- read_acttrust("recordings/P001.txt")
#' issues <- check_consistency(rec)
#' issues
#' }
check_consistency <- function(x, gap_s = 120, datetime_col = "datetime") {
  if (!datetime_col %in% names(x)) {
    zeitr_abort(
      "Column {.val {datetime_col}} not found in {.arg x}."
    )
  }

  times  <- as.POSIXct(x[[datetime_col]])
  deltas <- as.numeric(diff(times), units = "secs")
  years  <- as.integer(format(times, "%Y"))

  # Gaps
  gap_idx <- which(!is.na(deltas) & deltas > gap_s) + 1L
  gap_rows <- if (length(gap_idx) > 0L)
    tibble::tibble(
      row      = gap_idx,
      datetime = times[gap_idx],
      issue    = "gap",
      detail   = sprintf("%.0f s gap before this epoch", deltas[gap_idx - 1L])
    )
  else NULL

  # Backward jumps
  bj_idx <- which(!is.na(deltas) & deltas < 0) + 1L
  bj_rows <- if (length(bj_idx) > 0L)
    tibble::tibble(
      row      = bj_idx,
      datetime = times[bj_idx],
      issue    = "backward_jump",
      detail   = sprintf("timestamp went back %.0f s", abs(deltas[bj_idx - 1L]))
    )
  else NULL

  # Year artefacts
  ya_idx <- which(!is.na(years) & years %in% c(1970L, 2000L))
  ya_rows <- if (length(ya_idx) > 0L)
    tibble::tibble(
      row      = ya_idx,
      datetime = times[ya_idx],
      issue    = "year_artefact",
      detail   = sprintf("suspicious year %d (likely firmware bug)", years[ya_idx])
    )
  else NULL

  out <- do.call(rbind, Filter(Negate(is.null), list(gap_rows, bj_rows, ya_rows)))

  if (is.null(out)) {
    return(tibble::tibble(
      row      = integer(),
      datetime = as.POSIXct(character()),
      issue    = character(),
      detail   = character()
    ))
  }

  out[order(out$row), ]
}
