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
  n      <- length(times)
  issues <- vector("list", n)
  k      <- 0L

  deltas <- as.numeric(diff(times), units = "secs")

  for (i in seq_len(n - 1L)) {
    d <- deltas[i]

    if (!is.na(d) && d > gap_s) {
      k <- k + 1L
      issues[[k]] <- list(
        row      = i + 1L,
        datetime = times[i + 1L],
        issue    = "gap",
        detail   = sprintf("%.0f s gap before this epoch", d)
      )
    }

    if (!is.na(d) && d < 0) {
      k <- k + 1L
      issues[[k]] <- list(
        row      = i + 1L,
        datetime = times[i + 1L],
        issue    = "backward_jump",
        detail   = sprintf("timestamp went back %.0f s", abs(d))
      )
    }
  }

  # Year artefacts
  years <- as.integer(format(times, "%Y"))
  for (i in seq_len(n)) {
    if (!is.na(years[i]) && years[i] %in% c(1970L, 2000L)) {
      k <- k + 1L
      issues[[k]] <- list(
        row      = i,
        datetime = times[i],
        issue    = "year_artefact",
        detail   = sprintf("suspicious year %d (likely firmware bug)", years[i])
      )
    }
  }

  if (k == 0L) {
    return(tibble::tibble(
      row      = integer(),
      datetime = as.POSIXct(character()),
      issue    = character(),
      detail   = character()
    ))
  }

  issues <- issues[seq_len(k)]
  out    <- tibble::tibble(
    row      = vapply(issues, `[[`, integer(1),  "row"),
    datetime = do.call(c, lapply(issues, `[[`, "datetime")),
    issue    = vapply(issues, `[[`, character(1), "issue"),
    detail   = vapply(issues, `[[`, character(1), "detail")
  )

  out[order(out$row), ]
}
