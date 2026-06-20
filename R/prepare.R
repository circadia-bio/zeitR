#' Prepare a raw actigraphy tibble for analysis
#'
#' Transforms the output of [read_acttrust()] (or any reader that returns a
#' compatible tibble) into the working data frame expected by all detection
#' functions. Specifically:
#'
#' 1. Clamps `int_temp` and `ext_temp` to the physiological range \[0, 42\] °C.
#' 2. Adds min-max scaled temperature columns `int_temp_` and `ext_temp_`
#'    (range \[0, 1\]) for plotting.
#' 3. Ensures `state`, `offwrist`, and `sleep` columns are present and
#'    initialised to `0`.
#'
#' The original tibble is never modified; a copy is returned.
#'
#' @param x A tibble as returned by [read_acttrust()], containing at minimum
#'   `int_temp` and `ext_temp` columns.
#'
#' @return A tibble of the same dimensions as `x` with additional or updated
#'   columns: `state`, `offwrist`, `sleep`, `int_temp_`, `ext_temp_`.
#'
#' @export
#'
#' @importFrom tidyr replace_na
#'
#' @examples
#' \dontrun{
#' rec  <- read_acttrust("recordings/P001.txt")
#' prep <- prepare_actigraphy(rec)
#' }
prepare_actigraphy <- function(x) {
  if (!inherits(x, "data.frame")) {
    zeitr_abort("{.arg x} must be a data frame or tibble.")
  }

  required <- c("int_temp", "ext_temp")
  missing  <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    zeitr_abort(
      "{.arg x} is missing required column(s): {.val {missing}}"
    )
  }

  x <- x

  # ── Temperature clamping ─────────────────────────────────────────────────────
  int_temp <- pmin(pmax(as.double(x$int_temp), 0), 42)
  ext_temp <- pmin(pmax(as.double(tidyr::replace_na(x$ext_temp, 0)), 0), 42)

  x$int_temp <- int_temp
  x$ext_temp <- ext_temp

  # ── Min-max scaling for plotting ─────────────────────────────────────────────
  scale_max <- max(max(int_temp, na.rm = TRUE), max(ext_temp, na.rm = TRUE))
  if (scale_max > 0) {
    x$int_temp_ <- int_temp / scale_max
    x$ext_temp_ <- ext_temp / scale_max
  } else {
    x$int_temp_ <- int_temp
    x$ext_temp_ <- ext_temp
  }

  # ── State columns ─────────────────────────────────────────────────────────────
  if (!"state"    %in% names(x)) x$state    <- 0
  if (!"offwrist" %in% names(x)) x$offwrist <- 0
  if (!"sleep"    %in% names(x)) x$sleep    <- 0

  x
}
