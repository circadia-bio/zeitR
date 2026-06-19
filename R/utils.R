# ── Internal utilities ────────────────────────────────────────────────────────
# Shared helper functions used across zeitR.
# Not exported.

#' Abort with a zeitR-prefixed message
#' @noRd
zeitr_abort <- function(msg, ...) {
  cli::cli_abort(msg, ...)
}

#' Warn with a zeitR-prefixed message
#' @noRd
zeitr_warn <- function(msg, ...) {
  cli::cli_warn(msg, ...)
}

#' Inform with a zeitR-prefixed message
#' @noRd
zeitr_inform <- function(msg, ...) {
  cli::cli_inform(msg, ...)
}
