# ── Circular statistics for clock-time variables ─────────────────────────────
# Fix 20 (unit-circle circular mean) and Fix 24 (circular SD) from the
# JRSV condor pipeline (pipeline_functions_fix27.py).

#' Circular mean of clock times in decimal hours
#'
#' Computes the mean of times on the 0-24 h clock using the unit-circle
#' method: convert each time to a unit-vector angle, average the sin and
#' cos components, and back-project via `atan2`. This correctly handles
#' midnight wrap (e.g. times at 23:30 and 00:30 average to 00:00, not 12:00).
#'
#' The older threshold approach (shift values >= 12 h by -24 h before
#' averaging) produces badly biased estimates when times straddle noon --
#' critically for sleep offset in late sleepers -- and is superseded by
#' this function (Fix 20 in the JRSV pipeline, ICC for sleep offset
#' 0.502 -> 0.897 on N=404).
#'
#' @param x `numeric`. Decimal hours in \[0, 24). `NA` values are silently
#'   dropped.
#'
#' @return A single `numeric` in \[0, 24), or `NA_real_` if `x` is empty
#'   after `NA` removal.
#'
#' @references
#' Pewsey, A., Neuhauser, M., & Ruxton, G. D. (2013). *Circular Statistics
#' in R*. Oxford University Press.
#'
#' @export
#'
#' @examples
#' circ_mean_h(c(23.5, 0.5))   # midnight wrap -> 0
#' circ_mean_h(c(23.0, 1.0))   # -> 0
#' circ_mean_h(c(7.0, 9.0))    # -> 8
circ_mean_h <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_real_)
  theta <- x * (2.0 * pi / 24.0)
  ang   <- atan2(mean(sin(theta)), mean(cos(theta)))
  (ang * 24.0 / (2.0 * pi) + 24.0) %% 24.0
}

#' Circular standard deviation of clock times in decimal hours
#'
#' Computes the circular SD using the mean resultant length formula:
#' `sqrt(-2 * log(R_bar)) * 24 / (2*pi)`, where `R_bar` is the mean
#' resultant length of the unit-circle representation of the times.
#' Linear SD inflates when times straddle midnight; this measure is
#' invariant to the wrap point (Fix 24 in the JRSV pipeline).
#'
#' @param x `numeric`. Decimal hours in \[0, 24). `NA` values are silently
#'   dropped. Returns `NA_real_` for fewer than 2 values.
#'
#' @return A single non-negative `numeric` (hours), or `NA_real_`.
#'
#' @references
#' Mardia, K. V., & Jupp, P. E. (2000). *Directional Statistics* (2nd ed.).
#' Wiley.
#'
#' @export
#'
#' @examples
#' circ_sd_h(c(23.5, 0.5))   # small SD for times close to midnight
#' circ_sd_h(c(6.0, 18.0))   # large SD for times 12 h apart
circ_sd_h <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2L) return(NA_real_)
  theta <- x * (2.0 * pi / 24.0)
  R_bar <- sqrt(mean(cos(theta))^2 + mean(sin(theta))^2)
  R_bar <- min(1.0, max(0.0, R_bar))   # numerical clamp
  sqrt(-2.0 * log(R_bar)) * 24.0 / (2.0 * pi)
}
