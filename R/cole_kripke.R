#' Score actigraphy epochs as wake or sleep using the Cole-Kripke algorithm
#'
#' Applies the Cole-Kripke algorithm to a vector of zero-crossing mode (ZCM)
#' activity counts, scoring each epoch as wake (`1`) or sleep (`0`) using a
#' weighted sum of activity in a surrounding window.
#'
#' Each epoch's score is computed as:
#'
#' \deqn{D_i = P \sum_{j=1}^{9} W_j^{-} \cdot A_{i-j} +
#'             P \sum_{j=1}^{8} W_j^{+} \cdot A_{i+j}}
#'
#' where \eqn{A_i} is the ZCM count at epoch \eqn{i}, \eqn{W^{-}} and
#' \eqn{W^{+}} are the before and after weight vectors from Cole et al. (1992),
#' and \eqn{P = 0.000464}. Epochs with \eqn{D_i \ge 1} are scored as wake.
#'
#' @param zcm `numeric` vector of ZCM activity counts, one value per epoch.
#' @param P `numeric(1)`. Scaling factor. Default is `0.000464` (Cole et al.,
#'   1992).
#' @param weights_before `numeric(9)`. Weights applied to the 9 epochs
#'   *before* the current epoch. Defaults to the values from Cole et al.
#'   (1992), Table 2.
#' @param weights_after `numeric(8)`. Weights applied to the 8 epochs *after*
#'   the current epoch. Defaults to the values from Cole et al. (1992),
#'   Table 2.
#'
#' @return An integer vector the same length as `zcm`, with `1` indicating
#'   wake and `0` indicating sleep.
#'
#' @references
#' Cole, R. J., Kripke, D. F., Gruen, W., Mullaney, D. J., & Gillin, J. C.
#' (1992). Automatic sleep/wake identification from wrist activity.
#' *Sleep*, 15(5), 461–469. \doi{10.1093/sleep/15.5.461}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rec    <- read_acttrust("recordings/P001.txt")
#' scores <- score_epochs_cole_kripke(rec$ZCMn)
#' table(scores)  # 0 = sleep, 1 = wake
#' }
score_epochs_cole_kripke <- function(
    zcm,
    P              = 0.000464,
    weights_before = c(34.5, 133.0, 529.0, 375.0, 408.0,
                       400.5, 1074.0, 2048.5, 2424.5),
    weights_after  = c(1920.0, 149.5, 257.5, 125.0,
                       111.5,  120.0,  69.0,  40.5)
) {
  zcm <- as.double(zcm)
  n   <- length(zcm)
  nb  <- length(weights_before)
  na_ <- length(weights_after)

  if (n < 2L) {
    zeitr_warn("ZCM vector has fewer than 2 epochs; returning all zeros.")
    return(integer(n))
  }

  scores <- numeric(n)

  # Contributions from epochs *before* the current epoch
  for (i in seq_along(weights_before)) {
    offset <- nb - i + 1L
    idx_to   <- seq(offset + 1L, n)
    idx_from <- seq(1L, n - offset)
    scores[idx_to] <- scores[idx_to] + weights_before[i] * zcm[idx_from]
  }

  # Contributions from epochs *after* the current epoch
  for (i in seq_along(weights_after)) {
    offset <- i
    idx_to   <- seq(1L, n - offset)
    idx_from <- seq(offset + 1L, n)
    scores[idx_to] <- scores[idx_to] + weights_after[i] * zcm[idx_from]
  }

  scores <- scores * P
  as.integer(scores >= 1.0)
}
