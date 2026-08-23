#' Cosinor rhythmometry (Cornelissen 2014)
#'
#' Fits a single-harmonic cosine model to an activity/light/temperature
#' signal with a FIXED period (default 24h), returning the acrophase (peak
#' time), MESOR (rhythm-adjusted mean), and amplitude.
#'
#' @details
#' Ports pyActigraphy's actual `Cosinor` class (`pyActigraphy/analysis/
#' cosinor.py`) as it is really used in the reference notebook (Cell 16 of
#' `vs_condor_py_pipeline_fix30_jrsv.ipynb`'s `_fit_cosinor()`), not
#' reconstructed from documentation. The model is
#' \deqn{y = M + A \cos(\omega x + \phi)}
#' with \eqn{\omega = 2\pi / T}, \eqn{x} the integer epoch position
#' (`0, 1, 2, ...`, matching `Cosinor._convert_timestamp_to_index()`'s
#' `(index - index[0]) / index.freq`, not clock time or elapsed seconds),
#' \eqn{M} the MESOR, \eqn{A} the amplitude, \eqn{\phi} the acrophase.
#'
#' Cell 16's wrapper explicitly locks `Period` (`vary=False`) rather than
#' fitting it -- its own comment explains why: leaving `Period` free let
#' the optimizer converge anywhere from 0.23h to 91.75h, making the
#' acrophase/period meaningless. With the period fixed, the model is
#' linear in disguise: expanding
#' \eqn{A\cos(\omega x+\phi) = A\cos\phi\cos(\omega x) - A\sin\phi\sin(\omega x)}
#' and substituting \eqn{\beta_1 = A\cos\phi}, \eqn{\beta_2 = -A\sin\phi}
#' gives \eqn{y = M + \beta_1\cos(\omega x) + \beta_2\sin(\omega x)}, an
#' ordinary linear regression in \eqn{(M, \beta_1, \beta_2)} -- solved here
#' via `qr.solve()`, not a nonlinear optimizer. Verified by direct
#' execution that this closed-form OLS solution matches pyActigraphy's
#' actual `lmfit` (Levenberg-Marquardt) fit to within floating-point
#' precision when the period is locked, which it always is in the
#' notebook's actual usage -- this is not an approximation.
#'
#' `acrophase_time`/`acrophase_time_neg` reproduce the notebook's own
#' sign-flip fix exactly: the model's peak occurs at
#' `t_peak = -phi / omega = -phi * T / (2*pi)`, minutes since the first
#' epoch -- the notebook's own comment documents a prior bug where the
#' unnegated `+phi * T/(2*pi)` gave the ANTI-peak (trough), offset ~12h
#' from the true peak. `t_peak` is then wrapped to a 24h clock time via
#' `round(t_peak) %% 1440` -- hardcoded to 1440 regardless of the actual
#' `period_min` used for the fit, matching the source exactly (the
#' wrap-to-clock-time step assumes a 24h day, since it's meant to report a
#' wall-clock peak time, not a literal fit-period position).
#' `acrophase_time_neg` maps this to `[-12, 12]` (values `>= 12:00` shifted
#' by `-24`), matching the `_om10_neg`/`_ol5_neg` convention used elsewhere
#' in the same notebook for circular-safe comparison of clock times.
#'
#' @param x A `zeitr_recording`/`zeitr_result`, or a data frame / tibble
#'   with at least a `datetime` column and the column named in `col`.
#' @param col `character(1)`. Name of the signal column to fit (e.g.
#'   `"activity"`, `"light"`, `"int_temp"`). Default `"activity"`.
#' @param period_min `numeric(1)`. Fixed period in minutes. Default `1440`
#'   (24h), matching the reference notebook's locked-period usage.
#'
#' @return A tibble with columns `participant_id`, `acrophase_time`
#'   (`"HH:MM"` string), `acrophase_time_neg` (numeric, `[-12, 12]`),
#'   `MESOR`, `amplitude`, `period_min`, `n_epochs`.
#'
#' @references
#' Cornelissen, G. (2014). Cosinor-based rhythmometry. *Theoretical
#' Biology and Medical Modelling*, 11(1), 16.
#' \doi{10.1186/1742-4682-11-16}
#'
#' @seealso [compute_npcra()], [compute_sri()]
#'
#' @export
#'
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' result <- run_pipeline_native("recordings/P001.txt", tz = "America/Sao_Paulo")
#' compute_cosinor(result)
#' compute_cosinor(result, col = "light")
#' }
compute_cosinor <- function(x, col = "activity", period_min = 1440) {

  # ── Extract epochs tibble and participant_id ──────────────────────────────
  if (inherits(x, "zeitr_result")) {
    epochs         <- x$data
    participant_id <- x$subject_id %||% NA_character_
  } else if (inherits(x, "zeitr_recording")) {
    epochs         <- x$epochs
    participant_id <- x$metadata$participant_id %||% NA_character_
  } else if (is.data.frame(x)) {
    epochs         <- x
    participant_id <- NA_character_
  } else {
    zeitr_abort("{.arg x} must be a {.cls zeitr_result}, {.cls zeitr_recording}, or a data frame.")
  }

  required <- c("datetime", col)
  missing  <- setdiff(required, names(epochs))
  if (length(missing) > 0L) {
    zeitr_abort("Missing required column(s): {.val {missing}}")
  }

  datetimes <- as.POSIXct(epochs$datetime)
  signal    <- as.double(epochs[[col]])

  ord <- order(datetimes)
  signal <- signal[ord]

  keep   <- !is.na(signal)
  signal <- signal[keep]
  n      <- length(signal)

  if (n < 3L) {
    zeitr_warn("Fewer than 3 valid (non-NA) epochs; {.fn compute_cosinor} returns NA.")
    return(tibble::tibble(
      participant_id     = participant_id,
      acrophase_time     = NA_character_,
      acrophase_time_neg = NA_real_,
      MESOR              = NA_real_,
      amplitude          = NA_real_,
      period_min         = period_min,
      n_epochs           = n
    ))
  }

  fit <- .cosinor_fit(signal, period_min = period_min)

  tibble::tibble(
    participant_id     = participant_id,
    acrophase_time     = fit$acrophase_time,
    acrophase_time_neg = fit$acrophase_time_neg,
    MESOR              = fit$MESOR,
    amplitude          = fit$amplitude,
    period_min         = fit$period_min,
    n_epochs           = n
  )
}

#' Cosinor fit with a locked period, via closed-form OLS
#'
#' Ports pyActigraphy's actual `Cosinor` class as used in the reference
#' notebook (period always locked -- see `compute_cosinor()`'s Details for
#' the full derivation and verification against real `lmfit` execution).
#'
#' @param x numeric vector, the signal to fit (already `NA`-free).
#' @param period_min `numeric(1)`. Fixed period in minutes.
#' @return A list with `acrophase_time`, `acrophase_time_neg`, `MESOR`,
#'   `amplitude`, `period_min`.
#' @noRd
.cosinor_fit <- function(x, period_min) {
  n     <- length(x)
  t_idx <- 0:(n - 1L)
  omega <- 2 * pi / period_min

  X    <- cbind(1, cos(omega * t_idx), sin(omega * t_idx))
  beta <- qr.solve(X, x)

  mesor <- beta[1L]
  b1    <- beta[2L]
  b2    <- beta[3L]
  amplitude <- sqrt(b1^2 + b2^2)
  acrophase <- atan2(-b2, b1)

  # Peak time (minutes since the first epoch): t_peak = -phi / omega =
  # -phi * period_min / (2*pi). Wrapped to a 24h clock time via %% 1440,
  # hardcoded regardless of period_min, matching the source exactly.
  acro_min  <- -acrophase * period_min / (2 * pi)
  total_min <- round(acro_min) %% 1440
  hh <- total_min %/% 60L
  mm <- total_min %% 60L
  acrophase_time <- sprintf("%02d:%02d", hh, mm)

  acro_val <- hh + mm / 60
  acrophase_time_neg <- if (acro_val >= 12) acro_val - 24 else acro_val

  list(
    acrophase_time     = acrophase_time,
    acrophase_time_neg = acrophase_time_neg,
    MESOR              = mesor,
    amplitude          = amplitude,
    period_min         = period_min
  )
}
