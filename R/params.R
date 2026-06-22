#' ActTrust device parameter preset
#'
#' Returns a named list of all device- and study-specific parameter values used
#' by the zeitR pipeline when processing ActTrust recordings. Passing this list
#' (or a modified copy) to [run_pipeline()] via the `params` argument fully
#' controls which values each detector stage receives, without having to touch
#' the individual function calls.
#'
#' The list is organised into four named sections that map directly to the four
#' pipeline stages:
#'
#' \describe{
#'   \item{`offwrist`}{Parameters for [detect_offwrist_bimodal()]. These values
#'     were validated against the Condor circadiaBase pipeline using ActTrust
#'     hardware. The refinement stage (`.bimodal_refine_acttrust`) is
#'     device-specific and currently only validated for ActTrust.}
#'   \item{`sleep`}{Parameters for [detect_sleep_crespo()]. `sleep_quantile`
#'     (0.365) is the ActTrust CSPD wrapper value; the original Crespo (2012)
#'     algorithm uses 1/3.}
#'   \item{`nap`}{Parameters for [detect_naps_crespo()], passed as the `params`
#'     argument. Equivalent to `.cspd_nap_params()`.}
#'   \item{`waso`}{Parameters for [compute_waso()].}
#' }
#'
#' To adapt the pipeline for a different device, copy the list, modify the
#' relevant values, and pass it to [run_pipeline()]:
#'
#' ```r
#' p <- acttrust_params()
#' p$sleep$sleep_quantile <- 1/3   # use the original Crespo threshold
#' result <- run_pipeline("recording.txt", params = p)
#' ```
#'
#' @return A named list with elements `offwrist`, `sleep`, `nap`, and `waso`.
#'
#' @seealso [run_pipeline()] for the pipeline entry point.
#'
#' @export
acttrust_params <- function() {
  list(

    # ── Off-wrist detection ─────────────────────────────────────────────────
    # detect_offwrist_bimodal() parameters. Validated against the Condor
    # circadiaBase pipeline using ActTrust hardware. Note: the border
    # refinement stage (.bimodal_refine_acttrust) is device-specific.
    offwrist = list(
      hws                 = 10L,
      activity_quantile   = 0.15,
      min_norm_activity   = 0.015,
      nbins               = 100L,
      min_offwrist_length = 10L,
      min_temp_threshold  = 0.35
    ),

    # ── Main sleep detection ────────────────────────────────────────────────
    # detect_sleep_crespo() parameters. sleep_quantile = 0.365 is the
    # ActTrust cspd_wrapper value; the original Crespo (2012) algorithm uses
    # 1/3. All other values are CSPD defaults.
    sleep = list(
      median_filter_h      = 8,
      pad_h                = 1,
      sleep_quantile       = 0.365,
      morph_size           = 61L,
      consec_zeros_thr     = 15L,
      awake_zeros_thr      = 2L,
      sleep_zeros_thr      = 30L,
      zero_mitigation_q    = 0.33,
      min_short_window_thr = 1.0
    ),

    # ── Nap detection ───────────────────────────────────────────────────────
    # detect_naps_crespo() params list. This is the full .cspd_nap_params()
    # set, reproduced here so it is visible and modifiable in one place.
    nap = .cspd_nap_params(),

    # ── WASO ────────────────────────────────────────────────────────────────
    # compute_waso() parameters.
    waso = list(
      wake_thresh = 60L
    )

  )
}
