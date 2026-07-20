# ── Physical activity intensity from activity counts ─────────────────────────
# Ports the published MET-estimation equations and count-based cut-points
# from a single controlled-treadmill validation study:
#
#   Batista ES, Basilio Silva Gomes SR, Bruno de Morais Ferreira A, França LGS,
#   Fontenele Araújo J, Mortatti AL, Leocadio-Miguel MA (2026). From movement
#   to METs: A validation of ActTrust(R) for energy expenditure estimation and
#   physical activity classification in young adults. PLoS ONE 21(5):e0348631.
#   https://doi.org/10.1371/journal.pone.0348631
#   Data/code: https://github.com/circadia-bio/ACTT_validation_study
#
# What is (and isn't) ported:
# The original study fits sqrt(MET) ~ sqrt(activity) * device_placement via
# lm(), then derives cut-points and their 95% CIs with msm::deltamethod().
# That model needs indirect-calorimetry data (breath-by-breath VO2) that
# essentially no zeitR user will have. What's reproduced here is the
# *published* coefficient/cut-point table (their Table 2) -- fixed constants,
# applied at inference time to new ActTrust or GT3X+ counts. No caret,
# MuMIn, msm, or pROC dependency is introduced; those packages were needed to
# validate the equations, not to apply them.
#
# Equation (their Eq. 2): sqrt(MET) = b0 + b1 * sqrt(activity_counts)
#   => MET = (b0 + b1 * sqrt(activity_counts))^2

# ── Coefficient / cut-point table ─────────────────────────────────────────────

#' Published activity-count equations for estimating METs and PA intensity
#'
#' Returns the regression coefficients and count-based cut-points for
#' estimating energy expenditure (METs) and classifying physical activity (PA)
#' intensity from ActTrust(R) (Condor Instruments) or ActiGraph(R) GT3X+
#' activity counts, hip or wrist placement.
#'
#' @details
#' # Important caveats
#' These equations come from a **single** controlled-laboratory validation
#' study (Batista et al. 2026; N = 56 healthy adults aged 18-35; treadmill
#' walking/running at 3-9 km/h only). Treat them as one available equation
#' set, not a universal standard:
#'
#' * The paper's own Discussion section compares its GT3X+ (hip) cut-points
#'   against two other published GT3X+ studies (Sasaki et al. 2011;
#'   Santos-Lozano et al. 2013) and finds differences of 2-65% depending on
#'   the MET threshold -- and those two reference studies differ from *each
#'   other* by 16-39%. The paper attributes that spread mainly to
#'   sample-level characteristics rather than device or methodological
#'   artefacts, but that reading applies to the comparison *between* Sasaki
#'   and Santos-Lozano -- it doesn't fully carry over to comparisons against
#'   this paper's own equations, since all three studies used different
#'   modelling approaches (Sasaki et al.: ActiGraph's two-regression model;
#'   Santos-Lozano et al.: an artificial neural network; this paper: a
#'   simple sqrt-transformed linear model with a device x placement
#'   interaction). Cross-study cut-point differences therefore reflect model
#'   choice as well as sample, not sample alone -- and that's checkable to
#'   different degrees: this paper's linear coefficients (`b0`/`b1` below)
#'   are transparent enough to compare term-by-term against Sasaki et al.'s
#'   two-regression model, so a discrepancy there is at least diagnosable.
#'   Santos-Lozano et al.'s cut-points come from an ANN, which has no
#'   inspectable coefficients -- there is no way to say *why* it disagrees
#'   with the equations here, only *that* it does.
#' * The ACTT (hip)/ACTT (wrist) equations are the first published cut-points
#'   for ActTrust(R) at all, so there is nothing yet to cross-check them
#'   against.
#' * Validated only for laboratory treadmill walking/running in healthy young
#'   adults. Applying these equations to free-living data, other age groups
#'   (children, older adults), clinical populations, or other activity types
#'   is an extrapolation the source study explicitly flags as untested.
#'
#' `equation_set` is exposed as an explicit argument (rather than hard-coding
#' a single table) so a future validation study covering a different
#' population or device can be added as an alternative set without changing
#' the [estimate_ee()] / [classify_pa_intensity()] API.
#'
#' @param equation_set `character(1)`. Currently only `"batista2026"` is
#'   available.
#'
#' @return A tibble with one row per device x placement combination and
#'   columns:
#'   \describe{
#'     \item{`device`}{`"ACTT"` or `"GT3X+"`.}
#'     \item{`placement`}{`"hip"` or `"wrist"`.}
#'     \item{`b0`, `b1`}{Intercept and slope of
#'       `sqrt(MET) = b0 + b1 * sqrt(activity_counts)`.}
#'     \item{`cut3`, `cut6`, `cut9`}{Published count/min cut-points
#'       (95% CI midpoints) at the 3, 6, and 9 MET thresholds -- i.e. the
#'       count value at which the fitted equation crosses that MET level.
#'       Provided for reference; [classify_pa_intensity()] classifies on
#'       estimated METs directly rather than re-deriving these.}
#'   }
#'
#' @seealso [estimate_ee()], [classify_pa_intensity()]
#'
#' @export
#'
#' @examples
#' pa_equations()
pa_equations <- function(equation_set = "batista2026") {

  if (!identical(equation_set, "batista2026")) {
    zeitr_abort(c(
      "Unknown {.arg equation_set}: {.val {equation_set}}.",
      "i" = 'Only {.val batista2026} is currently available.'
    ))
  }

  tibble::tibble(
    device    = c("GT3X+", "ACTT",  "ACTT",   "GT3X+"),
    placement = c("hip",   "hip",   "wrist",  "wrist"),
    b0        = c(1.062,   1.107,   1.233,    1.207),
    b1        = c(0.0199,  0.0088,  0.0081,   0.0127),
    cut3      = c(1132,    5057,    3761,     1698),
    cut6      = c(4853,    23339,   22368,    9503),
    cut9      = c(9468,    46410,   47203,    19787)
  )
}


# ── Energy expenditure estimation ─────────────────────────────────────────────

#' Estimate energy expenditure (METs) from activity counts
#'
#' Applies the published ActTrust(R)/GT3X+ regression equation
#' (`sqrt(MET) = b0 + b1 * sqrt(activity_counts)`) from
#' [pa_equations()] to convert raw activity counts into estimated metabolic
#' equivalents (METs).
#'
#' @inherit pa_equations details
#'
#' @param counts Numeric vector of activity counts (counts/min, over the same
#'   epoch length used by the source equation set -- 1-minute for
#'   `"batista2026"`). Negative values are treated as `0` with a warning, since
#'   the underlying square-root transform is undefined below zero.
#' @param device `character(1)`. `"ACTT"` (default) or `"GT3X+"`.
#' @param placement `character(1)`. `"hip"` (default) or `"wrist"`.
#' @param equation_set `character(1)`, passed to [pa_equations()]. Default
#'   `"batista2026"`.
#'
#' @return Numeric vector of estimated METs, same length as `counts`.
#'
#' @seealso [pa_equations()] for the underlying coefficients and their
#'   limitations, [classify_pa_intensity()] to convert METs into intensity
#'   bands.
#'
#' @export
#'
#' @examples
#' estimate_ee(c(0, 5000, 25000, 50000), device = "ACTT", placement = "hip")
estimate_ee <- function(counts,
                         device       = c("ACTT", "GT3X+"),
                         placement    = c("hip", "wrist"),
                         equation_set = "batista2026") {
  device    <- match.arg(device)
  placement <- match.arg(placement)

  if (any(counts < 0, na.rm = TRUE)) {
    zeitr_warn("{sum(counts < 0, na.rm = TRUE)} negative count value(s) treated as 0.")
    counts <- pmax(counts, 0)
  }

  eq  <- pa_equations(equation_set)
  row <- eq[eq$device == device & eq$placement == placement, ]

  if (nrow(row) == 0L) {
    zeitr_abort("No equation found for device = {.val {device}}, placement = {.val {placement}}.")
  }

  (row$b0 + row$b1 * sqrt(counts))^2
}


# ── PA intensity classification ───────────────────────────────────────────────

#' Classify physical activity intensity from METs
#'
#' Bins metabolic equivalent (MET) values into the four PA intensity classes
#' used throughout Batista et al. (2026) and consistent with WHO PA intensity
#' conventions: light, moderate, vigorous, very vigorous.
#'
#' @param mets Numeric vector of MET values, e.g. from [estimate_ee()] or
#'   measured directly (indirect calorimetry, published equations, etc.).
#'   MET-based, not device-specific -- any source of METs can be classified.
#'
#' @return An ordered factor with levels `c("light", "moderate", "vigorous",
#'   "very_vigorous")`, same length as `mets`. Bands: `[0,3)` light,
#'   `[3,6)` moderate, `[6,9)` vigorous, `[9,Inf)` very vigorous.
#'
#' @seealso [estimate_ee()], [pa_equations()]
#'
#' @export
#'
#' @examples
#' mets <- estimate_ee(c(0, 5000, 25000, 50000), device = "ACTT", placement = "hip")
#' classify_pa_intensity(mets)
classify_pa_intensity <- function(mets) {
  factor(
    cut(mets,
        breaks = c(-Inf, 3, 6, 9, Inf),
        labels = c("light", "moderate", "vigorous", "very_vigorous"),
        right  = FALSE),
    levels  = c("light", "moderate", "vigorous", "very_vigorous"),
    ordered = TRUE
  )
}


#' Classify physical activity intensity directly from activity counts
#'
#' Convenience wrapper combining [estimate_ee()] and [classify_pa_intensity()]
#' in one call: converts raw ActTrust(R)/GT3X+ activity counts straight to a PA
#' intensity factor via the published `"batista2026"` (or future) equation
#' set, without needing to handle the intermediate MET values.
#'
#' @inheritParams estimate_ee
#'
#' @return An ordered factor, same length as `counts` and same levels as
#'   [classify_pa_intensity()].
#'
#' @seealso [estimate_ee()], [classify_pa_intensity()], [pa_equations()] for
#'   the coefficients and their generalisability caveats.
#'
#' @export
#'
#' @examples
#' classify_pa_counts(c(0, 5000, 25000, 50000), device = "ACTT", placement = "hip")
classify_pa_counts <- function(counts,
                                device       = c("ACTT", "GT3X+"),
                                placement    = c("hip", "wrist"),
                                equation_set = "batista2026") {
  device    <- match.arg(device)
  placement <- match.arg(placement)
  classify_pa_intensity(
    estimate_ee(counts, device = device, placement = placement, equation_set = equation_set)
  )
}
