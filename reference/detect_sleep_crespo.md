# Detect main sleep periods using the Crespo algorithm

Identifies the main sleep period(s) in an actigraphy recording using the
algorithm described in Crespo et al. (2012). The method applies an
adaptive median filter to the activity signal, mitigates spuriously long
zero runs, and thresholds the result at a quantile of the filtered
signal. Morphological closing and opening operations are then used to
smooth the binary sleep/wake estimate.

## Usage

``` r
detect_sleep_crespo(
  x,
  epoch_h = NULL,
  median_filter_h = 8,
  pad_h = 1,
  sleep_quantile = 0.365,
  morph_size = 61L,
  consec_zeros_thr = 15L,
  awake_zeros_thr = 2L,
  sleep_zeros_thr = 30L,
  zero_mitigation_q = 0.33,
  min_short_window_thr = 1,
  refine = TRUE,
  condition = 0L
)
```

## Arguments

- x:

  A tibble as returned by
  [`detect_offwrist_bimodal()`](https://zeitr.circadia-lab.uk/reference/detect_offwrist_bimodal.md)
  (or
  [`prepare_actigraphy()`](https://zeitr.circadia-lab.uk/reference/prepare_actigraphy.md)
  if off-wrist detection is skipped), containing columns `datetime`,
  `activity`, and `state`. The detector runs on the on-wrist subset
  (`state != 4`) only, mirroring the Python `cspd_wrapper`.

- epoch_h:

  `numeric(1)`. Number of epochs per hour. If `NULL` (default), derived
  from the epoch duration (the mode of the on-wrist inter-epoch
  interval), as `3600 / duration`.

- median_filter_h:

  `numeric(1)`. Length of the preprocessing median filter window in
  hours. Default is `8`.

- pad_h:

  `numeric(1)`. Padding length in hours added before the adaptive median
  filter. Default is `1`.

- sleep_quantile:

  `numeric(1)`. Quantile of the filtered activity used as the MSP
  sleep/wake threshold. Default is `0.365` (the ActTrust CSPD value used
  by `cspd_wrapper`; the standalone Crespo algorithm uses `1/3`).

- morph_size:

  `integer(1)`. Size of the structuring element used in morphological
  closing/opening. Default is `61` epochs.

- consec_zeros_thr:

  `integer(1)`. Runs of zeros longer than this threshold are treated as
  invalid (zero mitigation). Default is `15`.

- awake_zeros_thr:

  `integer(1)`. Threshold for consecutive zeros within wake periods.
  Default is `2`.

- sleep_zeros_thr:

  `integer(1)`. Threshold for consecutive zeros within sleep periods.
  Default is `30`.

- zero_mitigation_q:

  `numeric(1)`. Quantile of activity used to determine the mitigation
  level for invalid zero runs. Default is `0.33`.

- min_short_window_thr:

  `numeric(1)`. Minimum value of the adaptive median threshold; if the
  fitted quantile falls below this, the threshold is clamped here.
  Default is `1.0`.

- refine:

  `logical(1)`. If `TRUE` (default), the MSP detection is refined into
  final sleep periods by the CSPD bed-time / get-up-time refiners
  (`.cspd_refine_periods`), reproducing the Python `refined_output`. If
  `FALSE`, the raw MSP detection is used directly.

- condition:

  `integer(1)`. Initial condition flag (default `0`). The MSP stage
  bumps it to `2` when its activity-median threshold clamps to
  `min_short_window_thr`; the refiner uses that effective condition.
  Affects only `refine = TRUE`.

## Value

The input tibble `x` with `state` and `sleep` columns updated. Sleep
epochs have `state == 1` and `sleep == 1`; off-wrist epochs
(`state == 4`) are preserved and excluded from the sleep column. With
`refine = TRUE` the sleep epochs delimit the refined sleep PERIODS (the
Python `refined_output`); per-epoch wake/sleep within them is scored
later by
[`compute_waso()`](https://zeitr.circadia-lab.uk/reference/compute_waso.md).

## References

Crespo, C., Aboy, M., Fernández, J. R., & Mojón, A. (2012). Automatic
identification of activity-rest periods based on actigraphy. *Journal of
Medical and Biological Engineering*, 32(4), 249–256.
[doi:10.5405/jmbe.1033](https://doi.org/10.5405/jmbe.1033)

## See also

[`detect_naps_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_naps_crespo.md)
for secondary sleep period (nap) detection.

## Examples

``` r
if (FALSE) { # \dontrun{
rec  <- read_acttrust("recordings/P001.txt")
prep <- prepare_actigraphy(rec)
prep <- detect_offwrist_bimodal(prep)
prep <- detect_sleep_crespo(prep)
} # }
```
