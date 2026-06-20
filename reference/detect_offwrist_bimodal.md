# Off-wrist detection using the Condor bimodal activity/temperature model

Detects periods where the actigraph was not worn using the bimodal
algorithm developed by Condor Instruments. The algorithm proceeds in
three stages:

## Usage

``` r
detect_offwrist_bimodal(
  x,
  hws = 10L,
  activity_quantile = 0.15,
  min_norm_activity = 0.015,
  nbins = 100L,
  min_offwrist_length = 10L,
  min_temp_threshold = 0.35
)
```

## Arguments

- x:

  A tibble as returned by
  [`prepare_actigraphy()`](https://zeitr.circadia-lab.uk/reference/prepare_actigraphy.md),
  containing columns `datetime`, `activity` (PIM), `int_temp`,
  `ext_temp`, `state`, and `offwrist`.

- hws:

  `integer(1)`. Half-window size (in epochs) for rolling feature
  extraction. Default is `10` (matching the Python pipeline).

- activity_quantile:

  `numeric(1)`. Quantile used to define "low activity". Default is
  `0.15`.

- min_norm_activity:

  `numeric(1)`. Minimum normalised activity threshold below which the
  low-activity cutoff is clamped. Default is `0.015`.

- nbins:

  `integer(1)`. Number of histogram bins used when fitting the GMM
  threshold. Default is `100`.

- min_offwrist_length:

  `integer(1)`. Minimum number of consecutive off-wrist epochs required
  to retain a detected period. Shorter runs are discarded. Default is
  `10`.

- min_temp_threshold:

  `numeric(1)`. Minimum normalised temperature threshold; the fitted GMM
  threshold is clamped to this value if it falls below it. Default is
  `0.35`.

## Value

The input tibble `x` with `state` and `offwrist` columns updated.
Off-wrist epochs have `state == 4` and `offwrist == 0.25`.

## Details

1.  **Feature extraction** — rolling median of PIM activity and rolling
    signal/variance of internal temperature are computed over a
    symmetric window of half-width `hws` epochs.

2.  **Bimodal temperature threshold** — among epochs with low activity
    median (below the `activity_quantile` quantile of normalised
    activity), a 2-component Gaussian Mixture Model is fitted to the
    normalised temperature distribution. The threshold between the two
    components is taken as the minimum of the GMM density between the
    two means (or the minimum of the smoothed histogram if that is
    lower). Ashman's D is computed as a bimodality quality metric.

3.  **Initial classification & refinement** — epochs simultaneously
    exhibiting low activity median AND low temperature are marked as
    off-wrist. Short spurious off-wrist runs shorter than
    `min_offwrist_length` epochs are removed.

Off-wrist epochs are encoded as `state == 4` (matching the Python
pipeline convention). The `offwrist` column is set to `0.25` for
off-wrist epochs for actogram overlay plotting.

## References

The bimodal off-wrist algorithm was developed by Julius A. P. P. de
Paula at Condor Instruments (2023). It is not published in peer-reviewed
literature but the source code is available in the circadiaBase pipeline
repository. The Ashman D statistic is described in:

Ashman, K. M., Bird, C. M., & Zepf, S. E. (1994). Detecting bimodality
in astronomical datasets. *The Astronomical Journal*, 108, 2348.
[doi:10.1086/117248](https://doi.org/10.1086/117248)

## Examples

``` r
if (FALSE) { # \dontrun{
rec  <- read_acttrust("recordings/P001.txt")
prep <- prepare_actigraphy(rec)
prep <- detect_offwrist_bimodal(prep)
sum(prep$state == 4)  # number of off-wrist epochs
} # }
```
