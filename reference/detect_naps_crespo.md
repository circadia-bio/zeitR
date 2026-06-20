# Detect secondary sleep periods (naps) using the Crespo algorithm

Identifies naps — secondary sleep periods — in an actigraphy recording
using the nap variant of the Crespo algorithm. This function should be
run *after*
[`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md).

## Usage

``` r
detect_naps_crespo(
  x,
  epoch_h = NULL,
  median_filter_h = 8,
  pad_h = 1,
  nap_median_thr = 2,
  nap_zero_prop_thr = 0.5,
  nap_zero_prop_hws = 5L,
  use_and = FALSE
)
```

## Arguments

- x:

  A tibble as returned by
  [`detect_offwrist_bimodal()`](https://zeitr.circadia-lab.uk/reference/detect_offwrist_bimodal.md)
  (or
  [`prepare_actigraphy()`](https://zeitr.circadia-lab.uk/reference/prepare_actigraphy.md)
  if off-wrist detection is skipped), containing columns `datetime`,
  `activity`, and `state`.

- epoch_h:

  `numeric(1)`. Number of epochs per hour. If `NULL` (default),
  estimated automatically from the median inter-epoch interval in
  `datetime`.

- median_filter_h:

  `numeric(1)`. Length of the preprocessing median filter window in
  hours. Default is `8`.

- pad_h:

  `numeric(1)`. Padding length in hours added before the adaptive median
  filter. Default is `1`.

- nap_median_thr:

  `numeric(1)`. Epochs with a rolling median activity below this value
  may be scored as nap sleep. Default is `2.0`.

- nap_zero_prop_thr:

  `numeric(1)`. Epochs with a rolling zero-activity proportion above
  this threshold may be scored as nap sleep. Default is `0.5`.

- nap_zero_prop_hws:

  `integer(1)`. Half-window size (epochs) for the rolling
  zero-proportion filter. Default is `5L`.

- use_and:

  `logical(1)`. If `TRUE`, both the median activity AND zero-proportion
  criteria must be met. If `FALSE` (default), either criterion is
  sufficient.

## Value

The input tibble `x` with `state` and `sleep` columns updated. Nap
epochs have `state == 7` and are merged into `sleep` as value `1`.

## Details

The nap algorithm combines two criteria: a low rolling median activity
threshold and a high zero-activity proportion around each epoch. Epochs
satisfying either (or both, if `use_and = TRUE`) criteria are scored as
nap sleep.

## References

Crespo, C., Aboy, M., Fernández, J. R., & Mojón, A. (2012). Automatic
identification of activity-rest periods based on actigraphy. *Journal of
Medical and Biological Engineering*, 32(4), 249–256.
[doi:10.5405/jmbe.1033](https://doi.org/10.5405/jmbe.1033)

## See also

[`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md)
for main sleep period detection.

## Examples

``` r
if (FALSE) { # \dontrun{
rec  <- read_acttrust("recordings/P001.txt")
prep <- prepare_actigraphy(rec)
prep <- detect_offwrist_bimodal(prep)
prep <- detect_sleep_crespo(prep)
prep <- detect_naps_crespo(prep)
} # }
```
