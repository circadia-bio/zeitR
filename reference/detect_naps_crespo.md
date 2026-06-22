# Detect secondary sleep periods (naps) using the Crespo nap algorithm

Faithful port of the Python `nap_wrapper`: runs the full CSPD model in
nap mode (`detect_naps = TRUE`) on the currently-awake epochs
(`state == 0`) and merges the detected naps into the sleep state. Must
be run *after*
[`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md).

## Usage

``` r
detect_naps_crespo(x, epoch_h = NULL, params = .cspd_nap_params())
```

## Arguments

- x:

  A tibble as returned by
  [`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md),
  containing columns `datetime`, `activity`, and `state`. Nap detection
  runs on the wake (`state == 0`) subset only, mirroring the Python
  `nap_bool` mask.

- epoch_h:

  `numeric(1)`. Number of epochs per hour. If `NULL` (default), derived
  from the wake-subsequence epoch duration as `3600 / duration`.

- params:

  CSPD nap configuration list (default `.cspd_nap_params()`), the port
  of `nap_wrapper`'s parameter set.

## Value

The input tibble `x` with `state` and `sleep` columns updated. Nap
epochs become `state == 1` and `sleep == 1`; off-wrist (`state == 4`)
and existing main-sleep epochs are preserved.

## Details

Nap detection uses a nap-mode MSP (a high zero-proportion combined with
a low adaptive-median activity, `.crespo_nap_msp()`) followed by the
same bed-time / get-up-time refiners as the main sleep detection, with
the nap parameter set (`.cspd_nap_params()`) and nap-specific
minimum-length post-processing. Detected naps are written as
`state == 1` (merged into "sleep"), matching `nap_wrapper`, which
assigns `state[wake] = 1 - refined_output` (i.e. naps are *not* a
distinct state).

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
