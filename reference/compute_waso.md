# Compute WASO and nightly sleep statistics

Scores each epoch within detected sleep periods as wake or sleep using
[`score_epochs_cole_kripke()`](https://zeitr.circadia-lab.uk/reference/score_epochs_cole_kripke.md),
then computes per-night statistics. This is a faithful port of the
Python `condor_pipeline` `detect_waso`: night boundaries are built by
`.nights_df()` (the `search_gap = FALSE` path of the reference
`nights_df`), each night's ZCM is scored with Cole-Kripke, and the
epoch-level `state` is rebuilt from a zero (wake) base so that
within-night WASO-wake epochs and all epochs outside detected nights are
scored as wake.

## Usage

``` r
compute_waso(x, wake_thresh = 60L, search_gap = FALSE)
```

## Arguments

- x:

  A tibble as returned by
  [`detect_naps_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_naps_crespo.md)
  (or
  [`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md)
  if nap detection is skipped), containing columns `datetime`, `ZCMn`,
  and `state`.

- wake_thresh:

  `integer(1)`. Minimum duration in epochs of a wake bout required to
  delimit a new sleep period boundary. Default is `60`.

- search_gap:

  `logical(1)`. Reserved for API compatibility with the reference
  pipeline. The reference `detect_waso` always builds boundaries with
  `search_gap = FALSE`, so this argument currently has no effect.

## Value

A list with two elements:

- `nights`:

  A tibble with one row per detected night/nap and columns `night`,
  `is_nap`, `bed_time`, `get_up_time`, `tbt`, `tst`, `waso`, `sol`,
  `soi`, `nw`, `eff`.

- `data`:

  The input tibble `x` with `state` and `sleep` updated with epoch-level
  Cole-Kripke wake/sleep scores.

## Details

The following metrics are computed for each detected night / nap:

|  |  |
|----|----|
| Metric | Definition |
| `tbt` | Total Bed Time — epochs from bed time to get-up time |
| `tst` | Total Sleep Time — `tbt - waso - sol - soi` |
| `waso` | Wake After Sleep Onset — wake epochs between sleep onset and final wake |
| `sol` | Sleep Onset Latency — epochs from bed time to first sleep epoch |
| `soi` | Sleep Offset Inertia — trailing wake epochs at end of sleep period |
| `nw` | Number of awakenings — count of wake-onset transitions |
| `eff` | Sleep efficiency — `tst / tbt` |

## References

Cole, R. J., Kripke, D. F., Gruen, W., Mullaney, D. J., & Gillin, J. C.
(1992). Automatic sleep/wake identification from wrist activity.
*Sleep*, 15(5), 461–469.
[doi:10.1093/sleep/15.5.461](https://doi.org/10.1093/sleep/15.5.461)

## Examples

``` r
if (FALSE) { # \dontrun{
rec    <- read_acttrust("recordings/P001.txt")
prep   <- prepare_actigraphy(rec)
prep   <- detect_offwrist_bimodal(prep)
prep   <- detect_sleep_crespo(prep)
prep   <- detect_naps_crespo(prep)
result <- compute_waso(prep)
result$nights
} # }
```
