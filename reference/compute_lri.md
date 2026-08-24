# Light Regularity Index (LRI)

Computes the Light Regularity Index (Hand et al. 2023): a measure of
day-to-day consistency in light exposure timing, analogous to the Sleep
Regularity Index but applied to a binarized light-exposure signal.

## Usage

``` r
compute_lri(x, col = "light", threshold = NULL, log_transform = FALSE)
```

## Arguments

- x:

  A `zeitr_recording`/`zeitr_result`, or a data frame / tibble with at
  least a `datetime` column and the column named in `col`.

- col:

  `character(1)`. Name of the light channel column. Default `"light"`.

- threshold:

  `numeric(1)` or `NULL`. Threshold applied to `col` (`col > threshold`
  = "light exposure"). Default `NULL`, which computes all five of the
  reference notebook's default thresholds (`10`, `20`, `50`, `100`,
  `300`, applied as `1`, `1.301`, `1.699`, `2`, `2.4771` – see Details
  on why these numbers don't mean lux against raw data) and returns one
  row with columns `LRI_10`...`LRI_300`. Pass a single number for one
  custom threshold instead (returns a single `LRI` column).

- log_transform:

  `logical(1)`. If `TRUE`, applies `log10(col + 1)` before thresholding
  – the probably-originally-intended behaviour implied by the reference
  notebook's (incorrect, per Details) comment. Default `FALSE`, matching
  what the reference actually computes.

## Value

A tibble with `participant_id`, either `LRI_10`, `LRI_20`, `LRI_50`,
`LRI_100`, `LRI_300` (default, `threshold = NULL`) or a single `LRI`
column (custom `threshold`), and `n_epochs`.

## Details

Ports the `artvalencio/pyActigraphy` fork's actual `LRI()` method
(`pyActigraphy/light/light_metrics.py`) exactly – found there after an
earlier investigation confirmed it does not exist anywhere in the
official `ghammad/pyActigraphy` repository (any branch, any commit in
history). `LRI()`'s own inner `prob_stability()`/`lri_profile()`
functions are structurally identical to
[`compute_sri()`](https://zeitr.circadia-lab.uk/reference/compute_sri.md)'s
`"sadeh"`/ `"ck"`/`"scripps"`/`"roenneberg"` paths' aggregation
(`.sri_pyactigraphy()`) – verified by direct execution that `LRI()`'s
own code and a direct call to `sri()`
(`pyActigraphy/sleep/scoring/ sri.py`) give IDENTICAL results on the
same real light data, to full floating-point precision. `compute_lri()`
reuses `.sri_pyactigraphy()` directly rather than duplicating the same
formula under a different name.

**A real discrepancy worth knowing about, not a guess**: the reference
notebook's own comment (Cell 16 of
`vs_condor_py_pipeline_fix30_jrsv.ipynb`) claims "pyActigraphy stores
light data in log10-transformed lux internally," and picks its default
thresholds accordingly
(`10/20/50/100/300 lux -> log10 -> 1.0/1.301/1.699/2.0/2.4771`). This is
not what the code actually does: read directly, `RawATR`'s constructor
(`pyActigraphy/io/atr/atr.py`) passes the `LIGHT` column straight into
`LightRecording` with no [`log10()`](https://rdrr.io/r/base/Log.html)
call anywhere, and real light values from actual ActTrust recordings
confirm this empirically – raw `LIGHT` ranges from 0 to over 20,000 in
real data checked this session, nowhere near the ~0-5 range `log10(lux)`
would produce. So the reference Python's own `LRI_10`/`LRI_20`/etc.
columns almost certainly don't measure what their names claim –
thresholds of ~1-2.5 applied to raw, untransformed light amount to "any
detectable light vs. darkness," not "50 lux vs. below." This function
replicates the reference's actual behaviour (raw threshold, no
transform) for parity purposes – matching what the real pipeline
computes, not what its comments say it computes. `log_transform = TRUE`
is available if you want the – probably originally intended – log10(lux)
behaviour instead; it is NOT the default, since the default should match
the real reference output.

## References

Hand, A. J., Lane, J. M., Payne, A. C., ... & Barker, D. H. (2023).
Measuring light regularity: sleep regularity is associated with
regularity of light exposure in adolescents. *Sleep*, 46(8), zsad001.
[doi:10.1093/sleep/zsad001](https://doi.org/10.1093/sleep/zsad001)

## See also

[`compute_sri()`](https://zeitr.circadia-lab.uk/reference/compute_sri.md),
[`compute_cosinor()`](https://zeitr.circadia-lab.uk/reference/compute_cosinor.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_pipeline_native("recordings/P001.txt", tz = "America/Sao_Paulo")
compute_lri(result)
compute_lri(result, threshold = 50)
} # }
```
