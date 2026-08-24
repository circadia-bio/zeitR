# Cosinor rhythmometry (Cornelissen 2014)

Fits a single-harmonic cosine model to an activity/light/temperature
signal with a FIXED period (default 24h), returning the acrophase (peak
time), MESOR (rhythm-adjusted mean), and amplitude.

## Usage

``` r
compute_cosinor(x, col = "activity", period_min = 1440)
```

## Arguments

- x:

  A `zeitr_recording`/`zeitr_result`, or a data frame / tibble with at
  least a `datetime` column and the column named in `col`.

- col:

  `character(1)`. Name of the signal column to fit (e.g. `"activity"`,
  `"light"`, `"int_temp"`). Default `"activity"`.

- period_min:

  `numeric(1)`. Fixed period in minutes. Default `1440` (24h), matching
  the reference notebook's locked-period usage.

## Value

A tibble with columns `participant_id`, `acrophase_time` (`"HH:MM"`
string), `acrophase_time_neg` (numeric, `[-12, 12]`), `MESOR`,
`amplitude`, `period_min`, `n_epochs`.

## Details

Ports pyActigraphy's actual `Cosinor` class
(`pyActigraphy/analysis/ cosinor.py`) as it is really used in the
reference notebook (Cell 16 of
`vs_condor_py_pipeline_fix30_jrsv.ipynb`'s `_fit_cosinor()`), not
reconstructed from documentation. The model is \$\$y = M + A \cos(\omega
x + \phi)\$\$ with \\\omega = 2\pi / T\\, \\x\\ the integer epoch
position (`0, 1, 2, ...`, matching
`Cosinor._convert_timestamp_to_index()`'s
`(index - index[0]) / index.freq`, not clock time or elapsed seconds),
\\M\\ the MESOR, \\A\\ the amplitude, \\\phi\\ the acrophase.

Cell 16's wrapper explicitly locks `Period` (`vary=False`) rather than
fitting it – its own comment explains why: leaving `Period` free let the
optimizer converge anywhere from 0.23h to 91.75h, making the
acrophase/period meaningless. With the period fixed, the model is linear
in disguise: expanding \\A\cos(\omega x+\phi) = A\cos\phi\cos(\omega
x) - A\sin\phi\sin(\omega x)\\ and substituting \\\beta_1 = A\cos\phi\\,
\\\beta_2 = -A\sin\phi\\ gives \\y = M + \beta_1\cos(\omega x) +
\beta_2\sin(\omega x)\\, an ordinary linear regression in \\(M, \beta_1,
\beta_2)\\ – solved here via
[`qr.solve()`](https://rdrr.io/r/base/qr.html), not a nonlinear
optimizer. Verified by direct execution that this closed-form OLS
solution matches pyActigraphy's actual `lmfit` (Levenberg-Marquardt) fit
to within floating-point precision when the period is locked, which it
always is in the notebook's actual usage – this is not an approximation.

`acrophase_time`/`acrophase_time_neg` reproduce the notebook's own
sign-flip fix exactly: the model's peak occurs at
`t_peak = -phi / omega = -phi * T / (2*pi)`, minutes since the first
epoch – the notebook's own comment documents a prior bug where the
unnegated `+phi * T/(2*pi)` gave the ANTI-peak (trough), offset ~12h
from the true peak. `t_peak` is then wrapped to a 24h clock time via
`round(t_peak) %% 1440` – hardcoded to 1440 regardless of the actual
`period_min` used for the fit, matching the source exactly (the
wrap-to-clock-time step assumes a 24h day, since it's meant to report a
wall-clock peak time, not a literal fit-period position).
`acrophase_time_neg` maps this to `[-12, 12]` (values `>= 12:00` shifted
by `-24`), matching the `_om10_neg`/`_ol5_neg` convention used elsewhere
in the same notebook for circular-safe comparison of clock times.

## References

Cornelissen, G. (2014). Cosinor-based rhythmometry. *Theoretical Biology
and Medical Modelling*, 11(1), 16.
[doi:10.1186/1742-4682-11-16](https://doi.org/10.1186/1742-4682-11-16)

## See also

[`compute_npcra()`](https://zeitr.circadia-lab.uk/reference/compute_npcra.md),
[`compute_sri()`](https://zeitr.circadia-lab.uk/reference/compute_sri.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_pipeline_native("recordings/P001.txt", tz = "America/Sao_Paulo")
compute_cosinor(result)
compute_cosinor(result, col = "light")
} # }
```
