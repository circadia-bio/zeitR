# Sleep Regularity Index (SRI)

Computes the Sleep Regularity Index (Phillips et al. 2017): a measure of
day-to-day consistency in the sleep/wake pattern, based on the
probability that an epoch's sleep/wake state matches the state at the
same clock time exactly 24 h later (or earlier), averaged across the
whole recording. Ranges from -100 (perfectly inverted day-to-day) to
+100 (perfectly regular); 0 corresponds to chance-level agreement.

## Usage

``` r
compute_sri(x, epoch_s = NULL, max_gap_min = 30)
```

## Arguments

- x:

  A `zeitr_recording`/`zeitr_result`, or a data frame / tibble with at
  least `datetime` and `state` columns. `state == 1` or `state == 7` is
  treated as sleep, `state == 4` as off-wrist (missing), and any other
  value as wake – matching the coding already used across zeitR's
  pipeline output
  ([`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md),
  [`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md),
  [`export_hypnogram()`](https://zeitr.circadia-lab.uk/reference/export_hypnogram.md)).

- epoch_s:

  `numeric(1)`. Epoch duration in seconds. If `NULL` (default),
  estimated automatically from the median inter-epoch interval.

- max_gap_min:

  `numeric(1)`. Off-wrist gaps of this many minutes or less are
  interpolated rather than excluded. Default `30`, matching Fix 30 and
  Fix 14's 30-minute threshold elsewhere in the pipeline.

## Value

A tibble with columns `participant_id`, `sri`, `n_pairs` (number of
valid 24h-apart epoch comparisons used), and `n_epochs` (total epochs
after off-wrist gap interpolation, before the 24h-pairing step). `sri`
is `NA` if the recording is shorter than 24 h, or if no valid pairs
remain after off-wrist exclusion.

## Details

Ported from **Fix 30** of the Python reference pipeline (`SRI_vallim`):
rather than deriving sleep/wake from a pyActigraphy scoring algorithm
(Sadeh, Cole-Kripke, Roenneberg, Scripps – all of which showed
substantially worse agreement with manual reference scoring in Julia's
concordance analysis, ICC 0.19-0.67 vs 0.82 here), `compute_sri()`
derives sleep/wake directly from the epoch-level `state` column already
produced by zeitR's own pipelines
([`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
/
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md))
– the same classification
[`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
and
[`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
already treat as ground truth for this recording.

Off-wrist handling mirrors Fix 30 exactly: off-wrist epochs
(`state == 4`) are treated as missing. Gaps of `max_gap_min` minutes or
less are interpolated (forward-filled from the last valid epoch, or
back-filled from the next valid epoch when the gap starts at the very
beginning of the recording); longer gaps are left as missing and
excluded from the day-to-day comparison entirely, rather than being
counted as a non-match.

\$\$SRI = -100 + 200 \times \frac{1}{M}\sum\_{t} \Psi(t, t + 24h)\$\$

where \\\Psi(t, t+24h) = 1\\ if the sleep/wake state at epoch \\t\\
matches the state 24 h later, \\0\\ otherwise, and \\M\\ is the number
of epoch pairs where both epochs have a valid (non-missing) state.

## References

Phillips, A. J. K., Clerx, W. M., O'Brien, C. S., Sano, A., Barger, L.
K., Picard, R. W., Lockley, S. W., Klerman, E. B., & Czeisler, C. A.
(2017). Irregular sleep/wake patterns are associated with poorer
academic performance and delayed circadian and sleep/wake timing.
*Scientific Reports*, 7, 3216.
[doi:10.1038/s41598-017-03171-4](https://doi.org/10.1038/s41598-017-03171-4)

## See also

[`compute_npcra()`](https://zeitr.circadia-lab.uk/reference/compute_npcra.md),
[`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md),
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_pipeline_native("recordings/P001.txt", tz = "America/Sao_Paulo")
compute_sri(result)
} # }
```
