# Batch LIDS ultradian-rhythm summary across a study

Computes
[`compute_lids()`](https://zeitr.circadia-lab.uk/reference/compute_lids.md)
for every participant in a batch of pipeline results and summarises each
participant's quality-filtered bouts into a single row – the LIDS
counterpart to
[`study_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/study_sleep_metrics.md)
and
[`study_summary()`](https://zeitr.circadia-lab.uk/reference/study_summary.md),
making per-participant median cycle length, amplitude, offset, and slope
database-ready for tools like `syncR::sync()`.

## Usage

``` r
study_lids_metrics(results, min_bouts = 3, ...)
```

## Arguments

- results:

  A named list of `zeitr_result` objects, as returned by
  [`run_pipeline_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_batch.md)
  or
  [`run_pipeline_native_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native_batch.md).

- min_bouts:

  `integer(1)`. Minimum number of quality-filtered bouts required for a
  participant's summary metrics to be reported (rather than `NA`, with
  `n_bouts`/`n_bouts_passing` still populated). Default `3`.

- ...:

  Forwarded to
  [`compute_lids()`](https://zeitr.circadia-lab.uk/reference/compute_lids.md)
  for every participant.

## Value

A tibble with one row per participant: `participant_id`, `n_bouts`,
`n_bouts_passing`, and (computed only over
`passes_quality_filter == TRUE` bouts, and only when there are at least
`min_bouts` of them) `period_min_median`, `amplitude_median`,
`offset_median`, `slope_per_60min_median`, plus an `_iqr` variant of
each.

## See also

[`compute_lids()`](https://zeitr.circadia-lab.uk/reference/compute_lids.md),
[`study_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/study_sleep_metrics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
results <- run_pipeline_native_batch("recordings/", tz = "America/Sao_Paulo")
study_lids_metrics(results)
} # }
```
