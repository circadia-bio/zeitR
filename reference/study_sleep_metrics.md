# Batch sleep-timing and chronotype metrics across a study

Computes
[`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
and
[`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
for every participant in a batch of pipeline results and stacks them
into a single tibble with one row per participant – the
sleep-timing/chronotype counterpart to
[`study_summary()`](https://zeitr.circadia-lab.uk/reference/study_summary.md),
which covers NPCRA activity-rhythm variables instead.

## Usage

``` r
study_sleep_metrics(
  results,
  min_tib_h = 5,
  min_tib_eve_h = 3,
  tz = "UTC",
  holidays = NULL,
  free_days = NULL
)
```

## Arguments

- results:

  A named list of `zeitr_result` objects, as returned by
  [`run_pipeline_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_batch.md)
  or
  [`run_pipeline_native_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native_batch.md).
  `participant_id` is taken from each result's own `$subject_id`,
  falling back to the list name if that is unavailable.

- min_tib_h:

  `numeric(1)`. Minimum total in-bed time (hours) for a night to be
  included in
  [`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md).
  Default `5.0`.

- min_tib_eve_h:

  `numeric(1)`. Minimum TBT (hours) for a night to qualify as a
  free-day-eve night in
  [`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md).
  Default `3.0`.

- tz:

  `character(1)`. Time zone for extracting clock hours. Default `"UTC"`.

- holidays, free_days:

  Forwarded to both
  [`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
  and
  [`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
  for every participant. Default `NULL` for both, in which case each
  participant's own `result$holidays`/`result$free_days` (set when the
  pipeline was run) are used instead. Supplying either here overrides
  that per-participant default for the whole study.

## Value

A tibble with one row per participant: `participant_id`, all
[`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
columns (`n_overall`/`n_wd`/`n_fd` and the twelve sleep-timing metrics
with `_wd`/`_fd` suffixes), and all
[`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
columns (`n_nights_cpd`, `n_free_days`, `n_workdays`,
`msw_h`/`msf_h`/`msfsc_h` and their `_hms` forms, `sjl_h`, `sjla_h`,
`cpd_s`/`cpd_min`/`cpd_h`).

## Details

[`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
and
[`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
each return a single named list per participant with no participant
identifier, and there was previously no batch wrapper analogous to
[`study_summary()`](https://zeitr.circadia-lab.uk/reference/study_summary.md)
for them. This is that wrapper – intended to make these chronobiological
phenotyping metrics database-ready for tools like `syncR::sync()`, which
expect one row per participant with a shared `participant_id` column
across sources.

If either metric computation fails for a participant (e.g. no free days
found, or no nights pass the `min_tib_h` filter), that participant's row
is filled with `NA` for the affected metrics and a warning is emitted –
the rest of the study is unaffected.

## See also

[`study_summary()`](https://zeitr.circadia-lab.uk/reference/study_summary.md)
for the NPCRA (activity-rhythm) analogue,
[`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md),
[`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md),
[`run_pipeline_native_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native_batch.md)

## Examples

``` r
if (FALSE) { # \dontrun{
results <- run_pipeline_native_batch("recordings/", tz = "America/Sao_Paulo")
study_sleep_metrics(results)

# Feed straight into syncR::sync()
sync(zeit = study_sleep_metrics(results))
} # }
```
