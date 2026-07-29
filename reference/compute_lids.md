# Compute LIDS ultradian-rhythm parameters for every sleep bout

The main entry point for the LIDS pipeline: extracts sleep bouts,
applies
[`lids_transform()`](https://zeitr.circadia-lab.uk/reference/lids_transform.md)
and [`fit_lids()`](https://zeitr.circadia-lab.uk/reference/fit_lids.md)
to each, and returns one row per bout with cosine-fit parameters and a
quality-filter flag. Ports the full pipeline described in Winnebeck et
al. (2018) and Hammad et al. (2026).

## Usage

``` r
compute_lids(
  x,
  bout_source = c("auto", "state", "roenneberg"),
  activity_col = "activity",
  method = c("gaussian", "mva"),
  win_min = 30,
  sigma_min = 5,
  period_range = c(30, 180),
  period_step = 2,
  duration_range = c(3, 12),
  min_r = 0.4,
  max_p = 0.05,
  offset_bounds = c(1, 99),
  bout_args = list()
)
```

## Arguments

- x:

  A `zeitr_result`, `zeitr_recording`, or a data frame / tibble with at
  least `datetime` and `activity` columns (and, for
  `bout_source = "state"`/`"auto"`, a `state` column).

- bout_source:

  `character(1)`. `"auto"` (default), `"state"`, or `"roenneberg"`. See
  Details.

- activity_col:

  `character(1)`. Name of the activity column in `x`. Default
  `"activity"`.

- method, win_min, sigma_min:

  Forwarded to
  [`lids_transform()`](https://zeitr.circadia-lab.uk/reference/lids_transform.md).

- period_range, period_step:

  Forwarded to
  [`fit_lids()`](https://zeitr.circadia-lab.uk/reference/fit_lids.md).

- duration_range:

  `numeric(2)`, hours. Bout duration bounds (both `bout_source` paths).
  Default `c(3, 12)`.

- min_r, max_p, offset_bounds:

  Quality-filter thresholds; see Details.

- bout_args:

  Named list of additional arguments forwarded to
  [`detect_lids_bouts()`](https://zeitr.circadia-lab.uk/reference/detect_lids_bouts.md)
  when `bout_source = "roenneberg"` (e.g. `relative_threshold`,
  `main_window`). Default [`list()`](https://rdrr.io/r/base/list.html).

## Value

A tibble with one row per bout: `participant_id`, `bout_id`,
`bout_start`, `bout_end`, `duration_h`, `period_min`, `amplitude`,
`offset`, `slope_per_60min`, `phase_rad`, `pearson_r`, `p_value`, `mri`,
`passes_quality_filter`.

## Where bouts come from (`bout_source`)

- `"state"` – uses the epoch-level `state` column already produced by
  zeitR's own pipelines
  ([`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
  /
  [`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)):
  contiguous `state == 1` runs are treated as bouts, filtered by
  `duration_range`. Off-wrist (`state == 4`) epochs break a run rather
  than being bridged.

- `"roenneberg"` – ignores any existing `state` column and runs the
  independent
  [`detect_lids_bouts()`](https://zeitr.circadia-lab.uk/reference/detect_lids_bouts.md)
  relative-immobility detector directly on the activity signal. Use this
  for standalone recordings that haven't been run through zeitR's
  Crespo/Vallim pipelines, or to reproduce Winnebeck/Hammad's own
  bout-detection method rather than zeitR's.

- `"auto"` (the default) – `"state"` if a `state` column is present in
  `x`, otherwise `"roenneberg"`.

## Quality filtering

Following Winnebeck et al. (2018) and Hammad et al. (2026), a bout
*passes* quality filtering when all of:

- `pearson_r >= min_r` (default `0.4` – a soft data-quality threshold,
  not a hard significance test; ~75% of adult bouts cleared this bar in
  Winnebeck et al. 2018),

- `p_value <= max_p` (default `0.05`),

- `offset_bounds[1] < offset < offset_bounds[2]` (default
  `1 < offset < 99`), excluding spuriously flat bouts (e.g. a
  lost/removed device).

Bouts failing quality filtering are still returned (with
`passes_quality_filter = FALSE`) rather than dropped, so callers can
inspect what was excluded.

## References

Winnebeck, E. C., Fischer, D., Leise, T., & Roenneberg, T. (2018).
Dynamics and Ultradian Structure of Human Sleep in Real Life. *Current
Biology*, 28(1), 49-59.e5.
[doi:10.1016/j.cub.2017.11.063](https://doi.org/10.1016/j.cub.2017.11.063)

Hammad, G., Schoch, S. F., Engelmann, M., Spock, Z., Kurth, S., &
Winnebeck, E. C. (2026). Charting infant sleep cycle development using
actigraphy: Longitudinal evidence for ultradian cycle lengthening within
the first year of life. *SLEEP*.

## See also

[`lids_transform()`](https://zeitr.circadia-lab.uk/reference/lids_transform.md),
[`fit_lids()`](https://zeitr.circadia-lab.uk/reference/fit_lids.md),
[`detect_lids_bouts()`](https://zeitr.circadia-lab.uk/reference/detect_lids_bouts.md),
[`study_lids_metrics()`](https://zeitr.circadia-lab.uk/reference/study_lids_metrics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_pipeline_native("recordings/P001.txt", tz = "America/Sao_Paulo")
compute_lids(result)
} # }
```
