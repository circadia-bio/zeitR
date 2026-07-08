# Export a zeitR pipeline result as a hypnoR-compatible hypnogram

Converts the epoch-level `data` tibble from
[`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
or
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)
into the tidy hypnogram format expected by `hypnoR` metric functions.

## Usage

``` r
export_hypnogram(
  result,
  subject_id = NULL,
  source = "zeitR",
  drop_offwrist = FALSE,
  epoch_sec = NULL
)
```

## Arguments

- result:

  A `zeitr_result` list as returned by
  [`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
  or
  [`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md),
  or a tibble with at minimum the columns `datetime` and `state`.

- subject_id:

  `character(1)` or `NULL`. Written to the `subject_id` column and
  always takes precedence. When `NULL` (default), `result$subject_id` is
  used if present; the column is omitted entirely when no ID is
  available from either source.

- source:

  `character(1)`. Label written to the `source` column. Default
  `"zeitR"`.

- drop_offwrist:

  `logical(1)`. Remove off-wrist epochs (`state == 4`) from the output
  and re-index epochs. Default `FALSE`.

- epoch_sec:

  `numeric(1)`. When supplied, warns if the observed epoch duration
  differs from this value. Default `NULL`.

## Value

A tibble with columns:

- `epoch`:

  Integer epoch index, 1-based.

- `time`:

  `POSIXct` timestamp for the start of each epoch.

- `stage`:

  Ordered factor: `c("W", "Sleep", "Quiet sleep")`.

- `subject_id`:

  Character subject identifier (omitted when not available).

- `source`:

  Character scorer label.

## Details

The coarse (3-state) stage mapping used by `zeitR` is:

|          |        |                 |
|----------|--------|-----------------|
| `state`  | `ZCMn` | Stage           |
| `0`      | any    | `"W"`           |
| `1`, `7` | `> 0`  | `"Sleep"`       |
| `1`, `7` | `== 0` | `"Quiet sleep"` |
| `4`      | any    | `"W"`           |

Zero-count epochs within sleep (`ZCMn == 0`) are mapped to
`"Quiet sleep"` as the standard actigraphy proxy for quiet/deep sleep.
When `ZCMn` is not present in the data, all sleep epochs are mapped to
`"Sleep"`.

The `stage` column is an ordered factor with levels
`c("W", "Sleep", "Quiet sleep")`, matching `hypnoR`'s coarse resolution
contract. `"Quiet sleep"` is not produced by actigraphy but the level is
present so downstream `hypnoR` functions do not throw factor-level
errors.

`subject_id` is taken from `result$subject_id` when not explicitly
supplied. Both
[`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
and
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)
derive this from the input filename stem automatically, so in most cases
no manual override is needed.

## See also

[`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md),
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md),
[`label_states()`](https://zeitr.circadia-lab.uk/reference/label_states.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_pipeline("recordings/P001.txt", tz = "America/Sao_Paulo")

# subject_id inferred from filename automatically
hyp <- export_hypnogram(result)
hyp$subject_id  # "P001"

# Override when filename does not match study code
hyp <- export_hypnogram(result, subject_id = "STUDY_001")

# Batch: subject_id inferred for every file automatically
results <- run_pipeline_batch("recordings/", tz = "America/Sao_Paulo")
hyps    <- lapply(results, export_hypnogram)

# Drop off-wrist epochs
hyp_clean <- export_hypnogram(result, drop_offwrist = TRUE)
} # }
```
