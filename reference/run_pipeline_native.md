# Run the Vallim native actigraphy sleep analysis pipeline

Orchestrates the vendor-independent (Vallim) sleep analysis pipeline for
a single ActTrust recording. Shares the read / prepare / off-wrist /
CSPD epoch-scorer stages with
[`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md),
then replaces Condor's `nights_df` classification with the rule set
developed by Julia Ribeiro da Silva Vallim and ported from her condor
pipeline (Fixes 25, 26a/b/c, 27, 29).

## Usage

``` r
run_pipeline_native(
  path,
  tz = "UTC",
  gap_s = 120,
  params = acttrust_params(),
  offwrist_args = list(),
  sleep_args = list(),
  classify_args = list(),
  holidays = NULL,
  free_days = c("Saturday", "Sunday"),
  quiet = FALSE
)
```

## Arguments

- path:

  `character(1)`. Path to the ActTrust `.txt` file.

- tz:

  `character(1)`. Recording time zone. Default `"UTC"`.

- gap_s:

  `numeric(1)`. Gap threshold (seconds) for
  [`check_consistency()`](https://zeitr.circadia-lab.uk/reference/check_consistency.md).
  Default `120`.

- params:

  Device parameter preset, as returned by
  [`acttrust_params()`](https://zeitr.circadia-lab.uk/reference/acttrust_params.md).

- offwrist_args:

  `list`. Additional arguments for
  [`detect_offwrist_bimodal()`](https://zeitr.circadia-lab.uk/reference/detect_offwrist_bimodal.md).

- sleep_args:

  `list`. Additional arguments for
  [`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md).

- classify_args:

  `list`. Additional arguments for
  [`classify_sleep_episodes()`](https://zeitr.circadia-lab.uk/reference/classify_sleep_episodes.md).

- holidays:

  Holidays to treat as free days in addition to the days in `free_days`.
  Accepts three forms, which can be mixed in the same vector:

  - `Date` objects or `"YYYY-MM-DD"` strings for year-specific dates
    (e.g. `as.Date("2019-03-04")` for one Carnival day).

  - `"DD-MM"` strings for dates that recur every year (e.g. `"25-12"`
    for Christmas, `"07-09"` for Brazilian Independence Day). Stored in
    `result$holidays` and auto-forwarded by the `zeitr_result` S3
    methods of
    [`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
    and
    [`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md).
    Default `NULL`.

- free_days:

  A character vector of day names (`"Monday"` through `"Sunday"`,
  case-insensitive) or ISO integers (1 = Monday ... 7 = Sunday)
  identifying which days of the week are unconditionally treated as free
  days. Stored in `result$free_days` and auto-forwarded by the
  `zeitr_result` S3 methods of
  [`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
  and
  [`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md).
  Default `c("Saturday", "Sunday")`.

- quiet:

  `logical(1)`. Suppress timestamp warnings. Default `FALSE`.

## Value

A `zeitr_result` S3 object with the same structure as
[`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md),
except `nights` additionally contains:

- `sleep_type`:

  `character` – `"main"` or `"secondary"`.

- `bed_time`, `get_up_time`:

  POSIXct – `bed_time` is the first sleep epoch, `get_up_time` is the
  first wake epoch (matches the Python reference's `bts`/`gts`
  convention exactly).

The `is_nap` column is retained for backwards compatibility
(`is_nap == (sleep_type == "secondary")`).

## Details

Steps:

1.  **Read** –
    [`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md)

2.  **Consistency check** –
    [`check_consistency()`](https://zeitr.circadia-lab.uk/reference/check_consistency.md)

3.  **Prepare** –
    [`prepare_actigraphy()`](https://zeitr.circadia-lab.uk/reference/prepare_actigraphy.md)

4.  **Off-wrist detection** –
    [`detect_offwrist_bimodal()`](https://zeitr.circadia-lab.uk/reference/detect_offwrist_bimodal.md)

5.  **Epoch scoring** –
    [`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md)
    (CSPD; same scorer as
    [`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md))

6.  **Episode extraction** –
    [`extract_sleep_episodes()`](https://zeitr.circadia-lab.uk/reference/extract_sleep_episodes.md)

7.  **Episode classification** –
    [`classify_sleep_episodes()`](https://zeitr.circadia-lab.uk/reference/classify_sleep_episodes.md)
    (JRSV rules)

8.  **Fine-grained state** –
    [`compute_waso()`](https://zeitr.circadia-lab.uk/reference/compute_waso.md)
    (Cole-Kripke within periods, for `result$data`)

## See also

[`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
for the vendor (Condor) pipeline,
[`extract_sleep_episodes()`](https://zeitr.circadia-lab.uk/reference/extract_sleep_episodes.md),
[`classify_sleep_episodes()`](https://zeitr.circadia-lab.uk/reference/classify_sleep_episodes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_pipeline_native("recordings/P001.txt", tz = "America/Sao_Paulo")
result$nights
} # }
```
