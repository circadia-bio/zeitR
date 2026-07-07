# Compute sleep metrics split by day type

Calculates mean sleep metrics separately for all nights, weekday nights,
and weekend/holiday nights. Output column names and derived metrics
mirror `compute_sleep_metrics()` from Julia Ribeiro da Silva Vallim's
`pipeline_functions_fix27.py`.

## Usage

``` r
compute_sleep_metrics(nights, min_tib_h = 5, tz = "UTC", holidays = NULL)
```

## Arguments

- nights:

  A `tibble` of nightly sleep statistics as returned by
  [`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)
  or
  [`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md).
  Must contain at minimum the columns `is_nap`, `bed_time`,
  `get_up_time`, `tbt`, `tst`, `sol`, `soi`, `waso`, `eff`.

- min_tib_h:

  `numeric(1)`. Minimum total in-bed time (hours) for a night to be
  included. Default is `5.0` (matching the Python reference).

- tz:

  `character(1)`. Time zone for extracting clock hours from timestamps.
  Default is `"UTC"`.

- holidays:

  A `Date` vector of public holidays to treat as free days in addition
  to Saturdays and Sundays. Default is `NULL` (weekends only).

## Value

A named list with metrics for three groups (`overall`, `wd` = weekday,
`fd` = free day):

- `n_overall`, `n_wd`, `n_fd`:

  Night counts.

- `sleep_onset_h`, `sleep_offset_h`:

  Circular mean onset and arithmetic mean offset in decimal hours.

- `fpr_tib_h`:

  Mean TBT in hours.

- `fps_h`:

  Mean free period sleep (TBT − SOL − SOI) in hours.

- `tst_h`:

  Mean TST in hours.

- `latencia_min`, `inertia_min`:

  Mean SOL and SOI in minutes.

- `waso_min`:

  Mean WASO in minutes.

- `sleep_eff_pct`:

  Mean sleep efficiency in percent (0–100).

- `tst_24h_h`:

  Same as `tst_h` (24-h TST for main sleep only).

- `dp_midsleep_min`, `dp_tst_min`:

  SD of mid-sleep and TST in minutes.

Workday and free-day metrics carry the suffix `_wd` and `_fd`.

## Details

A night is assigned to the **weekend** group when the get-up date falls
on a Saturday, Sunday, or any date supplied in `holidays`. Weekday
nights are all remaining nights.

`sleep_onset_h` uses a circular mean (values \>= 12 h are shifted by −24
h before averaging, then wrapped to \[0, 24)). `sleep_offset_h` uses a
plain arithmetic mean.

`fps_h` (free period sleep) equals
`fpr_tib_h − (latencia_min + inertia_min) / 60` — TBT net of sleep onset
latency and sleep inertia.

`dp_midsleep_min` and `dp_tst_min` are standard deviations of per-night
mid-sleep (minutes) and TST (minutes), respectively.

## See also

[`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md),
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_pipeline_native("recordings/P001.txt",
                              tz = "America/Sao_Paulo")

sm <- compute_sleep_metrics(result$nights, tz = "America/Sao_Paulo")
sm$tst_h            # mean TST in hours
sm$sleep_onset_h    # mean sleep onset (circular, decimal hours)
sm$dp_midsleep_min  # within-person SD of mid-sleep in minutes
} # }
```
