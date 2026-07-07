# Compute CPD, MSF, MSW, MSFsc, SJL, and SJLa

Calculates chronobiological phenotyping metrics from classified nightly
sleep data. Mirrors `compute_cpd_metrics()` and `nights_to_cpd_df()`
from Julia Ribeiro da Silva Vallim's `pipeline_functions_fix27.py`.

## Usage

``` r
compute_cpd_metrics(
  nights,
  min_tib_h = 3,
  min_tib_eve_h = 3,
  tz = "UTC",
  holidays = NULL
)
```

## Arguments

- nights:

  A `tibble` of nightly sleep statistics as returned by
  [`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)
  or
  [`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md).

- min_tib_h:

  `numeric(1)`. Minimum TBT (hours) for inclusion. Default is `3.0`.

- min_tib_eve_h:

  `numeric(1)`. Minimum TBT (hours) for a night to qualify as a
  free-day-eve night. Default is `3.0`.

- tz:

  `character(1)`. Time zone for extracting clock hours. Default is
  `"UTC"`.

- holidays:

  A `Date` vector of public holidays to treat as free days. Default is
  `NULL` (weekends only).

## Value

A named list with `n_nights_cpd`, `n_free_days`, `n_workdays`, `msw_h`,
`msw_hms`, `msf_h`, `msf_hms`, `msfsc_h`, `msfsc_hms`, `sjl_h`,
`sjl_min`, `sjla_h`, `sjla_min`, `cpd_s`, `cpd_min`, `cpd_h`.

## Details

Mid-sleep is computed per night as: \$\$\text{MS} = \left(\text{SO} +
\frac{\text{offset} - \text{onset}}{2}\right) \bmod 24\$\$ where onset =
bts + SOL and offset = gts − SOI (both in decimal hours).

**MSW** and **MSF** use plain arithmetic means of per-night mid-sleep
values. **MSFsc** adjusts MSF by the free-day-eve sleep onset and the
weighted weekly mean sleep duration when free-day duration exceeds
weekday duration. **CPD** is the RMS distance of each night's mid-sleep
from MSFsc in the time × sequence plane.

## See also

[`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md),
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_pipeline_native("recordings/P001.txt",
                              tz = "America/Sao_Paulo")

cpd <- compute_cpd_metrics(result$nights, tz = "America/Sao_Paulo")
cpd$sjl_min    # social jet lag in minutes
cpd$msf_hms    # mid-sleep on free days as HH:MM
cpd$cpd_min    # CPD in minutes
} # }
```
