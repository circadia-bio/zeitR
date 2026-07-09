# Compute CPD, MSF, MSW, MSFsc, SJL, and SJLa

Calculates chronobiological phenotyping metrics from classified nightly
sleep data. Mirrors `compute_cpd_metrics()` and `nights_to_cpd_df()`
from Julia Ribeiro da Silva Vallim's `pipeline_functions_fix27.py`, with
two updates to match the fix29 notebook:

- **MSF and MSW** are computed with the circular mean to correctly
  handle mid-sleep times that wrap midnight.

- **Truncated episodes** starting after noon on the last recording day
  are excluded before metric computation.

## Usage

``` r
compute_cpd_metrics(x, ...)

# S3 method for class 'zeitr_result'
compute_cpd_metrics(
  x,
  min_tib_h = 3,
  min_tib_eve_h = 3,
  tz = "UTC",
  holidays = x$holidays,
  free_days = x$free_days,
  ...
)

# Default S3 method
compute_cpd_metrics(
  x,
  min_tib_h = 3,
  min_tib_eve_h = 3,
  tz = "UTC",
  holidays = NULL,
  free_days = c("Saturday", "Sunday"),
  ...
)
```

## Arguments

- x:

  A `zeitr_result` object **or** a `tibble` of nightly sleep statistics
  as returned by
  [`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)
  or
  [`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md).

- ...:

  Not used; reserved for forward compatibility with future methods.

- min_tib_h:

  `numeric(1)`. Minimum TBT (hours) for inclusion. Default is `3.0`.

- min_tib_eve_h:

  `numeric(1)`. Minimum TBT (hours) for a night to qualify as a
  free-day-eve night. Default is `3.0`.

- tz:

  `character(1)`. Time zone for extracting clock hours. Default is
  `"UTC"`.

- holidays:

  Holidays to treat as free days in addition to the days in `free_days`.
  Accepts three forms, which can be mixed in the same vector:

  - `Date` objects or `"YYYY-MM-DD"` strings for year-specific dates
    (e.g. `as.Date("2019-03-04")` for one Carnival day).

  - `"DD-MM"` strings for dates that recur every year (e.g. `"25-12"`
    for Christmas, `"01-01"` for New Year). Default is `NULL`. When `x`
    is a `zeitr_result`, defaults to `x$holidays` automatically. A
    warning is emitted when `NULL`; suppress with
    `options(zeitR.no_holidays_warn = FALSE)`.

- free_days:

  A character vector of day names (`"Monday"` through `"Sunday"`,
  case-insensitive) or ISO integers (1 = Monday ... 7 = Sunday)
  identifying which days of the week are unconditionally treated as free
  days. Default is `c("Saturday", "Sunday")`. When `x` is a
  `zeitr_result`, defaults to `x$free_days` automatically.

## Value

A named list with `n_nights_cpd`, `n_free_days`, `n_workdays`, `msw_h`,
`msw_hms`, `msf_h`, `msf_hms`, `msfsc_h`, `msfsc_hms`, `sjl_h`,
`sjl_min`, `sjla_h`, `sjla_min`, `cpd_s`, `cpd_min`, `cpd_h`.

## Details

Dispatches on the class of `x`:

- **`data.frame` / `tibble`** – treated as a nights table (same
  behaviour as the original function signature).

- **`zeitr_result`** – extracts `x$nights` and inherits `x$holidays` and
  `x$free_days` automatically.

Mid-sleep is computed per night as: \$\$\text{MS} = \left(\text{SO} +
\frac{\text{offset} - \text{onset}}{2}\right) \bmod 24\$\$ where onset =
bts + SOL and offset = gts - SOI (both in decimal hours).

**MSW** and **MSF** use the circular mean of per-night mid-sleep values
(values \>= 12 h shifted by -24 before averaging, then wrapped to \[0,
24)), matching `calculate_msf()` / `calculate_msw()` from the fix29
notebook. **MSFsc** adjusts MSF by the free-day-eve sleep onset and the
weighted weekly mean sleep duration when free-day duration exceeds
weekday duration. **CPD** is the RMS distance of each night's mid-sleep
from MSFsc in the time x sequence plane.

## See also

[`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md),
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)
