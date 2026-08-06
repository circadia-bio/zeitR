# Compute sleep metrics split by day type

Calculates mean sleep metrics separately for all nights, weekday nights,
and weekend/holiday nights. Output column names and derived metrics
mirror `compute_sleep_metrics()` from Julia Ribeiro da Silva Vallim's
`pipeline_functions_fix27.py`.

## Usage

``` r
compute_sleep_metrics(x, ...)

# S3 method for class 'zeitr_result'
compute_sleep_metrics(
  x,
  min_tib_h = 5,
  tz = "UTC",
  holidays = x$holidays,
  free_days = x$free_days,
  ...
)

# Default S3 method
compute_sleep_metrics(
  x,
  min_tib_h = 5,
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
  Must contain at minimum the columns `is_nap`, `bed_time`,
  `get_up_time`, `tbt`, `tst`, `sol`, `soi`, `waso`, `eff`.

- ...:

  Not used; reserved for forward compatibility with future methods.

- min_tib_h:

  `numeric(1)`. Minimum total in-bed time (hours) for a night to be
  included. Default is `5.0` (matching the Python reference).

- tz:

  `character(1)`. Time zone for extracting clock hours from timestamps.
  Default is `"UTC"`.

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

A named list with metrics for three groups (`overall`, `wd` = workday,
`fd` = free day):

- `n_overall`, `n_wd`, `n_fd`:

  Night counts.

- `sleep_onset_h`, `sleep_offset_h`:

  Circular mean onset and offset in decimal hours.

- `fpr_tib_h`:

  Mean TBT in hours.

- `fps_h`:

  Mean free period sleep (TBT - SOL - SOI) in hours.

- `tst_h`:

  Mean TST in hours.

- `latencia_min`, `inertia_min`:

  Mean SOL and SOI in minutes.

- `waso_min`:

  Mean WASO in minutes.

- `sleep_eff_pct`:

  Mean sleep efficiency in percent (0-100).

- `tst_24h_h`:

  Same as `tst_h` (24-h TST for main sleep only).

- `dp_midsleep_min`:

  Circular SD of mid-sleep, in minutes.

- `dp_tst_min`:

  Plain SD of TST, in minutes.

Workday and free-day metrics carry the suffix `_wd` and `_fd`.

## Details

Dispatches on the class of `x`:

- **`data.frame` / `tibble`** – treated as a nights table (same
  behaviour as the original function signature).

- **`zeitr_result`** – extracts `x$nights` and inherits `x$holidays` and
  `x$free_days` automatically.

A night is assigned to the **free-day** group when the get-up date falls
on one of the days in `free_days` or matches an entry in `holidays`.
Day-of-week is determined via the ISO 8601 weekday number
(`format(date, "%u")`) rather than
[`weekdays()`](https://rdrr.io/r/base/weekday.POSIXt.html), which is
locale-dependent. All remaining nights are **workday** nights.

Both `sleep_onset_h` and `sleep_offset_h` use a circular mean, treating
clock hours as angles on a 24-hour circle so values near midnight from
opposite sides average correctly (e.g. 23:30 and 00:30 average to 00:00,
not 12:00).

`fps_h` (free period sleep) equals
`fpr_tib_h - (latencia_min + inertia_min) / 60` – TBT net of sleep onset
latency and sleep inertia.

`dp_midsleep_min` is the circular SD (R-bar mean-resultant-length
formula) of per-night mid-sleep, in minutes – matching production
Python's `_std_circular_h()` exactly, not a plain SD. `dp_tst_min` is a
plain SD of per-night TST, in minutes.

## See also

[`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md),
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# From a zeitr_result: holidays and free_days forwarded automatically
result <- run_pipeline_native("recordings/P001.txt",
                              tz        = "America/Sao_Paulo",
                              holidays  = my_holidays,
                              free_days = c("Saturday", "Sunday"))
sm <- compute_sleep_metrics(result, tz = "America/Sao_Paulo")

# Non-standard schedule: Friday + Saturday as free days
sm <- compute_sleep_metrics(result$nights,
                            tz        = "America/Sao_Paulo",
                            holidays  = my_holidays,
                            free_days = c("Friday", "Saturday"))

sm$tst_h            # mean TST in hours
sm$sleep_onset_h    # mean sleep onset (circular, decimal hours)
sm$dp_midsleep_min  # within-person SD of mid-sleep in minutes
} # }
```
