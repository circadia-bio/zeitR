# Plot a double-plotted actogram with activity bars

Same double-plot layout as
[`plot_actogram_double()`](https://zeitr.circadia-lab.uk/reference/plot_actogram_double.md)
— each day appears in both the left and right column of adjacent rows —
but each epoch is rendered as a vertical bar whose height is
proportional to the raw activity count (default: `ZCMn`, the
zero-crossing mean from ActTrust). Bars are coloured by sleep/wake
state, so the plot simultaneously conveys activity intensity and state
classification.

## Usage

``` r
plot_actogram_activity(
  result,
  tz = NULL,
  title = NULL,
  colours = NULL,
  activity_col = "ZCMn",
  activity_cap_quantile = 0.99,
  date_label_every = 7L,
  epoch_min = 1,
  base_size = 13
)
```

## Arguments

- result:

  A `zeitr_result` list (from
  [`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
  or
  [`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)),
  or a tibble with at least `datetime` and `state` columns.

- tz:

  `character(1)` or `NULL`. Timezone for date and time-of-day
  extraction. `NULL` (default) auto-detects from the timezone attribute
  embedded in the `datetime` POSIXct column; falls back to `"UTC"` when
  absent.

- title:

  `character(1)` or `NULL`. Plot title. `NULL` constructs
  `"Actogram \u2014 <subject_id>"` from `result$subject_id` when
  available.

- colours:

  Named character vector mapping state labels (`"wake"`, `"sleep"`,
  `"nap"`, `"off-wrist"`) to hex colours. `NULL` uses
  [`actogram_colours()`](https://zeitr.circadia-lab.uk/reference/actogram_colours.md).

- activity_col:

  `character(1)`. Name of the column in `result$data` to use as the
  activity signal. Default `"ZCMn"` (zero-crossing mean; present in all
  ActTrust recordings processed by `zeitR`).

- activity_cap_quantile:

  `numeric(1)` in (0, 1\]. Quantile of non-zero activity values used to
  cap bar heights before normalising. Default `0.99` (top 1 % of active
  epochs are clipped to full bar height).

- date_label_every:

  `integer(1)`. Label every Nth row on the y-axis. Default `7` (weekly
  ticks).

- epoch_min:

  `numeric(1)`. Epoch duration in minutes. Used as the tile width in
  [`ggplot2::geom_tile()`](https://ggplot2.tidyverse.org/reference/geom_tile.html).
  Default `1` (ActTrust standard).

- base_size:

  `numeric(1)`. Base font size passed to
  [`ggplot2::theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html).
  Default `13`.

## Value

A `ggplot` object.

## Details

Bar heights are capped at the `activity_cap_quantile` quantile of
non-zero epochs to prevent outlier activity bursts from compressing the
visible range for the rest of the recording. A thin baseline stub (2 %
of row height) is drawn for zero-activity epochs so that sleep periods
and off-wrist blocks remain faintly visible.

## See also

[`plot_actogram()`](https://zeitr.circadia-lab.uk/reference/plot_actogram.md),
[`plot_actogram_double()`](https://zeitr.circadia-lab.uk/reference/plot_actogram_double.md),
[`actogram_colours()`](https://zeitr.circadia-lab.uk/reference/actogram_colours.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_pipeline("recordings/P001.txt", tz = "America/Sao_Paulo")
plot_actogram_activity(result)

# Cap at 95th percentile to highlight moderate activity
plot_actogram_activity(result, activity_cap_quantile = 0.95)
} # }
```
