# Plot a double-plotted actogram

The classic chronobiology double-plot format. Each recording day is
drawn twice: in the left column of its own row (x = 00:00 to 24:00) and
in the right column of the row above (x = 24:00 to 48:00). This means
consecutive pairs of days share a row, making circadian phase drift
visible as a diagonal band across rows.

## Usage

``` r
plot_actogram_double(
  result,
  tz = NULL,
  title = NULL,
  colours = NULL,
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

Row \\i\\ shows day \\i\\ on the left and day \\i+1\\ on the right.
Every day therefore appears twice in the plot (except the first, which
has no left-column predecessor, and the last, which has no right-column
successor). A dashed vertical line marks the 24 h boundary between the
two columns.

## See also

[`plot_actogram()`](https://zeitr.circadia-lab.uk/reference/plot_actogram.md),
[`plot_actogram_activity()`](https://zeitr.circadia-lab.uk/reference/plot_actogram_activity.md),
[`actogram_colours()`](https://zeitr.circadia-lab.uk/reference/actogram_colours.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_pipeline("recordings/P001.txt", tz = "America/Sao_Paulo")
plot_actogram_double(result)
} # }
```
