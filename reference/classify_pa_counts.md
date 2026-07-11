# Classify physical activity intensity directly from activity counts

Convenience wrapper combining
[`estimate_ee()`](https://zeitr.circadia-lab.uk/reference/estimate_ee.md)
and
[`classify_pa_intensity()`](https://zeitr.circadia-lab.uk/reference/classify_pa_intensity.md)
in one call: converts raw ActTrust(R)/GT3X+ activity counts straight to
a PA intensity factor via the published `"batista2026"` (or future)
equation set, without needing to handle the intermediate MET values.

## Usage

``` r
classify_pa_counts(
  counts,
  device = c("ACTT", "GT3X+"),
  placement = c("hip", "wrist"),
  equation_set = "batista2026"
)
```

## Arguments

- counts:

  Numeric vector of activity counts (counts/min, over the same epoch
  length used by the source equation set – 1-minute for
  `"batista2026"`). Negative values are treated as `0` with a warning,
  since the underlying square-root transform is undefined below zero.

- device:

  `character(1)`. `"ACTT"` (default) or `"GT3X+"`.

- placement:

  `character(1)`. `"hip"` (default) or `"wrist"`.

- equation_set:

  `character(1)`, passed to
  [`pa_equations()`](https://zeitr.circadia-lab.uk/reference/pa_equations.md).
  Default `"batista2026"`.

## Value

An ordered factor, same length as `counts` and same levels as
[`classify_pa_intensity()`](https://zeitr.circadia-lab.uk/reference/classify_pa_intensity.md).

## See also

[`estimate_ee()`](https://zeitr.circadia-lab.uk/reference/estimate_ee.md),
[`classify_pa_intensity()`](https://zeitr.circadia-lab.uk/reference/classify_pa_intensity.md),
[`pa_equations()`](https://zeitr.circadia-lab.uk/reference/pa_equations.md)
for the coefficients and their generalisability caveats.

## Examples

``` r
classify_pa_counts(c(0, 5000, 25000, 50000), device = "ACTT", placement = "hip")
#> [1] light         light         vigorous      very_vigorous
#> Levels: light < moderate < vigorous < very_vigorous
```
