# Estimate energy expenditure (METs) from activity counts

Applies the published ActTrust(R)/GT3X+ regression equation
(`sqrt(MET) = b0 + b1 * sqrt(activity_counts)`) from
[`pa_equations()`](https://zeitr.circadia-lab.uk/reference/pa_equations.md)
to convert raw activity counts into estimated metabolic equivalents
(METs).

## Usage

``` r
estimate_ee(
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

Numeric vector of estimated METs, same length as `counts`.

## See also

[`pa_equations()`](https://zeitr.circadia-lab.uk/reference/pa_equations.md)
for the underlying coefficients and their limitations,
[`classify_pa_intensity()`](https://zeitr.circadia-lab.uk/reference/classify_pa_intensity.md)
to convert METs into intensity bands.

## Examples

``` r
estimate_ee(c(0, 5000, 25000, 50000), device = "ACTT", placement = "hip")
#> [1] 1.225449 2.990319 6.242013 9.454025
```
