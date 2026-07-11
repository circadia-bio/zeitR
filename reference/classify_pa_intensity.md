# Classify physical activity intensity from METs

Bins metabolic equivalent (MET) values into the four PA intensity
classes used throughout Batista et al. (2026) and consistent with WHO PA
intensity conventions: light, moderate, vigorous, very vigorous.

## Usage

``` r
classify_pa_intensity(mets)
```

## Arguments

- mets:

  Numeric vector of MET values, e.g. from
  [`estimate_ee()`](https://zeitr.circadia-lab.uk/reference/estimate_ee.md)
  or measured directly (indirect calorimetry, published equations,
  etc.). MET-based, not device-specific – any source of METs can be
  classified.

## Value

An ordered factor with levels
`c("light", "moderate", "vigorous", "very_vigorous")`, same length as
`mets`. Bands: `[0,3)` light, `[3,6)` moderate, `[6,9)` vigorous,
`[9,Inf)` very vigorous.

## See also

[`estimate_ee()`](https://zeitr.circadia-lab.uk/reference/estimate_ee.md),
[`pa_equations()`](https://zeitr.circadia-lab.uk/reference/pa_equations.md)

## Examples

``` r
mets <- estimate_ee(c(0, 5000, 25000, 50000), device = "ACTT", placement = "hip")
classify_pa_intensity(mets)
#> [1] light         light         vigorous      very_vigorous
#> Levels: light < moderate < vigorous < very_vigorous
```
