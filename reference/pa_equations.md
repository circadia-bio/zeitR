# Published activity-count equations for estimating METs and PA intensity

Returns the regression coefficients and count-based cut-points for
estimating energy expenditure (METs) and classifying physical activity
(PA) intensity from ActTrust(R) (Condor Instruments) or ActiGraph(R)
GT3X+ activity counts, hip or wrist placement.

## Usage

``` r
pa_equations(equation_set = "batista2026")
```

## Arguments

- equation_set:

  `character(1)`. Currently only `"batista2026"` is available.

## Value

A tibble with one row per device x placement combination and columns:

- `device`:

  `"ACTT"` or `"GT3X+"`.

- `placement`:

  `"hip"` or `"wrist"`.

- `b0`, `b1`:

  Intercept and slope of `sqrt(MET) = b0 + b1 * sqrt(activity_counts)`.

- `cut3`, `cut6`, `cut9`:

  Published count/min cut-points (95% CI midpoints) at the 3, 6, and 9
  MET thresholds – i.e. the count value at which the fitted equation
  crosses that MET level. Provided for reference;
  [`classify_pa_intensity()`](https://zeitr.circadia-lab.uk/reference/classify_pa_intensity.md)
  classifies on estimated METs directly rather than re-deriving these.

## Important caveats

These equations come from a **single** controlled-laboratory validation
study (Batista et al. 2026; N = 56 healthy adults aged 18-35; treadmill
walking/running at 3-9 km/h only). Treat them as one available equation
set, not a universal standard:

- The paper's own Discussion section compares its GT3X+ (hip) cut-points
  against two other published GT3X+ studies (Sasaki et al. 2011;
  Santos-Lozano et al. 2013) and finds differences of 2-65% depending on
  the MET threshold – and those two reference studies differ from *each
  other* by 16-39%. The paper attributes that spread mainly to
  sample-level characteristics rather than device or methodological
  artefacts, but that reading applies to the comparison *between* Sasaki
  and Santos-Lozano – it doesn't fully carry over to comparisons against
  this paper's own equations, since all three studies used different
  modelling approaches (Sasaki et al.: ActiGraph's two-regression model;
  Santos-Lozano et al.: an artificial neural network; this paper: a
  simple sqrt-transformed linear model with a device x placement
  interaction). Cross-study cut-point differences therefore reflect
  model choice as well as sample, not sample alone – and that's
  checkable to different degrees: this paper's linear coefficients
  (`b0`/`b1` below) are transparent enough to compare term-by-term
  against Sasaki et al.'s two-regression model, so a discrepancy there
  is at least diagnosable. Santos-Lozano et al.'s cut-points come from
  an ANN, which has no inspectable coefficients – there is no way to say
  *why* it disagrees with the equations here, only *that* it does.

- The ACTT (hip)/ACTT (wrist) equations are the first published
  cut-points for ActTrust(R) at all, so there is nothing yet to
  cross-check them against.

- Validated only for laboratory treadmill walking/running in healthy
  young adults. Applying these equations to free-living data, other age
  groups (children, older adults), clinical populations, or other
  activity types is an extrapolation the source study explicitly flags
  as untested.

`equation_set` is exposed as an explicit argument (rather than
hard-coding a single table) so a future validation study covering a
different population or device can be added as an alternative set
without changing the
[`estimate_ee()`](https://zeitr.circadia-lab.uk/reference/estimate_ee.md)
/
[`classify_pa_intensity()`](https://zeitr.circadia-lab.uk/reference/classify_pa_intensity.md)
API.

## See also

[`estimate_ee()`](https://zeitr.circadia-lab.uk/reference/estimate_ee.md),
[`classify_pa_intensity()`](https://zeitr.circadia-lab.uk/reference/classify_pa_intensity.md)

## Examples

``` r
pa_equations()
#> # A tibble: 4 × 7
#>   device placement    b0     b1  cut3  cut6  cut9
#>   <chr>  <chr>     <dbl>  <dbl> <dbl> <dbl> <dbl>
#> 1 GT3X+  hip        1.06 0.0199  1132  4853  9468
#> 2 ACTT   hip        1.11 0.0088  5057 23339 46410
#> 3 ACTT   wrist      1.23 0.0081  3761 22368 47203
#> 4 GT3X+  wrist      1.21 0.0127  1698  9503 19787
```
