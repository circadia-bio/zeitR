# Circular standard deviation of clock times in decimal hours

Computes the circular SD using the mean resultant length formula:
`sqrt(-2 * log(R_bar)) * 24 / (2*pi)`, where `R_bar` is the mean
resultant length of the unit-circle representation of the times. Linear
SD inflates when times straddle midnight; this measure is invariant to
the wrap point (Fix 24 in the JRSV pipeline).

## Usage

``` r
circ_sd_h(x)
```

## Arguments

- x:

  `numeric`. Decimal hours in \[0, 24). `NA` values are silently
  dropped. Returns `NA_real_` for fewer than 2 values.

## Value

A single non-negative `numeric` (hours), or `NA_real_`.

## References

Mardia, K. V., & Jupp, P. E. (2000). *Directional Statistics* (2nd ed.).
Wiley.

## Examples

``` r
circ_sd_h(c(23.5, 0.5))   # small SD for times close to midnight
#> [1] 0.5007167
circ_sd_h(c(6.0, 18.0))   # large SD for times 12 h apart
#> [1] 33.00549
```
