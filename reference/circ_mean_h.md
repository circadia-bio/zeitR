# Circular mean of clock times in decimal hours

Computes the mean of times on the 0-24 h clock using the unit-circle
method: convert each time to a unit-vector angle, average the sin and
cos components, and back-project via `atan2`. This correctly handles
midnight wrap (e.g. times at 23:30 and 00:30 average to 00:00, not
12:00).

## Usage

``` r
circ_mean_h(x)
```

## Arguments

- x:

  `numeric`. Decimal hours in \[0, 24). `NA` values are silently
  dropped.

## Value

A single `numeric` in \[0, 24), or `NA_real_` if `x` is empty after `NA`
removal.

## Details

The older threshold approach (shift values \>= 12 h by -24 h before
averaging) produces badly biased estimates when times straddle noon –
critically for sleep offset in late sleepers – and is superseded by this
function (Fix 20 in the JRSV pipeline, ICC for sleep offset 0.502 -\>
0.897 on N=404).

## References

Pewsey, A., Neuhauser, M., & Ruxton, G. D. (2013). *Circular Statistics
in R*. Oxford University Press.

## Examples

``` r
circ_mean_h(c(23.5, 0.5))   # midnight wrap -> 0
#> [1] 24
circ_mean_h(c(23.0, 1.0))   # -> 0
#> [1] 0
circ_mean_h(c(7.0, 9.0))    # -> 8
#> [1] 8
```
