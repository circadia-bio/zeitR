# Fit a sloped cosine to a LIDS profile, scanning candidate periods

For a fixed period \\T\\, the sloped-cosine model \$\$f(t) =
\alpha\cos(2\pi t/T) + \beta\sin(2\pi t/T) + b + s\\t\$\$ is linear in
\\(\alpha,\beta,b,s)\\, so it is solved by ordinary least squares rather
than a non-linear optimiser (Hammad et al. 2026). Candidate periods are
scanned over `period_range` in steps of `period_step`, and for each one
the fit's Munich Rhythmicity Index (\\\text{MRI} = 2 \times
\text{amplitude} \times r\\, Winnebeck et al. 2018) is computed. The
period with the highest MRI is returned as the bout's estimated
ultradian cycle length – the same selection rule used by
`pyActigraphy.analysis.LIDS.lids_fit()`, generalised here to include the
linear slope term.

## Usage

``` r
fit_lids(lids, epoch_min = 1, period_range = c(30, 180), period_step = 2)
```

## Arguments

- lids:

  Numeric vector of (smoothed) LIDS values for one sleep bout, as
  returned by
  [`lids_transform()`](https://zeitr.circadia-lab.uk/reference/lids_transform.md).
  Must not contain `NA`.

- epoch_min:

  `numeric(1)`. Epoch duration in minutes. Default `1`.

- period_range:

  `numeric(2)`. Candidate period bounds in minutes. Default `c(30, 180)`
  (Hammad et al. 2026, tuned for infant/ultradian cycles); use
  `c(60, 180)` with `period_step = 5` to match Winnebeck et al. (2018)'s
  adult/adolescent scan.

- period_step:

  `numeric(1)`. Step size in minutes for the period scan. Default `2`.

## Value

A named list:

- `period_min`:

  Estimated cycle length (minutes) at peak MRI.

- `amplitude`:

  \\\sqrt{\alpha^2+\beta^2}\\.

- `phase_rad`:

  \\-\mathrm{atan2}(\beta,\alpha)\\, radians; `0` = LIDS peak at bout
  start.

- `offset`:

  Inactivity level at bout start (\\b\\).

- `slope_per_60min`:

  Linear trend, rescaled to LIDS units per hour.

- `pearson_r`:

  Correlation between fitted and observed LIDS.

- `p_value`:

  Two-sided p-value for `pearson_r`
  ([`stats::cor.test()`](https://rdrr.io/r/stats/cor.test.html)).

- `mri`:

  Munich Rhythmicity Index at the selected period.

## References

Winnebeck, E. C., Fischer, D., Leise, T., & Roenneberg, T. (2018).
Dynamics and Ultradian Structure of Human Sleep in Real Life. *Current
Biology*, 28(1), 49-59.e5.
[doi:10.1016/j.cub.2017.11.063](https://doi.org/10.1016/j.cub.2017.11.063)

Hammad, G., Schoch, S. F., Engelmann, M., Spock, Z., Kurth, S., &
Winnebeck, E. C. (2026). Charting infant sleep cycle development using
actigraphy. *SLEEP*.

## See also

[`lids_transform()`](https://zeitr.circadia-lab.uk/reference/lids_transform.md),
[`compute_lids()`](https://zeitr.circadia-lab.uk/reference/compute_lids.md)

## Examples

``` r
set.seed(1)
t <- seq(0, 300, by = 1)
lids <- 85 + 15 * cos(2 * pi * t / 60) - 0.05 * t + rnorm(length(t), sd = 2)
fit_lids(lids)
#> $period_min
#> [1] 60
#> 
#> $amplitude
#> [1] 15.18198
#> 
#> $phase_rad
#> [1] 0.007157598
#> 
#> $offset
#> [1] 85.16796
#> 
#> $slope_per_60min
#> [1] -3.038271
#> 
#> $pearson_r
#> [1] 0.9866464
#> 
#> $mri
#> [1] 29.9585
#> 
#> $p_value
#> [1] 1.036419e-237
#> 
```
