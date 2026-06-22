# Convert integer epoch states to a labelled factor

Converts the integer `state` column produced by the zeitR pipeline into
a human-readable factor. Useful for display, plotting, and export — the
internal `state` column always stays integer to preserve Python
reference parity.

## Usage

``` r
label_states(x)
```

## Arguments

- x:

  integer (or numeric) vector of epoch states, as found in
  `result$data$state`.

## Value

An ordered factor with levels `c("wake", "sleep", "nap", "off-wrist")`,
the same length as `x`.

## Details

|         |               |
|---------|---------------|
| Integer | Label         |
| `0`     | `"wake"`      |
| `1`     | `"sleep"`     |
| `4`     | `"off-wrist"` |
| `7`     | `"nap"`       |

Any value not in the table above is silently converted to `NA`.

## Examples

``` r
label_states(c(0L, 1L, 0L, 4L, 1L, 7L))
#>         0         1         0         4         1         7 
#>      wake     sleep      wake off-wrist     sleep       nap 
#> Levels: wake < sleep < nap < off-wrist
# [1] wake  sleep wake  off-wrist sleep nap
# Levels: wake < sleep < nap < off-wrist

if (FALSE) { # \dontrun{
result <- run_pipeline("recordings/P001.txt")
result$data$state_label <- label_states(result$data$state)
} # }
```
