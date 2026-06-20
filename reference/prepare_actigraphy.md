# Prepare a raw actigraphy tibble for analysis

Transforms the output of
[`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md)
(or any reader that returns a compatible tibble) into the working data
frame expected by all detection functions. Specifically:

## Usage

``` r
prepare_actigraphy(x)
```

## Arguments

- x:

  A tibble as returned by
  [`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md),
  containing at minimum `int_temp` and `ext_temp` columns.

## Value

A tibble of the same dimensions as `x` with additional or updated
columns: `state`, `offwrist`, `sleep`, `int_temp_`, `ext_temp_`.

## Details

1.  Clamps `int_temp` and `ext_temp` to the physiological range \[0,
    42\] °C.

2.  Adds min-max scaled temperature columns `int_temp_` and `ext_temp_`
    (range \[0, 1\]) for plotting.

3.  Ensures `state`, `offwrist`, and `sleep` columns are present and
    initialised to `0`.

The original tibble is never modified; a copy is returned.

## Examples

``` r
if (FALSE) { # \dontrun{
rec  <- read_acttrust("recordings/P001.txt")
prep <- prepare_actigraphy(rec)
} # }
```
