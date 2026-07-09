# Default state colour palette for actogram plots

Returns a named character vector of hex colours keyed on the four state
labels produced by
[`label_states()`](https://zeitr.circadia-lab.uk/reference/label_states.md).
Pass the output to the `colours` argument of any actogram function to
override individual colours.

## Usage

``` r
actogram_colours()
```

## Value

Named character vector with elements `"wake"`, `"sleep"`, `"nap"`, and
`"off-wrist"`.

## Examples

``` r
actogram_colours()
#>      wake     sleep       nap off-wrist 
#> "#D9C8A0" "#3B2F6B" "#F0A500" "#C25E2A" 
# wake        sleep        nap    off-wrist
# "#D9C8A0" "#3B2F6B" "#F0A500" "#C25E2A"
```
