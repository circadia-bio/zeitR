# Score actigraphy epochs as wake or sleep using the Cole-Kripke algorithm

Applies the Cole-Kripke algorithm to a vector of zero-crossing mode
(ZCM) activity counts, scoring each epoch as wake (`1`) or sleep (`0`)
using a weighted sum of activity in a surrounding window.

## Usage

``` r
score_epochs_cole_kripke(
  zcm,
  P = 0.000464,
  weights_before = c(34.5, 133, 529, 375, 408, 400.5, 1074, 2048.5, 2424.5),
  weights_after = c(1920, 149.5, 257.5, 125, 111.5, 120, 69, 40.5)
)
```

## Arguments

- zcm:

  `numeric` vector of ZCM activity counts, one value per epoch.

- P:

  `numeric(1)`. Scaling factor. Default is `0.000464` (Cole et al.,
  1992).

- weights_before:

  `numeric(9)`. Weights applied to the 9 epochs *before* the current
  epoch. Defaults to the values from Cole et al. (1992), Table 2.

- weights_after:

  `numeric(8)`. Weights applied to the 8 epochs *after* the current
  epoch. Defaults to the values from Cole et al. (1992), Table 2.

## Value

An integer vector the same length as `zcm`, with `1` indicating wake and
`0` indicating sleep.

## Details

Each epoch's score is computed as:

\$\$D_i = P \sum\_{j=1}^{9} W_j^{-} \cdot A\_{i-j} + P \sum\_{j=1}^{8}
W_j^{+} \cdot A\_{i+j}\$\$

where \\A_i\\ is the ZCM count at epoch \\i\\, \\W^{-}\\ and \\W^{+}\\
are the before and after weight vectors from Cole et al. (1992), and \\P
= 0.000464\\. Epochs with \\D_i \ge 1\\ are scored as wake.

## References

Cole, R. J., Kripke, D. F., Gruen, W., Mullaney, D. J., & Gillin, J. C.
(1992). Automatic sleep/wake identification from wrist activity.
*Sleep*, 15(5), 461–469.
[doi:10.1093/sleep/15.5.461](https://doi.org/10.1093/sleep/15.5.461)

## Examples

``` r
if (FALSE) { # \dontrun{
rec    <- read_acttrust("recordings/P001.txt")
scores <- score_epochs_cole_kripke(rec$ZCMn)
table(scores)  # 0 = sleep, 1 = wake
} # }
```
