# Non-parametric circadian rhythm analysis (NPCRA)

Computes the standard non-parametric circadian rhythm analysis variables
from an actigraphy recording. All variables are derived from the 24-hour
activity profile following Van Someren et al. (1999) and Marler et al.
(2006).

## Usage

``` r
compute_npcra(x, epoch_s = NULL, L5_hours = 5, M10_hours = 10)
```

## Arguments

- x:

  A `zeitr_recording` as returned by
  [`read_actigraphy()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy.md),
  or a data frame / tibble with at least `datetime` and `activity`
  columns.

- epoch_s:

  `numeric(1)`. Epoch duration in seconds. If `NULL` (default),
  estimated automatically from the median inter-epoch interval.

- L5_hours:

  `numeric(1)`. Width of the least-active window in hours. Default is
  `5`.

- M10_hours:

  `numeric(1)`. Width of the most-active window in hours. Default is
  `10`.

## Value

A one-row tibble with columns `participant_id`, `IS`, `IV`, `RA`, `L5`,
`L5_onset`, `M10`, `M10_onset`, `n_days`, `n_epochs`.

## Details

The following variables are computed:

- `IS`:

  **Interdaily stability** — consistency of the 24 h activity pattern
  across days (range 0–1; higher = more stable).

- `IV`:

  **Intradaily variability** — fragmentation of the rest-activity rhythm
  (\>= 0; higher = more fragmented).

- `RA`:

  **Relative amplitude** — contrast between the most active 10 h window
  (M10) and least active 5 h window (L5) (range 0–1).

- `L5`:

  Mean activity during the least active 5 h window.

- `L5_onset`:

  Clock time of the L5 window midpoint (hh:mm).

- `M10`:

  Mean activity during the most active 10 h window.

- `M10_onset`:

  Clock time of the M10 window midpoint (hh:mm).

## References

Van Someren, E. J. W., Lijzenga, C., Mirmiran, M., & Swaab, D. F.
(1997). Long-term fitness training improves the circadian rest-activity
rhythm in healthy elderly males. *Journal of Biological Rhythms*, 12(2),
146–156.
[doi:10.1177/074873049701200206](https://doi.org/10.1177/074873049701200206)

Marler, M. R., Gehrman, P., Martin, J. L., & Ancoli-Israel, S. (2006).
The sigmoidally transformed cosine curve: a mathematical model for
circadian rhythms with symmetric non-sinusoidal shapes. *Statistics in
Medicine*, 25(22), 3893–3904.
[doi:10.1002/sim.2466](https://doi.org/10.1002/sim.2466)

Van Someren, E. J. W., Swaab, D. F., Colenda, C. C., Cohen, W., McCall,
W. V., & Rosenquist, P. B. (1999). Bright light therapy: Improved
sensitivity to its effects on rest-activity rhythms in Alzheimer
patients by application of nonparametric methods. *Chronobiology
International*, 16(4), 505–518.
[doi:10.3109/07420529908998724](https://doi.org/10.3109/07420529908998724)

## Examples

``` r
if (FALSE) { # \dontrun{
rec   <- read_actigraphy("recordings/P001.txt")
npcra <- compute_npcra(rec)
npcra
} # }
```
