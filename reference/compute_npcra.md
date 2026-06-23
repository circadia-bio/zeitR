# Non-parametric circadian rhythm analysis (NPCRA)

Computes the standard non-parametric circadian rhythm analysis variables
from an actigraphy recording, following Gonçalves et al. (2014) and Van
Someren et al. (1999). All variables are derived from the 24-hour
average activity profile built from **hourly means** (p = 24).

## Usage

``` r
compute_npcra(
  x,
  epoch_s = NULL,
  L5_hours = 5,
  M10_hours = 10,
  window_days = NULL
)
```

## Arguments

- x:

  A `zeitr_recording` as returned by
  [`read_actigraphy()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy.md),
  or a data frame / tibble with at least `datetime` and `activity`
  columns. If a `state` column is present, off-wrist epochs
  (`state == 4`) are excluded before computing all NPCRA variables.

- epoch_s:

  `numeric(1)`. Epoch duration in seconds. If `NULL` (default),
  estimated automatically from the median inter-epoch interval.

- L5_hours:

  `numeric(1)`. Width of the least-active window in hours. Default is
  `5`.

- M10_hours:

  `numeric(1)`. Width of the most-active window in hours. Default is
  `10`.

- window_days:

  `numeric(1)` or `NULL`. If supplied, the recording is split into
  non-overlapping windows of this length (in days) and NPCRA variables
  are computed for each window. A `window_start` column is added to the
  output. Partial final windows (shorter than `window_days`) are
  included but flagged via a lower `n_days` value. Default `NULL`
  computes a single estimate over the full recording.

## Value

A tibble with columns `participant_id`, `window_start` (if `window_days`
is set), `IS`, `IV`, `RA`, `L5`, `L5_onset`, `M10`, `M10_onset`,
`n_days`, `n_epochs`.

## Details

The following variables are computed:

- `IS`:

  **Interdaily stability** — consistency of the 24 h rest-activity
  pattern across days (range 0–1; higher = more stable).

- `IV`:

  **Intradaily variability** — fragmentation of the rest-activity rhythm
  (\>= 0; higher = more fragmented).

- `RA`:

  **Relative amplitude** — contrast between the most active 10 h window
  (M10) and least active 5 h window (L5) (range 0–1).

- `L5`:

  Mean activity during the least active 5 consecutive hours.

- `L5_onset`:

  Clock time of the L5 window onset (hh:mm).

- `M10`:

  Mean activity during the most active 10 consecutive hours.

- `M10_onset`:

  Clock time of the M10 window onset (hh:mm).

## References

Gonçalves, B. S. B., Adamowicz, T., Louzada, F. M., Moreno, C. R., &
Araujo, J. F. (2014). A fresh look at the use of nonparametric analysis
in actimetry. *Sleep Medicine Reviews*, 20, 84–91.
[doi:10.1016/j.smrv.2014.06.002](https://doi.org/10.1016/j.smrv.2014.06.002)

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

# Single estimate over the full recording
compute_npcra(rec)

# Per-fortnight estimates
compute_npcra(rec, window_days = 14)
} # }
```
