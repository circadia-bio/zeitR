# Non-parametric circadian rhythm analysis

## Overview

Non-parametric circadian rhythm analysis (NPCRA) characterises the
24-hour rest-activity pattern without assuming a sinusoidal waveform. It
was introduced by Van Someren et al. (1997, 1999) and has become the
standard approach for actigraphy-based circadian assessment in both
healthy and clinical populations.

zeitR implements
[`compute_npcra()`](https://zeitr.circadia-lab.uk/reference/compute_npcra.md),
which returns five core variables from a single function call.

------------------------------------------------------------------------

## The five variables

### Interdaily stability (IS)

IS quantifies how consistent the 24-hour activity pattern is across
days. It is computed as the ratio of the variance of the average hourly
activity profile to the overall variance of the raw signal:

``` math
\text{IS} = \frac{n}{p} \cdot
\frac{\sum_{h=1}^{p} \left( \bar{x}_h - \bar{x} \right)^2}
     {\sum_{i=1}^{n} \left( x_i   - \bar{x} \right)^2}
```

where $`n`$ is the total number of epochs, $`p`$ is the number of epochs
per day, $`\bar{x}_h`$ is the mean activity at hour-of-day $`h`$, and
$`\bar{x}`$ is the grand mean.

**Range:** 0–1. Higher values indicate a more stable, clock-driven
rhythm. Typical values in healthy adults are 0.6–0.8.

### Intradaily variability (IV)

IV measures the fragmentation of the rest-activity rhythm — how often
and how abruptly the signal transitions between rest and activity:

``` math
\text{IV} = \frac{n \sum_{i=2}^{n} (x_i - x_{i-1})^2}
                 {(n-1) \sum_{i=1}^{n} (x_i - \bar{x})^2}
```

**Range:** $`\geq 0`$. Higher values indicate a more fragmented rhythm.
Values above 1.0 are often seen in older adults and patients with
neurological disorders.

### Relative amplitude (RA)

RA captures the contrast between the most active and least active
periods of the day:

``` math
\text{RA} = \frac{M10 - L5}{M10 + L5}
```

**Range:** 0–1. Higher values indicate a stronger amplitude of the
circadian rhythm.

### L5 and M10

**L5** is the mean activity level during the least active 5-hour window
in the 24-hour averaged profile. **M10** is the mean activity during the
most active 10-hour window. Both are computed over the 24-hour average
profile (not raw epochs), making them robust to day-to-day variability.

The onset times (`L5_onset`, `M10_onset`) give the start time of the
respective windows in `HH:MM` format.

------------------------------------------------------------------------

## Computing NPCRA

[`compute_npcra()`](https://zeitr.circadia-lab.uk/reference/compute_npcra.md)
accepts either a `zeitr_recording` (from
[`read_actigraphy()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy.md))
or any data frame with `datetime` and `activity` columns.

``` r

rec   <- read_actigraphy("recordings/P001.txt", tz = "America/Sao_Paulo")
npcra <- compute_npcra(rec)
npcra
```

    ## # A tibble: 1 x 10
    ##   participant_id    IS    IV    RA    L5 L5_onset   M10 M10_onset n_days n_epochs
    ##   <chr>          <dbl> <dbl> <dbl> <dbl> <chr>    <dbl> <chr>      <dbl>    <int>
    ## 1 P001            0.72  0.43  0.89  12.3 02:30    84.7  11:00       7.0     10080

### Interpretation example

For participant P001:

- **IS = 0.72** — a moderately stable rhythm, consistent with a healthy
  adult with regular working hours.
- **IV = 0.43** — low fragmentation; the rest-activity transitions are
  relatively smooth.
- **RA = 0.89** — strong circadian amplitude; there is a large contrast
  between the active and rest phases.
- **L5 onset = 02:30** — the least active window is centred around
  02:30, as expected for a typical sleeper.
- **M10 onset = 11:00** — the most active window peaks at mid-morning.

------------------------------------------------------------------------

## Customising window widths

The 5-hour and 10-hour windows are defaults. They can be adjusted:

``` r

# Narrower windows — more sensitive to sharp peaks
compute_npcra(rec, L5_hours = 4, M10_hours = 8)
```

Note that IS and IV are independent of the window parameters; only L5,
M10, and RA are affected.

------------------------------------------------------------------------

## Using raw data frames

If you have activity data in a plain data frame (e.g. imported from
another tool),
[`compute_npcra()`](https://zeitr.circadia-lab.uk/reference/compute_npcra.md)
works directly:

``` r

df <- data.frame(
  datetime = seq(
    as.POSIXct("2021-05-27 00:00:00", tz = "UTC"),
    by = 60,
    length.out = 10080
  ),
  activity = c(...)  # your activity counts
)

compute_npcra(df, epoch_s = 60)
```

The `epoch_s` argument is optional — if omitted it is estimated from the
median inter-epoch interval in `datetime`.

------------------------------------------------------------------------

## Batch NPCRA across a study

For multiple participants,
[`study_summary()`](https://zeitr.circadia-lab.uk/reference/study_summary.md)
calls
[`compute_npcra()`](https://zeitr.circadia-lab.uk/reference/compute_npcra.md)
for every recording in a `zeitr_study` and returns a single summary
tibble. See
[`vignette("study-analysis")`](https://zeitr.circadia-lab.uk/articles/study-analysis.md)
for a full walkthrough.

------------------------------------------------------------------------

## References

Van Someren, E. J. W., Lijzenga, C., Mirmiran, M., & Swaab, D. F.
(1997). Long-term fitness training improves the circadian rest-activity
rhythm in healthy elderly males. *Journal of Biological Rhythms*, 12(2),
146–156. <https://doi.org/10.1177/074873049701200206>

Van Someren, E. J. W., Swaab, D. F., Colenda, C. C., Cohen, W., McCall,
W. V., & Rosenquist, P. B. (1999). Bright light therapy: Improved
sensitivity to its effects on rest-activity rhythms in Alzheimer
patients by application of nonparametric methods. *Chronobiology
International*, 16(4), 505–518.
<https://doi.org/10.3109/07420529908998724>

Marler, M. R., Gehrman, P., Martin, J. L., & Ancoli-Israel, S. (2006).
The sigmoidally transformed cosine curve: a mathematical model for
circadian rhythms with symmetric non-sinusoidal shapes. *Statistics in
Medicine*, 25(22), 3893–3904. <https://doi.org/10.1002/sim.2466>
