# Compute PIM, TAT, and ZCM activity counts from raw triaxial acceleration

Reproduces the general wrist-actigraph processing chain – per-axis DC
removal and band-pass filter, vector norm, epoch-level integration –
that converts raw acceleration samples into the epoch-level activity
metrics used throughout the rest of zeitR (e.g.
[`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md)'s
`activity`/`ZCMn` columns) and the wider PA-intensity literature (see
[`pa_equations()`](https://zeitr.circadia-lab.uk/reference/pa_equations.md)).

## Usage

``` r
compute_activity_counts(
  x,
  y,
  z,
  sampling_rate,
  epoch_sec = 60,
  filter_low = 0.5,
  filter_high = 2.7,
  zcm_threshold = 0.01,
  tat_threshold = 0.05,
  metrics = c("PIM", "TAT", "ZCM")
)
```

## Arguments

- x, y, z:

  Numeric vectors of raw triaxial acceleration samples (same units –
  typically g – and length, uniformly sampled at `sampling_rate`).

- sampling_rate:

  `numeric(1)`. Sampling rate of `x`/`y`/`z` in Hz.

- epoch_sec:

  `numeric(1)`. Epoch length in seconds. Default `60` (1-minute epochs,
  the ActTrust standard used elsewhere in zeitR).

- filter_low, filter_high:

  `numeric(1)`. Band-pass cutoffs in Hz, passed to
  [`mrpheus::bandpass_filter()`](https://mrpheus.circadia-lab.uk/reference/bandpass_filter.html).
  Default `0.5`/`2.7` (ActTrust-style, per Batista et al. 2026); use
  `0.25`/`2.5` for GT3X+-style processing instead.

- zcm_threshold:

  `numeric(1)`. Dead-band around zero (same units as `x`/`y`/`z`) a
  filtered axis must clear to register a zero crossing. Default `0.01`.

- tat_threshold:

  `numeric(1)`. Amplitude threshold (same units) the filtered norm must
  exceed to count toward TAT. Default `0.05`.

- metrics:

  [`character()`](https://rdrr.io/r/base/character.html). Which metrics
  to compute. Default `c("PIM", "TAT", "ZCM")` (all three).

## Value

A tibble with one row per epoch: `epoch` (integer index, `1`-based) plus
one column per requested metric.

## What this is (and isn't)

Filtering reuses
[`mrpheus::remove_dc()`](https://mrpheus.circadia-lab.uk/reference/remove_dc.html)
and
[`mrpheus::bandpass_filter()`](https://mrpheus.circadia-lab.uk/reference/bandpass_filter.html)
– the same zero-phase Butterworth implementation already validated as
part of mrpheus's YASA-parity PSG pipeline – rather than a new,
unvalidated filter written from scratch. This makes **mrpheus a required
package for this function specifically** (a `Suggests` dependency of
zeitR, checked at runtime; the rest of zeitR does not need it).

Reusing a validated filter primitive does not make the whole function
validated, though: the PIM/TAT/ZCM epoch-aggregation logic below (what
happens *after* filtering) is still an original implementation with
nothing to check it against – no raw-to-counts converter exists to
compare it to, and Condor's/ActiGraph's exact onboard threshold
constants are proprietary and unpublished. Built from the general
processing description in Batista et al. (2026, PLoS ONE) and standard
actigraphy literature. Treat its output as a reasonable approximation of
device-computed counts, not a validated reproduction – and prefer
device-computed counts
([`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md))
over this function whenever they're available.

## Processing steps

1.  Each of `x`, `y`, `z` has its DC offset removed
    ([`mrpheus::remove_dc()`](https://mrpheus.circadia-lab.uk/reference/remove_dc.html))
    – important for the z-axis in particular, which typically carries a
    ~1 g gravity offset – then is independently band-pass filtered with
    a zero-phase
    ([`mrpheus::bandpass_filter()`](https://mrpheus.circadia-lab.uk/reference/bandpass_filter.html))
    Butterworth filter (`filter_low`-`filter_high` Hz).

2.  A vector norm is computed per sample: `sqrt(x^2 + y^2 + z^2)` on the
    *filtered* axes, matching the order described for ActTrust in
    Batista et al. (2026): filter first, then combine axes.

3.  Samples are grouped into non-overlapping epochs of
    `epoch_sec * sampling_rate` samples. Trailing samples that don't
    fill a complete final epoch are dropped, with a warning.

4.  Per epoch:

    - **PIM** (proportional integration mode) – sum of the absolute
      filtered norm across the epoch.

    - **TAT** (time above threshold) – seconds within the epoch where
      the absolute filtered norm exceeds `tat_threshold`.

    - **ZCM** (zero crossing mode) – number of sign changes beyond a
      `zcm_threshold` dead-band, summed across the three filtered axes,
      within the epoch.

`zcm_threshold` and `tat_threshold` have no published reference value –
they are exposed as tunable parameters rather than given a
false-precision default. The values shipped here are small, physically
motivated starting points (see Arguments), not calibrated constants.

## See also

[`pa_equations()`](https://zeitr.circadia-lab.uk/reference/pa_equations.md),
[`estimate_ee()`](https://zeitr.circadia-lab.uk/reference/estimate_ee.md),
[`classify_pa_counts()`](https://zeitr.circadia-lab.uk/reference/classify_pa_counts.md)
to go from PIM counts to METs and PA intensity bands;
[`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md)
for device-computed counts (preferred over this function when
available);
[`mrpheus::bandpass_filter()`](https://mrpheus.circadia-lab.uk/reference/bandpass_filter.html),
[`mrpheus::remove_dc()`](https://mrpheus.circadia-lab.uk/reference/remove_dc.html)
for the underlying filter primitives.

## Examples

``` r
set.seed(1)
n  <- 25 * 60 * 5  # 5 minutes at 25 Hz
x  <- rnorm(n, sd = 0.05)
y  <- rnorm(n, sd = 0.05)
z  <- rnorm(n, sd = 0.05) + 1  # gravity on the z-axis

if (FALSE) { # \dontrun{
compute_activity_counts(x, y, z, sampling_rate = 25, epoch_sec = 60)
} # }
```
