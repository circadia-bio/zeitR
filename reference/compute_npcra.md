# Non-parametric circadian rhythm analysis (NPCRA)

Computes the standard non-parametric circadian rhythm analysis variables
from an actigraphy recording. `IS` and `IV` are computed from the
hourly-mean activity profile (p = 24), matching pyActigraphy's actual
`_interdaily_stability()`/`_intradaily_variability()` (not the
population-variance formula in Gonçalves et al. 2014 or Van Someren et
al. 1999's own text – the real implementation uses sample variance, ddof
= 1; see Details). `L5`/`M10` and their onsets follow a different
convention – see below – matching the notebook this package's
Vallim-pipeline comparisons were validated against.

## Usage

``` r
compute_npcra(
  x,
  epoch_s = NULL,
  L5_hours = 5,
  M10_hours = 10,
  window_days = NULL,
  trim_to_d1 = TRUE,
  exclude_offwrist = FALSE
)
```

## Arguments

- x:

  A `zeitr_recording` as returned by
  [`read_actigraphy()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy.md),
  or a data frame / tibble with at least `datetime` and `activity`
  columns. If a `state` column is present, off-wrist epochs
  (`state == 4`) are used as-is by default – see `exclude_offwrist`.

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

- trim_to_d1:

  `logical(1)`. If `TRUE` (default), the recording is trimmed to start
  at 00:00 of D+1 – the first full calendar day after recording onset –
  before any NPCRA variable is computed, matching the Python reference
  pipeline's convention (it always starts its NPCRA window at D+1 00:00
  rather than spanning the raw, typically fractional, recording length).
  Set to `FALSE` for the full untrimmed recording (the pre-`trim_to_d1`
  behaviour). If trimming would leave fewer than 2 epochs, a warning is
  emitted and the untrimmed recording is used instead. Off-wrist
  exclusion (`state == 4`, if `exclude_offwrist = TRUE`) still applies
  either way; this does not replicate the Python pipeline's separate
  30-min-threshold off-wrist-run rule for the M10/L5 windows
  specifically (see `exclude_offwrist`) – only the D+1 window start.

- exclude_offwrist:

  `logical(1)`. If `TRUE`, off-wrist epochs (`state == 4`, when a
  `state` column is present) are deleted before computing any NPCRA
  variable. Default `FALSE` matches the actual production Python
  pipeline (Cell 16 of `vs_condor_py_pipeline_fix30_jrsv.ipynb`):
  `_nonparam_metrics()` is always called with `mask_series=None` for
  every variable (IS, IV, M10, L5, RA) – off-wrist periods' raw device
  readings are used as-is, not deleted. Set `TRUE` for the more
  conservative (but non-Python-matching) behaviour of excluding them.
  This is a blunter tool than Python's own off-wrist handling for M10/L5
  specifically (short runs zeroed and kept, long runs excluded via
  `NA` + `min_periods`), which this package does not replicate – `TRUE`
  here simply deletes every off-wrist epoch outright, changing the time
  index rather than leaving gaps.

## Value

A tibble with columns `participant_id`, `window_start` (if `window_days`
is set), `IS`, `IV`, `ISm`, `IVm`, `RA`, `L5`, `L5_onset`, `M10`,
`M10_onset`, `n_days`, `n_epochs`.

## Details

The following variables are computed:

- `IS`:

  **Interdaily stability** — consistency of the 24 h rest-activity
  pattern across days (range 0–1; higher = more stable). From the
  hourly-mean profile (p = 24).

- `IV`:

  **Intradaily variability** — fragmentation of the rest-activity rhythm
  (\>= 0; higher = more fragmented). From the hourly-mean profile (p =
  24).

- `ISm`:

  Mean of `IS` computed at every divisor of 1440 minutes between 1 and
  60 min (22 resolutions total: 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16,
  18, 20, 24, 30, 32, 36, 40, 45, 48, 60), matching Cell 16's
  `_ISm_IVm_FREQS` loop. Unlike `IS`, missing bins at each resolution
  are omitted (not zero-filled).

- `IVm`:

  As `ISm`, for `IV`.

- `RA`:

  **Relative amplitude** — contrast between the most active 10 h window
  (M10) and least active 5 h window (L5) (range 0–1).

- `L5`:

  Mean activity during the least active 5-hour window, found by a
  rolling mean over a 10-min-resampled series, searched globally across
  the whole recording (not the p = 24 hourly profile used for IS/IV).

- `L5_onset`:

  Wall-clock time ("HH:MM") of the *end* of the least-active window –
  the time of day, wrapped at 24 h regardless of which calendar day the
  window actually falls on.

- `M10`:

  As `L5`, for the most active 10-hour window.

- `M10_onset`:

  As `L5_onset`, for the most-active window.

`IS`/`IV` build a 1h-resampled series `X` first: missing hourly bins get
a real zero (matching Python's
`s_1h = s.resample('1h').mean().fillna(0)`), not silent omission. `X` is
grouped by hour-of-day into the p = 24 hourly profile `Xh`. Both
variables then use **sample variance** (divide by n - 1, matching
pandas' `.var()` default) rather than the population variance (divide by
n) that the classic Witting/Van Someren/Gonçalves formulas describe on
paper: \$\$IS =
\frac{\sum_h(\bar{X}\_h-\bar{X})^2/(p-1)}{\sum_i(X_i-\bar{X})^2/(N-1)}\$\$
\$\$IV =
\frac{\sum_i(X_i-X\_{i-1})^2/(N-1)}{\sum_i(X_i-\bar{X})^2/(N-1)}\$\$
with `N` the number of hourly bins in the (zero-filled) recording and
`p` the number of hour-of-day groups present (24 for any recording
spanning a full day). The two formulas share the same denominator,
matching pyActigraphy's `d_1h = data.var()` being computed once and
reused for both.

`ISm`/`IVm` repeat this same computation at 22 other bin widths (every
divisor of 1440 min between 1-60 min) and average the results – but
missing bins at those resolutions are omitted from `N` (matching
`.dropna()`), not zero-filled like the main `IS`/`IV`. A resolution is
only excluded from the average if it has 1 or fewer bins; if the IS/IV
formula itself degenerates (e.g. zero variance) at some other
resolution, that `NaN`/`Inf` is included in the average like any other
value, exactly matching the real Python (a `try/except` around the whole
per-frequency block, not a check on the computed value).

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
