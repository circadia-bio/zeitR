# Sleep Regularity Index (SRI)

Computes the Sleep Regularity Index (Phillips et al. 2017): a measure of
day-to-day consistency in the sleep/wake pattern, based on the
probability that an epoch's sleep/wake state matches the state at the
same clock time exactly 24 h later (or earlier), averaged across the
whole recording. Ranges from -100 (perfectly inverted day-to-day) to
+100 (perfectly regular); 0 corresponds to chance-level agreement.

## Usage

``` r
compute_sri(x, epoch_s = NULL, max_gap_min = 30, algo = "vallim")
```

## Arguments

- x:

  A `zeitr_recording`/`zeitr_result`, or a data frame / tibble with at
  least `datetime` and `state` columns (`algo = "vallim"`) or `datetime`
  and `activity` columns (`algo = "sadeh"`, `"ck"`, `"scripps"`, or
  `"roenneberg"`). For `state`: `state == 1` or `state == 7` is treated
  as sleep, `state == 4` as off-wrist (missing), and any other value as
  wake – matching the coding already used across zeitR's pipeline output
  ([`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md),
  [`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md),
  [`export_hypnogram()`](https://zeitr.circadia-lab.uk/reference/export_hypnogram.md)).

- epoch_s:

  `numeric(1)`. Epoch duration in seconds. If `NULL` (default),
  estimated automatically from the median inter-epoch interval.

- max_gap_min:

  `numeric(1)`. Off-wrist gaps of this many minutes or less are
  interpolated rather than excluded. Default `30`, matching Fix 30 and
  Fix 14's 30-minute threshold elsewhere in the pipeline. Only used when
  `algo = "vallim"`.

- algo:

  `character(1)`. Which sleep/wake source and SRI aggregation to use:
  `"vallim"` (default, uses the pipeline's own `state` column),
  `"sadeh"`, `"ck"`, `"scripps"`, or `"roenneberg"` (all four score raw
  `activity`). See Details for the precise differences beyond just the
  scoring algorithm.

## Value

A tibble with columns `participant_id`, `sri`, `n_pairs` (number of
valid 24h-apart epoch comparisons used; `NA` for `algo != "vallim"`,
whose aggregation isn't a single pooled pair count), and `n_epochs`.
`sri` is `NA` if the recording is shorter than 24 h, or (for
`algo = "vallim"`) if no valid pairs remain after off-wrist exclusion.

## Details

`algo = "vallim"` (default) ports **Fix 30** of the Python reference
pipeline (`SRI_vallim`): sleep/wake comes directly from the epoch-level
`state` column already produced by zeitR's own pipelines
([`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
/
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md))
– the same classification
[`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
and
[`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
already treat as ground truth for this recording. Off-wrist handling
mirrors Fix 30 exactly: off-wrist epochs (`state == 4`) are treated as
missing. Gaps of `max_gap_min` minutes or less are interpolated
(forward-filled from the last valid epoch, or back-filled from the next
valid epoch when the gap starts at the very beginning of the recording);
longer gaps are left as missing and excluded from the day-to-day
comparison entirely, rather than being counted as a non-match. The
aggregation is a single flat average over every valid 24h-apart epoch
pair, pooled across the whole recording: \$\$SRI = -100 + 200 \times
\frac{1}{M}\sum\_{t} \Psi(t, t + 24h)\$\$ where \\\Psi(t, t+24h) = 1\\
if the sleep/wake state at epoch \\t\\ matches the state 24 h later,
\\0\\ otherwise, and \\M\\ is the number of epoch pairs where both
epochs have a valid (non-missing) state. This matches Julia's own
`compute_sri_vallim()` replica (also a flat pooled ratio, not a
per-time-of-day average – see below).

`algo = "sadeh"` instead derives sleep/wake from the raw `activity`
signal via the Sadeh et al. (1994) algorithm, matching pyActigraphy's
actual `Sadeh()`/`SleepRegularityIndex()` (Cell 16 of
`vs_condor_py_pipeline_fix30_jrsv.ipynb`) exactly – both sourced from
pyActigraphy's real code (`pyActigraphy/sleep/scoring_base.py`,
`_sadeh()`; `pyActigraphy/sleep/scoring/sri.py`, `sri()`), not
reconstructed from documentation. Two precise details that differ from
the `"vallim"` path above:

- **Aggregation is a two-step average, not a flat pooled one.**
  pyActigraphy's `sri()` first groups epochs by time-of-day (hour,
  minute, second) across all days, averages day-to-day stability
  *within* each time-of-day slot, and only then averages *across* slots.
  This is mathematically different from a flat pooled average whenever
  slots have unequal numbers of valid day-pairs (e.g. a partial
  first/last day) – both `"vallim"` and `"sadeh"` are ported faithfully
  to their own respective real source, not reconciled to use the same
  aggregation.

- **No off-wrist handling.** Sadeh scores whatever raw activity it is
  given; there is no off-wrist masking/interpolation step in
  pyActigraphy's own code for this path. `max_gap_min` has no effect
  when `algo = "sadeh"`.

Sadeh's own edge behaviour is also reproduced exactly: `mean_W5`/`NAT`
need a full centered 11-epoch window (5 before, self, 5 after) and are
`NA` for the first/last 5 epochs; `sd_Last6` needs a full trailing
6-epoch window (self + 5 before) and is `NA` for the first 5 epochs;
`logAct` uses the *following* epoch's activity (`shift(-1)`) and is `NA`
for the last epoch. Where the resulting `PS` score is `NA`, pandas'
`NaN > threshold` evaluates to `False` (not propagated as missing), so
that epoch is scored wake (`0`), not excluded – reproduced here
explicitly, since R's `NA > threshold` gives `NA`, not `FALSE`.

`algo = "ck"` derives sleep/wake via pyActigraphy's native Cole-Kripke
implementation (`.CK()`, default `settings = "30sec_max_non_overlap"`) –
a DIFFERENT weight set from the Condor-native `ColeKripke` class used
elsewhere in zeitR's pipeline (`R/cole_kripke.R`); the two share an
algorithm family name but are otherwise unrelated. Shares `"sadeh"`'s
two-step SRI aggregation and lack of off-wrist handling (see above).
Despite the reference notebook resampling to 30-second bins before
calling this, that round-trip is a mathematical no-op on data that was
only ever 1-minute resolution (verified by direct execution – see
`?.ck_native_score`), so this operates directly on native-resolution
`activity` with no resampling needed. Also applies Webster's (1982)
rescoring rules afterward, matching pyActigraphy's default. Uses the
*opposite* threshold polarity from Sadeh (`D < threshold` = sleep here,
vs `PS > threshold` = sleep for Sadeh) – both faithful to their own
respective source.

`algo = "scripps"` derives sleep/wake via pyActigraphy's native Scripps
Clinic implementation (`Scripps()`/`_scripps()`). Structurally identical
to `algo = "ck"` – same centered rolling weighted dot product, same
`D < threshold` = sleep polarity, same two-step SRI aggregation and lack
of off-wrist handling – just a different scale/window/threshold, and no
Webster rescoring step (pyActigraphy's `Scripps()` doesn't call
`rescore()`, unlike `CK()`).

`algo = "roenneberg"` derives sleep/wake via pyActigraphy's native
Roenneberg et al. algorithm (`roenneberg()`) – by far the most involved
of the four: trend extraction (a 24h centered rolling mean, allowing a
partial window down to 12h at the recording's edges), a 15%-of-trend
threshold categorization, seed-finding (candidate sleep-onset runs at
least 30 min long), then an iterative correlation-based bout-cleaning
loop (each candidate onset is tested against a family of triangular
"sleep bout ending at position i" templates over the following 12h,
accepting the HIGHEST correlation peak as the bout's offset – not simply
the first qualifying one; see `?.find_highest_peak_idx` for a real
version mismatch this caught). Shares the two-step SRI aggregation and
lack of off-wrist handling with the other three raw-activity algorithms
above. No rescoring step (rescoring is specific to `CK()`).

## References

Phillips, A. J. K., Clerx, W. M., O'Brien, C. S., Sano, A., Barger, L.
K., Picard, R. W., Lockley, S. W., Klerman, E. B., & Czeisler, C. A.
(2017). Irregular sleep/wake patterns are associated with poorer
academic performance and delayed circadian and sleep/wake timing.
*Scientific Reports*, 7, 3216.
[doi:10.1038/s41598-017-03171-4](https://doi.org/10.1038/s41598-017-03171-4)

Sadeh, A., Sharkey, M., & Carskadon, M. A. (1994). Activity-Based
Sleep-Wake Identification: An Empirical Test of Methodological Issues.
*Sleep*, 17(3), 201-207.
[doi:10.1093/sleep/17.3.201](https://doi.org/10.1093/sleep/17.3.201)

## See also

[`compute_npcra()`](https://zeitr.circadia-lab.uk/reference/compute_npcra.md),
[`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md),
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_pipeline_native("recordings/P001.txt", tz = "America/Sao_Paulo")
compute_sri(result)
compute_sri(result, algo = "sadeh")
} # }
```
