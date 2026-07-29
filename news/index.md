# Changelog

## zeitR (development version)

### ✨ New features

- New LIDS (Locomotor Inactivity During Sleep) module –
  [`lids_transform()`](https://zeitr.circadia-lab.uk/reference/lids_transform.md),
  [`fit_lids()`](https://zeitr.circadia-lab.uk/reference/fit_lids.md),
  [`detect_lids_bouts()`](https://zeitr.circadia-lab.uk/reference/detect_lids_bouts.md),
  [`compute_lids()`](https://zeitr.circadia-lab.uk/reference/compute_lids.md),
  and
  [`study_lids_metrics()`](https://zeitr.circadia-lab.uk/reference/study_lids_metrics.md)
  – porting the ultradian-rhythm methodology of Winnebeck et al. (2018,
  *Current Biology*) and its infant extension in Hammad et al. (2026,
  *SLEEP*, <https://zenodo.org/records/18199381>).
  [`lids_transform()`](https://zeitr.circadia-lab.uk/reference/lids_transform.md)
  applies the `100/(1+x)` non-linear transform plus Gaussian (Hammad
  2026, default) or moving-average (Winnebeck 2018 / `pyActigraphy`)
  smoothing;
  [`fit_lids()`](https://zeitr.circadia-lab.uk/reference/fit_lids.md)
  scans candidate periods (30-180 min by default) with an OLS
  sloped-cosine fit, selecting the period with the highest Munich
  Rhythmicity Index.
  [`compute_lids()`](https://zeitr.circadia-lab.uk/reference/compute_lids.md)
  is the main entry point: it extracts sleep bouts either from an
  existing zeitR pipeline’s `state` column (`bout_source = "state"`) or
  via the new standalone
  [`detect_lids_bouts()`](https://zeitr.circadia-lab.uk/reference/detect_lids_bouts.md)
  Roenneberg relative-immobility detector (`bout_source = "roenneberg"`,
  for raw activity that hasn’t been run through
  [`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)/[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)),
  fits each bout, and applies the Winnebeck/Hammad quality filter
  (`pearson_r`, `p_value`, offset bounds).
  [`study_lids_metrics()`](https://zeitr.circadia-lab.uk/reference/study_lids_metrics.md)
  is the batch/`syncR::sync()`-ready counterpart, summarising each
  participant’s quality-filtered bouts (median +/- IQR period,
  amplitude, offset, slope) into one row, alongside
  [`study_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/study_sleep_metrics.md)
  and
  [`study_summary()`](https://zeitr.circadia-lab.uk/reference/study_summary.md).

  Ported from a prototype R notebook draft by Mario Leocadio-Miguel
  (itself adapted from an older MATLAB script); one bug in that draft is
  fixed here – its bout-fusing step was truncated mid-statement
  (`fused.a...`) and would not have run as written; see
  [`?detect_lids_bouts`](https://zeitr.circadia-lab.uk/reference/detect_lids_bouts.md)
  for the full reimplementation (`.fuse_bouts()`). Not yet validated
  against `pyActigraphy`’s `LIDS` class or an external reference dataset
  – treat results accordingly until a parity check is run.

## zeitR 0.1.6 (2026-07)

### ✨ New features

- New
  [`read_axivity()`](https://zeitr.circadia-lab.uk/reference/read_axivity.md)
  – bridges [axR](https://github.com/circadia-bio/axR)’s
  `axivity_read_cwa()` raw per-sample output (AX3/AX6 `.cwa` files) into
  the same 9-column epoch-level shape
  [`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md)
  produces, converting raw triaxial acceleration to PIM/TAT/ZCM via
  [`compute_activity_counts()`](https://zeitr.circadia-lab.uk/reference/compute_activity_counts.md)
  (GT3X+ filter preset, `0.25`-`2.5` Hz, since no validated
  Axivity-specific preset exists). Wired into
  `read_actigraphy(device = "axivity")`. `axR` added to `Suggests` now
  that it’s published on r-universe. Treat `activity`/`ZCMn` from this
  path as an unvalidated approximation – see
  [`?read_axivity`](https://zeitr.circadia-lab.uk/reference/read_axivity.md)
  for the full caveats (sampling-rate detection, time-zone re-labelling
  instead of conversion, and what specifically hasn’t been checked
  against a reference like GGIR).
- New
  [`compute_sri()`](https://zeitr.circadia-lab.uk/reference/compute_sri.md)
  — Sleep Regularity Index (Phillips et al. 2017), ported from **Fix
  30** of the Python reference pipeline (`SRI_vallim`). Derives
  sleep/wake directly from the epoch-level `state` column zeitR’s own
  pipelines already produce, rather than from a pyActigraphy scoring
  algorithm – Julia’s concordance analysis against manual reference
  scoring showed this is substantially more accurate than
  Sadeh/Cole-Kripke/ Roenneberg/Scripps (ICC 0.82 vs 0.19-0.67, N=404).
  Off-wrist gaps ≤30 min are interpolated; longer gaps are excluded from
  the day-to-day comparison rather than counted as a mismatch.
- [`compute_npcra()`](https://zeitr.circadia-lab.uk/reference/compute_npcra.md)
  gains a `trim_to_d1` argument (default `TRUE`): the recording is now
  trimmed to start at 00:00 of D+1 before computing IS, IV, RA, L5, and
  M10, matching the Python reference pipeline’s convention (which never
  spans the raw, typically fractional, recording length – the previous
  default here reported e.g. `n_days = 7.5` where Python reports a clean
  `7`). Set `trim_to_d1 = FALSE` for the previous behaviour. Does not
  replicate Python’s separate 30-min-threshold rule for the M10/L5
  windows specifically – only the D+1 window start.

### 🐛 Bug fixes

- **Fix 25 / Fix 26c interaction**: `.recover_fragmented_episodes()`
  ([`classify_sleep_episodes()`](https://zeitr.circadia-lab.uk/reference/classify_sleep_episodes.md)’s
  Fix 26c step) had no knowledge of the boundary Fix 25 uses to exclude
  episodes truncated by the end of the recording. When Fix 25 correctly
  excluded such an episode, its date became “uncovered”, and the Fix 26c
  recovery scan would reconstruct the *same* episode from the *same* raw
  epochs – silently undoing Fix 25 and producing a biologically
  implausible extra main night (e.g. 8 main nights on a 7-day
  recording). Root-caused while investigating Julia’s report of this
  exact symptom (matching the Python pipeline’s now-fixed Fix 29f).
  Fixed by giving `.recover_fragmented_episodes()` the same
  last-day-noon boundary and skipping recovery for any candidate that
  would itself start at/after noon on the recording’s last calendar day.

### 🧪 Tests

- New `test-read-axivity.R`:
  [`read_axivity()`](https://zeitr.circadia-lab.uk/reference/read_axivity.md)’s
  bridging logic (column shape, per-epoch light/int_temp averaging,
  epoch-start datetimes, tz re-labelling vs shifting, metadata assembly,
  dominant-sample-rate detection with an outlier warning, and an exact
  match against a direct
  [`compute_activity_counts()`](https://zeitr.circadia-lab.uk/reference/compute_activity_counts.md)
  call on the same input) – `axivity_read_cwa()` itself is mocked
  throughout via
  [`testthat::local_mocked_bindings()`](https://testthat.r-lib.org/reference/local_mocked_bindings.html),
  so these don’t need a real `.cwa` binary fixture. Skipped via
  `skip_if_not_installed()` for `axR`/`mrpheus`/`withr` where needed.
  `withr` added to `Suggests`.
- `test-fix26c.R`: regression test reproducing the Fix 25 / Fix 26c
  interaction above – a short evening sleep-like run right at the file’s
  end (mirroring the notebook’s ID_0138 case) is no longer recovered.
- `test-npcra.R`: `trim_to_d1` default behaviour (exact day removed,
  `n_days` unaffected numerically on the repeating fixture), and the
  \<2-epoch fallback-with-warning path.
- New `test-sri.R`:
  [`compute_sri()`](https://zeitr.circadia-lab.uk/reference/compute_sri.md)
  on perfectly regular (SRI = 100) and perfectly inverted (SRI = -100)
  synthetic patterns, off-wrist gap interpolation vs exclusion at the
  30-min boundary, `zeitr_result`/bare data frame input, the \<24h
  fallback-with-warning path, and `.interpolate_short_gaps()` directly
  (ffill, start-of-vector bfill, gap spanning the whole vector).

## zeitR 0.1.5 (2026-07)

### ✨ New features

- New
  [`pa_equations()`](https://zeitr.circadia-lab.uk/reference/pa_equations.md),
  [`estimate_ee()`](https://zeitr.circadia-lab.uk/reference/estimate_ee.md),
  [`classify_pa_intensity()`](https://zeitr.circadia-lab.uk/reference/classify_pa_intensity.md),
  and
  [`classify_pa_counts()`](https://zeitr.circadia-lab.uk/reference/classify_pa_counts.md)
  – physical activity intensity classification from ActTrust(R)/GT3X+
  activity counts, porting the published cut-points and MET-estimation
  equations from Batista et al. (2026, PLoS ONE
  <https://doi.org/10.1371/journal.pone.0348631>). Extends zeitR beyond
  sleep staging into the other half of the 24h rest-activity cycle:
  light, moderate, vigorous, and very vigorous PA bands from hip- or
  wrist-worn ActTrust(R) or GT3X+ counts. Only the published
  coefficient/cut-point table is ported (not the original study’s
  calorimetry-fitted [`lm()`](https://rdrr.io/r/stats/lm.html)/`msm`
  pipeline, which needs data no zeitR user will have) – see
  [`?pa_equations`](https://zeitr.circadia-lab.uk/reference/pa_equations.md)
  for the equation-set design and important generalisability caveats:
  single lab-treadmill validation study (N=56, healthy adults 18-35),
  GT3X+ (hip) cut-points differ 2-65% from prior published GT3X+ studies
  which themselves differ from each other by 16-39%, and that spread
  reflects differing modelling approaches (two-regression, ANN, and this
  paper’s linear model) as well as sample – diagnosable against Sasaki
  et al.’s transparent two-regression model, not diagnosable against
  Santos-Lozano et al.’s ANN, which has no inspectable coefficients to
  compare against.

- New
  [`compute_activity_counts()`](https://zeitr.circadia-lab.uk/reference/compute_activity_counts.md)
  – converts raw triaxial acceleration (`x`/`y`/`z` + `sampling_rate`)
  into epoch-level PIM/TAT/ZCM activity counts, for devices or pipelines
  that only provide raw samples rather than onboard-computed counts
  (e.g. [`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md)’s
  `activity`/`ZCMn` columns). Filtering reuses
  [`mrpheus::remove_dc()`](https://mrpheus.circadia-lab.uk/reference/remove_dc.html)
  /
  [`mrpheus::bandpass_filter()`](https://mrpheus.circadia-lab.uk/reference/bandpass_filter.html)
  – the same zero-phase Butterworth implementation already validated as
  part of mrpheus’s YASA-parity PSG pipeline – rather than a new,
  unvalidated filter written from scratch. This makes `mrpheus` a
  cross-package `Suggests` dependency for this function specifically
  (precedented by hypnoR, which already `Suggests` both `mrpheus` and
  `zeitR`); the rest of zeitR stays independent of it. The epoch-level
  PIM/TAT/ZCM aggregation logic itself has no reference implementation
  to check against – no raw-to-counts converter exists in
  `condor_pipeline`/`circadiaBase_Docker` (only already-epoched data),
  and Condor’s/ActiGraph’s exact onboard thresholds are proprietary –
  see
  [`?compute_activity_counts`](https://zeitr.circadia-lab.uk/reference/compute_activity_counts.md)
  for what’s built from the general processing description in Batista et
  al. (2026, PLoS ONE) vs. what’s an open/tunable parameter
  (`zcm_threshold`, `tat_threshold`). Default band-pass cutoffs are
  ActTrust-style (`0.5`-`2.7` Hz); pass `filter_low`/`filter_high` for
  GT3X+-style (`0.25`-`2.5` Hz) processing instead.

### 🚀 CI

- Fixed mrpheus dependency resolution in CI: added
  `Additional_repositories: https://circadia-bio.r-universe.dev` and a
  repo-root `.Rprofile` setting `options(repos = ...)` directly,
  matching hypnoR’s existing setup for the same mrpheus/zeitR pairing.
  `Additional_repositories` alone doesn’t automatically wire into
  `pak`’s dependency resolution – it’s mainly a
  documentation/NOTE-suppression field – so the `.Rprofile` is what
  actually makes `mrpheus` resolvable for
  [`compute_activity_counts()`](https://zeitr.circadia-lab.uk/reference/compute_activity_counts.md).

### 🐛 Bug fixes

- `.zero_crossing_indicator()` (internal,
  [`compute_activity_counts()`](https://zeitr.circadia-lab.uk/reference/compute_activity_counts.md))
  called `rep(FALSE, n - 1)` with `n = 0` for a zero-length input, and
  [`rep()`](https://rdrr.io/r/base/rep.html) rejects a negative `times`
  argument. Surfaced while writing the degenerate-length test case, not
  by any real recording. Fixed with an early
  `if (n < 2L) return(logical(0))` guard.

### 📚 Documentation

- New
  [`vignette("physical-activity")`](https://zeitr.circadia-lab.uk/articles/physical-activity.md)
  – walks through
  [`pa_equations()`](https://zeitr.circadia-lab.uk/reference/pa_equations.md),
  [`estimate_ee()`](https://zeitr.circadia-lab.uk/reference/estimate_ee.md),
  [`classify_pa_intensity()`](https://zeitr.circadia-lab.uk/reference/classify_pa_intensity.md),
  and
  [`classify_pa_counts()`](https://zeitr.circadia-lab.uk/reference/classify_pa_counts.md)
  using the bundled ActTrust recording’s real PIM counts. Flags the
  extrapolation from treadmill-fitted equations to free-living data
  inline rather than burying it in a caveats section, and points to
  [`?pa_equations`](https://zeitr.circadia-lab.uk/reference/pa_equations.md)
  for the full generalisability discussion instead of duplicating it.
- New
  [`vignette("raw-accelerometry")`](https://zeitr.circadia-lab.uk/articles/raw-accelerometry.md)
  –
  [`compute_activity_counts()`](https://zeitr.circadia-lab.uk/reference/compute_activity_counts.md)
  on a simulated raw triaxial recording (quiet segment, then clear 1 Hz
  movement), covering ActTrust- vs. GT3X+-style filter cutoffs, tuning
  `zcm_threshold`/`tat_threshold` against sensor noise floor, and
  closing the full loop into
  [`estimate_ee()`](https://zeitr.circadia-lab.uk/reference/estimate_ee.md)/[`classify_pa_intensity()`](https://zeitr.circadia-lab.uk/reference/classify_pa_intensity.md)
  so the raw signal and PA-intensity vignettes read as one continuous
  pipeline rather than two unrelated features.
- README: added the `r-universe` badge and recommended
  `install.packages(..., repos = c("https://circadia-bio.r-universe.dev", ...))`
  install path (GitHub `pak` install kept as the dev-version fallback,
  matching hypnoR); filled in ten previously undocumented `Features`
  bullets
  ([`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md),
  [`prepare_actigraphy()`](https://zeitr.circadia-lab.uk/reference/prepare_actigraphy.md),
  the four PA-intensity functions, the four actogram-plotting functions)
  that were already shipped and in the pkgdown reference index but never
  listed; added the PA-intensity MET-band table to Computed Variables;
  brought the Project Structure tree and Dependencies table (now with an
  Imports/Suggests `Type` column) up to date with everything actually in
  `DESCRIPTION`.

### 🧪 Tests

- `test-raw-accelerometry.R`:
  `.zero_crossing_indicator()`/`.epoch_rowsum()` directly, plus
  [`compute_activity_counts()`](https://zeitr.circadia-lab.uk/reference/compute_activity_counts.md)
  input validation, the trailing-incomplete-epoch warning, a flat
  zero-motion signal (all metrics exactly zero), a clean 1 Hz sinusoid
  (ZCM ~ 2 x frequency x `epoch_sec` on interior epochs), PIM/TAT
  increasing with amplitude, and `metrics` subsetting. Skipped via
  `skip_if_not_installed("mrpheus")` where that dependency is needed.

## zeitR 0.1.4 (2026-07)

### ✨ New features

- New
  [`study_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/study_sleep_metrics.md)
  – batch wrapper computing
  [`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
  and
  [`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
  across every participant in a
  [`run_pipeline_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_batch.md)/[`run_pipeline_native_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native_batch.md)
  result, stacked into one tibble with a `participant_id` column – the
  sleep-timing/chronotype (CPD, MSF/MSW, social jetlag) counterpart to
  [`study_summary()`](https://zeitr.circadia-lab.uk/reference/study_summary.md)
  (NPCRA/activity-rhythm variables). Closes a gap found while checking
  `syncR::sync()` compatibility:
  [`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
  and
  [`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
  each return a single named list per participant with no participant
  identifier and no batch equivalent, unlike
  [`study_summary()`](https://zeitr.circadia-lab.uk/reference/study_summary.md),
  so there was previously no way to get these metrics into the
  one-row-per-participant shape `sync()` expects without writing manual
  glue code per study.

### 📊 Visualisation

- [`actogram_colours()`](https://zeitr.circadia-lab.uk/reference/actogram_colours.md):
  swapped the default `"wake"` and `"off-wrist"` colours (wake is now
  the warm terracotta `#C25E2A`; off-wrist is now the neutral sand
  `#D9C8A0`). Affects the default palette used by
  [`plot_actogram()`](https://zeitr.circadia-lab.uk/reference/plot_actogram.md),
  [`plot_actogram_double()`](https://zeitr.circadia-lab.uk/reference/plot_actogram_double.md),
  and
  [`plot_actogram_activity()`](https://zeitr.circadia-lab.uk/reference/plot_actogram_activity.md)
  whenever `colours` is not supplied explicitly.
- [`plot_actogram_activity()`](https://zeitr.circadia-lab.uk/reference/plot_actogram_activity.md)
  gains a `log_scale` argument. When `TRUE`, applies a
  [`log1p()`](https://rdrr.io/r/base/Log.html) transform to the activity
  signal before capping and normalising bar heights, compressing the
  dynamic range so structure among low-to-moderate activity epochs is
  easier to see against a right-skewed raw signal (occasional high
  bursts no longer dominate the visible range). Default `FALSE`
  preserves the existing linear-scale behaviour exactly.

### 🚀 Performance

- Removed `.adaptive_median_filter()`, a dead pure-R fallback in
  `sleep_periods.R` that was superseded by
  `adaptive_median_filter_cpp()` and never called. No behaviour change.
- Removed four further dead internal helpers from `utils.R` –
  `zero_sequences()`, `quantile_filter()`, `max_filter()`,
  `min_filter()` – confirmed unreferenced anywhere in `R/`, `dev/`, or
  the test suite; all superseded by direct calls to the corresponding
  Rcpp functions. No behaviour change. (`rolling_apply()` was initially
  removed too, but is kept – it’s the R reference implementation used by
  the `rolling_max_cpp()`/`rolling_min_cpp()` parity tests in
  `test-crespo-cpp-parity.R`.)
- [`run_pipeline_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_batch.md)
  and
  [`run_pipeline_native_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native_batch.md)
  gain a `parallel` argument. When `TRUE`, files are processed
  concurrently via
  [`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html)
  under whatever
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  the caller has set
  (e.g. `future::plan(future::multisession(workers = 4))`). Falls back
  to sequential processing with a warning if `future.apply` is not
  installed. Default remains `FALSE` (sequential), so existing code is
  unaffected. `future` and `future.apply` added to `Suggests`.

### 🐛 Bug fixes

- `zeitr_abort()`, `zeitr_warn()`, and `zeitr_inform()` (internal
  message wrappers around
  [`cli::cli_abort()`](https://cli.r-lib.org/reference/cli_abort.html)/`cli_warn()`/`cli_inform()`)
  did not forward `.envir`, so [glue](https://glue.tidyverse.org/)-style
  interpolation of a variable local to the *calling* function
  (e.g. `{.val {missing_cols}}`) silently failed with “object not found”
  instead of producing the intended message. Fixed by defaulting
  `.envir = parent.frame()` in all three wrappers; no call sites needed
  to change.
- Five call sites (three in `plot_actogram.R`, two in `export.R`) passed
  their error/warning message as multiple separate comma-delimited
  string arguments instead of one string. R does not auto-concatenate
  adjacent string literals, so the extra arguments were passed through
  to
  [`cli::cli_abort()`](https://cli.r-lib.org/reference/cli_abort.html)/`cli_warn()`
  as unnamed condition data, which `rlang` rejects (“Conditions must
  have named data fields”). Fixed by merging each into a single message
  string. None of the five affected error paths had previously been
  exercised by a test.
- [`plot_actogram()`](https://zeitr.circadia-lab.uk/reference/plot_actogram.md)
  and
  [`plot_actogram_double()`](https://zeitr.circadia-lab.uk/reference/plot_actogram_double.md):
  the epoch at `mins_since_midnight = 0` (midnight, i.e. the first epoch
  of every calendar day) had its left half clipped by `geom_tile()`’s
  centred tile extending outside a hard `scale_x_continuous()` limit of
  exactly `0`, producing a silently dropped/incomplete tile at the start
  of every row (and a `ggplot2` “missing values” warning once actually
  rendered under test). Fixed by widening the x-axis limits by half an
  epoch on each side in both
  [`plot_actogram()`](https://zeitr.circadia-lab.uk/reference/plot_actogram.md)
  and the shared `.scale_x_double()` helper (also used by
  [`plot_actogram_activity()`](https://zeitr.circadia-lab.uk/reference/plot_actogram_activity.md),
  which was not affected by the clipping itself since it uses
  `geom_rect()` with explicit epoch boundaries rather than a centred
  tile).
- [`export_hypnogram()`](https://zeitr.circadia-lab.uk/reference/export_hypnogram.md):
  `is.list(result)` doesn’t exclude a bare tibble (tibbles are lists
  too), so `result$subject_id` on a bare tibble without that column
  triggered a spurious “Unknown or uninitialised column” warning. Fixed
  by adding the missing `!is.data.frame(result)` guard, matching the
  pattern already used two lines above it in the same function.
- [`export_hypnogram()`](https://zeitr.circadia-lab.uk/reference/export_hypnogram.md):
  when `ZCMn` is absent (a documented, valid use case), `zcm` is `NULL`
  and `zcm == 0` evaluates to `logical(0)`, which
  [`dplyr::case_when()`](https://dplyr.tidyverse.org/reference/case-and-replace-when.html)
  cannot recycle against the other length-n conditions – a hard error,
  not just a warning. Fixed by precomputing a proper full-length
  `zcm_is_zero` vector before the
  [`case_when()`](https://dplyr.tidyverse.org/reference/case-and-replace-when.html)
  call.
- [`compute_npcra()`](https://zeitr.circadia-lab.uk/reference/compute_npcra.md):
  the same `is.null(x$col)`-on-a-tibble pattern for the optional `state`
  column triggered a spurious “Unknown or uninitialised column” warning
  whenever `state` was absent (also a documented, valid use case). Fixed
  with `"state" %in% names(epochs)`.
- `plot_actogram.R`’s internal `.actogram_title()` had the same
  [`is.list()`](https://rdrr.io/r/base/list.html)/bare-tibble gap as
  [`export_hypnogram()`](https://zeitr.circadia-lab.uk/reference/export_hypnogram.md),
  not yet triggered by any existing test but the same latent risk. Fixed
  proactively for consistency with the pattern used elsewhere.

### 🧪 Tests

- `test-batch-helper.R`: `.run_pipeline_over_files()` – sequential
  success, partial-failure skip-with-warning, all-failing batch returns
  empty list, parallel dispatch via `future_lapply()` (skipped if
  `future.apply` is not installed), sequential fallback when
  `future.apply` is unavailable (skipped if it is installed), and a
  regression guard confirming both exported batch wrappers still default
  to `parallel = FALSE`.
- `test-actogram-snapshots.R`: visual regression snapshots (`vdiffr`,
  skipped if not installed) for
  [`plot_actogram()`](https://zeitr.circadia-lab.uk/reference/plot_actogram.md),
  [`plot_actogram_double()`](https://zeitr.circadia-lab.uk/reference/plot_actogram_double.md),
  [`plot_actogram_activity()`](https://zeitr.circadia-lab.uk/reference/plot_actogram_activity.md),
  and
  [`plot_actogram_activity()`](https://zeitr.circadia-lab.uk/reference/plot_actogram_activity.md)
  with a custom `activity_cap_quantile` and with `log_scale = TRUE`,
  using a deterministic synthetic 2-day fixture. Also covers
  (independently of `vdiffr`) missing-column errors for all three
  functions, the missing-`activity_col` error, acceptance of a
  `zeitr_result` list as well as a bare tibble, and that
  `log_scale = FALSE` is byte-identical to the pre-`log_scale`
  behaviour. `vdiffr` added to `Suggests`.
- New test files bringing five previously 0%-covered files up to full or
  near-full coverage: `test-circ-utils.R`, `test-export-hypnogram.R`,
  `test-npcra.R`, `test-study-summary.R`, `test-read-actigraphy.R`. Also
  `test-utils.R`, covering edge-case branches in `norm_01()`,
  `zero_prop()`, `ashman_d()`, and `%||%` not guaranteed to be hit by
  ordinary pipeline data. Overall coverage moved from 81.1% to 87.9%.
  `covr` added to `Suggests`.
- `test-study-sleep-metrics.R`: synthetic multi-participant coverage for
  [`study_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/study_sleep_metrics.md)
  – both metric sets present with correct `n_overall`/`n_wd`/`n_fd`
  counts, holiday forwarding shifting a night between the
  workday/free-day groups, per-participant `holidays`/`free_days`
  fallback vs a study-level override, a participant whose
  [`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
  call fails (no free days) while
  [`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
  still succeeds for the same participant (only the failing metric set
  is `NA`-filled), skipping non-`zeitr_result` entries, the
  empty/all-invalid-batch paths, and the `subject_id`-missing fallback
  to the list name.

## zeitR 0.1.3 (2026-07)

### 📊 Visualisation

- [`plot_actogram()`](https://zeitr.circadia-lab.uk/reference/plot_actogram.md)
  – single-column raster actogram. One row per calendar day, time-of-day
  on the x-axis, filled by sleep/wake state. Oldest day at the top,
  following standard chronobiology convention. Equivalent to the ad-hoc
  ggplot2 code in the single-recording vignette but packaged as a
  reusable function with consistent defaults.
- [`plot_actogram_double()`](https://zeitr.circadia-lab.uk/reference/plot_actogram_double.md)
  – classic double-plotted actogram. Each recording day appears twice:
  in the left column of its own row (x = 00:00 to 24:00) and in the
  right column of the row above (x = 24:00 to 48:00). Circadian phase
  drift is visible as a diagonal band across consecutive rows. A dashed
  vertical line marks the 24 h column boundary.
- [`plot_actogram_activity()`](https://zeitr.circadia-lab.uk/reference/plot_actogram_activity.md)
  – double-plotted actogram with activity bars. Same row structure as
  [`plot_actogram_double()`](https://zeitr.circadia-lab.uk/reference/plot_actogram_double.md)
  but each epoch is drawn as a vertical bar whose height is proportional
  to the raw ZCMn activity count. Bars are coloured by sleep/wake state
  so activity intensity and state classification are read
  simultaneously. A 99th-percentile cap on bar heights prevents outlier
  bursts from compressing the rest of the range; a thin baseline stub
  keeps zero-activity epochs (sleep, off-wrist) faintly visible.
- [`actogram_colours()`](https://zeitr.circadia-lab.uk/reference/actogram_colours.md)
  – exported helper returning the named hex colour vector used as the
  default palette across all three actogram functions. Pass the result
  to any `colours` argument to inspect or partially override defaults.
- All three functions accept a `zeitr_result` list or a bare tibble with
  `datetime` and `state` columns. `ggplot2` remains in `Suggests`; a
  clear error is thrown if it is not installed.

### 🚀 Performance

- `rolling_median_prepadded_cpp()` added to `src/rolling_filters.cpp`.
  Replaces the `RcppRoll` / `zoo` / `vapply` fallback chain in
  `.estimate_sleep_padded()` with a single direct Rcpp call. Off-wrist
  sleep estimation is now unconditionally fast (O(n \* win) in C++)
  regardless of which optional packages are installed. `zoo` removed
  from Imports; `RcppRoll` removed from Suggests.
- [`check_consistency()`](https://zeitr.circadia-lab.uk/reference/check_consistency.md)
  vectorised. Two O(n) R `for` loops replaced with
  [`which()`](https://rdrr.io/r/base/which.html) calls. No behaviour
  change.

### 📅 Free-day classification

- New `free_days` parameter on
  [`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md),
  [`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md),
  and
  [`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md).
  Replaces the hardcoded Saturday + Sunday with any combination of days
  (English names or ISO integers 1–7). Default is
  `c("Saturday", "Sunday")`. Enables non-standard schedules such as
  Friday–Saturday weekends or compressed work weeks.
- `holidays` now accepts `"DD-MM"` strings for fixed-date annual
  holidays (e.g. `"25-12"` for Christmas) in addition to `Date` objects
  and `"YYYY-MM-DD"` strings. All three forms can be mixed in the same
  vector.
- [`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
  and
  [`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
  are now S3 generics. Passing a `zeitr_result` directly auto-forwards
  `result$holidays` and `result$free_days` — no need to repeat them
  manually.
- A warning is emitted when `holidays = NULL`; suppress with
  `options(zeitR.no_holidays_warn = FALSE)`.

### 🐛 Bug fixes

- Free-day detection was broken on non-English locales
  ([`weekdays()`](https://rdrr.io/r/base/weekday.POSIXt.html) returns
  `"sabado"` on `pt_BR`). Fixed by using the locale-independent ISO 8601
  weekday number.
- MSF and MSW now use the circular mean, matching the fix29 notebook.
  Plain mean gives wrong results when mid-sleep wraps midnight.
- [`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
  now drops episodes starting after noon on the last recording day
  (truncated by end of file), matching fix29’s filter.

### 🧪 Tests

- Free-day classification tests (`test-free-days.R`):
  `.parse_free_days()` input validation (English names, ISO integers,
  case insensitivity, range errors), `.is_free_day()` locale-independent
  weekday detection (vectorised over a full week, default and custom
  schedules), all three holiday input forms (`Date`, `"YYYY-MM-DD"`,
  `"DD-MM"`) including mixed-form vectors, year-specificity of
  `"YYYY-MM-DD"` vs recurrence of `"DD-MM"`, the
  `zeitR.no_holidays_warn` option, and `zeitr_result` S3 dispatch
  forwarding `free_days` and `holidays` to both
  [`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
  and
  [`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md).

### 📚 Documentation

- New vignette
  [`vignette("actogram")`](https://zeitr.circadia-lab.uk/articles/actogram.md)
  – covers all three plot functions and
  [`actogram_colours()`](https://zeitr.circadia-lab.uk/reference/actogram_colours.md):
  single-column vs double-plotted vs activity-bar formats, colour
  customisation, `date_label_every`, extending the returned `ggplot`
  object with additional layers, and working with bare tibbles instead
  of a `zeitr_result`.

## zeitR 0.1.2 (2026-07)

### 🌙 Vallim native pipeline

- [`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)
  — full single-recording pipeline using the Vallim (JRSV) rule set
  developed by Julia Ribeiro da Silva Vallim. Replaces Condor’s
  `nights_df` classification with a 7-step adaptive rule set: Fix 25
  (truncated episode exclusion), Fix 26a (adaptive nocturnal window
  inferred from wrist temperature and ambient light), Fix 26b
  (sleep-date collision resolution), Fix 26c (fragmented episode
  recovery), Fix 29 (14–16 h episode splitting), Rules 3–5
  (main/secondary classification), Rule 6 (longest-main-per-date
  selection), Rule 7 (days-without-main exclusion). Validated against
  the Python reference pipeline: all 52 main nights on the ActTrust
  validation recording classified identically.
- [`run_pipeline_native_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native_batch.md)
  — directory-level wrapper for
  [`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md).
- [`extract_sleep_episodes()`](https://zeitr.circadia-lab.uk/reference/extract_sleep_episodes.md)
  — convert a CSPD-scored epoch table into a per-episode tibble with
  WASO, SOL, SOI, TST, NW, EFF derived via Cole-Kripke epoch scoring.
- [`classify_sleep_episodes()`](https://zeitr.circadia-lab.uk/reference/classify_sleep_episodes.md)
  — apply the full JRSV rule set to a raw episode table and return
  `sleep_type` (`"main"` / `"secondary"`).
- [`circ_mean_h()`](https://zeitr.circadia-lab.uk/reference/circ_mean_h.md)
  — unit-circle circular mean for clock-time variables (Fix 20; handles
  midnight wrap correctly).
- [`circ_sd_h()`](https://zeitr.circadia-lab.uk/reference/circ_sd_h.md)
  — circular SD using mean resultant length formula (Fix 24; invariant
  to the wrap point).
- Julia Ribeiro da Silva Vallim (ORCID 0000-0001-8708-8479) added as
  author in DESCRIPTION, `_pkgdown.yml`, and pipeline documentation.

### 🚀 Performance (Rcpp)

- Five rolling filters replaced with Rcpp implementations — **215×**
  speedup on a 40,000-epoch recording (12 GB → 312 KB memory):
  `rolling_median_cpp`, `rolling_mean_cpp`, `rolling_var_cpp`,
  `rolling_zero_prop_cpp`, `rolling_quantile_cpp`.
- `diff5()` five-point stencil derivative ported to `diff5_cpp` —
  replaces the interior R `for` loop; used in off-wrist temperature
  derivative.
- [`score_epochs_cole_kripke()`](https://zeitr.circadia-lab.uk/reference/score_epochs_cole_kripke.md)
  ported to `score_epochs_cole_kripke_cpp` — single-pass O(n)
  convolution replaces 17 vectorised R additions.
- `Rcpp (>= 1.0.0)` added to `Imports` and `LinkingTo`.

### 🔧 Other changes

- `mclust` moved from `Suggests` to `Imports`; the GMM fallback warning
  is no longer emitted for standard ActTrust recordings.
- `_pkgdown.yml` updated with all new exports and Julia Vallim
  authorship.

### 🧪 Tests

- Rcpp rolling filter parity tests (`test-rolling-filters-parity.R`):
  all five filters, `diff5_cpp`, and `score_epochs_cole_kripke_cpp`
  validated against R reference implementations.
- Vallim pipeline classification parity tests (`test-vallim-parity.R`):
  episode count, sleep-date coverage, and `sleep_type` classification
  locked against `inst/extdata/vallim_nights.csv` (generated by
  `dev/parity_vallim.py` on the Python reference pipeline).

### 📤 hypnoR export

- [`export_hypnogram()`](https://zeitr.circadia-lab.uk/reference/export_hypnogram.md)
  — converts a `zeitr_result` into the tidy hypnogram format expected by
  `hypnoR`. Stage mapping: `ZCMn == 0` within sleep epochs becomes
  `"Quiet sleep"`; non-zero sleep activity becomes `"Sleep"`; wake and
  off-wrist epochs become `"W"`. `subject_id` is inferred automatically
  from `result$subject_id` (set by the pipeline from the filename stem)
  and can be overridden with an explicit argument. Works in both
  single-file and batch contexts.

### 📋 Sleep summary metrics

- [`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
  — per-night sleep metrics split by day type (overall / workday / free
  day): SOL, TST, TBT, WASO, sleep efficiency, sleep onset, get-up time,
  mid-sleep, and within-person SDs. Column names and arithmetic mirror
  Julia Vallim’s
  [`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)
  in `pipeline_functions_fix27.py`.
- [`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
  — CPD, MSW, MSF, MSFsc, SJL, and signed SJLa. Ports
  `nights_to_cpd_df()` and
  [`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
  from the same reference. Both functions accept a `holidays` argument
  for country-specific public holidays beyond weekends.

### 🐛 Bug fixes

- **Fix 26c (fragment recovery)**: two bugs closed.
  1.  The Python reference pipeline silently disabled temperature- and
      light-based gap merging due to column name mismatches
      (`'TEMPERATURE'` / `'LIGHT'` vs the actual `'int_temp'` /
      `'light'` columns in ActTrust data). `pipeline_functions_fix27.py`
      patched; `inst/extdata/vallim_nights.csv` regenerated against the
      corrected Python output and re-verified: 52/52 main nights, all
      sleep dates and classifications match R.
  2.  R’s Fix 26c was using the period-level CSPD `state` column to
      detect sleep runs within the recovery window, which incorrectly
      treated entire 19+ h CSPD periods as a single sleep run. Fix 26c
      now uses Cole-Kripke epoch scoring on `ZCMn` to determine
      sleep/wake within the candidate window, matching the intended
      behaviour. R is the reference implementation for Fix 26c.
- **`offwrist_refiner.R`**: fixed scalar `FALSE` assignment to a 0-row
  data frame (`$valley_peak <- rep(FALSE, nrow(...))`) that caused a
  crash on recordings with no valid off-wrist candidates.

### 🚀 Performance (Rcpp) — continued

- **Crespo MSP hot paths** — five additional Rcpp ports eliminating the
  remaining R `for` loops and `vapply` calls in the main sleep detector:
  - `rolling_max_cpp` / `rolling_min_cpp` — replace
    `rolling_apply(max/min)` in `.morphological_open_close()`.
  - `zero_mitigation_cpp` — replaces the zero-run mitigation `for` loop
    (pass 1 of `.crespo_msp()`).
  - `mark_invalid_zeros_cpp` — replaces the invalid-zero marking `for`
    loop (pass 2 of `.crespo_msp()`).
  - `adaptive_median_filter_cpp` — replaces the variable-window adaptive
    median `for` loop in both `.crespo_msp()` and `.crespo_nap_msp()`;
    this was the single hottest loop in the pipeline.
  - The coarse median filter `vapply` in `.crespo_msp()` now reuses the
    existing `rolling_median_cpp` with constant padding.

### 🧪 Tests

- Crespo C++ parity tests (`test-crespo-cpp-parity.R`):
  `rolling_max_cpp`, `rolling_min_cpp`, `zero_mitigation_cpp`,
  `mark_invalid_zeros_cpp`, `adaptive_median_filter_cpp`, morphological
  close/open pair, and end-to-end epoch count lock on `input1.txt`.
- Fix 26c regression test (`test-fix26c.R`): synthetic 1-min epoch
  recording with a bloated 19 h CSPD `state = 1` period containing two
  Cole-Kripke sleep runs (3 h + 5.5 h) separated by a warm/dark wake
  gap. Asserts that `.recover_fragmented_episodes()` merges the
  CK-derived runs (TBT ~ 9 h) rather than the CSPD state period (TBT ~
  19 h). Runs on CI; no external data required.

------------------------------------------------------------------------

## zeitR 0.1.0 (2026-06)

### 🚀 Pipeline

- Full actigraphy pipeline validated epoch-for-epoch (`0 / 76,196`
  mismatches) against the Condor circadiaBase Python reference on an
  ActTrust recording:
  [`detect_offwrist_bimodal()`](https://zeitr.circadia-lab.uk/reference/detect_offwrist_bimodal.md),
  [`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md),
  [`detect_naps_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_naps_crespo.md)
  (faithful `nap_wrapper` port),
  [`compute_waso()`](https://zeitr.circadia-lab.uk/reference/compute_waso.md).
- [`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
  gains a `params` argument (default
  [`acttrust_params()`](https://zeitr.circadia-lab.uk/reference/acttrust_params.md));
  device-specific defaults are now consolidated in one place and can be
  overridden without touching individual detector calls. `wake_thresh`
  is removed from the
  [`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
  signature — it now lives in `params$waso$wake_thresh`.
- [`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
  gains a `quiet` argument to suppress the timestamp-issue warning
  (useful in batch and testing contexts).

### ✨ New functions

- [`acttrust_params()`](https://zeitr.circadia-lab.uk/reference/acttrust_params.md)
  — exported device parameter preset consolidating all ActTrust-specific
  defaults across off-wrist, sleep, nap, and WASO stages. Copy and
  modify to adapt the pipeline to other devices.
- [`label_states()`](https://zeitr.circadia-lab.uk/reference/label_states.md)
  — converts the integer `state` column to a human-readable ordered
  factor (`"wake"`, `"sleep"`, `"nap"`, `"off-wrist"`).

### 🧪 Tests

- End-to-end pipeline parity regression test (`test-pipeline-parity.R`):
  epoch-level state, per-layer counts, and nightly statistics locked
  against `python_output.csv` and `python_nights.csv`.
- CSPD refiner parity tests (`test-cspd-refiner-parity.R`): stage-1
  peak-valley length filter, stage-2 sleep-gap separation, full
  `.cspd_refine_periods` output, and bedtime/getuptime indices against
  Python intermediates.
- Sleep Crespo wiring test (`test-sleep-crespo-wiring.R`): isolates
  `detect_sleep_crespo(refine = TRUE)` against
  `cspd_refined_output.csv`.
- WASO parity tests (`test-waso-parity.R`): `.nights_df` boundaries,
  per-night statistics, and within-night epoch agreement on
  boundary-matched nights.

### 🌱 Initial release

- Full package scaffold:
  [`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md),
  [`read_actigraphy()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy.md),
  [`read_actigraphy_dir()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy_dir.md),
  [`prepare_actigraphy()`](https://zeitr.circadia-lab.uk/reference/prepare_actigraphy.md),
  [`check_consistency()`](https://zeitr.circadia-lab.uk/reference/check_consistency.md),
  [`score_epochs_cole_kripke()`](https://zeitr.circadia-lab.uk/reference/score_epochs_cole_kripke.md),
  [`compute_npcra()`](https://zeitr.circadia-lab.uk/reference/compute_npcra.md),
  [`study_summary()`](https://zeitr.circadia-lab.uk/reference/study_summary.md),
  [`run_pipeline_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_batch.md).
- Three vignettes: getting started, NPCRA, study-level analysis.
- pkgdown site with Bootstrap 5 and Circadia Lab branding.
