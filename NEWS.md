## zeitR 0.1.0  (2026-06)

### Pipeline

* Full actigraphy pipeline validated epoch-for-epoch (`0 / 76,196` mismatches)
  against the Condor circadiaBase Python reference on an ActTrust recording:
  `detect_offwrist_bimodal()`, `detect_sleep_crespo()`,
  `detect_naps_crespo()` (faithful `nap_wrapper` port), `compute_waso()`.
* `run_pipeline()` gains a `params` argument (default `acttrust_params()`);
  device-specific defaults are now consolidated in one place and can be
  overridden without touching individual detector calls. `wake_thresh` is
  removed from the `run_pipeline()` signature — it now lives in
  `params$waso$wake_thresh`.
* `run_pipeline()` gains a `quiet` argument to suppress the timestamp-issue
  warning (useful in batch and testing contexts).

### New functions

* `acttrust_params()` — exported device parameter preset consolidating all
  ActTrust-specific defaults across off-wrist, sleep, nap, and WASO stages.
  Copy and modify to adapt the pipeline to other devices.
* `label_states()` — converts the integer `state` column to a human-readable
  ordered factor (`"wake"`, `"sleep"`, `"nap"`, `"off-wrist"`).

### Tests

* End-to-end pipeline parity regression test (`test-pipeline-parity.R`):
  epoch-level state, per-layer counts, and nightly statistics locked against
  `python_output.csv` and `python_nights.csv`.
* CSPD refiner parity tests (`test-cspd-refiner-parity.R`): stage-1
  peak-valley length filter, stage-2 sleep-gap separation, full
  `.cspd_refine_periods` output, and bedtime/getuptime indices against Python
  intermediates.
* Sleep Crespo wiring test (`test-sleep-crespo-wiring.R`): isolates
  `detect_sleep_crespo(refine = TRUE)` against `cspd_refined_output.csv`.
* WASO parity tests (`test-waso-parity.R`): `.nights_df` boundaries,
  per-night statistics, and within-night epoch agreement on
  boundary-matched nights.

### Initial release

* Full package scaffold: `read_acttrust()`, `read_actigraphy()`,
  `read_actigraphy_dir()`, `prepare_actigraphy()`, `check_consistency()`,
  `score_epochs_cole_kripke()`, `compute_npcra()`, `study_summary()`,
  `run_pipeline_batch()`.
* Three vignettes: getting started, NPCRA, study-level analysis.
* pkgdown site with Bootstrap 5 and Circadia Lab branding.
