# Package index

## Import

Read raw actigraphy files into R.

- [`read_actigraphy()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy.md)
  : Read an actigraphy file into a zeitr_recording object
- [`read_actigraphy_dir()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy_dir.md)
  : Read all actigraphy files in a directory
- [`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md)
  : Read a Condor Instruments ActTrust actigraphy file

## Prepare & validate

Prepare recordings for analysis and check data quality.

- [`prepare_actigraphy()`](https://zeitr.circadia-lab.uk/reference/prepare_actigraphy.md)
  : Prepare a raw actigraphy tibble for analysis
- [`check_consistency()`](https://zeitr.circadia-lab.uk/reference/check_consistency.md)
  : Check actigraphy timestamps for consistency issues

## Detection

Off-wrist, sleep period, and epoch-level scoring algorithms.

- [`detect_offwrist_bimodal()`](https://zeitr.circadia-lab.uk/reference/detect_offwrist_bimodal.md)
  : Off-wrist detection using the Condor bimodal activity/temperature
  model
- [`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md)
  : Detect main sleep periods using the Crespo algorithm
- [`detect_naps_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_naps_crespo.md)
  : Detect secondary sleep periods (naps) using the Crespo nap algorithm
- [`score_epochs_cole_kripke()`](https://zeitr.circadia-lab.uk/reference/score_epochs_cole_kripke.md)
  : Score actigraphy epochs as wake or sleep using the Cole-Kripke
  algorithm
- [`compute_waso()`](https://zeitr.circadia-lab.uk/reference/compute_waso.md)
  : Compute WASO and nightly sleep statistics

## Circadian analysis

Non-parametric circadian rhythm variables.

- [`compute_npcra()`](https://zeitr.circadia-lab.uk/reference/compute_npcra.md)
  : Non-parametric circadian rhythm analysis (NPCRA)
- [`compute_waso()`](https://zeitr.circadia-lab.uk/reference/compute_waso.md)
  : Compute WASO and nightly sleep statistics

## Study-level

Summarise results across multiple participants.

- [`study_summary()`](https://zeitr.circadia-lab.uk/reference/study_summary.md)
  : Summarise a zeitr_study across participants

## Pipeline

Run the full analysis pipeline on one file or a directory.

- [`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
  : Run the full actigraphy sleep analysis pipeline
- [`run_pipeline_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_batch.md)
  : Run the pipeline on all files in a directory
- [`acttrust_params()`](https://zeitr.circadia-lab.uk/reference/acttrust_params.md)
  : ActTrust device parameter preset

## Utilities

Helper functions for working with pipeline output.

- [`label_states()`](https://zeitr.circadia-lab.uk/reference/label_states.md)
  : Convert integer epoch states to a labelled factor

## Package

Package-level documentation.

- [`zeitR`](https://zeitr.circadia-lab.uk/reference/zeitR-package.md)
  [`zeitR-package`](https://zeitr.circadia-lab.uk/reference/zeitR-package.md)
  : zeitR: Actigraphy Data Parsing and Analysis for R
