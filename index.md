# ⌚️ zeitR

**Actigraphy data parsing and analysis for R.**

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://zeitr.circadia-lab.uk/LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3)](https://www.r-project.org/)
[![R CMD
CHECK](https://github.com/circadia-bio/zeitR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/circadia-bio/zeitR/actions/workflows/R-CMD-check.yaml)
[![Status](https://img.shields.io/badge/status-early%20development-orange)](https://github.com/circadia-bio/zeitR)
[![pkgdown](https://img.shields.io/badge/docs-zeitr.circadia--lab.uk-F0A500)](https://zeitr.circadia-lab.uk)

------------------------------------------------------------------------

> \[!WARNING\] **zeitR is in early development and has not been formally
> validated.** The pipeline has been validated epoch-for-epoch against
> the Condor circadiaBase Python reference on an ActTrust recording, but
> the package has not undergone formal peer review. Verify outputs
> independently before using in any research context.

------------------------------------------------------------------------

## 📖 What is zeitR?

zeitR is an R package for importing, parsing, and analysing raw
actigraphy recordings from wrist-worn devices. It runs a full
rest-activity pipeline — off-wrist detection, sleep period
identification, WASO computation — and computes standard non-parametric
circadian rhythm variables (IS, IV, RA, L5, M10), returning tidy data
frames ready for downstream chronobiological analysis.

zeitR is designed to complement
[slumbR](https://github.com/circadia-bio/slumbR) in the Circadia Lab
ecosystem: slumbR handles sleep diary and questionnaire data, zeitR
handles the actigraphy side of a study, and both speak the same tidy,
pipeline-friendly R idioms.

------------------------------------------------------------------------

## ✨ Features

- 📥
  **[`read_actigraphy()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy.md)**
  — parse a raw device file into a `zeitr_recording` object with
  `$epochs` and `$metadata`
- 📂
  **[`read_actigraphy_dir()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy_dir.md)**
  — batch-read a whole directory of files into a `zeitr_study`
- 🔍
  **[`check_consistency()`](https://zeitr.circadia-lab.uk/reference/check_consistency.md)**
  — flag timestamp gaps, backward jumps, and firmware year artefacts
- 🦾
  **[`detect_offwrist_bimodal()`](https://zeitr.circadia-lab.uk/reference/detect_offwrist_bimodal.md)**
  — Condor bimodal activity/temperature off-wrist detection (three-stage
  refiner)
- 😴
  **[`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md)**
  — main sleep period detection (Crespo et al., 2012)
- 💤
  **[`detect_naps_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_naps_crespo.md)**
  — secondary sleep period (nap) detection (Crespo et al., 2012)
- ⏱️
  **[`score_epochs_cole_kripke()`](https://zeitr.circadia-lab.uk/reference/score_epochs_cole_kripke.md)**
  — epoch-level wake/sleep scoring (Cole & Kripke, 1992)
- 📊
  **[`compute_waso()`](https://zeitr.circadia-lab.uk/reference/compute_waso.md)**
  — nightly TBT, TST, WASO, SOL, SOI, number of awakenings, and sleep
  efficiency
- 📐
  **[`compute_npcra()`](https://zeitr.circadia-lab.uk/reference/compute_npcra.md)**
  — non-parametric circadian rhythm analysis: IS, IV, RA, L5, M10
- 🗂️
  **[`study_summary()`](https://zeitr.circadia-lab.uk/reference/study_summary.md)**
  — participant-level NPCRA summary across a whole study
- 🚀
  **[`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)**
  — run the complete pipeline on a single file in one call
- 🗃️
  **[`run_pipeline_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_batch.md)**
  — run the complete pipeline across a directory of files
- 🏷️
  **[`label_states()`](https://zeitr.circadia-lab.uk/reference/label_states.md)**
  — convert integer epoch states to a human-readable factor (`"wake"`,
  `"sleep"`, `"nap"`, `"off-wrist"`)
- ⚙️
  **[`acttrust_params()`](https://zeitr.circadia-lab.uk/reference/acttrust_params.md)**
  — device parameter preset; copy and modify to adapt the pipeline to
  other devices

------------------------------------------------------------------------

## 🚀 Getting Started

### Installation

``` r

# Install from GitHub
# install.packages("pak")
pak::pak("circadia-bio/zeitR")
```

### Single recording

``` r

library(zeitR)

rec <- read_actigraphy("recordings/P001.txt", tz = "America/Sao_Paulo")
rec$epochs    # tidy epoch-level tibble
rec$metadata  # device info, firmware, epoch length
```

### Full pipeline

``` r

result <- run_pipeline("recordings/P001.txt", tz = "America/Sao_Paulo")

result$nights  # nightly sleep statistics
result$data    # epoch-level tibble with state, sleep, offwrist columns
result$issues  # timestamp consistency issues (0 rows if none)
```

### State labels

``` r

# Convert integer states to readable labels for display or plotting
result$data$state_label <- label_states(result$data$state)

table(result$data$state_label)
#>      wake     sleep       nap off-wrist
#>     48231     24603       892      2470
```

### Device configuration

``` r

# The pipeline ships with ActTrust-validated defaults via acttrust_params().
# To adapt for a different device, copy and modify:
p <- acttrust_params()
p$sleep$sleep_quantile <- 1/3   # original Crespo (2012) threshold

result <- run_pipeline("recordings/P001.txt", params = p)
```

### Non-parametric circadian rhythm analysis

``` r

npcra <- compute_npcra(rec)
npcra
#>   participant_id    IS    IV    RA    L5 L5_onset   M10 M10_onset n_days n_epochs
#>            P001  0.72  0.43  0.89  12.3    02:30  84.7     11:00    7.0    10080
```

### Whole study

``` r

study <- read_actigraphy_dir("recordings/", tz = "America/Sao_Paulo")
study_summary(study)
#>   participant_id n_epochs n_days    IS    IV    RA    L5 L5_onset   M10 M10_onset
#>            P001    10080   7.00  0.72  0.43  0.89  12.3    02:30  84.7     11:00
#>            P002     9950   6.91  0.68  0.51  0.85  10.1    03:00  79.4     10:30
```

------------------------------------------------------------------------

## 📐 Computed variables

### NPCRA (`compute_npcra()`)

| Variable | Definition |
|----|----|
| `IS` | Interdaily stability — consistency of the 24 h rhythm across days (0–1) |
| `IV` | Intradaily variability — fragmentation of the rest-activity rhythm (≥ 0) |
| `RA` | Relative amplitude — contrast between M10 and L5 (0–1) |
| `L5` / `L5_onset` | Mean activity and onset time of the least active 5 h window |
| `M10` / `M10_onset` | Mean activity and onset time of the most active 10 h window |

### Nightly sleep statistics (`compute_waso()`)

| Variable | Definition                      |
|----------|---------------------------------|
| `tbt`    | Total Bed Time (epochs)         |
| `tst`    | Total Sleep Time (epochs)       |
| `waso`   | Wake After Sleep Onset (epochs) |
| `sol`    | Sleep Onset Latency (epochs)    |
| `soi`    | Sleep Offset Inertia (epochs)   |
| `nw`     | Number of awakenings            |
| `eff`    | Sleep efficiency — TST / TBT    |

------------------------------------------------------------------------

## 🔬 Algorithms

| Step | Algorithm | Reference | Validated |
|----|----|----|----|
| Off-wrist detection | Condor bimodal activity/temperature model | Condor Instruments | ActTrust ✓ |
| Sleep period detection | Crespo adaptive median filter | Crespo et al. (2012) | ActTrust ✓ |
| Nap detection | Crespo zero-proportion filter | Crespo et al. (2012) | ActTrust ✓ |
| Epoch scoring | Cole-Kripke weighted ZCM sum | Cole & Kripke (1992) | ActTrust ✓ |

All four stages have been validated epoch-for-epoch (`0 / 76196`
mismatches) against the Condor circadiaBase Python reference pipeline on
an ActTrust recording. Default parameters in
[`acttrust_params()`](https://zeitr.circadia-lab.uk/reference/acttrust_params.md)
reproduce this reference exactly.

------------------------------------------------------------------------

## 🗂️ Project Structure

    zeitR/
    ├── R/
    │   ├── zeitR-package.R       # package-level docs
    │   ├── read_acttrust.R       # ActTrust file parser
    │   ├── read_actigraphy.R     # device-agnostic wrapper, zeitr_study
    │   ├── prepare.R             # temperature clamping, state columns
    │   ├── consistency.R         # timestamp quality checks
    │   ├── offwrist.R            # detect_offwrist_bimodal()
    │   ├── offwrist_refiner.R    # three-stage BimodalOffwristRefiner port
    │   ├── sleep_periods.R       # detect_sleep_crespo(), detect_naps_crespo()
    │   ├── cole_kripke.R         # score_epochs_cole_kripke()
    │   ├── waso.R                # compute_waso()
    │   ├── npcra.R               # compute_npcra()
    │   ├── study_summary.R       # study_summary()
    │   ├── params.R              # acttrust_params() device preset
    │   ├── pipeline.R            # run_pipeline(), run_pipeline_batch()
    │   └── utils.R               # label_states() + internal helpers
    ├── man/figures/
    │   ├── logo.svg              # hex sticker
    │   └── favicon.svg           # favicon
    ├── vignettes/
    │   ├── getting-started.Rmd
    │   ├── npcra.Rmd
    │   └── study-analysis.Rmd
    ├── tests/testthat/
    ├── .github/workflows/
    │   ├── R-CMD-check.yaml
    │   └── pkgdown.yaml
    ├── DESCRIPTION
    ├── NEWS.md
    └── zeitR.Rproj

------------------------------------------------------------------------

## 📦 Dependencies

| Package   | Version | Purpose                |
|-----------|---------|------------------------|
| cli       | ≥ 3.6.0 | Messages and progress  |
| lubridate | ≥ 1.9.0 | Date/time handling     |
| tibble    | ≥ 3.0.0 | Tidy data frames       |
| tidyr     | ≥ 1.3.0 | Pivoting and reshaping |

------------------------------------------------------------------------

## 👥 Authors

| Role | Name | Affiliation |
|----|----|----|
| Author, maintainer | Lucas França | Northumbria University, Circadia Lab |
| Author | Mario Leocadio-Miguel | Northumbria University, Circadia Lab |

------------------------------------------------------------------------

## 🤝 Related Tools

- 🌙 [**slumbR**](https://github.com/circadia-bio/slumbR) — R companion
  for Sleep Diaries exports (sleep variables, questionnaire scoring)
- 🧮 [**tallieR**](https://github.com/circadia-bio/tallieR) — R
  companion for ScoreMe questionnaire exports
- 🔬 [**circadia-bio**](https://github.com/circadia-bio) — the Circadia
  Lab GitHub organisation

------------------------------------------------------------------------

## 📄 Licence

![](inst/logo.png)

Released under the [MIT License](https://zeitr.circadia-lab.uk/LICENSE).

Copyright © Lucas França, Mario Leocadio-Miguel, 2026
