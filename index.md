# ⌚️ zeitR

**Actigraphy data parsing and analysis for R.**

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://zeitr.circadia-lab.uk/LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3)](https://www.r-project.org/)
[![R CMD
CHECK](https://github.com/circadia-bio/zeitR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/circadia-bio/zeitR/actions/workflows/R-CMD-check.yaml)
[![Coverage](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/circadia-bio/zeitR/gh-pages/badges/coverage.json)](https://github.com/circadia-bio/zeitR/actions/workflows/pkgdown.yaml)
[![Status](https://img.shields.io/badge/status-early%20development-orange)](https://github.com/circadia-bio/zeitR)
[![pkgdown](https://img.shields.io/badge/docs-zeitr.circadia--lab.uk-F0A500)](https://zeitr.circadia-lab.uk)

------------------------------------------------------------------------

> \[!WARNING\] **zeitR is in early development and has not been formally
> validated.** The CSPD pipeline has been validated epoch-for-epoch
> against the Condor circadiaBase Python reference on an ActTrust
> recording. The Vallim pipeline has been validated at the
> classification level against Julia Vallim’s Python reference notebook.
> Neither pipeline has undergone formal peer review. Verify outputs
> independently before using in any research context.

------------------------------------------------------------------------

## 📖 What is zeitR?

zeitR is an R package for importing, parsing, and analysing raw
actigraphy recordings from wrist-worn devices. It runs a full
rest-activity pipeline — off-wrist detection, sleep period
identification, WASO computation — and computes standard non-parametric
circadian rhythm variables (IS, IV, RA, L5, M10), returning tidy data
frames ready for downstream chronobiological analysis.

zeitR ships two end-to-end pipelines:

- **[`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)**
  — the Condor CSPD pipeline, validated epoch-for-epoch against the
  Condor circadiaBase Python reference.
- **[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)**
  — the **Vallim pipeline**, a post-processing layer developed by Julia
  Ribeiro da Silva Vallim that replaces Condor’s classification logic
  with an adaptive 7-rule rule set (Fix 25, 26a/b/c, 27, 29, Rules 3–7).
  Handles edge cases including fragmented episodes, long sleep periods,
  date collisions, and adaptive nocturnal window inference from wrist
  temperature and ambient light.

zeitR is designed to complement
[slumbR](https://github.com/circadia-bio/slumbR) in the Circadia Lab
ecosystem: slumbR handles sleep diary and questionnaire data, zeitR
handles the actigraphy side of a study.

------------------------------------------------------------------------

## ✨ Features

- 📥
  **[`read_actigraphy()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy.md)**
  — parse a raw device file into a `zeitr_recording` object
- 📂
  **[`read_actigraphy_dir()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy_dir.md)**
  — batch-read a whole directory into a `zeitr_study`
- 🔍
  **[`check_consistency()`](https://zeitr.circadia-lab.uk/reference/check_consistency.md)**
  — flag timestamp gaps, backward jumps, and firmware artefacts
- 🦾
  **[`detect_offwrist_bimodal()`](https://zeitr.circadia-lab.uk/reference/detect_offwrist_bimodal.md)**
  — Condor bimodal activity/temperature off-wrist detection
- 😴
  **[`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md)**
  — main sleep period detection (Crespo et al., 2012)
- 💤
  **[`detect_naps_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_naps_crespo.md)**
  — secondary sleep period detection (Crespo et al., 2012)
- ⏱️
  **[`score_epochs_cole_kripke()`](https://zeitr.circadia-lab.uk/reference/score_epochs_cole_kripke.md)**
  — epoch-level wake/sleep scoring (Cole & Kripke, 1992)
- 📊
  **[`compute_waso()`](https://zeitr.circadia-lab.uk/reference/compute_waso.md)**
  — nightly TBT, TST, WASO, SOL, SOI, awakenings, sleep efficiency
- 📐
  **[`compute_npcra()`](https://zeitr.circadia-lab.uk/reference/compute_npcra.md)**
  — non-parametric circadian rhythm analysis (IS, IV, RA, L5, M10)
- 🗂️
  **[`study_summary()`](https://zeitr.circadia-lab.uk/reference/study_summary.md)**
  — participant-level NPCRA summary across a whole study
- 🗂️
  **[`study_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/study_sleep_metrics.md)**
  — participant-level sleep-timing/chronotype summary (CPD, MSF/MSW,
  SJL) across a whole study
- 📋
  **[`compute_sleep_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md)**
  — per-night sleep metrics split by day type (overall / workday / free
  day)
- 📋
  **[`compute_cpd_metrics()`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)**
  — CPD, MSW, MSF, MSFsc, social jet lag (SJL, SJLa)
- 🚀
  **[`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)**
  — full CSPD pipeline on a single file
- 🗃️
  **[`run_pipeline_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_batch.md)**
  — CSPD pipeline across a directory
- 🌙
  **[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)**
  — full Vallim pipeline on a single file
- 🌙
  **[`run_pipeline_native_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native_batch.md)**
  — Vallim pipeline across a directory
- 📤
  **[`export_hypnogram()`](https://zeitr.circadia-lab.uk/reference/export_hypnogram.md)**
  — export to `hypnoR` format (`W` / `Sleep` / `Quiet sleep`;
  `subject_id` auto-inferred)
- 🔬
  **[`extract_sleep_episodes()`](https://zeitr.circadia-lab.uk/reference/extract_sleep_episodes.md)**
  — extract per-episode statistics from a CSPD-scored table
- 🔬
  **[`classify_sleep_episodes()`](https://zeitr.circadia-lab.uk/reference/classify_sleep_episodes.md)**
  — apply the JRSV rule set to classify episodes
- 📏
  **[`circ_mean_h()`](https://zeitr.circadia-lab.uk/reference/circ_mean_h.md)**
  — circular mean for clock-time variables (handles midnight wrap)
- 📏
  **[`circ_sd_h()`](https://zeitr.circadia-lab.uk/reference/circ_sd_h.md)**
  — circular SD for clock-time variables
- 🏷️
  **[`label_states()`](https://zeitr.circadia-lab.uk/reference/label_states.md)**
  — convert integer epoch states to a human-readable factor
- ⚙️
  **[`acttrust_params()`](https://zeitr.circadia-lab.uk/reference/acttrust_params.md)**
  — device parameter preset for the ActTrust actigraph

------------------------------------------------------------------------

## 🚀 Getting Started

### Installation

``` r

install.packages("pak")
pak::pak("circadia-bio/zeitR")
```

### CSPD pipeline

``` r

library(zeitR)

result <- run_pipeline("recordings/P001.txt", tz = "America/Sao_Paulo")

result$nights  # nightly sleep statistics
result$data    # epoch-level tibble with state column
result$issues  # timestamp consistency flags
```

### Vallim pipeline

``` r

result <- run_pipeline_native("recordings/P001.txt", tz = "America/Sao_Paulo")

# nights has an additional sleep_type column ("main" / "secondary")
result$nights |> dplyr::filter(sleep_type == "main")
```

### State labels

``` r

result$data$state_label <- label_states(result$data$state)

table(result$data$state_label)
#>      wake     sleep       nap off-wrist
#>     48231     24603       892      2470
```

### Circular statistics

``` r

main_nights <- result$nights |> dplyr::filter(sleep_type == "main")
onset_h     <- as.numeric(format(main_nights$bed_time, "%H")) +
               as.numeric(format(main_nights$bed_time, "%M")) / 60

circ_mean_h(onset_h)  # mean sleep onset (handles midnight wrap)
circ_sd_h(onset_h)    # within-person variability in sleep onset
```

### Non-parametric circadian rhythm analysis

``` r

rec   <- read_acttrust("recordings/P001.txt", tz = "America/Sao_Paulo")
npcra <- compute_npcra(rec)
npcra
#>   IS    IV    RA    L5 L5_onset   M10 M10_onset n_days
#>   0.72  0.43  0.89  12.3    02:30  84.7     11:00    7.0
```

### Device configuration

``` r

p <- acttrust_params()
p$sleep$sleep_quantile <- 1/3   # original Crespo (2012) threshold

result <- run_pipeline("recordings/P001.txt", params = p)
```

------------------------------------------------------------------------

## 📐 Computed variables

### NPCRA (`compute_npcra()`)

| Variable | Definition |
|----|----|
| `IS` | Interdaily stability — consistency of the 24 h rhythm across days (0–1) |
| `IV` | Intradaily variability — fragmentation of the rest-activity rhythm (≥ 0) |
| `RA` | Relative amplitude — contrast between M10 and L5 (0–1) |
| `L5` / `L5_onset` | Mean activity and onset of the least active 5 h window |
| `M10` / `M10_onset` | Mean activity and onset of the most active 10 h window |

### Nightly sleep statistics

| Variable     | Definition                                       |
|--------------|--------------------------------------------------|
| `tbt`        | Total Bed Time (minutes)                         |
| `tst`        | Total Sleep Time (minutes)                       |
| `waso`       | Wake After Sleep Onset (minutes)                 |
| `sol`        | Sleep Onset Latency (minutes)                    |
| `soi`        | Sleep Offset Inertia (minutes)                   |
| `nw`         | Number of awakenings                             |
| `eff`        | Sleep efficiency — TST / TBT                     |
| `sleep_type` | `"main"` or `"secondary"` (Vallim pipeline only) |

------------------------------------------------------------------------

## 🔬 Algorithms

| Step | Algorithm | Reference | Validated |
|----|----|----|----|
| Off-wrist detection | Condor bimodal activity/temperature model | Condor Instruments | ActTrust ✓ |
| Sleep period detection | Crespo adaptive median filter | Crespo et al. (2012) | ActTrust ✓ |
| Nap detection | Crespo zero-proportion filter | Crespo et al. (2012) | ActTrust ✓ |
| Epoch scoring | Cole-Kripke weighted ZCM sum | Cole & Kripke (1992) | ActTrust ✓ |
| Episode classification | Vallim JRSV rule set (Fixes 25, 26a/b/c, 27, 29) | Vallim (2024) | ActTrust ✓ |
| Sleep summary | Day-type metric split (overall / workday / free day) | Vallim (2024) | ActTrust ✓ |
| Chronotype | CPD, MSW, MSF, MSFsc, SJL | Roenneberg et al. | ActTrust ✓ |

The CSPD pipeline has been validated epoch-for-epoch (0 / 76,196
mismatches) against the Condor circadiaBase Python reference. The Vallim
pipeline has been validated at the classification level: all 52 main
nights on the ActTrust validation recording classified identically to
Julia Vallim’s Python reference notebook. R is now the reference
implementation for Fix 26c (fragment recovery), which correctly uses
Cole-Kripke epoch scoring and proper temperature/light column names that
were mismatched in the Python original.

------------------------------------------------------------------------

## 🗂️ Project Structure

    zeitR/
    ├── R/
    │   ├── zeitR-package.R       # package-level docs and Rcpp registration
    │   ├── read_acttrust.R       # ActTrust file parser
    │   ├── read_actigraphy.R     # device-agnostic wrapper, zeitr_study
    │   ├── prepare.R             # temperature clamping, state column init
    │   ├── consistency.R         # timestamp quality checks
    │   ├── offwrist.R            # detect_offwrist_bimodal()
    │   ├── offwrist_refiner.R    # three-stage BimodalOffwristRefiner port
    │   ├── sleep_periods.R       # detect_sleep_crespo(), detect_naps_crespo()
    │   ├── sleep_classify.R      # Vallim pipeline: extract + classify episodes
    │   ├── sleep_metrics.R       # compute_sleep_metrics(), compute_cpd_metrics()
    │   ├── cole_kripke.R         # score_epochs_cole_kripke()
    │   ├── waso.R                # compute_waso()
    │   ├── npcra.R               # compute_npcra()
    │   ├── study_summary.R       # study_summary()
    │   ├── study_sleep_metrics.R # study_sleep_metrics()
    │   ├── circ_utils.R          # circ_mean_h(), circ_sd_h()
    │   ├── params.R              # acttrust_params()
    │   ├── pipeline.R            # run_pipeline*(), run_pipeline_native*()
    │   ├── export.R              # export_hypnogram()
    │   └── utils.R               # label_states() + Rcpp wrappers + helpers
    ├── src/
    │   └── rolling_filters.cpp   # Rcpp: rolling filters, diff5, Cole-Kripke
    ├── man/figures/
    │   ├── logo.svg
    │   └── favicon.svg
    ├── vignettes/
    │   ├── getting-started.Rmd
    │   ├── npcra.Rmd
    │   ├── sleep-analysis.Rmd    # CSPD pipeline walkthrough
    │   ├── study-analysis.Rmd
    │   └── vallim-pipeline.Rmd   # Vallim pipeline walkthrough
    ├── tests/testthat/
    ├── inst/extdata/             # validation fixtures
    ├── DESCRIPTION
    ├── NEWS.md
    └── zeitR.Rproj

------------------------------------------------------------------------

## 📦 Dependencies

| Package   | Version | Purpose                               |
|-----------|---------|---------------------------------------|
| cli       | ≥ 3.6.0 | Messages and progress                 |
| lubridate | ≥ 1.9.0 | Date/time handling                    |
| mclust    | any     | Bimodal GMM for off-wrist detection   |
| Rcpp      | ≥ 1.0.0 | C++ rolling filters and epoch scoring |
| tibble    | ≥ 3.0.0 | Tidy data frames                      |
| tidyr     | ≥ 1.3.0 | Pivoting and reshaping                |

------------------------------------------------------------------------

## 👥 Authors

| Role | Name | Affiliation |
|----|----|----|
| Author, maintainer | Lucas França | Northumbria University, Circadia Lab |
| Author | Mario Leocadio-Miguel | Northumbria University, Circadia Lab |
| Author | Julia Ribeiro da Silva Vallim | Universidade Federal de São Paulo |

------------------------------------------------------------------------

## 🤝 Related Tools

- 🌙 [**slumbR**](https://github.com/circadia-bio/slumbR) — sleep diary
  processing and circadian metrics
- 🧮 [**tallieR**](https://github.com/circadia-bio/tallieR) —
  sociodemographic and questionnaire scoring
- 🔄 [**syncR**](https://github.com/circadia-bio/syncR) — unified
  participant-indexed database for the Circadia ecosystem
- 🔬 [**circadia-bio**](https://github.com/circadia-bio) — the Circadia
  Lab GitHub organisation

------------------------------------------------------------------------

## 📄 Licence

Released under the [MIT License](https://zeitr.circadia-lab.uk/LICENSE).

Copyright © Lucas França, Mario Leocadio-Miguel, 2026
