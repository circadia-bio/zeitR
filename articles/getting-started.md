# Getting started with zeitR

## What is zeitR?

**zeitR** is an R package for importing, parsing, and analysing raw
actigraphy recordings from wrist-worn devices. It reads device output
files, runs a full rest-activity pipeline, and returns tidy data frames
ready for chronobiological analysis.

zeitR is part of the [Circadia Lab](https://circadia-lab.uk) R
ecosystem:

| Package | Purpose |
|----|----|
| **zeitR** | Actigraphy parsing and rest-activity analysis |
| [slumbR](https://github.com/circadia-bio/slumbR) | Sleep diary and questionnaire scoring |
| [tallieR](https://github.com/circadia-bio/tallieR) | ScoreMe questionnaire exports |

------------------------------------------------------------------------

## Installation

``` r

# Install from GitHub
# install.packages("pak")
pak::pak("circadia-bio/zeitR")
```

------------------------------------------------------------------------

## Reading a single recording

The main entry point is
[`read_actigraphy()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy.md).
It accepts a path to a raw device file and returns a `zeitr_recording`
object with two slots: `$epochs` (the epoch-level tibble) and
`$metadata` (device and subject information from the file header).

``` r

rec <- read_actigraphy("recordings/P001.txt", tz = "America/Sao_Paulo")
rec
```

    ## zeitr_recording: P001
    ## * Device:   ActTrust2 (ID 828)
    ## * Firmware: 1.4
    ## * Interval: 60 s
    ## * Epochs:   10080
    ## * From:     2021-05-27 11:10:15
    ## * To:       2021-06-03 11:09:15

### The `$epochs` tibble

``` r

rec$epochs
```

    ## # A tibble: 10,080 x 8
    ##    datetime            activity int_temp ext_temp  ZCMn state offwrist sleep
    ##    <dttm>                 <dbl>    <dbl>    <dbl> <dbl> <dbl>    <dbl> <dbl>
    ##  1 2021-05-27 11:10:15     4856     24.2     23.9  0        0        0     0
    ##  2 2021-05-27 11:11:15     4483     24.4     24.1  1.5      0        0     0
    ##  3 2021-05-27 11:12:15      425     24.4     23.9  0.05     0        0     0
    ## ...

Each row is one epoch (60 seconds by default). The key columns are:

| Column     | Description                                                  |
|------------|--------------------------------------------------------------|
| `datetime` | Epoch timestamp                                              |
| `activity` | PIM activity count                                           |
| `int_temp` | Internal (on-body) temperature, °C                           |
| `ext_temp` | External (ambient) temperature, °C                           |
| `ZCMn`     | Normalised zero-crossing mode count                          |
| `state`    | Pipeline state (0 = wake, 1 = sleep, 4 = off-wrist, 7 = nap) |
| `offwrist` | Off-wrist indicator (0.25 = off-wrist)                       |
| `sleep`    | Binary sleep/wake (1 = sleep)                                |

### The `$metadata` list

``` r

rec$metadata
```

    ## $subject
    ## [1] "Julius"
    ## $device_model
    ## [1] "ActTrust2"
    ## $interval_s
    ## [1] 60
    ## ...

------------------------------------------------------------------------

## Running the full pipeline

[`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
chains all analysis steps in one call: timestamp consistency check →
prepare → off-wrist detection → sleep period detection (Crespo) → nap
detection (Crespo) → WASO computation (Cole-Kripke).

``` r

result <- run_pipeline("recordings/P001.txt", tz = "America/Sao_Paulo")
result
```

    ## zeitr_result: P001
    ## * Source:  .../recordings/P001.txt
    ## * Epochs:  10080
    ## * Nights:  7
    ## * Issues:  0

### Nightly sleep statistics

``` r

result$nights
```

    ## # A tibble: 7 x 11
    ##   night is_nap bed_time            get_up_time          tbt   tst  waso   sol   soi    nw   eff
    ##   <int> <lgl>  <dttm>              <dttm>             <int> <dbl> <dbl> <int> <int> <int> <dbl>
    ## 1     1 FALSE  2021-05-27 22:04:15 2021-05-28 07:18:15  494   432    48    12     2    18  0.874
    ## ...

The metrics are:

| Column | Definition                      |
|--------|---------------------------------|
| `tbt`  | Total Bed Time (epochs)         |
| `tst`  | Total Sleep Time (epochs)       |
| `waso` | Wake After Sleep Onset (epochs) |
| `sol`  | Sleep Onset Latency (epochs)    |
| `soi`  | Sleep Offset Inertia (epochs)   |
| `nw`   | Number of awakenings            |
| `eff`  | Sleep efficiency (TST / TBT)    |

> **Tip:** Multiply epoch counts by
> `result$data$metadata$interval_s / 60` to convert to minutes.

### Timestamp issues

``` r

result$issues
```

If no issues were found, this is a zero-row tibble. Issues are flagged
for gaps \> 120 s, backward jumps, and year artefacts (1970 or 2000).

------------------------------------------------------------------------

## Running pipeline steps individually

For more control, each step can be called separately. This is useful
when you want to inspect intermediate outputs or customise algorithm
parameters.

``` r

# 1. Read
rec  <- read_actigraphy("recordings/P001.txt", tz = "America/Sao_Paulo")

# 2. Check timestamps
issues <- check_consistency(rec$epochs)

# 3. Prepare (temperature clamping, state columns)
prep <- prepare_actigraphy(rec$epochs)

# 4. Off-wrist detection (Condor bimodal algorithm)
prep <- detect_offwrist_bimodal(prep)
sum(prep$state == 4)  # off-wrist epochs

# 5. Main sleep periods (Crespo, 2012)
prep <- detect_sleep_crespo(prep)

# 6. Naps (Crespo, 2012)
prep <- detect_naps_crespo(prep)

# 7. WASO + nightly stats (Cole-Kripke, 1992)
out <- compute_waso(prep)
out$nights
```

------------------------------------------------------------------------

## Next steps

- [`vignette("npcra")`](https://zeitr.circadia-lab.uk/articles/npcra.md)
  — computing IS, IV, RA, L5, and M10 from a recording
- [`vignette("study-analysis")`](https://zeitr.circadia-lab.uk/articles/study-analysis.md)
  — batch reading and participant-level summaries
- [`?run_pipeline`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
  — full documentation of all pipeline parameters
- [`?detect_sleep_crespo`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md)
  — Crespo algorithm parameters
- [`?score_epochs_cole_kripke`](https://zeitr.circadia-lab.uk/reference/score_epochs_cole_kripke.md)
  — Cole-Kripke scoring details
