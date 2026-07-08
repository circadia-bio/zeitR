# The Vallim pipeline: native sleep classification

This vignette introduces
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md),
the **Vallim pipeline** — a sleep classification layer developed by
Julia Ribeiro da Silva Vallim that addresses limitations of the base
CSPD pipeline in detecting and classifying sleep episodes from
actigraphy data.

------------------------------------------------------------------------

## Background

The Condor CSPD pipeline
([`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md))
detects sleep periods using an adaptive median filter and returns one
row per period. It works well for clear-cut recordings but can
misclassify episodes in several edge cases:

- Daytime naps adjacent to a long main sleep (date collision)
- Recordings ending mid-sleep (truncated episodes)
- Sleep periods fragmented across multiple short blocks
- Unusually long in-bed periods (TBT \> 14 h) that should be split

The Vallim pipeline keeps the same off-wrist detection and CSPD epoch
scorer as
[`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
but replaces the final classification step with a 7-rule adaptive system
ported from Julia Vallim’s Python notebook
(`vs_condor_py_pipeline_fix29_jrsv.ipynb`).

------------------------------------------------------------------------

## The rule set

The classification follows this fixed execution order (Fix 27):

| Step | Fix/Rule | Description |
|----|----|----|
| 1 | Fix 25 | Exclude truncated episodes at the recording end |
| 2 | Fix 26a | Infer adaptive nocturnal window from wrist temperature and light |
| 3 | Fix 29 / Rule 2 | Split TBT 14–16 h episodes; exclude TBT \> 16 h |
| 4 | Fix 26c | Recover fragmented sleep nights missed by the scorer |
| 5 | Rules 3–5 | Classify as main or secondary using the nocturnal window |
| 6 | Fix 26b | Resolve sleep-date collisions from the noon threshold |
| 7 | Rule 6 | Keep the longest main episode per sleep date |
| 8 | Rule 7 | Exclude all episodes on days with no main sleep |

The nocturnal window (Fix 26a) is inferred adaptively per participant
rather than fixed at 18:00–06:00. Episodes where wrist temperature ≥ 28
°C and ambient light ≤ 5 lux are used as candidates; the circular mean
of their onsets becomes the window anchor (±4 h).

------------------------------------------------------------------------

## Basic usage

``` r

FILE <- system.file("extdata", "input1.txt", package = "zeitR")
TZ   <- "America/Sao_Paulo"

result <- run_pipeline_native(FILE, tz = TZ, quiet = TRUE)
#> ℹ Reading input1.txt ...
#> ✔ [input1] Done. 52 main night(s), 0 secondary episode(s).
result
#> 
#> ── zeitr_result: input1 ────────────────────────────────────────────────────────
#> • Source: /home/runner/work/_temp/Library/zeitR/extdata/input1.txt
#> • Epochs: 76196
#> • Nights: 52
#> • Issues: 5
```

The `nights` table has the same schema as
[`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
but adds a `sleep_type` column:

``` r

result$nights |>
  mutate(
    bed_time    = format(bed_time,    "%Y-%m-%d %H:%M"),
    get_up_time = format(get_up_time, "%Y-%m-%d %H:%M"),
    tbt_h       = round(tbt / 60, 1),
    tst_h       = round(tst / 60, 1),
    eff_pct     = round(eff * 100, 1)
  ) |>
  select(night, sleep_type, bed_time, get_up_time, tbt_h, tst_h, waso, sol, eff_pct) |>
  head(10) |>
  knitr::kable(
    col.names = c("Night", "Type", "Bed time", "Get-up time",
                  "TBT (h)", "TST (h)", "WASO (min)", "SOL (min)", "Eff (%)"),
    align = "clllrrrrr"
  )
```

| Night | Type | Bed time | Get-up time | TBT (h) | TST (h) | WASO (min) | SOL (min) | Eff (%) |
|:--:|:---|:---|:---|---:|---:|---:|---:|---:|
| 1 | main | 2021-05-27 23:42 | 2021-05-28 06:53 | 7.2 | 6.9 | 18 | 0 | 95.8 |
| 2 | main | 2021-05-28 23:30 | 2021-05-29 07:15 | 7.8 | 7.3 | 25 | 0 | 94.4 |
| 3 | main | 2021-05-30 22:18 | 2021-05-31 07:57 | 9.7 | 9.2 | 27 | 0 | 95.3 |
| 4 | main | 2021-05-31 23:24 | 2021-06-01 07:13 | 7.8 | 7.6 | 12 | 0 | 97.0 |
| 5 | main | 2021-06-02 00:10 | 2021-06-02 07:18 | 7.2 | 6.8 | 23 | 0 | 94.6 |
| 6 | main | 2021-06-02 23:05 | 2021-06-03 07:46 | 8.7 | 8.2 | 30 | 0 | 94.3 |
| 7 | main | 2021-06-04 00:01 | 2021-06-04 08:28 | 8.5 | 7.8 | 39 | 0 | 92.3 |
| 8 | main | 2021-06-04 23:36 | 2021-06-05 05:19 | 5.7 | 5.3 | 4 | 9 | 93.3 |
| 9 | main | 2021-06-06 00:34 | 2021-06-06 06:19 | 5.8 | 5.7 | 7 | 0 | 98.0 |
| 10 | main | 2021-06-07 00:13 | 2021-06-07 06:58 | 6.8 | 6.7 | 6 | 0 | 98.5 |

------------------------------------------------------------------------

## Comparing the two pipelines

Running both pipelines on the same recording lets you see exactly what
the Vallim rule set changes:

``` r

cspd   <- run_pipeline(FILE, tz = TZ, quiet = TRUE)
#> ℹ Reading input1.txt ...
#> ✔ [input1] Done. 55 night(s) detected.
native <- run_pipeline_native(FILE, tz = TZ, quiet = TRUE)
#> ℹ Reading input1.txt ...
#> ✔ [input1] Done. 52 main night(s), 0 secondary episode(s).

cat("CSPD nights:  ", nrow(cspd$nights),   "\n")
#> CSPD nights:   55
cat("Vallim nights:", nrow(native$nights), "\n")
#> Vallim nights: 52
cat("  main:      ", sum(native$nights$sleep_type == "main"),      "\n")
#>   main:       52
cat("  secondary: ", sum(native$nights$sleep_type == "secondary"), "\n")
#>   secondary:  0
```

The Vallim pipeline consolidates duplicate-date episodes (Rule 6) and
discards dates where no valid main sleep exists (Rule 7), resulting in a
cleaner set of main nights.

------------------------------------------------------------------------

## Sleep onset and offset statistics

[`circ_mean_h()`](https://zeitr.circadia-lab.uk/reference/circ_mean_h.md)
and
[`circ_sd_h()`](https://zeitr.circadia-lab.uk/reference/circ_sd_h.md)
compute circular statistics for clock-time variables. Using standard
arithmetic mean for times like sleep onset (which often straddles
midnight) inflates variability estimates and produces biased means; the
circular mean handles the wrap correctly.

``` r

main <- native$nights |> filter(sleep_type == "main")

# Decimal hours of bed time and get-up time
onset_h  <- as.numeric(format(main$bed_time,    "%H", tz = TZ)) +
            as.numeric(format(main$bed_time,    "%M", tz = TZ)) / 60
offset_h <- as.numeric(format(main$get_up_time, "%H", tz = TZ)) +
            as.numeric(format(main$get_up_time, "%M", tz = TZ)) / 60

cat(sprintf("Mean sleep onset:  %05.2f h  (SD = %.2f h)\n",
            circ_mean_h(onset_h),  circ_sd_h(onset_h)))
#> Mean sleep onset:  23.80 h  (SD = 0.89 h)
cat(sprintf("Mean sleep offset: %05.2f h  (SD = %.2f h)\n",
            circ_mean_h(offset_h), circ_sd_h(offset_h)))
#> Mean sleep offset: 07.73 h  (SD = 0.77 h)
```

The
[`circ_mean_h()`](https://zeitr.circadia-lab.uk/reference/circ_mean_h.md)
function returns decimal hours in \[0, 24). To display as HH:MM:

``` r

fmt_h <- function(h) sprintf("%02d:%02d", as.integer(h) %% 24L,
                              as.integer((h %% 1) * 60))

cat("Mean sleep onset: ", fmt_h(circ_mean_h(onset_h)),  "\n")
#> Mean sleep onset:  23:47
cat("Mean get-up time: ", fmt_h(circ_mean_h(offset_h)), "\n")
#> Mean get-up time:  07:43
```

------------------------------------------------------------------------

## Advanced: running each step directly

For custom workflows,
[`extract_sleep_episodes()`](https://zeitr.circadia-lab.uk/reference/extract_sleep_episodes.md)
and
[`classify_sleep_episodes()`](https://zeitr.circadia-lab.uk/reference/classify_sleep_episodes.md)
expose the two stages independently.

``` r

# Step 1: read, prepare, detect off-wrist and CSPD-score epochs
rec  <- read_acttrust(FILE, tz = TZ)
prep <- prepare_actigraphy(rec)
prep <- detect_offwrist_bimodal(prep)
prep <- detect_sleep_crespo(prep)

# Step 2: extract per-episode statistics from the CSPD state
episodes_raw <- extract_sleep_episodes(prep)

# Step 3: apply the Vallim rule set
episodes <- classify_sleep_episodes(
  episodes = episodes_raw,
  data     = prep,
  verbose  = TRUE           # print each fix as it fires
)

episodes |> filter(sleep_type == "main")
```

The `verbose = TRUE` argument prints a trace of each rule:

    [Fix 26a] anchor=23.8 h -> nocturnal window 19.8 h-3.8 h
    [Fix 26c] No fragmented episodes recovered.
    [Fix 26b] No sleep-date collisions detected.

------------------------------------------------------------------------

## Batch processing

[`run_pipeline_native_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native_batch.md)
works identically to
[`run_pipeline_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_batch.md):

``` r

results <- run_pipeline_native_batch(
  "recordings/",
  tz    = "America/Sao_Paulo",
  quiet = TRUE
)

# Combine all main nights
main_all <- bind_rows(
  lapply(names(results), function(id) {
    results[[id]]$nights |>
      filter(sleep_type == "main") |>
      mutate(subject_id = id)
  })
)
```

------------------------------------------------------------------------

## Next steps

- [`vignette("sleep-analysis")`](https://zeitr.circadia-lab.uk/articles/sleep-analysis.md)
  — CSPD pipeline walkthrough with actogram
- [`?export_hypnogram`](https://zeitr.circadia-lab.uk/reference/export_hypnogram.md)
  — export results to `hypnoR` for sleep architecture metrics
- [`?compute_sleep_metrics`](https://zeitr.circadia-lab.uk/reference/compute_sleep_metrics.md),
  [`?compute_cpd_metrics`](https://zeitr.circadia-lab.uk/reference/compute_cpd_metrics.md)
  — day-type sleep summary and chronotype
- [`?classify_sleep_episodes`](https://zeitr.circadia-lab.uk/reference/classify_sleep_episodes.md)
  — full parameter reference for the rule set
- [`?circ_mean_h`](https://zeitr.circadia-lab.uk/reference/circ_mean_h.md),
  [`?circ_sd_h`](https://zeitr.circadia-lab.uk/reference/circ_sd_h.md) —
  circular statistics for clock-time variables
- [`?acttrust_params`](https://zeitr.circadia-lab.uk/reference/acttrust_params.md)
  — device parameter defaults
