# Extract sleep episodes from a CSPD-scored epoch table

Converts the epoch-level `state` vector produced by
[`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md)
into a per-episode tibble analogous to the Python SleepPipeline's
`results.nights`. Contiguous on-wrist sleep periods are delimited by
`.nights_df()`; per-epoch wake/sleep within each period is scored by
[`score_epochs_cole_kripke()`](https://zeitr.circadia-lab.uk/reference/score_epochs_cole_kripke.md)
to derive WASO, SOL, SOI, TST, EFF, and NW. All time metrics are
returned in **minutes**.

## Usage

``` r
extract_sleep_episodes(data, wake_thresh = 60L)
```

## Arguments

- data:

  A tibble as returned by
  [`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md),
  containing columns `datetime`, `ZCMn`, and `state`.

- wake_thresh:

  `integer(1)`. Minimum wake-bout length (epochs) used by `.nights_df()`
  to separate episode boundaries. Default `60L`.

## Value

A tibble with one row per detected episode and columns `bts` (POSIXct),
`gts` (POSIXct), `tbt`, `tst`, `sol`, `soi`, `waso` (minutes), `nw`
(integer), `eff` (0-1), `nap` (logical, always `FALSE` – classification
into main/secondary happens in
[`classify_sleep_episodes()`](https://zeitr.circadia-lab.uk/reference/classify_sleep_episodes.md)).

## See also

[`classify_sleep_episodes()`](https://zeitr.circadia-lab.uk/reference/classify_sleep_episodes.md),
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)
