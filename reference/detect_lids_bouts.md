# Detect nighttime sleep bouts via the Roenneberg relative-immobility method

Standalone sleep-bout detector for raw actigraphy, independent of
zeitR's main Crespo-based pipeline
([`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md))
or the Vallim native pipeline
([`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)).
Intended for LIDS analysis on activity recordings that haven't been run
through either of those. If you already have a `zeitr_result`, pass it
straight to
[`compute_lids()`](https://zeitr.circadia-lab.uk/reference/compute_lids.md)
instead (`bout_source = "state"`, the default there) and skip this
function.

## Usage

``` r
detect_lids_bouts(
  datetime,
  activity,
  relative_threshold = 0.15,
  ma_window_h = 24,
  bridge_min = 5,
  min_bout_min = 30,
  max_gap_min = 15,
  main_window = c("18:00", "08:00"),
  duration_range = c(3, 12),
  one_per_night = TRUE
)
```

## Arguments

- datetime:

  `POSIXct` vector of epoch timestamps, regularly spaced.

- activity:

  Numeric vector of activity counts, same length as `datetime`.

- relative_threshold:

  `numeric(1)`. Default `0.15`.

- ma_window_h:

  `numeric(1)`. Moving-average window in hours. Default `24`.

- bridge_min:

  `numeric(1)`. Maximum active-blip length (minutes) to bridge during
  consolidation. Default `5`.

- min_bout_min:

  `numeric(1)`. Minimum consolidated-run length (minutes) to keep as a
  candidate bout. Default `30`.

- max_gap_min:

  `numeric(1)`. Maximum gap (minutes) between candidate bouts to fuse
  them into one. Default `15`.

- main_window:

  `character(2)`. `c(start, end)` clock times (`"HH:MM"`) defining the
  window a bout must *start* within. Default `c("18:00", "08:00")`
  (wraps midnight).

- duration_range:

  `numeric(2)`. Final bout duration bounds, in hours. Default `c(3, 12)`
  (Winnebeck et al. 2018).

- one_per_night:

  `logical(1)`. If `TRUE` (default), keep only the longest surviving
  bout per calendar night (the "main sleep episode").

## Value

A tibble with one row per detected bout: `bout_id`, `bout_start`,
`bout_end`, `duration_h`, `n_epochs`.

## Details

Ported (with one bug fixed – see below) from a prototype R notebook
translation of Mario Leocadio-Miguel's method, itself adapted from the
relative-immobility algorithm used in Winnebeck et al. (2018) and Hammad
et al. (2026):

1.  Compute a `ma_window_h`-centered moving average of activity (the
    recording's own baseline).

2.  Flag epochs where activity \< `relative_threshold` times that moving
    average as *candidate sleep*.

3.  **Consolidate**: brief active blips (\<= `bridge_min` minutes),
    surrounded on both sides by candidate-sleep epochs, are relabelled
    as sleep.

4.  Keep only consolidated runs lasting \>= `min_bout_min`.

5.  **Fuse** bouts separated by a gap \<= `max_gap_min`. (This is the
    step truncated mid-statement in the source notebook, `fused.a...` –
    reimplemented in full here as the internal `.fuse_bouts()`. Unlike
    the notebook, which pads fused gaps with `NaN` to preserve a regular
    time axis, the fused bout here simply spans `start` to `end` over
    the original activity values – fine for bout *timing*, but if you
    need the gap epochs excluded from the LIDS fit itself, mask them to
    `NA` in `activity` beforehand.)

6.  Restrict to bouts starting within `main_window` (default
    18:00-08:00), excluding daytime naps / recording artefacts.

7.  Filter by total duration (`duration_range`) and, if
    `one_per_night = TRUE`, keep only the longest bout per calendar
    night.

The moving-average step uses zeitR's border-replicated
`rolling_mean_cpp()` rather than pandas' `min_periods=1` edge behaviour
(which shrinks the window near the recording boundary instead of
replicating the edge value) – a minor difference confined to the first/
last `ma_window_h`/2 hours of the recording.

## References

Winnebeck, E. C., Fischer, D., Leise, T., & Roenneberg, T. (2018).
Dynamics and Ultradian Structure of Human Sleep in Real Life. *Current
Biology*, 28(1), 49-59.e5.
[doi:10.1016/j.cub.2017.11.063](https://doi.org/10.1016/j.cub.2017.11.063)

## See also

[`compute_lids()`](https://zeitr.circadia-lab.uk/reference/compute_lids.md),
[`lids_transform()`](https://zeitr.circadia-lab.uk/reference/lids_transform.md),
[`fit_lids()`](https://zeitr.circadia-lab.uk/reference/fit_lids.md)
