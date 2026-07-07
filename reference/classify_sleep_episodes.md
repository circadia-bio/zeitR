# Classify sleep episodes as main or secondary (native pipeline)

Applies the JRSV classification rule set to a raw episode table from
[`extract_sleep_episodes()`](https://zeitr.circadia-lab.uk/reference/extract_sleep_episodes.md),
assigning each episode a `sleep_type` of `"main"` or `"secondary"`. The
execution order follows Fix 27:

## Usage

``` r
classify_sleep_episodes(
  episodes,
  data,
  max_tib_h = 16,
  max_main_tib_h = 14,
  min_main_tib_h = 4,
  nocturnal_onset_start = 18,
  nocturnal_onset_end = 6,
  temp_thresh = 28,
  light_thresh_window = 5,
  light_thresh_recovery = 10,
  min_fragment_tib_h = 1,
  rolling_window_min = 15L,
  max_split_iterations = 3L,
  collision_gap_h = 4,
  verbose = FALSE
)
```

## Arguments

- episodes:

  A tibble as returned by
  [`extract_sleep_episodes()`](https://zeitr.circadia-lab.uk/reference/extract_sleep_episodes.md).

- data:

  A tibble as returned by
  [`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md),
  containing at minimum `datetime`, `state`, `activity`. The ActTrust
  channels `int_temp` (wrist temperature) and `light` (ambient lux) must
  be present for Fix 26a/c; the function falls back to defaults when
  absent.

- max_tib_h:

  `numeric(1)`. Episodes with TBT \> this value are excluded directly.
  Default `16`.

- max_main_tib_h:

  `numeric(1)`. Episodes with TBT in `(max_main_tib_h, max_tib_h]` are
  split first. Default `14`.

- min_main_tib_h:

  `numeric(1)`. Minimum TBT (hours) for a main episode. Default `4`.

- nocturnal_onset_start:

  `numeric(1)`. Default nocturnal window start (decimal hours).
  Overridden by Fix 26a. Default `18`.

- nocturnal_onset_end:

  `numeric(1)`. Default nocturnal window end. Overridden by Fix 26a.
  Default `6`.

- temp_thresh:

  `numeric(1)`. Minimum wrist temperature (degC) for candidate episodes
  in Fix 26a/c. Default `28`.

- light_thresh_window:

  `numeric(1)`. Maximum ambient light (lux) for the nocturnal window
  inference (Fix 26a). Default `5`.

- light_thresh_recovery:

  `numeric(1)`. Maximum ambient light (lux) for gap merging in fragment
  recovery (Fix 26c). Default `10`.

- min_fragment_tib_h:

  `numeric(1)`. Minimum TBT (hours) for each fragment produced by an
  episode split (Fix 29). Default `1`.

- rolling_window_min:

  `integer(1)`. Smoothing window (epochs) for the activity signal used
  to locate episode split points. Default `15L`.

- max_split_iterations:

  `integer(1)`. Maximum recursive splits per episode. Default `3L`.

- collision_gap_h:

  `numeric(1)`. Minimum gap (hours) between two same-date main episodes
  to trigger sleep-date reassignment (Fix 26b). Default `4`.

- verbose:

  `logical(1)`. Print step-by-step diagnostics. Default `FALSE`.

## Value

A tibble with the same columns as `episodes` plus `sleep_type` (`"main"`
or `"secondary"`) and `is_nap` (logical, `TRUE` when
`sleep_type == "secondary"`, for backwards compatibility with the
`zeitr_result$nights` schema).

## Details

1.  **Fix 25** – exclude truncated episodes at the recording end.

2.  **Fix 26a** – infer adaptive nocturnal window from `int_temp` and
    `light`.

3.  **Fix 29** – split TBT 14-16 h episodes first; exclude TBT \> 16 h
    directly without attempting a split.

4.  **Fix 26c** – recover fragmented sleep nights missed by the scorer.

5.  **Rules 3-5** – classify as main or secondary using the nocturnal
    window and `min_main_tib_h`.

6.  **Fix 26b** – resolve sleep-date collisions from the noon threshold.

7.  **Rule 6** – keep the longest main episode per sleep date.

8.  **Rule 7** – exclude all episodes on days with no main sleep.

## See also

[`extract_sleep_episodes()`](https://zeitr.circadia-lab.uk/reference/extract_sleep_episodes.md),
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)
