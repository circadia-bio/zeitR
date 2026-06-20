# Summarise a zeitr_study across participants

Computes per-participant summary statistics from a `zeitr_study` object
(as returned by
[`read_actigraphy_dir()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy_dir.md)).
For each recording, the function computes NPCRA variables (IS, IV, RA,
L5, M10) and basic recording quality metrics.

## Usage

``` r
study_summary(study, epoch_s = NULL, L5_hours = 5, M10_hours = 10)
```

## Arguments

- study:

  A `zeitr_study` object as returned by
  [`read_actigraphy_dir()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy_dir.md),
  or a named list of `zeitr_recording` objects.

- epoch_s:

  `numeric(1)`. Epoch duration in seconds. If `NULL` (default),
  estimated separately for each recording.

- L5_hours:

  `numeric(1)`. Width of the L5 window in hours. Default is `5`.

- M10_hours:

  `numeric(1)`. Width of the M10 window in hours. Default is `10`.

## Value

A tibble with one row per participant and columns:

- `participant_id`:

  Participant identifier (filename stem).

- `n_epochs`:

  Total number of epochs in the recording.

- `n_days`:

  Recording duration in days.

- `start`:

  `POSIXct` — first epoch timestamp.

- `end`:

  `POSIXct` — last epoch timestamp.

- `IS`:

  Interdaily stability.

- `IV`:

  Intradaily variability.

- `RA`:

  Relative amplitude.

- `L5`:

  Mean activity in the least active 5 h window.

- `L5_onset`:

  Clock time of L5 midpoint.

- `M10`:

  Mean activity in the most active 10 h window.

- `M10_onset`:

  Clock time of M10 midpoint.

## See also

[`compute_npcra()`](https://zeitr.circadia-lab.uk/reference/compute_npcra.md)
for single-recording NPCRA,
[`read_actigraphy_dir()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy_dir.md)
to create a `zeitr_study`.

## Examples

``` r
if (FALSE) { # \dontrun{
study <- read_actigraphy_dir("recordings/", tz = "America/Sao_Paulo")
study_summary(study)
} # }
```
