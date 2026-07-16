# Read an Axivity AX3/AX6 .cwa file into a zeitR-standard epoch tibble

Bridges
[`axR::axivity_read_cwa()`](https://axr.circadia-lab.uk/reference/axivity_read_cwa.html)'s
raw per-sample output (triaxial acceleration at the device's native
sampling rate) into the same epoch-level, 9-column shape
[`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md)
produces – so an Axivity recording can flow through the rest of zeitR
([`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md),
[`compute_npcra()`](https://zeitr.circadia-lab.uk/reference/compute_npcra.md),
[`compute_sri()`](https://zeitr.circadia-lab.uk/reference/compute_sri.md),
etc.) exactly like an ActTrust one.

## Usage

``` r
read_axivity(
  path,
  tz = "UTC",
  epoch_sec = 60,
  filter_low = 0.25,
  filter_high = 2.5,
  zcm_threshold = 0.01,
  tat_threshold = 0.05
)
```

## Arguments

- path:

  `character(1)`. Path to a `.cwa`/AX6 file, forwarded to
  [`axR::axivity_read_cwa()`](https://axr.circadia-lab.uk/reference/axivity_read_cwa.html).

- tz:

  `character(1)`. Time zone the device's clock was set to. Default
  `"UTC"`.

- epoch_sec:

  `numeric(1)`. Epoch length in seconds, forwarded to
  [`compute_activity_counts()`](https://zeitr.circadia-lab.uk/reference/compute_activity_counts.md).
  Default `60`, matching the rest of zeitR.

- filter_low, filter_high:

  `numeric(1)`. Band-pass cutoffs in Hz, forwarded to
  [`compute_activity_counts()`](https://zeitr.circadia-lab.uk/reference/compute_activity_counts.md).
  Default `0.25`/`2.5` (GT3X+-style preset – see Details; no validated
  Axivity-specific preset exists).

- zcm_threshold, tat_threshold:

  `numeric(1)`. Forwarded to
  [`compute_activity_counts()`](https://zeitr.circadia-lab.uk/reference/compute_activity_counts.md).
  Defaults `0.01`/`0.05`.

## Value

A tibble with one row per epoch and the same columns as
[`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md)
(`datetime`, `activity`, `int_temp`, `ext_temp`, `ZCMn`, `light`,
`state`, `offwrist`, `sleep`), plus one extra column not present there:
`TAT` (time above threshold, seconds/epoch, from
[`compute_activity_counts()`](https://zeitr.circadia-lab.uk/reference/compute_activity_counts.md)).
`activity` is PIM. `ext_temp` is always `NA` – Axivity devices have a
single on-body temperature sensor, no separate ambient sensor. The
tibble carries a `"zeitr_axivity"` class and a `metadata` attribute (a
named list with `device_id`, `session_id`, `sample_rate`, `epoch_sec`,
`filter_low`, `filter_high`, `cwa_metadata` (axR's raw device metadata
string), `source_file`).

## What this does (and doesn't) validate

The raw-to-counts conversion itself is
[`compute_activity_counts()`](https://zeitr.circadia-lab.uk/reference/compute_activity_counts.md)
– already documented there as an unvalidated approximation of onboard
device processing (no reference converter exists to check it against).
Axivity devices additionally have **no published or validated
filter/threshold preset** at all (unlike ActTrust and GT3X+, which at
least have a documented processing description to approximate). The
`filter_low`/`filter_high` defaults here (`0.25`/`2.5` Hz) are
[`compute_activity_counts()`](https://zeitr.circadia-lab.uk/reference/compute_activity_counts.md)'s
GT3X+-style preset, reused because it's the closer starting point of the
two existing options for a research-grade wrist accelerometer like the
AX3 – **not** because it has been checked against real Axivity output.
Treat `activity`/`ZCMn` from this function as a rough approximation
only; validate against a reference (e.g. GGIR) before relying on it for
any published analysis.

`ZCMn` is
[`compute_activity_counts()`](https://zeitr.circadia-lab.uk/reference/compute_activity_counts.md)'s
raw `ZCM` count with no additional normalisation applied – named `ZCMn`
only for column-name compatibility with
[`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md)'s
CK-scoring input, not because a normalisation step has actually been
performed.

## Sampling rate

`axivity_read_cwa()` reports `sample_rate` per sample (block-level, as
stored in the `.cwa` file). This function takes the single most common
value across the whole recording and uses it for the entire conversion;
if any samples report a different rate (a genuine rate change mid
recording, or a corrupt block), a warning names the discrepancy but the
dominant rate is still used throughout. Epoch boundaries and grouping
for `light`/`int_temp` averaging are derived from this same dominant
rate, matching
[`compute_activity_counts()`](https://zeitr.circadia-lab.uk/reference/compute_activity_counts.md)'s
own epoch grouping exactly (including which trailing samples get
dropped).

## Time zone

`axivity_read_cwa()` tags `timestamp` as UTC by convention (the device's
own real-time clock, not a true UTC source) – exactly like
[`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md)
does for ActTrust's `DATE/TIME` column. `tz` here re-labels the same
clock reading under the recording's actual local time zone (via
[`lubridate::force_tz()`](https://lubridate.tidyverse.org/reference/force_tz.html);
no shift in wall-clock value), rather than converting it – set it to the
time zone the device's clock was actually set to for correct circadian
alignment downstream.

## See also

[`read_actigraphy()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy.md)
(`device = "axivity"`),
[`compute_activity_counts()`](https://zeitr.circadia-lab.uk/reference/compute_activity_counts.md)
for the underlying raw-to-counts conversion and its validation caveats,
[`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md)
for the column shape this matches.

## Examples

``` r
if (FALSE) { # \dontrun{
rec <- read_axivity("recordings/P001.cwa", tz = "America/Sao_Paulo")
rec
attr(rec, "metadata")
} # }
```
