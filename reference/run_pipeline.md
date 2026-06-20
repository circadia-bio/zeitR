# Run the full actigraphy sleep analysis pipeline

Orchestrates the complete pipeline for a single ActTrust recording:

## Usage

``` r
run_pipeline(
  path,
  tz = "UTC",
  wake_thresh = 60L,
  gap_s = 120,
  offwrist_args = list(),
  sleep_args = list(),
  nap_args = list()
)
```

## Arguments

- path:

  `character(1)`. Path to the ActTrust `.txt` file.

- tz:

  `character(1)`. Recording time zone. Passed to
  [`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md).
  Default is `"UTC"`.

- wake_thresh:

  `integer(1)`. Minimum wake-bout length (epochs) used to separate sleep
  periods in
  [`compute_waso()`](https://zeitr.circadia-lab.uk/reference/compute_waso.md).
  Default is `60`.

- gap_s:

  `numeric(1)`. Gap threshold (seconds) for
  [`check_consistency()`](https://zeitr.circadia-lab.uk/reference/check_consistency.md).
  Default is `120`.

- offwrist_args:

  `list`. Additional arguments passed to
  [`detect_offwrist_bimodal()`](https://zeitr.circadia-lab.uk/reference/detect_offwrist_bimodal.md).
  Default is an empty list.

- sleep_args:

  `list`. Additional arguments passed to
  [`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md).
  Default is an empty list.

- nap_args:

  `list`. Additional arguments passed to
  [`detect_naps_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_naps_crespo.md).
  Default is an empty list.

## Value

A `zeitr_result` S3 object — a named list with:

- `subject_id`:

  `character` — derived from the input filename stem.

- `source_file`:

  `character` — absolute path to the input file.

- `data`:

  `tibble` — final epoch-level data frame with all state columns
  populated.

- `nights`:

  `tibble` — per-night sleep statistics.

- `issues`:

  `tibble` — timestamp consistency issues (0 rows if none).

- `metadata`:

  `list` — device and subject metadata from the file header.

## Details

1.  **Read** —
    [`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md)

2.  **Consistency check** —
    [`check_consistency()`](https://zeitr.circadia-lab.uk/reference/check_consistency.md)

3.  **Prepare** —
    [`prepare_actigraphy()`](https://zeitr.circadia-lab.uk/reference/prepare_actigraphy.md)

4.  **Off-wrist detection** —
    [`detect_offwrist_bimodal()`](https://zeitr.circadia-lab.uk/reference/detect_offwrist_bimodal.md)

5.  **Main sleep period detection** —
    [`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md)

6.  **Nap detection** —
    [`detect_naps_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_naps_crespo.md)

7.  **WASO + nightly statistics** —
    [`compute_waso()`](https://zeitr.circadia-lab.uk/reference/compute_waso.md)

## See also

[`run_pipeline_batch()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_batch.md)
for processing a directory of files.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_pipeline("recordings/P001.txt", tz = "America/Sao_Paulo")
result$nights
result$data
} # }
```
