# Read an actigraphy file into a zeitr_recording object

A device-agnostic wrapper that reads a raw actigraphy file and returns a
`zeitr_recording` object with `$epochs` (a tidy tibble) and `$metadata`
(a named list of device and recording information).

## Usage

``` r
read_actigraphy(path, device = "acttrust", tz = "UTC", ...)
```

## Arguments

- path:

  `character(1)`. Path to the raw actigraphy file.

- device:

  `character(1)`. Device type. One of `"acttrust"` (default). Additional
  devices will be added in future versions.

- tz:

  `character(1)`. Recording time zone passed to the underlying reader.
  Default is `"UTC"`.

- ...:

  Additional arguments forwarded to the device-specific reader (e.g.
  `encoding` for
  [`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md)).

## Value

A `zeitr_recording` S3 object — a named list with:

- `$epochs`:

  A tibble with one row per epoch and columns `datetime`, `activity`,
  `int_temp`, `ext_temp`, `ZCMn`, `state`, `offwrist`, `sleep`.

- `$metadata`:

  A named list with `subject`, `device_id`, `device_model`,
  `firmware_version`, `interval_s`, `source_file`, and `participant_id`
  (derived from the filename stem).

## Details

Currently supported devices:

- `"acttrust"` — Condor Instruments ActTrust / ActTrust2 (`.txt`)

## See also

[`read_actigraphy_dir()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy_dir.md)
to read a whole directory at once.

## Examples

``` r
if (FALSE) { # \dontrun{
rec <- read_actigraphy("recordings/P001.txt")
rec$epochs
rec$metadata
} # }
```
