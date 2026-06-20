# Read a Condor Instruments ActTrust actigraphy file

Parses a Condor ActTrust `.txt` export into a tidy tibble. The file
format consists of a variable-length key-value header block followed by
semicolon-delimited epoch rows. The header ends at the line beginning
with `DATE/TIME`.

## Usage

``` r
read_acttrust(path, tz = "UTC", encoding = "latin1")
```

## Arguments

- path:

  `character(1)` or
  [`fs::path`](https://fs.r-lib.org/reference/path.html). Path to the
  ActTrust `.txt` file.

- tz:

  `character(1)`. Time zone string passed to
  [`lubridate::parse_date_time()`](https://lubridate.tidyverse.org/reference/parse_date_time.html).
  Defaults to `"UTC"`. Set to the local recording time zone for correct
  circadian alignment.

- encoding:

  `character(1)`. File encoding. Defaults to `"latin1"`, which matches
  Condor's default export encoding.

## Value

A tibble with one row per epoch and the following columns:

- `datetime`:

  `POSIXct` — epoch timestamp.

- `activity`:

  `double` — PIM activity count.

- `int_temp`:

  `double` — internal (on-body) temperature, °C.

- `ext_temp`:

  `double` — external (ambient) temperature, °C. `NA` if unavailable.

- `ZCMn`:

  `double` — normalised zero-crossing mode count. `NA` if unavailable.

- `state`:

  `double` — state column, initialised to `0`.

- `offwrist`:

  `double` — off-wrist indicator, initialised to `0`.

- `sleep`:

  `double` — sleep indicator, initialised to `0`.

The tibble carries a `"zeitr_recording"` class and a `metadata`
attribute (a named list with `subject`, `device_id`, `device_model`,
`firmware_version`, `interval_s`, `source_file`).

## Examples

``` r
if (FALSE) { # \dontrun{
rec <- read_acttrust("recordings/P001.txt")
rec
attr(rec, "metadata")
} # }
```
