# Check actigraphy timestamps for consistency issues

Scans a recording tibble for three classes of timestamp problem:

## Usage

``` r
check_consistency(x, gap_s = 120, datetime_col = "datetime")
```

## Arguments

- x:

  A tibble as returned by
  [`read_acttrust()`](https://zeitr.circadia-lab.uk/reference/read_acttrust.md)
  or
  [`prepare_actigraphy()`](https://zeitr.circadia-lab.uk/reference/prepare_actigraphy.md),
  containing a `datetime` column.

- gap_s:

  `numeric(1)`. Gap threshold in seconds. Intervals longer than this are
  flagged. Default is `120` (2 minutes).

- datetime_col:

  `character(1)`. Name of the datetime column. Default is `"datetime"`.

## Value

A tibble with one row per detected issue and columns:

- `row`:

  `integer` — row index in `x` where the issue occurs.

- `datetime`:

  `POSIXct` — timestamp at that row.

- `issue`:

  `character` — one of `"gap"`, `"backward_jump"`, or `"year_artefact"`.

- `detail`:

  `character` — human-readable description.

Returns a zero-row tibble if no issues are found.

## Details

- **Gaps** — intervals between consecutive epochs longer than `gap_s`
  seconds.

- **Backward jumps** — timestamps that go backwards in time.

- **Year artefacts** — timestamps in the years 1970 or 2000, which
  typically indicate firmware epoch-counter rollover bugs.

## Examples

``` r
if (FALSE) { # \dontrun{
rec    <- read_acttrust("recordings/P001.txt")
issues <- check_consistency(rec)
issues
} # }
```
