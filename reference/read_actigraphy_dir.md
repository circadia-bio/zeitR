# Read all actigraphy files in a directory

Applies
[`read_actigraphy()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy.md)
to every file matching `pattern` in `folder`. Returns a `zeitr_study`
object — a named list of `zeitr_recording` objects, one per file. Files
that fail to parse are skipped with a warning.

## Usage

``` r
read_actigraphy_dir(
  folder,
  device = "acttrust",
  pattern = "*.txt",
  tz = "UTC",
  ...
)
```

## Arguments

- folder:

  `character(1)`. Path to a directory containing actigraphy files.

- device:

  `character(1)`. Device type passed to
  [`read_actigraphy()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy.md).
  Default is `"acttrust"`.

- pattern:

  `character(1)`. Glob pattern for file discovery. Default is `"*.txt"`.

- tz:

  `character(1)`. Recording time zone. Default is `"UTC"`.

- ...:

  Additional arguments forwarded to
  [`read_actigraphy()`](https://zeitr.circadia-lab.uk/reference/read_actigraphy.md).

## Value

A `zeitr_study` S3 object — a named list of `zeitr_recording` objects.
Names are participant IDs (filename stems).

## See also

[`study_summary()`](https://zeitr.circadia-lab.uk/reference/study_summary.md)
to summarise a `zeitr_study`.

## Examples

``` r
if (FALSE) { # \dontrun{
study <- read_actigraphy_dir("recordings/", tz = "America/Sao_Paulo")
study_summary(study)
} # }
```
