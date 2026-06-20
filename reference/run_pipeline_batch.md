# Run the pipeline on all files in a directory

Applies
[`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
to every file matching `pattern` in `folder`, returning a list of
`zeitr_result` objects. Files that fail are skipped with a warning.

## Usage

``` r
run_pipeline_batch(folder, pattern = "*.txt", ...)
```

## Arguments

- folder:

  `character(1)`. Path to a directory containing ActTrust files.

- pattern:

  `character(1)`. Glob pattern for file discovery. Default is `"*.txt"`.

- ...:

  Additional arguments forwarded to
  [`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md).

## Value

A named list of `zeitr_result` objects, one per successfully processed
file. Names are the file stem (subject IDs).

## Examples

``` r
if (FALSE) { # \dontrun{
results <- run_pipeline_batch("recordings/", tz = "America/Sao_Paulo")
lapply(results, function(r) r$nights)
} # }
```
