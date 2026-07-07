# Run the native pipeline on all files in a directory

Applies
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)
to every file matching `pattern` in `folder`, returning a named list of
`zeitr_result` objects. Files that fail are skipped with a warning.

## Usage

``` r
run_pipeline_native_batch(folder, pattern = "*.txt", ...)
```

## Arguments

- folder:

  `character(1)`. Path to a directory of ActTrust files.

- pattern:

  `character(1)`. Glob pattern. Default `"*.txt"`.

- ...:

  Additional arguments forwarded to
  [`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md).

## Value

A named list of `zeitr_result` objects.

## Examples

``` r
if (FALSE) { # \dontrun{
results <- run_pipeline_native_batch("recordings/", tz = "America/Sao_Paulo")
lapply(results, function(r) r$nights)
} # }
```
