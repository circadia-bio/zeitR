# Run the native pipeline on all files in a directory

Applies
[`run_pipeline_native()`](https://zeitr.circadia-lab.uk/reference/run_pipeline_native.md)
to every file matching `pattern` in `folder`, returning a named list of
`zeitr_result` objects. Files that fail are skipped with a warning.

## Usage

``` r
run_pipeline_native_batch(folder, pattern = "*.txt", parallel = FALSE, ...)
```

## Arguments

- folder:

  `character(1)`. Path to a directory of ActTrust files.

- pattern:

  `character(1)`. Glob pattern. Default `"*.txt"`.

- parallel:

  `logical(1)`. If `TRUE`, processes files in parallel using
  [`future.apply::future_lapply()`](https://rdrr.io/pkg/future.apply/man/future_lapply.html)
  (a `Suggests`-only dependency). Requires a
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  to be set beforehand, e.g.
  `future::plan(future::multisession(workers = 4))`; falls back to
  sequential processing with a warning if `future.apply` is not
  installed. Default `FALSE`.

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

# Parallel: process a large cohort across 4 workers
future::plan(future::multisession(workers = 4))
results <- run_pipeline_native_batch("recordings/", tz = "America/Sao_Paulo",
                                     parallel = TRUE)
} # }
```
