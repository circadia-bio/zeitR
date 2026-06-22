# ActTrust device parameter preset

Returns a named list of all device- and study-specific parameter values
used by the zeitR pipeline when processing ActTrust recordings. Passing
this list (or a modified copy) to
[`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
via the `params` argument fully controls which values each detector
stage receives, without having to touch the individual function calls.

## Usage

``` r
acttrust_params()
```

## Value

A named list with elements `offwrist`, `sleep`, `nap`, and `waso`.

## Details

The list is organised into four named sections that map directly to the
four pipeline stages:

- `offwrist`:

  Parameters for
  [`detect_offwrist_bimodal()`](https://zeitr.circadia-lab.uk/reference/detect_offwrist_bimodal.md).
  These values were validated against the Condor circadiaBase pipeline
  using ActTrust hardware. The refinement stage
  (`.bimodal_refine_acttrust`) is device-specific and currently only
  validated for ActTrust.

- `sleep`:

  Parameters for
  [`detect_sleep_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_sleep_crespo.md).
  `sleep_quantile` (0.365) is the ActTrust CSPD wrapper value; the
  original Crespo (2012) algorithm uses 1/3.

- `nap`:

  Parameters for
  [`detect_naps_crespo()`](https://zeitr.circadia-lab.uk/reference/detect_naps_crespo.md),
  passed as the `params` argument. Equivalent to `.cspd_nap_params()`.

- `waso`:

  Parameters for
  [`compute_waso()`](https://zeitr.circadia-lab.uk/reference/compute_waso.md).

To adapt the pipeline for a different device, copy the list, modify the
relevant values, and pass it to
[`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md):

    p <- acttrust_params()
    p$sleep$sleep_quantile <- 1/3   # use the original Crespo threshold
    result <- run_pipeline("recording.txt", params = p)

## See also

[`run_pipeline()`](https://zeitr.circadia-lab.uk/reference/run_pipeline.md)
for the pipeline entry point.
