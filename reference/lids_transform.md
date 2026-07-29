# Apply the LIDS (Locomotor Inactivity During Sleep) transform

Converts an activity signal into "inactivity" via the LIDS non-linear
transform, then smooths it – the first step of the LIDS ultradian-rhythm
pipeline (Winnebeck et al. 2018; Hammad et al. 2026). LIDS is only ever
computed *inside* an already-identified sleep bout; extract bouts first
via
[`compute_lids()`](https://zeitr.circadia-lab.uk/reference/compute_lids.md)
or
[`detect_lids_bouts()`](https://zeitr.circadia-lab.uk/reference/detect_lids_bouts.md).

## Usage

``` r
lids_transform(
  activity,
  epoch_min = 1,
  method = c("gaussian", "mva"),
  win_min = 30,
  sigma_min = 5
)
```

## Arguments

- activity:

  Numeric vector of activity counts for a single sleep bout, in
  chronological order at a constant epoch length.

- epoch_min:

  `numeric(1)`. Epoch duration in minutes. Default `1`.

- method:

  `character(1)`. `"gaussian"` (default) or `"mva"`.

- win_min:

  `numeric(1)`. Smoothing window width in minutes (full width for
  `"mva"`; +/- 3 sigma width for `"gaussian"`). Default `30`.

- sigma_min:

  `numeric(1)`. Gaussian kernel standard deviation in minutes, used only
  when `method = "gaussian"`. Default `5` (Hammad et al. 2026); ignored
  for `"mva"`.

## Value

Numeric vector of smoothed LIDS values, same length as `activity`.

## Details

\$\$\text{LIDS}\_i = \frac{100}{1+x_i}\$\$

where \\x_i\\ is the raw activity count at epoch \\i\\. A LIDS value of
100 means zero movement; it falls toward 0 as movement increases.

Two smoothing methods are available, both over a nominal 30-min window:

- `"gaussian"` (default) – Gaussian kernel, standard deviation
  `sigma_min` (default 5 min), truncated at +/- 3 sigma (Hammad et al.
  2026).

- `"mva"` – centered moving average (Winnebeck et al. 2018;
  `pyActigraphy`'s default). Uses zeitR's border-replicated internal
  `rolling_mean_cpp()`, which replicates the edge value rather than
  shrinking the window near the bout boundary (pandas' `min_periods=1`
  behaviour) – a minor difference confined to the first/ last ~15 min of
  each bout.

Any `NA` in `activity` (e.g. a brief off-wrist gap inside an otherwise
valid bout) is linearly interpolated first –
[`fit_lids()`](https://zeitr.circadia-lab.uk/reference/fit_lids.md)
cannot handle missing values.

## References

Winnebeck, E. C., Fischer, D., Leise, T., & Roenneberg, T. (2018).
Dynamics and Ultradian Structure of Human Sleep in Real Life. *Current
Biology*, 28(1), 49-59.e5.
[doi:10.1016/j.cub.2017.11.063](https://doi.org/10.1016/j.cub.2017.11.063)

Hammad, G., Schoch, S. F., Engelmann, M., Spock, Z., Kurth, S., &
Winnebeck, E. C. (2026). Charting infant sleep cycle development using
actigraphy: Longitudinal evidence for ultradian cycle lengthening within
the first year of life. *SLEEP*.

## See also

[`fit_lids()`](https://zeitr.circadia-lab.uk/reference/fit_lids.md),
[`detect_lids_bouts()`](https://zeitr.circadia-lab.uk/reference/detect_lids_bouts.md),
[`compute_lids()`](https://zeitr.circadia-lab.uk/reference/compute_lids.md)

## Examples

``` r
set.seed(1)
activity <- pmax(0, 20 + 15 * sin(seq(0, 6 * pi, length.out = 360)) +
                    rnorm(360, sd = 5))
lids <- lids_transform(activity, method = "gaussian")
```
