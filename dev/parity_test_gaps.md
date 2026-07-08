# Vallim parity test gaps — status

## Gap 1 — Regenerate `vallim_nights.csv` against corrected Python

**Status: DONE**

`pipeline_functions_fix27.py` was patched: `temp_col='TEMPERATURE'` ->
`'int_temp'` and `light_col='LIGHT'` -> `'light'` in both `compute_episode_stats`
and `recover_fragmented_episodes`. Fix 26c now fires correctly in Python.

`inst/extdata/vallim_nights.csv` regenerated. All 4 Vallim parity tests pass.

---

## Gap 2 — ID0003 fixture and test

**Status: DEFERRED (out of testthat; dev artefacts ready)**

ID0003.txt cannot be committed to the repo (privacy). A permanently-skipping
test in the suite is noise, so this gap is handled as a dev script rather than
a testthat test.

Artefacts ready in `dev/`:
- `dev/vallim_nights.csv` — 7-main-night fixture generated from ID0003.txt
  against the patched Python pipeline. Confirms Fix 26c recovery of Aug 11
  sleep date (bts 22:13, gts 09:59 next day, TBT = 707 min).
- `dev/vallim_estimates.csv` — CPD metrics fixture for ID0003.txt.

To re-run manually when ID0003.txt is available locally:
```r
# In zeitR project root
devtools::load_all()
result <- run_pipeline_native(
  "path/to/ID0003.txt",
  tz = "America/Sao_Paulo", quiet = FALSE
)
stopifnot(sum(result$nights$sleep_type == "main") == 7L)
aug11 <- as.Date("2018-08-11")
main_dates <- zeitR:::.noon_date(
  result$nights$bts[result$nights$sleep_type == "main"],
  tz = "America/Sao_Paulo"
)
stopifnot(aug11 %in% main_dates)
```

---

## Gap 3 — Fix 26c unit test

**Status: DONE** (`tests/testthat/test-fix26c.R`)

Synthetic test — no external data required, runs on CI.

Scenario: a 19-h bloated CSPD state=1 period contains two CK sleep runs
(3 h + 5.5 h) separated by a 30-min warm/dark wake gap. Fix 26c merges
the CK runs -> TBT ~ 9 h. If the bug were present (state used instead of
CK scoring), TBT would be ~ 19 h and the `expect_lt(tbt, 14 * 60)` assertion
would fail.

---

## Summary

| Gap | Action                             | Status  |
|-----|------------------------------------|---------|
| 1   | Regenerate `vallim_nights.csv`     | Done    |
| 2   | ID0003 fixture + test              | Deferred (dev script) |
| 3   | Fix 26c synthetic unit test        | Done    |
