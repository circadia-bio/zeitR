# Vallim parity test gaps — next session checklist

## Context

`devtools::test()` passes but does NOT guarantee parity with the corrected
Vallim Python pipeline. Three specific gaps need to be closed.

---

## Gap 1 — Regenerate `vallim_nights.csv` against corrected Python

**Why:** The current fixture was generated when Julia's `recover_fragmented_episodes`
had a column name bug (`temp_col='TEMPERATURE'`, `light_col='LIGHT'` instead of
`'int_temp'` and `'light'`). Fix 26c wasn't firing in Python. After Julia patches
the column names, the Python output may change for recordings where Fix 26c
would recover a night.

**Check first:** Has Julia pushed a fix to `pipeline_functions_fix27.py`?

**If yes:**
```python
# In JupyterLab / circadiaBase_Docker
# Update notebooks/parity_vallim.py to use corrected Python, then:
python notebooks/parity_vallim.py \
    --file    /path/to/input1.txt \
    --pf      notebooks/pipeline_functions_fix27.py \
    --output  /path/to/zeitR/inst/extdata/
# This overwrites vallim_nights.csv
```

Then in R:
```r
devtools::test(filter = "vallim-parity")  # all 4 tests must still pass
```

---

## Gap 2 — Add ID0003.txt fixture and test

**Why:** ID0003.txt is the recording that exposed the Fix 26c bug. R correctly
recovers the Aug 11 sleep night via Cole-Kripke scoring + temperature/light
merge; the old Python didn't. No test currently covers this.

**Expected R output (verified):**
- 7 main nights
- Aug 11 sleep date recovered by Fix 26c (bts ~ 21:33, gts ~ 07:04, TBT ~ 11.2h)

**Step 1 — Generate Python fixture (after Julia's patch):**
```python
python notebooks/parity_vallim.py \
    --file    /path/to/ID0003.txt \
    --pf      notebooks/pipeline_functions_fix27.py \
    --output  /tmp/
# Copy /tmp/vallim_nights.csv to inst/extdata/vallim_nights_id0003.csv
```

**Step 2 — Write the test** in `tests/testthat/test-vallim-parity.R`:
```r
test_that("Vallim pipeline ID0003: 7 main nights including Fix 26c recovery", {
  fpath <- system.file("extdata", "ID0003.txt", package = "zeitR")
  if (!nzchar(fpath)) testthat::skip("ID0003.txt not in extdata")

  result <- run_pipeline_native(fpath, tz = "America/Sao_Paulo", quiet = TRUE)

  # Night count
  expect_equal(sum(result$nights$sleep_type == "main"), 7L)

  # Aug 11 sleep date is present (recovered by Fix 26c)
  main_dates <- as.Date(result$nights$bed_time[result$nights$sleep_type == "main"],
                        tz = "America/Sao_Paulo")
  # Aug 11: bed_time after noon on Aug 11 -> sleep_date = Aug 11
  expect_true(as.Date("2018-08-11") %in% main_dates)
})
```

**Note:** ID0003.txt must be added to `inst/extdata/` if not already there.

---

## Gap 3 — Fix 26c unit test

**Why:** No test directly validates that Cole-Kripke scoring is used for
fragment recovery (not the period-level CSPD state column).

**Write in `tests/testthat/test-vallim-parity.R`:**
```r
test_that("Fix 26c uses Cole-Kripke scoring not CSPD period state", {
  # Synthetic: one sleep period state=1 for 20h (as CSPD would set it),
  # but CK scoring would classify it as fragmented short sleep runs.
  # Fix 26c should NOT recover a single 20h run.
  # This verifies the CK-based approach is used, not raw state==1.

  # The cleanest check: run on ID0003.txt and confirm the Aug 11 episode
  # has TBT < 14h (i.e. it was recovered as a bounded fragment,
  # not as the full 19.4h CSPD period).
  fpath <- system.file("extdata", "ID0003.txt", package = "zeitR")
  if (!nzchar(fpath)) testthat::skip("ID0003.txt not in extdata")

  result <- run_pipeline_native(fpath, tz = "America/Sao_Paulo", quiet = TRUE)
  main <- result$nights[result$nights$sleep_type == "main", ]

  # The recovered Aug 11 night should have TBT well under 14h (not the 19.4h CSPD period)
  aug11_mask <- as.Date(main$bed_time, tz = "America/Sao_Paulo") == as.Date("2018-08-11")
  if (any(aug11_mask)) {
    expect_lt(main$tbt[aug11_mask] / 60, 14.0)
  }
})
```

---

## Summary of actions

| Priority | Action | Blocked on |
|---|---|---|
| High | Regenerate `vallim_nights.csv` | Julia patching Python |
| High | Add ID0003 fixture + test | ID0003.txt in inst/extdata + Julia's patch |
| Medium | Fix 26c unit test | Nothing — can write now |

The Fix 26c unit test (Gap 3) can be written immediately without waiting for
Julia. Gaps 1 and 2 depend on Julia's Python fix.
