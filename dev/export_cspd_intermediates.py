#!/usr/bin/env python3
"""
dev/export_cspd_intermediates.py

Export CSPD refiner stage intermediates for validating the R port
(R/sleep_refine.R) against the Python reference.

Run inside the circadiaBase_Docker Python environment (the jupyter Dockerfile
pins: numpy 1.26.4, scipy 1.14.1, pandas 2.2.2). It

  1. reconstructs the sleep-detection input from inst/extdata/python_output.csv
     (activity + datetime + state; the on-wrist subset is state != 4, exactly
     as cspd_wrapper does internally),
  2. runs the real CSPD pipeline with light instrumentation that dumps the
     detection vector at each refinement-stage boundary, and
  3. writes fixtures into inst/extdata/.

The vendor source is NOT modified: a copy is instrumented in a temp dir and
imported from there.

Fixtures (all on the on-wrist subset, length == sum(state != 4)):
  cspd_stage1_in.csv       detection entering the peak-valley length filter
  cspd_stage1_out.csv      detection after the peak-valley length filter
  cspd_stage2_out.csv      detection after sleep-gap separation
  cspd_pre_transitions.csv raw transitions (index, direction) pre-refinement
  cspd_refined_output.csv  final refiner output (1 = wake, 0 = sleep period)
  cspd_refined_sleep.csv   per-night bedtime/getuptime indices (post-refine)

Usage
-----
    python dev/export_cspd_intermediates.py

Override the source locations with environment variables if your checkout
differs (e.g. when running inside the Docker container):
    CONDOR_VENDOR_DIR=/path/to/vendor/condor \
    EXTDATA_DIR=/path/to/zeitR/inst/extdata \
        python dev/export_cspd_intermediates.py
"""

import os, sys, shutil, tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import warnings; warnings.filterwarnings("ignore")

# ── CONFIG — adjust paths if your checkout differs ───────────────────────────
HERE    = Path(__file__).resolve().parent                 # .../zeitR/dev
REPO    = HERE.parent                                      # .../zeitR
EXTDATA = Path(os.environ.get("EXTDATA_DIR", REPO / "inst" / "extdata"))
VENDOR  = Path(os.environ.get(
    "CONDOR_VENDOR_DIR",
    "/Users/lucas/Documents/GitHub/circadiaBase_Docker/"
    "condor_pipeline/algorithms/vendor/condor",
))
EXTDATA.mkdir(parents=True, exist_ok=True)

# ── 1. instrument a copy of the vendor sources in a temp dir ─────────────────
tmp = Path(tempfile.mkdtemp(prefix="cspd_instr_"))
for py in VENDOR.glob("*.py"):
    shutil.copy(py, tmp / py.name)

src_path = tmp / "cspd_without_prints.py"
src = src_path.read_text()
o = str(EXTDATA).replace("\\", "/")

def ins_before(text, anchor, snippet, count=1):
    if text.count(anchor) < count:
        raise RuntimeError(f"anchor not found ({text.count(anchor)}x): {anchor!r}")
    return text.replace(anchor, snippet + anchor, count)

def ins_after(text, anchor, snippet, count=1):
    if text.count(anchor) < count:
        raise RuntimeError(f"anchor not found ({text.count(anchor)}x): {anchor!r}")
    return text.replace(anchor, anchor + snippet, count)

src = ins_before(src,
    "            if self.do_peak_valley_length_filter:\n",
    f'            np.savetxt(r"{o}/cspd_stage1_in.csv", np.asarray(final_sleep_detection), fmt="%d")\n')
src = ins_before(src,
    "            datetime_difference = datetime_diff(self.datetime_stamps)\n",
    f'            np.savetxt(r"{o}/cspd_stage1_out.csv", np.asarray(final_sleep_detection), fmt="%d")\n')
src = ins_before(src,
    "            sleep_period_borders = np.diff(final_sleep_detection)\n",
    f'            np.savetxt(r"{o}/cspd_stage2_out.csv", np.asarray(final_sleep_detection), fmt="%d")\n',
    count=1)  # first occurrence == stage-2 boundary
src = ins_before(src,
    "            refined_transitions = transitions.copy()\n",
    f'            np.savetxt(r"{o}/cspd_pre_transitions.csv", np.asarray(transitions, dtype=float), '
    f'fmt="%g", delimiter=",", header="index,direction", comments="")\n')
# final assignment is the last line of the file (no trailing newline) -> append
src = ins_after(src,
    "        self.refined_output = refined_output",
    f'\n        np.savetxt(r"{o}/cspd_refined_output.csv", np.asarray(refined_output), fmt="%d")'
    f'\n        if getattr(self, "refined_sleep_df", None) is not None and len(self.refined_sleep_df):'
    f'\n            self.refined_sleep_df.to_csv(r"{o}/cspd_refined_sleep.csv", index=False)\n')

src_path.write_text(src)

# ── 2. import instrumented wrapper from temp dir ─────────────────────────────
sys.path.insert(0, str(tmp))
from cspd_wrapper_without_prints import cspd_wrapper

# ── 3. reconstruct detection input and run ───────────────────────────────────
df = pd.read_csv(EXTDATA / "python_output.csv")
df["datetime"] = pd.to_datetime(df["datetime"])
print(f"input rows: {len(df)}   on-wrist (state!=4): {int((df['state']!=4).sum())}")

cspd_wrapper(df.copy())

print("\nfixtures written to", EXTDATA)
for f in ["cspd_stage1_in.csv","cspd_stage1_out.csv","cspd_stage2_out.csv",
          "cspd_pre_transitions.csv","cspd_refined_output.csv","cspd_refined_sleep.csv"]:
    p = EXTDATA / f
    print(f"  {'ok  ' if p.exists() else 'MISS'} {f}")

shutil.rmtree(tmp, ignore_errors=True)
