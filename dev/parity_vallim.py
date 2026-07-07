#!/usr/bin/env python3
"""
parity_vallim.py
----------------
Runs the Vallim (JRSV) native sleep classification pipeline on a single
ActTrust file and exports two CSVs for comparison against zeitR's
run_pipeline_native() output:

  vallim_nights.csv    -- per-night: bts, gts, tbt, tst, sol, soi, waso,
                          nw, eff, sleep_type
  vallim_estimates.csv -- summary sleep metrics (overall, workday, weekend)
                          + CPD / SJL (if day-type split is available)

Usage
-----
    python parity_vallim.py \\
        --file   /path/to/recording.txt \\
        --pf     /path/to/pipeline_functions_fix27.py \\
        --output /path/to/output/dir

Pipeline steps implemented (mirrors notebook cells 3 -> 5 -> 7)
---------------------------------------------------------------
  Fix 25  : exclude truncated episodes at recording end
  Fix 26a : infer adaptive nocturnal window (TEMP + LIGHT)
  Fix 29  : split TBT 14-16 h episodes; exclude TBT > 16 h directly
  Fix 26c : recover fragmented episodes
  Rules 3-5: classify as main / secondary
  Fix 26b : resolve sleep-date collisions
  Rule 6  : keep longest main per sleep date
  Rule 7  : exclude days without main sleep

Known minor discrepancy vs R
-----------------------------
Python classification uses actual sleep ONSET time (bts + sol/60 h) for
the nocturnal window test; R uses bts directly (first epoch of the sleep
period). For episodes with short SOL (< 30 min) the difference is
negligible; longer SOL may cause isolated episodes to be classified
differently.
"""

import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd


# -- CLI -----------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('--file',   required=True,
                   help='Path to ActTrust .txt recording')
    p.add_argument('--pf',     required=True,
                   help='Path to pipeline_functions_fix27.py')
    p.add_argument('--output', default='.',
                   help='Output directory for CSV files (default: .)')
    p.add_argument('--device', default='acttrust',
                   help='Device type (default: acttrust)')
    return p.parse_args()


# -- Dynamic import of pipeline_functions --------------------------------------

def load_pipeline_functions(path):
    spec = importlib.util.spec_from_file_location('pipeline_functions_fix27', path)
    mod  = importlib.util.module_from_spec(spec)
    sys.modules['pipeline_functions_fix27'] = mod
    spec.loader.exec_module(mod)
    return mod


# -- Helpers (adapted from cell 3) ---------------------------------------------

def get_search_pct(bts, noc_start, noc_end):
    """Search window direction based on onset time (cell 3 logic)."""
    h = pd.Timestamp(bts).hour + pd.Timestamp(bts).minute / 60.0
    if h >= noc_start or h < noc_end:
        return (0.50, 0.90)   # nocturnal
    return (0.10, 0.50)       # diurnal


def sleep_date(bts):
    """Noon-threshold sleep date (cell 3 logic)."""
    t = pd.Timestamp(bts)
    return t.date() if t.hour >= 12 else (t - pd.Timedelta(days=1)).date()


def split_episode_by_activity(row, df_activity, min_frag_h,
                               rolling_min=15, search_pct=(0.10, 0.50)):
    """Split one episode at its activity peak (cell 3 logic)."""
    bts_ts    = pd.Timestamp(row['bts'])
    gts_ts    = pd.Timestamp(row['gts'])
    total_dur = (gts_ts - bts_ts).total_seconds()

    search_start = bts_ts + pd.Timedelta(seconds=total_dur * search_pct[0])
    search_end   = bts_ts + pd.Timedelta(seconds=total_dur * search_pct[1])

    mask   = (df_activity.index >= search_start) & (df_activity.index <= search_end)
    window = df_activity.loc[mask, 'activity']

    if window.empty:
        return [row]

    smoothed  = window.rolling(window=rolling_min, center=True, min_periods=1).mean()
    anchor_ts = smoothed.idxmax()

    tbt_before = (anchor_ts - bts_ts).total_seconds() / 60.0
    tbt_after  = (gts_ts - anchor_ts).total_seconds() / 60.0

    if tbt_before < min_frag_h * 60.0 or tbt_after < min_frag_h * 60.0:
        return [row]

    r_bef        = row.copy()
    r_bef['gts'] = anchor_ts
    r_bef['tbt'] = tbt_before
    r_bef['sol'] = 0.0
    r_bef['soi'] = 0.0

    r_aft        = row.copy()
    r_aft['bts'] = anchor_ts
    r_aft['tbt'] = tbt_after
    r_aft['sol'] = 0.0
    r_aft['soi'] = 0.0

    return [r_bef, r_aft]


def split_recursive(row, df_activity, min_frag_h, max_tib,
                    noc_start, noc_end, rolling_min=15, iteration=0, max_iter=3):
    if row['tbt'] / 60.0 <= max_tib or iteration >= max_iter:
        return [row]
    pct   = get_search_pct(row['bts'], noc_start, noc_end)
    frags = split_episode_by_activity(row, df_activity, min_frag_h,
                                      rolling_min=rolling_min, search_pct=pct)
    if len(frags) == 1:
        return frags
    result = []
    for f in frags:
        result.extend(split_recursive(f, df_activity, min_frag_h, max_tib,
                                      noc_start, noc_end, rolling_min,
                                      iteration + 1, max_iter))
    return result


def keep_longest_main(group):
    """Rule 6: keep the longest main episode per sleep date (cell 3 logic)."""
    mains = group[group['sleep_type'] == 'main']
    if len(mains) <= 1:
        return group
    onset_h  = pd.to_datetime(mains['bts']).dt.hour
    evening  = mains[onset_h >= 12]
    keep_idx = (evening['tbt'].idxmax() if len(evening) >= 1
                else mains['tbt'].idxmax())
    demote   = mains.index.difference([keep_idx])
    group.loc[demote, 'sleep_type'] = 'secondary'
    return group


# -- Main pipeline -------------------------------------------------------------

def run_vallim_pipeline(file_path, pf, device='acttrust'):
    """
    Run the full Vallim classification pipeline on one recording.
    Returns classified nights_data DataFrame.
    """
    from condor_pipeline.pipeline import SleepPipeline

    # Parameters (mirrors cell 3 defaults)
    max_tib_h          = 16.0
    max_main_tib_h     = 14.0
    min_main_tib_h     =  4.0
    noc_start          = 18.0   # overridden below by Fix 26a
    noc_end            =  6.0
    min_fragment_tib_h =  1.0
    rolling_window_min = 15
    max_split_iter     =  3

    print(f'Running SleepPipeline on {os.path.basename(file_path)} ...')
    pipe    = SleepPipeline(file_path, device=device)
    results = pipe.run()

    df          = results.df
    nights_data = results.nights.copy()
    print(f'  Initial episodes: {len(nights_data)}')

    # Fix 25: exclude truncated episodes at recording end
    last_day  = pd.to_datetime(df['datetime']).max().normalize()
    last_noon = last_day + pd.Timedelta(hours=12)
    trunc     = pd.to_datetime(nights_data['bts']) >= last_noon
    if trunc.any():
        print(f'  [Fix 25] Excluded {trunc.sum()} truncated episode(s) at recording end.')
        nights_data = nights_data[~trunc].reset_index(drop=True)

    # Fix 26a: infer adaptive nocturnal window
    noc_start, noc_end = pf.infer_nocturnal_window(df, nights_data, verbose=True)
    print(f'  [Fix 26a] Nocturnal window: {noc_start:.1f}h - {noc_end:.1f}h')

    # Fix 29 / Rule 1: split TBT 14-16 h episodes
    df_indexed = df.set_index('datetime') if 'datetime' in df.columns else df
    long_mask  = ((nights_data['tbt'] / 60.0 > max_main_tib_h) &
                  (nights_data['tbt'] / 60.0 <= max_tib_h))
    if long_mask.any():
        print(f'  [Fix 29] Splitting {long_mask.sum()} episode(s) with TBT 14-16 h ...')
        new_rows = []
        for _, row in nights_data.iterrows():
            if (row['tbt'] / 60.0 > max_main_tib_h) and (row['tbt'] / 60.0 <= max_tib_h):
                new_rows.extend(split_recursive(
                    row, df_indexed, min_fragment_tib_h, max_main_tib_h,
                    noc_start, noc_end, rolling_window_min, max_iter=max_split_iter))
            else:
                new_rows.append(row)
        nights_data = pd.DataFrame(new_rows).reset_index(drop=True)

    # Rule 1 cont: exclude still > max_main_tib_h after split
    before = len(nights_data)
    nights_data = nights_data[nights_data['tbt'] / 60.0 <= max_main_tib_h].reset_index(drop=True)
    if (before - len(nights_data)):
        print(f'  [Fix 29] Excluded {before - len(nights_data)} episode(s) '
              f'still > {max_main_tib_h} h after split.')

    # Rule 2: exclude > max_tib_h directly
    before = len(nights_data)
    nights_data = nights_data[nights_data['tbt'] / 60.0 <= max_tib_h].reset_index(drop=True)
    if (before - len(nights_data)):
        print(f'  [Rule 2] Excluded {before - len(nights_data)} episode(s) '
              f'with TBT > {max_tib_h} h.')

    # Sleep date column (needed for Fix 26c and Rule 6)
    nights_data['_sleep_date'] = pd.to_datetime(nights_data['bts']).apply(sleep_date)

    # Preliminary classification for Fix 26c (needs sleep_type to identify covered dates)
    bts_h_pre  = (pd.to_datetime(nights_data['bts']).dt.hour +
                  pd.to_datetime(nights_data['bts']).dt.minute / 60.0 +
                  nights_data['sol'] / 60.0)
    in_noc_pre = (bts_h_pre >= noc_start) | (bts_h_pre < noc_end)
    nights_data['sleep_type'] = 'secondary'
    nights_data.loc[in_noc_pre & (nights_data['tbt'] / 60.0 >= min_main_tib_h),
                    'sleep_type'] = 'main'

    # Fix 26c: recover fragmented episodes
    nights_data = pf.recover_fragmented_episodes(
        df, nights_data,
        temp_thresh=28.0, light_thresh=10.0, min_tib_h=3.0,
        nocturnal_start=noc_start, nocturnal_end=noc_end,
        verbose=True,
    )
    nights_data['bts'] = pd.to_datetime(nights_data['bts'])
    nights_data['gts'] = pd.to_datetime(nights_data['gts'])

    # Rules 3-5: final classification
    # NOTE: uses sleep ONSET time (bts + sol/60) — minor discrepancy vs R (uses bts directly)
    bts_h  = (pd.to_datetime(nights_data['bts']).dt.hour +
              pd.to_datetime(nights_data['bts']).dt.minute / 60.0 +
              nights_data['sol'] / 60.0)
    in_noc = (bts_h >= noc_start) | (bts_h < noc_end)
    nights_data['sleep_type'] = 'secondary'
    nights_data.loc[in_noc & (nights_data['tbt'] / 60.0 >= min_main_tib_h),
                    'sleep_type'] = 'main'

    # Re-stamp sleep_date after recovery
    nights_data['_sleep_date'] = pd.to_datetime(nights_data['bts']).apply(sleep_date)

    # Fix 26b: resolve sleep-date collisions
    nights_data = pf.fix_sleep_date_collision(nights_data, gap_h=4.0, verbose=True)

    # Rule 6: keep longest main per sleep date
    nights_data = (nights_data
                   .groupby('_sleep_date', group_keys=False)
                   .apply(keep_longest_main))

    # Rule 7: exclude days without main sleep
    days_with_main = nights_data[nights_data['sleep_type'] == 'main']['_sleep_date'].unique()
    excluded       = nights_data[
        ~nights_data['_sleep_date'].isin(days_with_main)]['_sleep_date'].unique()
    if len(excluded):
        print(f'  [Rule 7] Excluded {len(excluded)} sleep date(s) with no valid main episode.')
    nights_data = nights_data[nights_data['_sleep_date'].isin(days_with_main)]

    nights_data['nap'] = nights_data['sleep_type'] != 'main'
    nights_data = nights_data.drop(columns=['_sleep_date']).reset_index(drop=True)

    n_main = (nights_data['sleep_type'] == 'main').sum()
    n_sec  = (nights_data['sleep_type'] == 'secondary').sum()
    print(f'  Done. {n_main} main night(s), {n_sec} secondary episode(s).')
    return nights_data


# -- Main ----------------------------------------------------------------------

def main():
    args = parse_args()
    os.makedirs(args.output, exist_ok=True)

    print(f'Loading pipeline_functions from {args.pf} ...')
    pf = load_pipeline_functions(args.pf)

    nights_data = run_vallim_pipeline(args.file, pf, device=args.device)

    # Export vallim_nights.csv
    nights_out = nights_data[[
        'bts', 'gts', 'tbt', 'tst', 'sol', 'soi', 'waso', 'nw', 'eff', 'sleep_type', 'nap'
    ]].copy()
    nights_path = os.path.join(args.output, 'vallim_nights.csv')
    nights_out.to_csv(nights_path, index=False)
    print(f'\nSaved {nights_path}  ({len(nights_out)} rows)')

    # Compute and export vallim_estimates.csv
    print('\nComputing sleep estimates ...')
    estimates = {}

    sm = pf.compute_sleep_metrics(nights_data)
    estimates.update(sm)

    try:
        cpd = pf.compute_cpd_metrics(nights_data)
        estimates.update(cpd)
    except Exception as e:
        print(f'  [CPD] Skipped: {e}')

    est_df   = pd.DataFrame([estimates])
    est_path = os.path.join(args.output, 'vallim_estimates.csv')
    est_df.to_csv(est_path, index=False)
    print(f'Saved {est_path}')
    print('\nDone.')


if __name__ == '__main__':
    main()
