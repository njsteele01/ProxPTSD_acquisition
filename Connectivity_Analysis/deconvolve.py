#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 21:37:27 2026
Updated to save original timeseries with simplified names.
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import gamma

# ================= Configuration =================
#BASE_VOLUME = '/mnt/munin' 
BASE_VOLUME = '/Volumes' 

# Subject List
if BASE_VOLUME == '/mnt/munin':
    subjects = [sys.argv[1]]
elif BASE_VOLUME == '/Volumes':
    subjects = [f"sub-{i:02d}" for i in range(1, 51)]
else:
    print("Error: Base volume path is incorrect or attempting to run on a new system.")
    sys.exit(1)
    
base_dir = f"{BASE_VOLUME}/Morey/ProxPTSD.01/Analysis/derivatives"

# Input Directory Pattern
timeseries_dir = os.path.join(base_dir, "rois/{subj}/acq/timeseries")
# Output Directory Pattern
output_dir = os.path.join(base_dir, "rois/gppi_seeds/gPPI_regressors")

# ================= ROI Mapping =================
# Key:   The suffix used to find the INPUT file (e.g., sub-01_SUFFIX)
# Value: The simplified name used for the OUTPUT files
roi_mapping = {
    "ventral_PFC_mask_ts_demeaned":             "ventral_PFC",
    "PAG_mask_ts_demeaned":                     "PAG",
    "post_hippo_aal_bilateral_mask_ts_demeaned":"post_hippo",
    "ant_hippo_aal_bilateral_mask_ts_demeaned": "ant_hippo",
    "amyg_bilateral_aal_ts_demeaned":           "amygdala"
}

# Master Events Path
events_dir = f"{BASE_VOLUME}/Morey/ProxPTSD.01/Analysis/behav/acquisition"

os.makedirs(output_dir, exist_ok=True)

# Parameters
TR = 2.0             
NT_bins = 16         # Microtime resolution

# Condition files
condition_files = [
    "CSminus_close_early",
    "CSminus_close_late",
    "CSminus_far_early",
    "CSminus_far_late",
    "CSplus_close_early",
    "CSplus_close_late",
    "CSplus_far_early",
    "CSplus_far_late",
    "shock"
]

# ================= User Provided Deconvolution Functions =================
def dctmtx_numpy_vect(N: int, K: int) -> np.ndarray:
    n = np.arange(N)
    C = np.zeros((N, K))
    C[:, 0] = 1 / np.sqrt(N)
    k = np.arange(1, K)
    C[:, 1:K] = np.sqrt(2 / N) * np.cos(np.pi * (2 * n[:, np.newaxis]) * k / (2 * N))
    return C

def compute_xb_Hxb(N: int, NT: int, TR: float) -> np.ndarray:
    dt = TR / NT
    k = np.arange(0, N * NT, NT)
    t = np.arange(0, 32 + dt, dt)
    hrf = gamma.pdf(t, 6) - gamma.pdf(t, NT) / 6
    hrf = hrf / np.sum(hrf)

    M = N * NT + 128
    xb = dctmtx_numpy_vect(M, N)
    Hxb = np.zeros((N, N))
    for i in range(N):
        Hx = np.convolve(xb[:, i], hrf, mode='full')
        Hxb[:, i] = Hx[k + 128]
    xb = xb[128:, :]
    return xb, Hxb

def ridge_regress_deconvolution(BOLD: np.ndarray, TR: float, alpha: float = 0.005, NT: int = 16, xb: np.ndarray = None, Hxb: np.ndarray = None) -> np.ndarray:
    if len(BOLD) == 0: raise ValueError('BOLD signal is empty.')
    if TR == 0: raise ZeroDivisionError('TR should be more than 0')
    if xb is None:
        N = len(BOLD)
        xb, Hxb = compute_xb_Hxb(N, NT, TR)
    N = len(BOLD)
    C = np.linalg.solve(Hxb.T @ Hxb + alpha * np.eye(N), Hxb.T @ BOLD)
    neuro = xb @ C
    return neuro.flatten()

# ================= Wrapper for gPPI Logic =================
def get_hrf_for_reconvolution(TR, NT):
    dt = TR / NT
    t = np.arange(0, 32 + dt, dt)
    hrf = gamma.pdf(t, 6) - gamma.pdf(t, NT) / 6
    return hrf / np.sum(hrf)

def generate_gppi_vector(roi_timeseries, onsets, durations, tr, nt=16):
    n_scans = len(roi_timeseries)
    dt = tr / nt
    
    # 1. Deconvolution
    neural_micro = ridge_regress_deconvolution(roi_timeseries, tr, NT=nt)
    
    # 2. Create Task Vector
    task_micro = np.zeros(len(neural_micro))
    for onset, dur in zip(onsets, durations):
        eff_dur = max(dur, dt) 
        
        start_idx = int(np.round(onset / dt))
        end_idx = int(np.round((onset + eff_dur) / dt))
        
        if start_idx < len(task_micro):
            end_idx = min(end_idx, len(task_micro))
            task_micro[start_idx:end_idx] = 1.0

    # 3. Interaction & Reconvolution
    ppi_neural = neural_micro * task_micro
    hrf = get_hrf_for_reconvolution(tr, nt)
    ppi_bold_micro = np.convolve(ppi_neural, hrf, mode='full')
    
    # 4. Downsample & Clean
    ppi_final = ppi_bold_micro[::nt]
    if len(ppi_final) > n_scans:
        ppi_final = ppi_final[:n_scans]
    elif len(ppi_final) < n_scans:
        ppi_final = np.pad(ppi_final, (0, n_scans - len(ppi_final)))
        
    return ppi_final

# ================= Load Master Events Once =================
print("Loading Master Event Files...")
master_events = {}

for cond in condition_files:
    ev_path = os.path.join(events_dir, f"{cond}.txt")
    if os.path.exists(ev_path):
        try:
            d = pd.read_csv(ev_path, sep=r'\s+', header=None, engine='python')
            master_events[cond] = {'onsets': d[0].values, 'durations': d[1].values}
            print(f"  Loaded {cond}: {len(d)} events")
        except Exception as e:
            print(f"  [!] Error reading {cond}: {e}")
    else:
        print(f"  [!] Missing file: {ev_path}")

if not master_events:
    raise ValueError("No event files loaded! Check paths.")

# ================= Main Subject Loop =================
for subj in subjects:
    print(f"Processing {subj}...", flush=True)
    
    subj_ts_dir = timeseries_dir.format(subj=subj)
    subj_out_dir = os.path.join(output_dir, subj)
    os.makedirs(subj_out_dir, exist_ok=True)
    
    for roi_suffix, roi_short_name in roi_mapping.items():
        
        # 1. Input Filename Construction
        ts_filename = f"{subj}_{roi_suffix}.txt"
        ts_path = os.path.join(subj_ts_dir, ts_filename)
        
        if not os.path.exists(ts_path):
            print(f"  [!] Missing ROI file: {ts_filename}. Skipping.")
            continue
            
        try:
            # Load Data
            roi_signal = np.loadtxt(ts_path)
            
            # 2. SAVE ORIGINAL TIMESERIES (Physio Regressor)
            # Saved as: {roi_short_name}.txt (e.g., PAG.txt)
            physio_out_fname = f"{roi_short_name}.txt"
            np.savetxt(os.path.join(subj_out_dir, physio_out_fname), roi_signal, fmt='%.6f')
            
        except Exception as e:
            print(f"  [!] Error loading/saving {ts_filename}: {e}")
            continue

        # 3. Generate PPI Interactions
        for cond_name, timing in master_events.items():
            try:
                ppi_vec = generate_gppi_vector(
                    roi_timeseries=roi_signal,
                    onsets=timing['onsets'],
                    durations=timing['durations'],
                    tr=TR,
                    nt=NT_bins
                )
                
                # Saved as: PPI_{roi_short_name}_{condition}.txt
                fname = f"PPI_{roi_short_name}_{cond_name}.txt"
                np.savetxt(os.path.join(subj_out_dir, fname), ppi_vec, fmt='%.6f')
                
            except Exception as e:
                print(f"  [!] Error calculating {roi_short_name} x {cond_name}: {e}")

    print(f"  -> Finished {subj}")

print("Done.")
