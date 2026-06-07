import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tdt
from matplotlib.patches import Ellipse
from matplotlib.gridspec import GridSpec
from scipy.signal import butter, filtfilt, savgol_filter, medfilt, detrend, welch, iirnotch
from scipy.integrate import trapezoid
from scipy.stats import ttest_ind
from sklearn.neighbors import KDTree

# ============================================
# 1. USER CONFIGURATION
# Personal sessions: copy local_config.example.py → local_config.py and fill it in.
# local_config.py is gitignored and overrides every default below.
# ============================================
block_path     = "PASTE_YOUR_PATH_HERE"
start_time_s   = 0
end_time_s     = 10000
stim_intervals = []
store          = "EEGw"
channels       = [1, 2, 3, 4]  # Chest=1–2, Ear=3–4
outdir         = None           # set below after possible override

try:
    from local_config import *  # noqa: F401, F403
except ImportError:
    pass

if outdir is None:
    outdir = f"{block_path}/analysis_output_Rodent_I"
os.makedirs(outdir, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 600,
    "savefig.dpi": 600,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.size": 20
})

HRV_UNITS = {
    "mean_rr":"ms","mean_hr":"bpm","SDNN":"ms","RMSSD":"ms",
    "pnn3":"%","pnn5":"%","LF":"ms^2","HF":"ms^2",
    "total_power":"ms^2","LF_HF":"ratio",
    "sd1":"ms","sd2":"ms","sd1_sd2_ratio":"ratio",
    "sample_entropy":"dim","dfa_alpha_I":"dim",
    "dfa_alpha_II":"dim","dfa_alpha_ratio":"ratio"
}

# ============================================
# 2. HELPER FUNCTIONS
# ============================================

def save_plot(fig, path):
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)

def auto_invert(data, fs, channels):
    data_adj = data.copy()
    for i in range(data.shape[0]):
        pos = np.sum(data_adj[i][data_adj[i] > 0] ** 2)
        neg = np.sum(data_adj[i][data_adj[i] < 0] ** 2)
        if abs(neg) > abs(pos):
            data_adj[i] = -data_adj[i]
    
    rms1, rms2 = np.sqrt(np.mean(data_adj[0]**2)), np.sqrt(np.mean(data_adj[1]**2))
    master_idx = 0 if rms1 > rms2 else 1
    master = data_adj[master_idx]

    for i in range(data.shape[0]):
        if i != master_idx and np.corrcoef(master, data_adj[i])[0,1] < 0:
            data_adj[i] = -data_adj[i]
    return data_adj

def clean_ecg_signal(data, fs):
    b_notch, a_notch = iirnotch(60.0, 30.0, fs)
    data_notched = filtfilt(b_notch, a_notch, data)
    b_hi, a_hi = butter(2, 0.5/(fs/2), btype='high')
    return filtfilt(b_hi, a_hi, data_notched)

def clean_rr_for_hrv(rr, low=0.10, high=0.35, max_rel_dev=0.20, reference_rr=None):
    """Filter RR intervals to a physiological band, then reject outliers relative to
    a central rate. By default the central rate is the local median (self-referenced).

    When reference_rr is supplied (the simultaneous chest median for this window),
    it anchors the deviation test instead — this rejects ear T-wave false detections
    that the local median would otherwise lock onto, using the trusted lead's known
    heart rate as a physiological prior. A wider band is allowed in this mode to
    preserve genuine beat-to-beat variation around the anchored rate."""
    rr = np.asarray(rr)
    rr_f = rr[(rr >= low) & (rr <= high)]
    if len(rr_f) < 5:
        return rr_f
    if reference_rr is not None:
        med = reference_rr
        max_rel_dev = max(max_rel_dev, 0.30)
    else:
        med = np.median(rr_f)
    return rr_f[np.abs(rr_f - med) / med <= max_rel_dev]

def dfa_alpha(rr):
    x = rr - np.mean(rr)
    y = np.cumsum(x)
    N = len(x)
    n_vals = np.unique(np.logspace(np.log10(4), np.log10(max(N // 4, 5)), 12).astype(int))
    F = []
    valid_n = []
    for n in n_vals:
        if n < 3:
            continue
        segs = N // n
        if segs < 2:
            continue
        rms = []
        for i in range(segs):
            seg = y[i*n:(i+1)*n]
            t = np.arange(n)
            p = np.polyfit(t, seg, 1)
            rms.append(np.sqrt(np.mean((seg - (p[0]*t + p[1]))**2)))
        F.append(np.mean(rms))
        valid_n.append(n)
    n_vals = np.array(valid_n)
    F = np.array(F)
    if len(F) < 3:
        return np.nan, np.nan
    short = n_vals <= 16
    long_ = n_vals > 16
    a1 = np.polyfit(np.log(n_vals[short]), np.log(F[short]), 1)[0] if short.sum() > 2 else np.nan
    a2 = np.polyfit(np.log(n_vals[long_]), np.log(F[long_]), 1)[0] if long_.sum() > 2 else np.nan
    return a1, a2

def qrs_filter(ecg, fs, bp_band=(0.1,200), env_smooth_ms=4, mad_k=1.5, refractory_ms=80, min_rr_ms=40, max_rr_ms=400):
    x = ecg - np.mean(ecg)
    nyq = fs/2
    b,a = butter(2,[bp_band[0]/nyq, bp_band[1]/nyq], btype="bandpass")
    y = filtfilt(b,a,x)
    y /= (np.percentile(np.abs(y),97) + 1e-12)

    env = np.abs(y)
    win = max(int((env_smooth_ms/1000)*fs),3)
    env = np.convolve(env, np.ones(win)/win, mode="same")

    win_med = int(0.15*fs) | 1
    medv = medfilt(env, win_med)
    madv = medfilt(np.abs(env-medv), win_med)
    thr  = medv + mad_k*madv

    diff = np.diff((env > thr).astype(int))
    starts = np.where(diff==1)[0]+1
    cand = []
    for s in starts:
        idx = np.argmax(y[s:s+int(0.01*fs)]) + s
        if y[idx] > 0: cand.append(idx)

    Rl = []
    last = -1e9
    min_samples = int((refractory_ms/1000)*fs)
    for c in sorted(cand):
        if c - last >= min_samples:
            Rl.append(c)
            last = c
        elif abs(y[c]) > abs(y[Rl[-1]]):
            Rl[-1] = c
            last = c
            
    tpk = np.array(Rl)/fs
    if len(tpk) < 2: return np.array([]), np.array([]), np.array(Rl), tpk
    rr = np.diff(tpk)
    keep = np.r_[True, (rr*1000 >= min_rr_ms) & (rr*1000 <= max_rr_ms)]
    return rr[keep[1:]], 60/rr[keep[1:]], np.array(Rl)[keep], tpk[keep]

# ============================================
# 3. HRV METRICS
# ============================================

def hrv_time_domain(rr):
    rr_ms = rr*1000
    if len(rr_ms)<3: return {k: np.nan for k in ["mean_rr","mean_hr","SDNN","RMSSD","pnn3","pnn5"]}
    diff = np.diff(rr_ms)
    return dict(mean_rr=np.mean(rr_ms), mean_hr=60000/np.mean(rr_ms), SDNN=np.std(rr_ms,ddof=1),
                RMSSD=np.sqrt(np.mean(diff**2)), pnn3=100*np.mean(np.abs(diff)>3), pnn5=100*np.mean(np.abs(diff)>5))

def hrv_frequency_domain(rr, fs_interp=4.0):
    if len(rr)<4: return {"LF":np.nan,"HF":np.nan,"total_power":np.nan,"LF_HF":np.nan,"f":np.array([]), "psd":np.array([])}
    t = np.cumsum(rr) - rr[0]
    rr_interp = np.interp(np.arange(0, t[-1], 1/fs_interp), t, rr)
    f, px = welch(detrend(rr_interp-np.mean(rr_interp)), fs=fs_interp, nperseg=min(256,len(rr_interp)))
    px *= 1e6
    lf_m = (f>=0.2)&(f<=0.75)
    hf_m = (f>=0.75)&(f<=3.0)
    LF, HF = trapezoid(px[lf_m], f[lf_m]), trapezoid(px[hf_m], f[hf_m])
    return dict(LF=LF, HF=HF, total_power=trapezoid(px[(f>=0.2)&(f<=3.0)], f[(f>=0.2)&(f<=3.0)]), LF_HF=LF/HF if HF>0 else np.nan, f=f, psd=px)

def sample_entropy(x, m=2, r=None):
    x = np.asarray(x)
    if len(x)<m+2: return np.nan
    r = 0.2*np.std(x) if r is None else r
    Xm, Xm1 = [np.array([x[i:i+k] for i in range(len(x)-k)]) for k in [m, m+1]]
    Cm = np.sum(KDTree(Xm, metric="chebyshev").query_radius(Xm, r, count_only=True)-1)
    Cm1 = np.sum(KDTree(Xm1, metric="chebyshev").query_radius(Xm1, r, count_only=True)-1)
    return -np.log(Cm1/Cm) if Cm*Cm1 > 0 else np.nan

def compute_hrv(rr):
    out = hrv_time_domain(rr)
    fd = hrv_frequency_domain(rr)
    out.update({k: fd[k] for k in ["LF","HF","total_power","LF_HF"]})
    diff_ms = np.diff(rr*1000)
    out["sd1"] = np.sqrt(0.5)*np.std(diff_ms, ddof=1)
    out["sd2"] = np.sqrt(2*np.std(rr*1000, ddof=1)**2 - 0.5*np.std(diff_ms, ddof=1)**2)
    out["sd1_sd2_ratio"] = out["sd1"] / out["sd2"] if out["sd2"] else np.nan
    out["sample_entropy"] = sample_entropy(rr)
    a1, a2 = dfa_alpha(rr)
    out["dfa_alpha_I"]    = a1
    out["dfa_alpha_II"]   = a2
    out["dfa_alpha_ratio"] = a1 / a2 if (a1 and a2) else np.nan
    return out

# ============================================
# 3b. SEQUENTIAL PRE/POST (FLANKING-REST) ANALYSIS
# ============================================
# For each stimulation window, compute HRV over the rest period immediately
# *before* stim onset and immediately *after* stim offset. Each slice is taken
# from the flanking non-stim (rest/washout) window: up to SEQ_SLICE_S seconds
# ending at onset, and the same span starting at offset. Comparing the two
# isolates the carryover effect of stimulation between two quiescent baselines,
# with no during-stim artifact in either slice. Because each slice is drawn from
# one end of a rest window, the "after" slice of one stim and the "before" slice
# of the next come from opposite ends of the same washout — so consecutive
# brackets share no data and the paired differences stay independent.

SEQ_SLICE_S   = 300.0   # flanking rest slice length (5 minutes)
SEQ_MIN_S     = 120.0   # discard a slice shorter than this (rest window too brief)
SEQ_MIN_BEATS = 30      # discard a slice with too few cleaned beats

def _slice_hrv(tpk, win, ref=None):
    """HRV over peaks within the [start, end] time slice, with RR cleaning
    (ear anchored to the chest reference when ref is supplied)."""
    m = (tpk >= win[0]) & (tpk <= win[1])
    if np.sum(m) <= SEQ_MIN_BEATS:
        return None, int(np.sum(m))
    rr = clean_rr_for_hrv(np.diff(tpk[m]), reference_rr=ref)
    if len(rr) < SEQ_MIN_BEATS:
        return None, len(rr)
    return compute_hrv(rr), len(rr)

def sequential_prepost_rows(stim_intervals, tpk_c, tpk_e):
    rows = []
    for i, ev in enumerate(stim_intervals):
        if "_On" not in ev["label"]:
            continue
        prev_ev = stim_intervals[i - 1] if i > 0 else None
        next_ev = stim_intervals[i + 1] if i + 1 < len(stim_intervals) else None

        slices = {}
        if prev_ev is not None and "_On" not in prev_ev["label"]:
            b0 = max(ev["start"] - SEQ_SLICE_S, prev_ev["start"])
            if ev["start"] - b0 >= SEQ_MIN_S:
                slices["before"] = (b0, ev["start"])
        if next_ev is not None and "_On" not in next_ev["label"]:
            a1 = min(ev["end"] + SEQ_SLICE_S, next_ev["end"])
            if a1 - ev["end"] >= SEQ_MIN_S:
                slices["after"] = (ev["end"], a1)

        for phase, win in slices.items():
            # Chest first: its cleaned median anchors ear cleaning for this slice.
            mc = (tpk_c >= win[0]) & (tpk_c <= win[1])
            chest_rr = clean_rr_for_hrv(np.diff(tpk_c[mc])) if np.sum(mc) > 5 else np.array([])
            chest_ref = float(np.median(chest_rr)) if len(chest_rr) >= 5 else None
            for site, tpk in [("CHEST", tpk_c), ("EAR", tpk_e)]:
                hrv, nb = _slice_hrv(tpk, win, ref=chest_ref if site == "EAR" else None)
                if hrv is None:
                    continue
                rows.append(dict(
                    stim=ev["label"], phase=phase, site=site,
                    slice_start=round(win[0], 1), slice_end=round(win[1], 1),
                    slice_dur_s=round(win[1] - win[0], 1), n_beats=nb,
                    **{k: hrv.get(k, np.nan) for k in HRV_UNITS}))
    return rows

# ============================================
# 4. EXECUTION PIPELINE
# ============================================

def run_trial(block_path, start_time_s=0, end_time_s=10000, stim_intervals=None,
              store="EEGw", channels=None, outdir=None, label=None):
    if stim_intervals is None:
        stim_intervals = []
    if channels is None:
        channels = [1, 2, 3, 4]
    if outdir is None:
        outdir = f"{block_path}/analysis_output_Rodent_I"
    if label is None:
        label = os.path.basename(block_path)

    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)

    print(f"\n=== {label} ===")
    print("Loading data...")
    blk = tdt.read_block(block_path, store=store)
    fs = int(blk.streams[store].fs)
    X_adj = auto_invert(np.asarray(blk.streams[store].data)[[i-1 for i in channels]], fs, channels)

    t_all = np.arange(X_adj.shape[1]) / fs
    mask  = (t_all >= start_time_s) & (t_all <= end_time_s)
    CHEST = clean_ecg_signal(np.mean(X_adj[0:2, mask], axis=0), fs)
    EAR   = clean_ecg_signal(np.mean(X_adj[2:4, mask], axis=0), fs)

    print("Detecting peaks...")
    rr_c, hr_c, Rl_c, tpk_c = qrs_filter(CHEST, fs)
    rr_e, hr_e, Rl_e, tpk_e = qrs_filter(EAR,   fs)

    # ── Diagnostic peak-detection panel ──────────────────────────────────────
    # 6 rows: for each of Chest and Ear — raw ECG + bandpass-filtered signal
    #         with detected peaks, shown in two windows: early in recording and
    #         at first stim interval. Then RR timeseries and histogram.
    t_sig = np.arange(len(CHEST)) / fs

    # Compute the bandpass-filtered signal the QRS detector actually uses.
    def _bandpass(sig, fs, lo=0.1, hi=200):
        from scipy.signal import butter, filtfilt
        nyq = fs / 2
        b, a = butter(2, [lo/nyq, hi/nyq], btype="bandpass")
        y = filtfilt(b, a, sig - np.mean(sig))
        y /= (np.percentile(np.abs(y), 97) + 1e-12)
        return y

    bp_c = _bandpass(CHEST, fs)
    bp_e = _bandpass(EAR,   fs)

    # Two 10-second windows: start of recording and start of first stim interval.
    early_start = start_time_s + 5.0
    stim_start  = stim_intervals[0]["start"] if stim_intervals else early_start + 60.0
    windows = [
        (early_start,       early_start + 10.0, "Early baseline"),
        (stim_start,        stim_start  + 10.0, "First stim onset"),
    ]

    fig = plt.figure(figsize=(22, 26))
    fig.suptitle(f"Peak Detection Diagnostic — {label}", fontsize=18)
    gs = fig.add_gridspec(6, 2, hspace=0.5, wspace=0.35)

    # Mark which peaks survive clean_rr_for_hrv so the diagnostic colours match the
    # cleaner applied before HRV computation. Ear is anchored to the chest median
    # (as in the per-window table logic), using a global chest median here for the
    # whole-recording overview.
    def _cleaned_peak_mask(tpk, reference_rr=None, low=0.10, high=0.35, max_rel_dev=0.20):
        if len(tpk) < 2:
            return np.zeros(len(tpk), dtype=bool)
        rr = np.diff(tpk)
        in_bounds = (rr >= low) & (rr <= high)
        if reference_rr is not None:
            med = reference_rr
            max_rel_dev = max(max_rel_dev, 0.30)
        else:
            med = np.median(rr[in_bounds]) if in_bounds.sum() >= 5 else None
        if med is not None:
            clean = in_bounds & (np.abs(rr - med) / med <= max_rel_dev)
        else:
            clean = in_bounds
        used = np.zeros(len(tpk), dtype=bool)
        for i in np.where(clean)[0]:
            used[i]     = True
            used[i + 1] = True
        return used

    chest_global_ref = float(np.median(clean_rr_for_hrv(rr_c))) if len(rr_c) > 10 else None
    used_c = _cleaned_peak_mask(tpk_c)
    used_e = _cleaned_peak_mask(tpk_e, reference_rr=chest_global_ref)

    for row, (site, raw_sig, bp_sig, tpk, used, color) in enumerate([
        ("Chest", CHEST, bp_c, tpk_c, used_c, "steelblue"),
        ("Ear",   EAR,   bp_e, tpk_e, used_e, "firebrick"),
    ]):
        for col, (ws, we, wlabel) in enumerate(windows):
            ax  = fig.add_subplot(gs[row*2,     col])
            ax2 = fig.add_subplot(gs[row*2 + 1, col], sharex=ax)

            wm  = (t_sig >= ws) & (t_sig <= we)
            pm  = (tpk   >= ws) & (tpk   <= we)
            pm_used     = pm & used
            pm_filtered = pm & ~used

            for ax_i, sig_i in [(ax, raw_sig), (ax2, bp_sig)]:
                ax_i.plot(t_sig[wm], sig_i[wm], color=color, lw=0.8, alpha=0.7)
                if pm_filtered.any():
                    filt_idx = (tpk[pm_filtered] * fs).astype(int).clip(0, len(sig_i)-1)
                    ax_i.plot(tpk[pm_filtered], sig_i[filt_idx],
                              'x', color='gray', ms=7, mew=1.5, zorder=4, label="Filtered out")
                if pm_used.any():
                    used_idx = (tpk[pm_used] * fs).astype(int).clip(0, len(sig_i)-1)
                    ax_i.plot(tpk[pm_used], sig_i[used_idx],
                              'v', color='black', ms=7, zorder=5, label="Used for HRV")

            rr_w = np.diff(tpk[pm_used]) * 1000
            n_used = pm_used.sum()
            n_filt = pm_filtered.sum()
            mrr = f"{rr_w.mean():.0f} ms ({60000/rr_w.mean():.0f} bpm)" if len(rr_w) else "—"
            ax.set_title(f"{site} — {wlabel} (raw ECG)  |  ▼ used: {n_used}  ✗ filtered: {n_filt}", fontsize=9)
            ax.set_ylabel("Amplitude")
            ax.legend(fontsize=8, loc="upper right")
            plt.setp(ax.get_xticklabels(), visible=False)
            ax2.set_title(f"Bandpass  |  mean RR (used) = {mrr}", fontsize=9)
            ax2.set_xlabel("Time (s)")
            ax2.set_ylabel("Norm. amplitude")

    # RR interval timeseries (row 4)
    ax = fig.add_subplot(gs[4, :])
    if len(rr_c): ax.plot(tpk_c[1:], rr_c*1000, color="steelblue", lw=0.6, alpha=0.7, label="Chest")
    if len(rr_e): ax.plot(tpk_e[1:], rr_e*1000, color="firebrick", lw=0.6, alpha=0.7, label="Ear")
    for event in stim_intervals:
        ax.axvspan(event["start"], event["end"],
                   color="red" if "_On" in event["label"] else "gray", alpha=0.08)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("RR interval (ms)")
    ax.set_title("RR Interval Timeseries — full recording")
    ax.legend(loc="upper right", fontsize=11)

    # RR histogram (row 5)
    ax = fig.add_subplot(gs[5, :])
    bins = np.linspace(50, 600, 80)
    if len(rr_c): ax.hist(rr_c*1000, bins=bins, color="steelblue", alpha=0.6, label="Chest")
    if len(rr_e): ax.hist(rr_e*1000, bins=bins, color="firebrick", alpha=0.6, label="Ear")
    ax.axvspan(100, 350, color="green", alpha=0.08, label="clean_rr_for_hrv bounds")
    ax.set_xlabel("RR interval (ms)")
    ax.set_ylabel("Count")
    ax.set_title("RR Interval Distribution")
    ax.legend(loc="upper right", fontsize=11)

    save_plot(fig, f"{outdir}/Diagnostic_Peak_Detection.png")
    print("  Diagnostic panel saved.")
    # ─────────────────────────────────────────────────────────────────────────

    rows = []
    for event in stim_intervals:
        # Chest first: its cleaned median is this window's heart-rate reference,
        # used to anchor ear cleaning (rejects ear T-wave double-detections).
        mc = (tpk_c >= event["start"]) & (tpk_c <= event["end"])
        chest_rr = clean_rr_for_hrv(np.diff(tpk_c[mc])) if np.sum(mc) > 5 else np.array([])
        chest_ref = float(np.median(chest_rr)) if len(chest_rr) >= 5 else None

        for site, tpk in [("CHEST", tpk_c), ("EAR", tpk_e)]:
            m = (tpk >= event["start"]) & (tpk <= event["end"])
            if np.sum(m) > 5:
                ref = chest_ref if site == "EAR" else None
                hrv = compute_hrv(clean_rr_for_hrv(np.diff(tpk[m]), reference_rr=ref))
                rows.append(dict(window=event["label"], site=site, **{k: hrv.get(k, np.nan) for k in HRV_UNITS}))

    pd.DataFrame(rows).to_csv(f"{outdir}/Stimulation_Analysis_Summary.csv", index=False)

    # Sequential pre/post analysis: HRV of the 5-min rest slices flanking each stim.
    seq_rows = sequential_prepost_rows(stim_intervals, tpk_c, tpk_e)
    pd.DataFrame(seq_rows).to_csv(f"{outdir}/Sequential_PrePost_Summary.csv", index=False)
    print(f"  Sequential pre/post slices: {len(seq_rows)} rows.")

    fig, ax = plt.subplots(figsize=(15, 5))
    ax.plot(tpk_c[1:], hr_c, 'k', alpha=0.3, label="Chest HR")
    for event in stim_intervals:
        ax.axvspan(event["start"], event["end"],
                   color='red' if "On" in event["label"] else 'gray', alpha=0.15)
    ax.set_title(f"Global Heart Rate & Stimulation Overlays — {label}")
    save_plot(fig, f"{outdir}/Global_HR_Summary.png")

    print(f"Results saved to: {outdir}")


# If local_config defines a `trials` list, run each; otherwise fall back to single-block defaults.
try:
    trials
except NameError:
    trials = [dict(block_path=block_path, start_time_s=start_time_s, end_time_s=end_time_s,
                   stim_intervals=stim_intervals, store=store, channels=channels, outdir=outdir)]

for trial in trials:
    run_trial(**trial)

print("\nAll trials complete.")




