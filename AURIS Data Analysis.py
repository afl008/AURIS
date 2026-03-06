import os
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
# 1. USER CONFIGURATION (Replace with your session info)
# ============================================
block_path   = "PASTE_YOUR_PATH_HERE"
start_time_s = 0
end_time_s   = 10000 
stim_intervals = [
    # Example format: {"label": "Baseline", "start": 0, "end": 60},
]

# Analysis settings
store      = "EEGw"
channels   = [1, 2, 3, 4]  # Chest=1–2, Ear=3–4
outdir     = f"{block_path}/analysis_output_Rodent_I"
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
    out["sample_entropy"] = sample_entropy(rr)
    return out

# ============================================
# 4. EXECUTION PIPELINE
# ============================================

print("Loading data...")
blk = tdt.read_block(block_path, store=store)
fs = int(blk.streams[store].fs)
X_adj = auto_invert(np.asarray(blk.streams[store].data)[[i-1 for i in channels]], fs, channels)

t_all = np.arange(X_adj.shape[1])/fs
mask = (t_all >= start_time_s) & (t_all <= end_time_s)
CHEST, EAR = clean_ecg_signal(np.mean(X_adj[0:2, mask], axis=0), fs), clean_ecg_signal(np.mean(X_adj[2:4, mask], axis=0), fs)
t = t_all[mask] - t_all[mask][0]

print("Detecting peaks...")
rr_c, hr_c, Rl_c, tpk_c = qrs_filter(CHEST, fs)
rr_e, hr_e, Rl_e, tpk_e = qrs_filter(EAR, fs)

# Build Analysis Table
rows = []
for event in stim_intervals:
    for site, tpk in [("CHEST", tpk_c), ("EAR", tpk_e)]:
        m = (tpk >= event["start"]) & (tpk <= event["end"])
        if np.sum(m) > 5:
            hrv = compute_hrv(np.diff(tpk[m]))
            rows.append(dict(window=event["label"], site=site, **{k: hrv.get(k, np.nan) for k in HRV_UNITS}))

df_stim = pd.DataFrame(rows)
df_stim.to_csv(f"{outdir}/Stimulation_Analysis_Summary.csv", index=False)

# Quick Plot: Global Heart Rate
fig, ax = plt.subplots(figsize=(15, 5))
ax.plot(tpk_c[1:], hr_c, 'k', alpha=0.3, label="Chest HR")
for event in stim_intervals:
    ax.axvspan(event["start"], event["end"], color='red' if "On" in event["label"] else 'gray', alpha=0.15)
ax.set_title("Global Heart Rate & Stimulation Overlays")
save_plot(fig, f"{outdir}/Global_HR_Summary.png")

print(f"Analysis complete. Results saved to: {outdir}")




