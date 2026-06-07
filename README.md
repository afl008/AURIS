# AURIS — ECG/HRV Analysis Pipeline

Analysis pipeline for rodent cardiac data acquired via Tucker-Davis Technologies (TDT) hardware. Processes simultaneous chest and ear ECG recordings, detects R-peaks, and computes a comprehensive suite of HRV metrics across user-defined stimulation windows.

## What it does

- Loads a TDT block and extracts the EEG/ECG store
- Auto-corrects signal polarity across channels
- Applies a 60 Hz notch filter and 0.5 Hz high-pass to clean each channel
- Detects R-peaks with an adaptive MAD threshold + refractory-period filter
- Cleans the RR series (physiological band + outlier rejection), with the ear channel anchored to the simultaneous chest heart rate per window — see [RR cleaning & chest anchoring](#rr-cleaning--chest-anchoring)
- Computes **time-domain**, **frequency-domain**, and **nonlinear** HRV metrics per stimulation interval
- Saves results to CSV, a heart-rate overview plot, and a per-trial peak-detection diagnostic panel

### HRV metrics computed

| Domain | Metrics |
|---|---|
| Time | Mean RR, Mean HR, SDNN, RMSSD, pNN3, pNN5 |
| Frequency | LF, HF, Total Power, LF/HF ratio |
| Nonlinear | SD1, SD2, SD1/SD2, Sample Entropy, DFA α1, DFA α2, DFA α Ratio |

LF, HF, and Total Power are reported in ms².

## Setup

```bash
# Create and activate a virtual environment
python3.10 -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

> **TDT SDK note:** `tdt` requires a TDT system license. Install from [tdt.com](https://www.tdt.com/docs/sdk/offline-data-analysis/offline-data-python/) or via `pip install tdt`.

## Configuration

Session configuration lives in a **gitignored** `local_config.py` so your personal paths and intervals never appear in version history.

```bash
cp local_config.example.py local_config.py
# edit local_config.py with your paths and stim intervals
```

### Multi-trial mode (recommended)

Define a `trials` list and the script runs each block in sequence, clearing and recreating its output directory each time:

```python
trials = [
    {
        "label":          "FUS-VNS 2025-10-24",
        "block_path":     "/path/to/TDT/block_1",
        "start_time_s":   0,
        "end_time_s":     10000,
        "stim_intervals": [
            {"label": "Baseline", "start":   0, "end":  60},
            {"label": "Stim On",  "start":  60, "end": 120},
            {"label": "Recovery", "start": 120, "end": 180},
        ],
        # Optional per-trial overrides: "store", "channels", "outdir"
    },
    {
        "label":          "FUS-VNS 2025-11-03",
        "block_path":     "/path/to/TDT/block_2",
        "start_time_s":   0,
        "end_time_s":     10000,
        "stim_intervals": [],
    },
]
```

### Single-trial mode

Omit `trials` and set top-level variables instead:

```python
block_path     = "/path/to/your/TDT/block"
start_time_s   = 0
end_time_s     = 10000
stim_intervals = [
    {"label": "Baseline", "start":   0, "end":  60},
    {"label": "Stim On",  "start":  60, "end": 120},
]
# store    = "EEGw"
# channels = [1, 2, 3, 4]
# outdir   = f"{block_path}/my_custom_output"
```

If `local_config.py` is absent, the script runs with template defaults.

## Running

```bash
python "AURIS Data Analysis.py"
```

Outputs are written to `<block_path>/analysis_output_Rodent_I/`:

| File | Description |
|---|---|
| `Stimulation_Analysis_Summary.csv` | HRV metrics per interval × site |
| `Sequential_PrePost_Summary.csv` | HRV of the 5-min rest slices flanking each stim (before onset / after offset) × site — see [Sequential pre/post analysis](#sequential-prepost-analysis) |
| `Global_HR_Summary.png` | Heart rate trace with stimulation overlays |
| `Diagnostic_Peak_Detection.png` | Per-trial QC panel: raw + bandpass ECG with detected peaks (early-baseline and first-stim windows), RR timeseries, and RR histogram. Peaks used for HRV are marked ▼; peaks rejected by RR cleaning are marked ✗. |

## Methods

### RR cleaning & chest anchoring

Detected R-peaks are converted to RR intervals, then cleaned before any HRV metric is
computed (`clean_rr_for_hrv`):

1. **Physiological band** — intervals outside 100–350 ms (≈170–600 bpm for rodents) are dropped.
2. **Outlier rejection** — for the **chest** channel, intervals deviating more than ±20 % from
   the window's median are dropped (self-referenced).
3. **Chest anchoring (ear only)** — the **ear** channel is cleaned against the *simultaneous
   chest median for that same window* rather than its own median, with a ±30 % band.

Anchoring is necessary because the two ear electrodes are weakly coherent
(inter-electrode r ≈ 0 on some trials vs. r ≈ 0.5–0.85 for chest), so the ear detector
intermittently locks onto T-waves, producing doubled RR intervals. A self-referenced median
can land on the doubled rate and report half the true heart rate. Using the trusted chest
lead's per-window heart rate as a physiological prior rejects these artifacts. The anchor is
computed **per window** (not globally) so it tracks heart-rate drift across the recording.

### Sequential pre/post analysis

To assess the effect of a stimulation block, each stim is compared against the *rest
periods that flank it* rather than against the stim itself (`Sequential_PrePost_Summary.csv`):

- **Before** — up to 5 minutes of the preceding washout, ending at stim onset.
- **After** — up to 5 minutes of the following washout, starting at stim offset.

The paired contrast is `after − before` (post-stim rest vs. pre-stim rest), so both slices
are quiescent baselines and no during-stim artifact enters either side. Because each slice is
drawn from one end of a washout, the *after* slice of one stim and the *before* slice of the
next come from opposite ends of the same ≥10-min washout — so consecutive brackets share no
samples and the paired differences remain statistically independent. A stim contributes a pair
only when both flanking slices exist and span ≥2 minutes.

### Chest/ear sensor agreement

With per-window chest anchoring, agreement between the two sensors falls into three tiers:

| Tier | Metrics | Agreement |
|---|---|---|
| Rate | Mean RR, Mean HR | excellent (≈0–4 %) |
| Normalized / nonlinear | LF/HF, SD1/SD2, Sample Entropy, DFA α1/α2/Ratio | good (≈2–19 %) |
| Absolute variance | SDNN, RMSSD, LF, HF, Total Power, SD1, SD2 | ear runs ≈25–250 % higher |

The ear reliably recovers heart **rate** and **normalized** HRV dynamics. Absolute time-domain
variability is inflated on the ear by residual beat-detection jitter; this is a sensor
limitation, not a processing artifact — Savitzky–Golay smoothing rescales both channels
proportionally but does **not** close the gap, so it is intentionally **not** applied.

## Dependencies

| Package | Purpose |
|---|---|
| `numpy` | Numerical arrays |
| `pandas` | CSV output |
| `matplotlib` | Plotting |
| `scipy` | Filtering, spectral analysis, statistics |
| `scikit-learn` | KD-tree for sample entropy |
| `tdt` | TDT block I/O |

## Citation

If you use this pipeline, please cite the associated ASME paper (citation forthcoming).
