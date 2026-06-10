# AURIS — ECG/HRV Analysis Pipeline

Analysis pipeline for rodent cardiac data acquired via Tucker-Davis Technologies (TDT) hardware. It processes simultaneous chest and ear ECG recordings, detects R-peaks, and computes HRV metrics across user-defined stimulation windows.

## What it does

- Loads a TDT block and extracts the ECG store
- Auto-corrects signal polarity across channels
- Applies a notch filter and high-pass filter to clean each channel
- Detects R-peaks with an adaptive MAD threshold and refractory-period filter
- Cleans RR intervals using physiological band filtering and outlier rejection, with the ear channel anchored to the simultaneous chest heart rate per window
- Computes **time-domain**, **frequency-domain**, and **nonlinear** HRV metrics per stimulation interval
- Saves results to CSV along with a heart-rate overview plot and per-trial peak-detection diagnostics

### HRV metrics computed

| Domain | Metrics |
|---|---|
| Time | Mean RR, Mean HR, SDNN, RMSSD, pNN3, pNN5 |
| Frequency | LF, HF, Total Power, LF/HF ratio |
| Nonlinear | SD1, SD2, SD1/SD2, Sample Entropy, DFA α1, DFA α2, DFA α Ratio |

## Setup

```bash
python3.10 -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

> **TDT SDK note:** `tdt` requires a TDT system license. Install from [tdt.com](https://www.tdt.com/docs/sdk/offline-data-analysis/offline-data-python/) or via `pip install tdt`.

## Configuration

Session configuration lives in a **gitignored** `local_config.py` so paths and intervals stay out of version history.

```bash
cp local_config.example.py local_config.py
# edit local_config.py with your paths and stim intervals
```

### Multi-trial mode (recommended)

Define a `trials` list to run multiple blocks in sequence:

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
| `Sequential_PrePost_Summary.csv` | HRV of rest slices flanking each stim (before onset / after offset) × site |
| `Global_HR_Summary.png` | Heart rate trace with stimulation overlays |
| `Diagnostic_Peak_Detection.png` | Per-trial QC panel with detected peaks, RR timeseries, and RR histogram |

## Methods

### RR cleaning & chest anchoring

After R-peak detection, RR intervals are cleaned before HRV metrics are computed:

1. **Physiological band** — intervals outside the expected rodent heart rate range are dropped.
2. **Outlier rejection** — intervals deviating too far from the window median are removed.
3. **Chest anchoring (ear only)** — the ear channel is cleaned against the simultaneous chest median for each window rather than its own, accounting for the weaker signal coherence of the ear electrodes.

Anchoring is computed per window so it tracks heart-rate drift across the recording.

### Sequential pre/post analysis

Each stimulation is compared against the rest periods flanking it rather than the stim itself. Up to 5 minutes of preceding washout ("before") and following washout ("after") are extracted and contrasted as paired baselines. This keeps any during-stim artifact out of the comparison, and consecutive brackets share no overlapping samples.

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

If you use this pipeline, please cite:

> Bohluli R. et al., "Monitoring Autonomic Tone During Spinal Cord Neuromodulation Using Wearable AURIS Sensor." https://doi.org/10.64898/2026.03.07.709943
