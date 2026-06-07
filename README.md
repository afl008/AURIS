# AURIS — ECG/HRV Analysis Pipeline

Analysis pipeline for rodent cardiac data acquired via Tucker-Davis Technologies (TDT) hardware. Processes simultaneous chest and ear ECG recordings, detects R-peaks, and computes a comprehensive suite of HRV metrics across user-defined stimulation windows.

## What it does

- Loads a TDT block and extracts the EEG/ECG store
- Auto-corrects signal polarity across channels
- Applies a 60 Hz notch filter and 0.5 Hz high-pass to clean each channel
- Detects R-peaks with an adaptive MAD threshold + refractory-period filter
- Computes **time-domain**, **frequency-domain**, and **nonlinear** HRV metrics per stimulation interval
- Saves results to CSV and a heart rate overview plot

### HRV metrics computed

| Domain | Metrics |
|---|---|
| Time | Mean RR, Mean HR, SDNN, RMSSD, pNN3, pNN5 |
| Frequency | LF, HF, Total Power, LF/HF ratio |
| Nonlinear | SD1, SD2, SD1/SD2, Sample Entropy |

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

Open `AURIS Data Analysis.py` and edit the **User Configuration** block at the top:

```python
block_path   = "/path/to/your/TDT/block"   # path to the .Tbk block folder
start_time_s = 0                            # analysis window start (seconds)
end_time_s   = 10000                        # analysis window end (seconds)

stim_intervals = [
    {"label": "Baseline",   "start":   0, "end":  60},
    {"label": "Stim On",    "start":  60, "end": 120},
    {"label": "Recovery",   "start": 120, "end": 180},
]

store    = "EEGw"       # TDT store name
channels = [1, 2, 3, 4] # Chest = 1–2, Ear = 3–4
```

## Running

```bash
python "AURIS Data Analysis.py"
```

Outputs are written to `<block_path>/analysis_output_Rodent_I/`:

| File | Description |
|---|---|
| `Stimulation_Analysis_Summary.csv` | HRV metrics per interval × site |
| `Global_HR_Summary.png` | Heart rate trace with stimulation overlays |

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
