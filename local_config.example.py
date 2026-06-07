# Copy this file to local_config.py and fill in your session details.
# local_config.py is gitignored — it will never be committed.

block_path   = "/path/to/your/TDT/block"
start_time_s = 0
end_time_s   = 10000

stim_intervals = [
    {"label": "Baseline", "start":   0, "end":  60},
    {"label": "Stim On",  "start":  60, "end": 120},
    {"label": "Recovery", "start": 120, "end": 180},
]

# Optional overrides (defaults shown — omit lines you don't need to change)
# store    = "EEGw"
# channels = [1, 2, 3, 4]   # Chest=1–2, Ear=3–4
# outdir   = f"{block_path}/analysis_output_Rodent_I"
