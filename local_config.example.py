# Copy this file to local_config.py and fill in your session details.
# local_config.py is gitignored — it will never be committed.
#
# MULTI-TRIAL MODE: define a `trials` list and the script will run each in sequence,
# clearing and recreating the output directory for each trial.
#
# SINGLE-TRIAL MODE: omit `trials` and set the top-level variables instead.

# --- Multi-trial (recommended) ---
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
        # Optional per-trial overrides:
        # "store":    "EEGw",
        # "channels": [1, 2, 3, 4],
        # "outdir":   "/custom/output/path",
    },
    {
        "label":          "FUS-VNS 2025-11-03",
        "block_path":     "/path/to/TDT/block_2",
        "start_time_s":   0,
        "end_time_s":     10000,
        "stim_intervals": [],
    },
]

# --- Single-trial (alternative) ---
# block_path     = "/path/to/your/TDT/block"
# start_time_s   = 0
# end_time_s     = 10000
# stim_intervals = [
#     {"label": "Baseline", "start":   0, "end":  60},
#     {"label": "Stim On",  "start":  60, "end": 120},
# ]
# store    = "EEGw"
# channels = [1, 2, 3, 4]
# outdir   = f"{block_path}/my_custom_output"
