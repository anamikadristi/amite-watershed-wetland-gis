"""
run_all.py: run the analysis scripts in order (02 to 09).

Run 01_download_data.py first (once). Then:
    python scripts/run_all.py
Add --download to run the download step as well.
"""

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEPS = [
    "02_prepare_boundaries.py",
    "03_landcover_nlcd.py",
    "04_wetlands_nwi.py",
    "05_water_quality.py",
    "06_site_basins.py",
    "07_relationships.py",
    "08_discharge_context.py",
    "09_maps.py",
]

if "--download" in sys.argv:
    STEPS.insert(0, "01_download_data.py")

for step in STEPS:
    print(f"\n===== {step} =====", flush=True)
    subprocess.run([sys.executable, str(HERE / step)], check=True)
print("\nAll steps finished.")
