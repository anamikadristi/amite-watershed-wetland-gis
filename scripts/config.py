"""
config.py: shared settings for every numbered script.

Folder paths, coordinate reference systems, study-area codes and analysis
thresholds live here so they are defined once.

Large downloads go to DATA_RAW, which is kept out of Git. By default this is
<repo>/data_raw/. Set the environment variable AMITE_DATA_RAW to use a
different folder, for example a scratch disk.
"""

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA_RAW = Path(os.environ.get("AMITE_DATA_RAW", ROOT / "data_raw"))
DATA_DIR = ROOT / "data"            # small derived public files kept in Git
OUT_DIR = ROOT / "outputs"
FIG_DIR = OUT_DIR / "figures"
TAB_DIR = OUT_DIR / "tables"
GPKG = DATA_DIR / "amite_maurepas_layers.gpkg"

for d in (DATA_RAW, DATA_DIR, FIG_DIR, TAB_DIR):
    d.mkdir(parents=True, exist_ok=True)

# Coordinate reference systems
CRS_GEO = "EPSG:4326"     # WGS84 lat/lon, used for web service queries
CRS_AREA = "EPSG:5070"    # NAD83 CONUS Albers equal area, the NLCD grid

# Study area: three HUC8 units that drain to Lake Maurepas
HUC8 = {
    "08070202": "Amite",
    "08070203": "Tickfaw",
    "08070204": "Lake Maurepas",
}

NLCD_YEARS = [2001, 2021]

# Water Quality Portal query window
WQ_START = "01-01-2000"   # MM-DD-YYYY, the format WQP expects
WQ_END = "12-31-2025"

# USGS gauges used for discharge context
GAUGES = {
    "07378500": "Amite River near Denham Springs, LA",
    "07377500": "Comite River near Comite, LA",
    "07380120": "Amite River at Port Vincent, LA",
}

# A site needs at least this many samples of a parameter to get a median
MIN_SAMPLES_PER_SITE = 5
# A parameter needs at least this many sites before a correlation is run
MIN_SITES_FOR_TEST = 6

# NLCD class codes and the groups used in this project
NLCD_CLASSES = {
    11: "Open water", 12: "Perennial ice/snow",
    21: "Developed, open space", 22: "Developed, low intensity",
    23: "Developed, medium intensity", 24: "Developed, high intensity",
    31: "Barren land", 41: "Deciduous forest", 42: "Evergreen forest",
    43: "Mixed forest", 52: "Shrub/scrub", 71: "Grassland/herbaceous",
    81: "Pasture/hay", 82: "Cultivated crops", 90: "Woody wetlands",
    95: "Emergent herbaceous wetlands",
}
NLCD_GROUPS = {
    "water": [11, 12],
    "developed": [21, 22, 23, 24],
    "barren": [31],
    "forest": [41, 42, 43],
    "shrub_grass": [52, 71],
    "agriculture": [81, 82],
    "wetland": [90, 95],
}
