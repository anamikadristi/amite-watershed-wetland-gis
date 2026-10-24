"""
02_prepare_boundaries.py: build the HUC12 study-area layer and its
upstream network.

Steps
  1. Read WBD HUC12 and HUC8 polygons and NHDPlus V2 flowlines/waterbodies.
  2. Compute HUC12 area in an equal-area projection (EPSG:5070).
  3. Use the WBD "tohuc" field to find, for every HUC12, all HUC12 units
     upstream of it (including itself). Later scripts sum land cover over
     this set to describe the area that drains to a sampling site.

Outputs
  data/amite_maurepas_layers.gpkg   layers: huc8, huc12, flowlines, waterbodies
  data/huc12_upstream_members.csv   long table: huc12, upstream_huc12
  outputs/tables/huc12_attributes.csv
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import geopandas as gpd
import pandas as pd

import config


def upstream_sets(huc12, tohuc):
    """Return {huc: set of HUC12s that drain through huc, itself included}."""
    down = dict(zip(huc12, tohuc))
    members = {h: {h} for h in huc12}
    for h in huc12:
        seen = set()
        nxt = down.get(h)
        while nxt in members and nxt not in seen:  # guard against loops
            members[nxt].add(h)
            seen.add(nxt)
            nxt = down.get(nxt)
    return members


def main():
    raw = config.DATA_RAW
    h12 = gpd.read_file(raw / "wbd_huc12.geojson")
    h8 = gpd.read_file(raw / "wbd_huc8.geojson")
    h12["huc8"] = h12["huc12"].str[:8]
    h12["huc8_name"] = h12["huc8"].map(config.HUC8)
    h12["area_km2"] = h12.to_crs(config.CRS_AREA).area / 1e6
    h8["area_km2"] = h8.to_crs(config.CRS_AREA).area / 1e6

    members = upstream_sets(h12["huc12"].tolist(), h12["tohuc"].tolist())
    area = dict(zip(h12["huc12"], h12["area_km2"]))
    h12["n_upstream_huc12"] = h12["huc12"].map(lambda h: len(members[h]))
    h12["upstream_area_km2"] = h12["huc12"].map(lambda h: sum(area[m] for m in members[h]))

    long = pd.DataFrame(
        [(h, m) for h, ms in members.items() for m in sorted(ms)],
        columns=["huc12", "upstream_huc12"],
    )
    long.to_csv(config.DATA_DIR / "huc12_upstream_members.csv", index=False)

    flw = gpd.read_file(raw / "nhdplus_flowlines_order3plus.geojson")
    wb = gpd.read_file(raw / "nhdplus_waterbodies_1km2.geojson")
    study = h12.dissolve()[["geometry"]]
    flw = gpd.clip(flw.to_crs(study.crs), study)
    wb = gpd.clip(wb.to_crs(study.crs), study)
    flw["geometry"] = flw.geometry.simplify(0.0002)   # about 20 m, for mapping only
    wb["geometry"] = wb.geometry.simplify(0.0002)

    gpkg = config.GPKG
    if gpkg.exists():
        gpkg.unlink()
    h8.to_file(gpkg, layer="huc8", driver="GPKG")
    h12.to_file(gpkg, layer="huc12", driver="GPKG")
    flw.to_file(gpkg, layer="flowlines", driver="GPKG")
    wb.to_file(gpkg, layer="waterbodies", driver="GPKG")

    cols = ["huc12", "name", "huc8", "huc8_name", "tohuc", "area_km2",
            "n_upstream_huc12", "upstream_area_km2"]
    h12[cols].sort_values("huc12").to_csv(
        config.TAB_DIR / "huc12_attributes.csv", index=False, float_format="%.2f")

    print(f"HUC12 units: {len(h12)}")
    print(h12.groupby("huc8_name")["area_km2"].agg(["count", "sum"]).round(1).to_string())
    print(f"Study area: {h12['area_km2'].sum():,.1f} km2")
    print(f"Largest upstream area: {h12['upstream_area_km2'].max():,.1f} km2")
    print(f"Flowline segments (order 3+): {len(flw)}, waterbodies >= 1 km2: {len(wb)}")


if __name__ == "__main__":
    main()
