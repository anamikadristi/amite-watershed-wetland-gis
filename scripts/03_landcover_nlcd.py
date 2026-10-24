"""
03_landcover_nlcd.py: NLCD 2001 and 2021 land cover per HUC12, upstream
totals, and 2001 to 2021 change.

Method
  The HUC12 polygons are rasterized onto the NLCD 30 m grid (EPSG:5070).
  A cell belongs to a HUC12 when its centre falls inside the polygon. Cells
  are then counted by HUC12 and NLCD class with one numpy bincount, which is
  exact and fast. Percentages use all classified cells (water included).

Outputs (outputs/tables/)
  nlcd_class_area_study_area.csv     class areas for both years
  nlcd_groups_huc12_local.csv        group % per HUC12, both years, and change
  nlcd_groups_huc12_upstream.csv     group % for the whole area upstream of each HUC12
  nlcd_transition_2001_2021.csv      group-to-group transition matrix (km2)
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize

import config

CELL_KM2 = 30 * 30 / 1e6


def zone_raster(h12, src):
    """Rasterize HUC12 index (1..n) onto the grid of an open rasterio dataset."""
    shapes = ((g, i + 1) for i, g in enumerate(h12.to_crs(src.crs).geometry))
    return rasterize(shapes, out_shape=(src.height, src.width), transform=src.transform,
                     fill=0, dtype="int32")


def counts_by_zone(zones, lc, n_zones):
    """Matrix [zone, class] of cell counts (classes 0..255)."""
    idx = zones.astype(np.int64) * 256 + lc.astype(np.int64)
    c = np.bincount(idx.ravel(), minlength=(n_zones + 1) * 256)
    return c.reshape(n_zones + 1, 256)


def group_pct(counts_df):
    """Convert a [unit x class] count table to % by NLCD group."""
    out = pd.DataFrame(index=counts_df.index)
    total = counts_df.sum(axis=1)
    for g, codes in config.NLCD_GROUPS.items():
        cols = [c for c in codes if c in counts_df.columns]
        out[f"pct_{g}"] = 100 * counts_df[cols].sum(axis=1) / total
    return out


def main():
    h12 = gpd.read_file(config.GPKG, layer="huc12").sort_values("huc12").reset_index(drop=True)
    n = len(h12)
    codes = sorted(config.NLCD_CLASSES)
    per_year = {}
    grids = {}
    for year in config.NLCD_YEARS:
        with rasterio.open(config.DATA_RAW / f"nlcd{year}_study_area.tif") as src:
            lc = src.read(1)
            zones = zone_raster(h12, src)
        grids[year] = lc
        m = counts_by_zone(zones, lc, n)[1:, :]          # drop zone 0 (outside)
        df = pd.DataFrame(m[:, codes], columns=codes, index=h12["huc12"])
        per_year[year] = df
        print(f"NLCD {year}: {int(df.values.sum()):,} cells inside HUC12 units")

    # Study-area class table
    rows = []
    for year, df in per_year.items():
        tot = df.values.sum()
        for c in codes:
            rows.append({"year": year, "nlcd_code": c, "class_name": config.NLCD_CLASSES[c],
                         "area_km2": df[c].sum() * CELL_KM2,
                         "percent": 100 * df[c].sum() / tot})
    cls = pd.DataFrame(rows)
    cls = cls[cls.groupby("nlcd_code")["area_km2"].transform("sum") > 0]
    cls.to_csv(config.TAB_DIR / "nlcd_class_area_study_area.csv", index=False, float_format="%.3f")

    # Local (per HUC12) group percentages and change
    y0, y1 = config.NLCD_YEARS
    loc0 = group_pct(per_year[y0]).add_suffix(f"_{y0}")
    loc1 = group_pct(per_year[y1]).add_suffix(f"_{y1}")
    local = loc0.join(loc1)
    for g in config.NLCD_GROUPS:
        local[f"chg_{g}_pts"] = local[f"pct_{g}_{y1}"] - local[f"pct_{g}_{y0}"]
    local = h12.set_index("huc12")[["name", "huc8_name", "area_km2"]].join(local)
    local.to_csv(config.TAB_DIR / "nlcd_groups_huc12_local.csv", float_format="%.3f")

    # Upstream totals: sum counts over every HUC12 that drains through a unit
    members = pd.read_csv(config.DATA_DIR / "huc12_upstream_members.csv", dtype=str)
    up_tabs = []
    for year, df in per_year.items():
        up = members.merge(df, left_on="upstream_huc12", right_index=True)
        up = up.groupby("huc12")[codes].sum()
        up_tabs.append(group_pct(up).add_suffix(f"_{year}"))
    upstream = h12.set_index("huc12")[["name", "upstream_area_km2"]].join(pd.concat(up_tabs, axis=1))
    upstream.to_csv(config.TAB_DIR / "nlcd_groups_huc12_upstream.csv", float_format="%.3f")

    # Transition matrix 2001 -> 2021 by group, cells inside the study area only
    with rasterio.open(config.DATA_RAW / f"nlcd{y1}_study_area.tif") as src:
        inside = zone_raster(h12, src) > 0
    lut = np.zeros(256, dtype=np.int16) - 1
    gnames = list(config.NLCD_GROUPS)
    for gi, g in enumerate(gnames):
        lut[config.NLCD_GROUPS[g]] = gi
    a = lut[grids[y0][inside]]
    b = lut[grids[y1][inside]]
    ok = (a >= 0) & (b >= 0)
    k = len(gnames)
    tm = np.bincount(a[ok] * k + b[ok], minlength=k * k).reshape(k, k) * CELL_KM2
    tmat = pd.DataFrame(tm, index=[f"{g}_{y0}" for g in gnames], columns=[f"{g}_{y1}" for g in gnames])
    tmat.to_csv(config.TAB_DIR / "nlcd_transition_2001_2021.csv", float_format="%.3f")
    changed = 100 * (tm.sum() - np.trace(tm)) / tm.sum()

    # Console summary
    tot = {y: per_year[y].values.sum() for y in per_year}
    print("\nStudy-area group percentages")
    for g, cds in config.NLCD_GROUPS.items():
        p0 = 100 * per_year[y0][[c for c in cds if c in codes]].values.sum() / tot[y0]
        p1 = 100 * per_year[y1][[c for c in cds if c in codes]].values.sum() / tot[y1]
        print(f"  {g:12s} {y0}: {p0:6.2f}%  {y1}: {p1:6.2f}%  change {p1 - p0:+.2f} pts")
    print(f"Cells that changed group between {y0} and {y1}: {changed:.2f}%")
    print(f"HUC12 with largest developed gain: "
          f"{local['chg_developed_pts'].idxmax()} ({local['chg_developed_pts'].max():+.2f} pts)")


if __name__ == "__main__":
    main()
