"""
common.py: small helpers shared by several scripts.
"""

import geopandas as gpd
import numpy as np
import pandas as pd

import config

# NWI WETLAND_TYPE -> analysis group
NWI_GROUP = {
    "Freshwater Forested/Shrub Wetland": "vegetated_wetland",
    "Freshwater Emergent Wetland": "vegetated_wetland",
    "Estuarine and Marine Wetland": "vegetated_wetland",
    "Freshwater Pond": "pond",
    "Lake": "deepwater",
    "Riverine": "deepwater",
    "Estuarine and Marine Deepwater": "deepwater",
    "Other": "deepwater",
}


def canonical_site(org, site):
    """LDEQ lists the same site under LADEQWPD and LADEQWPD_WQX, sometimes with
    and sometimes without leading zeros (for example -43 and -0043); merge them."""
    if org in ("LADEQWPD", "LADEQWPD_WQX"):
        num = site.split("-", 1)[1]
        return "LDEQ-" + ((num.lstrip("0") or "0") if num.isdigit() else num)
    return site


def read_nwi(verbose=True):
    """Read the three NWI HU8 wetland shapefiles straight from their zip files."""
    parts = []
    for code in config.HUC8:
        z = config.DATA_RAW / f"NWI_HU8_{code}_Watershed.zip"
        path = f"zip://{z}!HU8_{code}_Watershed/HU8_{code}_Wetlands.shp"
        # The shapefile truncates the field name WETLAND_TYPE to WETLAND_TY
        g = gpd.read_file(path, columns=["ATTRIBUTE", "WETLAND_TY"])
        g = g.rename(columns={"WETLAND_TY": "WETLAND_TYPE"})
        g["source_hu8"] = code
        parts.append(g.to_crs(config.CRS_AREA))
        if verbose:
            print(f"  NWI HU8 {code}: {len(g):,} polygons")
    nwi = pd.concat(parts, ignore_index=True)
    nwi["geometry"] = nwi.geometry.make_valid()
    nwi["group"] = nwi["WETLAND_TYPE"].map(NWI_GROUP).fillna("deepwater")
    return gpd.GeoDataFrame(nwi, crs=config.CRS_AREA)


def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values (q-values); NaN stays NaN."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1)
    q[ok] = out
    return q
