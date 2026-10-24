"""
04_wetlands_nwi.py: USFWS National Wetlands Inventory (NWI) wetland area per
HUC12, local and upstream.

NWI polygons are read directly from the HU8 zip files, projected to EPSG:5070
and intersected with the HUC12 polygons. NWI WETLAND_TYPE values are grouped:

  vegetated wetland : Freshwater Forested/Shrub Wetland,
                      Freshwater Emergent Wetland, Estuarine and Marine Wetland
  pond              : Freshwater Pond
  deepwater         : Lake, Riverine, Estuarine and Marine Deepwater, Other

Outputs
  outputs/tables/nwi_type_area_study_area.csv
  outputs/tables/nwi_huc12_local.csv
  outputs/tables/nwi_huc12_upstream.csv
  outputs/figures/fig03_nwi_wetlands.png
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch

import config
from common import read_nwi

TYPE_COLORS = {
    "Freshwater Forested/Shrub Wetland": "#008300",
    "Freshwater Emergent Wetland": "#1baf7a",
    "Estuarine and Marine Wetland": "#4a3aa7",
    "Freshwater Pond": "#2a78d6",
    "Lake": "#8fb8e8",
    "Riverine": "#8fb8e8",
    "Estuarine and Marine Deepwater": "#8fb8e8",
    "Other": "#8fb8e8",
}


def main():
    h12 = gpd.read_file(config.GPKG, layer="huc12").to_crs(config.CRS_AREA)
    print("Reading NWI")
    nwi = read_nwi()
    print("Intersecting NWI with HUC12 polygons")
    inter = gpd.overlay(nwi[["WETLAND_TYPE", "group", "geometry"]],
                        h12[["huc12", "geometry"]], how="intersection", keep_geom_type=True)
    inter["area_km2"] = inter.area / 1e6

    by_type = inter.groupby("WETLAND_TYPE")["area_km2"].sum().sort_values(ascending=False)
    study_km2 = h12.area.sum() / 1e6
    t = by_type.rename("area_km2").to_frame()
    t["pct_of_study_area"] = 100 * t["area_km2"] / study_km2
    t.to_csv(config.TAB_DIR / "nwi_type_area_study_area.csv", float_format="%.3f")

    piv = inter.pivot_table(index="huc12", columns="group", values="area_km2",
                            aggfunc="sum", fill_value=0)
    for col in ("vegetated_wetland", "pond", "deepwater"):
        if col not in piv:
            piv[col] = 0.0
    area = h12.set_index("huc12").area / 1e6
    piv = piv.reindex(area.index, fill_value=0)
    local = pd.DataFrame({
        "area_km2": area,
        "nwi_vegwet_km2": piv["vegetated_wetland"],
        "pct_nwi_vegwet": 100 * piv["vegetated_wetland"] / area,
        "pct_nwi_pond": 100 * piv["pond"] / area,
        "pct_nwi_deepwater": 100 * piv["deepwater"] / area,
    })
    local.to_csv(config.TAB_DIR / "nwi_huc12_local.csv", float_format="%.3f")

    members = pd.read_csv(config.DATA_DIR / "huc12_upstream_members.csv", dtype=str)
    m = members.merge(local[["area_km2", "nwi_vegwet_km2"]], left_on="upstream_huc12",
                      right_index=True).groupby("huc12")[["area_km2", "nwi_vegwet_km2"]].sum()
    m["pct_nwi_vegwet_upstream"] = 100 * m["nwi_vegwet_km2"] / m["area_km2"]
    m.to_csv(config.TAB_DIR / "nwi_huc12_upstream.csv", float_format="%.3f")

    print("\nNWI area by type inside the study area")
    print(t.round(2).to_string())
    veg = piv["vegetated_wetland"].sum()
    print(f"Vegetated NWI wetland: {veg:,.1f} km2 = {100 * veg / study_km2:.2f}% of {study_km2:,.1f} km2")
    print(f"HUC12 vegetated wetland %: median {local['pct_nwi_vegwet'].median():.1f}, "
          f"range {local['pct_nwi_vegwet'].min():.1f} to {local['pct_nwi_vegwet'].max():.1f}")

    # Map
    fig, ax = plt.subplots(figsize=(7.5, 8.5))
    h12.boundary.plot(ax=ax, color="#9a9a96", linewidth=0.3)
    draw = inter.copy()
    draw["geometry"] = draw.geometry.simplify(15)
    for wt, col in TYPE_COLORS.items():
        sub = draw[draw["WETLAND_TYPE"] == wt]
        if len(sub):
            sub.plot(ax=ax, color=col, linewidth=0)
    gpd.read_file(config.GPKG, layer="huc8").to_crs(config.CRS_AREA).boundary.plot(
        ax=ax, color="#0b0b0b", linewidth=0.8)
    handles = [Patch(color=c, label=k) for k, c in TYPE_COLORS.items()
               if k in ("Freshwater Forested/Shrub Wetland", "Freshwater Emergent Wetland",
                        "Estuarine and Marine Wetland", "Freshwater Pond")]
    handles.append(Patch(color="#8fb8e8", label="Lake, riverine and other deepwater"))
    ax.legend(handles=handles, loc="lower left", fontsize=7.5, frameon=True)
    ax.set_title("USFWS National Wetlands Inventory, Amite, Tickfaw and Lake Maurepas HUC8 units",
                 fontsize=9.5)
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(config.FIG_DIR / "fig03_nwi_wetlands.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
