"""
09_maps.py: maps and summary figures.

  fig01_study_area.png          HUC8/HUC12 units, NHDPlus streams, WQ sites, gauges
  fig02_nlcd2021.png            NLCD 2021 land cover (standard NLCD colours)
  fig04_huc12_cover_2021.png    wetland, forest and developed % per HUC12
  fig05_huc12_change.png        2001 to 2021 change in those groups per HUC12
  fig06_study_area_groups.png   study-area group % in 2001 and 2021
  fig10_site_basins.png         NHDPlus upstream basins used to link sites to land cover
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import io

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from rasterio.features import rasterize

import config

NLCD_COLORS = {
    11: "#466b9f", 21: "#dec5c5", 22: "#d99282", 23: "#eb0000", 24: "#ab0000",
    31: "#b3ac9f", 41: "#68ab5f", 42: "#1c5f2c", 43: "#b5c58f", 52: "#ccb879",
    71: "#dfdfc2", 81: "#dcd939", 82: "#ab6c28", 90: "#b8d9eb", 95: "#6c9fb8",
}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
ORG_STYLE = {"LADEQWPD_WQX": ("#2a78d6", "o", "LDEQ ambient network"),
             "LDEQ319": ("#eb6834", "s", "LDEQ 319 nonpoint program"),
             "USGS-MS": ("#1baf7a", "^", "USGS Mississippi")}


def layers():
    L = {k: gpd.read_file(config.GPKG, layer=k).to_crs(config.CRS_AREA)
         for k in ("huc8", "huc12", "flowlines", "waterbodies", "wq_sites")}
    return L


def base(ax, L, h12_lw=0.3):
    L["waterbodies"].plot(ax=ax, color="#b8d4ee", linewidth=0)
    L["huc12"].boundary.plot(ax=ax, color="#9a9a96", linewidth=h12_lw)
    L["huc8"].boundary.plot(ax=ax, color=INK, linewidth=0.9)
    ax.set_axis_off()


def gauges():
    path = config.DATA_RAW / "nwis_gauge_sites.rdb"
    lines = [ln for ln in path.read_text().splitlines() if not ln.startswith("#")]
    df = pd.read_csv(io.StringIO("\n".join([lines[0]] + lines[2:])), sep="\t", dtype=str)
    g = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.dec_long_va.astype(float),
                                                         df.dec_lat_va.astype(float)), crs="EPSG:4269")
    return g.to_crs(config.CRS_AREA)


def fig01(L):
    fig, ax = plt.subplots(figsize=(7.5, 9))
    base(ax, L)
    f = L["flowlines"]
    for order in sorted(f["streamorde"].dropna().unique()):
        lw = {3: 0.35, 4: 0.6, 5: 0.9, 6: 1.3}.get(int(order), 1.7)
        f[f["streamorde"] == order].plot(ax=ax, color="#2a78d6", linewidth=lw)
    s = L["wq_sites"]
    s = s[s["used_in_medians"].astype(bool)]
    handles = []
    for org, (col, m, lab) in ORG_STYLE.items():
        sub = s[s["OrganizationIdentifier"] == org]
        if len(sub):
            sub.plot(ax=ax, color=col, marker=m, markersize=26, edgecolor="white", linewidth=0.6, zorder=4)
            handles.append(Line2D([], [], color=col, marker=m, linestyle="", markersize=6,
                                  markeredgecolor="white", label=f"{lab} (n={len(sub)})"))
    g = gauges()
    g.plot(ax=ax, color=INK, marker="*", markersize=110, edgecolor="white", linewidth=0.6, zorder=5)
    for _, r in g.iterrows():
        ax.annotate(f"USGS {r['site_no']}", (r.geometry.x, r.geometry.y), xytext=(6, 4),
                    textcoords="offset points", fontsize=7, color=INK)
    handles.append(Line2D([], [], color=INK, marker="*", linestyle="", markersize=10, label="USGS gauge"))
    handles.append(Line2D([], [], color="#2a78d6", linewidth=1.2, label="NHDPlus V2 stream, order 3+"))
    handles.append(Patch(facecolor="none", edgecolor=INK, label="HUC8 boundary"))
    handles.append(Patch(facecolor="none", edgecolor="#9a9a96", label="HUC12 boundary"))
    for _, r in L["huc8"].iterrows():
        c = r.geometry.representative_point()
        ax.text(c.x, c.y, r["name"], fontsize=9, color=INK2, ha="center", style="italic",
                bbox=dict(facecolor="white", alpha=0.6, edgecolor="none", pad=1))
    ax.legend(handles=handles, loc="lower left", fontsize=7, frameon=True)
    ax.set_title("Study area: Amite, Tickfaw and Lake Maurepas HUC8 units (LA and MS)", fontsize=10)
    fig.tight_layout()
    fig.savefig(config.FIG_DIR / "fig01_study_area.png", dpi=170)
    plt.close(fig)


def fig02(L):
    with rasterio.open(config.DATA_RAW / f"nlcd{config.NLCD_YEARS[-1]}_study_area.tif") as src:
        lc = src.read(1)
        inside = rasterize(((g, 1) for g in L["huc12"].geometry), out_shape=lc.shape,
                           transform=src.transform, fill=0, dtype="uint8")
        b = src.bounds
    lc = np.where(inside == 1, lc, 0)[::2, ::2]   # 60 m display resolution
    codes = sorted(NLCD_COLORS)
    lut = np.zeros(256, dtype=np.uint8)
    for i, c in enumerate(codes):
        lut[c] = i + 1
    idx = lut[lc]
    cmap = ListedColormap(["#ffffff"] + [NLCD_COLORS[c] for c in codes])
    fig, ax = plt.subplots(figsize=(9.5, 9))
    ax.imshow(np.ma.masked_equal(idx, 0), cmap=cmap, vmin=0, vmax=len(codes),
              extent=(b.left, b.right, b.bottom, b.top), interpolation="nearest")
    L["huc8"].boundary.plot(ax=ax, color=INK, linewidth=0.8)
    present = [c for c in codes if (lc == c).any()]
    ax.legend(handles=[Patch(color=NLCD_COLORS[c], label=config.NLCD_CLASSES[c]) for c in present],
              loc="lower left", bbox_to_anchor=(1.0, 0.0), fontsize=7, frameon=False, ncol=1)
    ax.set_axis_off()
    ax.set_title(f"NLCD {config.NLCD_YEARS[-1]} land cover", fontsize=10)
    fig.tight_layout()
    fig.savefig(config.FIG_DIR / "fig02_nlcd2021.png", dpi=170)
    plt.close(fig)


def choropleth_panels(L, df, cols, titles, fname, suptitle, diverging=False):
    h = L["huc12"].merge(df, left_on="huc12", right_index=True)
    fig, axes = plt.subplots(1, len(cols), figsize=(4.0 * len(cols), 5.6))
    for ax, col, t in zip(axes, cols, titles):
        if diverging:
            lim = max(abs(h[col].min()), abs(h[col].max()), 0.5)
            norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
            h.plot(ax=ax, column=col, cmap="PuOr_r", norm=norm, edgecolor="white", linewidth=0.3)
            sm = plt.cm.ScalarMappable(norm=norm, cmap="PuOr_r")
        else:
            bounds = [0, 5, 10, 20, 30, 40, 50, 60, 80, 100]
            cmap = plt.get_cmap({"wetland": "Blues", "forest": "Greens", "developed": "Reds"}
                                .get(next((k for k in ("wetland", "forest", "developed") if k in col), ""), "Greys"))
            norm = BoundaryNorm(bounds, cmap.N)
            h.plot(ax=ax, column=col, cmap=cmap, norm=norm, edgecolor="white", linewidth=0.3)
            sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        L["huc8"].boundary.plot(ax=ax, color=INK, linewidth=0.7)
        cb = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.045, pad=0.02)
        cb.ax.tick_params(labelsize=7)
        cb.set_label("percentage points" if diverging else "% of HUC12 area", fontsize=7.5)
        ax.set_title(t, fontsize=9)
        ax.set_axis_off()
    fig.suptitle(suptitle, fontsize=10)
    fig.tight_layout()
    fig.savefig(config.FIG_DIR / fname, dpi=170)
    plt.close(fig)


def fig06():
    t = pd.read_csv(config.TAB_DIR / "nlcd_class_area_study_area.csv")
    inv = {c: g for g, cs in config.NLCD_GROUPS.items() for c in cs}
    t["group"] = t["nlcd_code"].map(inv)
    g = t.groupby(["group", "year"])["percent"].sum().unstack()
    order = ["wetland", "forest", "developed", "agriculture", "shrub_grass", "water", "barren"]
    g = g.reindex([o for o in order if o in g.index])
    y0, y1 = config.NLCD_YEARS
    fig, ax = plt.subplots(figsize=(7.5, 3.8))
    x = np.arange(len(g))
    w = 0.38
    ax.bar(x - w / 2 - 0.01, g[y0], w, color="#8fb8e8", label=str(y0))
    ax.bar(x + w / 2 + 0.01, g[y1], w, color="#2a78d6", label=str(y1))
    for i, (a, b) in enumerate(zip(g[y0], g[y1])):
        ax.text(i + w / 2 + 0.01, b + 0.6, f"{b - a:+.1f}", ha="center", fontsize=7.5, color=INK2)
    ax.set_xticks(x, [s.replace("_", "/") for s in g.index], fontsize=8)
    ax.set_ylabel("% of study area", fontsize=8)
    ax.tick_params(labelsize=7.5)
    ax.grid(axis="y", color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=8, frameon=False)
    ax.set_title(f"NLCD land-cover groups, {y0} and {y1} (labels: change in percentage points)", fontsize=9)
    fig.tight_layout()
    fig.savefig(config.FIG_DIR / "fig06_study_area_groups.png", dpi=170)
    plt.close(fig)


def fig10(L):
    b = gpd.read_file(config.GPKG, layer="site_basins").to_crs(config.CRS_AREA)
    b = b.sort_values("basin_area_km2", ascending=False)
    fig, ax = plt.subplots(figsize=(7.5, 9))
    base(ax, L, h12_lw=0.2)
    nb = b[b["link_method"] == "nhdplus_basin"]
    nb.boundary.plot(ax=ax, color="#2a78d6", linewidth=0.6, alpha=0.8)
    s = L["wq_sites"]
    s = s[s["used_in_medians"].astype(bool)]
    s.plot(ax=ax, color="#eb6834", markersize=14, edgecolor="white", linewidth=0.5, zorder=4)
    ax.legend(handles=[Line2D([], [], color="#2a78d6", label=f"NHDPlus upstream basin (n={len(nb)})"),
                       Line2D([], [], color="#eb6834", marker="o", linestyle="", label="Water-quality site")],
              loc="lower left", fontsize=7.5)
    ax.set_title("Upstream basins used to link sites to land cover", fontsize=10)
    fig.tight_layout()
    fig.savefig(config.FIG_DIR / "fig10_site_basins.png", dpi=170)
    plt.close(fig)


def main():
    L = layers()
    fig01(L)
    fig02(L)
    loc = pd.read_csv(config.TAB_DIR / "nlcd_groups_huc12_local.csv", dtype={"huc12": str}).set_index("huc12")
    y0, y1 = config.NLCD_YEARS
    choropleth_panels(L, loc, [f"pct_wetland_{y1}", f"pct_forest_{y1}", f"pct_developed_{y1}"],
                      ["Wetland (NLCD 90, 95)", "Forest (NLCD 41-43)", "Developed (NLCD 21-24)"],
                      "fig04_huc12_cover_2021.png", f"HUC12 land cover, NLCD {y1}")
    choropleth_panels(L, loc, ["chg_wetland_pts", "chg_forest_pts", "chg_developed_pts"],
                      ["Wetland", "Forest", "Developed"],
                      "fig05_huc12_change.png", f"Change in HUC12 cover, NLCD {y0} to {y1}",
                      diverging=True)
    fig06()
    fig10(L)
    print("Figures written to outputs/figures/")


if __name__ == "__main__":
    main()
