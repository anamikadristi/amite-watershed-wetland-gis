"""
06_site_basins.py: land cover and wetland cover in the area that drains to
each water-quality site.

Point-to-watershed linking
  1. Preferred: a local NHDPlus V2 delineation. The site is placed in the
     NHDPlus catchment that contains it (or the nearest catchment within
     SNAP_M metres). The network is then walked upstream with the
     fromnode/tonode fields and all upstream catchments are merged. This is
     the same logic the USGS NLDI basin service uses (without splitting the
     local catchment at the site), and a few NLDI basins are compared with
     it as a check.
  2. Fallback: if no catchment is found, the site's own HUC12 polygon
     (local HUC12 composition). These sites are flagged in link_method and
     a sensitivity run in 07_relationships.py drops them.

For each basin we compute NLCD 2001 and 2021 group percentages, NWI
vegetated wetland %, and the share of the basin covered by the downloaded
NLCD and NWI extents. Sites with NLCD coverage below 95% are flagged and
NWI % is left blank where NWI covers less than 95% of the basin.

Outputs
  outputs/tables/site_basin_landcover.csv
  data/amite_maurepas_layers.gpkg  layer: site_basins (simplified)
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask

import config
from common import canonical_site, read_nwi

CELL_KM2 = 30 * 30 / 1e6


SNAP_M = 300


def read_pages(stem):
    files = sorted(config.DATA_RAW.glob(f"{stem}_p*.geojson"))
    return pd.concat([gpd.read_file(f) for f in files], ignore_index=True)


def delineate(sites):
    """Upstream NHDPlus basin for each site; returns GeoDataFrame in EPSG:5070."""
    net = read_pages("nhdplus_network").drop_duplicates("comid")
    cat = read_pages("nhdplus_catchments").drop_duplicates("featureid")
    cat = gpd.GeoDataFrame(cat, geometry="geometry").to_crs(config.CRS_AREA)
    cat["featureid"] = cat["featureid"].astype(int)
    net["comid"] = net["comid"].astype(int)
    print(f"NHDPlus flowlines: {len(net):,}; catchments: {len(cat):,}")

    # upstream adjacency: flowline u is upstream of d when u.tonode == d.fromnode
    by_tonode = net.groupby("tonode")["comid"].apply(list).to_dict()
    fromnode = dict(zip(net["comid"], net["fromnode"]))
    totda = dict(zip(net["comid"], net["totdasqkm"]))
    geom = dict(zip(cat["featureid"], cat.geometry))

    def upstream(c):
        seen, stack = set(), [c]
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            stack.extend(by_tonode.get(fromnode.get(x), []))
        return seen

    pts = gpd.GeoDataFrame(sites[["site_id"]].copy(),
                           geometry=gpd.points_from_xy(sites.lon, sites.lat), crs=config.CRS_GEO
                           ).to_crs(config.CRS_AREA)
    hit = gpd.sjoin_nearest(pts, cat[["featureid", "geometry"]], how="left",
                            max_distance=SNAP_M, distance_col="snap_m")
    hit = hit.drop_duplicates("site_id")
    rows = []
    for _, r in hit.iterrows():
        if pd.isna(r["featureid"]) or int(r["featureid"]) not in fromnode:
            continue
        c = int(r["featureid"])
        ups = upstream(c)
        polys = [geom[u] for u in ups if u in geom]
        rows.append({"site_id": r["site_id"], "comid": c, "snap_m": r["snap_m"],
                     "n_catchments": len(polys), "nhdplus_totdasqkm": totda.get(c, np.nan),
                     "geometry": gpd.GeoSeries(polys, crs=config.CRS_AREA).union_all()})
    return gpd.GeoDataFrame(rows, geometry="geometry", crs=config.CRS_AREA)


def nldi_check(basins):
    """Compare local basins with NLDI basins downloaded with --nldi, if any."""
    folder = config.DATA_RAW / "nldi_basins"
    rows = []
    for f in sorted(folder.glob("*.geojson")) if folder.exists() else []:
        g = gpd.read_file(f)
        if g.empty:
            continue
        mlid = f.stem
        sid = canonical_site(mlid.split("-", 1)[0], mlid)
        loc = basins[basins["site_id"] == sid]
        if loc.empty:
            continue
        ng = g.to_crs(config.CRS_AREA).union_all()
        lg = loc.geometry.iloc[0]
        inter = ng.intersection(lg).area
        rows.append({"site_id": sid, "nldi_km2": ng.area / 1e6, "local_km2": lg.area / 1e6,
                     "overlap_iou": inter / ng.union(lg).area})
    t = pd.DataFrame(rows)
    if len(t):
        t.to_csv(config.TAB_DIR / "basin_check_nldi_vs_local.csv", index=False, float_format="%.3f")
        print("Local delineation vs NLDI basins:")
        print(t.round(3).to_string(index=False))


def nlcd_pct(src, geom):
    data, _ = mask(src, [geom], crop=True, filled=True, nodata=0)
    v = data[0]
    v = v[v > 0]
    counts = np.bincount(v, minlength=256)
    total = counts.sum()
    out = {"_cells": int(total)}
    for g, codes in config.NLCD_GROUPS.items():
        out[f"pct_{g}"] = 100 * counts[codes].sum() / total if total else np.nan
    return out


def main():
    sites = pd.read_csv(config.DATA_DIR / "wq_sites.csv", dtype={"huc12": str})
    sites = sites[sites["used_in_medians"]]
    basins = delineate(sites)
    print(f"Sites used in medians: {len(sites)}; NHDPlus basins delineated: {len(basins)}")
    nldi_check(basins)

    h12 = gpd.read_file(config.GPKG, layer="huc12").to_crs(config.CRS_AREA)
    study = h12.union_all()
    link = sites[["site_id", "MonitoringLocationName", "huc12"]].merge(
        basins[["site_id", "comid", "snap_m", "nhdplus_totdasqkm", "geometry"]], on="site_id", how="left")
    link["link_method"] = np.where(link["geometry"].notna(), "nhdplus_basin", "huc12_local")
    h12_geom = h12.set_index("huc12").geometry
    link.loc[link["geometry"].isna(), "geometry"] = link.loc[
        link["geometry"].isna(), "huc12"].map(h12_geom)
    link = gpd.GeoDataFrame(link, geometry="geometry", crs=config.CRS_AREA)
    link["basin_area_km2"] = link.area / 1e6
    link["frac_in_study_area"] = link.geometry.intersection(study).area / link.area

    rows = []
    for year in config.NLCD_YEARS:
        with rasterio.open(config.DATA_RAW / f"nlcd{year}_study_area.tif") as src:
            for i, geom in enumerate(link.geometry):
                r = nlcd_pct(src, geom)
                r = {f"{k}_{year}": v for k, v in r.items()}
                r["i"] = i
                rows.append(r)
    lc = pd.DataFrame(rows).groupby("i").first()
    link = link.join(lc)
    link["nlcd_coverage"] = link[f"_cells_{config.NLCD_YEARS[-1]}"] * CELL_KM2 / link["basin_area_km2"]

    print("Overlaying NWI with site basins")
    nwi = read_nwi(verbose=False)
    veg = nwi[nwi["group"] == "vegetated_wetland"][["geometry"]]
    lb = link[["site_id", "geometry"]].copy()
    inter = gpd.overlay(veg, lb, how="intersection", keep_geom_type=True)
    vw = (inter.area / 1e6).groupby(inter["site_id"]).sum()
    link["nwi_vegwet_km2"] = link["site_id"].map(vw).fillna(0.0)
    link["pct_nwi_vegwet"] = 100 * link["nwi_vegwet_km2"] / (link["basin_area_km2"] * link["frac_in_study_area"])
    # NLCD was downloaded for a padded box, so it covers whole basins. NWI was
    # downloaded only for the three HUC8 units, so NWI % is left blank when
    # less than 95% of a basin lies inside them.
    link["coverage_ok"] = link["nlcd_coverage"] >= 0.95
    link.loc[link["frac_in_study_area"] < 0.95, "pct_nwi_vegwet"] = np.nan

    drop = [c for c in link.columns if c.startswith("_cells")]
    out = link.drop(columns=drop + ["geometry"])
    out.to_csv(config.TAB_DIR / "site_basin_landcover.csv", index=False, float_format="%.3f")
    simp = link[["site_id", "link_method", "basin_area_km2", "geometry"]].copy()
    simp["geometry"] = simp.geometry.simplify(60)
    simp.to_file(config.GPKG, layer="site_basins", driver="GPKG")

    print(out["link_method"].value_counts().to_string())
    print(f"NLCD coverage OK: {out['coverage_ok'].sum()} of {len(out)}; "
          f"NWI % available: {out['pct_nwi_vegwet'].notna().sum()}")
    print(f"Basin area km2: median {out['basin_area_km2'].median():.1f}, "
          f"range {out['basin_area_km2'].min():.2f} to {out['basin_area_km2'].max():,.1f}")
    # Check basin areas against drainage areas listed in WQP station metadata
    sta = pd.read_csv(config.DATA_RAW / "wqp_stations.zip", dtype=str)
    da = pd.to_numeric(sta["DrainageAreaMeasure/MeasureValue"], errors="coerce")
    unit = sta["DrainageAreaMeasure/MeasureUnitCode"].fillna("").str.lower()
    sta["listed_km2"] = np.where(unit.str.contains("mi"), da * 2.58999, np.nan)
    chk = out[out["link_method"] == "nhdplus_basin"].merge(
        sta[["MonitoringLocationIdentifier", "listed_km2"]].dropna(),
        left_on="site_id", right_on="MonitoringLocationIdentifier")
    chk = chk[["site_id", "basin_area_km2", "nhdplus_totdasqkm", "listed_km2"]]
    chk.to_csv(config.TAB_DIR / "basin_area_check_wqp_metadata.csv", index=False, float_format="%.2f")
    print("Local basin area vs drainage area listed in WQP metadata:")
    print(chk.round(1).to_string(index=False))

if __name__ == "__main__":
    main()
