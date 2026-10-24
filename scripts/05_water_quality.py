"""
05_water_quality.py: clean Water Quality Portal (WQP) results, link sites
to HUC12 units, and compute per-site medians.

Cleaning rules (all applied in code below)
  * Keep water samples; drop quality-control activities (replicates, blanks).
  * Keep numeric results; a "Not Detected" result is set to half of its
    reported detection limit, or dropped if no limit is given.
  * Parameters and unit handling:
      SC_uScm    Specific conductance; umho/cm and uS/cm kept. Values labelled
                 mS/cm are multiplied by 1000 only when below 10; values of 10
                 or more labelled mS/cm (all from one LDEQ project, range
                 about 18 to 620) are fresh-water uS/cm magnitudes with a unit
                 label error and are kept as uS/cm.
      DO_mgL     Dissolved oxygen, mg/L
      pH         pH, standard units
      Temp_C     Water temperature, deg C
      NOx_N_mgL  Nitrate plus nitrite as N. If a sample has no nitrate plus
                 nitrite value, nitrate as N is used. Results reported as NO3
                 or in ueq/L are not used.
      TP_mgL     Total phosphorus (Phosphorus, Total fraction), mg/L as P;
                 ug/L divided by 1000
      Turb_NTU   Turbidity in NTU, FNU or NTRU (treated as comparable)
      DOC_mgL    Organic carbon, Dissolved fraction, mg/L
  * One value per site, parameter and activity (the first by priority).
  * Plausibility screens: pH 3 to 11, SC > 0, DO 0 to 25, Temp -1 to 40.
  * LDEQ reports the same site under two organization codes (LADEQWPD and
    LADEQWPD_WQX with the same number). These are merged into one site.

Medians use samples from WQ_MEDIAN_START to the end of the download window,
the years closest to NLCD 2021. A site needs MIN_SAMPLES_PER_SITE samples.

Outputs
  data/wq_results_clean.csv          cleaned long table (public WQP data)
  data/wq_sites.csv                  one row per site with HUC12 assignment
  outputs/tables/wq_site_medians.csv site x parameter medians and counts
  outputs/tables/wq_parameter_summary.csv
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import geopandas as gpd
import numpy as np
import pandas as pd

import config
from common import canonical_site

WQ_MEDIAN_START = 2011

# (parameter, WQP CharacteristicName, allowed fractions or None, unit -> factor, priority)
RULES = [
    ("SC_uScm", "Specific conductance", None,
     {"umho/cm": 1, "us/cm": 1, "us/cm @25c": 1, "ms/cm": 1000}, 1),
    ("DO_mgL", "Dissolved oxygen (DO)", None, {"mg/l": 1}, 1),
    ("pH", "pH", None, {"std units": 1, "": 1}, 1),
    ("Temp_C", "Temperature, water", None, {"deg c": 1}, 1),
    ("NOx_N_mgL", "Nitrate + Nitrite", None, {"mg/l": 1, "mg/l as n": 1}, 1),
    ("NOx_N_mgL", "Inorganic nitrogen (nitrate and nitrite)", None, {"mg/l as n": 1}, 1),
    ("NOx_N_mgL", "Nitrate", None, {"mg/l as n": 1, "mg/l": 1}, 2),
    ("TP_mgL", "Phosphorus", {"Total"}, {"mg/l": 1, "mg/l as p": 1, "ppm": 1}, 1),
    ("TP_mgL", "Total Phosphorus, mixed forms", None, {"mg/l": 1, "ug/l": 0.001}, 2),
    ("Turb_NTU", "Turbidity", None, {"ntu": 1, "fnu": 1, "ntru": 1}, 1),
    ("DOC_mgL", "Organic carbon", {"Dissolved"}, {"mg/l": 1}, 1),
]
LIMITS = {"pH": (3, 11), "SC_uScm": (5, 60000), "DO_mgL": (0, 25), "Temp_C": (-1, 40)}
PARAM_LABELS = {
    "SC_uScm": "Specific conductance (uS/cm)",
    "DO_mgL": "Dissolved oxygen (mg/L)",
    "pH": "pH",
    "Temp_C": "Water temperature (deg C)",
    "NOx_N_mgL": "Nitrate plus nitrite (mg/L as N)",
    "TP_mgL": "Total phosphorus (mg/L as P)",
    "Turb_NTU": "Turbidity (NTU)",
    "DOC_mgL": "Dissolved organic carbon (mg/L)",
}


def clean(res):
    res = res[res["ActivityMediaName"].eq("Water")]
    res = res[~res["ActivityTypeCode"].str.contains("Quality Control", na=False)].copy()
    res["unit"] = res["ResultMeasure/MeasureUnitCode"].fillna("").str.strip().str.lower()
    res["value"] = pd.to_numeric(res["ResultMeasureValue"], errors="coerce")
    nd = res["ResultDetectionConditionText"].eq("Not Detected")
    dl = pd.to_numeric(res["DetectionQuantitationLimitMeasure/MeasureValue"], errors="coerce")
    dl_unit = res["DetectionQuantitationLimitMeasure/MeasureUnitCode"].fillna("").str.lower()
    res.loc[nd, "value"] = dl[nd] / 2
    res.loc[nd & (res["unit"] == ""), "unit"] = dl_unit[nd]
    res["censored"] = nd

    out = []
    for param, char, fracs, units, prio in RULES:
        sub = res[res["CharacteristicName"].eq(char)]
        if fracs is not None:
            sub = sub[sub["ResultSampleFractionText"].isin(fracs)]
        sub = sub[sub["unit"].isin(units)]
        factor = sub["unit"].map(units)
        if param == "SC_uScm":
            mislabel = sub["unit"].eq("ms/cm") & (sub["value"] >= 10)
            factor = factor.where(~mislabel, 1)
        sub = sub.assign(param=param, priority=prio, value=sub["value"] * factor)
        out.append(sub)
    df = pd.concat(out)
    df = df.dropna(subset=["value"])
    for p, (lo, hi) in LIMITS.items():
        bad = df["param"].eq(p) & ~df["value"].between(lo, hi)
        df = df[~bad]
    df["site_id"] = [canonical_site(o, s) for o, s in
                     zip(df["OrganizationIdentifier"], df["MonitoringLocationIdentifier"])]
    df["date"] = pd.to_datetime(df["ActivityStartDate"])
    df = (df.sort_values(["priority"])
            .drop_duplicates(["site_id", "param", "ActivityIdentifier"])
            .drop_duplicates(["site_id", "param", "date", "value"]))
    keep = ["site_id", "MonitoringLocationIdentifier", "OrganizationIdentifier", "date",
            "param", "value", "censored", "CharacteristicName", "ResultSampleFractionText",
            "unit", "ActivityTypeCode"]
    return df[keep].sort_values(["site_id", "param", "date"]).reset_index(drop=True)


def main():
    res = pd.read_csv(config.DATA_RAW / "wqp_results.zip", low_memory=False, dtype=str)
    sta = pd.read_csv(config.DATA_RAW / "wqp_stations.zip", low_memory=False, dtype=str)
    print(f"WQP rows downloaded: {len(res):,}; stations: {len(sta)}")
    df = clean(res)
    print(f"Rows after cleaning: {len(df):,}")
    df.to_csv(config.DATA_DIR / "wq_results_clean.csv", index=False, float_format="%.6g")

    sta["site_id"] = [canonical_site(o, s) for o, s in
                      zip(sta["OrganizationIdentifier"], sta["MonitoringLocationIdentifier"])]
    sta["lat"] = pd.to_numeric(sta["LatitudeMeasure"], errors="coerce")
    sta["lon"] = pd.to_numeric(sta["LongitudeMeasure"], errors="coerce")
    sites = (sta.dropna(subset=["lat", "lon"])
               .sort_values("OrganizationIdentifier", ascending=False)  # prefer *_WQX names
               .drop_duplicates("site_id")
               [["site_id", "OrganizationIdentifier", "MonitoringLocationName",
                 "MonitoringLocationTypeName", "StateCode", "lat", "lon"]])
    sites = sites[sites["site_id"].isin(df["site_id"])]
    g = gpd.GeoDataFrame(sites, geometry=gpd.points_from_xy(sites.lon, sites.lat), crs=config.CRS_GEO)
    h12 = gpd.read_file(config.GPKG, layer="huc12")[["huc12", "name", "upstream_area_km2", "geometry"]]
    g = gpd.sjoin(g, h12.to_crs(config.CRS_GEO), how="left", predicate="within").drop(columns="index_right")
    g = g.rename(columns={"name": "huc12_name"})
    print(f"Sites with data: {len(g)}; inside a study-area HUC12: {g['huc12'].notna().sum()}")

    win = df[df["date"].dt.year >= WQ_MEDIAN_START]
    stats = (win.groupby(["site_id", "param"])
                .agg(n=("value", "size"), median=("value", "median"),
                     n_censored=("censored", "sum"),
                     first=("date", "min"), last=("date", "max"))
                .reset_index())
    stats = stats[stats["n"] >= config.MIN_SAMPLES_PER_SITE]
    stats["first_year"] = stats["first"].dt.year
    stats["last_year"] = stats["last"].dt.year
    stats = stats.drop(columns=["first", "last"])
    stats = stats.merge(g[["site_id", "MonitoringLocationName", "huc12", "huc12_name"]], on="site_id")
    stats = stats[stats["huc12"].notna()]
    stats.to_csv(config.TAB_DIR / "wq_site_medians.csv", index=False, float_format="%.4g")

    n_sites = stats.groupby("site_id").size()
    g["used_in_medians"] = g["site_id"].isin(n_sites.index)
    g.drop(columns="geometry").to_csv(config.DATA_DIR / "wq_sites.csv", index=False)
    g.to_file(config.GPKG, layer="wq_sites", driver="GPKG")

    summ = (stats.groupby("param")
                 .agg(n_sites=("site_id", "nunique"), n_samples=("n", "sum"),
                      median_of_site_medians=("median", "median"),
                      min_site_median=("median", "min"), max_site_median=("median", "max"),
                      n_huc12=("huc12", "nunique")))
    summ["label"] = summ.index.map(PARAM_LABELS)
    summ.to_csv(config.TAB_DIR / "wq_parameter_summary.csv", float_format="%.4g")
    print(f"\nSite medians, {WQ_MEDIAN_START} onward, sites with >= {config.MIN_SAMPLES_PER_SITE} samples")
    print(summ.drop(columns="label").to_string())
    print(f"Sites used: {stats['site_id'].nunique()}, in {stats['huc12'].nunique()} HUC12 units")
    print("Organizations:", g[g.used_in_medians]["OrganizationIdentifier"].value_counts().to_dict())


if __name__ == "__main__":
    main()
