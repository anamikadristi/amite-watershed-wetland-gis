"""
07_relationships.py: relationships between basin land cover and site
water-quality medians.

For each parameter and each land-cover predictor:
  * Spearman rank correlation (rho, p) across sites
  * simple least-squares regression of the site median on the predictor
    (log10 of the median for skewed parameters: SC, NOx, TP, turbidity)
  * Benjamini-Hochberg q-values across all tests within a site subset

Site subsets
  all          every site with a basin link and >= 95% data coverage
  no_maurepas  drops sites in the Lake Maurepas HUC8 (08070204), where tidal
               exchange and backwater from the lake, not upstream land cover,
               are likely to control salinity and nutrients
  no_ldeq319   drops the LDEQ nonpoint-source (319) project sites, which are
               clustered in the urbanizing Comite and Baton Rouge area and
               were sampled over shorter periods than the ambient network

These are observational, cross-sectional associations among a small number
of sites that are partly nested (downstream basins contain upstream ones).
They are not causal estimates.

Outputs
  outputs/tables/wq_landcover_tests.csv
  outputs/tables/wq_site_medians_with_landcover.csv
  outputs/figures/fig07_spearman_heatmap.png
  outputs/figures/fig08_scatter_landcover_wq.png
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import config
from common import bh_fdr

PARAMS = ["SC_uScm", "DO_mgL", "pH", "Temp_C", "NOx_N_mgL", "TP_mgL", "Turb_NTU", "DOC_mgL"]
LOG_PARAMS = {"SC_uScm", "NOx_N_mgL", "TP_mgL", "Turb_NTU", "DOC_mgL"}
LABELS = {
    "SC_uScm": "Specific conductance", "DO_mgL": "Dissolved oxygen", "pH": "pH",
    "Temp_C": "Water temperature", "NOx_N_mgL": "Nitrate + nitrite N",
    "TP_mgL": "Total phosphorus", "Turb_NTU": "Turbidity", "DOC_mgL": "Dissolved organic carbon",
}
UNITS = {"SC_uScm": "uS/cm", "DO_mgL": "mg/L", "pH": "", "Temp_C": "deg C",
         "NOx_N_mgL": "mg/L", "TP_mgL": "mg/L", "Turb_NTU": "NTU", "DOC_mgL": "mg/L"}
Y = config.NLCD_YEARS[-1]
PREDICTORS = {
    f"pct_wetland_{Y}": "NLCD wetland %",
    "pct_nwi_vegwet": "NWI vegetated wetland %",
    f"pct_forest_{Y}": "NLCD forest %",
    f"pct_developed_{Y}": "NLCD developed %",
    f"pct_agriculture_{Y}": "NLCD agriculture %",
    "log10_basin_km2": "log10 basin area",
}
HUC8_COLORS = {"Amite": "#2a78d6", "Tickfaw": "#eb6834", "Lake Maurepas": "#1baf7a"}


def run_tests(df, subset):
    rows = []
    for p in PARAMS:
        d = df[df["param"] == p]
        for pred, plab in PREDICTORS.items():
            x = d[pred].to_numpy(float)
            y = d["median"].to_numpy(float)
            ok = ~np.isnan(x) & ~np.isnan(y)
            x, y = x[ok], y[ok]
            n = len(x)
            r = {"subset": subset, "param": p, "predictor": pred, "predictor_label": plab, "n_sites": n}
            if n >= config.MIN_SITES_FOR_TEST and np.ptp(x) > 0:
                rho, pr = stats.spearmanr(x, y)
                yy = np.log10(y) if p in LOG_PARAMS else y
                lr = stats.linregress(x, yy)
                r.update(spearman_rho=rho, spearman_p=pr,
                         ols_response="log10(median)" if p in LOG_PARAMS else "median",
                         ols_slope=lr.slope, ols_intercept=lr.intercept,
                         ols_r2=lr.rvalue ** 2, ols_p=lr.pvalue)
            rows.append(r)
    t = pd.DataFrame(rows)
    if "spearman_p" in t:
        t["spearman_q_bh"] = bh_fdr(t["spearman_p"])
    return t


def main():
    med = pd.read_csv(config.TAB_DIR / "wq_site_medians.csv", dtype={"huc12": str})
    lc = pd.read_csv(config.TAB_DIR / "site_basin_landcover.csv", dtype={"huc12": str})
    lc["log10_basin_km2"] = np.log10(lc["basin_area_km2"])
    lc["huc8_name"] = lc["huc12"].str[:8].map(config.HUC8)
    df = med.merge(lc.drop(columns=["MonitoringLocationName", "huc12"]), on="site_id", how="inner")
    df = df[df["coverage_ok"]]
    df.to_csv(config.TAB_DIR / "wq_site_medians_with_landcover.csv", index=False, float_format="%.4g")
    print(f"Site x parameter rows with land cover: {len(df)}; sites: {df['site_id'].nunique()}")

    subsets = {
        "all": df,
        "no_maurepas": df[df["huc8_name"] != "Lake Maurepas"],
        "no_ldeq319": df[~df["site_id"].str.startswith("LDEQ319")],
    }
    tests = pd.concat([run_tests(d, k) for k, d in subsets.items()], ignore_index=True)
    tests.to_csv(config.TAB_DIR / "wq_landcover_tests.csv", index=False, float_format="%.4g")

    for k in subsets:
        t = tests[(tests["subset"] == k) & tests["spearman_rho"].notna()]
        print(f"\nSubset {k}: {len(t)} tests, sites {subsets[k]['site_id'].nunique()}, "
              f"p<0.05: {(t['spearman_p'] < 0.05).sum()}, q<0.05: {(t['spearman_q_bh'] < 0.05).sum()}")
        show = t[t["spearman_q_bh"] < 0.05].sort_values("spearman_p")
        print(show[["param", "predictor_label", "n_sites", "spearman_rho", "spearman_p",
                    "spearman_q_bh", "ols_r2"]].round(4).to_string(index=False))

    # Heatmap, subset "all"
    t = tests[tests["subset"] == "all"]
    params = [p for p in PARAMS if t[(t["param"] == p)]["spearman_rho"].notna().any()]
    preds = list(PREDICTORS)
    M = np.full((len(params), len(preds)), np.nan)
    Q = np.full_like(M, np.nan)
    N = np.zeros_like(M)
    for i, p in enumerate(params):
        for j, pr in enumerate(preds):
            row = t[(t["param"] == p) & (t["predictor"] == pr)].iloc[0]
            M[i, j] = row.get("spearman_rho", np.nan)
            Q[i, j] = row.get("spearman_q_bh", np.nan)
            N[i, j] = row["n_sites"]
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    im = ax.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    for i in range(len(params)):
        for j in range(len(preds)):
            if np.isnan(M[i, j]):
                continue
            mark = "*" if Q[i, j] < 0.05 else ""
            ax.text(j, i, f"{M[i, j]:+.2f}{mark}", ha="center", va="center", fontsize=8,
                    color="white" if abs(M[i, j]) > 0.6 else "#0b0b0b")
    ax.set_xticks(range(len(preds)), [PREDICTORS[p] for p in preds], rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(len(params)),
                  [f"{LABELS[p]} (n={int(N[i].max())})" for i, p in enumerate(params)], fontsize=8)
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("Spearman rho", fontsize=8)
    ax.set_title("Site median water quality vs. upstream basin cover (NLCD 2021, NWI)\n"
                 "* Benjamini-Hochberg q < 0.05 across all tests shown", fontsize=9)
    fig.tight_layout()
    fig.savefig(config.FIG_DIR / "fig07_spearman_heatmap.png", dpi=170)
    plt.close(fig)

    # Scatter small multiples: fixed set of parameters and predictors
    sp = [p for p in ["SC_uScm", "DO_mgL", "TP_mgL", "NOx_N_mgL"] if p in params]
    spred = [f"pct_wetland_{Y}", f"pct_forest_{Y}", f"pct_developed_{Y}"]
    fig, axes = plt.subplots(len(sp), len(spred), figsize=(9, 2.4 * len(sp)), squeeze=False)
    d_all = subsets["all"]
    for i, p in enumerate(sp):
        d = d_all[d_all["param"] == p]
        for j, pr in enumerate(spred):
            ax = axes[i, j]
            for h8, col in HUC8_COLORS.items():
                s = d[d["huc8_name"] == h8]
                ax.scatter(s[pr], s["median"], s=22, color=col, edgecolor="white", linewidth=0.6,
                           label=h8, zorder=3)
            if p in LOG_PARAMS:
                ax.set_yscale("log")
            row = t[(t["param"] == p) & (t["predictor"] == pr)].iloc[0]
            if not np.isnan(row.get("spearman_rho", np.nan)):
                ax.set_title(f"rho = {row['spearman_rho']:+.2f}, n = {int(row['n_sites'])}",
                             fontsize=7.5, color="#52514e", loc="left")
            ax.grid(color="#e6e5e0", linewidth=0.5)
            ax.tick_params(labelsize=7)
            if i == len(sp) - 1:
                ax.set_xlabel(PREDICTORS[pr], fontsize=8)
            if j == 0:
                ax.set_ylabel(f"{LABELS[p]}\n({UNITS[p]})" if UNITS[p] else LABELS[p], fontsize=8)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=8, frameon=False, title="HUC8 unit",
               title_fontsize=8)
    fig.suptitle("Site medians (2011 onward) against basin land cover, NLCD 2021", fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(config.FIG_DIR / "fig08_scatter_landcover_wq.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
