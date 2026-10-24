"""
08_discharge_context.py: streamflow context for the water-quality samples.

In the WQP download used here (stream sites, 2000 onward) there are no USGS
discrete samples at the three named gauges (07378500, 07377500, 07380120);
their WQP records for the parameters checked predate 2000. Instead, state
and federal samples on the Amite River main stem are paired with same-day
mean discharge at the nearest Amite main-stem gauge.

Steps
  1. Read NWIS daily mean discharge (parameter 00060, ft3/s) for 2000 to 2025.
  2. Pick sites whose name starts with "Amite River" and that lie within
     CQ_MAX_KM of an Amite main-stem gauge (07378500 or 07380120). Each site
     is paired with its nearest gauge.
  3. For SC, turbidity and TP, per gauge: log10(C) vs log10(Q) regression
     and Spearman rho, using all samples from 2000 onward.
  4. Compare the discharge on sampling days with the full daily record, to
     show whether routine sampling over- or under-represents high flow.

Outputs
  outputs/tables/discharge_summary.csv
  outputs/tables/cq_amite_main_stem.csv
  outputs/figures/fig09_discharge_cq.png
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import io

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import config

GAUGE = "07378500"           # gauge shown in the hydrograph panel
AMITE_GAUGES = ["07378500", "07380120"]
CQ_MAX_KM = 5
CQ_PARAMS = {"SC_uScm": "Specific conductance (uS/cm)", "Turb_NTU": "Turbidity (NTU)",
             "TP_mgL": "Total phosphorus (mg/L)"}


def read_rdb(path):
    lines = [ln for ln in Path(path).read_text().splitlines() if not ln.startswith("#")]
    if len(lines) < 2:
        return pd.DataFrame()
    body = "\n".join([lines[0]] + lines[2:])          # skip the column-format row
    return pd.read_csv(io.StringIO(body), sep="\t", dtype=str)


def read_dv():
    df = read_rdb(config.DATA_RAW / "nwis_dv_discharge.rdb")
    # One RDB file holds several sites; each block repeats the header row
    df = df[df["agency_cd"] == "USGS"].copy()
    qcol = [c for c in df.columns if c.endswith("_00060_00003")]
    df["q_cfs"] = pd.to_numeric(df[qcol].bfill(axis=1).iloc[:, 0], errors="coerce")
    df["date"] = pd.to_datetime(df["datetime"])
    return df[["site_no", "date", "q_cfs"]].dropna()


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = p2 - p1
    dl = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * r * np.arcsin(np.sqrt(a))


def main():
    dv = read_dv()
    summ = (dv.groupby("site_no")["q_cfs"]
              .agg(n_days="size", median_cfs="median",
                   p10_cfs=lambda s: s.quantile(0.10), p90_cfs=lambda s: s.quantile(0.90),
                   max_cfs="max"))
    summ = summ.join(dv.groupby("site_no")["date"].agg(first="min", last="max"))
    summ["name"] = summ.index.map(config.GAUGES)
    summ.to_csv(config.TAB_DIR / "discharge_summary.csv", float_format="%.1f")
    print(summ.to_string())

    gs = read_rdb(config.DATA_RAW / "nwis_gauge_sites.rdb")
    gs = gs[gs["site_no"].isin(AMITE_GAUGES)]
    sites = pd.read_csv(config.DATA_DIR / "wq_sites.csv")
    sites = sites[sites["MonitoringLocationName"].str.startswith("Amite River")].copy()
    for _, g in gs.iterrows():
        sites[f"km_{g['site_no']}"] = haversine_km(sites["lat"], sites["lon"],
                                                   float(g["dec_lat_va"]), float(g["dec_long_va"]))
    kcols = [f"km_{g}" for g in gs["site_no"]]
    sites["gauge"] = sites[kcols].idxmin(axis=1).str[3:]
    sites["km_to_gauge"] = sites[kcols].min(axis=1)
    near = sites[sites["km_to_gauge"] <= CQ_MAX_KM]
    print(f"\nAmite main-stem sites within {CQ_MAX_KM} km of an Amite gauge:")
    print(near[["site_id", "MonitoringLocationName", "gauge", "km_to_gauge"]].round(1).to_string(index=False))

    wq = pd.read_csv(config.DATA_DIR / "wq_results_clean.csv", parse_dates=["date"])
    wq = wq[wq["site_id"].isin(near["site_id"]) & wq["param"].isin(CQ_PARAMS)]
    wq = wq.merge(near[["site_id", "gauge"]], on="site_id")
    pair = wq.merge(dv.rename(columns={"site_no": "gauge"}), on=["gauge", "date"], how="inner")
    pair = pair[(pair["q_cfs"] > 0) & (pair["value"] > 0)]

    rows = []
    for (gauge, p), d in pair.groupby(["gauge", "param"]):
        if len(d) < 10:
            continue
        lr = stats.linregress(np.log10(d["q_cfs"]), np.log10(d["value"]))
        rho, pr = stats.spearmanr(d["q_cfs"], d["value"])
        allq = dv[(dv["site_no"] == gauge) & dv["date"].between(d["date"].min(), d["date"].max())]["q_cfs"]
        rows.append({"gauge": gauge, "param": p, "n_samples": len(d), "n_sites": d["site_id"].nunique(),
                     "sites": ";".join(sorted(d["site_id"].unique())),
                     "first": d["date"].min().date(), "last": d["date"].max().date(),
                     "loglog_slope": lr.slope, "loglog_r2": lr.rvalue ** 2, "loglog_p": lr.pvalue,
                     "spearman_rho": rho, "spearman_p": pr,
                     "median_q_sample_days_cfs": d["q_cfs"].median(),
                     "median_q_all_days_same_period_cfs": allq.median(),
                     "sample_day_median_q_percentile": 100 * (allq < d["q_cfs"].median()).mean()})
    cq = pd.DataFrame(rows)
    cq.to_csv(config.TAB_DIR / "cq_amite_main_stem.csv", index=False, float_format="%.4g")
    print("\nConcentration-discharge, Amite River main stem")
    print(cq.drop(columns="sites").round(3).to_string(index=False))
    q = dv[dv["site_no"] == GAUGE][["date", "q_cfs"]]
    pair_g = pair[pair["gauge"] == GAUGE]
    cq_g = cq[cq["gauge"] == GAUGE]

    # Figure
    fig = plt.figure(figsize=(10, 6.4))
    ax0 = fig.add_axes([0.07, 0.58, 0.9, 0.36])
    qq = q[q["date"].dt.year >= 2000]
    ax0.plot(qq["date"], qq["q_cfs"], color="#2a78d6", linewidth=0.6)
    sd = pair_g[pair_g["param"] == "SC_uScm"].drop_duplicates("date")
    ax0.scatter(sd["date"], sd["q_cfs"], s=10, color="#eb6834", zorder=3, label="Sampling day (SC sample)")
    ax0.set_yscale("log")
    ax0.set_ylabel("Daily mean Q (ft3/s)", fontsize=8)
    ax0.set_title(f"USGS {GAUGE} {config.GAUGES[GAUGE]}: daily discharge and sampling days",
                  fontsize=9)
    ax0.legend(fontsize=7.5, loc="upper left")
    ax0.grid(color="#e6e5e0", linewidth=0.5)
    ax0.tick_params(labelsize=7)
    for k, p in enumerate(cq_g["param"]):
        ax = fig.add_axes([0.09 + k * 0.32, 0.08, 0.24, 0.36])
        d = pair_g[pair_g["param"] == p]
        ax.scatter(d["q_cfs"], d["value"], s=12, color="#2a78d6", edgecolor="white", linewidth=0.4)
        r = cq_g[cq_g["param"] == p].iloc[0]
        xs = np.logspace(np.log10(d["q_cfs"].min()), np.log10(d["q_cfs"].max()), 50)
        lr = stats.linregress(np.log10(d["q_cfs"]), np.log10(d["value"]))
        ax.plot(xs, 10 ** (lr.intercept + lr.slope * np.log10(xs)), color="#0b0b0b", linewidth=1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(f"Same-day Q at {GAUGE} (ft3/s)", fontsize=8)
        ax.set_ylabel(CQ_PARAMS[p], fontsize=8)
        ax.text(0.04, 0.05, f"slope {r['loglog_slope']:+.2f}, R2 {r['loglog_r2']:.2f}\n"
                f"n = {int(r['n_samples'])}", transform=ax.transAxes, fontsize=7.5, color="#52514e")
        ax.grid(color="#e6e5e0", linewidth=0.5)
        ax.tick_params(which="both", labelsize=7)
    fig.savefig(config.FIG_DIR / "fig09_discharge_cq.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
