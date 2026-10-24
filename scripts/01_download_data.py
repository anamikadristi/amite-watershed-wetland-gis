"""
01_download_data.py: download every public input used by the project.

Sources (all public, no API key needed):
  1. USGS Watershed Boundary Dataset (WBD), HUC8 and HUC12 polygons,
     The National Map ArcGIS REST service.
  2. NLCD 2001 and 2021 land cover (MRLC GeoServer WCS), clipped to the
     study-area bounding box on the server.
  3. USFWS National Wetlands Inventory (NWI), HU8 watershed downloads.
  4. NHDPlus V2 flowlines (stream order 3 and up) and large waterbodies,
     USGS WMA GeoServer WFS.
  5. Water Quality Portal (WQP) discrete results and station metadata for
     stream sites in the three HUC8 units.
  6. USGS NWIS daily mean discharge for three gauges.
  7. NHDPlus V2 network attributes (all flowlines) and catchment polygons
     for the study-area bounding box, used to delineate the upstream basin
     of each water-quality site locally (script 06).
  8. Optional (--nldi): USGS NLDI basins for a few WQP sites, used only as
     an independent check of the local delineation. The NLDI service can
     take minutes per basin.

Files are written to config.DATA_RAW (default: data_raw/, ignored by Git).
A manifest with URL, access date and size is written to
data/download_manifest.csv.

Run:  python scripts/01_download_data.py          (skips files already present)
      python scripts/01_download_data.py --force  (downloads again)
      python scripts/01_download_data.py --nldi   (also fetch NLDI check basins)
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import csv
import datetime as dt
import time

import requests

import config

FORCE = "--force" in sys.argv
RAW = config.DATA_RAW
SESSION = requests.Session()
SESSION.headers["User-Agent"] = "amite-watershed-wetland-gis (public research script)"
MANIFEST = []


def fetch(url, out, params=None, timeout=600, tries=3):
    """Download url to out (Path), with retries. Records the file in MANIFEST."""
    out = Path(out)
    if out.exists() and out.stat().st_size > 0 and not FORCE:
        print(f"  exists, skipping: {out.name}")
    else:
        for attempt in range(1, tries + 1):
            try:
                with SESSION.get(url, params=params, timeout=timeout, stream=True) as r:
                    r.raise_for_status()
                    tmp = out.with_suffix(out.suffix + ".part")
                    with open(tmp, "wb") as f:
                        for chunk in r.iter_content(1 << 20):
                            f.write(chunk)
                    tmp.rename(out)
                break
            except Exception as err:
                print(f"  attempt {attempt} failed for {out.name}: {str(err)[:150]}")
                if attempt == tries:
                    raise
                time.sleep(10 * attempt)
        print(f"  saved {out.name} ({out.stat().st_size / 1e6:.2f} MB)")
    full = requests.Request("GET", url, params=params).prepare().url
    MANIFEST.append({
        "file": out.name,
        "url": full,
        "accessed": dt.date.fromtimestamp(out.stat().st_mtime).isoformat(),
        "size_mb": round(out.stat().st_size / 1e6, 3),
    })
    return out


# --------------------------------------------------------------------------
# 1. WBD
# --------------------------------------------------------------------------
WBD = "https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer"


def download_wbd():
    print("WBD HUC12 and HUC8 polygons")
    codes = list(config.HUC8)
    where12 = " OR ".join(f"huc12 LIKE '{c}%'" for c in codes)
    fetch(f"{WBD}/6/query", RAW / "wbd_huc12.geojson", params={
        "where": where12, "outFields": "huc12,name,tohuc,areasqkm,hutype",
        "outSR": 4326, "f": "geojson", "returnGeometry": "true",
    })
    where8 = " OR ".join(f"huc8 = '{c}'" for c in codes)
    fetch(f"{WBD}/4/query", RAW / "wbd_huc8.geojson", params={
        "where": where8, "outFields": "huc8,name,areasqkm",
        "outSR": 4326, "f": "geojson", "returnGeometry": "true",
    })


def study_bbox_albers(pad_m=2000):
    """Bounding box of the HUC12 layer in EPSG:5070, snapped to the NLCD grid."""
    import geopandas as gpd
    h = gpd.read_file(RAW / "wbd_huc12.geojson").to_crs(config.CRS_AREA)
    minx, miny, maxx, maxy = h.total_bounds
    # NLCD L48 grid origin is x = -2493045, y = 177285 with 30 m cells
    def snap(v, origin, down):
        k = (v - origin) / 30.0
        k = int(k) if down else int(k) + 1
        return origin + 30 * k
    return (snap(minx - pad_m, -2493045, True), snap(miny - pad_m, 177285, True),
            snap(maxx + pad_m, -2493045, False), snap(maxy + pad_m, 177285, False))


# --------------------------------------------------------------------------
# 2. NLCD via WCS (tiled, then mosaicked)
# --------------------------------------------------------------------------
WCS = "https://www.mrlc.gov/geoserver/mrlc_download/wcs"


def download_nlcd():
    import rasterio
    from rasterio.merge import merge

    minx, miny, maxx, maxy = study_bbox_albers()
    tile = 1500 * 30  # 45 km tiles
    for year in config.NLCD_YEARS:
        final = RAW / f"nlcd{year}_study_area.tif"
        cov = f"mrlc_download__NLCD_{year}_Land_Cover_L48"
        print(f"NLCD {year}")
        parts = []
        x = minx
        i = 0
        while x < maxx:
            y = miny
            while y < maxy:
                x2, y2 = min(x + tile, maxx), min(y + tile, maxy)
                part = RAW / f"nlcd{year}_tile{i:02d}.tif"
                fetch(WCS, part, params=[
                    ("service", "WCS"), ("version", "2.0.1"), ("request", "GetCoverage"),
                    ("coverageId", cov), ("subset", f"X({x},{x2})"),
                    ("subset", f"Y({y},{y2})"), ("format", "image/geotiff"),
                ])
                parts.append(part)
                i += 1
                y = y2
            x = x2
        if final.exists() and not FORCE:
            print(f"  exists, skipping mosaic: {final.name}")
            continue
        srcs = [rasterio.open(p) for p in parts]
        mosaic, transform = merge(srcs, nodata=0)
        meta = srcs[0].meta.copy()
        meta.update(height=mosaic.shape[1], width=mosaic.shape[2], transform=transform,
                    nodata=0, compress="lzw", driver="GTiff")
        for s in srcs:
            s.close()
        with rasterio.open(final, "w", **meta) as dst:
            dst.write(mosaic)
        print(f"  mosaic written: {final.name} {mosaic.shape[2]} x {mosaic.shape[1]} cells")


# --------------------------------------------------------------------------
# 3. NWI HU8 downloads
# --------------------------------------------------------------------------
NWI = "https://documentst.ecosphere.fws.gov/wetlands/downloads/watershed/HU8_{}_Watershed.zip"


def download_nwi():
    print("NWI wetlands by HU8")
    for code in config.HUC8:
        fetch(NWI.format(code), RAW / f"NWI_HU8_{code}_Watershed.zip", timeout=1200)


# --------------------------------------------------------------------------
# 4. NHDPlus V2 flowlines and waterbodies
# --------------------------------------------------------------------------
WFS = "https://api.water.usgs.gov/geoserver/wmadata/ows"


def download_nhd():
    import geopandas as gpd
    print("NHDPlus V2 flowlines and waterbodies")
    h = gpd.read_file(RAW / "wbd_huc12.geojson")
    x0, y0, x1, y1 = h.total_bounds
    box = f"BBOX(the_geom,{x0:.4f},{y0:.4f},{x1:.4f},{y1:.4f},'EPSG:4326')"
    common = {"service": "WFS", "version": "1.1.0", "request": "GetFeature",
              "outputFormat": "application/json"}
    fetch(WFS, RAW / "nhdplus_flowlines_order3plus.geojson", params={
        **common, "typeName": "wmadata:nhdflowline_network",
        "propertyName": "comid,gnis_name,streamorde,ftype,the_geom",
        "CQL_FILTER": f"streamorde >= 3 AND {box}",
    })
    fetch(WFS, RAW / "nhdplus_waterbodies_1km2.geojson", params={
        **common, "typeName": "wmadata:nhdwaterbody",
        "propertyName": "comid,gnis_name,areasqkm,ftype,the_geom",
        "CQL_FILTER": f"areasqkm >= 1 AND {box}",
    })


# --------------------------------------------------------------------------
# 5. Water Quality Portal
# --------------------------------------------------------------------------
WQP = "https://www.waterqualitydata.us/data"
CHARACTERISTICS = [
    "Organic carbon", "Specific conductance", "Dissolved oxygen (DO)", "pH",
    "Temperature, water", "Nitrate", "Nitrate + Nitrite",
    "Inorganic nitrogen (nitrate and nitrite)", "Phosphorus",
    "Total Phosphorus, mixed forms", "Turbidity",
]


def wqp_params():
    p = [("huc", ";".join(config.HUC8)), ("siteType", "Stream"),
         ("startDateLo", config.WQ_START), ("startDateHi", config.WQ_END),
         ("mimeType", "csv"), ("zip", "yes")]
    p += [("characteristicName", c) for c in CHARACTERISTICS]
    return p


def download_wqp():
    print("Water Quality Portal results and stations")
    fetch(f"{WQP}/Result/search", RAW / "wqp_results.zip", params=wqp_params(), timeout=1800)
    fetch(f"{WQP}/Station/search", RAW / "wqp_stations.zip", params=wqp_params(), timeout=900)


# --------------------------------------------------------------------------
# 6. NWIS daily discharge
# --------------------------------------------------------------------------
def download_nwis_dv():
    print("NWIS daily discharge")
    fetch("https://waterservices.usgs.gov/nwis/dv/", RAW / "nwis_dv_discharge.rdb", params={
        "sites": ",".join(config.GAUGES), "parameterCd": "00060",
        "startDT": "2000-01-01", "endDT": "2025-12-31", "format": "rdb",
    })
    fetch("https://waterservices.usgs.gov/nwis/site/", RAW / "nwis_gauge_sites.rdb", params={
        "sites": ",".join(config.GAUGES), "siteOutput": "expanded", "format": "rdb",
    })


# --------------------------------------------------------------------------
# 7. NHDPlus V2 network and catchments (paged WFS requests)
# --------------------------------------------------------------------------
NET_FIELDS = "comid,fromnode,tonode,hydroseq,dnhydroseq,divergence,streamorde,totdasqkm,gnis_name,ftype"


def download_nhd_network(page=2000):
    import geopandas as gpd
    print("NHDPlus V2 network attributes and catchments")
    h = gpd.read_file(RAW / "wbd_huc12.geojson")
    x0, y0, x1, y1 = h.total_bounds
    box = f"BBOX(the_geom,{x0 - 0.02:.4f},{y0 - 0.02:.4f},{x1 + 0.02:.4f},{y1 + 0.02:.4f},'EPSG:4326')"
    common = {"service": "WFS", "version": "1.1.0", "request": "GetFeature",
              "outputFormat": "application/json", "CQL_FILTER": box, "maxFeatures": page}
    jobs = [("nhdplus_network", "wmadata:nhdflowline_network", "comid", NET_FIELDS + ",the_geom"),
            ("nhdplus_catchments", "wmadata:catchmentsp", "featureid", "featureid,areasqkm,the_geom")]
    for stem, layer, sort, fields in jobs:
        hits = SESSION.get(WFS, params={**common, "typeName": layer, "resultType": "hits",
                                         "outputFormat": None, "maxFeatures": None}, timeout=300).text
        n = int(hits.split('numberOfFeatures="')[1].split('"')[0])
        print(f"  {layer}: {n} features")
        for k, start in enumerate(range(0, n, page)):
            fetch(WFS, RAW / f"{stem}_p{k:02d}.geojson", params={
                **common, "typeName": layer, "propertyName": fields,
                "sortBy": sort, "startIndex": start})


# --------------------------------------------------------------------------
# 8. Optional NLDI basins for a few sites (independent check)
# --------------------------------------------------------------------------
NLDI = "https://api.water.usgs.gov/nldi/linked-data"
NLDI_CHECK_SITES = ["LADEQWPD_WQX-0036", "LADEQWPD_WQX-0043", "LADEQWPD_WQX-0116",
                    "LADEQWPD_WQX-0228", "LADEQWPD_WQX-0233"]


def download_nldi_basins():
    out_dir = RAW / "nldi_basins"
    out_dir.mkdir(exist_ok=True)
    print(f"NLDI check basins for {len(NLDI_CHECK_SITES)} sites (slow service)")
    for site in NLDI_CHECK_SITES:
        try:
            fetch(f"{NLDI}/wqp/{site}/basin", out_dir / f"{site}.geojson", timeout=900, tries=2)
        except Exception as err:
            print(f"  NLDI failed for {site}: {str(err)[:100]}")


def main():
    print(f"Download folder: {RAW}")
    download_wbd()
    download_nlcd()
    download_nwi()
    download_nhd()
    download_wqp()
    download_nwis_dv()
    download_nhd_network()
    if "--nldi" in sys.argv:
        download_nldi_basins()
    # Keep one manifest row per file (tiles included), newest wins
    rows = {r["file"]: r for r in MANIFEST}
    out = config.DATA_DIR / "download_manifest.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["file", "url", "accessed", "size_mb"])
        w.writeheader()
        w.writerows(rows.values())
    print(f"Manifest written: {out.relative_to(config.ROOT)}")


if __name__ == "__main__":
    main()
