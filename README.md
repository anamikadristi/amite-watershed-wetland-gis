# Land Cover, Wetlands and Water Quality in the Amite River Basin, Louisiana

Anamika Dristi, Ph.D.

## Summary

This project links land cover, mapped wetlands and stream water quality in
three HUC8 units that drain to Lake Maurepas in southeast Louisiana and
southwest Mississippi: the Amite (08070202), Tickfaw (08070203) and Lake
Maurepas (08070204) units. It uses only public data from USGS, MRLC, USFWS
and the Water Quality Portal. A Python pipeline downloads the data, summarizes
NLCD 2001 and 2021 land cover and NWI wetlands for 82 HUC12 subwatersheds,
delineates the upstream NHDPlus basin of 74 water-quality sites, and tests
how site median water quality relates to basin land cover.

Main results from the run in this repository:

* The study area covers 8,606.4 km2. In NLCD 2021 it is 33.15% wetland,
  29.77% forest, 13.68% developed and 12.46% agriculture.
* From 2001 to 2021, developed land rose from 12.05% to 13.68% of the area
  (+1.64 percentage points) and rose in all 82 HUC12 units. The largest gain
  was +12.99 points in Grand Goudine Bayou-New River (080702040103).
* NWI maps 2,060.7 km2 of vegetated wetland (23.94% of the study area),
  mostly freshwater forested/shrub wetland.
* Across 74 sites, median specific conductance rose with developed cover
  (Spearman rho = +0.85) and fell with forest cover (rho = -0.83). Total
  phosphorus rose with developed cover (rho = +0.81, n = 33 sites).
  Wetland cover showed no clear relationship with dissolved oxygen, nitrate
  plus nitrite or total phosphorus.

These are cross-sectional associations among a limited set of sites, not
causal estimates. See Key findings and Limitations.

## Independence statement

This is an independent analysis of public data. It is separate from the
author's LSU dissertation research. No dissertation field measurements or
laboratory data are used, copied or referenced in this repository.

## Related published work (context only)

The author's peer-reviewed work on this region gives scientific context. None
of its data are used here.

* Dristi, A., and Xu, Y. J. (2023). Water, 15(14), 2621.
* Dristi, A., and Xu, Y. J. (2025). Water, 17(17), 2589.
* Dristi, A., et al. (2025). Estuarine, Coastal and Shelf Science, 324, 109465.

## Data sources

All data were accessed on 2026-10-23. The exact request URLs, file sizes and
access dates are listed in `data/download_manifest.csv`.

| Data set | Provider | Access point | Use |
|---|---|---|---|
| Watershed Boundary Dataset, HUC8 and HUC12 | USGS | [The National Map WBD service](https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer) | Study units, HUC12 routing (tohuc) |
| NLCD 2001 and 2021 land cover, 30 m | MRLC / USGS | [MRLC](https://www.mrlc.gov/data), [GeoServer WCS](https://www.mrlc.gov/geoserver/mrlc_download/wcs) | Land-cover composition and change |
| National Wetlands Inventory, HU8 downloads | USFWS | [NWI](https://www.fws.gov/program/national-wetlands-inventory/data-download), [HU8 files](https://documentst.ecosphere.fws.gov/wetlands/downloads/watershed/) | Mapped wetland area |
| NHDPlus V2 flowlines, waterbodies, catchments | USGS / U.S. EPA | [USGS WMA GeoServer](https://api.water.usgs.gov/geoserver/wmadata/ows) | Streams, site basin delineation |
| Discrete water-quality results and stations | Water Quality Portal (LDEQ, USGS, EPA NARS and others) | [waterqualitydata.us](https://www.waterqualitydata.us/) | Site medians |
| Daily mean discharge (00060) | USGS NWIS | [NWIS Water Services](https://waterservices.usgs.gov/) | Flow context |
| Upstream basins for 5 sites (check only) | USGS NLDI | [NLDI](https://api.water.usgs.gov/nldi/) | Independent check of the local delineation |

Water Quality Portal query: stream sites in the three HUC8 units, 2000-01-01
to 2025-12-31, characteristics: organic carbon, specific conductance,
dissolved oxygen, pH, water temperature, nitrate, nitrate plus nitrite,
inorganic nitrogen (nitrate and nitrite), phosphorus, total phosphorus mixed
forms, and turbidity.

## Methods

1. **Study units** (`02_prepare_boundaries.py`). WBD HUC12 polygons for the
   three HUC8 units (82 units). Areas are computed in NAD83 CONUS Albers
   (EPSG:5070). The WBD `tohuc` field is used to list every HUC12 upstream of
   each unit.
2. **Land cover** (`03_landcover_nlcd.py`). NLCD 2001 and 2021 rasters are
   clipped on the server to the study-area box. HUC12 polygons are
   rasterized onto the 30 m grid and cells are counted by class. Classes are
   grouped as water (11), developed (21 to 24), barren (31), forest (41 to
   43), shrub/grass (52, 71), agriculture (81, 82) and wetland (90, 95).
   Outputs include per-HUC12 percentages, upstream totals, 2001 to 2021
   change and a group transition matrix.
3. **Wetlands** (`04_wetlands_nwi.py`). NWI polygons are intersected with
   HUC12 polygons. "Vegetated wetland" is the sum of freshwater
   forested/shrub, freshwater emergent, and estuarine and marine wetland.
   Ponds, lakes, riverine and deepwater classes are reported separately.
4. **Water quality** (`05_water_quality.py`). Quality-control samples are
   removed, units are harmonized, non-detects are set to half the reported
   detection limit (348 of 15,917 cleaned values), and LDEQ sites listed under
   two organization codes are merged. Specific conductance values labelled
   mS/cm but in the fresh-water uS/cm range (all from one LDEQ project) are
   treated as a unit label error and kept as uS/cm. Site medians use samples
   from 2011 onward, the years closest to NLCD 2021, and need at least 5
   samples.
5. **Point to watershed linking** (`06_site_basins.py`). Each site is placed
   in the NHDPlus V2 catchment that contains it, and the network is walked
   upstream with the fromnode/tonode fields. The upstream catchments form the
   site basin, from which NLCD and NWI percentages are computed. Each site is
   also tagged with the HUC12 that contains it.
6. **Relationships** (`07_relationships.py`). For each parameter and
   predictor (NLCD wetland, NWI vegetated wetland, forest, developed and
   agriculture percentages, and log10 basin area): Spearman rank correlation
   and a simple least-squares regression (log10 of the median for
   conductance, nutrients and turbidity). Benjamini-Hochberg q-values are
   computed across the 42 tests in each site subset. Two sensitivity subsets
   drop the Lake Maurepas HUC8 sites (tidal and backwater influence) or the
   LDEQ nonpoint-source (319) project sites (spatially clustered).
7. **Flow context** (`08_discharge_context.py`). Amite River main-stem
   samples within 5 km of USGS 07378500 (Denham Springs) or 07380120 (Port
   Vincent) are paired with same-day mean discharge for log-log
   concentration-discharge fits.
8. **Maps and figures** (`09_maps.py`, plus figures from scripts 04, 07, 08).

## Repository structure

```
amite-watershed-wetland-gis/
  README.md  LICENSE  requirements.txt  .gitignore
  scripts/
    config.py                 paths, CRS, study codes, thresholds
    common.py                 shared helpers (site IDs, NWI reader, FDR)
    01_download_data.py       downloads all public inputs to data_raw/
    02_prepare_boundaries.py  HUC8/HUC12, streams, HUC12 upstream lists
    03_landcover_nlcd.py      NLCD 2001/2021 per HUC12, change, transitions
    04_wetlands_nwi.py        NWI wetland area per HUC12, NWI map
    05_water_quality.py       WQP cleaning, site medians
    06_site_basins.py         NHDPlus upstream basins, basin land cover
    07_relationships.py       Spearman and regression tests, figures
    08_discharge_context.py   discharge and concentration-discharge
    09_maps.py                maps and summary figures
    run_all.py                runs scripts 02 to 09
  data/                       small derived public files
    amite_maurepas_layers.gpkg   huc8, huc12, flowlines, waterbodies, wq_sites, site_basins
    wq_results_clean.csv         cleaned WQP results (2000 to 2025)
    wq_sites.csv                 site list with HUC12 assignment
    huc12_upstream_members.csv   HUC12 upstream lists
    download_manifest.csv        URLs, sizes and access dates
  outputs/
    figures/                  fig01 to fig10 (PNG)
    tables/                   CSV results
```

Raw downloads (about 120 MB of NLCD tiles, NWI zip files, NHDPlus pages and
WQP files) go to `data_raw/`, which is not tracked by Git.

## How to run

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python scripts/01_download_data.py         # about 5 minutes; add --nldi for the NLDI check basins (slow)
python scripts/run_all.py                  # about 1 minute
```

Set the environment variable `AMITE_DATA_RAW` to keep downloads in another
folder. No API keys or accounts are needed. Web services change over time, so
a later download can give slightly different numbers.

## Key findings

All numbers come from the tables in `outputs/tables/`.

### Land cover and wetlands (HUC12 level)

| Group | NLCD 2001 (%) | NLCD 2021 (%) | Change (points) |
|---|---|---|---|
| Wetland | 33.40 | 33.15 | -0.25 |
| Forest | 25.73 | 29.77 | +4.04 |
| Developed | 12.05 | 13.68 | +1.64 |
| Agriculture | 15.04 | 12.46 | -2.58 |
| Shrub/grass | 9.71 | 6.79 | -2.92 |
| Water | 3.64 | 3.68 | +0.04 |
| Barren | 0.43 | 0.46 | +0.03 |

* HUC12 medians in 2021: wetland 28.6%, forest 31.4%, developed 5.6%.
  The Lake Maurepas unit is mostly wetland (HUC12 median 68.7%), while the
  Amite and Tickfaw units are mostly forest (medians 44.9% and 37.5%).
* 14.92% of cells changed group between 2001 and 2021. Much of this is
  exchange between forest and shrub/grass (331.4 km2 forest to shrub/grass,
  573.7 km2 the other way), which is consistent with timber harvest and
  regrowth.
  The forest gain should be read with that in mind.
* About 29.5 km2 of 2001 wetland and 54.9 km2 of 2001 agriculture were
  developed by 2021. The largest developed gains were in Grand Goudine
  Bayou-New River (33.06% to 46.05%), Bayou Fountain-Bayou Manchac (37.02% to
  46.96%) and West Colyell Creek (14.45% to 22.04%), all in the western
  part of the study area near Baton Rouge.
* NWI vegetated wetland covers 23.94% of the study area (forested/shrub
  1,944.7 km2, emergent 68.0 km2, estuarine 47.9 km2). This is lower than
  NLCD wetland (33.15%) because the two products use different mapping
  methods, dates and definitions. HUC12 NWI vegetated wetland ranges from
  0.7% to 97.4% (median 14.0%).

### Site basins

* 74 sites met the sample rules (32 LDEQ ambient network, 29 LDEQ
  nonpoint-source project, 13 USGS Mississippi sites) in 39 HUC12 units.
  Basin areas range from 1.09 to 8,543.4 km2 (median 131.2 km2).
* The local NHDPlus delineation matched USGS NLDI basins for 5 check sites
  (overlap, intersection over union, 0.990 to 0.998) and was within 7% of
  the drainage areas listed in WQP metadata for 3 USGS sites (471.5 vs 463.6,
  244.0 vs 261.6 and 73.8 vs 74.6 km2).

### Water quality and basin land cover (Spearman rho, all sites)

| Parameter | Sites | Developed % | Forest % | NLCD wetland % | NWI veg. wetland % |
|---|---|---|---|---|---|
| Specific conductance | 74 | +0.85* | -0.83* | +0.45* | +0.27 |
| Total phosphorus | 33 | +0.81* | -0.65* | +0.13 | -0.08 |
| Nitrate + nitrite N | 32 | +0.48* | -0.25 | -0.19 | -0.34 |
| Turbidity | 46 | +0.54* | -0.48* | +0.51* | +0.28 |
| Dissolved oxygen | 61 | -0.36* | +0.42* | -0.15 | -0.06 |
| pH | 74 | +0.45* | -0.46* | +0.31* | +0.31* |

\* Benjamini-Hochberg q < 0.05 across 42 tests. NWI columns have fewer sites
(basins outside the NWI download extent are left blank). Full results,
including regression slopes and R2, are in
`outputs/tables/wq_landcover_tests.csv`.

* The clearest pattern is a gradient from forested headwater basins (low
  conductance, low phosphorus, higher dissolved oxygen) to developed basins
  in the western part of the basin near Baton Rouge (higher conductance,
  phosphorus and nitrate plus nitrite). This pattern held when the Lake
  Maurepas sites were dropped (64 sites; conductance vs developed rho =
  +0.88; phosphorus vs developed rho = +0.85, n = 23) and when the LDEQ 319
  project sites were dropped (45 sites; conductance vs forest rho = -0.88).
* Total phosphorus and nitrate plus nitrite were lower at sites with larger
  basins (rho = -0.77 and -0.59), which is consistent with dilution and
  in-stream processing along the main stems but is also confounded with
  where the urban sites are.
* Wetland cover was not a clear predictor of nutrients or dissolved oxygen
  at this scale. The positive links between NLCD wetland and conductance or
  turbidity probably reflect location (low-gradient lower-basin sites), not
  a wetland effect.
* Forest and developed cover are strongly related across site basins
  (Spearman rho = -0.71 between them), so their separate effects cannot be
  told apart with these data.

### Flow context

* Daily mean discharge at USGS 07378500 (2000 to 2025) had a median of 798
  ft3/s and a maximum of 202,000 ft3/s.
* At Denham Springs (LDEQ site, 2017 to 2021, about 21 to 23 samples),
  conductance fell with flow (log-log slope -0.25, R2 0.52), and turbidity
  (slope +0.83, R2 0.60) and total phosphorus (slope +0.43, R2 0.31) rose
  with flow. At Port Vincent (2001 to 2014, 40 to 45 samples), turbidity rose
  with flow (slope +0.46, R2 0.51) and the conductance and phosphorus fits
  were weak (R2 0.15 and 0.00).
* The median flow on sampling days sat at the 48th to 55th percentile of
  daily flows in the same period, so routine samples mostly describe
  ordinary flow, not floods.

## Limitations

* **Small and uneven samples.** Phosphorus and nitrate plus nitrite have
  only 33 and 32 sites. Sites come from different programs, periods and
  sampling seasons, and medians do not adjust for season or flow.
* **Not independent.** Many basins are nested along the same rivers, and the
  LDEQ 319 sites are clustered. P-values therefore overstate the evidence.
  No spatial or mixed model was fitted.
* **Confounding.** Land-cover groups are correlated with each other and with
  position in the basin, soils and wastewater inputs, none of which are
  modelled here.
* **Dissolved organic carbon.** Only 43 DOC values at 13 sites were found
  for 2000 onward, mostly from 2000, and no site met the sample rule, so DOC
  could not be analysed.
* **USGS gauges.** WQP discrete samples at 07378500, 07377500 and 07380120
  for the parameters checked predate 2000, so these gauges serve as flow
  references only. Discharge at 07377500 and 07380120 ends in 2017 and 2015.
* **Data currency.** The legacy WQP service used here does not include USGS
  discrete data added after March 2024. LDEQ data run to 2025.
* **Land cover.** NLCD has classification error, and its wetland classes
  differ from NWI. NWI mapping dates vary across the area. The NLCD 2001 layer was taken from
  the MRLC download service and its release version was not checked
  separately. Comparisons between years are indicative only.
* **Unit and coding issues.** Unit harmonization follows the rules in
  `05_water_quality.py`; turbidity units (NTU, FNU, NTRU) are treated as
  comparable; one LDEQ conductance unit label was corrected by value range.
* **Tidal sites.** Pass Manchac and other Lake Maurepas sites are affected by
  lake and estuary exchange, which upstream land cover does not describe.

## Contact

Anamika Dristi, Ph.D. (anamikadristi@gmail.com)
