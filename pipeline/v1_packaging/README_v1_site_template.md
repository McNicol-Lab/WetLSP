# WetLSP v1.0 site outputs

Per-site outputs of the WetLSP (Wetland Land Surface Phenology) pipeline, version 1.0, derived from PlanetScope surface-reflectance imagery.
Source PlanetScope imagery itself is not redistributed here.

## Provenance

- Code: https://github.com/McNicol-Lab/WetLSP, branch `qc-outlier-fixes` (PR #2), commit `36817b8`; `pipelineVersion` "1.0".
  Includes c2a902a (GDAL NUM_THREADS=1 worker-write fix) and 36817b8 (per-scene UDM2 matching, `udm2MissingAction` = "drop").
- Steps 01–05 run from the raw PlanetScope scenes (`rawImage-unzip/<SITE>/`) with the branch's `pipeline/wetlsp-parameters.json`
  (only paths changed; `writeCSV` false). Vegetation index: EVI.
- Raw scenes were read through a flat symlink view of the local raw copy (`raw_flat/<SITE>/<scene files>`), because step 01 takes the
  acquisition date from the first 8-digit run in the relative path; with the Planet order layout `data/<order-uuid>/PSScene/...`,
  order UUIDs containing 8 digits would give wrong/NA dates. Outputs are unaffected by this workaround.
- Step 03 ran with a harness copy that only changes the number of worker cores.
- Processing location/dates: see `history`/`comment` attributes in each NetCDF. TODO: final processing summary and DOI.

## Files

- `<SITE>-wetlsp-<YEAR>.nc` (2021–2024): annual phenometrics written by step 05 and annotated with
  `pipeline/netcdf-metadata-populate.R` + `pipeline/nc-patch.R`; `product_version` = "v1.0".
  Same variables, units and fill values as the published v0.1 files.
- `<SITE_>_pixels_{geom,meta,timeseries}.parquet`: fetch circle at the published radius
  (`final_fetch_r` in `wetlsp_cyverse_site_catalog_final.csv`), built with `pipeline/extract_fetch_pixels.R`.
  For extra sites not in the published package the radius is set as described under "Fetch radius for extra sites" below;
  the radius used is recorded in the meta table key `radius_m`.
  Same layout as the v0.1 `README_parquet.md`, with these differences:
  * new column `qc_raw` (int8; null for `series = 'spline'`): v1.0 per-observation QC code for raw EVI
    (`QC_CODES` in `pipeline/wetlsp-functions.R`: 0 ok, 1 invalid (<=0, fill, EVI outside [-0.2, 1]), 2 despike,
    3 second-spline QC, 4 in-season below dormant floor, 5 floored off-season);
  * raw rows whose EVI was set to NA by v1.0 are kept with `evi = null` so their `qc_raw` is visible;
  * meta table has an extra key `product_version` = "v1.0".
- `compare_*.png`, `compare_summary.json` (package sites only): comparison against the published v0.1 package.

## Fetch radius for extra sites

Extra sites that lack a published fetch radius were assigned 300 m if the site latitude is extra-tropical
(|lat| ≥ 23.436°) and 600 m if tropical (|lat| < 23.436°). Sites with a published or catalog `final_fetch_r` used that value
instead. These radii define the fetch circle used to extract the pixel time-series Parquet files. Extras have no
published-package reference, so v0.1-vs-v1.0 comparison figures (`compare_*`) were not generated for them.
Latitudes are from `wetlsp_cyverse_site_catalog_final.csv` or `metadata/wetlsp-site-metadata.csv` (decision 2026-11-09).

- Catalog/metadata `final_fetch_r`: ES-Pdu 500, SE-Deg 300, RU-Ch2 500, RU-Fy2 500, RU-Vrk 500, MY-SBW 1600, PE-QFR 2200,
  PE-TNR 2200, US-NC4 2200, NA-Mym 600 (metadata row BD-Mym), NA-Amp 3000 (PE-Amp), NA-Los 3700 (PE-AMG),
  NA-Mam 1300 (BR-Mam).
- Assumed 600 m (tropical): HK-MPM, ID-DTP, ID-ITP.
- Assumed 300 m (extra-tropical): FI-Lom, FI-Si2, US-A03, US-A10, US-Beo, US-Bes, US-HRA, US-NGB, US-Tw4, US-xBA.

## Site ID aliases (download ID → final ID)

Four extra sites were processed from Planet/download folders that still use older IDs.
Per-site outputs are written under the final site IDs used in `wetlsp_cyverse_site_catalog_final.csv`:

- `NA-Amp` → `PE-Amp` (Ampiyacu, Peru)
- `NA-Mam` → `BR-Mam` (Mamirauá, Brazil)
- `NA-Mym` → `BD-Mym` (Mymensingh, Bangladesh)
- `NA-Los` → `PE-AMG` (Los Amigos, Peru)

Raw imagery remains under `rawImage-unzip/<download ID>/`. The catalog uses the final IDs and paths.

## Contact

Gavin McNicol, University of Illinois Chicago (gmcnicol@uic.edu). TODO: license, citation and acknowledgements text.
