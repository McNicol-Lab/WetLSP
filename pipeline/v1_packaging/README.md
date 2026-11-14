# WetLSP v1.0 packaging scripts (record of the v1.0 reprocess)

These scripts turned pipeline steps 01–05 output into the v1.0 site products (annual NetCDF +
`<SITE_>_pixels_{geom,meta,timeseries}.parquet` + READMEs) for the 108-site WetLSP v1.0 package,
November 2026. They are committed as an exact record of what was run. They contain hard-coded paths from the
processing Mac (`~/wetlsp_v1_test2/...`) and are not yet a general tool (see the release plan; a configurable
runner is future work).

| File | Role |
|---|---|
| `mkparams.py` | Per-site copy of `pipeline/wetlsp-parameters.json` (paths only; `writeCSV` false; core count) |
| `post.R` | Annotates step-05 NetCDFs (`netcdf-metadata-populate.R`, `nc-patch.R`, `product_version` v1.0) and extracts the fetch-circle pixel time series with `extract_fetch_pixels.R`, keeping `qc_raw` |
| `post_nc.R` | NetCDF-only variant of `post.R` |
| `finalize.py` | Consolidates per-site parquet into single sorted zstd files with the published write settings and geometry CRS metadata |
| `compare.py` | v1.0 vs published-package comparison figures and `compare_summary.json` (published-package sites only) |
| `make_extras_readmes.py`, `make_package_readmes.py` | Write v1.0 `README_parquet.md/.json` per site (extras / published-package sites) and remove `candidate_extras` wording from `README_v1.md` |
| `README_v1_site_template.md` | Per-site `README_v1.md` (provenance) as shipped |
| `extras_fetch_radii.csv` | Fetch radius and source for the 26 sites not in the published package |
| `catalog/wetlsp_cyverse_site_catalog_final.{csv,json}` | 108-site v1.0 catalog (copy of the NAS `wetlsp-v1` root) |

Orchestration (staging, retries, watchdog, verified NAS push, power schedule) lived in the processing harness and
is archived with its logs on the lab NAS (`wlsp/wetlsp-v1_harness_backup_20261113`).
