# WetLSP v1.0

**Land surface phenology for wetland flux-tower sites from PlanetScope imagery**

WetLSP turns PlanetScope surface reflectance into pixel-wise daily EVI (observed and gap-filled) and annual phenometrics around eddy-covariance tower sites. It was built for wetland methane research, starting with FLUXNET-CH4 sites. It is a fork of the PLSP algorithms of Moon et al. (2022), extended to wetlands and with new outlier and quality-control steps.

The data were developed for, and first shared with, the McNicol ecoϕlab (UIC), the ESIIL AI for Natural Methane Working Group and the USGS Powell Center wetlands synthesis.

## What's new in v1.0
Compared with WetLSP v0.1 (`c50f58c`):
- **Invalid-value screening.** Band values ≤ 0, non-finite EVI and EVI outside [-0.2, 1] are set to NA before despiking.
- **Corrected spike test.** The MAD in `madTest` is now computed on the double difference, not on the EVI values. This bug is inherited from PLSP. `MADspikeThresh` goes from 2 to 4.
- **Second-spline dip filter** (`SecondSplineQC`). Removes short dips that return to baseline, while keeping real sustained drops such as mowing, flooding or fire.
- **Season-aware dormant floor.** Low values are raised to the dormant floor only outside the growing season; inside it they become NA.
- **Upper-envelope spline weighting.** Points below the fitted curve are down-weighted and the spline is refit.
- **Per-observation QC flags** (`qc_raw`): 0 ok, 1 invalid, 2 despike, 3 second-spline, 4 in-season below floor, 5 floored.
- **UTM-zone boundary fix** for sites near a zone edge.
- **Version tagging.** NetCDF files carry `product_version = "v1.0"`.

Full details, parameter defaults and a step-by-step comparison with PLSP are in [pipeline/QC_CHANGES.md](pipeline/QC_CHANGES.md).

## Data products
The v1.0 package covers 108 sites: the 82 sites of the original (v0.1) package plus 26 additional sites. For each site it contains:
- **Phenometrics:** one NetCDF per site per year.
- **Time series:** daily pixel-wise EVI, observed and spline-filled, with QC flags, as parquet clipped to the tower fetch area.
- **Geometry and metadata tables,** plus per-site `README_v1.md` and `README_parquet.md/.json`.

The site list and fetch radii are in [pipeline/v1_packaging/catalog](pipeline/v1_packaging/catalog).

**Archive:** ESS-DIVE, DOI to be added.

## Repository layout
```
pipeline/
  code/                 00–05 processing steps (download → NetCDF)
  wetlsp-functions.R    core phenology, QC and spline functions
  wetlsp-parameters.json  run settings, QC and phenology parameters
  run_all.R             runs steps 01–05 for one or more sites
  write_chunks_to_parquet_ds.R, extract_fetch_pixels.R, combine_fetch_sf*.R
                        fetch-area time-series extraction
  v1_packaging/         scripts, catalog and README generators used to build the v1.0 package
  tests/                synthetic QC test, UDM2 matching test
  QC_CHANGES.md         v1.0 methods and differences from PLSP
  validation/, timeseries/  exploratory evaluation and visualization (not required)
workflows/              R Markdown review and site-selection notebooks
geojson/AMFLX/          site AOI polygons
metadata/               site lists and product-layer table
```

## Requirements
- R ≥ 4 with `terra`, `raster`, `sf`, `ncdf4`, `foreach`, `doParallel`, `rjson`, `jsonlite`, `geojsonR` and `gdalUtilities`, plus `arrow` for parquet output.
- GDAL, PROJ, GEOS and udunits system libraries.
- Python 3 with `requests` and `pandas` for the image download.
- A Planet account with API access.

A `Dockerfile` based on the CyVerse geospatial RStudio image is included.

## Running the pipeline
1. **Clone** the repo and work from its root; every script reads `pipeline/wetlsp-parameters.json` by relative path.
2. **Configure** the `setup` block of `pipeline/wetlsp-parameters.json`. Set `workDir`, `geojsonDir`, `dataDir` and `outDir`, the year range (`phenStartYr`, `phenEndYr`), and `numChunks` and `numCores`.
3. **Add site polygons** as GeoJSON in `geojson/AMFLX/`. The pipeline expects that subfolder.
4. **Download imagery.** Set your Planet API key in `pipeline/code/00_img_download.py`, then run it to order and download PlanetScope surface reflectance and UDM2 into `dataDir`.
5. **Process** one or more sites by site index:
   ```bash
   Rscript pipeline/run_all.R 12          # one site
   Rscript pipeline/run_all.R 12 13 14    # several sites
   ```
   `run_all.R` runs these steps in order and writes logs to `logs/site_<N>/`:
   - `01-wetlsp-process.R`: clip, mask and mosaic the images
   - `02-wetlsp-chunks.R`: split the scene into pixel chunks
   - `03-wetlsp-phenometrics.R`: despike, QC, spline and phenometrics
   - `04-wetlsp-generate.R`: assemble the chunks
   - `05-wetlsp-format.R`: write the annual NetCDF

   Each step can also be run alone, for example `Rscript pipeline/code/03-wetlsp-phenometrics.R 12`.
6. **Extract the fetch-area time series** to parquet with `pipeline/write_chunks_to_parquet_ds.R` and the `extract_fetch_pixels.R` / `combine_fetch_sf*.R` helpers.

### Reproducing the v1.0 package
`pipeline/v1_packaging/` holds the scripts that built the original (v0.1) v1.0 site products from steps 01–05: `post.R`, `post_nc.R`, `finalize.py`, `compare.py`, `mkparams.py` and the README generators. They keep the paths of the original processing machine and are a record of what was run, not yet a general tool. See [pipeline/v1_packaging/README.md](pipeline/v1_packaging/README.md).

### Testing
```bash
Rscript pipeline/tests/test_qc_synthetic.R
```
Builds a synthetic 2023–2025 series with injected dips and zero-valued observations, then checks that the v1.0 QC removes them without biasing the summer spline.

## Downstream workflows (`workflows/`)
- `nc-check-populate.Rmd`: standardize and populate the per-site, per-year phenometric NetCDFs.
- `core-site-selection.Rmd`: find sites with enough coverage for analysis.
- `review-final-wetlsp.Rmd`: inspect phenometrics by site-year and the observed and gap-filled EVI.
- `wetlsp-timeseries.Rmd`: build the fetch perimeter and extract pixel-wise EVI to `.rds`.
- `review-fetch-rds.Rmd`: inspect the fetch perimeter and explore the `.rds` output.

## Known limitations
- PLSP's open-water mask is disabled. WetLSP uses EVI and UDM2 only, not EVI2 with the legacy UDM.
- Stricter UDM2 masking (requiring clear cloud, shadow, haze and snow bands) is not yet implemented because it requires rerunning the image step.
- Hard-coded paths remain in some scripts. Check `wetlsp-parameters.json` and the `v1_packaging` scripts before running elsewhere.

## Versions
- **WetLSP v0.1:** the original fork (`c50f58c`).
- **WetLSP v1.0:** QC update, zone fix and packaging record (this release).

## Citation
If you use WetLSP, please cite the data package (ESS-DIVE DOI to be added) and the original method:

Moon, M., Richardson, A.D., Milliman, T. and Friedl, M.A. (2022). A high spatial resolution land surface phenology dataset for AmeriFlux and NEON sites. *Scientific Data*, 9, 448. https://doi.org/10.1038/s41597-022-01570-5

Moon, M., Richardson, A.D., Milliman, T. and Friedl, M.A. (2023). Land Surface Phenology, Eddy Covariance Tower Sites, North America, 2017–2021. ORNL DAAC. https://doi.org/10.3334/ORNLDAAC/2033

## Acknowledgments
Forked from the BU-LCSC PLSP algorithms (Moon et al., 2022). The image download script is by Dennis Milechin, adapted by Minkyu Moon. Developed in the McNicol ecoϕlab at the University of Illinois Chicago.

## License
To be added.
