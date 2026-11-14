#!/usr/bin/env Rscript
# Post-pipeline packaging for the v1.0 test: annotated NetCDFs (workflows/nc-check-populate.Rmd functions)
# and fetch-circle time series (pipeline/extract_fetch_pixels.R, as workflows/wetlsp-timeseries.Rmd), with qc_raw kept.
suppressPackageStartupMessages({library(sf); library(dplyr); library(tidyr); library(arrow); library(ncdf4); library(terra); library(purrr); library(stringr)})
a <- commandArgs(TRUE); site <- a[1]; radius_m <- as.numeric(a[2])
R <- path.expand("~/wetlsp_v1_test2"); C <- file.path(R, "code"); W <- file.path(R, "work", site)
source(file.path(C, "pipeline/netcdf-metadata-populate.R")); source(file.path(C, "pipeline/nc-patch.R"))
source(file.path(C, "pipeline/extract_fetch_pixels.R"))
site_dir <- file.path(W, "output", site); out <- file.path(R, "final", site); dir.create(out, recursive = TRUE, showWarnings = FALSE)
years <- 2021:2024
for (yr in years) {
  nc_in <- file.path(W, "output", "Product_netCDF", site, paste0("WetLSP_", yr, ".nc"))
  if (!file.exists(nc_in)) { message("MISSING nc ", yr); next }
  nc_out <- file.path(out, paste0(site, "-wetlsp-", yr, ".nc"))
  annotate_wetlsp_netcdf(nc_path = nc_in, out_path = nc_out, site_id = site, year = yr, creator_email = "gmcnicol@uic.edu",
                         pipeline_version = "v1.0", years_included = paste(range(years), collapse = "–"))
  patch_cf_coords_and_grid_mapping(nc_out)
  nc <- nc_open(nc_out, write = TRUE); ncatt_put(nc, 0, "product_version", "v1.0")
  ncatt_put(nc, 0, "comment", paste0("WetLSP v1.0, pipeline commit 36817b8 (branch qc-outlier-fixes), processed ", format(Sys.Date()), ". This file contains phenometrics only (no time dimension). Companion time-series artifacts are distributed separately."))
  nc_close(nc)
}
