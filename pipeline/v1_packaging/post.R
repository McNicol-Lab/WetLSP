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
# Streaming fetch filter: same logic as extract_fetch_pixels() (make_site_buffer + filter_chunk_by_buffer),
# but holds only one chunk in memory (the repo function accumulates kept_list -> 24 GB limit at 2000 m).
extract_fetch_stream <- function(site_dir, site_id, geojson_path, radius_m, in_subdir = "tables_sf", out_subdir, ...) {
  in_dir <- file.path(site_dir, in_subdir); out_dir <- file.path(site_dir, out_subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cf <- list.files(in_dir, pattern = "^chunk_\\d{3}_evi_sf\\.rds$", full.names = TRUE)
  if (!length(cf)) stop("No chunk sf RDS files found in: ", in_dir)
  f1 <- readRDS(cf[1]); buf <- make_site_buffer(geojson_path, st_crs(f1), radius_m); rm(f1)
  for (i in seq_along(cf)) {
    f_out <- file.path(out_dir, basename(cf[i])); if (file.exists(f_out)) unlink(f_out)
    k <- filter_chunk_by_buffer(readRDS(cf[i]), buf)
    if (nrow(k) > 0) saveRDS(k, f_out) else message("No pixels intersect buffer in ", basename(cf[i]))
    rm(k); if (i %% 10 == 0) gc(FALSE)
    if (i %% 25 == 0) message("Processed ", i, "/", length(cf))
  }
  invisible(TRUE)
}
if (!nzchar(Sys.getenv("POST_TEST"))) res <- extract_fetch_stream(site_dir = site_dir, site_id = site, geojson_path = file.path(W, "geojson/AMFLX", paste0(site, ".geojson")),
                            radius_m = radius_m, out_subdir = paste0("tables_sf_fetch_", radius_m, "m"), overwrite = TRUE, return_objects = FALSE)
fdir <- file.path(site_dir, paste0("tables_sf_fetch_", radius_m, "m"))
files <- sort(list.files(fdir, pattern = "^chunk_\\d{3}_evi_sf\\.rds$", full.names = TRUE))
if (nzchar(Sys.getenv("POST_TEST"))) files <- sort(list.files(file.path(site_dir, "tables_sf"), pattern = "^chunk_\\d{3}_evi_sf\\.rds$", full.names = TRUE))[1:2]
if (!length(files)) stop("no fetch chunk files in ", fdir)
tmp <- file.path(R, "work", site, "pq_ds"); unlink(tmp, recursive = TRUE); dir.create(file.path(tmp, "ts"), recursive = TRUE); dir.create(file.path(tmp, "geom"))
crs_wkt <- NULL; npx <- 0
for (i in seq_along(files)) {
  s <- readRDS(files[i]); if (!nrow(s)) next
  if (is.null(crs_wkt)) crs_wkt <- st_crs(s)$wkt
  xy <- st_coordinates(s$geometry); pid <- as.integer(s$cell); npx <- npx + nrow(s)
  write_parquet(tibble::tibble(pixel_id = pid, cell = as.integer(s$cell), x = xy[,1], y = xy[,2]), file.path(tmp, "geom", sprintf("g_%03d.parquet", i)))
  d <- st_drop_geometry(s)
  sp <- tibble::tibble(pixel_id = rep(pid, lengths(d$evi_spline)), date = as.Date(unlist(d$dates_spline), origin = "1970-01-01"), evi = unlist(d$evi_spline))
  sp <- sp[!is.na(sp$evi), ]; sp$series <- "spline"; sp$qc_raw <- NA_integer_
  qc <- if ("qc_raw" %in% names(d)) d$qc_raw else vector("list", nrow(d))
  rw <- tibble::tibble(pixel_id = rep(pid, lengths(d$evi_raw)), date = as.Date(unlist(d$dates_raw), origin = "1970-01-01"), evi = unlist(d$evi_raw))
  qv <- unlist(qc); rw$qc_raw <- if (length(qv) == nrow(rw)) as.integer(qv) else NA_integer_
  if (length(qv) != nrow(rw)) message("WARN qc_raw length mismatch in ", basename(files[i]))
  rw$series <- "raw"   # raw rows kept even when evi is NA (v1.0 invalid -> NA) so qc_raw is visible
  ts <- bind_rows(sp, rw) |> transmute(pixel_id, series, date, year = as.integer(format(date, "%Y")), evi = as.numeric(evi), qc_raw = as.integer(qc_raw))
  write_parquet(ts, file.path(tmp, "ts", sprintf("t_%03d.parquet", i)))
  rm(s, d, sp, rw, ts, qc, qv, xy); gc(FALSE)
}
write_parquet(tibble::tibble(key = c("crs_wkt", "site_id", "radius_m", "pixel_id_source", "product_version"),
                             value = c(crs_wkt, site, as.character(radius_m), "cell", "v1.0")), file.path(tmp, "meta.parquet"))
message("POST_DONE ", site, " fetch pixels=", npx)
