library(sf)
library(terra)

# ---- Resolve NetCDF path for a given site/year ----
resolve_wetlsp_nc <- function(site_dir, site_id, year, nc_subdir = "netcdf") {
  nc_dir <- file.path(site_dir, nc_subdir)
  
  # List .nc files with full paths
  nc_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
  
  if (length(nc_files) == 0) {
    stop("No .nc files found in: ", nc_dir)
  }
  
  # Prefer files that contain BOTH site_id and year anywhere in name
  hits <- nc_files[
    grepl(site_id, basename(nc_files), fixed = TRUE) &
      grepl(as.character(year), basename(nc_files), fixed = TRUE)
  ]
  
  # If none match site_id, match by year only (common if naming varies)
  if (length(hits) == 0) {
    hits <- nc_files[grepl(as.character(year), basename(nc_files), fixed = TRUE)]
  }
  
  if (length(hits) == 0) {
    stop("No NetCDF matched year ", year, " in: ", nc_dir,
         "\nFiles seen:\n", paste(basename(nc_files), collapse = "\n"))
  }
  
  # If multiple matches, pick the most recent
  hits <- hits[order(file.info(hits)$mtime, decreasing = TRUE)]
  hits[[1]]
}

# ---- Your inputs ----
site <- "AT-Nsd"
year <- 2023
site_dir <- file.path("output", site)

nc_path <- resolve_wetlsp_nc(site_dir, site, year)
message("Using NetCDF: ", nc_path)

geojson_path <- "geojson/AMFLX/AT-Nsd.geojson"
radius_m <- 100

site_dir <- file.path("output/", site)
sf_all <- readRDS(file.path(site_dir, "AT-Nsd_evi_timeseries_fetch_sf_100m.rds"))

# ---- A) Ensure spline dates are Date ONCE (not inside make_day_raster) ----
sf_all$dates_spline <- lapply(sf_all$dates_spline, function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1970-01-01"))
  as.Date(x)
})

# ---- B) Build buffer in sf_all CRS ----
target_crs <- st_crs(sf_all)

aoi <- st_read(geojson_path, quiet = TRUE)
site_pt <- st_centroid(st_union(aoi))
site_pt <- st_transform(site_pt, target_crs)

buf   <- st_buffer(site_pt, dist = radius_m)
buf_v <- terra::vect(buf)

# ---- C) Template raster: read, assign CRS from sf_all (no reprojection), crop ----
# nc_path must exist in your environment
tmpl_full <- terra::rast(nc_path)[[1]]
tmpl_full <- tmpl_full * NA
terra::crs(tmpl_full) <- target_crs$wkt

tmpl <- terra::crop(tmpl_full, terra::ext(buf_v))
tmpl <- tmpl * NA

# Convert sf to SpatVector once
sv_geom <- terra::vect(sf_all)

# ---- D) Day raster function (robust match) ----
make_day_raster <- function(day) {
  day_i <- as.integer(day)
  
  vals <- mapply(
    FUN = function(dts, evi) {
      if (length(dts) == 0L || length(evi) == 0L) return(NA_real_)
      j <- match(day_i, as.integer(dts))   # robust Date/POSIX conversion
      if (is.na(j)) NA_real_ else as.numeric(evi[[j]])
    },
    dts = sf_all$dates_spline,
    evi = sf_all$evi_spline,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
  
  # Optional: quick debug—how many pixels had a value this day?
  n_hit <- sum(!is.na(vals))
  if (n_hit == 0) return(NULL)
  
  sv <- sv_geom
  sv$evi <- vals
  
  r <- terra::rasterize(sv, tmpl, field = "evi", touches = FALSE)
  r <- terra::mask(r, buf_v)
  
  r
}

# ---- E) Write outputs ----
out_dir <- file.path(site_dir, paste0("evi_daily_", year, "_r", radius_m, "m"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

thin <- 3
days <- seq(as.Date(paste0(year, "-01-01")), as.Date(paste0(year, "-12-31")), by = "day")
days_use <- days[seq(1, length(days), by = thin)]

for (i in seq_along(days_use)) {
  day <- days_use[i]
  f_out <- file.path(out_dir, paste0("EVI_", format(day, "%Y%m%d"), ".tif"))
  
  r <- make_day_raster(day)
  
  if (is.null(r)) {
    message("No matched pixels for ", day, " (skipping)")
    next
  }
  
  # Sanity: skip if completely empty after rasterize/mask
  if (terra::global(!is.na(r), "sum", na.rm = TRUE)[1,1] == 0) {
    message("All NA after rasterize/mask for ", day, " (skipping)")
    next
  }
  
  terra::writeRaster(r, f_out, overwrite = TRUE)
  
  if (i %% 25 == 0) message("Wrote ", i, "/", length(days_use))
}