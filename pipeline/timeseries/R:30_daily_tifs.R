# R/30_daily_tifs_parquet.R

write_daily_evi_tifs_from_parquet <- function(
    site_id,
    year,
    radius_m,
    geojson_path,
    site_dir = file.path("output", site_id),
    nc_path,
    geom_dir,       # e.g., output/AT-Nsd/tables_sf_fetch_600m/pixels_geom_ds
    ts_dir,         # e.g., output/AT-Nsd/tables_sf_fetch_600m/pixels_timeseries_ds
    meta_path = NULL,   # optional: parquet with key/value including crs_wkt
    crs_wkt = NULL,     # optional override
    thin = 3,
    overwrite = TRUE
) {
  stopifnot(dir.exists(geom_dir), dir.exists(ts_dir), file.exists(nc_path))
  
  requireNamespace("sf", quietly = TRUE)
  requireNamespace("terra", quietly = TRUE)
  requireNamespace("arrow", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  
  # ---- CRS (from meta if not provided) ----
  if (is.null(crs_wkt) || !nzchar(crs_wkt)) {
    if (!is.null(meta_path) && file.exists(meta_path)) {
      meta <- arrow::read_parquet(meta_path)
      crs_wkt <- meta$value[match("crs_wkt", meta$key)]
    }
    if (is.null(crs_wkt) || is.na(crs_wkt) || !nzchar(crs_wkt)) {
      stop("Need CRS WKT via crs_wkt= or meta_path with key=crs_wkt.")
    }
  }
  
  # ---- Buffer (in target CRS) ----
  aoi <- sf::st_read(geojson_path, quiet = TRUE)
  site_pt <- sf::st_centroid(sf::st_union(aoi))
  site_pt <- sf::st_transform(site_pt, sf::st_crs(crs_wkt))
  buf <- sf::st_buffer(site_pt, dist = radius_m)
  buf_v <- terra::vect(buf)
  
  # ---- Template raster (CROP FIRST; never * NA on full raster) ----
  r0 <- terra::rast(nc_path, lyrs = 1)
  if (is.na(terra::crs(r0)) || !nzchar(terra::crs(r0))) terra::crs(r0) <- crs_wkt
  r_crop <- terra::crop(r0, buf_v, snap = "out")
  tmpl <- terra::rast(r_crop)
  terra::values(tmpl) <- NA_real_
  
  # ---- Geometry: read once (small) ----
  geom_ds <- arrow::open_dataset(geom_dir, format = "parquet")
  geom_tbl <- geom_ds |>
    dplyr::select(pixel_id, x, y) |>
    dplyr::collect()
  
  if (!all(c("pixel_id","x","y") %in% names(geom_tbl))) {
    stop("geom_dir must contain columns: pixel_id, x, y")
  }
  
  # Build a SpatVector once (no sf needed)
  geom_tbl$pixel_id <- as.integer(geom_tbl$pixel_id)
  geom_tbl$x <- as.numeric(geom_tbl$x)
  geom_tbl$y <- as.numeric(geom_tbl$y)
  
  sv_geom <- terra::vect(
    geom_tbl,
    geom = c("x","y"),
    crs = crs_wkt,
    keepgeom = TRUE
  )
  
  # ---- Timeseries dataset (lazy) ----
  ts_ds <- arrow::open_dataset(ts_dir, format = "parquet")

  
  # ---- Output folder ----
  out_dir <- file.path(site_dir, paste0("evi_daily_", year, "_r", radius_m, "m"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  days <- seq(as.Date(paste0(year, "-01-01")), as.Date(paste0(year, "-12-31")), by = "day")
  days_use <- days[seq(1, length(days), by = thin)]
  
  for (i in seq_along(days_use)) {
    day <- days_use[i]
    f_out <- file.path(out_dir, paste0("EVI_", format(day, "%Y%m%d"), ".tif"))
    if (file.exists(f_out) && !overwrite) next
    
    # Pull ONLY this day's spline values
    ts_day <- ts_ds |>
      dplyr::filter(series == "spline", year == !!year) |>
      dplyr::select(pixel_id, date, evi) |>
      dplyr::collect() |>
      dplyr::mutate(date = as.Date(date)) |>
      dplyr::filter(date == day) |>
      dplyr::select(pixel_id, evi)
    
    if (nrow(ts_day) == 0) {
      message("No values for ", day, " (skipping)")
      next
    }
    
    ts_day$pixel_id <- as.integer(ts_day$pixel_id)
    ts_day$evi <- as.numeric(ts_day$evi)
    
    # Join onto geometry table via pixel_id (match is fast + low memory)
    idx <- match(sv_geom$pixel_id, ts_day$pixel_id)
    vals <- ts_day$evi[idx]  # aligned to sv_geom order; NAs where missing
    
    if (all(is.na(vals))) {
      message("All NA after join for ", day, " (skipping)")
      next
    }
    
    sv <- sv_geom
    sv$evi <- vals
    
    r <- terra::rasterize(sv, tmpl, field = "evi", touches = FALSE)
    r <- terra::mask(r, buf_v)
    
    if (terra::global(!is.na(r), "sum", na.rm = TRUE)[1,1] == 0) {
      message("All NA after rasterize/mask for ", day, " (skipping)")
      next
    }
    
    terra::writeRaster(r, f_out, overwrite = TRUE)
    
    if (i %% 25 == 0) message("Wrote ", i, "/", length(days_use))
  }
  
  invisible(list(out_dir = out_dir))
}