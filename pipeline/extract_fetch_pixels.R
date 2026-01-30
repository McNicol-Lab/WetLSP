library(sf)
library(terra)
library(purrr)
library(stringr)

# ---- helpers ----
### make template
make_template_from_netcdf_or_sf <- function(nc_path = NULL, sf_ref, buffer_geom = NULL, pad_m = 0) {
  stopifnot(inherits(sf_ref, "sf"))
  
  # 1) Start with NetCDF grid if provided
  if (!is.null(nc_path) && file.exists(nc_path)) {
    tmpl <- terra::rast(nc_path)[[1]]
    
    # If netcdf crs is empty, copy WKT from sf_ref
    if (is.na(terra::crs(tmpl)) || terra::crs(tmpl) == "") {
      terra::crs(tmpl) <- sf::st_crs(sf_ref)$wkt
    }
    
    tmpl <- tmpl * NA
  } else {
    # 2) Fallback: template from sf geometry (less ideal; but consistent CRS)
    sv <- terra::vect(sf_ref)
    e  <- terra::ext(sv)
    tmpl <- terra::rast(e, resolution = 3, crs = sf::st_crs(sf_ref)$wkt) * NA
  }
  
  # Optional zoom template around buffer
  tmpl_zoom <- NULL
  if (!is.null(buffer_geom)) {
    buf_sv <- terra::vect(buffer_geom)
    e_zoom <- terra::ext(buf_sv)
    if (pad_m > 0) e_zoom <- e_zoom + c(-pad_m, pad_m, -pad_m, pad_m)
    tmpl_zoom <- terra::crop(tmpl, e_zoom)
  }
  
  list(template = tmpl, template_zoom = tmpl_zoom)
}

# Create a radius buffer around a site point read from a GeoJSON.
# - geojson_path: e.g. "geojson/AMFLX/AT-Nsd.geojson"
# - target_crs: optional; if NULL, will use the CRS of the chunk sf objects
# - radius_m: buffer radius in meters (e.g., 1000 for 1 km)
make_site_buffer <- function(geojson_path, target_crs = NULL, radius_m = 1000) {
  aoi <- st_read(geojson_path, quiet = TRUE)
  
  # If the GeoJSON is a polygon AOI, use its centroid as the "site point"
  # If it is already a point, centroid() will return the point itself.
  site_pt <- st_centroid(st_union(aoi))
  
  if (!is.null(target_crs)) {
    site_pt <- st_transform(site_pt, target_crs)
  } else {
    # keep as-is for now; we'll transform later once we know chunk CRS
  }
  
  # Buffer expects units consistent with CRS. In projected CRS (meters), radius_m is meters.
  st_buffer(site_pt, dist = radius_m)
}

# Filter one chunk sf object to only features intersecting the buffer
filter_chunk_by_buffer <- function(chunk_sf, buffer_geom) {
  # Ensure same CRS
  if (st_crs(chunk_sf) != st_crs(buffer_geom)) {
    buffer_geom <- st_transform(buffer_geom, st_crs(chunk_sf))
  }
  
  idx <- st_intersects(chunk_sf, buffer_geom, sparse = FALSE)[, 1]
  chunk_sf[idx, , drop = FALSE]
}

# ---- main driver ----
extract_fetch_pixels <- function(
    site_dir,
    site_id,
    geojson_path,
    radius_m = 1000,
    in_subdir = "tables_sf",
    out_subdir = "tables_sf_fetch",
    overwrite = FALSE,
    write_combined = FALSE,
    return_objects = TRUE,
    export_global = FALSE,
    global_prefix = NULL
) {
  library(sf)
  
  in_dir  <- file.path(site_dir, in_subdir)
  out_dir <- file.path(site_dir, out_subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  chunk_files <- list.files(
    in_dir,
    pattern = "^chunk_\\d{3}_evi_sf\\.rds$",
    full.names = TRUE
  )
  if (length(chunk_files) == 0) {
    stop("No chunk sf RDS files found in: ", in_dir)
  }
  
  first_chunk <- readRDS(chunk_files[1])
  target_crs  <- st_crs(first_chunk)
  
  buffer_geom <- make_site_buffer(
    geojson_path = geojson_path,
    target_crs   = target_crs,
    radius_m     = radius_m
  )
  
  kept_list <- vector("list", length(chunk_files))
  names(kept_list) <- basename(chunk_files)
  
  for (i in seq_along(chunk_files)) {
    f_in  <- chunk_files[i]
    f_out <- file.path(out_dir, basename(f_in))
    
    if (file.exists(f_out) && !overwrite) {
      message("Skipping existing: ", basename(f_out))
      kept_list[[i]] <- NULL
      next
    }
    
    chunk_sf   <- readRDS(f_in)
    chunk_kept <- filter_chunk_by_buffer(chunk_sf, buffer_geom)
    
    if (nrow(chunk_kept) > 0) {
      saveRDS(chunk_kept, f_out)
    } else {
      message("No pixels intersect buffer in ", basename(f_in))
    }
    
    kept_list[[i]] <- chunk_kept
    if (i %% 25 == 0) message("Processed ", i, "/", length(chunk_files))
  }
  
  combined_path <- NULL
  combined_sf   <- NULL
  
  if (write_combined) {
    combined <- kept_list[!vapply(kept_list, is.null, logical(1))]
    combined <- combined[vapply(combined, nrow, integer(1)) > 0]
    
    if (length(combined) > 0) {
      combined_sf <- do.call(rbind, combined)
      combined_path <- file.path(site_dir, paste0(site_id, "_evi_timeseries_fetch_sf.rds"))
      saveRDS(combined_sf, combined_path)
      message("Wrote combined fetch sf: ", combined_path)
    } else {
      message("No intersecting pixels found across all chunks; no combined file written.")
    }
  }
  
  # ---- optionally export key objects to the global environment ----
  if (export_global) {
    pref <- if (is.null(global_prefix) || !nzchar(global_prefix)) "" else paste0(global_prefix, "_")
    
    assign(paste0(pref, "kept_list"), kept_list, envir = .GlobalEnv)
    assign(paste0(pref, "site_dir"),  site_dir,  envir = .GlobalEnv)
    assign(paste0(pref, "site_id"),   site_id,   envir = .GlobalEnv)
    
    # handy extras, if you want them
    assign(paste0(pref, "fetch_in_dir"),  in_dir,  envir = .GlobalEnv)
    assign(paste0(pref, "fetch_out_dir"), out_dir, envir = .GlobalEnv)
    assign(paste0(pref, "fetch_buffer"),  buffer_geom, envir = .GlobalEnv)
  }
  
  if (!return_objects) return(invisible(TRUE))
  
  list(
    site_dir = site_dir,
    site_id  = site_id,
    in_dir   = in_dir,
    out_dir  = out_dir,
    kept_list = kept_list,
    buffer_geom = buffer_geom,
    combined_path = combined_path,
    combined_sf = combined_sf
  )
}