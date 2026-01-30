# R/20_geometry_template.R
# Build buffer + template raster in sf_all CRS (no reprojection)

build_site_buffer <- function(geojson_path, target_crs, radius_m) {
  aoi <- sf::st_read(geojson_path, quiet = TRUE)
  site_pt <- sf::st_centroid(sf::st_union(aoi))
  site_pt <- sf::st_transform(site_pt, target_crs)
  sf::st_buffer(site_pt, dist = radius_m)
}

build_template_from_nc <- function(nc_path, target_crs_wkt, mask_vect) {
  # Read ONE layer as a file-backed SpatRaster (no values loaded yet)
  r0 <- terra::rast(nc_path, lyrs = 1)
  
  # NetCDF may not have CRS set; assign (do NOT project)
  if (is.na(terra::crs(r0)) || !nzchar(terra::crs(r0))) {
    terra::crs(r0) <- target_crs_wkt
  }
  
  # Crop using the vector directly (avoids ext() coercions; stays file-backed)
  r_crop <- terra::crop(r0, mask_vect, snap = "out")
  
  # Create a SMALL in-memory empty template on the cropped grid
  tmpl <- terra::rast(r_crop)
  terra::values(tmpl) <- NA_real_
  
  tmpl
}