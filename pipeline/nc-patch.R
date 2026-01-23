# Patch CF-ish coord attrs + ensure grid_mapping on all data variables
# Requires: ncdf4

patch_cf_coords_and_grid_mapping <- function(nc_path, gridmap_var = "transverse_mercator") {
  stopifnot(file.exists(nc_path))
  if (!requireNamespace("ncdf4", quietly = TRUE)) install.packages("ncdf4")
  
  nc <- ncdf4::nc_open(nc_path, write = TRUE)
  on.exit(try(ncdf4::nc_close(nc), silent = TRUE), add = TRUE)
  
  # ---- 1) Fix x/y attributes ----
  for (d in c("x", "y")) {
    if (d %in% names(nc$dim)) {
      # Remove non-CF 'axis: projection_x_coordinate' style if present
      old_axis <- try(ncdf4::ncatt_get(nc, d, "axis")$value, silent = TRUE)
      if (!inherits(old_axis, "try-error") && is.character(old_axis) &&
          grepl("^projection_[xy]_coordinate$", old_axis)) {
        try(ncdf4::ncatt_del(nc, d, "axis"), silent = TRUE)
      }
      
      if (d == "x") {
        ncdf4::ncatt_put(nc, d, "standard_name", "projection_x_coordinate")
        ncdf4::ncatt_put(nc, d, "axis", "X")
        ncdf4::ncatt_put(nc, d, "long_name", "x coordinate of projection")
      } else {
        ncdf4::ncatt_put(nc, d, "standard_name", "projection_y_coordinate")
        ncdf4::ncatt_put(nc, d, "axis", "Y")
        ncdf4::ncatt_put(nc, d, "long_name", "y coordinate of projection")
      }
    }
  }
  
  # ---- 2) Add grid_mapping to all data variables (except dims + projection var) ----
  var_names <- names(nc$var)
  for (v in var_names) {
    if (identical(v, gridmap_var)) next
    
    # Skip 1D coordinate-like vars if present (rare in your case)
    if (v %in% c("x", "y")) next
    
    # Add/overwrite grid_mapping
    ncdf4::ncatt_put(nc, v, "grid_mapping", gridmap_var)
    
    # Optional: help some readers (safe even if redundant)
    # ncdf4::ncatt_put(nc, v, "coordinates", "x y")
  }
  
  invisible(TRUE)
}

# Example:
# patch_cf_coords_and_grid_mapping("output/AT-Nsd/netcdf/WetLSP_2021_annotated.nc")