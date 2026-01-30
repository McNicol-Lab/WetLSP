# R/10_io_resolve.R
# IO helpers: resolve NetCDF + load sf_all + normalize dates

resolve_wetlsp_nc <- function(site_dir, site_id, year, nc_subdir = "netcdf") {
  nc_dir <- file.path(site_dir, nc_subdir)
  nc_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
  
  if (length(nc_files) == 0) {
    stop("No .nc files found in: ", nc_dir)
  }
  
  hits <- nc_files[
    grepl(site_id, basename(nc_files), fixed = TRUE) &
      grepl(as.character(year), basename(nc_files), fixed = TRUE)
  ]
  
  if (length(hits) == 0) {
    hits <- nc_files[grepl(as.character(year), basename(nc_files), fixed = TRUE)]
  }
  
  if (length(hits) == 0) {
    stop(
      "No NetCDF matched year ", year, " in: ", nc_dir,
      "\nFiles seen:\n", paste(basename(nc_files), collapse = "\n")
    )
  }
  
  hits <- hits[order(file.info(hits)$mtime, decreasing = TRUE)]
  hits[[1]]
}

read_sf_all_fetch <- function(site_dir, site_id, radius_m,
                              suffix = "_evi_timeseries_fetch_sf_100m.rds") {
  # keep flexible: you can change naming later
  # If you prefer to key by radius, pass a different suffix
  path <- file.path(site_dir, paste0(site_id, suffix))
  if (!file.exists(path)) stop("sf_all file not found: ", path)
  readRDS(path)
}

coerce_dates_list_to_Date <- function(x_list) {
  lapply(x_list, function(x) {
    if (inherits(x, "Date")) return(x)
    if (inherits(x, "POSIXt")) return(as.Date(x))
    if (is.numeric(x)) return(as.Date(x, origin = "1970-01-01"))
    as.Date(x) # character fallback
  })
}