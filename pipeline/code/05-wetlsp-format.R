#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# A High Spatial Resolution Land Surface Phenology Dataset for AmeriFlux & NEON Sites
# 05: Export yearly phenology GeoTIFFs → CF-compliant NetCDF
# Author: Minkyu Moon; modernized by Gavin McNicol
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# ----------------------------- Dependencies -----------------------------------
required_packages <- c("terra","rjson","geojsonR","ncdf4")
install_if_missing <- function(pkg){
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
  }
  library(pkg, character.only = TRUE)
}
invisible(lapply(required_packages, install_if_missing))
library(terra)
library(ncdf4)

# ---------------------------- Arguments & Params -------------------------------
args <- commandArgs(trailingOnly = TRUE)
numSite <- as.numeric(args[1])
if (is.na(numSite)) stop("❌ Provide numSite as first argument")
message("📍 Running site index: ", numSite)

params_path <- "pipeline/wetlsp-parameters.json"
if (!file.exists(params_path)) params_path <- "PLSP_Parameters_gm_v1.json"
params <- fromJSON(file = params_path)
source(params$setup$rFunctions)

# ------------------------------- Site Lookup -----------------------------------
geojsonDir <- params$setup$geojsonDir
siteInfo   <- GetSiteInfo(numSite, geojsonDir, params)
strSite <- siteInfo[[2]]
message("🗺️  Site: ", strSite)

# ------------------------------ Directories ------------------------------------
inBase  <- file.path(params$setup$outDir, "Product_GeoTiff", strSite)
outBase <- file.path(params$setup$outDir, "Product_netCDF", strSite)
dir.create(outBase, recursive = TRUE, showWarnings = FALSE)

# ------------------------------ Inputs -----------------------------------------
productTable <- read.csv(params$setup$productTable, stringsAsFactors = FALSE)
# make sure numeric columns are numeric
num_cols <- c("fill_value","scale","offset","valid_min","valid_max")
for (cc in num_cols) if (cc %in% names(productTable)) productTable[[cc]] <- as.numeric(productTable[[cc]])

phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr

base_path <- file.path(params$setup$outDir, strSite, "base_image.tif")
if (!file.exists(base_path)) stop("❌ Base image missing: ", base_path)
baseR <- terra::rast(base_path)

# ------------------------------ Define dimensions safely -----------------------
nx <- ncol(baseR); ny <- nrow(baseR)
resxy <- terra::res(baseR)  # c(resx, resy)
xminv <- terra::xmin(baseR); xmaxv <- terra::xmax(baseR)
yminv <- terra::ymin(baseR); ymaxv <- terra::ymax(baseR)

# build by length.out to avoid rounding drifts
x_vals <- as.numeric(seq(from = xminv + resxy[1]/2, by = resxy[1], length.out = nx))
y_vals <- as.numeric(seq(from = yminv + resxy[2]/2, by = resxy[2], length.out = ny))

if (length(x_vals) != nx || length(y_vals) != ny) {
  stop("❌ x/y coordinate vector lengths do not match raster dimensions")
}

dimx <- ncdim_def(
  name = "x", longname = "x coordinate", units = "m",
  vals = x_vals, create_dimvar = TRUE
)
dimy <- ncdim_def(
  name = "y", longname = "y coordinate", units = "m",
  vals = rev(y_vals), create_dimvar = TRUE  # CF expects y decreasing north→south
)

# ------------------------------ Yearly Loop ------------------------------------
for (yy in phenYrs) {
  yearDir <- file.path(inBase, as.character(yy))
  if (!dir.exists(yearDir)) { message("⚠️  Missing GeoTIFFs for ", yy); next }

  tif_files <- list.files(yearDir, pattern="\\.tif$", full.names=TRUE)
  if (!length(tif_files)) { message("⚠️  No TIFFs for ", yy); next }

  outFile <- file.path(outBase, paste0("WetLSP_", yy, ".nc"))
  if (file.exists(outFile)) { message("⏭️  NetCDF already exists for ", yy); next }

  message("💾 Creating NetCDF for ", strSite, " — ", yy)

  # Define variables ------------------------------------------------------------
  results <- vector("list", nrow(productTable) + 1)
  results[[1]] <- ncvar_def("transverse_mercator","",list(),prec="char")

  for (i in seq_len(nrow(productTable))) {
    lyr  <- productTable[i,]
    prec <- ifelse(isTRUE(lyr$data_type == "Int16"), "short", "float")
    results[[i+1]] <- ncvar_def(
      name     = lyr$short_name,
      units    = lyr$units,
      dim      = list(dimx, dimy),      # order: x, y
      missval  = as.numeric(lyr$fill_value),
      longname = lyr$long_name,
      prec     = prec,
      compression = 4
    )
  }

  # Check dim classes for ncdf4 (ncdim4)
  if (!inherits(dimx, "ncdim4") || !inherits(dimy, "ncdim4")) {
    stop("❌ Invalid NetCDF dimensions — check x/y definitions.")
  }

  if (file.exists(outFile)) file.remove(outFile)
  ncout <- nc_create(outFile, results, force_v4 = TRUE)

  # Write each product layer ----------------------------------------------------
  for (i in seq_len(nrow(productTable))) {
    lyr <- productTable[i,]
    # Step 4 names: "%02d_%Y_%short.tif" → match on _shortname.tif
    tif <- grep(paste0("_", lyr$short_name, "\\.tif$"), tif_files, value=TRUE)
    if (!length(tif)) { message("⚠️  Missing layer ", lyr$short_name, " for ", yy); next }

    r <- terra::rast(tif[1])

    # ✅ Always resample to match base grid (ensures consistent geometry)
    if (!isTRUE(terra::compareGeom(r, baseR, stopOnError = FALSE))) {
      message("↩️  Resampling ", lyr$short_name, " to base grid")
      r <- terra::resample(r, baseR, method = "bilinear")
    }
    
    # Flip vertically (so north is up) and transpose to [x, y] order
    m_yx <- as.matrix(terra::flip(r, "vertical"), wide = TRUE)
    m_xy <- t(m_yx)

    # apply valid range and fill
    vmin <- as.numeric(lyr$valid_min); vmax <- as.numeric(lyr$valid_max)
    fill <- as.numeric(lyr$fill_value)
    m_xy[m_xy < vmin | m_xy > vmax | is.na(m_xy)] <- fill

    ncvar_put(ncout, results[[i+1]], m_xy)

    # per-variable attributes
    if ("scale" %in% names(lyr))   ncatt_put(ncout, lyr$short_name, "scale",     as.numeric(lyr$scale))
    if ("offset" %in% names(lyr))  ncatt_put(ncout, lyr$short_name, "offset",    as.numeric(lyr$offset))
    if ("data_type" %in% names(lyr)) ncatt_put(ncout, lyr$short_name, "data_type", as.character(lyr$data_type))
    ncatt_put(ncout, lyr$short_name, "valid_min", vmin)
    ncatt_put(ncout, lyr$short_name, "valid_max", vmax)
  }

  # ------------------------- Projection attributes -----------------------------
  crs_wkt <- terra::crs(baseR)
  ncatt_put(ncout, "transverse_mercator", "long_name",           "CRS definition")
  ncatt_put(ncout, "transverse_mercator", "grid_mapping_name",   "transverse_mercator")
  ncatt_put(ncout, "transverse_mercator", "spatial_ref",         gsub("\\", "", crs_wkt, fixed=TRUE))
  ncatt_put(ncout, "transverse_mercator", "GeoTransform",
            paste(xminv, resxy[1], 0, ymaxv, 0, -resxy[2]))

  # ------------------------- Global attributes ---------------------------------
  ncatt_put(ncout, 0, "title",               "Wetland Surface Phenology from PlanetScope (WetLSP)")
  ncatt_put(ncout, 0, "summary",             "High-resolution Land-Surface Phenology for FLUXNET-CH4 wetlands")
  ncatt_put(ncout, 0, "product_version",     "v001")
  ncatt_put(ncout, 0, "software_repository", "https://github.com/McNicol-Lab/WetLSP")

  ncatt_put(ncout, 0, "origin_creator_name", "Land Cover & Surface Climate Group, Boston University")
  ncatt_put(ncout, 0, "fork_creator_name",   "Ecosystem and Planetary Health Integration Lab, UIC")
  ncatt_put(ncout, 0, "origin_contributor_name",
            "Minkyu Moon, A.R. Richardson, T. Milliman, M.A. Friedl")
  ncatt_put(ncout, 0, "fork_contributor_name",
            "Gavin McNicol, Aarin Bothra, Youmi Oh, Lichen Liu, Sparkle Malone, Qing Zhu")

  ncatt_put(ncout, "x", "axis", "projection_x_coordinate")
  ncatt_put(ncout, "y", "axis", "projection_y_coordinate")

  nc_close(ncout)
  message("✅ Saved NetCDF: ", outFile)
}

message("🎉 Step 5 complete for site ", strSite)