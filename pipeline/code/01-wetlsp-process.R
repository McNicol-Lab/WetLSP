#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# A High Spatial Resolution Land Surface Phenology Dataset for AmeriFlux and NEON Sites
#
# 01: A script for PlanetScope image process
# 
# Author: Minkyu Moon; moon.minkyu@gmail.com
# Revised: Gavin McNicol; gmcnicol@uic.edu
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Submitted to the job for each site with shell script:
# #!/bin/bash
# echo Submitting $1
# R --vanilla < ~/01_img_process.R $1
#
# example submission command using default parameters:
# qsub -V -pe omp 28 -l h_rt=12:00:00 run_01.sh numSite
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 01-image-process.R  —  Efficient, CyVerse-friendly refactor
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

required_packages <- c(
  "raster", "gdalUtilities", "terra", "rjson", "geojsonR",
  "dplyr", "foreach", "doParallel"
)

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
  }
  library(pkg, character.only = TRUE)
}
invisible(lapply(required_packages, install_if_missing))

library(terra)
library(raster)
library(foreach)
library(doParallel)

#---------------------------------------------
# Command-line argument
args <- commandArgs(trailingOnly = TRUE)
numSite <- as.numeric(args[1])
message("📍 Running site index: ", numSite)

#---------------------------------------------
# Load parameters and functions
params <- fromJSON(file = "pipeline/wetlsp-parameters.json")
source(params$setup$rFunctions)

#---------------------------------------------
# Get site info
geojsonDir <- params$setup$geojsonDir
siteInfo <- GetSiteInfo(numSite, geojsonDir, params)

imgDir  <- siteInfo[[1]]
strSite <- siteInfo[[2]]
cLong   <- siteInfo[[3]]
cLat    <- siteInfo[[4]]

message("🗺️  Site: ", strSite)
message("📂 Image dir: ", imgDir)

#---------------------------------------------
# Collect PlanetScope files
dfileSR   <- list.files(path = imgDir, pattern = glob2rx('*MS_SR*.tif'), recursive = TRUE)
fileSR    <- list.files(path = imgDir, pattern = glob2rx('*MS_SR*.tif'), recursive = TRUE, full.names = TRUE)
dfileUDM2 <- list.files(path = imgDir, pattern = glob2rx('*_udm2*.tif'), recursive = TRUE)
fileUDM2  <- list.files(path = imgDir, pattern = glob2rx('*_udm2*.tif'), recursive = TRUE, full.names = TRUE)

# Extract acquisition dates
sr_dates <- regmatches(dfileSR, regexpr("\\d{8}", dfileSR))
sr_dates_all <- as.Date(sr_dates, format = "%Y%m%d")
valid_sr <- !is.na(sr_dates_all)

file_table <- data.frame(
  file_name = dfileSR[valid_sr],
  full_path = fileSR[valid_sr],
  date = sr_dates_all[valid_sr],
  stringsAsFactors = FALSE
) %>% arrange(date)

dates <- unique(file_table$date)

#---------------------------------------------
# Prepare output and reference geometry
outDir <- file.path(params$setup$outDir, strSite)
if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

siteWin <- GetSiteShp(fileSR, cLong, cLat)

# Get base image and the exact path it was saved to
base_res <- GetBaseImg(fileSR, siteWin, outDir, save = TRUE)
imgBase  <- base_res$rast
base_path <- base_res$path  # may be in /home/jovyan/work if /data-store isn't writable

# Prefer opening inside each worker (lighter than exporting a big raster)
message("📄 Base image path for workers: ", ifelse(is.null(base_path), "<in-memory>", base_path))

# reopen base image if needed
if (!inherits(imgBase, "SpatRaster")) {
  imgBase <- terra::rast(base_path)
}

#---------------------------------------------
# Choose a fast, writable temp directory for terra
# Priority: /home/jovyan/work → /home/jovyan → /tmp
local_tmp_candidates <- c(
  "/home/jovyan/work",
  "/home/jovyan",
  "/tmp"
)

# pick the first that exists or can be created
local_tmp <- NULL
for (p in local_tmp_candidates) {
  if (dir.exists(p) || dir.create(p, recursive = TRUE, showWarnings = FALSE)) {
    local_tmp <- file.path(p, "terra_tmp")
    if (!dir.exists(local_tmp))
      dir.create(local_tmp, recursive = TRUE, showWarnings = FALSE)
    break
  }
}

if (is.null(local_tmp)) stop("❌ No writable temp directory found!")

# Clear existing terra temp directory manually if it exists
if (dir.exists(local_tmp)) {
  old_tmp <- list.files(local_tmp, full.names = TRUE)
  if (length(old_tmp) > 0) unlink(old_tmp, recursive = TRUE, force = TRUE)
}

terraOptions(
  tempdir = local_tmp,
  memfrac = 0.8,
  threads = 4
)

message("🌐 terra tempdir: ", terraOptions()$tempdir)
message("🌐 terra threads: ", terraOptions()$threads)


#---------------------------------------------
# Parallel backend setup (safe version)
num_cores <- min(8, parallel::detectCores() - 1)   # safer upper bound; adjust if you have more RAM
message("🔁 Setting up parallel backend with ", num_cores, " workers")

if (.Platform$OS.type == "unix") {
  # On CyVerse (Linux): shared-memory backend avoids serialization errors
  if (!requireNamespace("doMC", quietly = TRUE)) install.packages("doMC", repos = "https://cran.rstudio.com/")
  library(doMC)
  registerDoMC(cores = num_cores)
  message("✅ Using doMC (shared-memory) backend")
} else {
  # Fallback for Windows or non-Unix environments
  cl <- parallel::makeCluster(num_cores, timeout = 3600)
  doParallel::registerDoParallel(cl)
  message("✅ Using PSOCK cluster backend")
}


#---------------------------------------------
# Output directory for mosaics
outDirMosaic <- file.path(outDir, "mosaic")
if (!dir.exists(outDirMosaic)) dir.create(outDirMosaic, recursive = TRUE)

#---------------------------------------------
# Identify which mosaics already exist
existing_files <- list.files(outDirMosaic, pattern = "_clipped_mosaic\\.tif$", full.names = TRUE)
existing_dates <- as.Date(regmatches(existing_files, regexpr("\\d{8}", existing_files)), "%Y%m%d")

# Filter only missing dates
dates_to_process <- setdiff(dates, existing_dates)

if (length(dates_to_process) == 0) {
  message("✅ All mosaics already processed — nothing to do.")
  quit(save = "no")
}

message("🧭 Found ", length(existing_files), " existing mosaics. ",
        "Processing remaining ", length(dates_to_process), " dates.")

#---------------------------------------------
# Parallel loop across dates
foreach(
  dd = seq_along(dates_to_process),
  .packages = c("terra", "raster"),
  .export   = c("file_table", "dfileUDM2", "fileUDM2", "siteWin", "outDirMosaic")
) %dopar% {

  current_date <- dates_to_process[dd]
  message("🔥 [PID ", Sys.getpid(), "] Processing date ", current_date)

  files_today <- file_table[file_table$date == current_date, ]

  message("🔥 [PID ", Sys.getpid(), "] Processing date ", current_date)

  files_today <- file_table[file_table$date == current_date, ]
  imgValid <- integer(0)

  # --- Quick raster validity check
  for (mm in seq_len(nrow(files_today))) {
    r <- try({
      img <- terra::rast(files_today$full_path[mm])
      img <- terra::crop(img, siteWin)
      img
    }, silent = TRUE)
    if (!inherits(r, "try-error") && nlyr(r) >= 4) {
      imgValid <- c(imgValid, mm)
    }
  }

  if (length(imgValid) == 0) {
    message("⛔ No valid images for ", current_date)
    return(NULL)
  }

  # --- Process valid rasters in-memory
  imgB <- vector("list", length(imgValid))
  for (nn in seq_along(imgValid)) {
    idx <- imgValid[nn]
    base_file <- files_today$file_name[idx]
    date_str <- regmatches(base_file, regexpr("\\d{8}", base_file))
    message("🧩 Reading ", base_file)

    imgT <- terra::rast(files_today$full_path[idx])
    imgT <- terra::crop(imgT, siteWin)

    # --- Optional UDM2 masking ---
    udm2_match <- grep(date_str, dfileUDM2)
    if (length(udm2_match) > 0) {
      message("☁️  Applying UDM2 mask")
      udm2T <- try(terra::rast(fileUDM2[udm2_match[1]]), silent = TRUE)
      if (!inherits(udm2T, "try-error")) {
        udm2T <- terra::crop(udm2T, siteWin)

        # --- Align mask geometry to image before masking ---
        if (!isTRUE(terra::compareGeom(imgT, udm2T, stopOnError = FALSE))) {
        message("⚙️  Aligning UDM2 mask geometry to image for ", base_file)
        udm2T <- terra::project(udm2T, imgT)  # reprojects/resamples mask to match image grid
}

mask_clear <- udm2T[[1]] == 1
imgT <- try(terra::mask(imgT, mask_clear, maskvalue = 0), silent = TRUE)

if (inherits(imgT, "try-error")) {
  message("⚠️  Masking failed even after alignment, skipping mask for ", base_file)
  imgT <- terra::crop(terra::rast(files_today$full_path[idx]), siteWin)
}
      }
    } else {
      message("⚠️  No UDM2 match for ", date_str)
    }
    imgB[[nn]] <- imgT
  }

  # --- Mosaic in-memory (safe version) ---
valid_nonempty <- vapply(imgB, function(x) inherits(x, "SpatRaster") && nlyr(x) > 0, logical(1))
imgB <- imgB[valid_nonempty]

if (length(imgB) == 0) {
  message("⛔ No usable rasters for ", dates[dd], " after masking.")
  return(NULL)
}

if (length(imgB) == 1) {
  Rast <- imgB[[1]]
  message("ℹ️  Only one valid raster for ", dates[dd], " — skipping mosaic.")
} else {
  template <- imgB[[1]]
  imgB_aligned <- lapply(imgB, function(r) {
    if (!isTRUE(terra::compareGeom(r, template, stopOnError = FALSE))) {
      terra::project(r, template)
    } else r
  })
  Rast <- try(do.call(terra::mosaic, imgB_aligned), silent = TRUE)
  if (inherits(Rast, "try-error")) {
    message("⚠️  Mosaic failed for ", dates[dd], ", returning first raster instead.")
    Rast <- template
  }
}

# --- Write mosaic or single raster ---
outFile <- file.path(outDirMosaic, paste0(format(dates[dd], "%Y%m%d"), "_clipped_mosaic.tif"))
terra::writeRaster(Rast, outFile, overwrite = TRUE, filetype = "GTiff")
message("✅ Saved mosaic: ", outFile)

return(outFile)

}  # <-- this closes the foreach loop completely

#---------------------------------------------
# Clean up cluster if used
if (exists("cl")) parallel::stopCluster(cl)

message("✅ Mosaic process complete.")
message(format(Sys.time()), " - Total mosaics written: ",
        length(list.files(outDirMosaic, pattern = "\\.tif$")))