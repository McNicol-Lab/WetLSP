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
siteInfo   <- GetSiteInfo(numSite, geojsonDir, params)

imgDir  <- siteInfo[[1]]
strSite <- siteInfo[[2]]
cLong   <- siteInfo[[3]]
cLat    <- siteInfo[[4]]

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

#---------------------------------------------
# Match each SR scene to its own UDM2 by scene ID (filename stem before "_3B_"),
# not by date: dates with several scenes would otherwise share one mask.
udm2MissingAction <- params$setup$udm2MissingAction
if (is.null(udm2MissingAction)) udm2MissingAction <- "drop"
udm2MissingAction <- match.arg(udm2MissingAction, c("drop", "keep"))
udm2_idx <- MatchUDM2(file_table$file_name, dfileUDM2)
file_table$udm2_path <- ifelse(is.na(udm2_idx), NA_character_, fileUDM2[udm2_idx])
n_nomask <- sum(is.na(file_table$udm2_path))
n_multi  <- sum(table(file_table$date) > 1)
message("☁️  UDM2 matched by scene ID: ", nrow(file_table) - n_nomask, "/", nrow(file_table),
        " SR scenes; ", n_multi, " dates with >1 scene; udm2MissingAction = '", udm2MissingAction, "'")
if (n_nomask > 0) {
  message("⚠️  ", n_nomask, " SR scene(s) have no matching UDM2 (",
          ifelse(udm2MissingAction == "drop", "dropped", "kept UNMASKED"), "): ",
          paste(head(file_table$file_name[is.na(file_table$udm2_path)], 10), collapse = ", "),
          ifelse(n_nomask > 10, ", ...", ""))
}

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

# 🔧 ensure this really is a Date vector
dates_to_process <- sort(unique(as.Date(dates_to_process, origin = "1970-01-01")))

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
  .export   = c("file_table", "dfileUDM2", "fileUDM2", "siteWin", "outDirMosaic", "dates_to_process", "udm2MissingAction")
) %dopar% {
  
  # Wrap ENTIRE per-date body so one failure doesn't kill the job
  result <- tryCatch({
    
    current_date <- dates_to_process[dd]
    
    # extra safety: if the Date class has been dropped, restore it
    if (!inherits(current_date, "Date")) {
      current_date <- as.Date(current_date, origin = "1970-01-01")
    }
    
    message("🔥 [PID ", Sys.getpid(), "] Processing date ", current_date)
    
    files_today <- file_table[file_table$date == current_date, ]
    imgValid <- integer(0)
    
    # --- Quick raster validity check (SR images only) ---
    for (mm in seq_len(nrow(files_today))) {
      r <- try({
        img <- terra::rast(files_today$full_path[mm])
        img <- terra::crop(img, siteWin)  # may fail if extents don't overlap
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
    
    # --- Helper: simple extent-overlap test (no terra::relate) ---
    has_overlap <- function(a, b) {
      ea <- terra::ext(a); eb <- terra::ext(b)
      !(ea$xmax <= eb$xmin ||
          ea$xmin >= eb$xmax ||
          ea$ymax <= eb$ymin ||
          ea$ymin >= eb$ymax)
    }
    
    # --- Process valid rasters in-memory ---
    imgB <- vector("list", length(imgValid))
    
    for (nn in seq_along(imgValid)) {
      idx <- imgValid[nn]
      base_file <- files_today$file_name[idx]
      full_file <- files_today$full_path[idx]
      date_str  <- regmatches(base_file, regexpr("\\d{8}", base_file))
      
      message("🧩 Reading ", base_file)
      
      imgT <- try(terra::rast(full_file), silent = TRUE)
      if (inherits(imgT, "try-error")) {
        message("⚠️  Failed to read image for ", base_file, " – skipping.")
        next
      }
      
      # Crop to site window (already CRS-aligned at template stage)
      imgT <- try(terra::crop(imgT, siteWin), silent = TRUE)
      if (inherits(imgT, "try-error")) {
        message("⚠️  Crop to siteWin failed for ", base_file, " – skipping.")
        next
      }
      
      # --- Optional UDM2 masking (defensive) ---
      # UDM2 is pre-matched to this scene by scene ID (see MatchUDM2).
      masked <- FALSE
      udm2_path <- files_today$udm2_path[idx]
      if (!is.na(udm2_path)) {
        message("☁️  Applying UDM2 mask")
        
        udm2T <- try(terra::rast(udm2_path), silent = TRUE)
        if (!inherits(udm2T, "try-error")) {
          
          # Ensure CRS compatibility
          if (!isTRUE(terra::compareGeom(imgT, udm2T, crs = TRUE,
                                         ext = FALSE, rowcol = FALSE,
                                         stopOnError = FALSE))) {
            message("⚙️  Reprojecting UDM2 to match image CRS/grid for ", base_file)
            udm2T <- try(terra::project(udm2T, imgT), silent = TRUE)
            if (inherits(udm2T, "try-error")) {
              message("⚠️  UDM2 reprojection failed for ", base_file, " – skipping mask.")
              udm2T <- NULL
            }
          }
          
          if (!is.null(udm2T)) {
            # Check overlap BEFORE any crop/mask
            if (!has_overlap(udm2T, imgT)) {
              message("⚠️  UDM2 extent does not overlap image for ", base_file,
                      " – skipping mask.")
            } else {
              # Align extents roughly to imgT
              udm2T <- try(terra::crop(udm2T, imgT, snap = "out"), silent = TRUE)
              if (!inherits(udm2T, "try-error")) {
                mask_clear <- udm2T[[1]] == 1
                imgT_masked <- try(terra::mask(imgT, mask_clear, maskvalue = 0),
                                   silent = TRUE)
                if (!inherits(imgT_masked, "try-error")) {
                  imgT <- imgT_masked
                  masked <- TRUE
                } else {
                  message("⚠️  Masking failed for ", base_file, " – keeping unmasked image.")
                }
              } else {
                message("⚠️  UDM2 crop failed for ", base_file, " – skipping mask.")
              }
            }
          }
        } else {
          message("⚠️  Failed to read UDM2 for ", base_file, " – skipping mask.")
        }
      } else {
        message("⚠️  No UDM2 for scene ", base_file)
      }
      if (!masked && udm2MissingAction == "drop") {
        message("🚫 Dropping unmasked scene ", base_file, " (udm2MissingAction = 'drop')")
        next
      }
      
      imgB[[nn]] <- imgT
    }
    
    # --- Mosaic in-memory (safe version) ---
    valid_nonempty <- vapply(
      imgB,
      function(x) inherits(x, "SpatRaster") && nlyr(x) > 0,
      logical(1)
    )
    imgB <- imgB[valid_nonempty]
    
    if (length(imgB) == 0) {
      message("⛔ No usable rasters for ", current_date, " after masking.")
      return(NULL)
    }
    
    if (length(imgB) == 1) {
      Rast <- imgB[[1]]
      message("ℹ️  Only one valid raster for ", current_date, " — skipping mosaic.")
    } else {
      template <- imgB[[1]]
      imgB_aligned <- lapply(imgB, function(r) {
        if (!isTRUE(terra::compareGeom(r, template, stopOnError = FALSE))) {
          terra::project(r, template)
        } else {
          r
        }
      })
      
      Rast <- try(do.call(terra::mosaic, imgB_aligned), silent = TRUE)
      if (inherits(Rast, "try-error")) {
        message("⚠️  Mosaic failed for ", current_date, ", returning first raster instead.")
        Rast <- template
      }
    }
    
    # --- Write mosaic or single raster ---
    date_str <- strftime(current_date, "%Y%m%d")
    outFile <- file.path(outDirMosaic, paste0(date_str, "_clipped_mosaic.tif"))
    # NUM_THREADS=1: forked (doMC) workers deadlock inside GDAL's multithreaded
    # GeoTIFF compression once the parent process has written a compressed GeoTIFF
    # (GetBaseImg). Keep worker-side writes single-threaded.
    terra::writeRaster(Rast, outFile, overwrite = TRUE, filetype = "GTiff",
                       gdal = c("COMPRESS=LZW", "NUM_THREADS=1"))
    message("✅ Saved mosaic: ", outFile)
    
    return(outFile)
    
  }, error = function(e) {
    # This prevents foreach from bombing out and lets you see *which* date failed
    message("❌ ERROR on date ", dates_to_process[dd], " (index ", dd, "): ", conditionMessage(e))
    NULL
  })
  
  result
}

#---------------------------------------------
# Clean up cluster if used
if (exists("cl")) parallel::stopCluster(cl)

message("✅ Mosaic process complete.")
message(format(Sys.time()), " - Total mosaics written: ",
        length(list.files(outDirMosaic, pattern = "\\.tif$")))