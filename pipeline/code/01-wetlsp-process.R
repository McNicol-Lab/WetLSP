#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# A High Spatial Resolution Land Surface Phenology Dataset for AmeriFlux and NEON Sites
#
# 01: A script for PlanetScope image process
# 
# Author: Minkyu Moon; moon.minkyu@gmail.com
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

required_packages <- c(
  "raster",
  # "rgdal",            # Deprecated; included with raster or sf now
  "gdalUtilities",
  # "rgeos",            # Deprecated; merged into sf/terra
  "terra",
  "rjson",
  "geojsonR",
  "dplyr",
  "foreach",
  "doParallel"
)

# Install any packages that are not already installed
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
    library(pkg, character.only = TRUE)
  } else {
    library(pkg, character.only = TRUE)
  }
}

invisible(lapply(required_packages, install_if_missing))

library(raster)
# library(rgdal)
library(gdalUtilities)
# library(rgeos)
library(terra)
library(rjson)
library(geojsonR)
library(dplyr)
library(foreach)
library(doParallel)

########################################
args <- commandArgs(trailingOnly = TRUE)
numSite <- as.numeric(args[1])
print(args)

########################################
## Load parameters
params <- fromJSON(file='pipeline/wetlsp-parameters.json')
source(params$setup$rFunctions)

########################################
## Get site name, image directory and coordinate
geojsonDir <- params$setup$geojsonDir

print(list.files(geojsonDir, pattern = "*.geojson"))
print(numSite)

siteInfo <- GetSiteInfo(numSite,geojsonDir,params)

imgDir <- siteInfo[[1]]; strSite <- siteInfo[[2]]
print(paste(strSite,';',imgDir))

cLong <- siteInfo[[3]];cLat <- siteInfo[[4]]
print(paste(cLong,';',cLat))



########################################
## Get list of files
dfileSR   <- list.files(path=imgDir, pattern=glob2rx('*MS_SR*.tif'), recursive=TRUE)
# dfileUDM  <- list.files(path=imgDir, pattern=glob2rx('*_DN_udm*.tif'), recursive=TRUE)
dfileUDM2 <- list.files(path=imgDir, pattern=glob2rx('*_udm2*.tif'), recursive=TRUE)

fileSR    <- list.files(path=imgDir, pattern=glob2rx('*MS_SR*.tif'), recursive=TRUE, full.names=TRUE)
# fileUDM   <- list.files(path=imgDir, pattern=glob2rx('*_DN_udm*.tif'), recursive=TRUE, full.names=TRUE)
fileUDM2  <- list.files(path=imgDir, pattern=glob2rx('*_udm2*.tif'), recursive=TRUE, full.names=TRUE)

# Build lookup table for SR data
sr_dates <- regmatches(dfileSR, regexpr("\\d{8}", dfileSR))
sr_dates_all <- as.Date(sr_dates, format = "%Y%m%d")
valid_sr <- !is.na(sr_dates_all)

file_table <- data.frame(
  file_name = dfileSR[valid_sr],
  full_path = fileSR[valid_sr],
  date = sr_dates_all[valid_sr],
  stringsAsFactors = FALSE
) %>% arrange(date)

# Unique valid dates
dates <- unique(file_table$date)

########################################
## Image process
outDir <- file.path(params$setup$outDir, strSite)
if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

siteWin <- GetSiteShp(fileSR, cLong, cLat)
imgBase <- GetBaseImg(fileSR, siteWin, outDir, save=TRUE)

num_cores <- min(params$setup$numCores, parallel::detectCores())
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat("🔁 Parallel backend registered with", num_cores, "workers\n")

outDir <- file.path(params$setup$outDir, strSite, "mosaic")
if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

# for(dd in 1:length(dates)){
  foreach(dd=1:length(dates)) %dopar% {
  cat(paste0("🔥 Started dd = ", dd, ", PID = ", Sys.getpid(), "\n"))
  
  files_today <- file_table[file_table$date == dates[dd], ]
  imgValid <- c()
  
  for (mm in seq_len(nrow(files_today))) {
    log <- try({
      img <- raster(files_today$full_path[mm]) 
      img <- crop(img, siteWin) 
    }, silent=TRUE)
    
    if (!inherits(log, 'try-error')) {
      if (nbands(raster(files_today$full_path[mm])) >= 4) {
        imgValid <- c(imgValid, mm)
      }
    }
  }

  if (length(imgValid) == 0) {
    cat("⛔ No valid images for date", dates[dd], "\n")
  }
  
  if (length(imgValid) > 0) {
    imgB <- vector('list', length(imgValid))
    
    for (nn in seq_along(imgValid)) {
      idx <- imgValid[nn]
      base_file <- files_today$file_name[idx]
      img <- brick(files_today$full_path[idx])
      numBand <- nbands(img)
      cat(paste0("🔢 Image has ", numBand, " bands\n"))
      
        imgT <- rast(files_today$full_path[idx])
        imgT <- crop(imgT, siteWin)
        
        cat("🔍 Checking for UDM2 match...\n")
        
        date_str <- regmatches(base_file, regexpr("\\d{8}", base_file))
        udm2_match <- grep(date_str, dfileUDM2)
        
        if (length(udm2_match) == length(imgValid)) {
          cat("☁️☁️ Found UDM2 — applying mask\n")
          log <- try({
            
            udm2T <- rast(fileUDM2[udm2_match[nn]])
            udm2T <- crop(udm2T, siteWin) 
            
            if (compareRaster(imgT, udm2T, extent=FALSE, rowcol=FALSE, crs=TRUE, stopiffalse=FALSE)) {
                  
              usable_mask <- udm2T[[1]] == 1 # Band 1 = clear 
              imgT <- raster::mask(imgT, usable_mask, maskvalue = 0)
              
            } else {
              cat("⚠️  Raster and UDM2 mismatch — skipped masking\n")
            }
          }, silent=TRUE)
        } else {
          cat("⚠️  No UDM2 match — skipping masking\n")
        }
        imgB[[nn]] <- imgT
      cat("✅ Finished processing image ", nn, "\n")
    }
    
    # Mosaic valid images
    rast_list <- vector("list", numBand)
    
    for (b in 1:numBand) {
      temp <- list( raster( lapply(imgB, function(x) x[[b]])[[1]] ) )
      temp[[length(temp)+1]] <- imgBase
      temp <- c(temp, fun = mean, na.rm = TRUE)
      rast_list[[b]] <- do.call(raster::mosaic, temp)
      rast_list[[b]] <- rast(rast_list[[b]])
    }
    
      Rast <- c(rast_list[[1]], 
                rast_list[[2]],
                rast_list[[3]], 
                rast_list[[4]],
                rast_list[[5]], 
                rast_list[[6]],
                rast_list[[7]], 
                rast_list[[8]])
      
    names(Rast) <- c("1","2","3","4","5","6","7","8")
    outFile <- file.path(outDir, paste0(format(dates[dd], "%Y%m%d"), "_cliped_mosaic.tif"))

    try({
      raster::writeRaster(Rast, filename = outFile, overwrite = TRUE)
      cat("✅ Saved to: ", outFile, "\n")
    }, silent = TRUE)
    
  }
}

# Summary
cat("✅ Mosaic process complete\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Starting date:", dates[dd], "\n")
cat("Number of output mosaics: ", length(list.files(path=outDir, pattern="\\.tif$")), "\n")
cat("Total dates processed: ", dd, "\n")