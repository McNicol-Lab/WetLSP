#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# A High Spatial Resolution Land Surface Phenology Dataset for AmeriFlux and NEON Sites
#
# 03: A script for estimating phenometrics
# 
# Author: Minkyu Moon; moon.minkyu@gmail.com
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Submitted to the job for a single chunk for a site with shell script:
# #!/bin/bash
# echo Submitting $1
# R --vanilla < ~/03_LSP_script.R $1
#
# example submission command using default parameters:
# qsub -V -l h_rt=12:00:00 run_03.sh numSite chunk
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

library(raster)
# library(rgdal)
library(gdalUtilities)

library(rjson)
library(geojsonR)

library(doMC)
library(doParallel)

library(sf)
library(dplyr)
library(tidyr)


########################################
args <- commandArgs(trailingOnly = TRUE)
numSite <- as.numeric(args[1])
print(paste("Running site", numSite))

# numSite <- 1


########################################
## Load parameters
params <- fromJSON(file='PLSP_Parameters_gm_v1.json')
source(params$setup$rFunctions)


########################################
## Get site name, image directory and coordinate
strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

ckDir <- paste0(params$setup$outDir,strSite,'/chunk')
print(ckDir)

imgBase <- raster::raster(file.path(params$setup$outDir, strSite, "base_image.tif"))

## Load all chunks 
chunk_files <- list.files(path=ckDir,pattern=glob2rx(paste0('*.rda')),full.names=T)

# ## Load water mask
# waterRater <- raster(paste0(params$setup$outDir,strSite,'/water_mask_30_1.tif'))
# 
# numCk <- params$setup$numChunks
# chunk <- length(waterRater)%/%numCk
# if(cc==numCk){
#   chunks <- c((chunk*(cc-1)+1):length(waterRater))
# }else{
#   chunks <- c((chunk*(cc-1)+1):(chunk*cc))
# }
# waterMask <- values(waterRater)[chunks]

for(f in chunk_files) {
  
  cat("📦 Processing chunk file:", basename(f), "\n")
  load(f)
  
  ckNum <- gsub(".*chunk_|\\.rda", "", basename(f))  # extract '###'
  
  numChunks   <- params$setup$numChunks
  numPixSite  <- length(imgBase)
  chunkSize   <- numPixSite %/% numChunks
  ckIdx       <- as.integer(ckNum)
  
  if (ckIdx == numChunks) {
    cells <- (chunkSize * (ckIdx - 1) + 1):numPixSite
  } else {
    cells <- (chunkSize * (ckIdx - 1) + 1):(chunkSize * ckIdx)
  }
  coords <- xyFromCell(imgBase, cells)
  
  # per-chunk accumulators for list-columns
  n <- nrow(coords)
  ts_acc <- new.env(parent = emptyenv())
  ts_acc$raw_dates    <- vector("list", n)
  ts_acc$raw_evi      <- vector("list", n)
  ts_acc$smooth_dates <- vector("list", n)
  ts_acc$smooth_evi   <- vector("list", n)
  
  # callback sink: DoPhenologyPlanet will call this
  ts_sink <- function(pix_meta, dates_raw, evi_raw, pred_dates, evi_spline) {
    idx <- pix_meta$chunk_row  # 1..n within this chunk
    if (!is.null(dates_raw)) {
      ts_acc$raw_dates[[idx]] <<- as.Date(dates_raw)
      ts_acc$raw_evi[[idx]]   <<- as.numeric(evi_raw)
    }
    if (!is.null(pred_dates)) {
      if (is.null(ts_acc$smooth_dates[[idx]])) {
        ts_acc$smooth_dates[[idx]] <<- as.Date(pred_dates)
        ts_acc$smooth_evi[[idx]]   <<- as.numeric(evi_spline)
      } else {
        ts_acc$smooth_dates[[idx]] <<- c(ts_acc$smooth_dates[[idx]], as.Date(pred_dates))
        ts_acc$smooth_evi[[idx]]   <<- c(ts_acc$smooth_evi[[idx]],   as.numeric(evi_spline))
      }
    }
  }

  numPix <- dim(band1)[1]
  phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr
  pheno_mat <- matrix(NA, numPix, 24 * length(phenYrs))
  
  for (i in 1:numPix) {
    
    # if (i == 1 || i %% 1000 == 0) {
    #   print(paste("Checking pixel", i))
    #   print(summary(band2[i, ]))
    #   print(summary(band3[i, ]))
    #   print(summary(band5[i, ]))
    #   print(summary(band7[i, ]))
    # }
    
    
      pix_meta <- list(chunk_row = i, cell = cells[i], xy = coords[i, ])
      pheno_mat[i,] <- DoPhenologyPlanet(band2[i,], band4[i,],
                                         band6[i,], band8[i,],
                                         dates, phenYrs, params, waterMask = 0,
                                         ts_sink = ts_sink, pix_meta = pix_meta)
    
  
    #pheno_mat[i,] <- DoPhenologyPlanet(band2[i,], band4[i,],
                        #band6[i,], band8[i,],
                        #dates, phenYrs, params, waterMask = 0)
    
    # if (i) {
    #   if (all(res %in% c(NA, 4))) {
    #     cat("Pixel", i, "returned fallback phenology vector (mostly NA, with 4s).\n")
    #   } else {
    #     cat("Pixel", i, "returned real phenology metrics.\n")
    #     print(summary(res))
    #   }
    # }
    
   # cat("   ...pixel", i, "of", numPix, "\n")
  }
  
  cat("Raw time series stored:", sum(lengths(ts_acc$raw_evi) > 0), "\n")
  cat("Smoothed series stored:", sum(lengths(ts_acc$smooth_evi) > 0), "\n")
    
    sf_obj <- st_as_sf(
      data.frame(
        cell = cells,
        x    = coords[,1],
        y    = coords[,2]
      ),
      coords = c("x","y"),
      crs    = raster::crs(imgBase)
    ) %>%
      
    dplyr::mutate(
      dates_raw    = ts_acc$raw_dates,
      evi_raw      = ts_acc$raw_evi,
      dates_spline = ts_acc$smooth_dates,
      evi_spline   = ts_acc$smooth_evi
    )
  
  out_sf_dir <- file.path(params$setup$outDir, strSite, "tables_sf")
  dir.create(out_sf_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(sf_obj, file = file.path(out_sf_dir, paste0("chunk_", ckNum, "_evi_sf.rds")))
  cat("✅ Saved sf:", file.path(out_sf_dir, paste0("chunk_", ckNum, "_evi_sf.rds")), "\n")
  
  ckPheDir <- file.path(params$setup$outDir, strSite, 'chunk_phe')
  if (!dir.exists(ckPheDir)) dir.create(ckPheDir)
  save(pheno_mat, file = file.path(ckPheDir, paste0("chunk_phe_", ckNum, ".rda")))
  cat("✅ Saved:", paste0("chunk_phe_", ckNum, ".rda"), "\n")
  
  #save to csv
  # --- add coords back into the sf object ---
  coords_df <- st_coordinates(sf_obj)
  
  sf_obj2 <- sf_obj %>%
    st_drop_geometry() %>%
    mutate(x = coords_df[,1], 
           y = coords_df[,2])
  
  # --- unnest raw ---
  ts_raw <- sf_obj2 %>%
    select(cell, x, y, date = dates_raw, raw_value = evi_raw) %>%
    tidyr::unnest(c(date, raw_value))
  
  # --- unnest smoothed ---
  ts_smooth <- sf_obj2 %>%
    select(cell, x, y, date = dates_spline, smooth_value = evi_spline) %>%
    tidyr::unnest(c(date, smooth_value))
  
  # --- merge raw + smooth by cell,x,y,date ---
  ts_wide <- full_join(ts_raw, ts_smooth,
                       by = c("cell", "x", "y", "date")) %>%
    arrange(cell, date)
  
  # --- save CSV ---
  out_csv <- file.path(out_sf_dir, paste0("chunk_", ckNum, "_timeseries.csv"))
  write.csv(ts_wide, out_csv, row.names = FALSE)
  
  cat("✅ Saved per-pixel time series CSV:", out_csv, "\n")
}
