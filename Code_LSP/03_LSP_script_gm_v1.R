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
    
    
      pheno_mat[i,] <- DoPhenologyPlanet(band2[i,], band4[i,],
                        band6[i,], band8[i,],
                        dates, phenYrs, params, waterMask = 0)
    
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
  
  ckPheDir <- file.path(params$setup$outDir, strSite, 'chunk_phe')
  if (!dir.exists(ckPheDir)) dir.create(ckPheDir)
  save(pheno_mat, file = file.path(ckPheDir, paste0("chunk_phe_", ckNum, ".rda")))
  cat("✅ Saved:", paste0("chunk_phe_", ckNum, ".rda"), "\n")
  
}
