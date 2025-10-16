#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# A High Spatial Resolution Land Surface Phenology Dataset for AmeriFlux and NEON Sites
# 02: Process mosaiced PlanetScope images into chunked .rda files
# Author: Minkyu Moon; Revised: Gavin McNicol
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#---------------------------------------------
# Load dependencies
required_packages <- c("terra", "rjson", "geojsonR", "foreach", "doParallel")
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
  }
  library(pkg, character.only = TRUE)
}
invisible(lapply(required_packages, install_if_missing))

library(terra)
library(foreach)
library(doParallel)

#---------------------------------------------
# Command-line argument
args <- commandArgs(trailingOnly = TRUE)
numSite <- as.numeric(args[1])
message("📍 Running site index: ", numSite)

#---------------------------------------------
# Load parameters and helper functions
params <- fromJSON(file = "pipeline/wetlsp-parameters.json")
source(params$setup$rFunctions)

#---------------------------------------------
# Get site info (consistent with Script 01)
geojsonDir <- params$setup$geojsonDir
siteInfo <- GetSiteInfo(numSite, geojsonDir, params)

imgDir  <- siteInfo[[1]]
strSite <- siteInfo[[2]]
cLong   <- siteInfo[[3]]
cLat    <- siteInfo[[4]]

message("🗺️  Site: ", strSite)
message("📂 Image dir: ", imgDir)

# Mosaic output directory for this site
imgDir <- file.path(params$setup$outDir, strSite, "mosaic")
if (!dir.exists(imgDir)) stop("❌ Mosaic directory not found: ", imgDir)

#---------------------------------------------
# Gather mosaiced image files
dfiles <- list.files(path = imgDir, pattern = "_clipped_mosaic\\.tif$", full.names = FALSE)
files  <- list.files(path = imgDir, pattern = "_clipped_mosaic\\.tif$", full.names = TRUE)
if (length(files) == 0) stop("❌ No mosaiced .tif files found.")

# Extract valid dates from filenames (YYYYMMDD)
dates <- as.Date(regmatches(dfiles, regexpr("\\d{8}", dfiles)), "%Y%m%d")
message("🧭 Found ", length(dates), " mosaics spanning ", min(dates), " to ", max(dates))

#---------------------------------------------
# Base raster for reference
base_path <- file.path(params$setup$outDir, strSite, "base_image.tif")
if (!file.exists(base_path)) stop("❌ Base image missing: ", base_path)
imgBase <- terra::rast(base_path)

numCk <- params$setup$numChunks
total_cells <- ncell(imgBase)
chunk_size <- ceiling(total_cells / numCk)
message("📦 Base raster has ", total_cells, " cells → ", numCk, " chunks (~", chunk_size, " per chunk)")

#---------------------------------------------
# Prepare output directories
ckDir <- file.path(params$setup$outDir, strSite, "chunk")
ckDirTemp <- file.path(ckDir, "temp")

dir.create(ckDir, recursive = TRUE, showWarnings = FALSE)
dir.create(ckDirTemp, recursive = TRUE, showWarnings = FALSE)

#---------------------------------------------
# Parallel backend
num_cores <- min(params$setup$numCores, parallel::detectCores() - 1)
message("🔁 Using ", num_cores, " parallel workers")
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)
on.exit(parallel::stopCluster(cl), add = TRUE)

#---------------------------------------------
# STEP 1 — Divide mosaiced images into chunks
message("✂️  Step 1: Dividing mosaics into chunks and saving temporary .rda files")

foreach(i = seq_along(dates), .packages = "terra") %dopar% {
  date_str <- format(dates[i], "%Y%m%d")
  tif_path <- files[i]

  r <- try(terra::rast(tif_path), silent = TRUE)
  if (inherits(r, "try-error")) {
    message("⚠️  Skipping corrupt raster: ", tif_path)
    return(NULL)
  }

  vals <- terra::values(r)
  total_cells_i <- nrow(vals)
  chunk_size_i  <- ceiling(total_cells_i / numCk)

  message("📏 [", date_str, "] ", total_cells_i, " valid cells across ", numCk, " chunks")

  for (cc in seq_len(numCk)) {
    start_idx <- (cc - 1) * chunk_size_i + 1
    end_idx   <- min(cc * chunk_size_i, total_cells_i)
    if (start_idx > total_cells_i) next

    ckNum <- sprintf("%03d", cc)
    dirTemp <- file.path(ckDirTemp, ckNum)
    if (!dir.exists(dirTemp)) dir.create(dirTemp, recursive = TRUE, showWarnings = FALSE)

    chunk_vals <- vals[start_idx:end_idx, , drop = FALSE]
    b1 <- chunk_vals[, 1]; b2 <- chunk_vals[, 2]; b3 <- chunk_vals[, 3]; b4 <- chunk_vals[, 4]
    b5 <- chunk_vals[, 5]; b6 <- chunk_vals[, 6]; b7 <- chunk_vals[, 7]; b8 <- chunk_vals[, 8]

    save_path <- file.path(dirTemp, paste0(date_str, ".rda"))
    save(b1, b2, b3, b4, b5, b6, b7, b8, file = save_path)
  }

  message("✅ [", date_str, "] Split into ", numCk, " chunks")
}

#---------------------------------------------
# STEP 2 — Merge chunk time series
message("🧩 Step 2: Merging temporal chunks across dates")

foreach(cc = seq_len(numCk)) %dopar% {
  ckNum <- sprintf("%03d", cc)
  dirTemp <- file.path(ckDirTemp, ckNum)
  chunk_files <- list.files(dirTemp, full.names = TRUE)

  if (length(chunk_files) != length(dates)) {
    message("⚠️  Chunk ", ckNum, " incomplete (", length(chunk_files), "/", length(dates), ") — skipping")
    return(NULL)
  }

  load_rda <- function(f) {
    e <- new.env()
    load(f, envir = e)
    as.list(e)
  }

  loaded <- lapply(chunk_files, load_rda)

  # --- Find common valid length across all dates ---
  lengths_per_date <- vapply(loaded, function(x) length(x$b1), numeric(1))
  n_common <- min(lengths_per_date)
  T <- length(loaded)
  message("🧮 Chunk ", ckNum, ": using ", n_common, " pixels × ", T, " dates")

  # --- Initialize band matrices ---
  band_list <- lapply(1:8, function(b) matrix(NA_real_, nrow = n_common, ncol = T))

  # --- Fill matrices safely (truncate/pad as needed) ---
  for (t in seq_along(loaded)) {
    for (b in 1:8) {
      v <- loaded[[t]][[paste0("b", b)]]
      if (length(v) > n_common) v <- v[seq_len(n_common)]
      if (length(v) < n_common) v <- c(v, rep(NA_real_, n_common - length(v)))
      band_list[[b]][, t] <- v
    }
  }

  # --- Save merged chunk ---
  save_path <- file.path(ckDir, paste0("chunk_", ckNum, ".rda"))
  save(list = c(paste0("band", 1:8), "dates"),
       file = save_path,
       envir = list2env(c(setNames(band_list, paste0("band", 1:8)),
                          list(dates = dates))))
  message("💾 Saved merged chunk: ", save_path)
}

#---------------------------------------------
# STEP 3 — Cleanup temporary files
message("🧹 Cleaning up temporary directories...")
unlink(ckDirTemp, recursive = TRUE, force = TRUE)
message("✅ Chunk creation complete for site ", strSite)