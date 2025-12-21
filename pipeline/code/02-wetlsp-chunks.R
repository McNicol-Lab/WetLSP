#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# A High Spatial Resolution Land Surface Phenology Dataset for AmeriFlux & NEON
# 02: Process mosaiced PlanetScope images into base-aligned chunked .rda files
# Author: Minkyu Moon; Revised: Gavin McNicol; Hardened by ChatGPT :)
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
if (is.na(numSite)) stop("❌ Provide numSite as first argument")
message("📍 Running site index: ", numSite)

#---------------------------------------------
# Load parameters and helper functions
params <- fromJSON(file = "pipeline/wetlsp-parameters.json")
source(params$setup$rFunctions)

#---------------------------------------------
# Get site info (consistent with Script 01)
geojsonDir <- params$setup$geojsonDir
siteInfo   <- GetSiteInfo(numSite, geojsonDir, params)

imgDirSrc  <- siteInfo[[1]]
strSite    <- siteInfo[[2]]
cLong      <- siteInfo[[3]]
cLat       <- siteInfo[[4]]

message("🗺️  Site: ", strSite)
message("📂 Mosaic dir (input): ", file.path(params$setup$outDir, strSite, "mosaic"))

# Mosaic input directory (from Step 1)
mosaicDir <- file.path(params$setup$outDir, strSite, "mosaic")
if (!dir.exists(mosaicDir)) stop("❌ Mosaic directory not found: ", mosaicDir)

#---------------------------------------------
# Gather mosaiced image files
dfiles <- list.files(path = mosaicDir, pattern = "_clipped_mosaic\\.tif$", full.names = FALSE)
files  <- list.files(path = mosaicDir, pattern = "_clipped_mosaic\\.tif$", full.names = TRUE)
if (length(files) == 0) stop("❌ No mosaiced .tif files found in: ", mosaicDir)

# Extract valid dates from filenames (YYYYMMDD)
dates <- as.Date(regmatches(dfiles, regexpr("\\d{8}", dfiles)), "%Y%m%d")
ord   <- order(dates)
dates <- dates[ord]
files <- files[ord]
message("🧭 Found ", length(dates), " mosaics spanning ", min(dates), " → ", max(dates))

#---------------------------------------------
# Base raster for reference (defines grid for all chunks)
base_path <- file.path(params$setup$outDir, strSite, "base_image.tif")
if (!file.exists(base_path)) stop("❌ Base image missing: ", base_path)
imgBase <- terra::rast(base_path)

# Chunking plan (fixed by base grid)
numCk        <- params$setup$numChunks
total_cells  <- ncell(imgBase)
chunk_size   <- ceiling(total_cells / numCk)
message("📦 Base grid: ", ncol(imgBase), "×", nrow(imgBase), " (", total_cells,
        " cells) → ", numCk, " chunks (~", chunk_size, " cells/chunk)")

# Precompute fixed base chunk indices (spatially consistent across all dates)
base_chunks <- lapply(seq_len(numCk), function(cc){
  s <- (cc - 1) * chunk_size + 1
  e <- min(cc * chunk_size, total_cells)
  if (s > e) integer(0) else seq.int(s, e)
})

#---------------------------------------------
# Prepare output directories (skip if chunks already exist)
ckDir     <- file.path(params$setup$outDir, strSite, "chunk")
ckDirTemp <- file.path(ckDir, "temp")

numCk <- params$setup$numChunks

# Optional flag in JSON: "overwriteChunks": true/false
overwrite_chunks <- isTRUE(params$setup$overwriteChunks)

# If chunk dir exists and all expected chunk_XXX.rda files are present, skip work
if (dir.exists(ckDir) && !overwrite_chunks) {
  existing_chunks <- list.files(
    ckDir,
    pattern = "^chunk_\\d{3}\\.rda$",
    full.names = TRUE
  )
  
  if (length(existing_chunks) >= numCk) {
    message("✅ All ", numCk, " chunk files already exist in: ", ckDir)
    message("   Skipping chunk creation for site ", strSite,
            " (set 'overwriteChunks': true in wetlsp-parameters.json to force rebuild).")
    quit(save = "no", status = 0)
  } else {
    message("ℹ️  Chunk dir exists but only ", length(existing_chunks),
            " / ", numCk, " chunks found → will (re)build.")
  }
}

# Only get here if we need to (re)build chunks
dir.create(ckDir, recursive = TRUE, showWarnings = FALSE)
dir.create(ckDirTemp, recursive = TRUE, showWarnings = FALSE)

#---------------------------------------------
# Parallel backend setup (safe for terra)
num_cores <- min(params$setup$numCores, max(1, parallel::detectCores() - 1))
message("🔁 Using ", num_cores, " parallel workers")

if (.Platform$OS.type == "unix") {
  if (!requireNamespace("doMC", quietly = TRUE)) {
    install.packages("doMC", repos = "https://cran.rstudio.com/")
  }
  library(doMC)
  registerDoMC(cores = num_cores)
  message("✅ Using doMC shared-memory backend (no socket transfer)")
} else {
  cl <- parallel::makeCluster(num_cores, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  message("✅ Using PSOCK cluster backend (Windows fallback)")
}

# Keep terra from spinning up its own threads inside foreach
terra::terraOptions(threads = 1, memfrac = 0.75)

# ---------------------------------------------
# STEP 1 — Divide mosaics into base-aligned chunks
message("✂️  Step 1: Dividing mosaics into chunks and saving temporary .rda files")

is_readable_tif <- function(f) {
  r <- try(terra::rast(f), silent = TRUE)
  if (inherits(r, "try-error")) return(FALSE)
  
  # Force real pixel reads (like your plot/values call does)
  ok <- try(terra::global(r[[1]], "mean", na.rm = TRUE), silent = TRUE)
  !inherits(ok, "try-error")
}

good <- vapply(files, is_readable_tif, logical(1))
bad_files <- files[!good]

if (length(bad_files)) {
  message("⚠️ Skipping corrupted mosaics:")
  message(paste(basename(bad_files), collapse = ", "))
}

files <- files[good]
dates <- dates[good]

foreach(
  i = seq_along(dates),
  .packages = c("terra"),
  .export   = c("base_path", "ckDirTemp", "numCk", "base_chunks", "total_cells")
) %dopar% {
  date_str <- format(dates[i], "%Y%m%d")
  tif_path <- files[i]
  message("🔥 [PID ", Sys.getpid(), "] ", date_str, " → ", basename(tif_path))

  if (!file.exists(tif_path)) {
    message("⚠️  Missing file on disk: ", tif_path)
    return(NULL)
  }

  r <- try(terra::rast(tif_path), silent = TRUE)
  if (inherits(r, "try-error")) {
    message("⚠️  Failed to read raster for ", date_str)
    return(NULL)
  }
  
  # Force a read to catch corrupted strips/tiles early
  ok <- try(terra::global(r[[1]], "mean", na.rm = TRUE), silent = TRUE)
  if (inherits(ok, "try-error")) {
    message("⚠️  Corrupt/unreadable pixel blocks for ", date_str, " — skipping.  ", ok)
    return(NULL)
  }

  # Re-open base inside worker (don’t serialize SpatRaster)
  baseR <- terra::rast(base_path)

  # ✅ Force identical geometry to base grid (prevents striping later)
  if (!isTRUE(terra::compareGeom(r, baseR, stopOnError = FALSE))) {
    r <- tryCatch(
      {
        terra::resample(r, baseR, method = "bilinear")  # or "near"
      },
      error = function(e) {
        message("❌ [resample] failed for ", basename(tif_path), " — skipping.  ",
                conditionMessage(e))
        NULL
      }
    )
    if (is.null(r)) return(NULL)
  }

  vals <- terra::values(r, mat = TRUE)  # [ncell(base) x nbands]
  if (is.null(vals) || nrow(vals) != total_cells) {
    message("⚠️  Value matrix mismatch for ", date_str, " (", nrow(vals), " vs ", total_cells, ")")
    return(NULL)
  }

  # Split by precomputed base chunks
  for (cc in seq_len(numCk)) {
    idx <- base_chunks[[cc]]
    if (!length(idx)) next

    ckNum   <- sprintf("%03d", cc)
    dirTemp <- file.path(ckDirTemp, ckNum)
    if (!dir.exists(dirTemp)) dir.create(dirTemp, recursive = TRUE, showWarnings = FALSE)

    chunk_vals <- vals[idx, , drop = FALSE]  # rows by base cells
    # save as 8 vectors (b1..b8) for Step-3 compatibility
    b1 <- chunk_vals[,1]; b2 <- chunk_vals[,2]; b3 <- chunk_vals[,3]; b4 <- chunk_vals[,4]
    b5 <- chunk_vals[,5]; b6 <- chunk_vals[,6]; b7 <- chunk_vals[,7]; b8 <- chunk_vals[,8]

    save(file = file.path(dirTemp, paste0(date_str, ".rda")),
         list = c("b1","b2","b3","b4","b5","b6","b7","b8"))
  }

  rm(r, baseR, vals); gc()
  message("✅ ", date_str, " → ", numCk, " chunks written")
}

#---------------------------------------------
# STEP 2 — Merge chunk time series (robust to missing dates)
message("🧩 Step 2: Merging temporal chunks across dates (NA-fill missing)")

foreach(
  cc = seq_len(numCk),
  .packages = character(),
  .export   = c("ckDirTemp", "ckDir", "dates", "base_chunks")
) %dopar% {
  ckNum   <- sprintf("%03d", cc)
  dirTemp <- file.path(ckDirTemp, ckNum)
  
  if (!dir.exists(dirTemp)) {
    message("⚠️  Missing temp dir for chunk ", ckNum, " — skipping")
    return(NULL)
  }
  
  # Expected date strings (chronological, matches 'dates')
  expected_dates <- format(dates, "%Y%m%d")
  expected_paths <- file.path(dirTemp, paste0(expected_dates, ".rda"))
  
  exists_vec <- file.exists(expected_paths)
  n_found    <- sum(exists_vec)
  n_expected <- length(expected_paths)
  
  if (n_found == 0) {
    message("⚠️  Chunk ", ckNum, " has 0 / ", n_expected, " temp files — skipping")
    return(NULL)
  }
  
  if (n_found < n_expected) {
    miss <- expected_dates[!exists_vec]
    message("⚠️  Chunk ", ckNum, " missing ", (n_expected - n_found),
            " / ", n_expected, " dates → NA-fill. e.g. ",
            paste(head(miss, 5), collapse = ", "))
  }
  
  n_pix <- length(base_chunks[[cc]])   # rows for this chunk
  T     <- n_expected                  # number of time steps (dates)
  
  # Initialize band matrices [n_pix × T] with NA
  band_list <- lapply(1:8, function(b) matrix(NA_real_, nrow = n_pix, ncol = T))
  
  # For each expected date: load if present, otherwise keep NA
  for (t in seq_len(T)) {
    if (!exists_vec[t]) next
    
    e <- new.env(parent = emptyenv())
    ok <- try(load(expected_paths[t], envir = e), silent = TRUE)
    if (inherits(ok, "try-error")) {
      # If a temp .rda is itself corrupted/unreadable, treat as missing
      next
    }
    
    for (b in 1:8) {
      nm <- paste0("b", b)
      v  <- if (exists(nm, envir = e, inherits = FALSE)) get(nm, envir = e) else NULL
      
      if (!is.null(v)) {
        # guard length mismatch (should be exact, but keep safe)
        if (length(v) != n_pix) {
          len <- min(length(v), n_pix)
          tmp <- rep(NA_real_, n_pix)
          tmp[seq_len(len)] <- v[seq_len(len)]
          v <- tmp
        }
        band_list[[b]][, t] <- v
      }
    }
    
    rm(e)
  }
  
  # Save merged chunk
  save_path <- file.path(ckDir, paste0("chunk_", ckNum, ".rda"))
  
  save(
    list = c(paste0("band", 1:8), "dates"),
    file = save_path,
    envir = list2env(
      c(setNames(band_list, paste0("band", 1:8)),
        list(dates = dates)),
      parent = emptyenv()
    )
  )
  
  message("💾 Saved merged chunk: ", save_path, " (", n_found, "/", n_expected, " dates present)")
}

#---------------------------------------------
# STEP 3 — Cleanup temporary files
message("🧹 Cleaning up temporary directories…")
unlink(ckDirTemp, recursive = TRUE, force = TRUE)
message("✅ Chunk creation complete for site ", strSite)