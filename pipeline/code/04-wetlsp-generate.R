#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# A High Spatial Resolution Land Surface Phenology Dataset for AmeriFlux and NEON Sites
# 04: Save phenology layers to GeoTIFF (per year, per product)
# Author: Minkyu Moon; Modernized: Gavin McNicol
# Updated: Always overwrite existing GeoTIFFs (no skipping)
#          FIX: Each output GeoTIFF is a single-layer phenology product
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# ----------------------------- Dependencies -----------------------------------
required_packages <- c("terra","rjson","geojsonR")
install_if_missing <- function(pkg){
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
  }
  library(pkg, character.only = TRUE)
}
invisible(lapply(required_packages, install_if_missing))
library(terra)

# ---------------------------- Arguments & Params -------------------------------
args <- commandArgs(trailingOnly = TRUE)
numSite <- as.numeric(args[1])
if (is.na(numSite)) stop("❌ Provide numSite as first argument")
message("📍 Running site index: ", numSite)

# Prefer the unified parameters file; fall back to legacy if needed
params_path <- "pipeline/wetlsp-parameters.json"
if (!file.exists(params_path)) params_path <- "PLSP_Parameters_gm_v1.json"
params <- fromJSON(file = params_path)
source(params$setup$rFunctions)

# Product table
if (!file.exists(params$setup$productTable)) {
  stop("❌ productTable not found: ", params$setup$productTable)
}
productTable <- read.csv(params$setup$productTable, header = TRUE, stringsAsFactors = FALSE)
if (nrow(productTable) < 24 || !"short_name" %in% names(productTable)) {
  stop("❌ productTable must have at least 24 rows and a 'short_name' column.")
}
phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr

# ------------------------------- Site Lookup -----------------------------------
geojsonDir <- params$setup$geojsonDir
siteInfo   <- GetSiteInfo(numSite, geojsonDir, params)

imgDir  <- siteInfo[[1]]
strSite <- siteInfo[[2]]
message("🗺️  Site: ", strSite)
message("📂 Image dir: ", imgDir)

# ----------------------- Base geometry & chunk inputs --------------------------
base_path <- file.path(params$setup$outDir, strSite, "base_image.tif")
if (!file.exists(base_path)) stop("❌ Base image missing: ", base_path)
baseR <- terra::rast(base_path)

# ✅ Use a single-layer template (one band) for all phenology outputs
template_r <- baseR[[1]]

ckPheDir <- file.path(params$setup$outDir, strSite, "chunk_phe")
if (!dir.exists(ckPheDir)) stop("❌ chunk_phe directory not found: ", ckPheDir)

# All chunk phenology files produced in Step 3
phe_files <- list.files(ckPheDir, pattern = "^chunk_phe_\\d{3}\\.rda$", full.names = TRUE)
if (!length(phe_files)) stop("❌ No chunk phenology files found in: ", ckPheDir)
message("🧭 Found ", length(phe_files), " chunk phenology files")

# Output root
pheDirRoot <- file.path(params$setup$outDir, "Product_GeoTiff", strSite)
dir.create(pheDirRoot, recursive = TRUE, showWarnings = FALSE)

# ----------------------------- terra options -----------------------------------
local_tmp_candidates <- c("/home/jovyan/work", "/home/jovyan", "/tmp")
local_tmp <- NULL
for (p in local_tmp_candidates) {
  if (dir.exists(p) || dir.create(p, recursive = TRUE, showWarnings = FALSE)) {
    local_tmp <- file.path(p, "terra_tmp")
    if (!dir.exists(local_tmp)) dir.create(local_tmp, recursive = TRUE, showWarnings = FALSE)
    break
  }
}
if (is.null(local_tmp)) stop("❌ No writable temp directory found!")

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

# --------------------------- Helpers & indexing --------------------------------
numPix     <- ncell(baseR)
numChunks  <- params$setup$numChunks
chunk_size <- ceiling(numPix / numChunks)

chunk_cells_from_base <- function(base_rast, ckNum, numChunks){
  total_cells <- ncell(base_rast)
  ch_size     <- ceiling(total_cells / numChunks)
  start_idx   <- (ckNum - 1) * ch_size + 1
  end_idx     <- min(ckNum * ch_size, total_cells)
  if (start_idx > end_idx) integer(0) else seq(start_idx, end_idx)
}

# ------------------------------ Main per-year loop -----------------------------
for (y in phenYrs) {
  ylab <- as.character(y)
  message("📅 Year ", ylab, " — assembling 24 products")
  
  pheDirYear <- file.path(pheDirRoot, ylab)
  dir.create(pheDirYear, recursive = TRUE, showWarnings = FALSE)
  
  # Always rebuild — overwrite existing GeoTIFFs
  vals_year <- matrix(NA_real_, nrow = numPix, ncol = 24)
  
  # Fill from each chunk file
  for (cf in phe_files) {
    ckNum_str <- sub("^chunk_phe_(\\d{3})\\.rda$", "\\1", basename(cf))
    if (is.na(suppressWarnings(as.integer(ckNum_str)))) {
      message("⚠️  Skipping malformed chunk file: ", cf)
      next
    }
    ckNum <- as.integer(ckNum_str)
    
    e <- new.env()
    log <- try(load(cf, envir = e), silent = TRUE)
    if (inherits(log, "try-error") || !"pheno_mat" %in% ls(e)) {
      message("⚠️  Could not load pheno_mat from: ", cf)
      next
    }
    pm <- get("pheno_mat", envir = e)
    
    if (!is.matrix(pm) || ncol(pm) %% 24 != 0) {
      message("⚠️  pheno_mat shape unexpected in ", cf, " — skipping")
      next
    }
    
    nYears_pm <- ncol(pm) / 24
    years_seq <- params$setup$phenStartYr:params$setup$phenEndYr
    if (length(years_seq) != nYears_pm) {
      years_seq <- seq_len(nYears_pm) + params$setup$phenStartYr - 1L
    }
    
    y_idx <- which(years_seq == y)
    if (length(y_idx) != 1) {
      message("⚠️  Year ", ylab, " not present in chunk ", ckNum_str, " — skipping file")
      next
    }
    col_start <- (y_idx - 1L) * 24L + 1L
    col_end   <- y_idx * 24L
    pm_year   <- pm[, col_start:col_end, drop = FALSE]
    
    cells_chunk <- chunk_cells_from_base(baseR, ckNum, numChunks)
    if (!length(cells_chunk)) {
      message("⚠️  Empty cell range for chunk ", ckNum_str)
      next
    }
    
    n_target <- length(cells_chunk)
    n_src    <- nrow(pm_year)
    n_copy   <- min(n_target, n_src)
    if (n_copy <= 0) next
    
    vals_year[cells_chunk[seq_len(n_copy)], ] <- pm_year[seq_len(n_copy), ]
  }
  
  # --------------------- Write 24 GeoTIFFs for this year -----------------------
  message("💾 Writing (overwriting) 24 single-layer GeoTIFFs for ", ylab)
  for (k in 1:24) {
    out_name <- file.path(
      pheDirYear,
      sprintf("%02d_%s_%s.tif", k, ylab, productTable$short_name[k])
    )
    
    # Single-layer raster with correct geometry
    r_k <- template_r
    values(r_k) <- vals_year[, k]
    
    terra::writeRaster(
      r_k, out_name, overwrite = TRUE, filetype = "GTiff"
    )
  }
  
  message("✅ Finished year ", ylab)
}

# Optional cleanup of intermediates
if (isTRUE(params$setup$cleanupIntermediates)) {
  site_root <- file.path(params$setup$outDir, strSite)
  dirs_to_remove <- c("mosaic", "chunk", "chunk_phe", "tables_sf")
  for (d in dirs_to_remove) {
    p <- file.path(site_root, d)
    if (dir.exists(p)) {
      message("🧹 Removing intermediate directory: ", p)
      unlink(p, recursive = TRUE, force = TRUE)
    }
  }
  message("🧼 Intermediate cleanup complete for site ", strSite)
}

message("🎉 Step 4 complete for site ", strSite)