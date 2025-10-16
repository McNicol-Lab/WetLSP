#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# A High Spatial Resolution Land Surface Phenology Dataset for AmeriFlux and NEON Sites
# 03: Estimate phenometrics from chunked mosaics
# Author: Minkyu Moon; Revised: Gavin McNicol (modernized)
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# ----------------------------- Dependencies -----------------------------------
required_packages <- c("terra","rjson","geojsonR","foreach","doParallel","sf","dplyr","tidyr")
install_if_missing <- function(pkg){
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
  }
  library(pkg, character.only = TRUE)
}
invisible(lapply(required_packages, install_if_missing))

library(terra)
library(foreach)
library(doParallel)
library(sf)
library(dplyr)
library(tidyr)

# ---------------------------- Arguments & Params -------------------------------
args <- commandArgs(trailingOnly = TRUE)
numSite <- as.numeric(args[1])
if (is.na(numSite)) stop("❌ Provide numSite as first argument")
message("📍 Running site index: ", numSite)

params <- fromJSON(file = "pipeline/wetlsp-parameters.json")
source(params$setup$rFunctions)

# ------------------------------- Site Lookup -----------------------------------
geojsonDir <- params$setup$geojsonDir
siteInfo   <- GetSiteInfo(numSite, geojsonDir, params)

imgDir  <- siteInfo[[1]]
strSite <- siteInfo[[2]]
cLong   <- siteInfo[[3]]
cLat    <- siteInfo[[4]]

message("🗺️  Site: ", strSite)
message("📂 Image dir: ", imgDir)

# ------------------------------ Paths & Inputs ---------------------------------
base_path <- file.path(params$setup$outDir, strSite, "base_image.tif")
if (!file.exists(base_path)) stop("❌ Base image missing: ", base_path)
baseR <- terra::rast(base_path)

ckDir     <- file.path(params$setup$outDir, strSite, "chunk")
tablesDir <- file.path(params$setup$outDir, strSite, "tables_sf")
pheDir    <- file.path(params$setup$outDir, strSite, "chunk_phe")

dir.create(tablesDir, recursive = TRUE, showWarnings = FALSE)
dir.create(pheDir,    recursive = TRUE, showWarnings = FALSE)

chunk_files <- list.files(ckDir, pattern = "^chunk_\\d{3}\\.rda$", full.names = TRUE)
if (!length(chunk_files)) stop("❌ No chunk files found in: ", ckDir)
message("🧭 Found ", length(chunk_files), " chunk files")

# ---------------------------- Execution Options --------------------------------
# Toggle per-pixel parallelism inside a chunk (OFF by default—CPU heavy, high overhead)
use_inner_parallel <- FALSE
num_cores <- min(params$setup$numCores, parallel::detectCores() - 1)
if (use_inner_parallel) {
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  message("🔁 Inner parallelism enabled with ", num_cores, " workers")
} else {
  message("🔁 Inner parallelism disabled (sequential within each chunk)")
}

# --------------------------- Helper: coords by chunk ---------------------------
# Compute base-raster cell indices for a given chunk number (1..numChunks)
chunk_cells_from_base <- function(base_rast, ckNum, numChunks){
  total_cells <- ncell(base_rast)
  chunk_size  <- ceiling(total_cells / numChunks)
  start_idx   <- (ckNum - 1) * chunk_size + 1
  end_idx     <- min(ckNum * chunk_size, total_cells)
  if (start_idx > end_idx) return(integer(0))
  seq(start_idx, end_idx)
}

# ------------------------------ Process Chunks ---------------------------------
phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr
phen_vec_len <- 24 * length(phenYrs)  # expected length returned by DoPhenologyPlanet

for (f in chunk_files) {
  ckNum_str <- sub("^chunk_(\\d{3})\\.rda$", "\\1", basename(f))
  ckNum     <- as.integer(ckNum_str)
  if (is.na(ckNum)) {
    message("⚠️  Skipping malformed chunk filename: ", f)
    next
  }

  # Skip if outputs exist (resume-safe)
  phe_out_rda  <- file.path(pheDir,    paste0("chunk_phe_", ckNum_str, ".rda"))
  timeseries_csv <- file.path(tablesDir, paste0("chunk_", ckNum_str, "_timeseries.csv"))
  if (file.exists(phe_out_rda) && file.exists(timeseries_csv)) {
    message("⏭️  Skipping chunk ", ckNum_str, " (outputs already exist)")
    next
  }

  message("📦 Processing chunk ", ckNum_str, "  (", basename(f), ")")
  e <- new.env()
  load(f, envir = e)

  # Expect these objects from Step 2: band1..band8, dates
  needed <- c(paste0("band", 1:8), "dates")
  if (!all(needed %in% ls(e))) {
    message("⚠️  Chunk ", ckNum_str, " missing expected objects; found: ", paste(ls(e), collapse = ", "))
    next
  }

  # Matrices [rows = pixels (n_common), cols = dates]
  B <- lapply(1:8, function(b) get(paste0("band", b), envir = e))
  dates <- get("dates", envir = e)

  # Make sure dims line up
  n_pix  <- nrow(B[[1]])
  n_time <- ncol(B[[1]])
  if (!all(vapply(B, function(m) nrow(m) == n_pix && ncol(m) == n_time, logical(1)))) {
    message("⚠️  Band dimensions differ in chunk ", ckNum_str, " — skipping")
    next
  }

  # Get base cell ids and coordinates for this chunk, then align to n_pix (truncate if needed)
  cells_chunk <- chunk_cells_from_base(baseR, ckNum, params$setup$numChunks)
  if (!length(cells_chunk)) {
    message("⚠️  Empty cell range for chunk ", ckNum_str, " — skipping")
    next
  }
  # If merged chunk (Step 2) is shorter than base chunk, truncate cells/coords to n_pix
  if (length(cells_chunk) > n_pix) cells_chunk <- cells_chunk[seq_len(n_pix)]
  xy <- terra::xyFromCell(baseR, cells_chunk)

  # Accumulators for saving time series per pixel
  ts_acc <- new.env(parent = emptyenv())
  ts_acc$raw_dates    <- vector("list", n_pix)
  ts_acc$raw_evi      <- vector("list", n_pix)
  ts_acc$smooth_dates <- vector("list", n_pix)
  ts_acc$smooth_evi   <- vector("list", n_pix)

  # Callback sink for DoPhenologyPlanet
  ts_sink <- function(pix_meta, dates_raw, evi_raw, pred_dates, evi_spline) {
    idx <- pix_meta$chunk_row  # 1..n_pix
    if (!is.null(dates_raw)) {
      ts_acc$raw_dates[[idx]] <<- as.Date(dates_raw)
      ts_acc$raw_evi[[idx]]   <<- as.numeric(evi_raw)
    }
    if (!is.null(pred_dates)) {
      # Store or append smoothed series
      if (is.null(ts_acc$smooth_dates[[idx]])) {
        ts_acc$smooth_dates[[idx]] <<- as.Date(pred_dates)
        ts_acc$smooth_evi[[idx]]   <<- as.numeric(evi_spline)
      } else {
        ts_acc$smooth_dates[[idx]] <<- c(ts_acc$smooth_dates[[idx]], as.Date(pred_dates))
        ts_acc$smooth_evi[[idx]]   <<- c(ts_acc$smooth_evi[[idx]],   as.numeric(evi_spline))
      }
    }
  }

  # Pre-allocate phenology output
  pheno_mat <- matrix(NA_real_, nrow = n_pix, ncol = phen_vec_len)

  # Choose sequential or inner-parallel per-pixel
  if (use_inner_parallel) {
    pheno_mat <- foreach(i = seq_len(n_pix), .combine = rbind, .packages = character()) %dopar% {
      pix_meta <- list(chunk_row = i, cell = cells_chunk[i], xy = xy[i, ])
      res <- DoPhenologyPlanet(
        B[[2]][i,], B[[4]][i,], B[[6]][i,], B[[8]][i,],
        dates, phenYrs, params, waterMask = 0,
        ts_sink = ts_sink, pix_meta = pix_meta
      )
      if (length(res) != phen_vec_len) {
        # pad or truncate safely
        len <- min(length(res), phen_vec_len)
        out <- rep(NA_real_, phen_vec_len)
        out[seq_len(len)] <- res[seq_len(len)]
        out
      } else res
    }
  } else {
    for (i in seq_len(n_pix)) {
      if (i %% 2000 == 0 || i == 1) message("   … pixel ", i, "/", n_pix)
      pix_meta <- list(chunk_row = i, cell = cells_chunk[i], xy = xy[i, ])
      res <- DoPhenologyPlanet(
        B[[2]][i,], B[[4]][i,], B[[6]][i,], B[[8]][i,],
        dates, phenYrs, params, waterMask = 0,
        ts_sink = ts_sink, pix_meta = pix_meta
      )
      if (length(res) != phen_vec_len) {
        len <- min(length(res), phen_vec_len)
        pheno_mat[i,] <- c(res[seq_len(len)], rep(NA_real_, phen_vec_len - len))
      } else {
        pheno_mat[i,] <- res
      }
    }
  }

  message("🧾 Time series captured: raw=", sum(lengths(ts_acc$raw_evi) > 0),
          ", smooth=", sum(lengths(ts_acc$smooth_evi) > 0))

  # ----------------------------- Save SF --------------------------------------
  sf_obj <- st_as_sf(
    data.frame(cell = cells_chunk, x = xy[,1], y = xy[,2]),
    coords = c("x","y"),
    crs    = terra::crs(baseR)
  ) |>
    dplyr::mutate(
      dates_raw    = ts_acc$raw_dates,
      evi_raw      = ts_acc$raw_evi,
      dates_spline = ts_acc$smooth_dates,
      evi_spline   = ts_acc$smooth_evi
    )

  saveRDS(sf_obj, file = file.path(tablesDir, paste0("chunk_", ckNum_str, "_evi_sf.rds")))
  message("💾 Saved sf: ", file.path(tablesDir, paste0("chunk_", ckNum_str, "_evi_sf.rds")))

  # ----------------------------- Save Phenology -------------------------------
  save(pheno_mat, file = phe_out_rda)
  message("💾 Saved phenology matrix: ", phe_out_rda)

  # ----------------------------- Save CSV (wide) ------------------------------
  coords_df <- sf::st_coordinates(sf_obj)
  sf_obj2 <- sf_obj |>
    st_drop_geometry() |>
    mutate(x = coords_df[,1], y = coords_df[,2])

  ts_raw <- sf_obj2 |>
    select(cell, x, y, date = dates_raw, raw_value = evi_raw) |>
    tidyr::unnest(c(date, raw_value))

  ts_smooth <- sf_obj2 |>
    select(cell, x, y, date = dates_spline, smooth_value = evi_spline) |>
    tidyr::unnest(c(date, smooth_value))

  ts_wide <- dplyr::full_join(ts_raw, ts_smooth, by = c("cell","x","y","date")) |>
    arrange(cell, date)

  out_csv <- file.path(tablesDir, paste0("chunk_", ckNum_str, "_timeseries.csv"))
  write.csv(ts_wide, out_csv, row.names = FALSE)
  message("💾 Saved per-pixel time series CSV: ", out_csv)
}

message("✅ Step 3 complete for site ", strSite)