#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# A High Spatial Resolution Land Surface Phenology Dataset for AmeriFlux & NEON
# 03: Estimate phenometrics from chunked mosaics (fast, base-aligned)
# Author: Minkyu Moon; Revised: Gavin McNicol; Parallel + smooth fix by ChatGPT
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# ----------------------------- Dependencies -----------------------------------
required_packages <- c(
  "terra","rjson","geojsonR","foreach","doParallel",
  "sf","dplyr","tidyr"
)
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

# ---------------------------------------------
# Base & I/O paths
base_path <- file.path(params$setup$outDir, strSite, "base_image.tif")
if (!file.exists(base_path)) stop("❌ Base image missing: ", base_path)
baseR <- terra::rast(base_path)

ckDir      <- file.path(params$setup$outDir, strSite, "chunk")
tablesDir  <- file.path(params$setup$outDir, strSite, "tables_sf")
pheDir     <- file.path(params$setup$outDir, strSite, "chunk_phe")
siteTblDir <- file.path(params$setup$outDir, strSite, "tables_site")

dir.create(tablesDir,  recursive = TRUE, showWarnings = FALSE)
dir.create(pheDir,     recursive = TRUE, showWarnings = FALSE)
dir.create(siteTblDir, recursive = TRUE, showWarnings = FALSE)

# Aggregated site-level CSV
aggFile <- file.path(siteTblDir, paste0(strSite, "_timeseries_all.csv"))
if (file.exists(aggFile)) file.remove(aggFile)

# Chunks from Step-02
chunk_files <- list.files(ckDir, pattern = "^chunk_\\d{3}\\.rda$", full.names = TRUE)
if (!length(chunk_files)) stop("❌ No chunk files found in: ", ckDir)
chunk_files <- sort(chunk_files)
message("🧭 Found ", length(chunk_files), " chunk files")

# ---------------------------- Execution Options --------------------------------
use_inner_parallel <- TRUE
n_cores_inner <- min(6, max(1, parallel::detectCores() - 1))  # tune as needed

if (use_inner_parallel) {
  if (.Platform$OS.type == "unix") {
    if (!requireNamespace("doMC", quietly = TRUE)) {
      install.packages("doMC", repos = "https://cran.rstudio.com/")
    }
    library(doMC)
    registerDoMC(cores = n_cores_inner)
    message("🔁 Using inner parallelism (", n_cores_inner, " cores per chunk)")
  } else {
    cl_inner <- parallel::makeCluster(n_cores_inner, type = "PSOCK")
    doParallel::registerDoParallel(cl_inner)
    on.exit(parallel::stopCluster(cl_inner), add = TRUE)
    message("🔁 Using inner parallelism via PSOCK (", n_cores_inner, " cores)")
  }
} else {
  message("🔁 Inner parallelism disabled (sequential within each chunk)")
}

terra::terraOptions(threads = 1, memfrac = 0.75)

# --------------------------- Helper: cells for chunk ---------------------------
chunk_cells_from_base <- function(base_rast, ckNum, numChunks){
  total_cells <- terra::ncell(base_rast)
  chunk_size  <- ceiling(total_cells / numChunks)
  s <- (ckNum - 1) * chunk_size + 1
  e <- min(ckNum * chunk_size, total_cells)
  if (s > e) integer(0) else seq.int(s, e)
}

# ------------------------------ Process Chunks ---------------------------------
phenYrs      <- params$setup$phenStartYr:params$setup$phenEndYr
phen_vec_len <- 24 * length(phenYrs)  # DoPhenologyPlanet length

for (f in chunk_files) {
  ckNum_str <- sub("^chunk_(\\d{3})\\.rda$", "\\1", basename(f))
  ckNum     <- as.integer(ckNum_str)
  if (is.na(ckNum)) { message("⚠️  Bad chunk name: ", f); next }

  phe_out_rda   <- file.path(pheDir,    paste0("chunk_phe_", ckNum_str, ".rda"))
  out_csv_chunk <- file.path(tablesDir, paste0("chunk_", ckNum_str, "_timeseries.csv"))

  if (file.exists(phe_out_rda) && file.exists(out_csv_chunk)) {
    message("⏭️  Skipping chunk ", ckNum_str, " (outputs already exist)")
    next
  }

  message("📦 Processing chunk ", ckNum_str, "  (", basename(f), ")")

  e <- new.env()
  load(f, envir = e)
  needed <- c(paste0("band", 1:8), "dates")
  if (!all(needed %in% ls(e))) {
    message("⚠️  Missing objects in chunk ", ckNum_str, ": ", paste(ls(e), collapse=", "))
    next
  }

  # Matrices [n_pix × n_time]
  B     <- lapply(1:8, function(b) get(paste0("band", b), envir = e))
  dates <- get("dates", envir = e)

  n_pix  <- nrow(B[[1]])
  n_time <- ncol(B[[1]])
  if (!all(vapply(B, function(m) nrow(m) == n_pix && ncol(m) == n_time, logical(1)))) {
    message("⚠️  Band dims differ in chunk ", ckNum_str, " — skipping")
    next
  }

  # Base-aligned pixel indices & coords for this chunk
  cells_chunk <- chunk_cells_from_base(baseR, ckNum, params$setup$numChunks)
  if (!length(cells_chunk)) { message("⚠️  Empty cell range for ", ckNum_str); next }
  if (length(cells_chunk) > n_pix) cells_chunk <- cells_chunk[seq_len(n_pix)]
  xy <- terra::xyFromCell(baseR, cells_chunk)

  # ---------------- Parallel pixel processing with local sink ----------------
  message("🧮 Processing ", n_pix, " pixels (", n_cores_inner, " cores)...")

  pheno_list <- if (use_inner_parallel) {
    foreach(
      i = seq_len(n_pix),
      .packages = character(),
      .export   = c("DoPhenologyPlanet","dates","phenYrs","params",
                    "cells_chunk","xy","phen_vec_len")
    ) %dopar% {
      pix_meta <- list(chunk_row = i, cell = cells_chunk[i], xy = xy[i, ])
      local_acc <- new.env(parent = emptyenv())
      local_acc$raw_dates <- NULL
      local_acc$raw_evi <- NULL
      local_acc$smooth_dates <- NULL
      local_acc$smooth_evi <- NULL

      local_sink <- function(pix_meta, dates_raw, evi_raw, pred_dates, evi_spline) {
        if (!is.null(dates_raw)) {
          local_acc$raw_dates <- as.Date(dates_raw)
          local_acc$raw_evi   <- as.numeric(evi_raw)
        }
        if (!is.null(pred_dates)) {
          local_acc$smooth_dates <- as.Date(pred_dates)
          local_acc$smooth_evi   <- as.numeric(evi_spline)
        }
      }

      res <- DoPhenologyPlanet(
        B[[2]][i,], B[[4]][i,], B[[6]][i,], B[[8]][i,],
        dates, phenYrs, params, waterMask = 0,
        ts_sink = local_sink, pix_meta = pix_meta
      )

      list(
        pheno = if (length(res) == phen_vec_len) res else rep(NA_real_, phen_vec_len),
        raw_dates = local_acc$raw_dates,
        raw_evi = local_acc$raw_evi,
        smooth_dates = local_acc$smooth_dates,
        smooth_evi = local_acc$smooth_evi
      )
    }
  } else {
    # Sequential fallback
    lapply(seq_len(n_pix), function(i) {
      pix_meta <- list(chunk_row = i, cell = cells_chunk[i], xy = xy[i, ])
      local_acc <- new.env(parent = emptyenv())
      local_sink <- function(pix_meta, dates_raw, evi_raw, pred_dates, evi_spline) {
        if (!is.null(dates_raw)) {
          local_acc$raw_dates <- as.Date(dates_raw)
          local_acc$raw_evi   <- as.numeric(evi_raw)
        }
        if (!is.null(pred_dates)) {
          local_acc$smooth_dates <- as.Date(pred_dates)
          local_acc$smooth_evi   <- as.numeric(evi_spline)
        }
      }

      res <- DoPhenologyPlanet(
        B[[2]][i,], B[[4]][i,], B[[6]][i,], B[[8]][i,],
        dates, phenYrs, params, waterMask = 0,
        ts_sink = local_sink, pix_meta = pix_meta
      )

      list(
        pheno = if (length(res) == phen_vec_len) res else rep(NA_real_, phen_vec_len),
        raw_dates = local_acc$raw_dates,
        raw_evi = local_acc$raw_evi,
        smooth_dates = local_acc$smooth_dates,
        smooth_evi = local_acc$smooth_evi
      )
    })
  }

  # Collect results back
  pheno_mat <- do.call(rbind, lapply(pheno_list, `[[`, "pheno"))
  ts_raw_dates    <- lapply(pheno_list, `[[`, "raw_dates")
  ts_raw_evi      <- lapply(pheno_list, `[[`, "raw_evi")
  ts_smooth_dates <- lapply(pheno_list, `[[`, "smooth_dates")
  ts_smooth_evi   <- lapply(pheno_list, `[[`, "smooth_evi")

  n_raw    <- sum(lengths(ts_raw_evi) > 0)
  n_smooth <- sum(lengths(ts_smooth_evi) > 0)
  message("🧾 Time series captured: raw=", n_raw, ", smooth=", n_smooth)

  # ----------------------------- Build SF + CSV data ---------------------------
  sf_obj <- st_as_sf(
    data.frame(cell = cells_chunk, x = xy[,1], y = xy[,2]),
    coords = c("x","y"),
    crs    = terra::crs(baseR)
  ) |>
    dplyr::mutate(
      dates_raw    = ts_raw_dates,
      evi_raw      = ts_raw_evi,
      dates_spline = ts_smooth_dates,
      evi_spline   = ts_smooth_evi
    )

  # Fallback for raw EVI if all missing
  if (n_raw == 0) {
    message("⚠️  No raw series from callback — computing directly from bands (2,6,8)")
    if (median(B[[6]], na.rm = TRUE) > 2) for (b in 1:8) B[[b]] <- B[[b]] / 10000
    evi_mat <- 2.5 * (B[[8]] - B[[6]]) / (B[[8]] + 6*B[[6]] - 7.5*B[[2]] + 1)
    evi_mat[!is.finite(evi_mat)] <- NA_real_
    evi_mat[evi_mat < -1] <- -1; evi_mat[evi_mat > 1] <- 1
    ts_raw_dates <- replicate(n_pix, dates, simplify = FALSE)
    ts_raw_evi   <- split(evi_mat, row(evi_mat))
    sf_obj$dates_raw <- ts_raw_dates
    sf_obj$evi_raw   <- ts_raw_evi
    n_raw <- n_pix
  }

  # ----------------------------- Save SF --------------------------------------
  saveRDS(sf_obj, file = file.path(tablesDir, paste0("chunk_", ckNum_str, "_evi_sf.rds")))
  message("💾 Saved sf: ", file.path(tablesDir, paste0("chunk_", ckNum_str, "_evi_sf.rds")))

  # ----------------------------- Save Phenology -------------------------------
  save(pheno_mat, file = phe_out_rda)
  message("💾 Saved phenology matrix: ", phe_out_rda)

 # ----------------------------- Save CSVs ------------------------------------
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
    
    # join and keep all available values
    ts_wide <- full_join(ts_raw, ts_smooth, by = c("cell","x","y","date")) |>
      arrange(cell, date)
    
    # --- FILTER LOGIC FIX ---
    # Keep all smooth dates within the configured phenology period,
    # and any raw values at their native dates
    if ("phenStartYr" %in% names(params$setup) && "phenEndYr" %in% names(params$setup)) {
      ts_wide <- ts_wide |>
        filter(lubridate::year(date) >= params$setup$phenStartYr &
               lubridate::year(date) <= params$setup$phenEndYr)
    }
    
    # drop empty rows
    ts_wide <- ts_wide %>%
      filter(!is.na(raw_value) | !is.na(smooth_value))
    
    if (nrow(ts_wide) == 0) {
      message("⚠️  No valid EVI rows to append for chunk ", ckNum_str)
      next
    }
    
    write.csv(ts_wide, out_csv_chunk, row.names = FALSE)
    message("💾 Saved per-chunk CSV: ", out_csv_chunk, " (", nrow(ts_wide), " rows)")
    
    # Append to aggregated site CSV
    if (file.exists(aggFile)) {
      write.table(ts_wide, aggFile, sep = ",", row.names = FALSE,
                  col.names = FALSE, append = TRUE)
    } else {
      write.table(ts_wide, aggFile, sep = ",", row.names = FALSE,
                  col.names = TRUE,  append = FALSE)
    }
    message("📈 Appended ", nrow(ts_wide), " rows → ", basename(aggFile))
}

message("✅ Step 3 complete for site ", strSite)

if (file.exists(aggFile)) {
  n_lines <- length(readLines(aggFile))
  message("🧮 Site-wide CSV: ", aggFile,
          " (", format(n_lines, big.mark=","), " lines)")
}