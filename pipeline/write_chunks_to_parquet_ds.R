write_chunks_to_parquet_ds <- function(
    chunk_files,
    out_geom,
    out_ts,
    out_meta,
    site_id,
    radius_m,
    batch_k = 10L,
    overwrite = TRUE,
    prefer_pixel_id = c("cell", "row_number")  # "cell" strongly recommended
) {
  prefer_pixel_id <- match.arg(prefer_pixel_id)
  
  dir.create(out_geom, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_ts,   recursive = TRUE, showWarnings = FALSE)
  dir.create(out_meta, recursive = TRUE, showWarnings = FALSE)
  
  if (overwrite) {
    unlink(out_geom, recursive = TRUE, force = TRUE)
    unlink(out_ts,   recursive = TRUE, force = TRUE)
    unlink(out_meta, recursive = TRUE, force = TRUE)
    dir.create(out_geom, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_ts,   recursive = TRUE, showWarnings = FALSE)
    dir.create(out_meta, recursive = TRUE, showWarnings = FALSE)
  }
  
  requireNamespace("sf", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("arrow", quietly = TRUE)
  
  # ---- capture CRS once ----
  crs_wkt <- NULL
  
  # ---- batching ----
  geom_batch <- list()
  ts_batch   <- list()
  b <- 0L
  batch_id <- 0L
  total_pix <- 0L
  
  # helper: longify
  longify <- function(df, dates_col, evi_col, series_name) {
    df |>
      sf::st_drop_geometry() |>
      dplyr::select(pixel_id, dates = dplyr::all_of(dates_col), evi = dplyr::all_of(evi_col)) |>
      tidyr::unnest(c(dates, evi)) |>
      dplyr::filter(!is.na(.data$evi)) |>
      dplyr::transmute(
        pixel_id = as.integer(.data$pixel_id),
        series   = series_name,
        date     = as.Date(.data$dates),
        year     = as.integer(format(.data$date, "%Y")),
        evi      = as.numeric(.data$evi)
      )
  }
  
  for (i in seq_along(chunk_files)) {
    f <- chunk_files[i]
    sf_i <- readRDS(f)
    if (nrow(sf_i) == 0) next
    
    # CRS
    if (is.null(crs_wkt)) crs_wkt <- sf::st_crs(sf_i)$wkt
    
    # pixel_id
    if (prefer_pixel_id == "cell") {
      if (!("cell" %in% names(sf_i))) {
        stop("Expected column 'cell' but it was not found in: ", basename(f))
      }
      sf_i$pixel_id <- as.integer(sf_i$cell)
    } else {
      sf_i$pixel_id <- seq_len(nrow(sf_i))
    }
    
    # geometry table (store x/y; keep cell for exact round-trip)
    xy <- sf::st_coordinates(sf_i$geometry)
    geom_batch[[length(geom_batch) + 1]] <- tibble::tibble(
      pixel_id = as.integer(sf_i$pixel_id),
      cell     = as.integer(sf_i$cell),
      x        = as.numeric(xy[, 1]),
      y        = as.numeric(xy[, 2])
    )
    
    # timeseries: spline + (optional) raw
    if (!all(c("dates_spline","evi_spline") %in% names(sf_i))) {
      stop("Missing spline cols in: ", basename(f))
    }
    ts_spline <- longify(sf_i, "dates_spline", "evi_spline", "spline")
    
    if (all(c("dates_raw","evi_raw") %in% names(sf_i))) {
      ts_raw <- longify(sf_i, "dates_raw", "evi_raw", "raw")
      ts_i <- dplyr::bind_rows(ts_spline, ts_raw)
    } else {
      ts_i <- ts_spline
    }
    
    ts_batch[[length(ts_batch) + 1]] <- ts_i
    
    b <- b + 1L
    total_pix <- total_pix + nrow(sf_i)
    
    if (b %% batch_k == 0L || i == length(chunk_files)) {
      batch_id <- batch_id + 1L
      message("Flushing batch ", batch_id, " at file ", i, "/", length(chunk_files),
              " | cumulative pixels: ", total_pix)
      
      geom_out <- dplyr::bind_rows(geom_batch) |>
        dplyr::distinct(.data$pixel_id, .keep_all = TRUE)
      
      ts_out <- dplyr::bind_rows(ts_batch)
      
      arrow::write_parquet(
        geom_out,
        sink = file.path(out_geom, sprintf("geom_batch_%03d.parquet", batch_id)),
        compression = "zstd"
      )
      
      arrow::write_parquet(
        ts_out,
        sink = file.path(out_ts, sprintf("ts_batch_%03d.parquet", batch_id)),
        compression = "zstd"
      )
      
      geom_batch <- list()
      ts_batch   <- list()
    }
  }
  
  if (is.null(crs_wkt)) stop("No non-empty inputs; CRS not captured.")
  
  meta <- tibble::tibble(
    key   = c("crs_wkt", "site_id", "radius_m", "pixel_id_source"),
    value = c(crs_wkt, site_id, as.character(radius_m), prefer_pixel_id)
  )
  
  arrow::write_parquet(meta, file.path(out_meta, "meta.parquet"), compression = "zstd")
  
  invisible(list(out_geom = out_geom, out_ts = out_ts, out_meta = out_meta))
}