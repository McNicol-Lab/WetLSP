rebuild_sf_from_parquet <- function(
    geom_dir,
    ts_dir,
    meta_path = NULL,   # e.g. file.path(out_meta, "meta.parquet")
    crs_wkt = NULL,     # optional override
    keep_pixel_id = TRUE,
    keep_xy = FALSE,
    keep_long_ts = FALSE,
    sort_by = c("cell", "pixel_id", "none"),
    quiet = FALSE
) {
  stopifnot(dir.exists(geom_dir), dir.exists(ts_dir))
  
  requireNamespace("arrow", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("sf", quietly = TRUE)
  
  sort_by <- match.arg(sort_by)
  
  msg <- function(...) if (!quiet) message(...)
  
  # ---- 1) Read metadata (CRS) ----
  if (is.null(crs_wkt) || !nzchar(crs_wkt)) {
    if (!is.null(meta_path) && file.exists(meta_path)) {
      meta <- arrow::read_parquet(meta_path)
      if (!all(c("key", "value") %in% names(meta))) {
        stop("meta_path must point to a parquet with columns {key, value}.")
      }
      crs_wkt <- meta$value[match("crs_wkt", meta$key)]
      if (length(crs_wkt) == 0 || is.na(crs_wkt) || !nzchar(crs_wkt)) {
        stop("Could not find non-empty crs_wkt in meta parquet at: ", meta_path)
      }
      msg("Loaded CRS WKT from meta.")
    } else {
      stop("CRS not provided. Supply either meta_path (recommended) or crs_wkt.")
    }
  }
  
  crs_obj <- sf::st_crs(crs_wkt)
  if (is.na(crs_obj)) stop("Invalid CRS WKT supplied/read from meta.")
  
  # ---- 2) Geometry: read and build sf points ----
  geom_ds  <- arrow::open_dataset(geom_dir, format = "parquet")
  geom_tbl <- dplyr::collect(geom_ds)
  
  req_geom <- c("pixel_id", "cell", "x", "y")
  missing_geom <- setdiff(req_geom, names(geom_tbl))
  if (length(missing_geom) > 0) {
    stop("Geometry dataset missing required columns: ", paste(missing_geom, collapse = ", "))
  }
  
  # coerce + de-dup (pixel_id is stable; cell should match pixel_id in your writer)
  geom_tbl <- geom_tbl %>%
    dplyr::mutate(
      pixel_id = as.integer(.data$pixel_id),
      cell     = as.integer(.data$cell),
      x        = as.numeric(.data$x),
      y        = as.numeric(.data$y)
    ) %>%
    dplyr::distinct(.data$pixel_id, .keep_all = TRUE)
  
  geom_sf <- sf::st_as_sf(
    geom_tbl,
    coords = c("x", "y"),
    crs = crs_obj
  )
  
  if (keep_xy) {
    # st_as_sf drops x/y columns; restore from geometry
    xy <- sf::st_coordinates(geom_sf$geometry)
    geom_sf$x <- as.numeric(xy[, 1])
    geom_sf$y <- as.numeric(xy[, 2])
  }
  
  # ---- 3) Timeseries: read and pivot to list-columns ----
  ts_ds  <- arrow::open_dataset(ts_dir, format = "parquet")
  ts_tbl <- dplyr::collect(ts_ds)
  
  req_ts <- c("pixel_id", "series", "date", "evi")
  missing_ts <- setdiff(req_ts, names(ts_tbl))
  if (length(missing_ts) > 0) {
    stop("Timeseries dataset missing required columns: ", paste(missing_ts, collapse = ", "))
  }
  
  ts_tbl <- ts_tbl %>%
    dplyr::mutate(
      pixel_id = as.integer(.data$pixel_id),
      series   = as.character(.data$series),
      date     = as.Date(.data$date),
      evi      = as.numeric(.data$evi)
    ) %>%
    dplyr::filter(!is.na(.data$date), !is.na(.data$evi))
  
  # Nest vectors per pixel_id x series
  ts_nested <- ts_tbl %>%
    dplyr::arrange(.data$pixel_id, .data$series, .data$date) %>%
    dplyr::group_by(.data$pixel_id, .data$series) %>%
    dplyr::summarise(
      dates = list(.data$date),
      evi   = list(.data$evi),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from  = .data$series,
      values_from = c(.data$dates, .data$evi),
      names_sep   = "_"
    )
  
  # We require spline series
  if (!all(c("dates_spline", "evi_spline") %in% names(ts_nested))) {
    stop("Timeseries parquet must include series == 'spline' (dates_spline/evi_spline missing).")
  }
  
  # Ensure raw list-columns exist even if not written
  if (!("dates_raw" %in% names(ts_nested))) {
    ts_nested$dates_raw <- replicate(nrow(ts_nested), list(as.Date(character(0))), simplify = FALSE)
  }
  if (!("evi_raw" %in% names(ts_nested))) {
    ts_nested$evi_raw <- replicate(nrow(ts_nested), list(numeric(0)), simplify = FALSE)
  }
  
  # ---- 4) Join geometry + timeseries ----
  out <- geom_sf %>%
    dplyr::left_join(ts_nested, by = "pixel_id")
  
  # Fill missing ts (should be rare)
  fill_dates <- function(x) if (is.null(x) || length(x) == 0) as.Date(character(0)) else x
  fill_num   <- function(x) if (is.null(x) || length(x) == 0) numeric(0) else x
  
  out$dates_spline <- lapply(out$dates_spline, fill_dates)
  out$evi_spline   <- lapply(out$evi_spline,   fill_num)
  out$dates_raw    <- lapply(out$dates_raw,    fill_dates)
  out$evi_raw      <- lapply(out$evi_raw,      fill_num)
  
  # ---- 5) Match original-style columns/order ----
  # Original: cell, dates_raw, evi_raw, dates_spline, evi_spline, geometry
  # We can keep pixel_id and/or x/y optionally.
  base_cols <- c("cell", "dates_raw", "evi_raw", "dates_spline", "evi_spline", "geometry")
  
  if (keep_pixel_id) {
    if (keep_xy) {
      out <- out %>% dplyr::select(pixel_id, dplyr::all_of(base_cols), x, y)
    } else {
      out <- out %>% dplyr::select(pixel_id, dplyr::all_of(base_cols))
    }
  } else {
    if (keep_xy) {
      out <- out %>% dplyr::select(dplyr::all_of(base_cols), x, y)
    } else {
      out <- out %>% dplyr::select(dplyr::all_of(base_cols))
    }
  }
  
  # ---- 6) Optional sorting ----
  if (sort_by == "cell") {
    out <- out %>% dplyr::arrange(.data$cell)
  } else if (sort_by == "pixel_id" && keep_pixel_id) {
    out <- out %>% dplyr::arrange(.data$pixel_id)
  }
  
  msg("Rebuilt sf: ", nrow(out), " pixels | cols: ", ncol(out))
  
  if (keep_long_ts) {
    return(list(sf = out, ts_long = ts_tbl))
  }
  
  out
}