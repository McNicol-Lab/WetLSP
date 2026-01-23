combine_fetch_sf_progress <- function(
    fetch_dir,
    site_dir,
    site_id,
    pattern = "^chunk_\\d{3}_evi_sf\\.rds$",
    progress_every = 10
) {
  
  files <- list.files(fetch_dir, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No fetch chunk files found in: ", fetch_dir)
  }
  
  cat("Found ", length(files), " fetch chunks\n")
  flush.console()
  
  cat("Combining with progress reporting…")
  flush.console()
  
  combined_sf <- NULL
  kept <- 0
  
  for (i in seq_along(files)) {
    sf_i <- readRDS(files[i])
    
    if (nrow(sf_i) == 0) next
    
    if (is.null(combined_sf)) {
      combined_sf <- sf_i
    } else {
      combined_sf <- rbind(combined_sf, sf_i)
    }
    
    kept <- kept + nrow(sf_i)
    
    if (i %% progress_every == 0 || i == length(files)) {
      cat(
        sprintf(
          "  [%3d / %3d] chunks processed | %6d pixels retained",
          i, length(files), kept
        )
      )
    }
  }
  flush.console()
  
  if (is.null(combined_sf)) {
    cat("No intersecting pixels found — nothing written.")
    flush.console()
    return(invisible(NULL))
  }
  
  out_path <- file.path(
    site_dir,
    paste0(site_id, "_evi_timeseries_fetch_sf.rds")
  )
  
  saveRDS(combined_sf, out_path)
  cat("✅ Wrote combined fetch sf: ", out_path)
  flush.console()
    
  invisible(combined_sf)
}