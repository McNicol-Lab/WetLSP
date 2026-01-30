combine_fetch_sf_progress_split4 <- function(
    fetch_dir,
    site_dir,
    site_id,
    radius_m,
    pattern = "^chunk_\\d{3}_evi_sf\\.rds$",
    n_parts = 4L,
    progress_every = 10L
) {
  files <- sort(list.files(fetch_dir, pattern = pattern, full.names = TRUE))
  if (length(files) == 0) stop("No fetch chunk files found in: ", fetch_dir)
  
  message("Found ", length(files), " fetch chunks")
  part_ids <- (seq_along(files) - 1L) %% n_parts
  
  out_paths <- character(n_parts)
  
  for (p in 0:(n_parts - 1L)) {
    message("\n--- Part ", p + 1L, "/", n_parts, " ---")
    files_p <- files[part_ids == p]
    
    combined_sf <- NULL
    kept <- 0L
    
    for (i in seq_along(files_p)) {
      sf_i <- readRDS(files_p[i])
      if (nrow(sf_i) == 0) next
      
      if (is.null(combined_sf)) combined_sf <- sf_i else combined_sf <- rbind(combined_sf, sf_i)
      kept <- kept + nrow(sf_i)
      
      if (i %% progress_every == 0L || i == length(files_p)) {
        message(sprintf("  [%3d/%3d] files | %8d pixels", i, length(files_p), kept))
      }
    }
    
    out_path <- file.path(
      fetch_dir,
      sprintf("%s-evi-timeseries-fetch-sf-%sm-part%02dof%02d.rds",
              site_id, radius_m, p + 1L, n_parts)
    )
    saveRDS(combined_sf, out_path)
    message("✅ Wrote part: ", out_path)
    
    out_paths[p + 1L] <- out_path
    rm(combined_sf); gc()
  }
  
  invisible(out_paths)
}