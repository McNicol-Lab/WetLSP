stitch_fetch_sf_parts <- function(
    part_paths,
    out_path,
    progress = TRUE
) {
  stopifnot(length(part_paths) > 1)
  
  if (progress) {
    message("Stitching ", length(part_paths), " parts:")
    for (p in part_paths) message("  - ", basename(p))
  }
  
  sf_all <- NULL
  
  for (i in seq_along(part_paths)) {
    sf_i <- readRDS(part_paths[i])
    
    if (is.null(sf_all)) {
      sf_all <- sf_i
    } else {
      sf_all <- rbind(sf_all, sf_i)
    }
    
    if (progress) {
      message("  merged ", i, "/", length(part_paths),
              " | rows so far: ", nrow(sf_all))
    }
  }
  
  saveRDS(sf_all, out_path)
  message("✅ Wrote stitched sf to: ", out_path)
  
  invisible(out_path)
}