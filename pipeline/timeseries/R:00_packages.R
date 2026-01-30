# R/00_packages.R
# Centralized package loading (so your Rmd stays clean)

load_wetlsp_packages <- function() {
  pkgs <- c(
    "sf", "terra", "dplyr", "magick",
    "tidyverse", "scales", "patchwork", "zoo",
    "rjson", "geojsonR", "ncdf4", "arrow"
  )
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Missing package: ", p, ". Install it first.")
    }
  }
  invisible(TRUE)
}