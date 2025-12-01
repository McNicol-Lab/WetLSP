## ✅ Helper 1: Simple pixel-time-series viewer (choose pixel by index)


ViewPixelTimeseries <- function(site_sf, cell_vec) {
  if (!("cell" %in% names(site_sf)))
    stop("sf object must have a 'cell' column")
  
  # Filter only pixels that exist
  pix <- site_sf[site_sf$cell %in% cell_vec, ]
  if (nrow(pix) == 0) stop("None of the requested cells were found")
  
  # Build combined dataframe for all cells
  df_list <- lapply(seq_len(nrow(pix)), function(i) {
    tibble::tibble(
      cell  = pix$cell[i],
      date  = c(pix$dates_raw[[i]],    pix$dates_spline[[i]]),
      evi   = c(pix$evi_raw[[i]],       pix$evi_spline[[i]]),
      type  = c(rep("raw",    length(pix$dates_raw[[i]])),
                rep("smooth", length(pix$dates_spline[[i]])))
    )
  })
  
  df <- dplyr::bind_rows(df_list)
  
  # Nice contrasting colors for cells
  n_cells <- length(unique(df$cell))
  cell_colors <- scales::hue_pal()(n_cells)
  
  ggplot2::ggplot(df, ggplot2::aes(x = as.Date(date), y = evi,
                                   color = as.factor(cell),
                                   shape = type,
                                   linetype = type)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    # ggplot2::geom_line(data = df[df$type == "smooth", ],
                       # size = 0.8, alpha = 0.9) +
    ggplot2::scale_color_manual(values = cell_colors) +
    ggplot2::scale_shape_manual(values = c(raw = 16, smooth = 1)) +
    # ggplot2::scale_linetype_manual(values = c(raw = "blank", smooth = "solid")) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = "EVI Time Series for Cells: ",
      color = "Cell #",
      shape = "Type",
      linetype = "Type",
      x = "Date",
      y = "EVI"
    )
}

## ✅ Helper 2: Interactive click-to-select pixel (optional)

PickPixelInteractive <- function(site_sf) {
  library(mapview)
  library(mapedit)
  
  m <- mapview::mapview(site_sf, zcol = "cell", legend = FALSE)
  sel <- mapedit::editMap(m)
  
  if (length(sel$features) == 0) {
    message("No pixel selected.")
    return(NULL)
  }
  
  chosen_id <- sel$features[[1]]$properties$cell
  message("Selected pixel: ", chosen_id)
  
  ViewPixelTimeseries(site_sf, chosen_id)
}