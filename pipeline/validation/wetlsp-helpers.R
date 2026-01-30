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
                rep("gap-filled", length(pix$dates_spline[[i]])))
    )
  })
  
  df <- dplyr::bind_rows(df_list)
  
  # Nice contrasting colors for cells
  n_cells <- length(unique(df$cell))
  cell_colors <- scales::hue_pal()(n_cells)
  
  ggplot2::ggplot(df, ggplot2::aes(x = as.Date(date), y = evi,
                                   color = as.factor(cell),
                                   # shape = as.factor(cell),
                                   linetype = type,
                                   alpha = type)) +
    ggplot2::geom_point() +
    # ggplot2::geom_line(data = df[df$type == "smooth", ],
                       # size = 0.8, alpha = 0.9) +
    ggplot2::scale_x_date(limits = c(ymd("2021-01-01"), ("2024-12-31"))) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_color_manual(values = c(cell_colors)) +
    ggplot2::scale_alpha_manual(values = c(0.1, 1)) +
    # ggplot2::scale_shape_manual(values = c(raw = 16, smooth = 1)) +
    # ggplot2::scale_linetype_manual(values = c(raw = "blank", smooth = "solid")) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = "EVI Time Series for Cells: ",
      color = "Type",
      shape = "Type",
      alpha = "Type",
      linetype = "Type",
      x = "Date",
      y = "EVI"
    )
}

## Helper 2: Display map of pixels in fetch with option to pick and highlight a cell

PickPixelInteractive <- function(site_sf, cell_id = NULL) {
  stopifnot("cell" %in% names(site_sf))
  
  library(mapview)
  library(mapedit)
  library(sf)
  
  # Base layer: all pixels
  base <- mapview(
    site_sf,
    zcol   = "cell",
    legend = FALSE,
    alpha  = 0.25,
    cex    = 0.08
  )
  
  # Optional highlighted pixel
  if (!is.null(cell_id)) {
    cell_id <- as.integer(cell_id)
    hit <- site_sf[site_sf$cell == cell_id, ]
    
    if (nrow(hit) == 0) {
      stop("cell_id not found in site_sf: ", cell_id)
    }
    
    hi <- mapview(
      hit,
      col.regions = "red",
      alpha = 1,
      cex   = 0.9,
      legend = FALSE,
      bbox = st_bbox(hit)   # <-- THIS is the zoom
    )
    
    m <- base + hi
  } else {
    m <- base
  }
  
  sel <- mapedit::editMap(m)
  
  if (length(sel$features) == 0) {
    message("No pixel selected.")
    return(invisible(NULL))
  }
  
  chosen_id <- sel$features[[1]]$properties$cell
  message("Selected pixel: ", chosen_id)
  
  ViewPixelTimeseries(site_sf, chosen_id)
}

## Helper 3: Violin plots

ViewPixelViolin <- function(site_sf, cell_vec) {
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
                rep("gap-filled", length(pix$dates_spline[[i]])))
    )
  })
  
  df <- dplyr::bind_rows(df_list)
  
  # Nice contrasting colors for cells
  n_cells <- length(unique(df$cell))
  cell_colors <- scales::hue_pal()(n_cells)
  
  ggplot2::ggplot(df %>% filter(type == "gap-filled"), 
                  ggplot2::aes(x = cell, y = evi,
                                   fill = as.factor(cell),
                                   # shape = as.factor(cell),
                                   # linetype = type,
                                   # alpha = type
                               )) +
    # ggplot2::geom_point() +
    ggplot2::geom_violin() +
    # ggplot2::geom_line(data = df[df$type == "smooth", ],
    # size = 0.8, alpha = 0.9) +
    # ggplot2::scale_x_date(limits = c(ymd("2021-01-01"), ("2024-12-31"))) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_color_manual(values = c(cell_colors)) +
    # ggplot2::scale_alpha_manual(values = c(0.1, 1)) +
    # ggplot2::scale_shape_manual(values = c(raw = 16, smooth = 1)) +
    # ggplot2::scale_linetype_manual(values = c(raw = "blank", smooth = "solid")) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = "EVI Time Series for Cells: ",
      color = "Type",
      shape = "Type",
      alpha = "Type",
      linetype = "Type",
      x = "Cell",
      y = "EVI"
    )
}

## Helper 4: Boxplots

ViewPixelBoxplot <- function(site_sf, cell_vec) {
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
                rep("gap-filled", length(pix$dates_spline[[i]])))
    )
  })
  
  df <- dplyr::bind_rows(df_list)
  
  # Nice contrasting colors for cells
  n_cells <- length(unique(df$cell))
  cell_colors <- scales::hue_pal()(n_cells)
  
  ggplot2::ggplot(df %>% filter(type == "gap-filled"), 
                  ggplot2::aes(x = as.factor(cell), y = evi,
                               fill = as.factor(cell),
                               # shape = as.factor(cell),
                               # linetype = type,
                               # alpha = type
                  )) +
    # ggplot2::geom_point() +
    ggplot2::geom_boxplot() +
    # ggplot2::geom_line(data = df[df$type == "smooth", ],
    # size = 0.8, alpha = 0.9) +
    # ggplot2::scale_x_date(limits = c(ymd("2021-01-01"), ("2024-12-31"))) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_color_manual(values = c(cell_colors)) +
    # ggplot2::scale_alpha_manual(values = c(0.1, 1)) +
    # ggplot2::scale_shape_manual(values = c(raw = 16, smooth = 1)) +
    # ggplot2::scale_linetype_manual(values = c(raw = "blank", smooth = "solid")) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = "EVI Time Series for Cells: ",
      color = "Type",
      shape = "Type",
      alpha = "Type",
      linetype = "Type",
      x = "Cell",
      y = "EVI"
    )
}