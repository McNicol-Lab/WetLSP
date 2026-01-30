# R/40_gif_ggplot.R
# Render tifs -> ggplot frames -> animated gif
# Includes: (a) rolling max scale (±4 days by default), (b) 3-color ramp, (c) inset annual mean curve

build_evi_daily_mean_curve <- function(sf_all) {
  tibble::tibble(
    date = unlist(sf_all$dates_spline),
    evi  = unlist(sf_all$evi_spline)
  ) |>
    dplyr::filter(!is.na(evi)) |>
    dplyr::group_by(date) |>
    dplyr::summarise(mean_evi = mean(evi), .groups = "drop")
}

make_inset_plot <- function(evi_long, day) {
  ggplot2::ggplot(evi_long, ggplot2::aes(date, mean_evi)) +
    ggplot2::geom_line(color = "grey30", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = as.numeric(day), color = "red", linewidth = 0.5) +
    ggplot2::scale_y_continuous(limits = c(0, max(evi_long$mean_evi, na.rm = TRUE))) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
}

get_day_from_tif <- function(f) {
  dstr <- stringr::str_match(basename(f), "EVI_(\\d{8})")[, 2]
  as.Date(dstr, format = "%Y%m%d")
}

compute_rolling_upper <- function(tifs, days, half_window = 4, upper_floor = 0.12) {
  # per-frame max
  file_max <- vapply(tifs, function(f) {
    r <- terra::rast(f)
    terra::global(r, "max", na.rm = TRUE)[1, 1]
  }, numeric(1))
  
  width <- 2 * half_window + 1
  
  roll_max <- zoo::rollapply(
    file_max,
    width = width,
    FUN = max,
    align = "center",
    fill = NA_real_,
    na.rm = TRUE
  )
  # fill edges
  roll_max <- zoo::na.locf(roll_max, na.rm = FALSE)
  roll_max <- zoo::na.locf(roll_max, fromLast = TRUE, na.rm = FALSE)
  
  # floor prevents winter collapse
  roll_max <- pmax(roll_max, upper_floor)
  
  setNames(roll_max, as.character(days))
}

tif_to_gg_magick <- function(f, site_id, year, upper, mid, evi_long, inset = TRUE) {
  r <- terra::rast(f)
  df <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(df)[3] <- "evi"
  
  day <- get_day_from_tif(f)
  
  # clip low end only (keep outliers from breaking scale via 'upper')
  df$evi <- pmax(df$evi, 0)
  
  # 3-color seasonal ramp
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = evi)) +
    ggplot2::geom_raster() +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_gradientn(
      colours = c("#4B1D3F", "#F1E05A", "#1A9850"),
      limits  = c(0, upper),
      values  = scales::rescale(c(0, mid, upper), to = c(0, 1)),
      oob     = scales::squish,
      na.value = "transparent",
      name = "EVI"
    ) +
    ggplot2::labs(
      title = paste0(site_id, " EVI ", format(day, "%Y-%m-%d")),
      subtitle = paste0("Scale: 0–rolling max (±window). upper=", signif(upper, 3),
                        " mid=", signif(mid, 3)),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
  
  if (inset) {
    p_inset <- make_inset_plot(evi_long, day)
    p <- p + patchwork::inset_element(
      p_inset, left = 0.62, bottom = 0.62, right = 0.98, top = 0.98
    )
  }
  
  png_file <- tempfile(fileext = ".png")
  ggplot2::ggsave(png_file, plot = p, width = 9, height = 7, dpi = 120)
  magick::image_read(png_file)
}

render_evi_gif_ggplot <- function(
    site_id,
    year,
    thin,
    radius_m,
    site_dir = file.path("output", site_id),
    out_dir  = file.path(site_dir, paste0("evi_daily_", year, "_r", radius_m, "m")),
    sf_all_path = file.path(site_dir, paste0(site_id, "_evi_timeseries_fetch_sf_100m.rds")),
    half_window_days = 4,
    upper_floor = 0.12,
    fps = 5
) {
  load_wetlsp_packages()
  
  if (!file.exists(sf_all_path)) stop("sf_all_path not found: ", sf_all_path)
  sf_all <- readRDS(sf_all_path)
  sf_all$dates_spline <- coerce_dates_list_to_Date(sf_all$dates_spline)
  
  tifs <- sort(list.files(out_dir, pattern = "\\.tif$", full.names = TRUE))
  if (length(tifs) == 0) stop("No .tif files found in: ", out_dir)
  
  days <- as.Date(vapply(tifs, get_day_from_tif, as.Date("1970-01-01")))
  ord <- order(days)
  days <- days[ord]
  tifs <- tifs[ord]
  
  # rolling upper
  roll_lookup <- compute_rolling_upper(tifs, days, half_window = half_window_days, upper_floor = upper_floor)
  
  # annual mean curve for inset
  evi_long <- build_evi_daily_mean_curve(sf_all)
  
  # build frames
  imgs <- vector("list", length(tifs))
  for (i in seq_along(tifs)) {
    day <- days[i]
    upper <- as.numeric(roll_lookup[[as.character(day)]])
    # midpoint: median of values on that day but bounded (prevents weird behavior)
    mid <- 0.35 * upper
    
    imgs[[i]] <- tif_to_gg_magick(
      f = tifs[i],
      site_id = site_id,
      year = year,
      upper = upper,
      mid = mid,
      evi_long = evi_long,
      inset = TRUE
    )
    
    if (i %% 25 == 0) message("Rendered ", i, "/", length(tifs))
  }
  
  gif <- magick::image_animate(magick::image_join(imgs), fps = fps)
  gif_path <- file.path(site_dir, paste0(site_id, "_EVI_", year, "_thin", thin, "_ggplot.gif"))
  magick::image_write(gif, gif_path)
  
  gif_path
}