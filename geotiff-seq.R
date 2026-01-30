library(terra)
library(sf)
library(ggplot2)
library(magick)
library(stringr)
library(viridisLite)
library(scales)
library(dplyr)
library(patchwork)
library(zoo)  # you use zoo::rollapply and zoo::na.locf

make_inset <- function(day) {
  ggplot(evi_long, aes(date, mean_evi)) +
    geom_line(color = "grey30", linewidth = 0.4) +
    geom_vline(xintercept = day, color = "red", linewidth = 0.5) +   # <- Date, not numeric
    scale_x_date(expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, max(evi_long$mean_evi, na.rm = TRUE))) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 7) +
    theme(
      panel.grid  = element_blank(),
      axis.text   = element_blank(),
      axis.ticks  = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    )
}

tif_to_gg_magick <- function(f, site_id, year, roll_lookup, upper_floor = 0.12) {
  r <- terra::rast(f)
  
  df <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(df)[3] <- "evi"
  
  # Parse day
  dstr <- stringr::str_match(basename(f), "EVI_(\\d{8})")[, 2]
  day  <- as.Date(dstr, format = "%Y%m%d")
  
  # Clip EVI to [0, 1] first (keeps outliers from wrecking scaling)
  df$evi <- pmin(pmax(df$evi, 0), 1)
  
  # Rolling upper limit for this frame (± window defined upstream)
  upper <- unname(roll_lookup[as.character(day)])
  if (!is.finite(upper)) {
    # fallback if roll_lookup missing (edges / NA)
    upper <- terra::global(r, "max", na.rm = TRUE)[1, 1]
  }
  upper <- max(upper, upper_floor, na.rm = TRUE)
  
  # Midpoint: use median of "signal" pixels (avoid near-zero winter domination)
  vals <- df$evi[is.finite(df$evi) & !is.na(df$evi)]
  vals <- vals[vals <= upper]
  vals_sig <- vals[vals > 0.02]  # tune threshold if needed
  
  mid <- if (length(vals_sig) >= 10) {
    as.numeric(stats::median(vals_sig, na.rm = TRUE))
  } else if (length(vals) >= 10) {
    as.numeric(stats::median(vals, na.rm = TRUE))
  } else {
    0.35 * upper
  }
  mid <- max(0, min(mid, upper))
  
  p_map <- ggplot(df, aes(x = x, y = y, fill = evi)) +
    geom_raster() +
    coord_equal() +
    scale_fill_gradientn(
      colours = c(
        "#4B1D3F",  # dark brown–purple (low / dormant)
        "#F1E05A",  # yellow (transition)
        "#1A9850"   # green (peak growth)
      ),
      limits = c(0, upper),
      values = scales::rescale(c(0, mid, upper), to = c(0, 1)),
      oob = scales::squish,
      na.value = "transparent",
      name = "EVI"
    ) +
    labs(
      title    = paste0(site_id, " EVI ", format(day, "%Y-%m-%d")),
      subtitle = paste0("Scale: 0–rolling max (± window). Midpoint = ", signif(mid, 3)),
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  p_inset <- make_inset(day)
  
  p_final <- p_map +
    inset_element(
      p_inset,
      left = 0.75, bottom = 0.75, right = 0.98, top = 0.98
    )
  
  png_file <- tempfile(fileext = ".png")
  ggsave(png_file, plot = p_final, width = 5, height = 4, dpi = 120, bg = "white")  # <- p_final
  magick::image_read(png_file)
}

roll_lookup <- setNames(roll15_max, as.character(days))

imgs <- lapply(
  tifs,
  tif_to_gg_magick,
  site_id = site_id,
  year = year,
  roll_lookup = roll_lookup
)

gif <- magick::image_animate(magick::image_join(imgs), fps = 5)
gif_path <- file.path(site_dir, paste0(site_id, "_EVI_", year, "_thin", thin, "_ggplot.gif"))
magick::image_write(gif, gif_path)

gif_path