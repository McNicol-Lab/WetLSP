# Helper: CONUS map with Alaska inset (flux-tower ready)
# - df must have: lon, lat, site_id
# - returns a ggplot object

library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(patchwork)
library(ggrepel)

plot_us_with_ak_inset <- function(
    df,
    point_size = 2,
    point_alpha = 0.8,
    label_size = 4,
    label_nudge_y = 1,
    conus_xlim = c(-125, -66),
    conus_ylim = c(24, 50),
    ak_xlim = c(-170, -130),
    ak_ylim = c(51, 72),
    ak_inset_scale = 0.5,      # size of inset relative to full plot
    ak_inset_pos = c(0.02, 0.68) # (x,y) in npc units: bottom-left of inset
) {
  stopifnot(all(c("lon","lat","site_id") %in% names(df)))
  
  # US states (sf)
  us <- rnaturalearth::ne_states(
    country = "United States of America",
    returnclass = "sf"
  ) |>
    st_transform(4326)
  
  conus_states <- us |>
    filter(!name_en %in% c("Alaska", "Hawaii", "Puerto Rico"))
  
  ak_state <- us |>
    filter(name_en == "Alaska")
  
  # Split tower points
  conus_pts <- df |>
    filter(.data$lon >= conus_xlim[1], .data$lon <= conus_xlim[2],
           .data$lat >= conus_ylim[1], .data$lat <= conus_ylim[2])
  
  ak_pts <- df |>
    filter(.data$lon >= ak_xlim[1], .data$lon <= ak_xlim[2],
           .data$lat >= ak_ylim[1], .data$lat <= ak_ylim[2])
  
  # --- main CONUS plot ---
  p_conus <- ggplot() +
    geom_sf(data = conus_states, fill = "grey95", color = "grey70", linewidth = 0.25) +
    geom_point(
      data = conus_pts,
      aes(x = lon, y = lat, color = site_id),
      size = point_size,
      alpha = point_alpha
    ) +
    ggrepel::geom_text_repel(
      data = conus_pts,
      aes(x = lon, y = lat, label = site_id),
      size = label_size,
      box.padding = 0.3,
      point.padding = 0.2,
      force = 2,
      force_pull = 0.5,
      direction = "both",
      min.segment.length = 0,
      segment.color = NA,
      seed = 1
    ) +
    coord_sf(xlim = conus_xlim, ylim = conus_ylim, expand = FALSE) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank()
    ) + scale_color_brewer(palette = "Set1")
  
  # --- Alaska inset ---
  p_ak <- ggplot() +
    geom_sf(data = ak_state, fill = "grey95", color = "grey70", linewidth = 0.25) +
    geom_point(
      data = ak_pts,
      aes(x = lon, y = lat, color = site_id),
      size = point_size,
      shape = 16,
      alpha = point_alpha
    ) +
    ggrepel::geom_text_repel(
      data = ak_pts,
      aes(x = lon, y = lat, label = site_id),
      size = label_size,
      box.padding = 0.3,
      point.padding = 0.2,
      force = 2,
      force_pull = 0.5,
      direction = "both",
      min.segment.length = 0,
      segment.color = NA,
      seed = 1
    ) +
    coord_sf(xlim = ak_xlim, ylim = ak_ylim, expand = FALSE) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = "grey60", linewidth = 0.4)
    ) + scale_color_brewer(palette = "Set1")
  
  # Inset via patchwork (overlay)
  # We draw CONUS then place Alaska in the bottom-left (tunable)
  p_conus +
    inset_element(
      p_ak,
      left   = ak_inset_pos[1],
      bottom = ak_inset_pos[2],
      right  = ak_inset_pos[1] + ak_inset_scale,
      top    = ak_inset_pos[2] + ak_inset_scale
    )
}

# ---- example usage ----
# df <- tibble(site_id = c("US-XYZ","AK-ABC"), lon = c(-90, -150), lat = c(41, 64))
# plot_us_with_ak_inset(df)