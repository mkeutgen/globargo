library(tidyverse)
library(ggplot2)
library(sf)
library(viridis)

### --- 5° Grid Map --- ###
stipple_resolution <- 5


# 1. Aggregate Argo float locations to 5° bins for the season
argo_bins <- df_argo_clean %>%
  mutate(
    lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
    lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
  ) %>%
  group_by(lon_bin_stipple, lat_bin_stipple) %>%
  summarize(count = n(), .groups = "drop")

# Create a full grid using lower left corners
full_grid <- expand.grid(
  lon_bin_stipple = seq(-180, 180 - stipple_resolution, by = stipple_resolution),
  lat_bin_stipple = seq(-90, 90 - stipple_resolution, by = stipple_resolution)
) %>% as_tibble()

# Left join and replace NAs with 0; add discrete bins for the fill scale
argo_bins_full <- full_grid %>%
  left_join(argo_bins, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
  mutate(
    count = ifelse(is.na(count), 0, count),
    count_bin = cut(count,
                    breaks = c(-Inf, 1, 5, 25, 100,500, Inf),
                    labels = c("0", "1-5", "5-25", "25-100", "100-500","500+")),
    # Compute center coordinates for geom_tile (tiles are centered at these points)
    lon_center = lon_bin_stipple + stipple_resolution / 2,
    lat_center = lat_bin_stipple + stipple_resolution / 2
  )

# Identify sampled areas (cells with at least 1 profile)
covered <- argo_bins_full %>% filter(count > 0)

nrow(covered)/(71/100*nrow(argo_bins_full)) # 0.6509737
 
# Plotting with geom_tile, stippled points, and fixed axis scales
map_5degree <- ggplot() +
  # Tile layer with discrete fill
  geom_tile(data = argo_bins_full, 
            aes(x = lon_center, y = lat_center, fill = count_bin),
            width = stipple_resolution, height = stipple_resolution,
            alpha = 0.95) +
  scale_fill_viridis_d(
    option = "viridis",  
    name = "Number of Argo profiles",
    guide = guide_legend(reverse = TRUE, na.translate = FALSE)
  ) +
  # Overlay the world map (assuming 'world' is an sf object)
  geom_sf(data = world, fill = "white", color = "white", inherit.aes = FALSE) +
  # Fix the axes: properly label longitude and latitude with tick marks every 30°
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid = element_blank()) +
  labs(title = "126,591 of Argo Profiles Aggregated in 5° Grid",
       subtitle = "66 % of 5° gridcells (~ planetary scale) are sampled",x="Longitude",y="Latitude")

ggsave(map_5degree, filename = "figures/map5deg.png",width = 14, height = 10, dpi = 300)


### --- 1° Grid Map --- ###
library(tidyverse)
library(ggplot2)
library(sf)
library(viridis)

### --- 1° Grid Map --- ###
stipple_resolution <- 1

# 1. Aggregate Argo float locations to 5° bins for the season
argo_bins <- df_argo_clean %>%
  mutate(
    lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
    lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
  ) %>%
  group_by(lon_bin_stipple, lat_bin_stipple) %>%
  summarize(count = n(), .groups = "drop")

# Create a full grid using lower left corners
full_grid <- expand.grid(
  lon_bin_stipple = seq(-180, 180 - stipple_resolution, by = stipple_resolution),
  lat_bin_stipple = seq(-90, 90 - stipple_resolution, by = stipple_resolution)
) %>% as_tibble()

# Left join and replace NAs with 0; add discrete bins for the fill scale
argo_bins_full <- full_grid %>%
  left_join(argo_bins, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
  mutate(
    count = ifelse(is.na(count), 0, count),
    count_bin = cut(count,
                    breaks = c(-Inf, 1, 5,10, 25,100, Inf),
                    labels = c("0", "1-5","5-10", "10-25", "25-100","100+")),
    # Compute center coordinates for geom_tile (tiles are centered at these points)
    lon_center = lon_bin_stipple + stipple_resolution / 2,
    lat_center = lat_bin_stipple + stipple_resolution / 2
  )
# Identify sampled areas (cells with at least 1 profile)
covered <- argo_bins_full %>% filter(count > 0)

nrow(covered)/(71/100*nrow(argo_bins_full)) # 0.32



# Plotting with geom_tile, stippled points, and fixed axis scales
map_1degree <- ggplot() +
  # Tile layer with discrete fill
  geom_tile(data = argo_bins_full, 
            aes(x = lon_center, y = lat_center, fill = count_bin),
            width = stipple_resolution, height = stipple_resolution,
            alpha = 0.95) +
  scale_fill_viridis_d(
    option = "viridis",  
    name = "Number of Argo profiles",
    guide = guide_legend(reverse = TRUE, na.translate = FALSE)
  ) +
  # Overlay the world map (assuming 'world' is an sf object)
  geom_sf(data = world, fill = "white", color = "white", inherit.aes = FALSE) +
  # Fix the axes: properly label longitude and latitude with tick marks every 30°
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid = element_blank()) +
  labs(title = "126,591 Argo Profiles Aggregated in 1° Grid",x="Longitude",y="Latitude",
       subtitle = "33 % of 1° (~ mesoscale) gridcells are sampled")

ggsave(map_1degree, filename = "figures/map1deg.png",width = 14, height = 10, dpi = 300)

### --- 0.25° Grid Map --- ###
stipple_resolution <- 0.25

# 1. Aggregate Argo float locations to 5° bins for the season
argo_bins <- df_argo_clean %>%
  mutate(
    lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
    lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
  ) %>%
  group_by(lon_bin_stipple, lat_bin_stipple) %>%
  summarize(count = n(), .groups = "drop")

# Create a full grid using lower left corners
full_grid <- expand.grid(
  lon_bin_stipple = seq(-180, 180 - stipple_resolution, by = stipple_resolution),
  lat_bin_stipple = seq(-90, 90 - stipple_resolution, by = stipple_resolution)
) %>% as_tibble()

# Left join and replace NAs with 0; add discrete bins for the fill scale
argo_bins_full <- full_grid %>%
  left_join(argo_bins, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
  mutate(
    count = ifelse(is.na(count), 0, count),
    count_bin = cut(count,
                    breaks = c(-Inf, 1, 5, 25, 100,500, Inf),
                    labels = c("0", "1-5", "5-25", "25-100", "100-500","500+")),
    # Compute center coordinates for geom_tile (tiles are centered at these points)
    lon_center = lon_bin_stipple + stipple_resolution / 2,
    lat_center = lat_bin_stipple + stipple_resolution / 2
  )

# Identify sampled areas (cells with at least 1 profile)
covered <- argo_bins_full %>% filter(count > 0)

nrow(covered)/(71/100*nrow(argo_bins_full)) # 0.07


# Plotting with geom_tile, stippled points, and fixed axis scales
map_025degree <- ggplot() +
  # Tile layer with discrete fill
  geom_tile(data = argo_bins_full, 
            aes(x = lon_center, y = lat_center, fill = count_bin),
            width = stipple_resolution, height = stipple_resolution,
            alpha = 0.95) +
  scale_fill_viridis_d(
    option = "viridis",  
    name = "Number of Argo profiles",
    guide = guide_legend(reverse = TRUE, na.translate = FALSE)
  ) +
  # Overlay the world map (assuming 'world' is an sf object)
  geom_sf(data = world, fill = "white", color = "white", inherit.aes = FALSE) +
  # Fix the axes: properly label longitude and latitude with tick marks every 30°
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid = element_blank()) +
  labs(title = "126,591Argo Profiles Aggregated in 0.25° Grid",x="Longitude",y="Latitude",
       subtitle = "7 % of 0.25° (~ submesoscale) gridcells are sampled")

ggsave(map_025degree, filename = "figures/map025deg.png",width = 14, height = 10, dpi = 300)

