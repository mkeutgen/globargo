# Load necessary libraries
library(tidyverse)
library(robustbase)
library(gsw)
library(zoo)
library(oce)
library(ggpubr)
library(segmented)
library(dplyr)
library(conflicted)
library(pracma)
library(fs)
library(ggnewscale)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggspatial)
library(viridis)
library(sp)
library(spdep)
library(mgcv)
library(patchwork)
library(scales)
library(fitdistrplus)
library(poweRlaw)
library(ggplot2)
library(tidync)

# Resolve function conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Load Qiu's data
Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")
setwd("/data/GLOBARGO/src/")
df_argo_clean <- read_csv("data/df_argo_loc.csv")
df_argo_clean$TIME %>% min(na.rm = T) # "2010-05-30"
df_argo_clean$TIME %>% max(na.rm = T) # "2024-05-02"

df_complete_clean <- read_csv("data/df_eddy_subduction_anom.csv")
df_carbon_clean <- read_csv("data/df_carbon_subduction_anom.csv")
df_eke <- read_table("data/Qiu2018_data_fig5.asc")

df_eke$LON_num %>% unique() %>% sort() # 3° res

breaks_eke_bl <- c(0, 1e-4, 1e-3, 2.5e-3, 5e-3, 7.5e-3, 1e-2, 2.5e-2, 5e-2, 7.5e-2, 1e-1, 1)
breaks_eke_unbl <- c(0, 1e-4, 1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3, 1e-2, 2.5e-2, 1e-1)


# Produce maps of EKE : 
df_eke$LON <- df_eke$`lon(E)` 
df_eke$LAT <- df_eke$`lat(N)`
df_eke$EKE_BL <- df_eke$`eke(bl)(m/s)^2` %>%  as.numeric()
df_eke$EKE_UNBL <- df_eke$`eke(unbl)(m/s)^2` %>%  as.numeric()

df_eke <- df_eke %>% select(LON,LAT,EKE_BL,EKE_UNBL)



df_eke <- df_eke %>%
  mutate(
    LON = ifelse(as.numeric(LON) > 180, as.numeric(LON) - 360, as.numeric(LON))
  )

df_eke$EKE_BL 



df_eke <- df_eke %>%
  mutate(
    EKE_BL_bin = cut(EKE_BL,
                     breaks = breaks_eke_bl,
                     include.lowest = TRUE),
    EKE_UNBL_bin = cut(EKE_UNBL,
                      breaks = breaks_eke_unbl,
                      include.lowest = TRUE),
    EKE_UNBL_bin_bl = cut(EKE_UNBL,
                          breaks = breaks_eke_unbl,
                          include.lowest = TRUE)
  )





# 1) Balanced EKE map

eke_balanced_map <- ggplot() +
  geom_rect(
    aes(xmin = -180, xmax = 180, ymin = 53, ymax = 90),
    fill = "lightblue",
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin = -180, xmax = 180, ymin = -90, ymax = -67),
    fill = "lightblue",
    inherit.aes = FALSE
  )+
  geom_tile(
    data = df_eke,
    aes(x = as.numeric(LON), 
        y = as.numeric(LAT), 
        fill = EKE_BL_bin)
  ) +
  scale_fill_viridis_d() + 
  geom_sf(data = world, fill = "grey", color = "grey", inherit.aes = FALSE) +
  coord_sf(
    xlim = c(-180, 180),
    ylim = c(-90, 90),
    expand = FALSE
  ) +
  labs(
    title = "Balanced Eddy Kinetic Energy",
    x = "Longitude",
    y = "Latitude",
    fill = "EKE (BL)"
  ) +
  theme_minimal() 



# 2) Unbalanced EKE map
eke_unbalanced_map <- ggplot() +
  geom_rect(
    aes(xmin = -180, xmax = 180, ymin = 53, ymax = 90),
    fill = "lightblue",
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin = -180, xmax = 180, ymin = -90, ymax = -67),
    fill = "lightblue",
    inherit.aes = FALSE
  )+
  geom_tile(
    data = df_eke,
    aes(x = as.numeric(LON), 
        y = as.numeric(LAT), 
        fill = EKE_UNBL_bin)
  ) +
  geom_contour(
    data = df_eke,
    aes(x = as.numeric(LON), 
        y = as.numeric(LAT), 
        z = EKE_UNBL),
    color = "white", alpha = 0.3
  ) +
  scale_fill_viridis_d(
    option = "viridis") +
  geom_sf(data = world, fill = "grey", color = "grey", inherit.aes = FALSE) +
  coord_sf(
    xlim = c(-180, 180),
    ylim = c(-90, 90),
    expand = FALSE
  ) +
  labs(
    title = "Unbalanced Eddy Kinetic Energy",
    x = "Longitude",
    y = "Latitude",
    fill = "EKE (BL)"
  ) +
  theme_minimal() 


# Finally display
print(eke_balanced_map)
print(eke_unbalanced_map)

# Smoothing with GAMs
# Convert LON, LAT to numeric and ensure positivity of EKE
df_eke <- df_eke %>%
  mutate(
    LON_num = as.numeric(LON),
    LAT_num = as.numeric(LAT),
    EKE_BL_pos = ifelse(EKE_BL <= 0, NA_real_, EKE_BL),
    EKE_UNBL_pos = ifelse(EKE_UNBL <= 0, NA_real_, EKE_UNBL),
    logEKE_BL   = log(EKE_BL_pos),
    logEKE_UNBL = log(EKE_UNBL_pos)
  ) %>%
  # filter out any rows with NA (which includes any EKE <= 0)
  filter(!is.na(LON_num), !is.na(LAT_num),
         !is.na(logEKE_BL), !is.na(logEKE_UNBL))



# Gaussian on log(EKE_BL)
gam_bl_gaussian_log <- gam(
  logEKE_BL ~ s(LAT_num,LON_num, bs = "sos", k = 600),
  data   = df_eke,
  family = gaussian(),
  method = "REML"
)


gam_bl_gaussian_log %>% summary() #R-sq.(adj) =  0.913   Deviance explained = 92.3%

gam_unbl_gaussian_log <- gam(
  logEKE_UNBL ~ s(LAT_num,LON_num, bs = "sos", k = 600),
  data   = df_eke,
  family = gaussian(),
  method = "REML"
)


# Suppose we define a function to make a simple lat-lon grid:
make_prediction_grid <- function(df, lon_step = 0.25, lat_step = 0.25) {
  lon_seq <- seq(-180,180,
                 by = lon_step)
  lat_seq <- seq(-90,
                 90,
                 by = lat_step)
  expand.grid(LON_num = lon_seq, LAT_num = lat_seq)
}

pred_grid <- make_prediction_grid(df_eke)

pred_grid$EKE_bl_hat_gaussian_log <- exp(
  predict(gam_bl_gaussian_log, newdata = pred_grid, type = "link")
)

pred_grid$EKE_unbl_hat_gaussian_log <- exp(
  predict(gam_unbl_gaussian_log, newdata = pred_grid, type = "link")
)



pred_grid$EKE_bl_hat_binned <- cut(pred_grid$EKE_bl_hat_gaussian_log,
                                   breaks = breaks_eke_bl,
                                   include.lowest = TRUE)

pred_grid$EKE_unbl_hat_binned <- cut(pred_grid$EKE_unbl_hat_gaussian_log,
                                     breaks = breaks_eke_unbl,
                                     include.lowest = TRUE)


pred_grid$EKE_bl_hat_binned %>% levels()


pred_grid$EKE_unbl_hat_binned %>% levels()


saveRDS(pred_grid,"data/pred_grid_EKE_bl_fromQiu2018.Rds")
# Custom function to format legend labels
format_eke_labels <- function(bins) {
  labels <- gsub("[\\(\\)\\[\\]]", "", bins)  # Remove all types of brackets
  labels <- gsub(",", " - ", labels)          # Replace commas with " - "
  return(labels)
}

pred_grid <- pred_grid %>% na.omit() 

# Define global scales for probability contours (same as before)
global_contour_breaks <- c(0.05, 0.10, 0.15, 0.20, 0.25)
global_contour_labels <- as.character(round(global_contour_breaks * 100, 0))


# 1. Aggregate Argo float locations to 5° bins for the season
argo_bins <- df_argo_clean %>%
  mutate(
    lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
    lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
  ) %>%
  group_by(lon_bin_stipple, lat_bin_stipple) %>%
  summarize(count = n(), .groups = "drop")

full_grid <- expand.grid(
  lon_bin_stipple = seq(-180, 180, by = 5),
  lat_bin_stipple = seq(-90, 90, by = 5)
) %>%
  as_tibble()

# Left join the existing argo_bins to the full grid and replace NAs with 0
argo_bins_full <- full_grid %>%
  left_join(argo_bins, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
  mutate(count = ifelse(is.na(count), 0, count))


# 2. Identify undersampled areas (e.g., fewer than 5 profiles)
undersampled <- argo_bins_full %>% filter(count < 1)




undersampled_corners <- undersampled %>%
  rowwise() %>%
  mutate(corners = list(
    data.frame(
      LON = c(lon_bin_stipple,
              lon_bin_stipple + stipple_resolution,
              lon_bin_stipple,
              lon_bin_stipple + stipple_resolution),
      LAT = c(lat_bin_stipple,
              lat_bin_stipple,
              lat_bin_stipple + stipple_resolution,
              lat_bin_stipple + stipple_resolution)
    )
  )) %>%
  ungroup() %>%
  unnest(corners)

# # ----- Prepare the prediction grid for probability of subduction (yearly) -----
# # (Assumes pred_carb_year exists with columns: lon_bin, lat_bin, proportion)
# prediction_grid_year <- pred_carb_year %>%
#   mutate(
#     lon_bin_stipple = floor(lon_bin / stipple_resolution) * stipple_resolution,
#     lat_bin_stipple = floor(lat_bin / stipple_resolution) * stipple_resolution
#   ) %>%
#   left_join(argo_bins_year, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
#   mutate(undersampled = ifelse(is.na(count), FALSE, TRUE))
# 
# # ----- Clean and prepare the EKE prediction grid (yearly) -----
# # (Assumes pred_grid contains the balanced EKE data with columns LON_num, LAT_num, and EKE_bl_hat_binned)
# pred_grid <- pred_grid %>% na.omit()
# 
# # Build undersampled corners
# argo_bins <- df_argo_clean %>%
#   mutate(
#     lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
#     lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
#   ) %>%
#   group_by(lon_bin_stipple, lat_bin_stipple) %>%
#   summarize(count = n(), .groups = "drop")
# 
# # 2. Identify undersampled areas (e.g., fewer than 5 profiles)
# undersampled <- argo_bins %>% filter(count < 5)
# 
# # 3. Generate corner points for each undersampled cell (4 corners per cell)
# undersampled_corners <- undersampled %>%
#   rowwise() %>%
#   mutate(corners = list(
#     data.frame(
#       LON = c(lon_bin_stipple,
#               lon_bin_stipple + stipple_resolution,
#               lon_bin_stipple,
#               lon_bin_stipple + stipple_resolution),
#       LAT = c(lat_bin_stipple,
#               lat_bin_stipple,
#               lat_bin_stipple + stipple_resolution,
#               lat_bin_stipple + stipple_resolution)
#     )
#   )) %>%
#   ungroup() %>%
#   unnest(corners)
# 
# ggplot()+geom_contour(data = pred_full_subd,
#              aes(x = lon_bin, y = lat_bin, z = proportion,
#                  color = factor(round(after_stat(level) * 100, 0))),
#              breaks = global_contour_breaks,
#              alpha = 1,
#              inherit.aes = FALSE,
#              show.legend = TRUE)
# 
# # ----- Build the Plot -----
# Manually defined labels for the balanced EKE legend.
eke_labels_bl <- c(
  "< 1e-4",
  "(1e-4, 1e-3)",
  "(1e-3, 2.5e-3)",
  "(2.5e-3, 5e-3)",
  "(5e-3, 7.5e-3)",
  "(7.5e-3, 1e-2)",
  "(1e-2, 2.5e-2)",
  "(2.5e-2, 5e-2)",
  "(5e-2, 7.5e-2)",
  "(7.5e-2, 1e-1)",
  "(1e-1, 1)"
)

# Manually defined labels for the unbalanced EKE legend.
eke_labels_unbl <- c(
  "< 1e-4",
  "(1e-4, 1e-3)",
  "(1e-3, 2e-3)",
  "(2e-3, 3e-3)",
  "(3e-3, 4e-3)",
  "(4e-3, 5e-3)",
  "(5e-3, 6e-3)",
  "(6e-3, 7e-3)",
  "(7e-3, 8e-3)",
  "(8e-3, 9e-3)",
  "(9e-3, 1e-2)",
  "(1e-2, 2.5e-2)",
  "(2.5e-2, 1e-1)"
)

# Define a common theme without panel borders
common_theme <- theme_bw(base_size = 25) +
  theme(
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 20),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

map_eke_bl <- ggplot() +
  # 2. Rectangles for polar regions (optional)
  geom_rect(aes(xmin = -180, xmax = 180, ymin = 50, ymax = 90),
            fill = "lightblue", inherit.aes = FALSE) +
  geom_rect(aes(xmin = -180, xmax = 180, ymin = -90, ymax = -60),
            fill = "lightblue", inherit.aes = FALSE) +
  
  # 1. Tile plot for yearly EKE background
  geom_tile(data = pred_grid, 
            aes(x = LON_num, y = LAT_num, fill = EKE_bl_hat_binned),
            alpha = 0.95) +
  scale_fill_viridis_d(
    option = "viridis",  
    name = "Balanced EKE\n (m²/s²)",
    guide = guide_legend(reverse = TRUE, na.translate = FALSE),
    labels=eke_labels_bl,
  ) +
  # 3. Overlay probability of subduction contours in red
  new_scale_color() +
  geom_contour(data = pred_full_subd,
               aes(x = lon_bin, y = lat_bin, z = proportion,
                   color = factor(round(after_stat(level) * 100, 0))),
               breaks = global_contour_breaks,
               alpha = 0.75,
               linewidth = 1.5,
               inherit.aes = FALSE,
               show.legend = TRUE) +
  scale_color_brewer(
    palette = "Reds",
    name = "Probability of \nsubduction (%)",
    limits = global_contour_labels,
    breaks = global_contour_labels,
    labels = function(x) paste0(x, "%")
  ) +
  
  # 4. Add stippling for undersampled cells
  geom_point(data = undersampled_corners,
             aes(x = LON, y = LAT),
             color = "white", alpha = 0.6, size = 1.5, shape = 20) +
  
  # 5. Overlay the world map
  geom_sf(data = world,
          fill = "white", color = "white",
          inherit.aes = FALSE) +
  
  # 6. Coordinate system and labels
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  labs(
    title = paste("Balanced Eddy Kinetic Energy & Subduction Probability"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    panel.grid = element_blank()
  )

# Display the plot
print(map_eke_bl)
# Unbalanced EKE Map

map_eke_unbl <- ggplot() +
  # 2. Rectangles for polar regions (optional)
  geom_rect(aes(xmin = -180, xmax = 180, ymin = 50, ymax = 90),
            fill = "lightblue", inherit.aes = FALSE) +
  geom_rect(aes(xmin = -180, xmax = 180, ymin = -90, ymax = -60),
            fill = "lightblue", inherit.aes = FALSE) +
  
  # 1. Tile plot for yearly EKE background
  geom_tile(data = pred_grid, 
            aes(x = LON_num, y = LAT_num, fill = EKE_unbl_hat_binned),
            alpha = 0.75) +
  scale_fill_viridis_d(
    option = "viridis",  
    name = "Unbalanced EKE\n (m²/s²)",
    guide = guide_legend(reverse = TRUE, na.translate = FALSE),
    labels=eke_labels_unbl,
  ) +
  # 3. Overlay probability of subduction contours in red
  new_scale_color() +
  geom_contour(data = pred_full_subd,
               aes(x = lon_bin, y = lat_bin, z = proportion,
                   color = factor(round(after_stat(level) * 100, 0))),
               breaks = global_contour_breaks,
               alpha = 1,
               linewidth = 1.5,
               inherit.aes = FALSE,
               show.legend = TRUE) +
  scale_color_brewer(
    palette = "Reds",
    name = "Probability of \nsubduction (%)",
    limits = global_contour_labels,
    breaks = global_contour_labels,
    labels = function(x) paste0(x, "%")
  ) +
  
  # 4. Add stippling for undersampled cells
  geom_point(data = undersampled_corners,
             aes(x = LON, y = LAT),
             color = "white", alpha = 0.6, size = 1.5, shape = 20) +
  
  # 5. Overlay the world map
  geom_sf(data = world,
          fill = "white", color = "white",
          inherit.aes = FALSE) +
  
  # 6. Coordinate system and labels
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  labs(
    title = paste("Unbalanced Eddy Kinetic Energy & Subduction Probability"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    panel.grid = element_blank()
  )



ggsave(plot=map_eke_bl,filename = "figures/smoothed_map_eke_bl.png",width = 18,height = 10,dpi = 300)
ggsave(plot=map_eke_unbl,filename = "figures/smoothed_map_eke_unbl.png",width = 18,height = 10,dpi = 300)

combined_discrete_eke_map <- (eke_balanced_map + map_subd_full) / (eke_unbalanced_map + map_subd_full)


combined_smoothed_eke_map <- (map_eke_bl + map_subd_full) / (map_eke_unbl + map_subd_full)

ggsave("figures/combined_discrete_eke_map.png", combined_discrete_eke_map, width = 14, height = 10, dpi = 300)

ggsave("figures/combined_smoothed_eke_map.png", combined_smoothed_eke_map, width = 14, height = 10, dpi = 300)



# Confirm it statistically with a GLM 
# pred_full_sub has 1° resolution 
# pred_grid (EKE data) has 0.25° resolution, it needs to be downscaled to 1° resolution 
# so the two datasets can be merged an I can build a GLM of logit(pred_full_subd$proportion) ~ pred_grid$EKE
# Downscale pred_grid from 0.25° resolution to a 1° grid 
# with bin centers computed as floor(LON_num) + 0.5, floor(LAT_num) + 0.5
pred_grid_down <- pred_grid %>%
  mutate(
    lon_bin = floor(LON_num) + 0.5,
    lat_bin = floor(LAT_num) + 0.5
  ) %>%
  group_by(lon_bin, lat_bin) %>%
  summarize(
    EKE_bl_hat_gaussian_log = mean(EKE_bl_hat_gaussian_log, na.rm = TRUE),
    EKE_unbl_hat_gaussian_log = mean(EKE_unbl_hat_gaussian_log, na.rm = TRUE),
    .groups = "drop"
  )

# Check the resulting downscaled grid
print(head(pred_grid_down))

# Merge the downscaled pred_grid with the 1° resolution dataset (pred_full_subd)
merged_data <- left_join(pred_full_subd, pred_grid_down, by = c("lon_bin", "lat_bin"))
merged_data <- left_join(pred_full_subd_highres, pred_grid, by = c("LON_num", "LAT_num"))

pred_full_subd_highres$LON_num <- pred_full_subd_highres$lon_bin 
pred_full_subd_highres$LAT_num <- pred_full_subd_highres$lat_bin 

# Check the merged data
print(head(merged_data))

# Compute the logit of the proportion (ensure that 'proportion' is strictly between 0 and 1)
merged_data <- merged_data %>%
  mutate(logit_proportion = log(proportion / (1 - proportion)))

extra_tropics <- merged_data %>% filter(lat_bin >= 30 | lat_bin <= 30)
tropics <- merged_data %>% filter(abs(lat_bin) < 30)



# Fit a GLM using the downscaled EKE (using balanced EKE as an example)
glm_fit_et <- glm(logit_proportion ~ log(EKE_bl_hat_gaussian_log)+ log(EKE_unbl_hat_gaussian_log),
               data = extra_tropics,
               family = gaussian())
summary(glm_fit_et)


# Fit a GLM using the downscaled EKE (using balanced EKE as an example)
glm_fit_tp <- glm(logit_proportion ~ log(EKE_bl_hat_gaussian_log)+log(EKE_unbl_hat_gaussian_log),
                     data = tropics,
                     family = gaussian())
summary(glm_fit_tp)

# Display the summary of the GLM fit
glm_fit_tp_unbl <- glm(logit_proportion ~ log(EKE_unbl_hat_gaussian_log),
                       data = tropics,
                       family = gaussian())
summary(glm_fit_tp_unbl)


