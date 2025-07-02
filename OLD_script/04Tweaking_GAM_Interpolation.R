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

library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggspatial)
library(viridis)
library(sp)
library(spdep)
library(mgcv)
library(patchwork)

library(fitdistrplus)
library(poweRlaw)

library(viridisLite)
library(scales)



# Resolve function conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")
setwd("/data/GLOBARGO/src/")
df_argo_clean <- read_csv("data/df_argo_loc.csv")
df_complete_clean <- read_csv("data/df_eddy_subduction_anom.csv")
df_carbon_clean <- read_csv("data/df_carbon_subduction_anom.csv")

all_anomalies <- read_csv("/data/GLOBARGO/data/detected_events_abs_sal_var_v5.csv")

# Bin data
# Define bin size for longitude and latitude
bin_size <- 1
longitude_bins <- seq(floor(min(df_argo_clean$LONGITUDE)), ceiling(max(df_argo_clean$LONGITUDE)), by = bin_size)
latitude_bins <- seq(floor(min(df_argo_clean$LATITUDE)), ceiling(max(df_argo_clean$LATITUDE)), by = bin_size)

# Compute bin centers for labeling
lon_centers <- longitude_bins[-length(longitude_bins)] + bin_size / 2
lat_centers <- latitude_bins[-length(latitude_bins)] + bin_size / 2

# Assign bins to Argo profiles and anomalies for each season
assign_bins <- function(df) {
  df %>%
    mutate(
      lon_bin = cut(LONGITUDE, breaks = longitude_bins, include.lowest = TRUE, labels = lon_centers),
      lat_bin = cut(LATITUDE, breaks = latitude_bins, include.lowest = TRUE, labels = lat_centers)
    ) %>%
    mutate(
      lon_bin = as.numeric(as.character(lon_bin)),
      lat_bin = as.numeric(as.character(lat_bin))
    )
}



# Define months for the four distinct periods: DJF, MAM, JJA, SON
djf_months <- c(12, 1, 2)   # December, January, February
mam_months <- c(3, 4, 5)    # March, April, May
jja_months <- c(6, 7, 8)    # June, July, August
son_months <- c(9, 10, 11)  # September, October, November

# Filter data for DJF, MAM, JJA, and SON
df_argo_djf <- df_argo_clean %>% filter(month(TIME) %in% djf_months)
df_argo_mam <- df_argo_clean %>% filter(month(TIME) %in% mam_months)
df_argo_jja <- df_argo_clean %>% filter(month(TIME) %in% jja_months)
df_argo_son <- df_argo_clean %>% filter(month(TIME) %in% son_months)

df_complete_djf <- df_complete_clean %>% filter(month(TIME) %in% djf_months)
df_complete_mam <- df_complete_clean %>% filter(month(TIME) %in% mam_months)
df_complete_jja <- df_complete_clean %>% filter(month(TIME) %in% jja_months)
df_complete_son <- df_complete_clean %>% filter(month(TIME) %in% son_months)

df_carbon_djf <- df_carbon_clean %>% filter(month(TIME) %in% djf_months)
df_carbon_mam <- df_carbon_clean %>% filter(month(TIME) %in% mam_months)
df_carbon_jja <- df_carbon_clean %>% filter(month(TIME) %in% jja_months)
df_carbon_son <- df_carbon_clean %>% filter(month(TIME) %in% son_months)

df_carbon_poc_djf <- df_carbon_with_poc %>% filter(month(TIME) %in% djf_months)
df_carbon_poc_mam <- df_carbon_with_poc %>% filter(month(TIME) %in% mam_months)
df_carbon_poc_jja <- df_carbon_with_poc %>% filter(month(TIME) %in% jja_months)
df_carbon_poc_son <- df_carbon_with_poc %>% filter(month(TIME) %in% son_months)



# Assign bins to Argo profiles and anomalies

df_argo_full <- assign_bins(df_argo_clean)
df_argo_djf <- assign_bins(df_argo_djf)
df_argo_mam <- assign_bins(df_argo_mam)
df_argo_jja <- assign_bins(df_argo_jja)
df_argo_son <- assign_bins(df_argo_son)

df_complete_full <- assign_bins(df_complete_clean)
df_complete_djf <- assign_bins(df_complete_djf)
df_complete_mam <- assign_bins(df_complete_mam)
df_complete_jja <- assign_bins(df_complete_jja)
df_complete_son <- assign_bins(df_complete_son)

df_carbon_full <- assign_bins(df_carbon_clean)
df_carbon_djf <- assign_bins(df_carbon_djf)
df_carbon_mam <- assign_bins(df_carbon_mam)
df_carbon_jja <- assign_bins(df_carbon_jja)
df_carbon_son <- assign_bins(df_carbon_son)



df_with_fp_full <- assign_bins(all_anomalies)


# Compute counts of total profiles and anomalies for each season
compute_counts <- function(df_argo, df_complete) {
  total_counts <- df_argo %>%
    group_by(lon_bin, lat_bin) %>%
    summarize(count_total = n(), .groups = 'drop')
  
  anomaly_counts <- df_complete %>%
    group_by(lon_bin, lat_bin) %>%
    summarize(count_anomaly = n(), .groups = 'drop')
  
  # Merge counts and compute proportions
  merged_counts <- full_join(total_counts, anomaly_counts, by = c("lon_bin", "lat_bin")) %>%
    mutate(
      count_anomaly = ifelse(is.na(count_anomaly), 0, count_anomaly),
      proportion = count_anomaly / count_total
    ) %>%
    filter(count_total > 0) %>%
    filter(proportion < 0.8) %>%  # Filter out high proportions (anomalies > 80%)
    filter(!is.na(lat_bin))
  
  return(merged_counts)
}
# Compute merged_counts of SUBDUCTION for whole year and each season

merged_counts_full <- compute_counts(df_argo_full,df_complete_full)
merged_counts_djf <- compute_counts(df_argo_djf, df_complete_djf)
merged_counts_mam <- compute_counts(df_argo_mam, df_complete_mam)
merged_counts_jja <- compute_counts(df_argo_jja, df_complete_jja)
merged_counts_son <- compute_counts(df_argo_son, df_complete_son)

# Compute merged_counts of CARBON SUBDUCTION for each season
merged_carbon_counts_full <- compute_counts(df_argo_full,df_carbon_full)
merged_carbon_counts_djf <- compute_counts(df_argo_djf, df_carbon_djf)
merged_carbon_counts_mam <- compute_counts(df_argo_mam, df_carbon_mam)
merged_carbon_counts_jja <- compute_counts(df_argo_jja, df_carbon_jja)
merged_carbon_counts_son <- compute_counts(df_argo_son, df_carbon_son)

merged_counts_all_anomalies <- compute_counts(df_argo_full,df_with_fp_full)




# Create proportion maps of SUBDUCTION (w/o carbon) without saving them to files directly

prop_map_full <- ggplot(merged_counts_full, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_full$lon_bin), ylim = range(merged_counts_full$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (whole year)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()

prop_map_all_anom <- ggplot(merged_counts_all_anomalies, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_djf$lon_bin), ylim = range(merged_counts_djf$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (whole year)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()


# Save combined figure
ggsave(filename = "figures/TimeSpaceVar/4SEASONS/combined_proportion_maps.png", plot = combined_proportion_maps, width = 15, height = 12)

###############################################################################
# 1) GAM Fitting Function
###############################################################################
# This function expects a data frame (merged_counts) with:
#  - lon_bin, lat_bin (numeric)
#  - count_total
#  - count_anomaly
# and fits a GAM:   cbind(anomaly, total - anomaly) ~ s(lat, lon).
#
# We specify k=600 or some other suitable value. Feel free to adjust.
# We specify minimum Argo bin in the 5*5 gridcell to learn : 5
fit_gam_season <- function(merged_counts, k_value = 600,argo_min = 0) {
  # Filter out any invalid rows if necessary
  df <- merged_counts %>%
    filter(!is.na(lon_bin), !is.na(lat_bin), count_total > argo_min)
  
  gam_model <- gam(
    cbind(count_anomaly, count_total - count_anomaly) ~
      s(lat_bin, lon_bin, bs = "sos", k = k_value),
    family = binomial(link = "logit"),
    data = df,
    method = "REML"
  )
  return(gam_model)
}

###############################################################################
# 2) Create Prediction Grid Function
###############################################################################
# We'll sample lat_bin, lon_bin in 1° increments (adjust as desired).
# Make sure the range covers the data domain for that season.

create_prediction_grid <- function(merged_counts, step) {
  lon_seq <- seq(min(merged_counts$lon_bin, na.rm = TRUE),
                 max(merged_counts$lon_bin, na.rm = TRUE),by = step)
  
  
  lat_seq <- seq(min(merged_counts$lat_bin, na.rm = TRUE),
                 max(merged_counts$lat_bin, na.rm = TRUE),
                 by = step)
  expand.grid(lon_bin = lon_seq, lat_bin = lat_seq)
}

###############################################################################
# 3) Predict Smoothed Proportions
###############################################################################

predict_gam <- function(gam_model, merged_counts, step) {
  pred_grid <- create_prediction_grid(merged_counts, step)
  pred_grid$proportion <- stats::predict(gam_model, newdata = pred_grid, type = "response")
  pred_grid <- pred_grid %>% filter(!is.na(proportion))
  return(pred_grid)
}

###############################################################################
# 4) Plot a GAM Result (Lat-Lon) With Discrete Color Bins
###############################################################################
# We'll define the scale in a separate function so we can re-use it across plots.

make_discrete_scale <- function(prob_max) {
  # E.g. if prob_max is 0.6 and binwidth is 0.1 => breaks = seq(0, 0.6, 0.1)
  breaks_vec <- c(0,0.1,0.5,0.10,0.20,0.30,0.40,0.50,0.60)
  
  scale_fill_viridis_b(
    name = "Probability",
    breaks = breaks_vec,
    limits = c(0, 0.6),
    oob = scales::squish
  )
}



stipple_resolution <- 5  # degrees

plot_gam_map <- function(pred_grid, world_data, season_label, event_label,argo_months) {
  # 1. Aggregate Argo float locations to 5° bins for the season
  argo_bins <- df_argo_clean %>%
    filter(month(TIME) %in% argo_months) %>%
    mutate(
      lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
      lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
    ) %>%
    group_by(lon_bin_stipple, lat_bin_stipple) %>%
    summarize(count = n(), .groups = "drop")
  
  full_grid <- expand.grid(
    lon_bin_stipple = seq(-180, 175, by = 5),
    lat_bin_stipple = seq(-90, 90, by = 5)
  ) %>%
    as_tibble()
  
  # Left join the existing argo_bins to the full grid and replace NAs with 0
  argo_bins_full <- full_grid %>%
    left_join(argo_bins, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
    mutate(count = ifelse(is.na(count), 0, count))
  
  
  # 2. Identify undersampled areas (e.g., fewer than 5 profiles)
  undersampled <- argo_bins_full %>% filter(count < 1)
  
  # 3. Generate corner points for each undersampled cell (4 corners per cell)
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
  
  
  x_limits <- range(pred_grid$lon_bin, na.rm = TRUE)
  y_limits <- range(pred_grid$lat_bin, na.rm = TRUE)
  
  ggplot() +
    geom_tile(data = pred_grid,
              aes(x = lon_bin, y = lat_bin, fill = proportion)) +
    geom_contour(data = pred_grid,
                 aes(x = lon_bin, y = lat_bin, z = proportion),
                 color = "white", alpha = 0.3,breaks = c(0.001,0.01,0.05,0.10,0.20,0.30,0.40,0.60)) +
    # Add stippling: plot four corner points per undersampled grid cell
    geom_point(data = undersampled_corners,
               aes(x = LON, y = LAT,color="undersampled regions"),
               color = "red", alpha = 1, size = 1, shape = 20) +
    # Add continents
    geom_sf(data = world_data, fill = "lightgrey", color = "lightgrey", inherit.aes = FALSE) +
    coord_sf(xlim = c(-180,180), ylim = c(-90,90), expand = FALSE, crs = st_crs(4326)) +
    labs(
      title = paste0("Estimated ", event_label, " Probability ", season_label, ""),
      x = "Longitude", y = "Latitude"
    ) +
    scale_fill_viridis_b(
      name = "Probability",
      breaks = c(0.001,0.01,0.05,0.10,0.20,0.30,0.40,0.60),
      labels = c("0.1 %","1 %","5 %","10 %","20 %","30 %","40 %","60 %"),
      limits = c(0.001, 0.6),
      oob = scales::squish, na.value = "white") +    # apply the discrete color scale
    theme_minimal(base_size = 25) +  # Increased text size for readability
    theme(legend.title = element_blank(),legend.position = 'bottom', 
          panel.grid = element_blank(),
          legend.text = element_text(size = 20),  # Improve legend labels
          legend.direction = "horizontal", 
          legend.box = "horizontal",
          legend.key.height = unit(1,"lines"),
          legend.key.width = unit(8,"lines"))+ 
    guides(color = guide_legend(title.position = "top", 
                                # hjust = 0.5 centres the title horizontally
                                title.hjust = 0.5,
                                label.position = "bottom")) 
  
  
}

plot_carb_gam_map <- function(pred_grid, world_data, season_label, event_label,argo_months) {
  # 1. Aggregate Argo float locations to 5° bins for the season
  argo_bins <- df_argo_clean %>%
    filter(month(TIME) %in% argo_months) %>%
    mutate(
      lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
      lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
    ) %>%
    group_by(lon_bin_stipple, lat_bin_stipple) %>%
    summarize(count = n(), .groups = "drop")
  
  full_grid <- expand.grid(
    lon_bin_stipple = seq(-180, 175, by = 5),
    lat_bin_stipple = seq(-90, 90, by = 5)
  ) %>%
    as_tibble()
  
  # Left join the existing argo_bins to the full grid and replace NAs with 0
  argo_bins_full <- full_grid %>%
    left_join(argo_bins, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
    mutate(count = ifelse(is.na(count), 0, count))
  
  
  # 2. Identify undersampled areas (e.g., fewer than 5 profiles)
  undersampled <- argo_bins_full %>% filter(count < 1)
  
  # 3. Generate corner points for each undersampled cell (4 corners per cell)
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
  
  
  x_limits <- range(pred_grid$lon_bin, na.rm = TRUE)
  y_limits <- range(pred_grid$lat_bin, na.rm = TRUE)
  
  ggplot() +
    geom_tile(data = pred_grid,
              aes(x = lon_bin, y = lat_bin, fill = proportion)) +
    geom_contour(data = pred_grid,
                 aes(x = lon_bin, y = lat_bin, z = proportion),
                 color = "white", alpha = 0.3,breaks = c(0.0005,0.005,0.025,0.05,0.10,0.15,0.20,0.30) ) +
    # Add stippling: plot four corner points per undersampled grid cell
    geom_point(data = undersampled_corners,
               aes(x = LON, y = LAT,color="undersampled regions"),
               color = "red", alpha = 1, size = 1, shape = 20) +
    # Add continents
    geom_sf(data = world_data, fill = "lightgrey", color = "lightgrey", inherit.aes = FALSE) +
    coord_sf(xlim = c(-180,180), ylim = c(-90,90), expand = FALSE, crs = st_crs(4326)) +
    labs(
      title = paste0("Estimated ", event_label, " Probability ", season_label, ""),
      x = "Longitude", y = "Latitude"
    ) +
    scale_fill_viridis_b(
      name = "Probability",
      breaks = c(0.0005,0.005,0.025,0.05,0.10,0.15,0.20,0.30),
      labels = c("0.05 %","0.5 %","2.5 %","5 %","10 %","15 %","20 %","30 %"),
      limits = c(0.0005, 0.30),
      oob = scales::squish, na.value = "white") +    # apply the discrete color scale
    theme_minimal(base_size = 25) +  # Increased text size for readability
    theme(legend.title = element_blank(),legend.position = 'bottom', 
          panel.grid = element_blank(),
          legend.text = element_text(size = 20),  # Improve legend labels
          legend.direction = "horizontal", 
          legend.box = "horizontal",
          legend.key.height = unit(1,"lines"),
          legend.key.width = unit(8,"lines"))+ 
    guides(color = guide_legend(title.position = "top", 
                                # hjust = 0.5 centres the title horizontally
                                title.hjust = 0.5,
                                label.position = "bottom")) 
  
  
}


###############################################################################
# 5) Build and Plot Subduction GAM for Each Season
###############################################################################

# Fit GAMs

gam_full_subd <- fit_gam_season(merged_counts_full,  k_value = 800,argo_min = 5)

gam_full_all_anom <- fit_gam_season(merged_counts_all_anomalies,  k_value = 800,argo_min = 5)

# High k variant, k = 600 (originally it was 300)

# Predict at submesoscale resolution
pred_full_subd <- predict_gam(gam_full_subd, merged_counts_full, step = 0.25)


saveRDS(pred_full_subd,file = "data/pred_full_subd025.Rds")
# And the false positives as well
pred_full_subd_all_anom_highres <- predict_gam(gam_full_all_anom, merged_counts_full, step = 0.25)


# mask observations with less than 0.001 (0.1) 
pred_full_subd2 <- pred_full_subd
pred_full_subd2$proportion <- ifelse(pred_full_subd$proportion < 0.001,NA, pred_full_subd$proportion)


pred_full_subd_all_anom_highres2 <- pred_full_subd_all_anom_highres
pred_full_subd_all_anom_highres2$proportion <- ifelse(
  pred_full_subd_all_anom_highres2$proportion < 0.001,NA, pred_full_subd_all_anom_highres$proportion)


# Predict by season

pred_djf_subd <- predict_gam(gam_djf_subd, merged_counts_djf, step = 0.25)
pred_mam_subd <- predict_gam(gam_mam_subd, merged_counts_mam, step = 0.25)
pred_jja_subd <- predict_gam(gam_jja_subd, merged_counts_jja, step = 0.25)
pred_son_subd <- predict_gam(gam_son_subd, merged_counts_son, step = 0.25)

# Determine global min/max for subduction proportions
subd_max_full <- max(pred_full_subd2$proportion)
subd_full_scale <- make_discrete_scale(0.6)



# Plot each season with the common discrete scale

map_subd_full <- plot_gam_map(pred_full_subd, world, "(Whole Year)", "Subduction",
                                argo_months = c(1:12))
map_subd_full2 <- plot_gam_map(pred_full_subd2, world, "(Whole Year)", "Subduction",
                              argo_months = c(1:12))


ggsave("figures/map_subd_annual.png",
       map_subd_full, width = 18, height = 10)

ggsave("figures/TimeSpaceVar/4SEASONS/map_subd_annual.png",
       map_subd_full2, width = 18, height = 10)


map_subd_full_all_anom <- plot_gam_map(pred_full_subd_all_anom_highres,
                                       world, "Whole Year (high res + false pos)",
                                       "Subduction", subd_full_scale,argo_months = c(1:12))

map_subd_full_all_anom2 <- plot_gam_map(pred_full_subd_all_anom_highres2,
                                        world, "Whole Year (high res + false pos)",
                                       "Subduction", subd_full_scale,argo_months = c(1:12))
map_subd_full_all_anom2 <- ggarrange(map_subd_full_all_anom2,
                           common.legend = T,legend="bottom")


ggsave("figures/map_subd_full_all_anom2.png",
       map_subd_full_all_anom2, width = 18, height = 10)


map_subd_djf <- plot_gam_map(pred_djf_subd, world, "DJF", "Subduction", subd_scale,c(12,1,2))
map_subd_mam <- plot_gam_map(pred_mam_subd, world, "MAM", "Subduction", subd_scale,c(3:5))
map_subd_jja <- plot_gam_map(pred_jja_subd, world, "JJA", "Subduction", subd_scale,c(6:8))
map_subd_son <- plot_gam_map(pred_son_subd, world, "SON", "Subduction", subd_scale,c(9:11))

map_subd_full <- ggarrange(map_subd_full,
                           common.legend = T,legend="bottom")


combined_subd <- ggarrange(map_subd_djf,map_subd_mam,map_subd_jja,map_subd_son,ncol = 2,nrow=2,
                           common.legend = T,legend="bottom")
ggsave("figures/TimeSpaceVar/4SEASONS/gam_subduction_discrete_yearly.png",
       map_subd_full, width = 18, height = 10)


ggsave("figures/TimeSpaceVar/4SEASONS/gam_subduction_discrete_4seasons.png",
       combined_subd, width = 20, height = 15)

###############################################################################
# 6) Build and Plot Carbon Subduction GAM for Each Season
###############################################################################
gam_full_carb <- fit_gam_season(merged_carbon_counts_full,  k_value = 800,argo_min = 5)

gam_djf_carb <- fit_gam_season(merged_carbon_counts_djf,  k_value = 600,argo_min = 5)
gam_mam_carb <- fit_gam_season(merged_carbon_counts_mam,  k_value = 600,argo_min = 5)
gam_jja_carb <- fit_gam_season(merged_carbon_counts_jja,  k_value = 600,argo_min = 5)
gam_son_carb <- fit_gam_season(merged_carbon_counts_son,  k_value = 600,argo_min = 5)

# Predict
pred_full_carb <- predict_gam(gam_full_carb, merged_carbon_counts_full, step = 0.25)
saveRDS(pred_full_carb,file = "data/pred_full_carb025.Rds")



pred_djf_carb <- predict_gam(gam_djf_carb, merged_carbon_counts_djf, step = 0.25)
pred_mam_carb <- predict_gam(gam_mam_carb, merged_carbon_counts_mam, step = 0.25)
pred_jja_carb <- predict_gam(gam_jja_carb, merged_carbon_counts_jja, step = 0.25)
pred_son_carb <- predict_gam(gam_son_carb, merged_carbon_counts_son, step = 0.25)

pred_djf_carb$season <- "DJF"
pred_mam_carb$season <- "MAM"
pred_jja_carb$season <- "JJA"
pred_son_carb$season <- "SON"

pred_son_carb_season_combined <- rbind(pred_djf_carb,
                                       pred_mam_carb,
                                       pred_jja_carb,
                                       pred_son_carb)

saveRDS(pred_son_carb_season_combined,"data/pred_grid_carb_prob.Rds")

pred_full_carb2 <- pred_full_carb
pred_full_carb2$proportion <- ifelse(pred_full_carb$proportion < 0.0005,NA, pred_full_carb$proportion)


# Global min/max for carbon subduction
carb_values <- c(pred_djf_carb$proportion, pred_mam_carb$proportion,
                 pred_jja_carb$proportion, pred_son_carb$proportion)
carb_min <- 0
carb_max <- max(carb_values, na.rm = TRUE)

carb_max_full <- max(pred_full_carb$proportion)
carb_full_scale <- make_discrete_scale(0,0.25,binwidth = 0.05)

# Discrete color scale for carbon subduction
carb_scale <- make_discrete_scale(carb_min, carb_max, binwidth = 0.05)

map_carb_full <- plot_carb_gam_map(pred_full_carb,  world, "", "Carbon Subduction",argo_months = c(1:12))
map_carb_full <- plot_carb_gam_map(pred_full_carb2,  world, "", "Carbon Subduction",argo_months = c(1:12))

ggsave(plot = map_carb_full,
       filename = "figures/map_carb_full2.png",
       width = 23,height = 8)



map_carb_full_subd_scale <- plot_gam_map(pred_full_carb,  world, "", "Carbon Subduction", subd_scale,argo_months = c(1:12))




library(cowplot)



ggsave("figures/combined_4maps_cowplot.png", final_cow, width=20, height=15, dpi=300)



map_carb_subd_full <- ggarrange(map_subd_full,map_carb_full_subd_scale,common.legend = TRUE)

map_carb_djf <- plot_carb_gam_map(pred_djf_carb, world, "DJF", "Carbon Subduction",c(12,1,2))
map_carb_mam <- plot_carb_gam_map(pred_mam_carb, world, "MAM", "Carbon Subduction",c(3:5))
map_carb_jja <- plot_carb_gam_map(pred_jja_carb, world, "JJA", "Carbon Subduction",c(6:8))
map_carb_son <- plot_carb_gam_map(pred_son_carb, world, "SON", "Carbon Subduction",c(9:11))

combined_carb <- ggarrange(map_carb_djf,map_carb_mam,map_carb_jja,map_carb_son,ncol = 2,nrow=2,
                           common.legend = T,legend="bottom")

ggsave("figures/TimeSpaceVar/4SEASONS/gam_carbon_subduction_discrete_yearly.png",
       map_carb_full, width = 18, height = 10)


ggsave("figures/TimeSpaceVar/4SEASONS/gam_carbon_subduction_discrete_4seasons.png",
       combined_carb, width = 23, height = 15.3)

california <- df_complete_clean %>% 
     filter(LONGITUDE >= -140, LONGITUDE <= -120) %>%
     filter(LATITUDE >= 30, LATITUDE <= 50)


california_tot <- df_argo_clean %>% 
  filter(LONGITUDE >= -140, LONGITUDE <= -120) %>%
  filter(LATITUDE >= 30, LATITUDE <= 50)

california_tot$WMO %>% unique()
gulf_stream <- df_complete_clean %>% 
  filter(LONGITUDE >= -80, LONGITUDE <= -60) %>%
  filter(LATITUDE >= 30, LATITUDE <= 50)


gulf_stream_tot <- df_argo_clean %>% 
  filter(LONGITUDE >= -80, LONGITUDE <= -60) %>%
  filter(LATITUDE >= 30, LATITUDE <= 50)

gulf_stream_tot$month <- month(gulf_stream_tot$TIME)

gulf_stream_tot$month %>% hist()

kuroshio <- df_argo_clean %>% 
  filter(LONGITUDE >= 130, LONGITUDE <= 150) %>%
  filter(LATITUDE >= 25, LATITUDE <= 40)


kuroshio <- df_complete_clean %>% 
  filter(LONGITUDE >= 130, LONGITUDE <= 150) %>%
  filter(LATITUDE >= 25, LATITUDE <= 40)
