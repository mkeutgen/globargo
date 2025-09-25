#############################
# Loading Required Libraries
#############################
library(tidyverse)
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
library(readr)
library(scales)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")
setwd("scripts/")
df_argo_clean <- read_csv("../data/df_argo_loc.csv")
df_argo_clean$TIME %>% min(na.rm = T)
df_argo_clean$TIME %>% max(na.rm = T)

df_complete_clean <- read_csv("../data/manually_verified_physical_subd_events.csv")
df_carbon_clean <- read_csv("../data/df_carbon_subduction_anom.csv")

df_carbon_clean_poc <- read_csv("../data/df_carbon_subduction_anom_with_poc_fromgali.csv")

all_anomalies <- read_csv("../data/detected_physical_subd_events.csv")



df_carbon_with_poc <- read_csv("../data/df_carbon_subduction_anom_with_poc_fromgali.csv")


########################
## Figure 1 : 4 panels 
# A det subd 
# B det carb
# C gam subd
# D gam carb
########################

stipple_resolution <- 5

# GAM 
pred_full_subd <- readRDS("../data/pred_full_subd025.Rds")
pred_full_subd$proportion <- ifelse(pred_full_subd$proportion < 0.001,NA, pred_full_subd$proportion)

pred_full_carb <- readRDS("../data/pred_full_carb025.Rds")
pred_full_carb$proportion <- ifelse(pred_full_carb$proportion < 0.0005,NA, pred_full_carb$proportion)




# ---- 1. Shared objects & helpers -------------------------------------------
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# one shared minimalist theme – keeps type & panel handling identical
theme_map <- function(base_size = 16){
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid      = element_blank(),
      panel.border    = element_rect(colour = "black", fill = NA, size = 0.6),
      axis.ticks      = element_line(colour = "black", size = 0.3),
      axis.text       = element_text(size = base_size * 0.8),
      plot.title      = element_text(hjust = 0.5, face = "bold",
                                     size  = base_size * 1.05,
                                     margin = margin(b = 6)),
      legend.position = "bottom",
      legend.text     = element_text(size = base_size * 0.8),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width  = unit(2.2 ,  "cm"),
      plot.margin       = margin(2, 2, 2, 2)
    )
}


label_lon <- function(x){
  paste0(abs(x), "°",
         ifelse(x < 0, "W",
                ifelse(x > 0, "E", "")))
}
label_lat <- function(x){
  paste0(abs(x), "°",
         ifelse(x < 0, "S",
                ifelse(x > 0, "N", "")))
}





plot_gam_map <- function(pred_grid, world_data, title_label,argo_months) {
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
      title = paste0(title_label),
      x = "Longitude", y = "Latitude"
    ) +
    scale_fill_viridis_b(
      name = "",
      breaks = c(0.001,0.01,0.05,0.10,0.20,0.30,0.40,0.60),
      labels = c("0.1","1","5","10","20","30","40","60 (%)"),
      limits = c(0.001, 0.6),
      oob = scales::squish, na.value = "white") +    # apply the discrete color scale
    theme_minimal(base_size = 25) +  # Increased text size for readability
    theme(legend.title = element_blank(),legend.position = 'bottom', 
          panel.grid = element_blank(),
          legend.text = element_text(size = 20),  # Improve legend labels
          legend.direction = "horizontal", 
          legend.box = "horizontal",
          legend.key.height = unit(1,"lines"),
          legend.key.width = unit(12,"lines"))+ 
    guides(color = guide_legend(title.position = "top", 
                                # hjust = 0.5 centres the title horizontally
                                title.hjust = 0.5,
                                label.position = "bottom")) 
  
  
}

plot_carb_gam_map <- function(pred_grid, world_data, title_label,argo_months) {
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
  
  
  # 2. Identify undersampled areas (e.g., fewer than 1 profiles)
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
      title = paste0(title_label),
      x = "Longitude", y = NULL
    ) +
    scale_fill_viridis_b(
      name = "",
      breaks = c(0.0005,0.005,0.025,0.05,0.10,0.15,0.20,0.30),
      labels = c("0.05","0.5","2.5","5","10","15","20","30 (%)"),
      limits = c(0.0005, 0.30),
      oob = scales::squish, na.value = "white") +    # apply the discrete color scale
    theme_minimal(base_size = 25) +  # Increased text size for readability
    theme(legend.title = element_blank(),legend.position = 'bottom', 
          panel.grid = element_blank(),
          legend.text = element_text(size = 20),  # Improve legend labels
          legend.direction = "horizontal", 
          legend.box = "horizontal",
          legend.key.height = unit(1,"lines"),
          legend.key.width = unit(12,"lines"))+ 
    guides(color = guide_legend(title.position = "top", 
                                # hjust = 0.5 centres the title horizontally
                                title.hjust = 0.5,
                                label.position = "bottom")) 
  
  
}


# For discrete‑vs‑continuous palettes
pal_events   <- c("Argo floats"        = "grey60",
                  "Physical Subduction events"  = "#E3E418FF",
                  "Carbon subduction events" = "#47C16EFF")

# helper to reuse long/lat ticks
ticks_x <- seq(-180, 180, by = 60)
ticks_y <- seq( -90,  90, by = 30)

# ---- 2. Point‑distribution maps --------------------------------------------
plot_subduct_distrib <- ggplot() +
  geom_sf(data = world, fill = "black", colour = "black", linewidth = 0.2) +
  geom_point(data = df_argo_clean,
             aes(LONGITUDE, LATITUDE, colour = "Argo floats"),
             alpha = .5, size = 1) +
  geom_point(data = df_complete_clean,
             aes(LONGITUDE, LATITUDE, colour = "Physical Subduction events"),
             alpha = .5, size = 1) +
  scale_colour_manual("", values = pal_events) +
  coord_sf(expand = FALSE, xlim = c(-180,180), ylim = c(-90,90)) +
  scale_x_continuous(breaks = ticks_x, labels = label_lon) +
  scale_y_continuous(breaks = ticks_y, labels = label_lat) +  labs(title = "a • Physical Subduction (4,377 Profiles)",
                                                                   x = NULL, y = "Latitude") +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_map()

plot_carbon_distrib <- ggplot() +
  geom_sf(data = world, fill = "black", colour = "black", linewidth = 0.2) +
  geom_point(data = df_argo_clean,
             aes(LONGITUDE, LATITUDE, colour = "Argo floats"),
             alpha = .5, size = 1) +
  geom_point(data = df_carbon_clean,
             aes(LONGITUDE, LATITUDE, colour = "Carbon subduction events"),
             alpha = .5, size = 1) +
  scale_colour_manual("", values = pal_events) +
  coord_sf(expand = FALSE, xlim = c(-180,180), ylim = c(-90,90)) +
  scale_x_continuous(breaks = ticks_x, labels = label_lon) +
  scale_y_continuous(breaks = ticks_y, labels = label_lat) +
  labs(title = "b • Carbon Subduction (1,333 Profiles)",
       x = NULL, y = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme_map()

# ---- 3. GAM maps ------------------------------------------------------------
# NB: your two “plot_*_gam_map()” functions already return complete ggplot
#     objects, so here we only apply the same theme & drop redundant legend.

map_subd_full <- plot_gam_map( pred_full_subd, world,"c • Physical Subduction Probability", 1:12) +
  theme_map() + guides(colour = "none")+theme(legend.spacing.x = unit(0.8, "cm"))  +
  scale_y_continuous(breaks = ticks_y, labels = label_lat) 

map_carb_full <- plot_carb_gam_map(pred_full_carb, world, "d • Carbon Subduction Probability", 1:12) +
  theme_map() + guides(colour = "none")+theme(legend.spacing.x = unit(0.8, "cm")) +
  scale_y_continuous(breaks = ticks_y, labels = label_lat) 

## (a) point‑cloud row — collect its legend only
top_row <- (plot_subduct_distrib | plot_carbon_distrib) +
  theme(legend.position = "bottom")           # put the single legend under that row

## (b) GAM row — DO *NOT* collect guides, leave each plot’s legend intact
bottom_row <- map_subd_full | map_carb_full            # two distinct legends

## (c) stack the rows
fig_final <- top_row / bottom_row +
  plot_layout(heights = c(1, 1),axis_titles = "collect")            # equal height rows

ggsave(fig_final,filename = "../pubfig/figure1_subd_gam.png",width = 16.18,height = 10)


########################
## Figure 2 : seasonality,  
# see Python notebook.
########################

##################################
## Figure 3 : Depths of subduction
##################################

df_complete_clean <- read_csv("../data/df_eddy_subduction_anom.csv")
add_month_region <- function(df) {
  df %>%
    mutate(
      month = month(TIME, label = TRUE, abbr = TRUE),  # Extract month as a factor
      region = case_when(
        LATITUDE > 30 & LONGITUDE >= -100 & LONGITUDE <= 20 ~ "North Atlantic",  # Above 30°N
        LATITUDE > 30 & (LONGITUDE < -100 | LONGITUDE > 120) ~ "North Pacific",  # Above 30°N
        LATITUDE >= 0 & LATITUDE <= 30 ~ "Northern Tropics",                     # 0° to 30°N
        LATITUDE < 0 & LATITUDE >= -30 ~ "Southern Tropics",                     # 0° to 30°S
        LATITUDE < -30 ~ "Southern Ocean",                                       # Below 30°S
        TRUE ~ NA_character_  # Exclude undefined regions
      )
    )
}
df_argo_clean <- add_month_region(df_argo_clean)
df_complete_clean <- add_month_region(df_complete_clean)
df_carbon_clean <- add_month_region(df_carbon_clean)


pres_data <- df_complete_clean %>% filter(!is.na(region), !is.na(month))
pres_data$PRES_ADJUSTED %>% median()
#pres_data$PRES_ADJUSTED <- pres_data$PRES_ADJUSTED-200


carbon_pres_data <- df_carbon_clean %>% filter(!is.na(region), !is.na(month))
#carbon_pres_data$PRES_ADJUSTED <- carbon_pres_data$PRES_ADJUSTED-200
carbon_pres_data$PRES_ADJUSTED %>% median()

carbon_pres_data$PRES_ADJUSTED %>% quantile()
pres_data$PRES_ADJUSTED %>% quantile(0.8)

pres_data$Type <- "Subduction"
carbon_pres_data$Type <- "Carbon Subduction"
# Combine data and create a new variable to distinguish datasets
combined_data <- bind_rows(
  pres_data %>% mutate(Type = "Physical Subduction Events"),
  carbon_pres_data %>% mutate(Type = "Carbon Subduction Events")
)

p_gt500 <- combined_data %>%                                        # your tibble
  group_by(region, Type) %>%                             # Region × event class
  summarise(
    n_total  = n(),                                      # denominator
    n_ge500  = sum(PRES_ADJUSTED >= 500),                # numerator
    p_ge500  = mean(PRES_ADJUSTED >= 500),               # ≡ n_ge500 / n_total
    .groups = "drop"
  ) %>% 
  arrange(region, desc(p_ge500))                         # optional

combined_data_copy <- combined_data
combined_data_copy$region <- "Global"

combined_data_global <- rbind(combined_data_copy,combined_data) 
combined_data_global$region <- factor(combined_data_global$region,levels=c("Global","North Atlantic","North Pacific","Northern Tropics",
                                                                           "Southern Tropics","Southern Ocean"))
region_levels <- levels(combined_data_global$region)
region_labs <- setNames(
  paste0(letters[seq_along(region_levels)], " \u2022 ", region_levels),
  region_levels
)

subd_distrib_plot <- ggplot(combined_data_global,
                            aes(x = PRES_ADJUSTED,
                                fill = Type, colour = Type)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(
    . ~ region,
    ncol   = 2,
    labeller = labeller(region = region_labs)   # ← HERE
  ) +
  theme_minimal() +
  theme(
    strip.text      = element_text(size = 15, face = "bold"),
    axis.text.x     = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y     = element_text(size = 15),
    axis.title.y    = element_text(size = 16),
    title           = element_text(size = 15),
    legend.text  =  element_text(size = 15),
    legend.position = "bottom"
  ) +
  labs(
    x     = "Depth (m)",
<<<<<<< HEAD
    y     = "Density"  ) +
=======
    y     = "Density",
    title = ""
  ) +
>>>>>>> 9e67aa73f649f0aca85c1e01172ba8d6c7997b52
  scale_x_continuous(
    breaks  = seq(0, 1000, by = 200),
    labels  = seq(0, 1000, by = 200)
  ) +
  guides(fill = "none")

# PDFs in mean seq time : 
mean_seq_df <- read_csv("../data/df_meanseq.csv")
mean_seq_df <- mean_seq_df %>% select(!depth) %>%
  na.omit() %>% filter(mean_seq_time > 0)

bin_size <- 1


# Create 3D coordinate matrices
coords_combined <- combined_data %>% 
  select(LATITUDE, LONGITUDE, PRES_ADJUSTED) %>% 
  as.matrix()

coords_meanseq <- mean_seq_df %>% 
  select(latitude, longitude, DEPTH) %>% 
  as.matrix()

# Find nearest neighbors
library(RANN)
nn <- nn2(coords_meanseq, coords_combined, k = 1)

# 4. Add the matched fseq values
combined_data$mean_seq <- mean_seq_df$mean_seq_time[nn$nn.idx]

subd_distrib_plot_time_seq <- ggplot(combined_data, aes(x = mean_seq, fill = Type,color = Type)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(. ~ region, scales = "free_y", ncol = 2,axes="all_x") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 16),
    legend.title = element_blank(),
    legend.text =  element_blank()
  ) +
  labs(
    x = "Depth (m)",
    y = "Density"  )+  scale_x_continuous(
    breaks = seq(0, 800, by = 200),
    labels = seq(200, 1000, by = 200)
  )+ 
  theme_minimal() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15,angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 16),
    title = element_text(size=15),
    legend.position = "bottom"
  )+  
  guides(fill = "none")

combined_data %>% filter(combined_data$Type == "Carbon Subduction Events") %>% pull(mean_seq) %>% summary()
combined_data %>% filter(combined_data$Type == "Carbon Subduction Events") %>% pull(mean_seq) %>% quantile(x = 0.80)

combined_data %>% filter(combined_data$Type == "Physical Subduction Events") %>% pull(mean_seq) %>% summary()

subd_distrib_plot_time_seq_combined <- ggplot(combined_data, aes(x = mean_seq, fill = Type,color = Type)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 16),
    legend.title = element_blank(),
    legend.text =  element_blank()
  ) +
  labs(
    x = "Mean Sequestration time (years)",
    y = "Density"  )+  scale_x_continuous(
    breaks = seq(0, 2000, by = 200),
    labels = seq(0, 2000, by = 200)
  )+ 
  theme_minimal() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15,angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 16),
    title = element_text(size=15),
    legend.position = "bottom"
  )+  
  guides(fill = "none")

# ------------------------------------------------------------------
# 1.  Single base‑theme so the two panels match perfectly
# ------------------------------------------------------------------
base_cdf_theme <- theme_minimal() %+replace%
  theme(
    strip.text      = element_text(size = 15, face = "bold"),
    axis.text.x     = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y     = element_text(size = 15,angle = 45),
    axis.title      = element_text(size = 16,margin = margin(t=40)),
    plot.title      = element_text(size = 15),
    legend.position = "bottom",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 15),
  )

# ------------------------------------------------------------------
# 2.  ECDF by sequestration time
# ------------------------------------------------------------------
subd_cdf_time <- ggplot(combined_data,
                        aes(mean_seq, colour = Type)) +
  stat_ecdf(size = 1.3, pad = FALSE) +
  scale_color_viridis_d() +
  labs(x = "Mean sequestration time (years)",
       y = "Cumulative probability",
       title = "a • Sequestration‑time of events CDF",
       color = NULL) +
  scale_x_continuous(breaks = seq(0, 1200, 200),
                     limits = c(0, 1200),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(breaks = seq(0, 1, 0.2),
                     labels = scales::percent_format(accuracy = 1))+ 
  base_cdf_theme+theme_map(base_size = 16)

# ------------------------------------------------------------------
# 3.  ECDF by depth
# ------------------------------------------------------------------
subd_cdf_depth <- ggplot(combined_data,
                         aes(PRES_ADJUSTED, colour = Type)) +
  stat_ecdf(size = 1.3, pad = FALSE) +
  scale_color_viridis_d() +
  labs(x = "Depth of events (m)",
       y = "Cumulative probability",
       title = "b • Depth of events CDF",color = NULL) +
  scale_x_continuous(breaks = seq(200, 1000, 200),
                     limits = c(200, 1000),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(breaks = seq(0, 1, 0.2),
                     labels = scales::percent_format(accuracy = 1)) +theme_map()+
  base_cdf_theme+theme_map(base_size = 16)

# ------------------------------------------------------------------
# 4.  Combine — one legend, two panels
# ------------------------------------------------------------------
cdf_combined <- (subd_cdf_time + subd_cdf_depth) +
  plot_layout(guides = "collect") &      # merge identical legend
  theme(legend.position = "bottom")      # keep legend centred




# ------------------------------------------------------------------
# 5.  Export
# ------------------------------------------------------------------
ggsave("../pubfig/fig3_cdf_time_vs_depth.png",
       cdf_combined,
       width  = 13,   # adjust to journal’s column width
       height = 10,
       dpi    = 300)


subd_pdf_time <- ggplot(combined_data,
                        aes(mean_seq, colour = Type, fill = Type)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x     = "Mean sequestration time (years)",
       y     = "Probability density",
       title = "a • Sequestration‑time of events PDF",
       color = NULL, fill = NULL) +
  scale_x_continuous(breaks = seq(0, 1200, 200),
                     limits = c(0, 1200),
                     expand = expansion(mult = c(0, 0.02))) +
  base_cdf_theme+theme_map(base_size = 16)

# ──────────────────────────────────────────────────────────────────
# 2.  PDF by depth
# ──────────────────────────────────────────────────────────────────
subd_pdf_depth <- ggplot(combined_data,
                         aes(PRES_ADJUSTED, colour = Type, fill = Type)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x     = "Depth of events (m)",
       y     = "Probability density",
       title = "b • Depth of events PDF",
       color = NULL, fill = NULL) +
  scale_x_continuous(breaks = seq(200, 1000, 200),
                     limits = c(200, 1000),
                     expand = expansion(mult = c(0, 0.02))) +
  base_cdf_theme+theme_map(base_size = 16)

# ──────────────────────────────────────────────────────────────────
# 3.  Combine ― share a single legend
# ──────────────────────────────────────────────────────────────────
pdf_combined <- subd_distrib_plot +theme(legend.position = "bottom")

# ──────────────────────────────────────────────────────────────────
# 4.  Save the figure
# ──────────────────────────────────────────────────────────────────
ggsave("../pubfig/fig3_pdf_time.png",
       subd_distrib_plot,
       width  = 13,   # same dimensions as CDF version
       height = 10,
       dpi    = 300)


# Figure 3 for SI 


subd_pdf_time_region <- subd_pdf_time + facet_wrap(. ~ region, axes = "all")+theme_minimal()+scale_x_continuous(breaks = seq(0, 1200, 400))
subd_pdf_depth_region <- subd_pdf_depth + facet_wrap(. ~ region, axes = "all")+theme_minimal()+scale_x_continuous(breaks = seq(0, 1000, 400))

pdf_combined_region <- (subd_pdf_time_region + subd_pdf_depth_region) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("../pubfig/fig3_SI_pdf_time_vs_depth_region.png",
       pdf_combined_region,
       width  = 13,   # same dimensions as CDF version
       height = 5.5,
       dpi    = 300)


subd_cdf_time_region <- subd_cdf_time + facet_wrap(. ~ region, axes = "all")+theme_minimal()+scale_x_continuous(breaks = seq(0, 1200, 400))
subd_cdf_depth_region <- subd_cdf_depth + facet_wrap(. ~ region, axes = "all")+theme_minimal()+scale_x_continuous(breaks = seq(0, 1000, 400))

cdf_combined_region <- (subd_cdf_time_region + subd_cdf_depth_region) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave("../pubfig/fig3_SI_cdf_time_vs_depth_region.png",
       cdf_combined_region,
       width  = 13,   # same dimensions as CDF version
       height = 5.5,
       dpi    = 300)

###############################################
## Figure 4 : EKE map combined with SI vs MLI 
###############################################
# data too big for github, can be downloaded on zenodo.
pred_grid <- readRDS("/data/GLOBARGO/src/data/pred_grid_EKE_bl_fromQiu2018.Rds")

format_eke_labels <- function(bins) {
  labels <- gsub("[\\(\\)\\[\\]]", "", bins)  # Remove all types of brackets
  labels <- gsub(",", " - ", labels)          # Replace commas with " - "
  return(labels)
}

breaks_eke_bl <- c(0, 1e-4, 1e-3, 2.5e-3, 5e-3, 7.5e-3, 1e-2, 2.5e-2, 5e-2, 7.5e-2, 1e-1, 1)

pred_grid$EKE_bl_hat_binned <- cut(pred_grid$EKE_bl_hat_gaussian_log,
                                   breaks = breaks_eke_bl,
                                   include.lowest = TRUE)

pred_grid$EKE_unbl_hat_binned <- cut(pred_grid$EKE_unbl_hat_gaussian_log,
                                     breaks = breaks_eke_bl,
                                     include.lowest = TRUE)



pred_grid <- pred_grid %>% na.omit() 

bin_levels <- cut(breaks_eke_bl[-1],               # cut the interval end-points
                  breaks_eke_bl,
                  include.lowest = TRUE) %>% 
  levels()

pred_grid <- pred_grid %>% 
  mutate(
    EKE_bl_hat_binned  = factor(EKE_bl_hat_binned,   levels = bin_levels),
    EKE_unbl_hat_binned= factor(EKE_unbl_hat_binned, levels = bin_levels)
  )

scale_EKE <- scale_fill_viridis_d(
  option = "viridis",
  name   = "EKE (m² s⁻²)",
  limits = bin_levels,          # guarantees the same palette order
  labels = eke_labels_bl,       # the list you already prepared
  guide  = guide_legend(reverse = TRUE, na.translate = FALSE)
)
# Define global scales for probability contours (same as before)
global_contour_breaks <- c(0.05,0.10,0.20,0.30,0.40,0.60)
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



eke_labels_bl <- list(
  expression("<"~10^-4),                                      #  < 10⁻⁴
  expression(10^-3),                           # 10⁻⁴ – 10⁻³
  expression(2.5%*%10^-3),                     # 10⁻³ – 2.5·10⁻³
  expression(5%*%10^-3),                 # 2.5·10⁻³ – 5·10⁻³
  expression(7.5%*%10^-3),                 # 5·10⁻³ – 7.5·10⁻³
  expression(10^-2),                     # 7.5·10⁻³ – 10⁻²
  expression(2.5%*%10^-2),                     # 10⁻² – 2.5·10⁻²
  expression(5%*%10^-2),                 # 2.5·10⁻² – 5·10⁻²
  expression(7.5%*%10^-2),                 # 5·10⁻² – 7.5·10⁻²
  expression(10^-1),                     # 7.5·10⁻² – 10⁻¹
  expression(">"~1)                                  # ≥ 10⁻¹
)

pred_grid$EKE_unbl_hat_binned %>% unique()


# eke_labels_unbl <- list(
#   expression("<"~10^-4),                                      #  < 10⁻⁴
#   expression(10^-3),                           # 10⁻⁴ – 10⁻³
#   expression(2%*%10^-3),                     # 10⁻³ – 2.5·10⁻³
#   expression(3%*%10^-3),                 # 2.5·10⁻³ – 5·10⁻³
#   expression(4%*%10^-3),                 # 5·10⁻³ – 7.5·10⁻³
#   expression(5%*%10^-3), 
#   expression(6%*%10^-3),                     # 10⁻³ – 2.5·10⁻³
#   expression(7%*%10^-3),                 # 2.5·10⁻³ – 5·10⁻³
#   expression(8%*%10^-3),                 # 5·10⁻³ – 7.5·10⁻³
#   expression(9%*%10^-3),                     # 7.5·10⁻³ – 10⁻²
#   expression(10^-2),                     # 7.5·10⁻³ – 10⁻²
#   expression(2.5%*%10^-2),                     # 7.5·10⁻³ – 10⁻²
#   expression(10^-1)
# )



map_eke_bl <- ggplot() +
  ## (optional) light‑blue polar rectangles
  geom_rect(aes(xmin = -180, xmax = 180, ymin =  50, ymax =  90),
            fill = "lightblue", colour = NA) +
  geom_rect(aes(xmin = -180, xmax = 180, ymin = -90, ymax = -60),
            fill = "lightblue", colour = NA) +
  ## 1‑B. Subduction probability contours
  new_scale_colour() +
  geom_contour(data  = pred_full_subd,
               aes(lon_bin, lat_bin, z = proportion,
                   colour = factor(round(after_stat(level) * 100))),
               breaks     = global_contour_breaks,
               linewidth  = 1.4, alpha = 0.75, show.legend = TRUE) +
  scale_colour_brewer(
    palette = "Reds",
    name    = "Physical Subduction\nprobability (%)",
    limits  = global_contour_labels,
    breaks  = global_contour_labels,
    labels  = paste0(global_contour_labels, "%"),
    guide   = guide_legend(
      override.aes = list(  # what the legend key actually shows
        fill      = NA,      # <- no rectangle
        linewidth = 2        # thicker line so it remains visible
      )
    )
  )+
  ## 1‑A. Balanced EKE background
  geom_tile(data = pred_grid,
            aes(LON_num, LAT_num, fill = EKE_bl_hat_binned),
            alpha = 0.95) +
  scale_fill_viridis_d(
    option = "viridis",
    name   = "Balanced EKE\n(m² s⁻²)",
    guide  = guide_legend(reverse = TRUE, na.translate = FALSE),
    labels = eke_labels_bl
  ) +
  geom_contour(data  = pred_full_subd,
               aes(lon_bin, lat_bin, z = proportion,
                   colour = factor(round(after_stat(level) * 100))),
               breaks     = global_contour_breaks,
               linewidth  = 1.4, alpha = 0.75, show.legend = FALSE)+
  
  ## 1‑C. Stippling undersampled 5°×5° cells
  geom_point(data = undersampled_corners,
             aes(LON, LAT), colour = "white",
             size = 1.4, alpha = 0.6, shape = 20) +
  
  ## 1‑D. Coastlines
  geom_sf(data = world, fill = "white", colour = "white", linewidth = 0.3) +
  
  ## 2.  Common axis / coordinate settings
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  scale_x_continuous(breaks = ticks_x, labels = label_lon) +
  scale_y_continuous(breaks = ticks_y, labels = label_lat) +
  
  ## 3.  Labels
  labs(title = "a • Balanced EKE & Subduction Probability",
       x = "Longitude", y = "Latitude") +
  
  ## 4.  Unified theme + legend style
  theme_map(base_size = 18) +
  theme(
    legend.position   = "right",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.key.width  = unit(2.2, "cm"),  # same as other plots
    legend.key.height = unit(0.5, "cm")
  ) +
  ## 4‑B.  Make the contour legend points larger
  guides(colour = guide_legend(override.aes = list(size = 4)))

# map unbalanced
map_eke_unbl <- ggplot() +
  ## (optional) light‑blue polar rectangles
  geom_rect(aes(xmin = -180, xmax = 180, ymin =  50, ymax =  90),
            fill = "lightblue", colour = NA) +
  geom_rect(aes(xmin = -180, xmax = 180, ymin = -90, ymax = -60),
            fill = "lightblue", colour = NA) +
  ## 1‑B. Subduction probability contours
  new_scale_colour() +
  geom_contour(data  = pred_full_subd,
               aes(lon_bin, lat_bin, z = proportion,
                   colour = factor(round(after_stat(level) * 100))),
               breaks     = global_contour_breaks,
               linewidth  = 1.4, alpha = 0.75, show.legend = TRUE) +
  scale_colour_brewer(
    palette = "Reds",
    name    = "Physical Subduction\nprobability (%)",
    limits  = global_contour_labels,
    breaks  = global_contour_labels,
    labels  = paste0(global_contour_labels, "%"),
    guide   = guide_legend(
      override.aes = list(  # what the legend key actually shows
        fill      = NA,      # <- no rectangle
        linewidth = 2        # thicker line so it remains visible
      )
    )
  )+
  ## 1‑A. Balanced EKE background
  geom_tile(data = pred_grid,
            aes(LON_num, LAT_num, fill = EKE_unbl_hat_binned),
            alpha = 0.95) +
  scale_fill_viridis_d(
    option = "viridis",
    name   = "Unbalanced EKE\n(m² s⁻²)",
    guide  = guide_legend(reverse = TRUE, na.translate = FALSE),
    labels = eke_labels_bl
  ) +
  geom_contour(data  = pred_full_subd,
               aes(lon_bin, lat_bin, z = proportion,
                   colour = factor(round(after_stat(level) * 100))),
               breaks     = global_contour_breaks,
               linewidth  = 1.4, alpha = 0.75, show.legend = FALSE)+
  
  ## 1‑C. Stippling undersampled 5°×5° cells
  geom_point(data = undersampled_corners,
             aes(LON, LAT), colour = "white",
             size = 1.4, alpha = 0.6, shape = 20) +
  
  ## 1‑D. Coastlines
  geom_sf(data = world, fill = "white", colour = "white", linewidth = 0.3) +
  
  ## 2.  Common axis / coordinate settings
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  scale_x_continuous(breaks = ticks_x, labels = label_lon) +
  scale_y_continuous(breaks = ticks_y, labels = label_lat) +
  
  ## 3.  Labels
  labs(title = "b • Unbalanced EKE & Subduction Probability",
       x = "Longitude", y = "Latitude") +
  
  ## 4.  Unified theme + legend style
  theme_map(base_size = 18) +
  theme(
    legend.position   = "right",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.key.width  = unit(2.2, "cm"),  # same as other plots
    legend.key.height = unit(0.5, "cm")
  ) +
  ## 4‑B.  Make the contour legend points larger
  guides(colour = guide_legend(override.aes = list(size = 4)))

















# Read data 
df_full <- read_csv("data/dataframe_full_mld_and_n2.csv") 

# Standardize the log-transformed predictors
df_full <- df_full %>%
  mutate(
    log_cleaned_mld_std = scale(log(cleaned_mld)),
    log_cleaned_N2_std = scale(log(cleaned_N2))
  )

df_full$Anomaly <- df_full$Anomaly %>% as_factor()
df_full$Anomaly_text <- ifelse(df_full$Anomaly == 1,"Subduction","No Subduction")

df_full %>% filter(Anomaly==0) %>% pull(cleaned_mld) %>% median()
df_full %>% filter(Anomaly==1) %>% pull(cleaned_mld) %>% median()

df_full %>% filter(Anomaly==0) %>% pull(cleaned_mld) %>% IQR()
df_full %>% filter(Anomaly==1) %>% pull(cleaned_mld) %>% IQR()


df_full %>% filter(Anomaly==0) %>% pull(cleaned_N2) %>% median(na.rm = T)
df_full %>% filter(Anomaly==1) %>% pull(cleaned_N2) %>% median(na.rm = T)

df_full %>% filter(Anomaly==0) %>% pull(cleaned_N2) %>% IQR(na.rm = T)
df_full %>% filter(Anomaly==1) %>% pull(cleaned_N2) %>% IQR(na.rm = T)


df_summary <- df_full %>%
  group_by(Anomaly) %>%
  summarise(
    median_cleaned_mld = median(cleaned_mld, na.rm = TRUE),
    mean_cleaned_mld   = mean(cleaned_mld, na.rm = TRUE),
    median_Max_N2      = median(Max_N2, na.rm = TRUE),
    mean_Max_N2        = mean(Max_N2, na.rm = TRUE),
    .groups = "drop"
  )

df_summary


# Density plot for Mixed Layer Depth (MLD)
density_mld <- ggplot(df_full, aes(x = log_cleaned_mld, fill = factor(Anomaly_text))) +
  geom_density(alpha = 0.7, color = "black", size = 1) +
  labs(
    title = "c • Mixed-Layer Depth PDF",
    x = "Log(Mixed-Layer Depth)",
    y = "Density",
    fill = ""
  ) +
  theme_bw(base_size = 25) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(face = "bold", size = 25, hjust = 0.5),
    plot.subtitle = element_text(size = 15, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 25),
    legend.text = element_text(size = 25)
  )+ scale_fill_viridis_d(begin=0.5)

# Density plot for Brunt-Väisälä Frequency (N²)
density_N2 <- ggplot(df_full, aes(x = log(cleaned_N2), fill = factor(Anomaly_text))) +
  geom_density(alpha = 0.7, color = "black", size = 1) +
  labs(
    title = "d • Brunt-Väisälä PDF",
    x = "Log(Brunt-Väisälä Frequency)",
    y = "Density",
    fill = ""
  ) +
  theme_bw(base_size = 25) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(face = "bold", size = 25, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 25),
    legend.text = element_text(size = 25),
  )+ scale_fill_viridis_d(begin=0.5)


density_mld  <- density_mld  + theme_map()
density_N2   <- density_N2   + theme_map()

# ──────────────────────────────────────────────────────────────────
# 2.  Compose: map on top, two PDFs below
#     (guides = "collect" merges identical legends into one row)
# ──────────────────────────────────────────────────────────────────
map_eke_bl_adj <- map_eke_bl + scale_EKE + theme_map() + theme(legend.position = "right")
map_eke_unbl_adj <- map_eke_unbl + scale_EKE + theme_map() + theme(legend.position = "none")

top_row <- (map_eke_bl_adj | map_eke_unbl_adj)+
  plot_layout(guides = "collect")

bottom_row <- (density_mld | density_N2) +
  plot_layout(guides = "collect")
combined_fig <- top_row /
  bottom_row +
  plot_layout(heights = c(1.2, 0.8)) 

# ──────────────────────────────────────────────────────────────────
# 3.  Save
# ──────────────────────────────────────────────────────────────────
ggsave("../pubfig/figure4_EKE_map_and_PDFs.png",
       combined_fig,
       width  = 18,   # adjust to journal’s column width
       height = 10,
       dpi    = 400)

# Figure 4 BIS, SI, KS distances,
# lines below are commented because the bootstrapping slows down code.

# 
# ks.test(df_full$cleaned_mld[df_full$Anomaly ==1],
#         df_full$cleaned_mld[df_full$Anomaly ==0],alternative="less")
# 
# ks.test(df_full$Max_N2[df_full$Anomaly ==1],
#         df_full$Max_N2[df_full$Anomaly ==0],alternative="greater")
# 
# set.seed(124)       # for reproducibility
# B <- 10000          # number of bootstrap iterations
# n <- nrow(df_full) # total number of observations
# 
# # Storage vectors for the KS statistics in each bootstrap iteration
# ks_mld <- numeric(B)
# ks_n2  <- numeric(B)
# 
# # Original (observed) KS stats for reference:
# ks_mld_obs <- ks.test(
#   df_full$cleaned_mld[df_full$Anomaly == 1],
#   df_full$cleaned_mld[df_full$Anomaly == 0],
#   alternative = "less"
# )$statistic
# 
# ks_n2_obs <- ks.test(
#   df_full$Max_N2[df_full$Anomaly == 1],
#   df_full$Max_N2[df_full$Anomaly == 0],
#   alternative = "greater"
# )$statistic
# 
# # Bootstrap loop
# for (b in seq_len(B)) {
#   # 1. Resample row indices (with replacement)
#   idx_boot <- sample(seq_len(n), size = n, replace = TRUE)
#   df_boot  <- df_full[idx_boot, ]
#   
#   # 2. Subset MLD in sub vs. no-sub
#   sub_mld    <- df_boot$cleaned_mld[df_boot$Anomaly == 1]
#   no_sub_mld <- df_boot$cleaned_mld[df_boot$Anomaly == 0]
#   
#   # 3. Subset N2 in sub vs. no-sub
#   sub_n2    <- df_boot$Max_N2[df_boot$Anomaly == 1]
#   no_sub_n2 <- df_boot$Max_N2[df_boot$Anomaly == 0]
#   
#   # 4. Compute the KS statistics in this bootstrap sample
#   #    (We'll just use the default two-sided test for the *statistic*,
#   #     because we mainly care about the *value* of the KS distance.)
#   ks_mld[b] <- suppressWarnings(
#     ks.test(sub_mld, no_sub_mld)$statistic
#   )
#   ks_n2[b] <- suppressWarnings(
#     ks.test(sub_n2, no_sub_n2)$statistic
#   )
# }
# 
# bootstrapped <- tibble(ks_mld, ks_n2) %>%
#   pivot_longer(
#     cols = everything(),
#     names_to = "Parameter",
#     values_to = "KS_Distance"
#   ) %>%
#   mutate(
#     Parameter = recode(
#       Parameter,
#       "ks_mld" = "Mixed-Layer Depth (MLD)",
#       "ks_n2"  = "Brunt-Väisälä Frequency (N²)"
#     )
#   ) %>%
#   ggplot() +
#   geom_density(aes(x = KS_Distance, fill = Parameter), alpha = 0.7) +
#   labs(
#     title = "Estimated Density of Simulated KS Distances",
#     x = "Bootstrapped Kolmogorov–Smirnov Distances",
#     y = "Density",
#     fill = NULL
#   ) + scale_fill_viridis_d(begin=0.5)+
#   theme_bw(base_size = 25) +
#   theme(
#     panel.border   = element_rect(color = "black", fill = NA, size = 0.5),
#     axis.ticks     = element_line(color = "black"),
#     plot.title     = element_text(face = "bold", size = 25, hjust = 0.5),
#     legend.position = "bottom",
#     legend.text    = element_text(size = 25)
#   )# For instance, compute the fraction of times that MLD's KS < N2's KS:
# 
# bootstrapped <- bootstrapped + theme_map()
# 
# ggsave("../pubfig/bootstrapped_ks.png",bootstrapped,width = 8, height = 8, dpi = 300)

# Figure 5 map of annual mean flux due to the ESP + characteristic velocity 
sensitivity_results <- readRDS(file = "../data/sensitivity_res.Rds")
sensitivity_region <- readRDS(file = "../data/sensitivity_region.Rds")

sensitivity_region <- sensitivity_region %>% pivot_longer(cols=!c(w,region))

sensitivity_region <- sensitivity_region %>% mutate(
  Metric = case_when(
    grepl("fifty_yrs", name) ~ "Export >50 yrs",
    grepl("total", name) ~ "Total export",
    TRUE ~ name
  ),
  Export_Type = case_when(
    grepl("_bl$", name) ~ "Baseline",
    grepl("_up$", name) ~ "Upper",
    TRUE ~ "Other"
  )
)

w200_df <- sensitivity_region %>% group_by(region,Metric) %>% filter(w == 200) %>%
  filter(Export_Type == "Baseline") %>% filter(!is.na(region))

tot_expw200 <- w200_df %>% filter(!is.na(region )) %>% filter(Export_Type == "Baseline") %>%
  filter(Metric == "Total export") %>% pull(value) %>% sum()
seq_expw200 <- w200_df %>% filter(!is.na(region )) %>% filter(Export_Type == "Baseline") %>%
  filter(Metric == "Export >50 yrs") %>% pull(value) %>% sum()

w200_df <- w200_df %>% mutate(
  percentage = case_when(
    Metric == "Total export" ~ 100* value/tot_expw200,
    Metric == "Export >50 yrs" ~ 100* value/seq_expw200,
    TRUE ~ NA_real_  # default case
  )
) %>% select(!Export_Type)

w200_df %>% filter(Metric == "Total export") %>% pull(percentage) %>% sum()
w200_df %>% filter(Metric == "Export >50 yrs") %>% pull(percentage) %>% sum()



# Extract NA values for each metric and export type combination
na_values <- w200_df %>%
  filter(is.na(region)) %>%
  select(Metric, Export_Type, na_value = value)


na_values$region <- "North Atlantic"

w200_df <- w200_df %>%
  left_join(na_values, by = c("region", "Metric", "Export_Type")) %>%
  mutate(
    value = case_when(
      !is.na(na_value) ~ value + na_value,  # Sum where na_value exists
      TRUE ~ value  # Keep original value otherwise
    )
  ) %>%
  select(-na_value)  # Remove the temporary na_value column

w200_df

sensitivity_region <- sensitivity_region %>% group_by(region,Metric) %>% mutate(min_export = min(value),
                                                        max_export = max(value),
                                                        mid_point = (min_export+max_export)/2) %>% ungroup()


sensitivity_region_sel <- sensitivity_region %>% select(region,Metric,min_export,max_export) %>%
  unique() %>% na.omit()

sensitivity_region_sel 

sensitivity_region_sel %>% group_by(Metric) %>% summarize(global_min = sum(min_export),global_max=sum(max_export))

na_values <- sensitivity_region_sel %>%
  filter(region == "Mediterranean") %>%
  select(Metric, region, min_export = min_export,max_export=max_export)


na_values$region <- "North Atlantic"

sensitivity_region_sel_updated <- sensitivity_region_sel %>%
  left_join(na_values, by = c("region", "Metric")) %>%
  mutate(
    min_export = case_when(
      !is.na(min_export.y) ~ min_export.x + min_export.y,  # Sum where na values exist
      TRUE ~ min_export.x  # Keep original value otherwise
    ),
    max_export = case_when(
      !is.na(max_export.y) ~ max_export.x + max_export.y,  # Sum where na values exist
      TRUE ~ max_export.x  # Keep original value otherwise
    ) 
  ) %>%
  select(-min_export.x, -min_export.y, -max_export.x, -max_export.y) %>%  # Remove temporary columns
  filter(region != "Mediterranean")  

sensitivity_region_sel_updated
sensitivity_region_sel %>% filter(Metric == "Export >50 yrs") %>% pull(min_export)  %>% sum()
sensitivity_region_sel %>% filter(Metric == "Export >50 yrs") %>% pull(max_export)  %>% sum()
sensitivity_region_sel %>% filter(Metric == "Total export") %>% pull(min_export)  %>% sum()
sensitivity_region_sel %>% filter(Metric == "Total export") %>% pull(max_export)  %>% sum()

global_data <- data.frame(
  region = rep("Global", 2),
  Metric = c("Export >50 yrs","Total export"),
  min_export = c(0.001434343,0.005042196),
  max_export = c(0.07034874, 0.279184)
)

combined_data <- bind_rows(sensitivity_region_sel_updated, global_data)
combined_data$value <- combined_data$min_export*10

pacific_tropics_combined <- combined_data %>%
  filter(region %in% c("North Pacific", "Northern Tropics", "Southern Tropics")) %>%
  group_by(Metric) %>%
  summarise(
    min_export = sum(min_export, na.rm = TRUE),
    max_export = sum(max_export, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(region = "North Pacific & Tropics")

# Keep other regions as they are
other_regions <- combined_data %>%
  filter(!region %in% c("North Pacific", "Northern Tropics", "Southern Tropics",NA))

# Combine the summed Pacific & Tropics with other regions
final_aggregated_data <- bind_rows(other_regions, pacific_tropics_combined)
final_aggregated_data$value <- final_aggregated_data$min_export*10


final_aggregated_data$Metric <- final_aggregated_data$Metric %>% as_factor()

final_aggregated_data

combined_data


esp_plot <- ggplot(
  final_aggregated_data %>% filter(region == "Global"),
  aes(
    x= Metric,           # x
    Mean,             # bar mid-point
    ymin = min_export,
    ymax = max_export,
    y  = value,
    color = Metric)
) + geom_point(position = position_dodge2(width = 0.8, preserve = "single"),size=6)+
  geom_errorbar(linewidth=2,
    width    = 0.8,
    position = position_dodge2(width = 0.6, preserve = "single")
  ) + scale_color_manual(
    values = c(
      "Total export"                       = "#404788FF",  # teal-ish
      "Export >50 yrs"     = "#55C667FF"),
    name = NULL
  )+
  ## SAME y-axis scale as esp_plot
  scale_y_continuous(
    breaks = c(0,0.1,0.2,0.3),limits=c(0,0.3)) +
  
  labs(
    title = "a • Global Export  ",
    x     = NULL,y = "Export Rate (Pg C year⁻¹)" ) +
  
  theme_map(base_size = 20)  + theme_bw(base_size = 20) +
  theme(
    panel.border    = element_rect(colour = "black", fill = NA, size = 0.6),
    axis.ticks      = element_line(colour = "black", size = 0.3),
    plot.title      = element_text(hjust = 0.5, face = "bold",
                                   size  = 20,
                                   margin = margin(b = 6)),
    legend.position = "none",
    legend.text     = element_text(size = 20 ),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width  = unit(2.2 ,  "cm"),
    plot.margin       = margin(2, 2, 2, 2),
    axis.text.y       = element_text(size = 18),
    axis.text.x       = element_text(size = 18,angle = 15, hjust = 1),  # tilt if names are long
  )

# Compute relative contributions
plot_data <- final_aggregated_data %>%
  group_by(Metric) %>%
  mutate(global_total = value[region == "Global"],
         rel_contrib  = value / global_total) %>%
  ungroup() %>%
  filter(region != "Global")  # exclude the global rows, only regions


# Compute relative contributions
plot_data_all <- combined_data %>%
  group_by(Metric) %>%
  mutate(global_total = value[region == "Global"],
         rel_contrib  = value / global_total) %>%
  ungroup() %>%
  filter(region != "Global")  # exclude the global rows, only regions



# Stacked barplot
esp_relative_plot <- ggplot(plot_data,
                            aes(x = Metric, y = rel_contrib, fill = region)) +
  geom_bar(stat = "identity", position = "stack",width=.5) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "b • Regional contributions",
    x     = NULL,
    y     = "Fraction of Global Export"
  ) +   scale_fill_viridis_d(option = "D",
                          begin  = 0,   # skip dark purple
                          end    = 0.7,   # stop before yellow
                          name   = NULL) +
  theme_bw(base_size = 20) +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, size = 0.6),
    axis.ticks        = element_line(colour = "black", size = 0.3),
    plot.title        = element_text(hjust = 0.5, face = "bold",
                                     size  = 20,
                                     margin = margin(b = 6)),
    legend.position   = "bottom",
    legend.text       = element_text(size = 20),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width  = unit(0.5 , "cm"),
    axis.text.x       = element_text(size = 18, angle = 15, hjust = 1),
    axis.text.y       = element_text(size = 18)
  )

plot_data_all$Metric <- plot_data_all$Metric %>% as_factor()

# Stacked barplot
esp_relative_plot_all <- ggplot(plot_data_all,
                            aes(x = Metric, y = rel_contrib, fill = region)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "b • Regional contributions to Global ESP Export",
    x     = NULL,
    y     = "Fraction of Global Export"
  ) +   scale_fill_viridis_d() +
  theme_bw(base_size = 18) +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, size = 0.6),
    axis.ticks        = element_line(colour = "black", size = 0.3),
    plot.title        = element_text(hjust = 0.5, face = "bold",
                                     size  = 25,
                                     margin = margin(b = 6)),
    legend.position   = "bottom",
    legend.text       = element_text(size = 18),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width  = unit(2.2 , "cm"),
    axis.text.x       = element_text(size = 18, angle = 15, hjust = 1),
    axis.text.y       = element_text(size = 18)
  )






# File too large for github
merged_df<- readRDS("/data/GLOBARGO/src/data/merged_dataset_poc_estim.Rds")
#merged_df<- readRDS("../../merged_dataset_poc_estim.Rds")

merged_df$LATITUDE %>% summary()


# POC flux assuming T = 4 days (mean speed of 50 m/s)
merged_df <- merged_df %>% mutate(poc_flux20w_bl = poc_bl * proportion * 20/200,  # lower bound
                                  poc_flux20w_up = poc_up * proportion * 20/200,
                                  poc_flux200w_bl = poc_bl * proportion * 200/200,
                                  poc_flux500w_bl = poc_bl * proportion * 500/200,
                                  poc_flux500w_up = poc_up * proportion * 500/200 # upper bound
)




merged_df$poc_flux200w_bl %>% summary()
merged_df$poc_flux500w_bl %>% summary()
# Get the daily mean fluxes

annual_mean_region <- merged_df %>%
  group_by(LONGITUDE, LATITUDE) %>%
  summarize(poc_flux_annual = mean(poc_flux200w_bl, na.rm = TRUE),
            poc_flux_annual_top = mean(poc_flux500w_up, na.rm = TRUE),
            poc_flux_annual_bot = mean(poc_flux20w_bl, na.rm = TRUE),
            poc_bl_annual = mean(poc_bl, na.rm = TRUE)
            )

stipple_resolution <- 5

argo_bins <- df_argo_clean %>%
  mutate(
    lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
    lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
  ) %>%
  group_by(lon_bin_stipple, lat_bin_stipple) %>%
  summarize(count = n(), .groups = "drop")

full_grid <- expand.grid(
  lon_bin_stipple = seq(-180, 175, by = stipple_resolution),
  lat_bin_stipple = seq(-90, 90, by = stipple_resolution)
) %>%
  as_tibble()

# Fill in any missing bins with count = 0
argo_bins_full <- full_grid %>%
  left_join(argo_bins, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
  mutate(count = ifelse(is.na(count), 0, count))

# 2. Identify undersampled areas (e.g., fewer than 1 profile)
undersampled <- argo_bins_full %>% filter(count < 1)

# 3. Generate corner points for each undersampled grid cell (4 corners per cell)
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

# 4. Filter the prediction grid to the desired season
world <- ne_countries(scale = "medium", returnclass = "sf")

annual_mean_region$poc_flux_annual_NA <- annual_mean_region$poc_flux_annual
annual_mean_region$poc_flux_annual_NA <- ifelse(annual_mean_region$poc_flux_annual <= 0.05,NA,annual_mean_region$poc_flux_annual)

annual_mean_region$poc_flux_annual_bot_NA <- annual_mean_region$poc_flux_annual_bot
annual_mean_region$poc_flux_annual_bot_NA <- ifelse(annual_mean_region$poc_flux_annual_bot <= 0.005,NA,annual_mean_region$poc_flux_annual_bot)

annual_mean_region$poc_flux_annual_top_NA <- annual_mean_region$poc_flux_annual_top
annual_mean_region$poc_flux_annual_top_NA <- ifelse(annual_mean_region$poc_flux_annual_top <= 0.125,NA,annual_mean_region$poc_flux_annual_top)


annual_mean_region$poc_flux_annual_top_NA %>% summary()
annual_mean_region$poc_flux_annual_bot_NA %>% summary()

poc_flux_map <- ggplot() +
  ## 1 ─────────────── POC‑flux background ──────────────────────────
  geom_tile(data = annual_mean_region,
            aes(LONGITUDE, LATITUDE, fill = poc_flux_annual_NA),
            alpha = 0.96) +
  scale_fill_viridis_b(
    name = "POC flux (mg C m⁻² day⁻¹)               ",
    breaks = c(0.05,0.1,0.5,1,2,5,7.5),
    limits = c(0.05,20 ),
    oob = scales::squish, na.value = "white") +                       # your viridis / magma scale object
  
  ## 3 ─────────────── Stippling for undersampled grid cells ────────
  geom_point(data = undersampled_corners,
             aes(LON, LAT),
             colour = "red", size = 1.3, alpha = 0.7, shape = 20) +
  
  ## 4 ─────────────── Coastlines  ──────────────────────────────────
  geom_sf(data = world, fill = "lightgrey", colour = "lightgrey",
          linewidth = 0.3) +
  
  ## 5 ─────────────── Axes, ticks, labels, theme  ──────────────────
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  scale_x_continuous(breaks = ticks_x, labels = label_lon) +
  scale_y_continuous(breaks = ticks_y, labels = label_lat) +
  labs(
    title = "c • Average annual POC flux", # (assuming W = 50 m day⁻¹)
    x = "Longitude", y = "Latitude", fill = "POC flux\n(mg C m⁻² day⁻¹)"
  ) +
  theme_map(base_size = 20) +
  theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.key.width  = unit(2.4, "cm"),
    legend.key.height = unit(0.5, "cm"))

poc_flux_map_top <- ggplot() +
  ## 1 ─────────────── POC‑flux background ──────────────────────────
  geom_tile(data = annual_mean_region,
            aes(LONGITUDE, LATITUDE, fill = poc_flux_annual_top_NA),
            alpha = 0.96) +
  scale_fill_viridis_b(
    name = "POC flux (mg C m⁻² day⁻¹)               ",
    oob = scales::squish, na.value = "white") +                       # your viridis / magma scale object
  
  ## 2 ─────────────── White contours (same variable) ───────────────
  geom_contour(data = annual_mean_region,
               aes(LONGITUDE, LATITUDE, z = poc_flux_annual_top_NA,
                   colour = after_stat(level)),
               colour  = "white",   # fixed colour, so legend disabled below
               alpha   = 0.30,
               binwidth = 2.5,
               show.legend = FALSE) +
  
  ## 3 ─────────────── Stippling for undersampled grid cells ────────
  geom_point(data = undersampled_corners,
             aes(LON, LAT),
             colour = "red", size = 1.3, alpha = 0.7, shape = 20) +
  
  ## 4 ─────────────── Coastlines  ──────────────────────────────────
  geom_sf(data = world, fill = "lightgrey", colour = "lightgrey",
          linewidth = 0.3) +
  
  ## 5 ─────────────── Axes, ticks, labels, theme  ──────────────────
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  scale_x_continuous(breaks = ticks_x, labels = label_lon) +
  scale_y_continuous(breaks = ticks_y, labels = label_lat) +
  labs(
    title = "c • Average annual POC flux", # (assuming W = 50 m day⁻¹)
    x = "Longitude", y = "Latitude", fill = "POC flux\n(mg C m⁻² day⁻¹)"
  ) +
  theme_map(base_size = 18) +
  theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.key.width  = unit(2.2, "cm"),
    legend.key.height = unit(0.5, "cm"))






## Seasonal subsets
djf_merged <- merged_df %>% filter(season == "DJF")
mam_merged <- merged_df %>% filter(season == "MAM")
jja_merged <- merged_df %>% filter(season == "JJA")
son_merged <- merged_df %>% filter(season == "SON")

## Compute global range for legend
poc_range <- range(merged_df$poc_bl, na.rm = TRUE)

## Helper function
make_poc_map <- function(data, season_label) {
  ggplot() +
    geom_tile(data = data,
              aes(LONGITUDE, LATITUDE, fill = poc_bl),
              alpha = 0.96) +
    scale_fill_viridis_b(
      n.breaks = 8,
      limits   = poc_range,        # <<< fixed scale across seasons
      oob      = scales::squish
    ) +
    geom_contour(data = data,
                 aes(LONGITUDE, LATITUDE, z = poc_bl,
                     colour = after_stat(level)),
                 colour  = "white",
                 alpha   = 0.30,
                 binwidth = 20,
                 show.legend = FALSE) +
    geom_point(data = undersampled_corners,
               aes(LON, LAT),
               colour = "red", size = 1.3, alpha = 0.7, shape = 20) +
    geom_sf(data = world, fill = "lightgrey", colour = "lightgrey",
            linewidth = 0.3) +
    coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
    scale_x_continuous(breaks = ticks_x, labels = label_lon) +
    scale_y_continuous(breaks = ticks_y, labels = label_lat) +
    labs(
      title = paste0(season_label, "  Average Integrated POC"),
      x = "Longitude", y = "Latitude", fill = "POC flux\n(mg C m⁻²)"
    ) +
    theme_map(base_size = 18) +
    theme(
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.box        = "horizontal",
      legend.key.width  = unit(2.2, "cm"),
      legend.key.height = unit(0.5, "cm"))
}

## Generate plots
poc_djf <- make_poc_map(djf_merged, "a • DJF")
poc_mam <- make_poc_map(mam_merged, "b • MAM")
poc_jja <- make_poc_map(jja_merged, "c • JJA")
poc_son <- make_poc_map(son_merged, "d • SON")

## Arrange and keep common legend
library(patchwork)
poc_all <- (poc_djf + poc_mam) / (poc_jja + poc_son) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom") +
  shrink_margins

ggsave(filename = "../pubfig/poc_in_si.png",poc_all,width = 17,height = 10,dpi = 400)



shrink_margins <- theme(plot.margin = margin(5, 5, 5, 5))   # 5 pts all round
esp_plot        <- esp_plot        + shrink_margins
poc_flux_map    <- poc_flux_map    + shrink_margins
esp_relative_plot <- esp_relative_plot + shrink_margins

(esp_plot | esp_relative_plot) / poc_flux_map + plot_layout(ncol = 2, widths = c(1, 2.2))

combined_fig <- (esp_plot | esp_relative_plot) / poc_flux_map + plot_layout(ncol = 2, widths = c(1, 2.2))


ggsave(plot = combined_fig,
       "../pubfig/fig5_ESP_export_rate_with_w_and_avg_POC_flux_map.png",
       width = 24,height = 8,dpi = 400)

# Sensitivity of annual export to W : 
df <- sensitivity_results
df$total_export_PgC_mean <- (df$total_export_PgC_bl + df$total_export_PgC_up)/2
df$total_export_PgC_50yrs_mean <- (df$total_export_PgC_50yrs_bl + df$total_export_PgC_50yrs_up)/2
# 1) Prepare the tidy data
df_plot <- bind_rows(
  df %>% 
    transmute(
      w_days,
      metric = "Total",
      min    = total_export_PgC_bl,
      mean   = total_export_PgC_mean,
      max    = total_export_PgC_up
    ),
  df %>% 
    transmute(
      w_days,
      metric = "Sequestered for >50 yrs",
      min    = total_export_PgC_50yrs_bl,
      mean   = total_export_PgC_50yrs_mean,
      max    = total_export_PgC_50yrs_up
    )
)


df_plot %>% group_by(metric) %>% summarise(mean=mean(mean))

sensitivity_plot <- ggplot(df_plot,
                           aes(x=w_days, y=mean,
                               ymin = min, ymax = max,
                               colour = metric)) +
  geom_errorbar(position = position_dodge(width = 0.6),
                fatten = 0.1, size = 1)+
  ## Viridis with the yellow–green range removed
  scale_colour_viridis_d(option = "D",
                         begin  = 0.15,   # skip dark purple
                         end    = 0.85,   # stop before yellow
                         name   = NULL) + scale_y_log10()+
  labs(
    title = "Sensitivity of Global Export of the Eddy Subduction Pump to Vertical Velocity",
    x = "Characteristic Vertical Velocity W (m day⁻¹)",
    y = "Global Export (Pg C year⁻¹)"
  ) +scale_x_continuous(breaks = c(20, 50, 100, 200, 300, 400, 500))+
  
  theme_map(base_size = 18) +
  theme(
    axis.text.x       = element_text(angle = 0, vjust = 0.5),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.key.width  = unit(2.2, "cm"),
    legend.key.height = unit(0.5, "cm")
  )

ggsave(plot = sensitivity_plot,filename = "../pubfig/fig5_SI_sensitivity_global_export.png",width = 15,height = 10)


############################
## Fig 7 nice profile plot
##########################

df <- readRDS(file = "/data/GLOBARGO/src/data/nice_profile_data.Rds") 

# ────────────────────────────────────────────────────────────────
# 0.  Common palette & theme  (same sizes as your maps)
# ────────────────────────────────────────────────────────────────
pal_prof <- viridis::viridis(2, option = "D", end = 0.8)
names(pal_prof) <- c("Observed", "Interpolated")

theme_profile <- theme_map(base_size = 18) %+replace%
  
  
  
  theme(legend.position   = "bottom",
        legend.direction  = "horizontal",
        legend.key.width  = unit(2.2, "cm"),
        legend.key.height = unit(0.5, "cm"),
        axis.text.x       = element_text(angle = 45, hjust = 1))

# helper for the x‑axis (pressure) so all four plots are identical
scale_pressure <- list(
  coord_flip(),
  scale_x_reverse(limits = c(1000, 0), expand = c(0, 0)),
  labs(x = "Adjusted pressure (dbar)")
)

# ────────────────────────────────────────────────────────────────
# 1.  Build each panel
# ────────────────────────────────────────────────────────────────
POC_plot <- ggplot(df, aes(PRES_ADJUSTED)) +  
  geom_vline(color="red",xintercept = 940,alpha=.1,size=20)+
  geom_line(aes(y = POC_OBS,        colour = "Observed"),     size = 1.4) +
  geom_line(aes(y = predicted_POC,  colour = "Interpolated"), size = 1.0) +
  scale_colour_manual(values = pal_prof, name = NULL) +
  scale_pressure +
  labs(y = "Particulate organic carbon (mg m⁻³)",
       title = paste0("d • POC profile — float ", wmo, " cycle ", cycle_number)) +
  theme_profile

bbp_plot <- ggplot(df, aes(PRES_ADJUSTED)) +
  geom_vline(color="red",xintercept = 940,alpha=.1,size=20)+
  geom_line(aes(y = BBP700_ADJUSTED, colour = "Observed"), size = 1.4) +
  scale_colour_manual(values = pal_prof["Observed"], drop = FALSE, name = NULL) +
  scale_pressure +
  labs(y = expression("bbp(700 nm) (m"^{-1}*" sr"^{-1}*")"),
       title = paste0("c • Backscatter — float ", wmo, " cycle ", cycle_number)) +
  theme_profile

AOU_plot <- ggplot(df, aes(PRES_ADJUSTED)) +
  geom_vline(color="red",xintercept = 940,alpha=.1,size=20)+
  geom_line(aes(y = AOU, colour = "Observed"), size = 1.4) +
  scale_colour_manual(values = pal_prof["Observed"], drop = FALSE, name = NULL) +
  scale_pressure +
  labs(y = "Apparent oxygen utilisation (µmol kg⁻¹)",
       title = paste0("a • AOU — float ", wmo, " cycle ", cycle_number)) +
  theme_profile

ABS_SAL_plot <- ggplot(df, aes(PRES_ADJUSTED)) +
  geom_vline(color="red",xintercept = 940,alpha=.1,size=20)+
  geom_line(aes(y = ABS_SAL, colour = "Observed"), size = 1.4) +
  scale_colour_manual(values = pal_prof["Observed"], drop = FALSE, name = NULL) +
  scale_pressure +
  labs(y = "Absolute salinity (g kg⁻¹)",
       title = paste0("b • Absolute salinity — float ", wmo, " cycle ", cycle_number)) +
  theme_profile

# ────────────────────────────────────────────────────────────────
# 2.  Combine & export
# ────────────────────────────────────────────────────────────────
profile_grid <- (AOU_plot | ABS_SAL_plot) /
  (bbp_plot | POC_plot) +
  plot_layout(guides = "collect") &   # single legend row
  theme(legend.position = "bottom")

profile_grid <- (AOU_plot | ABS_SAL_plot + POC_plot)+
  plot_layout(guides = "collect") &   # single legend row
  theme(legend.position = "bottom")

########################################
#########################################
##################### JUNK CODE
 #

# df_esp <- tibble(
#   Study  = c(
#     "This Study",
#     "This Study",
#     "(Boyd et al, 2019)"
#   ),
#   Metric = c(
#     "Total export",
#     "Export Sequestered for >50 yrs",
#     "Total export"
#   ),
#   Min   = c(0.005111526, 0.0015, 0),
#   Mean  = c(0.0921,      0.0239,      1),
#   Max   = c(0.282916640, 0.070641275, 2)
# )
# 
# # make sure Metric is a factor in the order you want
# df_esp$Metric <- factor(
#   df_esp$Metric,
#   levels = c(
#     "Total export",
#     "Export Sequestered for >50 yrs"
#   )
# )
# scale_boyd <- scale_y_continuous(limits = c(0, 2),
#                                  breaks = c(0,0.1,0.2,0.3,0.5,1,2),
#                                  labels=c("0","0.1","0.2","0.3","0.5","1","2"))
# 
# esp_plot_boyd <- ggplot(df_esp,
#                         aes(Metric, Mean,
#                             ymin = Min, ymax = Max,
#                             fill = Study)) +
#   geom_crossbar(                    # rectangle + a mid-bar
#     width   = 0.5,                # horizontal width of each bar
#     fatten  = 0,                   # kills the “mean” point
#     position = position_dodge2(width = 0.6,         # side-by-side spacing
#                                preserve = "single")
#   ) +
#   
#   ## Viridis with the yellow–green range removed
#   scale_fill_manual(
#     values = c(
#       "(Boyd et al, 2019)"                       = "#55C667FF",  # teal-ish
#       "This Study"     = "#404788FF"   # amber-ish
#     ))+  
#   ## y‑axis from 0 to max + headroom
#   scale_boyd+  
#   labs(
#     title = "a • Global Export",
#     x = NULL,
#     y = "Export (Pg C year⁻¹)"
#   ) +
#   
#   theme_map(base_size = 18) +
#   theme(
#     axis.text.x       = element_text(angle = 45, vjust = 0.5),
#     legend.position   = "bottom",
#     legend.direction  = "horizontal",
#     legend.key.width  = unit(2.2, "cm"),
#     legend.key.height = unit(0.5, "cm"),
#     legend.title = element_blank()
#   )
# 
# 
# esp_plot_log10 <- ggplot(df_esp,
#                          aes(Metric, Mean,
#                              ymin = Min, ymax = Max,
#                              fill = Study)) +
#   geom_crossbar(                    # rectangle + a mid-bar
#     width   = 0.5,                # horizontal width of each bar
#     fatten  = 0,                   # kills the “mean” point
#     position = position_dodge2(width = 0.6,         # side-by-side spacing
#                                preserve = "single")
#   ) +
#   
#   ## Viridis with the yellow–green range removed
#   scale_fill_viridis_d(option = "D",
#                        begin  = 0.15,   # skip dark purple
#                        end    = 0.85,   # stop before yellow
#                        name   = NULL) +
#   
#   ## y‑axis from 0 to max + headroom
#   scale_y_log10(breaks=c(0.0015,0.025,0.1,1,2),labels=c("0.0015","0.025","0.1","1","2")) +
#   
#   labs(
#     title = "a • Global Export of the Eddy Subduction Pump",
#     x = NULL,
#     y = "Global Export (Pg C year⁻¹)"
#   )+ 
#   
#   theme_map(base_size = 18) +
#   theme(
#     axis.text.x       = element_text(angle = 0, vjust = 0.5),
#     legend.position   = "bottom",
#     legend.direction  = "horizontal",
#     legend.key.width  = unit(2.2, "cm"),
#     legend.key.height = unit(0.5, "cm")
#   )
# df_esp$Metric
# 
# esp_plot <- ggplot(df_esp %>% filter(Study=="This Study"),
#                    aes(Metric, Mean,
#                        ymin = Min, ymax = Max,
#                        fill = Metric)) +
#   geom_crossbar(                    # rectangle + a mid-bar
#     width   = 0.5,                # horizontal width of each bar
#     fatten  = 0,                   # kills the “mean” point
#     position = position_dodge2(width = 0.6,         # side-by-side spacing
#                                preserve = "single")
#   ) +
#   
#   ## Viridis with the yellow–green range removed
#   scale_fill_manual(
#     values = c(
#       "Total export"                       = "#404788FF",  # teal-ish
#       "Export Sequestered for >50 yrs"     = "#55C667FF"   # amber-ish
#     ),
#     name = NULL
#   )+
#   labs(
#     title = "a • Global Export",
#     x = NULL,
#     y = "Export (Pg C year⁻¹)"
#   ) + scale_y_continuous(limits = c(0, 0.30),
#                          breaks = c(0.001,0.05,0.1,0.2,0.3),
#                          labels=c("0.001","0.05","0.1","0.2","0.3")) +
#   
#   theme_map(base_size = 18) +
#   theme(
#     axis.text.x       = element_text(angle = 0, vjust = 0.5),
#     legend.position   = "none",
#     legend.direction  = "horizontal",
#     legend.key.width  = unit(2.2, "cm"),
#     legend.key.height = unit(0.5, "cm")
#   )
# 
# 
# 
# 
# # Make sure Mean exists in the regional data frame
# sensitivity_region <- sensitivity_region %>%
#   mutate(Mean = 0.5 * (max_export + min_export))
# 
# # ─────────────────────────────────────────────────────────────
# # 1.  Build the regional plot with harmonised styling
# # ─────────────────────────────────────────────────────────────
# sensitivity_region
# # names of the bands you want to merge
# combo <- c("Northern Tropics", "Southern Tropics", "North Pacific")
# 
# sensitivity_region2 <- sensitivity_region %>% 
#   ## 1. build the new combined row
#   filter(region %in% combo) %>% 
#   summarise(
#     region      = "Tropics and North Pacific",
#     min_export  = sum(min_export,  na.rm = TRUE),
#     max_export  = sum(max_export,  na.rm = TRUE)
#   ) %>% 
#   mutate(Mean = 0.5 * (min_export + max_export)) %>%        # same formula
#   ## 2. stick it onto the original data
#   bind_rows(sensitivity_region) %>% 
#   ## 3. reorder if you like (optional)
#   relocate(region, min_export, max_export, Mean)
# 
# sensitivity_region2 <- sensitivity_region2 %>% filter(region %in% c("North Atlantic","Tropics and North Pacific","Southern Ocean"))
# 
# sensitivity_region2$Metric <- "Total export"
# sensitivity_region2$Study  <- "This Study"
# sensitivity_region2$Min <-  sensitivity_region2$min_export
# sensitivity_region2$Max <-  sensitivity_region2$max_export
# 
# df_esp_region <- bind_rows(df_esp,sensitivity_region2) %>% filter(Study == "This Study") %>% select(Metric,Min,Mean,Max,region)
# df_esp_region[1,5] <- "Global Export"
# df_esp_region[2,5] <- "Global Export Sequestered for >50 yrs"
# df_esp_region$region <- df_esp_region$region %>% as_factor()
# 
# 
# # re-factor with the new level order
# df_esp_region$region <- factor(
#   df_esp_region$region,
#   levels = c("Global Export","North Atlantic","Tropics and North Pacific","Southern Ocean","Global Export Sequestered for >50 yrs")
# )
# 
# esp_plot <- ggplot(
#   df_esp_region,
#   aes(
#     x= region,           # x
#     Mean,             # bar mid-point
#     ymin = Min,
#     ymax = Max,
#     fill = Metric)
# ) +
#   geom_crossbar(
#     width    = 0.5,
#     fatten   = 0,
#     position = position_dodge2(width = 0.6, preserve = "single")
#   ) + scale_fill_manual(
#     values = c(
#       "Total export"                       = "#404788FF",  # teal-ish
#       "Export Sequestered for >50 yrs"     = "#55C667FF"   # amber-ish
#     ),
#     name = NULL
#   )+
#   ## SAME y-axis scale as esp_plot
#   scale_y_continuous(
#     breaks = c(0,0.1,0.2,0.3),limits=c(0,0.3)) +
#   
#   labs(
#     title = "a • Export of the Eddy Subduction Pump",
#     x     = NULL,y = "Export Rate (Pg C year⁻¹)" ) +
#   
#   theme_map(base_size = 18) +
#   theme(
#     axis.text.x       = element_text(angle = 45, hjust = 1),  # tilt if names are long
#     legend.position   = "none",
#     legend.direction  = "horizontal",
#     legend.key.width  = unit(2.2, "cm"),
#     legend.key.height = unit(0.5, "cm")
#   )
# 
# esp_plot_log10 <- ggplot(
#   df_esp_region,
#   aes(
#     x= region,           # x
#     Mean,             # bar mid-point
#     ymin = Min,
#     ymax = Max,
#     fill = Metric)
# ) +
#   geom_crossbar(
#     width    = 0.5,
#     fatten   = 0,
#     position = position_dodge2(width = 0.6, preserve = "single")
#   ) + scale_fill_manual(
#     values = c(
#       "Total export"                       = "#404788FF",  # teal-ish
#       "Export Sequestered for >50 yrs"     = "#55C667FF"   # amber-ish
#     ),
#     name = NULL
#   )+
#   ## SAME y-axis scale as esp_plot
#   scale_y_log10() +
#   
#   labs(
#     title = "a • Export of the Eddy Subduction Pump",
#     x     = NULL,y = "Export Rate (Pg C year⁻¹)" ) +
#   
#   theme_map(base_size = 18) +
#   theme(
#     axis.text.x       = element_text(angle = 45, hjust = 1),  # tilt if names are long
#     legend.position   = "none",
#     legend.direction  = "horizontal",
#     legend.key.width  = unit(2.2, "cm"),
#     legend.key.height = unit(0.5, "cm"),
#     axis.ticks      = element_line(colour = "black", size = 0.3),
#     axis.text       = element_text(size = base_size * 0.8),
#     plot.title      = element_text(hjust = 0.5, face = "bold",
#                                    size  = base_size * 1.05,
#                                    margin = margin(b = 6)),
#   )
# 
# 
# # File too large for github
# merged_df<- readRDS("data/merged_dataset_poc_estim.Rds")
# 
# merged_df$LATITUDE %>% summary()
# 
# 
# # POC flux assuming T = 4 days (mean speed of 50 m/s)
# merged_df <- merged_df %>% mutate(poc_flux20w_bl = poc_bl * proportion * 20/200,  # lower bound
#                                   poc_flux20w_up = poc_up * proportion * 20/200,
#                                   poc_flux200w_bl = poc_bl * proportion * 200/200,
#                                   poc_flux500w_bl = poc_bl * proportion * 500/200,
#                                   poc_flux500w_up = poc_up * proportion * 500/200 # upper bound
# )
# 
# 
# 
# 
# merged_df$poc_flux200w_bl %>% summary()
# merged_df$poc_flux500w_bl %>% summary()
# # Get the daily mean fluxes
# annual_mean_region <- merged_df %>%
#   group_by(LONGITUDE, LATITUDE) %>%
#   summarize(poc_flux_annual_bl = mean(poc_flux20w_bl, na.rm = TRUE),
#             poc_flux_annual_up = mean(poc_flux500w_up, na.rm = TRUE))
# 
# annual_mean_region
# 
# stipple_resolution <- 5
# 
# argo_bins <- df_argo_clean %>%
#   mutate(
#     lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
#     lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
#   ) %>%
#   group_by(lon_bin_stipple, lat_bin_stipple) %>%
#   summarize(count = n(), .groups = "drop")
# 
# full_grid <- expand.grid(
#   lon_bin_stipple = seq(-180, 175, by = stipple_resolution),
#   lat_bin_stipple = seq(-90, 90, by = stipple_resolution)
# ) %>%
#   as_tibble()
# 
# # Fill in any missing bins with count = 0
# argo_bins_full <- full_grid %>%
#   left_join(argo_bins, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
#   mutate(count = ifelse(is.na(count), 0, count))
# 
# # 2. Identify undersampled areas (e.g., fewer than 1 profile)
# undersampled <- argo_bins_full %>% filter(count < 1)
# 
# # 3. Generate corner points for each undersampled grid cell (4 corners per cell)
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
# # 4. Filter the prediction grid to the desired season
# world <- ne_countries(scale = "medium", returnclass = "sf")
# 
# 
# annual_mean$poc_flux_annual_NA <- ifelse(df_plot$poc_flux_annual <= 0.05,NA,df_plot$poc_flux_annual)
# 
# poc_flux_map <- ggplot() +
#   ## 1 ─────────────── POC‑flux background ──────────────────────────
#   geom_tile(data = annual_mean,
#             aes(LONGITUDE, LATITUDE, fill = poc_flux_annual_NA),
#             alpha = 0.96) +
#   scale_fill_viridis_b(
#     name = "POC flux (mg C m⁻² day⁻¹)",
#     breaks = c(0.1,0.5,1,2,5,7.5),
#     limits = c(0.5,20 ),
#     oob = scales::squish, na.value = "white") +                       # your viridis / magma scale object
#   
#   ## 2 ─────────────── White contours (same variable) ───────────────
#   geom_contour(data = df_plot,
#                aes(LONGITUDE, LATITUDE, z = poc_flux_annual,
#                    colour = after_stat(level)),
#                colour  = "white",   # fixed colour, so legend disabled below
#                alpha   = 0.30,
#                binwidth = 2.5,
#                show.legend = FALSE) +
#   
#   ## 3 ─────────────── Stippling for undersampled grid cells ────────
#   geom_point(data = undersampled_corners,
#              aes(LON, LAT),
#              colour = "red", size = 1.3, alpha = 0.7, shape = 20) +
#   
#   ## 4 ─────────────── Coastlines  ──────────────────────────────────
#   geom_sf(data = world, fill = "lightgrey", colour = "lightgrey",
#           linewidth = 0.3) +
#   
#   ## 5 ─────────────── Axes, ticks, labels, theme  ──────────────────
#   coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
#   scale_x_continuous(breaks = ticks_x, labels = label_lon) +
#   scale_y_continuous(breaks = ticks_y, labels = label_lat) +
#   labs(
#     title = "b • Average annual POC flux", # (assuming W = 50 m day⁻¹)
#     x = "Longitude", y = "Latitude", fill = "POC flux\n(mg C m⁻² day⁻¹)"
#   ) +
#   theme_map(base_size = 18) +
#   theme(
#     legend.position   = "bottom",
#     legend.direction  = "horizontal",
#     legend.box        = "horizontal",
#     legend.key.width  = unit(2.2, "cm"),
#     legend.key.height = unit(0.5, "cm"))
# 
# shrink_margins <- theme(plot.margin = margin(5, 5, 5, 5))   # 5 pts all round
# esp_plot        <- esp_plot        + shrink_margins
# esp_region_plot <- esp_region_plot + shrink_margins
# poc_flux_map    <- poc_flux_map    + shrink_margins
# 
# 
# 
# combined_fig <- (esp_plot |  poc_flux_map) +
#   plot_layout(widths = c(1, 2),guides = "auto")
# 
# combined_fig_boyd <- (esp_plot_boyd |  poc_flux_map) +
#   plot_layout(widths = c(1, 2),guides = "auto")
# 
# 
# combined_fig_log <- (esp_plot_log10|  poc_flux_map )+
#   plot_layout(widths = c(1, 2),guides = "auto")
# 
# ggsave(plot = combined_fig,
#        "../pubfig/fig5_ESP_export_rate_with_w_and_avg_POC_flux_map.png",
#        width = 18,height = 10,dpi = 250)
# 
# ggsave(plot = combined_fig_log,
#        "../pubfig/fig5_ESP_export_rate_with_w_and_avg_POC_flux_map_log.png",
#        width = 25,height = 10,dpi = 250)
# 
# combined_fig <- esp_plot_log10 | poc_flux_map
# ggsave(plot = combined_fig,
#        "../pubfig/fig5_ESP_export_rate_with_w_and_avg_POC_flux_map_log10.png",
#        width = 19,height = 10,dpi = 400)
# 
# 
# 
# 