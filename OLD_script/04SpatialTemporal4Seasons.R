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
# Resolve function conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")
setwd("/data/GLOBARGO/src/")
df_argo_clean <- read_csv("data/df_argo_loc.csv")
df_complete_clean <- read_csv("data/df_eddy_subduction_anom.csv")
df_carbon_clean <- read_csv("data/df_carbon_subduction_anom.csv")

all_anomalies <- read_csv("/data/GLOBARGO/data/detected_events_abs_sal_var_v5.csv")



df_carbon_with_poc <- read_csv("/data/GLOBARGO/src/data/df_carbon_subduction_anom_with_poc_fromgali.csv")


# --- Compute argo density per 1 degree bin


# Now compute the average number of profiles per 1° grid cell for each region:

total_counts <- df_counts %>% group_by(region) %>% summarize(sum(count))
# Earth's radius in kilometers
R <- 6371

# Function: compute area (km²) for a spherical rectangle given lat & lon bounds in degrees
spherical_rectangle_area <- function(lat1, lat2, lon1, lon2, R = 6371) {
  # Convert degrees to radians
  lat1_rad <- lat1 * pi / 180
  lat2_rad <- lat2 * pi / 180
  lon1_rad <- lon1 * pi / 180
  lon2_rad <- lon2 * pi / 180
  
  # Compute area using the formula:
  # A = R^2 * |lon2 - lon1| * |sin(lat2) - sin(lat1)|
  area <- R^2 * abs(lon2_rad - lon1_rad) * abs(sin(lat2_rad) - sin(lat1_rad))
  return(area)
}

# Compute approximate areas for each region based on your definitions:
areas <- tribble(
  ~region,           ~area_km2,
  "North Atlantic",  spherical_rectangle_area(30, 90, -100, 20),
  "North Pacific",   spherical_rectangle_area(30, 90, -180, -100) + spherical_rectangle_area(30, 90, 120, 180),
  "Northern Tropics", spherical_rectangle_area(0, 30, -180, 180),
  "Southern Tropics", spherical_rectangle_area(-30, 0, -180, 180),
  "Southern Ocean",  spherical_rectangle_area(-90, -30, -180, 180)
)

print(areas)
# You will see (approximate) areas such as:
#  North Atlantic: ~4.26e7 km^2 (42.6 million km²)
#  North Pacific:  ~4.97e7 km^2 (49.7 million km²)
#  Northern Tropics: ~1.28e8 km^2 (127.9 million km²)
#  Southern Tropics: ~1.28e8 km^2 (127.9 million km²)
#  Southern Ocean:  ~1.28e8 km^2 (127.9 million km²)

# Suppose you have a tibble with total profile counts per region (as given):
region_counts <- tribble(
  ~region,           ~total_count,
  "North Atlantic",  23748,
  "North Pacific",    6609,
  "Northern Tropics", 26303,
  "Southern Ocean",   45138,
  "Southern Tropics", 16914
)

# Join counts with area estimates and compute density.
# Here we compute the number of profiles per million km².
density_df <- region_counts %>%
  left_join(areas, by = "region") %>%
  mutate(density_per_million_km2 = (total_count / area_km2) * 1e6)

# --- 1) Create WEEK and YEAR columns, then define bin edges ---
df_argo_clean <- df_argo_clean %>%
  mutate(
    WEEK = week(TIME),
    YEAR = year(TIME)
  )

bin_size <- 0.1  # ~10 km at mid-latitudes
longitude_bins <- seq(
  floor(min(df_argo_clean$LONGITUDE)),
  ceiling(max(df_argo_clean$LONGITUDE)),
  by = bin_size
)
latitude_bins <- seq(
  floor(min(df_argo_clean$LATITUDE)),
  ceiling(max(df_argo_clean$LATITUDE)),
  by = bin_size
)

# Compute bin centers for labeling
lon_centers <- longitude_bins[-length(longitude_bins)] + bin_size / 2
lat_centers <- latitude_bins[-length(latitude_bins)] + bin_size / 2

# --- 2) Assign each float to its bin ---
df_argo_clean <- df_argo_clean %>%
  mutate(
    lon_bin = cut(
      LONGITUDE,
      breaks = longitude_bins,
      include.lowest = TRUE,
      labels = lon_centers
    ),
    lat_bin = cut(
      LATITUDE,
      breaks = latitude_bins,
      include.lowest = TRUE,
      labels = lat_centers
    )
  )

# --- 3) Summarize Argo profiles by week, year, and bin ---
df_argo_submeso <- df_argo_clean %>%
  group_by(WEEK, YEAR, lon_bin, lat_bin) %>%
  summarize(count_total = n(), .groups = 'drop')

# --- 4) Filter for a particular week (e.g., Week 10 of 2024) ---
df_argo_randomweek <- df_argo_submeso %>%
  filter(WEEK == 10, YEAR == 2024) %>%
  mutate(
    lon_bin_num = as.numeric(as.character(lon_bin)),
    lat_bin_num = as.numeric(as.character(lat_bin))
  )

# --- 5) Plot: red squares where there's at least one float in that bin ---
world_map <- map_data("world")

ggplot() +
  # Draw world map polygons
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "grey80", 
    color = "black"
  ) + geom_tile(
    data = df_argo_randomweek,
    aes(x = lon_bin_num, y = lat_bin_num,fill=count_total),
    alpha = 1
  )+
  # Fix aspect ratio and set lat/lon bounds if desired
  coord_fixed(1.3, xlim = c(-180, 180), ylim = c(-90, 90)) +
  # Overlay bins (squares) with Argo data
  +
  labs(
    title = "Argo Float Distribution (Week 10, 2024) in ~10 km Bins",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal()



df_argo_submeso$count_total %>% as_factor() %>% summary()
ggplot() +
  geom_sf(data = wld.rob.sf, fill = "grey", color = "gray") +
  # Plot "not analyzed yet" points first (in black)
  geom_sf(data = df_argo_black, 
          aes(geometry = geometry), 
          color = "black", 
          alpha = 0.1, size = 2) +  labs(
    title = "0.2 gridcell that have at least 1 profile",
    subtitle = "Less than 15% of the Argo database has been analyzed for subduction signatures",
    color = "Previous studies"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(face = "bold", size = 30, hjust = 0.5),
    plot.subtitle = element_text(size = 25, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 25),
    panel.grid = element_blank()
  )


# Grid in 5 by 5 
stipple_resolution <- 5
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




argo_bins_full_sf <- argo_bins_full %>% filter(!is.na(lat_bin_stipple) & !is.na(lon_bin_stipple)) %>%
  st_as_sf(coords = c("lon_bin_stipple", "lat_bin_stipple"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

ggplot() +
  
  # 1. Tile plot for yearly EKE background
  geom_tile(data = argo_bins_full, 
            aes(x = lon_bin_stipple, y = lat_bin_stipple, fill = count),
            alpha = 0.95) +
  scale_fill_viridis_c(
    option = "viridis",  
    name = "Number of Argo floats",
    guide = guide_legend(reverse = TRUE, na.translate = FALSE)
  ) +
  # 4. Add stippling for undersampled cells
  geom_point(data = undersampled_corners,
             aes(x = LON, y = LAT),
             color = "white", alpha = 0.6, size = 1.5, shape = 20) +
  
  # 5. Overlay the world map
  geom_sf(data = world,
          fill = "white", color = "white",
          inherit.aes = TRUE)


# Find distinct WMO values for floats in the Kuroshio extension (assumed region)
kuroshio_floats <- df_argo_clean %>%
  filter(LONGITUDE >= 140, LONGITUDE <= 160,
         LATITUDE >= 30, LATITUDE <= 40) %>%
  distinct(WMO) %>%
  arrange(WMO)

southern_ocean_floats <- df_argo_clean %>%
  filter(LONGITUDE >= 151, LONGITUDE <= 180,
         LATITUDE >= -65, LATITUDE <= -55) %>%
  distinct(WMO) %>%
  arrange(WMO)

southern_ocean_floats

# Load world map data and transform to Robinson projection
world <- ne_countries(scale = "medium", returnclass = "sf")
meridian <- -180
wld.new <- st_break_antimeridian(world, lon_0 = meridian)

wld.rob.sf <- st_transform(wld.new, paste("+proj=robin +lon_0=", meridian, "+k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))



df_subset_llort <- df_argo_clean %>%
  # Filter by time: between 2012-01-01 and 2016-11-30
  filter(TIME >= as.Date("2012-01-01"),
         TIME <= as.Date("2016-11-30"),
         # Keep profiles in the Southern Ocean: between 65°S and 30°S.
         # Note: In the Southern Hemisphere, latitudes are negative.
         LATITUDE >= -65,
         LATITUDE <= -50)


df_subset_llort %>% filter(WMO %in% g)

# Define a lookup vector (names are the float numbers as strings)
study_map <- c("5906511" = "Chen et al., 2021", # true value : 5904034
               "2902878" = "Chen et al., 2021", # true value : 2901556
               "5904677" = "Llort et al., 2018",
               "5904479" = "Johnson & Omand, 2021",
               "5904105" = "Lacour et al., 2023", # true value : 6901516
               "5904179" = "Lacour et al., 2023", # true value : 6901480
               "5904672" = "Chen & Schofield, 2024",
               "5904673" = "Chen & Schofield, 2024",
               "5906206" = "Chen & Schofield, 2024")

# Artificial 
new_entries <- setNames(rep("Llort et al., 2018", length(g)), g)

study_map <- c(study_map, new_entries)
# Now add a new column to df_argo_clean.
# We first convert the WMO column to character (in case it's numeric)
# and then check if it matches one of the keys in study_map.
df_argo_clean <- df_argo_clean %>%
  mutate("Previous studies" = ifelse(as.character(WMO) %in% names(study_map),
                                     study_map[as.character(WMO)],
                                     "No subduction signatures detected"))

df_argo_clean$`Previous studies` %>% unique()




# Split the data: one for "No subduction signatures detected" and one for all other studies
df_argo_black <- df_argo_clean %>% filter(`Previous studies` == "No subduction signatures detected")
df_argo_other <- df_argo_clean %>% filter(`Previous studies` != "No subduction signatures detected")

ggplot() +
  geom_sf(data = world, fill = "grey", color = "gray")+coord_sf()

mapWorld <- borders("world",colour="gray70", fill="gray70") 
ggplot() + mapWorld + coord_sf(xlim=c(-180, 180),ylim=c(-22,40))




plot_argo_distrib <- ggplot() +
  # Plot world using an sf object in geographic coordinates (EPSG:4326)
  geom_sf(data = world, fill = "grey", color = "gray") +
  # Plot "not analyzed yet" points in black
  geom_point(data = df_argo_black, 
             aes(x = LONGITUDE, y = LATITUDE, color = 'No subduction signatures detected'),
             alpha = 0.1, size = 1.5) +
  # Then overlay points that have been analyzed, colored by study
  geom_point(data = df_argo_other, 
             aes(x = LONGITUDE, y = LATITUDE, color = `Previous studies`),
             alpha = 0.8, size = 1.5) +
  # Define the color scale with custom labels
  scale_color_manual(
    values = c("No subduction signatures detected" = "black",
               "Chen et al., 2021"    = "red",
               "Llort et al., 2018"    = "blue",
               "Johnson & Omand, 2021" = "green",
               "Lacour et al., 2023"   = "purple",
               "Chen & Schofield, 2024" = "pink"),
    labels = c("No subduction signatures detected" = "No subduction signatures detected",
               "Chen et al., 2021" = "Chen et al., 2021*",
               "Llort et al., 2018" = "Llort et al., 2018",
               "Johnson & Omand, 2021" = "Johnson & Omand, 2021",
               "Lacour et al., 2023" = "Lacour et al., 2023*",
               "Chen & Schofield, 2024" = "Chen & Schofield, 2024")
  ) +
  labs(
    title = "Argo Profile Distribution: 125,826 Profiles",
    subtitle = "Less than 15% of the Argo database has been analyzed for subduction signatures",
    color = "Previous studies",
    x = 'Longitude (°E)',y = 'Latitude (°N)'
  ) +
  # Set the coordinate system to geographic and define the limits
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  # Define breaks every 30° and add degree symbols to the labels
  scale_x_continuous(breaks = seq(-180, 180, 30), 
                     labels = function(x) paste0(x, "°")) +
  scale_y_continuous(breaks = seq(-90, 90, 30), 
                     labels = function(x) paste0(x, "°")) +
  # Switch to a theme that shows axes and ticks
  theme_bw(base_size = 25) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(face = "bold", size = 30, hjust = 0.5),
    plot.subtitle = element_text(size = 25, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 25)
  )

ggsave("figures/argo_distrib.png",plot = plot_argo_distrib,width = 20,height = 10,dpi = 300)

library(ggplot2)
library(ggpubr)  # for ggarrange


all_anomalies <- read_csv("/data/GLOBARGO/data/detected_events_abs_sal_var_v5.csv")

ggplot() +
  geom_sf(data = world, fill = "black", color = "black") +
  
  geom_point(
    data = df_argo_black,
    aes(x = LONGITUDE, y = LATITUDE, color = "Argo floats"),
    alpha = 0.5, size = 1.5
  ) +
  geom_point(
    data = all_anomalies,
    aes(x = LONGITUDE, y = LATITUDE, color = "All Subduction events including false positives"),
    alpha = 0.5, size = 1.5
  ) +
  
  scale_color_manual(
    name = "",
    values = c("Argo floats" = "grey", 
               "All Subduction events including false positives" = "orange")
  ) +
  
  # Increase the legend point size, and tweak spacing
  guides(
    color = guide_legend(
      override.aes = list(size = 6),  # Make the legend points larger
      keywidth = unit(3, "lines")   # Widen the key boxes
      # label.position = "right"      # Optionally, shift label to the right of the point
    )
  ) +
  
  labs(
    title = "All Subduction Anomalies Distribution: 9,979 Profiles",
    x = "Longitude (°E)",
    y = "Latitude (°N)"
  ) +
  
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  scale_x_continuous(breaks = seq(-180, 180, 30), labels = function(x) paste0(x, "°")) +
  scale_y_continuous(breaks = seq(-90, 90, 30),  labels = function(x) paste0(x, "°")) +
  
  theme_bw(base_size = 25) +
  theme(
    panel.border  = element_rect(color = "black", fill = NA, size = 0),
    axis.ticks    = element_line(color = "black"),
    plot.title    = element_text(size = 25, hjust = 0.5),
    legend.position = "bottom",
    legend.title  = element_text(face = "bold", size = 20),
    legend.text   = element_text(size = 25),
    
    # Add spacing between legend items themselves:
    legend.box = "horizontal",
    legend.spacing.x = unit(0.8, "cm")  # Adjust as needed
  )


plot_subduct_distrib <- ggplot() +
  geom_sf(data = world, fill = "black", color = "black") +
  
  geom_point(
    data = df_argo_black,
    aes(x = LONGITUDE, y = LATITUDE, color = "Argo floats"),
    alpha = 0.5, size = 1.5
  ) +
  geom_point(
    data = df_complete_clean,
    aes(x = LONGITUDE, y = LATITUDE, color = "Subduction events"),
    alpha = 0.5, size = 1.5
  ) +
  
  scale_color_manual(
    name = "",
    values = c("Argo floats" = "grey", 
               "Subduction events" = "blue")
  ) +
  
  # Increase the legend point size, and tweak spacing
  guides(
    color = guide_legend(
      override.aes = list(size = 6),  # Make the legend points larger
      keywidth = unit(3, "lines")   # Widen the key boxes
      # label.position = "right"      # Optionally, shift label to the right of the point
    )
  ) +
  
  labs(
    title = "Subduction Anomalies Distribution: 4,377 Profiles",
    x = "Longitude (°E)",
    y = "Latitude (°N)"
  ) +
  
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  scale_x_continuous(breaks = seq(-180, 180, 30), labels = function(x) paste0(x, "°")) +
  scale_y_continuous(breaks = seq(-90, 90, 30),  labels = function(x) paste0(x, "°")) +
  
  theme_bw(base_size = 25) +
  theme(
    panel.border  = element_rect(color = "black", fill = NA, size = 0),
    axis.ticks    = element_line(color = "black"),
    plot.title    = element_text(size = 25, hjust = 0.5),
    legend.position = "bottom",
    legend.title  = element_text(face = "bold", size = 20),
    legend.text   = element_text(size = 25),
    
    # Add spacing between legend items themselves:
    legend.box = "horizontal",
    legend.spacing.x = unit(0.8, "cm")  # Adjust as needed
  )

ggsave("figures/subduct_distrib.png",
       plot = plot_subduct_distrib,
       width = 20, height = 10, dpi = 300)


ggsave("figures/subduct_distrib.png",plot = plot_subduct_distrib,width = 20,height = 10,dpi = 300)
plot_carbon_distrib <- ggplot() +
  # Black continents
  geom_sf(data = world, fill = "black", color = "black") +
  
  # Argo floats in grey
  geom_point(
    data = df_argo_black,
    aes(x = LONGITUDE, y = LATITUDE, color = "Argo floats"),
    alpha = 0.5, size = 1.5
  ) +
  # Carbon subduction events in red
  geom_point(
    data = df_carbon_clean,
    aes(x = LONGITUDE, y = LATITUDE, color = "Carbon subduction events"),
    alpha = 0.5, size = 1.5
  ) +
  
  # Manual color mapping
  scale_color_manual(
    name = "",
    values = c("Argo floats" = "grey",
               "Carbon subduction events" = "red")
  ) +
  
  # Increase legend point size and add spacing
  guides(
    color = guide_legend(
      override.aes = list(size = 6),   # Larger points in the legend
      keywidth = unit(3, "lines")   # Wider legend key boxes
    )
  ) +
  
  labs(
    title = "Carbon Subduction Anomalies Distribution: 1,333 Profiles",
    x = "Longitude (°E)",
    y = "Latitude (°N)"
  ) +
  
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  
  # Add ticks every 30°
  scale_x_continuous(
    breaks = seq(-180, 180, 30),
    labels = function(x) paste0(x, "°")
  ) +
  scale_y_continuous(
    breaks = seq(-90, 90, 30),
    labels = function(x) paste0(x, "°")
  ) +
  
  theme_bw(base_size = 25) +
  theme(
    panel.border  = element_rect(color = "black", fill = NA, size = 0),
    axis.ticks    = element_line(color = "black"),
    plot.title    = element_text(size = 25, hjust = 0.5),
    legend.position = "bottom",
    legend.title  = element_text(face = "bold", size = 20),
    legend.text   = element_text(size = 25),
    legend.box    = "horizontal",
    legend.spacing.x = unit(0.8, "cm")   # Extra horizontal space between legend keys
  )



p <- ggarrange(plot_carbon_distrib,plot_subduct_distrib)
ggsave("figures/carbon_and_subd_distrib_combined.png",plot = p,width = 30,height = 10,dpi = 300)

ggsave("figures/carbon_distrib.png",plot = plot_carbon_distrib,width = 20,height = 10,dpi = 300)
ggsave("figures/carbon_distrib.pdf",plot = plot_carbon_distrib,width = 20,height = 10,dpi = 300)


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



# Convert data frames to spatial objects and transform to Robinson projection
df_argo_djf_sf <- df_argo_djf %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

df_argo_mam_sf <- df_argo_mam %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

df_argo_jja_sf <- df_argo_jja %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")

df_argo_son_sf <- df_argo_son %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE)) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = "+proj=robin")


# Function to plot Argo profiles for each season
plot_argo_profiles <- function(df_sf, title, output_file) {
  ggplot() +
    geom_sf(data = wld.rob.sf, fill = "grey", color = "gray") +
    geom_sf(data = df_sf, aes(geometry = geometry), color = "blue", alpha = 0.1, size = 0.1) +
    coord_sf(crs = "+proj=robin lon_0=180", datum = NA) +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(panel.grid.major = element_line(color = "gray90")) -> p
  
  # Save the plot
  ggsave(filename = output_file, plot = p, width = 10, height = 6)
}

# Generate maps for each season
plot_argo_profiles(df_argo_djf_sf, "DJF Argo Profiles", "figures/TimeSpaceVar/map_djf_argo_prof.jpg")
plot_argo_profiles(df_argo_mam_sf, "MAM Argo Profiles", "figures/TimeSpaceVar/map_mam_argo_prof.jpg")
plot_argo_profiles(df_argo_jja_sf, "JJA Argo Profiles", "figures/TimeSpaceVar/map_jja_argo_prof.jpg")
plot_argo_profiles(df_argo_son_sf, "SON Argo Profiles", "figures/TimeSpaceVar/map_son_argo_prof.jpg")
# Function for KDE analysis
kde_analysis <- function(df, n_grid) {
  MASS::kde2d(df$LONGITUDE, df$LATITUDE, n = n_grid)
}

# Run KDE for each season
kde_djf <- kde_analysis(df_argo_djf, 180)
kde_mam <- kde_analysis(df_argo_mam, 180)
kde_jja <- kde_analysis(df_argo_jja, 180)
kde_son <- kde_analysis(df_argo_son, 180)

# Convert KDE to data frame for plotting
kde_to_dataframe <- function(kde_result) {
  expand.grid(x = kde_result$x, y = kde_result$y) %>%
    mutate(z = as.vector(kde_result$z))
}

# Convert KDE results to data frames
kde_djf_df <- kde_to_dataframe(kde_djf)
kde_mam_df <- kde_to_dataframe(kde_mam)
kde_jja_df <- kde_to_dataframe(kde_jja)
kde_son_df <- kde_to_dataframe(kde_son)


# Function to plot KDE results
plot_kde <- function(kde_df, title, output_file) {
  ggplot() +
    geom_tile(data = kde_df, aes(x = x, y = y, fill = z)) +
    geom_contour(data = kde_df, aes(x = x, y = y, z = z), color = "white", alpha = 0.5) +
    coord_fixed() +
    geom_sf(data = world, fill = "gray80", color = "black") +
    scale_fill_viridis_c() +
    labs(title = title, x = "Longitude", y = "Latitude", fill = "Density") +
    theme_minimal() -> p
  
  # Save the plot
  ggsave(output_file, plot = p, width = 10, height = 6)
}

# Plot KDE maps for each season
plot_kde(kde_djf_df, "Density of DJF Argo Profiles", "figures/TimeSpaceVar/kde_djf_argo_profiles.png")
plot_kde(kde_mam_df, "Density of MAM Argo Profiles", "figures/TimeSpaceVar/kde_mam_argo_profiles.png")
plot_kde(kde_jja_df, "Density of JJA Argo Profiles", "figures/TimeSpaceVar/kde_jja_argo_profiles.png")
plot_kde(kde_son_df, "Density of SON Argo Profiles", "figures/TimeSpaceVar/kde_son_argo_profiles.png")

# Bin data
# Define bin size for longitude and latitude
bin_size <- 5
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


# Create proportion maps of SUBDUCTION (w/o carbon) without saving them to files directly

prop_map_full <- ggplot(merged_counts_full, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_djf$lon_bin), ylim = range(merged_counts_djf$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (whole year)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()


prop_map_djf <- ggplot(merged_counts_djf, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_djf$lon_bin), ylim = range(merged_counts_djf$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (DJF)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()

prop_map_mam <- ggplot(merged_counts_mam, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_mam$lon_bin), ylim = range(merged_counts_mam$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (MAM)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()

prop_map_jja <- ggplot(merged_counts_jja, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_jja$lon_bin), ylim = range(merged_counts_jja$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (JJA)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()

prop_map_son <- ggplot(merged_counts_son, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_son$lon_bin), ylim = range(merged_counts_son$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (SON)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()

# Combine proportion maps into a single figure using patchwork
combined_proportion_maps <- prop_map_djf + prop_map_mam + prop_map_jja + prop_map_son +
  plot_layout(ncol = 2, nrow = 2) +
  plot_annotation(title = "Proportion of Subduction Events Across Seasons")

# Save combined figure
ggsave(filename = "figures/TimeSpaceVar/4SEASONS/combined_proportion_maps.png", plot = combined_proportion_maps, width = 15, height = 12)

# Create proportion maps of CARBON SUBDUCTION without saving them to files directly
prop_map_full_carbon <- ggplot(merged_carbon_counts_full, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_djf$lon_bin), ylim = range(merged_counts_djf$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (whole year)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()




prop_map_djf_carbon <- ggplot(merged_carbon_counts_djf, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_djf$lon_bin), ylim = range(merged_counts_djf$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (DJF)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()

prop_map_mam_carbon <- ggplot(merged_carbon_counts_mam, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_mam$lon_bin), ylim = range(merged_counts_mam$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (MAM)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()

prop_map_jja_carbon <- ggplot(merged_carbon_counts_jja, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_jja$lon_bin), ylim = range(merged_counts_jja$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (JJA)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()

prop_map_son_carbon <- ggplot(merged_carbon_counts_son, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = proportion)) +
  geom_contour(aes(z = proportion), color = "white", alpha = 0.1) +
  geom_sf(data = world, fill = "gray80", color = "black", inherit.aes = FALSE) +
  coord_sf(xlim = range(merged_counts_son$lon_bin), ylim = range(merged_counts_son$lat_bin)) +
  scale_fill_viridis_c() +
  labs(title = "Proportion of Subduction Events (SON)", x = "Longitude", y = "Latitude", fill = "Proportion") +
  theme_minimal()

# Combine proportion maps into a single figure using patchwork
combined_proportion_maps_carbon <- prop_map_djf_carbon + prop_map_mam_carbon +
  prop_map_jja_carbon + prop_map_son_carbon +
  plot_layout(ncol = 2, nrow = 2) +
  plot_annotation(title = "Proportion of Carbon Subduction Events Across Seasons")

# Save combined figure
ggsave(filename = "figures/TimeSpaceVar/4SEASONS/combined_proportion_maps_carbon.png", plot = combined_proportion_maps_carbon, width = 15, height = 12)

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

create_prediction_grid <- function(merged_counts, step = 1) {
  lon_seq <- seq(min(merged_counts$lon_bin, na.rm = TRUE),
                 max(merged_counts$lon_bin, na.rm = TRUE))
                 
  
  lat_seq <- seq(min(merged_counts$lat_bin, na.rm = TRUE),
                 max(merged_counts$lat_bin, na.rm = TRUE),
                 by = step)
  expand.grid(lon_bin = lon_seq, lat_bin = lat_seq)
}

###############################################################################
# 3) Predict Smoothed Proportions
###############################################################################

predict_gam <- function(gam_model, merged_counts, step = 1) {
  pred_grid <- create_prediction_grid(merged_counts, step)
  pred_grid$proportion <- predict(gam_model, newdata = pred_grid, type = "response")
  pred_grid <- pred_grid %>% filter(!is.na(proportion))
  return(pred_grid)
}

###############################################################################
# 4) Plot a GAM Result (Lat-Lon) With Discrete Color Bins
###############################################################################
# We'll define the scale in a separate function so we can re-use it across plots.

make_discrete_scale <- function(prob_min, prob_max, binwidth = 0.1) {
  # E.g. if prob_max is 0.6 and binwidth is 0.1 => breaks = seq(0, 0.6, 0.1)
  breaks_vec <- seq(prob_min, prob_max, by = binwidth)
  
  scale_fill_viridis_b(
    name = "Probability",
    breaks = breaks_vec,
    limits = c(prob_min, prob_max),
    oob = scales::squish
  )
}
stipple_resolution <- 5  # degrees

plot_gam_map <- function(pred_grid, world_data, season_label, event_label, common_scale,argo_months) {
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
                 color = "white", alpha = 0.3,breaks = common_scale$breaks) +
    # Add stippling: plot four corner points per undersampled grid cell
    geom_point(data = undersampled_corners,
               aes(x = LON, y = LAT),
               color = "white", alpha = 1, size = 1, shape = 20) +
    # Add continents
    geom_sf(data = world_data, fill = "white", color = "white", inherit.aes = FALSE) +
    coord_sf(xlim = c(-180,180), ylim = c(-90,90), expand = FALSE, crs = st_crs(4326)) +
    labs(
      title = paste0("Estimated ", event_label, " Probability ", season_label, ""),
      x = "Longitude", y = "Latitude"
    ) +
    common_scale +    # apply the discrete color scale
    theme_minimal(base_size = 25) +  # Increased text size for readability
    theme(legend.title = element_blank(),legend.position = 'top', 
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

gam_full_subd <- fit_gam_season(merged_counts_full,  k_value = 600,argo_min = 25)
# High k variant, k = 600 (originally it was 300)
gam_djf_subd <- fit_gam_season(merged_counts_djf,  k_value = 600)
gam_mam_subd <- fit_gam_season(merged_counts_mam,  k_value = 600)
gam_jja_subd <- fit_gam_season(merged_counts_jja,  k_value = 600)
gam_son_subd <- fit_gam_season(merged_counts_son,  k_value = 600)

# Predict
pred_full_subd <- predict_gam(gam_full_subd, merged_counts_full, step = 1)
pred_full_subd_highres <- predict_gam(gam_full_subd, merged_counts_full, step = 0.25)


pred_djf_subd <- predict_gam(gam_djf_subd, merged_counts_djf, step = 1)
pred_mam_subd <- predict_gam(gam_mam_subd, merged_counts_mam, step = 1)
pred_jja_subd <- predict_gam(gam_jja_subd, merged_counts_jja, step = 1)
pred_son_subd <- predict_gam(gam_son_subd, merged_counts_son, step = 1)

# Determine global min/max for subduction proportions
subd_max_full <- max(pred_full_subd$proportion)
subd_full_scale <- make_discrete_scale(0,subd_max_full,binwidth = 0.05)

subd_values <- c(pred_djf_subd$proportion, pred_mam_subd$proportion,
                 pred_jja_subd$proportion, pred_son_subd$proportion)
subd_min <- 0  # typically 0
subd_max <- max(subd_values, na.rm = TRUE)

# Build a discrete color scale common to all subduction maps
# e.g. 0.1 step. Adjust if desired.
subd_scale <- make_discrete_scale(subd_min, subd_max, binwidth = 0.05)

# Plot each season with the common discrete scale

map_subd_full <- plot_gam_map(pred_full_subd, world, "Whole Year", "Subduction", subd_full_scale,argo_months = c(1:12))

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
gam_full_carb <- fit_gam_season(merged_carbon_counts_full,  k_value = 600,argo_min = 25)

gam_djf_carb <- fit_gam_season(merged_carbon_counts_djf,  k_value = 600,argo_min = 5)
gam_mam_carb <- fit_gam_season(merged_carbon_counts_mam,  k_value = 600,argo_min = 5)
gam_jja_carb <- fit_gam_season(merged_carbon_counts_jja,  k_value = 600,argo_min = 5)
gam_son_carb <- fit_gam_season(merged_carbon_counts_son,  k_value = 600,argo_min = 5)

# Predict
pred_full_carb <- predict_gam(gam_full_carb, merged_carbon_counts_full, step = 1)

pred_djf_carb <- predict_gam(gam_djf_carb, merged_carbon_counts_djf, step = 1)
pred_mam_carb <- predict_gam(gam_mam_carb, merged_carbon_counts_mam, step = 1)
pred_jja_carb <- predict_gam(gam_jja_carb, merged_carbon_counts_jja, step = 1)
pred_son_carb <- predict_gam(gam_son_carb, merged_carbon_counts_son, step = 1)

pred_djf_carb$season <- "DJF"
pred_mam_carb$season <- "MAM"
pred_jja_carb$season <- "JJA"
pred_son_carb$season <- "SON"

pred_son_carb_season_combined <- rbind(pred_djf_carb,
      pred_mam_carb,
      pred_jja_carb,
      pred_son_carb)

saveRDS(pred_son_carb_season_combined,"data/pred_grid_carb_prob.Rds")

# Global min/max for carbon subduction
carb_values <- c(pred_djf_carb$proportion, pred_mam_carb$proportion,
                 pred_jja_carb$proportion, pred_son_carb$proportion)
carb_min <- 0
carb_max <- max(carb_values, na.rm = TRUE)

carb_max_full <- max(pred_full_carb$proportion)
carb_full_scale <- make_discrete_scale(0,0.25,binwidth = 0.05)

# Discrete color scale for carbon subduction
carb_scale <- make_discrete_scale(carb_min, carb_max, binwidth = 0.05)

map_carb_full <- plot_gam_map(pred_full_carb,  world, "", "Carbon Subduction", carb_full_scale,argo_months = c(1:12))
map_carb_full <- ggarrange(map_carb_full,
                           common.legend = T,legend="bottom")

map_subd_full <- plot_gam_map(pred_full_subd,  world, "", "Subduction", subd_scale,argo_months = c(1:12))
map_subd_full <- ggarrange(map_subd_full,
                            common.legend = T,legend="bottom")
map_carb_full_subd_scale <- plot_gam_map(pred_full_carb,  world, "", "Carbon Subduction", subd_scale,argo_months = c(1:12))


ggarrange(plot_carbon_distrib,plot_subduct_distrib,map_carb_full,map_subd_full,nrow = 2,ncol = 2)

map_combined <- ggarrange(map_carb_full,map_subd_full,legend="bottom")

ggsave(plot = map_combined,
       filename = "figures/TimeSpaceVar/4SEASONS/gam_combined_discrete_yearly.png",
       width = 23,height = 8)

library(cowplot)



ggsave("figures/combined_4maps_cowplot.png", final_cow, width=20, height=15, dpi=300)



map_carb_subd_full <- ggarrange(map_subd_full,map_carb_full_subd_scale,common.legend = TRUE)

map_carb_djf <- plot_gam_map(pred_djf_carb, world, "DJF", "Carbon Subduction", carb_scale,c(12,1,2))
map_carb_mam <- plot_gam_map(pred_mam_carb, world, "MAM", "Carbon Subduction", carb_scale,c(3:5))
map_carb_jja <- plot_gam_map(pred_jja_carb, world, "JJA", "Carbon Subduction", carb_scale,c(6:8))
map_carb_son <- plot_gam_map(pred_son_carb, world, "SON", "Carbon Subduction", carb_scale,c(9:11))

combined_carb <- ggarrange(map_carb_djf,map_carb_mam,map_carb_jja,map_carb_son,ncol = 2,nrow=2,
                           common.legend = T,legend="bottom")

ggsave("figures/TimeSpaceVar/4SEASONS/gam_carbon_subduction_discrete_yearly.png",
       map_carb_full, width = 18, height = 10)


ggsave("figures/TimeSpaceVar/4SEASONS/gam_carbon_subduction_discrete_4seasons.png",
       combined_carb, width = 23, height = 15.3)


# Summary:
# - We create 8 total plots: 4 for subduction, 4 for carbon subduction,
#   each with a discrete color scale and a uniform range per dataset.
# - The scale is set by the global min (0) and max proportion among
#   the four seasons, so they're all directly comparable.
# - Adjust 'binwidth' in 'make_discrete_scale()' if I want finer or coarser bins.

###############################################################################
# Done!
###############################################################################

# Now let's do histograms faceted by region showing the probability of (carbon) subduction
# by months :

###############################################################################
# Seasonal Cycle of Subduction (Carbon and Non-Carbon) by Region and Month
###############################################################################

# Add Month and Region Information to Data
# Splitting Northern Tropics and Southern Tropics as a sanity check (should not be seasonal there)
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

# Not Splitting Northern Tropics and Southern Tropics as a sanity check (should not be seasonal there)
# add_month_region <- function(df) {
#   df %>%
#     mutate(
#       month = month(TIME, label = TRUE, abbr = TRUE),  # Extract month as a factor
#       region = case_when(
#         LATITUDE > 30 & LONGITUDE >= -100 & LONGITUDE <= 20 ~ "North Atlantic",  # Above 30°N
#         LATITUDE > 30 & (LONGITUDE < -100 | LONGITUDE > 120) ~ "North Pacific",  # Above 30°N
#         LATITUDE <= 30 & LATITUDE >= -30 ~ "Tropics",                     # 0° to 30°S
#         LATITUDE < -30 ~ "Southern Ocean",                                       # Below 30°S
#         TRUE ~ NA_character_  # Exclude undefined regions
#       )
#     )
# }
# 
# 
df_argo_clean <- add_month_region(df_argo_clean)
df_complete_clean <- add_month_region(df_complete_clean)
df_carbon_clean <- add_month_region(df_carbon_clean)
# 
# # Step 2: Compute Probabilities by Region and Month
compute_monthly_probabilities <- function(df_argo, df_events) {
   # Total Argo profiles by region and month
   total_counts <- df_argo %>%
     group_by(region, month) %>%
     summarize(count_total = n(), .groups = "drop")
   
   # Subduction events by region and month
   event_counts <- df_events %>%
     group_by(region, month) %>%
     summarize(count_event = n(), .groups = "drop")
   
   # Merge and compute proportions
   merged <- full_join(total_counts, event_counts, by = c("region", "month")) %>%
     mutate(
       count_event = replace_na(count_event, 0),
       proportion = ifelse(count_total > 0, count_event / count_total, NA)
     )
   return(merged)
 }

# Compute probabilities for subduction and carbon subduction
monthly_probs_subduction <- compute_monthly_probabilities(df_argo_clean, df_complete_clean)
monthly_probs_carbon <- compute_monthly_probabilities(df_argo_clean, df_carbon_clean)

monthly_probs_subduction$carbon <- "Subduction with/without carbon"
monthly_probs_carbon$carbon <- "Subduction with carbon"

monthly_probs_binded <- bind_rows(monthly_probs_carbon,monthly_probs_subduction) %>%
  filter(!is.na(region), !is.na(proportion), !is.na(month))  # Remove NA regions, proportions, and months %>%

ggplot(monthly_probs_binded, aes(x = month, y = proportion, fill = region)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.8) +
  facet_wrap(region~carbon, ncol = 2,axes="all",scales = "free") +
  scale_fill_viridis_d() +
  labs(
    title = title,
    x = "Month",
    y = "Probability of Subduction",
    fill = "Region"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15,angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 16),
    title = element_text(size=15)
  )+  
  guides(fill = "none")

# Step 3: Plot Monthly Histograms Faceted by Region
plot_monthly_histograms <- function(monthly_probs, title) {
  monthly_probs <- monthly_probs %>%
    filter(!is.na(region), !is.na(proportion), !is.na(month))  # Remove NA regions, proportions, and months
  
  ggplot(monthly_probs, aes(x = month, y = proportion, fill = region)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.8) +
    facet_wrap(~region, ncol = 2,axes="all") +
    scale_fill_viridis_d() +
    labs(
      title = title,
      x = "Month",
      y = "Probability of Subduction",
      fill = "Region"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 15,angle = 45, hjust = 1),
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 16),
      title = element_text(size=15)
    )+  
    guides(fill = "none")
}


# Plot subduction probabilities
plot_subduction <- plot_monthly_histograms(monthly_probs_subduction, "Probability that an Argo Float captures a subduction event")

# Plot carbon subduction probabilities
plot_carbon_subduction <- plot_monthly_histograms(monthly_probs_carbon, "Probability that an Argo Float captures a subduction event with carbon")

# Save the plots
ggsave("figures/seasonal_subduction_by_region.png", plot = plot_subduction, width = 10, height = 8)
ggsave("figures/seasonal_carbon_subduction_by_region.png", plot = plot_carbon_subduction, width = 10, height = 8)

### Add a statistical test to bring evidence against H0 that the monthly frequency of subduction and of carbon subduction
# sampled from a uniform distribution to test whether subduction is or not seasonal
# Remove NA months from both datasets
monthly_probs_subduction <- monthly_probs_subduction %>% filter(!is.na(month))
monthly_probs_carbon <- monthly_probs_carbon %>% filter(!is.na(month))

# Function to test uniformity for each region
test_seasonality_proportions <- function(monthly_probs) {
  monthly_probs %>%
    group_by(region) %>%
    summarize(
      chisq_stat = chisq.test(
        x = proportion * count_total,  # Observed weighted counts
        p = rep(1 / 12, 12),           # Expected uniform proportions
        rescale.p = TRUE               # Rescale expected to match total
      )$statistic,
      p_value = chisq.test(
        x = proportion * count_total,
        p = rep(1 / 12, 12),
        rescale.p = TRUE
      )$p.value,
      .groups = "drop"
    )
}


t <- monthly_probs_carbon %>% filter(region == "North Atlantic")
chisq.test(t$count_event)

# Test seasonality for subduction
seasonality_test_subduction <- test_seasonality_proportions(monthly_probs_subduction)

# Test for carbon subduction
seasonality_test_carbon <- test_seasonality_proportions(monthly_probs_carbon)

# Display results
print("Subduction Seasonality Test:")
print(seasonality_test_subduction)

print("Carbon Subduction Seasonality Test:")
print(seasonality_test_carbon)




# Run the post hoc test on the North Atlantic data
posthoc_results_carbon <- post_hoc_test_region(north_atlantic_carbon)
print("Post hoc test results for Carbon Subduction (North Atlantic):")
print(posthoc_results_carbon)

#############################################
#### SEASONALITY OF DEPTH ###################
#############################################

# 1. Summarize depth (pressure) by region and month
df_depth_summary <- df_complete_clean %>%
  filter(!is.na(region), !is.na(month)) %>%
  group_by(region, month) %>%
  summarize(
    n                = n(),
    median_pressure  = median(PRES_ADJUSTED, na.rm = TRUE),
    iqr_pressure     = IQR(PRES_ADJUSTED, na.rm = TRUE),
    .groups = "drop"
  )

df_carbon_depth_summary <- df_carbon_clean %>%
  filter(!is.na(region), !is.na(month)) %>%
  group_by(region, month) %>%
  summarize(
    n                = n(),
    median_pressure  = median(PRES_ADJUSTED, na.rm = TRUE),
    iqr_pressure     = IQR(PRES_ADJUSTED, na.rm = TRUE),
    .groups = "drop"
  )

# 1. Create a summary dataset of median depth per region and month
df_median_line <- df_complete_clean %>%
  filter(!is.na(region), !is.na(month)) %>%
  group_by(region, month) %>%
  summarize(
    median_pressure = median(PRES_ADJUSTED, na.rm = TRUE),
    .groups = "drop"
  )

df_median_carbon_line <- df_carbon_clean %>%
  filter(!is.na(region), !is.na(month)) %>%
  group_by(region, month) %>%
  summarize(
    median_pressure = median(PRES_ADJUSTED, na.rm = TRUE),
    .groups = "drop"
  )


# Print the summary table
print(df_depth_summary)
print(df_carbon_depth_summary)

# 2. Plot a boxplot of subduction depth by region and month
plot_subduction_depth <- ggplot(df_complete_clean %>% 
                                  filter(!is.na(region), !is.na(month)),
                                aes(x = month, 
                                    y = PRES_ADJUSTED, 
                                    fill = region)) +
  geom_boxplot() + # Median line (and optional points) for each region facet
  geom_line(
    data = df_median_line,
    aes(
      x = month, 
      y = median_pressure, 
      group = region
    ),
    color = "black", 
    size = 1
  ) +
  geom_point(
    data = df_median_line,
    aes(
      x = month, 
      y = median_pressure
    ),
    color = "black", 
    size = 2
  )+
  facet_wrap(~ region, ncol = 2,axes="all") +
  scale_fill_viridis_d() +
  labs(
    title = "Monthly Distribution of Subduction Depth by Region",
    x = "Month",
    y = "Pressure (dbar)",
    fill = "Region"
  ) +
  theme_minimal() +
  theme(
    strip.text     = element_text(size = 15, face = "bold"),
    axis.text.x    = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y    = element_text(size = 15),
    axis.title.y   = element_text(size = 16),
    title          = element_text(size = 15)
  ) +
  guides(fill = "none")+scale_y_reverse()

plot_carbon_subduction_depth <- ggplot(df_carbon_clean %>% 
                                  filter(!is.na(region), !is.na(month)),
                                aes(x = month, 
                                    y = PRES_ADJUSTED, 
                                    fill = region)) +
  geom_boxplot() +# Median line (and optional points) for each region facet
  geom_line(
    data = df_median_carbon_line,
    aes(
      x = month, 
      y = median_pressure, 
      group = region
    ),
    color = "black", 
    size = 1
  ) +
  geom_point(
    data = df_median_carbon_line,
    aes(
      x = month, 
      y = median_pressure
    ),
    color = "black", 
    size = 2
  )+
  facet_wrap(~ region, ncol = 2,axes="all") +
  scale_fill_viridis_d() +
  labs(
    title = "Monthly Distribution of Carbon Subduction Depth by Region",
    x = "Month",
    y = "Pressure (dbar)",
    fill = "Region"
  ) +
  theme_minimal() +
  theme(
    strip.text     = element_text(size = 15, face = "bold"),
    axis.text.x    = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y    = element_text(size = 15),
    axis.title.y   = element_text(size = 16),
    title          = element_text(size = 15)
  ) +
  guides(fill = "none")+scale_y_reverse()



# Display the plot
print(plot_subduction_depth)
plot_carbon_subduction_depth


df_carbon_clean$PRES_ADJUSTED


#  save the figure
ggsave("figures/seasonal_subduction_depth_by_region.png",
       plot   = plot_subduction_depth,
       width  = 10,
       height = 8)

ggsave("figures/seasonal_carbon_subduction_depth_by_region.png",
       plot   = plot_carbon_subduction_depth,
       width  = 10,
       height = 8)


# 3. Statistical test for monthly differences in subduction depth
#    Here we use a Kruskal-Wallis test for each region to see if 
#    the distribution of PRES_ADJUSTED differs by month
seasonality_test_depth <- df_complete_clean %>%
  filter(!is.na(region), !is.na(month)) %>%
  group_by(region) %>%
  do({
    kruskal_result <- kruskal.test(PRES_ADJUSTED ~ month, data = .)
    data.frame(
      kruskal_stat = kruskal_result$statistic,
      df           = kruskal_result$parameter,
      p_value      = kruskal_result$p.value
    )
  }) %>%
  ungroup()

seasonality_carbon_test_depth <- df_carbon_clean %>%
  filter(!is.na(region), !is.na(month)) %>%
  group_by(region) %>%
  do({
    kruskal_result <- kruskal.test(PRES_ADJUSTED ~ month, data = .)
    data.frame(
      kruskal_stat = kruskal_result$statistic,
      df           = kruskal_result$parameter,
      p_value      = kruskal_result$p.value
    )
  }) %>%
  ungroup()


# Print the test results
print("Seasonality Test for Subduction Depth (Kruskal-Wallis):")
print(seasonality_test_depth)
print(seasonality_carbon_test_depth)
anova_model <- aov(PRES_ADJUSTED ~ month, data = df_carbon_clean %>% filter(region == "Southern Ocean"))
summary(anova_model)

m <- lm(PRES_ADJUSTED ~ factor(month,ordered = FALSE), data = df_carbon_clean %>% filter(region == "Southern Ocean"))
m <- lm(PRES_ADJUSTED ~ region, data = df_carbon_clean)

m %>% summary()

anova_model <- aov(PRES_ADJUSTED ~ month , data = df_carbon_clean %>% filter(region == "North Atlantic"))

# Perform Tukey's HSD post hoc test to see which months differ significantly
tukey_results <- TukeyHSD(anova_model, conf.level = 0.95)
print(tukey_results)

# Seasonality of depth of subd Southern Ocean, seasonality, for southern tropics as well. 
# Seasonality of carbon depth of subd Southern Ocean, seasonality, for southern tropics as well.
##############################################################################
# 3) Fit the Exponential and Gamma distributions with fitdistrplus
##############################################################################

pres_data <- df_complete_clean %>% filter(!is.na(region), !is.na(month))
pres_data$PRES_ADJUSTED <- pres_data$PRES_ADJUSTED-200


carbon_pres_data <- df_carbon_clean %>% filter(!is.na(region), !is.na(month))
carbon_pres_data$PRES_ADJUSTED <- carbon_pres_data$PRES_ADJUSTED-200

##############################################################################
# 3) Define a helper function to fit Exp, Gamma, and Power Law, and 
#    return a data frame of x vs. PDF for each distribution.
##############################################################################

# Southern Ocean, seasonality, for southern tropics as well.

##############################################################################
# 3) Fit the Exponential and Gamma distributions with fitdistrplus
##############################################################################

pres_data <- df_complete_clean %>% filter(!is.na(region), !is.na(month))
pres_data$PRES_ADJUSTED <- pres_data$PRES_ADJUSTED-200


carbon_pres_data <- df_carbon_clean %>% filter(!is.na(region), !is.na(month))
carbon_pres_data$PRES_ADJUSTED <- carbon_pres_data$PRES_ADJUSTED-200
##############################################################################
# 3) Define a helper function to fit Exp, Gamma, LogNormal, Cauchy, and Power Law,
#    and return a data frame of x vs. PDF for each distribution.
##############################################################################
data.x <- carbon_pres_data %>% filter(region == "Southern Ocean") %>% pull(PRES_ADJUSTED)

m_pl = displ$new(data.x)
m_ln = dislnorm$new(data.x)

get_fitted_pdf_data <- function(sub_df, n_points = 500) {
  # sub_df is the subset for one region
  
  # 3a) Extract the numeric vector
  x_vals <- sub_df$PRES_ADJUSTED
  
  # 3b) Fit Exponential and Gamma via fitdistrplus
  fit_exp   <- fitdist(x_vals, "exp")
  fit_gamma <- fitdist(x_vals, "gamma")
  
  # 3c) Fit Log Normal via fitdistrplus
  fit_lognormal <- fitdist(x_vals, "lnorm")
  
  # 3d) Fit Cauchy via fitdistrplus
  fit_cauchy <- fitdist(x_vals, "cauchy",
                        start = list(location = mean(x_vals), scale = stats::sd(x_vals)))
  
  # 3e) Create a sequence of x-values spanning the data range
  x_seq <- seq(
    from = max(0, min(x_vals)),  # assume 0 as lower bound if needed
    to   = max(x_vals),
    length.out = n_points
  )
  
  # 3f) Compute PDFs from each fitted model
  exp_density       <- dexp(x_seq, rate  = fit_exp$estimate["rate"])
  gamma_density     <- dgamma(x_seq, shape = fit_gamma$estimate["shape"],
                              rate  = fit_gamma$estimate["rate"])
  lognormal_density <- dlnorm(x_seq, meanlog = fit_lognormal$estimate["meanlog"],
                              sdlog = fit_lognormal$estimate["sdlog"])
  cauchy_density    <- dcauchy(x_seq, location = fit_cauchy$estimate["location"],
                               scale    = fit_cauchy$estimate["scale"])
  
  # 3g) Create a long-format data frame of x vs. PDF for each distribution
  density_df <- data.frame(
    x = x_seq,
    Exp       = exp_density,
    Gamma     = gamma_density,
    LogNormal = lognormal_density,
    Cauchy    = cauchy_density
  ) %>%
    pivot_longer(
      cols = c("Exp", "Gamma", "LogNormal", "Cauchy"),
      names_to = "Distribution",
      values_to = "Density"
    )
  
  # 3h) Compute GOF measures using gofstat() for each fitted model
  gof_exp       <- gofstat(fit_exp)
  gof_gamma     <- gofstat(fit_gamma)
  gof_lognormal <- gofstat(fit_lognormal)
  gof_cauchy    <- gofstat(fit_cauchy)
  
  gof_df <- data.frame(
    Distribution = c("Exp", "Gamma", "LogNormal", "Cauchy"),
    AIC = c(gof_exp$aic, gof_gamma$aic, gof_lognormal$aic, gof_cauchy$aic),
    BIC = c(gof_exp$bic, gof_gamma$bic, gof_lognormal$bic, gof_cauchy$bic),
    KS  = c(gof_exp$ks, gof_gamma$ks, gof_lognormal$ks, gof_cauchy$ks),
    CvM = c(gof_exp$cvm, gof_gamma$cvm, gof_lognormal$cvm, gof_cauchy$cvm),
    AD  = c(gof_exp$ad, gof_gamma$ad, gof_lognormal$ad, gof_cauchy$ad)
  )
  
  # Return both the density data and GOF measures as a list
  return(list(density_df = density_df, gof_df = gof_df))
}

# Here each region returns a list with two elements: density_df and gof_df
fitted_list <- pres_data %>%
  group_by(region) %>%
  do(fit_out = get_fitted_pdf_data(.)[[1]])

gof_list <- pres_data %>%
  group_by(region) %>%
  do(fit_out = get_fitted_pdf_data(.)[[2]])

gof_list$fit_out[[1]]
gof_list$fit_out[[2]]
gof_list$fit_out[[3]]
gof_list$fit_out[[4]]



##############################################################################
# 4) Apply the fitting function by region, producing a data frame 
#    that contains x, Distribution, Density, region, etc.
##############################################################################
fitted_df <- pres_data %>%
  group_by(region) %>%
  do( get_fitted_pdf_data(.)[[1]] )  # for each region, produce the PDF data

# Make sure to keep region as a column (rather than a grouping var):
fitted_df <- ungroup(fitted_df)

##############################################################################
# 5) Plot faceted by region: histogram (or density) + fitted curves
##############################################################################
# We'll illustrate with a histogram, scaled to density on the y-axis.
# If you prefer a kernel density plot for the raw data, replace geom_histogram 
# with geom_density, but note that `bins = 30` applies to histograms, not densities.

# Comparing statistical distributions
x <- pres_data %>% filter(region == "Southern Ocean") %>% pull(PRES_ADJUSTED)



distrib_subduction <- ggplot() +
  # 5a) Histogram of the original data, faceted by region
  geom_histogram(
    data = pres_data,
    aes(x = PRES_ADJUSTED, y = ..density.., fill = region),
    bins  = 30,
    color = "black",
    fill  = "white",
    alpha = 0.8
  ) +
  # 5b) Overlaid lines from the fitted distributions
  geom_line(
    data = fitted_df ,
    aes(x = x, y = Density, color = Distribution),
    size = 1
  ) +
  # 5c) Facet by region
  facet_wrap(~ region, scales = "free_y") +
  # 5d) Color scale and labeling
  scale_color_manual(values = c("Exp" = "red", 
                                "Gamma" = "blue", 
                                "LogNormal" = "orange", 
                                "Cauchy" = "purple", 
                                "PL" = "green")) +
  labs(
    title = "Fitted Distributions by Region",
    x     = "Pressure Adjusted",
    y     = "Density"
  ) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  # 5e) Customizing x-axis ticks
  scale_x_continuous(
    breaks = seq(0, 800, by = 200),
    labels = seq(200, 1000, by = 200)
  )



fitted_df <- carbon_pres_data %>%
  group_by(region) %>%
  do( get_fitted_pdf_data(.)[[1]] )  # for each region, produce the PDF data

# Make sure to keep region as a column (rather than a grouping var):
fitted_df <- ungroup(fitted_df)


GOF_df <- carbon_pres_data %>%
  group_by(region) %>%
  do( get_fitted_pdf_data(.)[[2]] )

GOF_df
##############################################################################
# 6) Plot faceted by region: histogram (or density) + fitted curves
##############################################################################
# We'll illustrate with a histogram, scaled to density on the y-axis.
# If you prefer a kernel density plot for the raw data, replace geom_histogram 
# with geom_density, but note that `bins = 30` applies to histograms, not densities.
carbon_pres_data$region

distrib_carbon_subduction <- ggplot() +
  # 5a) Histogram of the original data, faceted by region
  geom_histogram(
    data = carbon_pres_data,
    aes(x = PRES_ADJUSTED, fill = region), #y=..density.. for density
    bins  = 30,
    color = "black",
    alpha = 0.8
  ) + scale_fill_viridis_d()+
  facet_wrap(~ region, scales = "free_y",ncol = 2,axes = "all_x") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  # 5e) Customizing x-axis ticks
  scale_x_continuous(
    breaks = seq(0, 800, by = 200),
    labels = seq(200, 1000, by = 200)
  )+ 
  theme_minimal() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15,angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 16),
    title = element_text(size=15)
  )+  
  guides(fill = "none")




ggsave("figures/distrib_carbon_subd.png",distrib_carbon_subduction,width = 10,height = 8)


distrib_subd_histogram <- ggplot() +
  # 5a) Histogram of the original data, faceted by region
  geom_histogram(
    data = pres_data,
    aes(x = PRES_ADJUSTED, fill = region),
    bins  = 30,
    color = "black",
    alpha = 0.8
  ) + scale_fill_viridis_d()+
  facet_wrap(~ region, scales = "free_y",ncol = 2,axes = "all_x") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  # 5e) Customizing x-axis ticks
  scale_x_continuous(
    breaks = seq(0, 800, by = 200),
    labels = seq(200, 1000, by = 200)
  )+ 
  theme_minimal() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15,angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 16),
    title = element_text(size=15)
  )+  
  guides(fill = "none")


pres_data$Type <- "Subduction"
carbon_pres_data$Type <- "Carbon Subduction"
# Combine data and create a new variable to distinguish datasets
combined_data <- bind_rows(
  pres_data %>% mutate(Type = "Subduction Events"),
  carbon_pres_data %>% mutate(Type = "Carbon Subduction Events")
)

# Plot overlaid density plots
subd_distrib_plot <- ggplot(combined_data, aes(x = PRES_ADJUSTED, fill = Type,color = Type)) +
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
    legend.text = element_text(size = 15)
  ) +
  labs(
    x = "Depth (m)",
    y = "Estimated Density",
    title = "Carbon subduction happens mostly between 200 and 500 meters"
  )+  scale_x_continuous(
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




ggsave(filename = "figures/subd_distrib_plot.png",plot = subd_distrib_plot,width = 10,height = 10)

# Compute proportion of events with PRES_ADJUSTED >= 300
proportions_df <- combined_data %>%
  group_by(Type, region) %>%
  summarize(
    total_events = n(),
    events_geq_300 = sum(PRES_ADJUSTED >= 300),
    proportion = events_geq_300 / total_events
  ) %>%
  arrange(region, Type)

# Display the results
print(proportions_df %>% filter(Type == "Carbon Subduction Events"))




ggsave("figures/distrib_subd_histogram.png",distrib_subd_histogram,width = 10,height = 8)


#  save the figure
ggsave("figures/seasonal_subduction_depth_distribution.png",
       plot   = distrib_subduction,
       width  = 10,
       height = 8)

ggsave("figures/seasonal_carbon_subduction_depth_distribution.png",
       plot   = distrib_carbon_subduction,
       width  = 10,
       height = 8)


# Plot ECDFs :

carbon_pres_data %>%
  ggplot(aes(x = PRES_ADJUSTED, color = region)) +
  stat_ecdf(size = 1.2) +
  labs(
    title = "Empirical CDF of PRES_ADJUSTED by Region",
    x = "PRES_ADJUSTED",
    y = "Empirical Cumulative Probability",
    color = "Region"
  ) +
  theme_minimal(base_size = 16) +  +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15,angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 16),
    title = element_text(size=15)
  )+  
  guides(fill = "none")




