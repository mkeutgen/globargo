# Supplementary Figure: Resolution Comparison
# This script creates a 6-panel figure comparing proportion maps at different resolutions
# Top row: 5 degree resolution (left: physical subduction, right: carbon subduction)
# Middle row: 1 degree resolution (left: physical subduction, right: carbon subduction)
# Bottom row: GAM 0.25 degree resolution (left: physical subduction, right: carbon subduction)

# Load necessary libraries
library(tidyverse)
library(rnaturalearth)
library(sf)
library(viridis)
library(patchwork)
library(lubridate)

# Load world map
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# theme_map function (from script_figures_generation.R)
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

# Load data
df_argo_clean <- read_csv("../data/df_argo_loc.csv")
df_complete_clean <- read_csv("../data/manually_verified_physical_subd_events.csv")
df_carbon_clean <- read_csv("../data/df_carbon_subduction_anom.csv")

# Load GAM predictions at 0.25 degree resolution
pred_full_subd <- readRDS("../data/pred_full_subd025.Rds")
pred_full_carb <- readRDS("../data/pred_full_carb025.Rds")

# Function to assign bins with variable bin size
# Use fixed global bins to ensure consistency across all datasets
assign_bins <- function(df, bin_size) {
  # Use global coverage to ensure all data uses same bins
  longitude_bins <- seq(-180, 180, by = bin_size)
  latitude_bins <- seq(-180, 180, by = bin_size)  # Changed from -90, 90 to match latitude range

  lon_centers <- longitude_bins[-length(longitude_bins)] + bin_size / 2
  lat_centers <- latitude_bins[-length(latitude_bins)] + bin_size / 2

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

# Function to compute counts and proportions
compute_counts <- function(df_argo, df_complete) {
  total_counts <- df_argo %>%
    group_by(lon_bin, lat_bin) %>%
    summarize(count_total = n(), .groups = 'drop')

  anomaly_counts <- df_complete %>%
    group_by(lon_bin, lat_bin) %>%
    summarize(count_anomaly = n(), .groups = 'drop')

  merged_counts <- full_join(total_counts, anomaly_counts, by = c("lon_bin", "lat_bin")) %>%
    mutate(
      count_anomaly = ifelse(is.na(count_anomaly), 0, count_anomaly),
      proportion = count_anomaly / count_total
    ) %>%
    filter(count_total > 0) %>%
    filter(proportion < 1) %>%  # Filter out proportions >= 1
    filter(!is.na(lat_bin))

  return(merged_counts)
}

# Create proportion maps at different resolutions
# 5 degree resolution
df_argo_5deg <- assign_bins(df_argo_clean, bin_size = 5)
df_complete_5deg <- assign_bins(df_complete_clean, bin_size = 5)
df_carbon_5deg <- assign_bins(df_carbon_clean, bin_size = 5)

merged_counts_5deg_subd <- compute_counts(df_argo_5deg, df_complete_5deg)
merged_counts_5deg_carb <- compute_counts(df_argo_5deg, df_carbon_5deg)


merged_counts_5deg_carb %>% View()
# 1 degree resolution
df_argo_1deg <- assign_bins(df_argo_clean, bin_size = 1)
df_complete_1deg <- assign_bins(df_complete_clean, bin_size = 1)
df_carbon_1deg <- assign_bins(df_carbon_clean, bin_size = 1)

merged_counts_1deg_subd <- compute_counts(df_argo_1deg, df_complete_1deg)
merged_counts_1deg_carb <- compute_counts(df_argo_1deg, df_carbon_1deg)

# Create plotting function for proportion maps
create_prop_map <- function(merged_counts, title, color_breaks, color_labels, color_limits) {
  ggplot(merged_counts, aes(x = lon_bin, y = lat_bin)) +
    geom_tile(aes(fill = proportion)) +
    geom_contour(aes(z = proportion), color = "white", alpha = 0.3, bins = 5) +
    geom_sf(data = world, fill = "lightgrey", color = "lightgrey", inherit.aes = FALSE) +
    coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE, crs = st_crs(4326)) +
    scale_fill_viridis_b(
      name = "",
      breaks = color_breaks,
      labels = color_labels,
      limits = color_limits,
      oob = scales::squish,
      na.value = "white"
    ) +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_map() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.width = unit(3, "cm")
    )
}

# Create the four panels
# Panel A: 5 degree - Physical Subduction
panel_5deg_subd <- create_prop_map(
  merged_counts_5deg_subd,
  "a • 5° Resolution - Physical Subduction",
  color_breaks = c(0.001, 0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.60),
  color_labels = c("0.1%", "1%", "5%", "10%", "20%", "30%", "40%", "60%"),
  color_limits = c(0.001, 0.6)
)

# Panel B: 5 degree - Carbon Subduction
panel_5deg_carb <- create_prop_map(
  merged_counts_5deg_carb,
  "b • 5° Resolution - Carbon Subduction",
  color_breaks = c(0.0005, 0.005, 0.025, 0.05, 0.10, 0.15, 0.20, 0.30),
  color_labels = c("0.05","0.5","2.5","5","10","15","20","30 (%)"),
  color_limits = c(0.0005, 0.30)
)

# Panel C: 1 degree - Physical Subduction
panel_1deg_subd <- create_prop_map(
  merged_counts_1deg_subd,
  "c • 1° Resolution - Physical Subduction",
  color_breaks = c(0.001, 0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.60),
  color_labels = c("0.1","1","5","10","20","30","40","60 (%)"),
  color_limits = c(0.001, 0.6)
)

# Panel D: 1 degree - Carbon Subduction
panel_1deg_carb <- create_prop_map(
  merged_counts_1deg_carb,
  "d • 1° Resolution - Carbon Subduction",
  color_breaks = c(0.0005, 0.005, 0.025, 0.05, 0.10, 0.15, 0.20, 0.30),
  color_labels = c("0.05","0.5","2.5","5","10","15","20","30 (%)"),
  color_limits = c(0.0005, 0.30)
)

# Panel E: GAM 0.25 degree - Physical Subduction

# Panel F: GAM 0.25 degree - Carbon Subduction
map_subd_full <- plot_gam_map( pred_full_subd, world,"e • Physical Subduction GAM Model", 1:12) +
  theme_map() + guides(colour = "none")+theme(legend.spacing.x = unit(0.8, "cm"))  +
  scale_y_continuous(breaks = ticks_y, labels = label_lat) 

map_carb_full <- plot_carb_gam_map(pred_full_carb, world, "f • Carbon Subduction GAM Model", 1:12) +
  theme_map() + guides(colour = "none")+theme(legend.spacing.x = unit(0.8, "cm")) +
  scale_y_continuous(breaks = ticks_y, labels = label_lat) 

# Combine all panels using patchwork
# Top row: 5 degree resolution
# Middle row: 1 degree resolution
# Bottom row: GAM 0.25 degree resolution
# Left column: Physical subduction
# Right column: Carbon subduction
combined_figure <- (panel_5deg_subd | panel_5deg_carb) /
                   (panel_1deg_subd | panel_1deg_carb) /
                   (map_subd_full | map_carb_full) +
  plot_annotation(
    title = "Comparison of Subduction Event Proportions Gridded at Different Grid Resolutions\n with the GAM Interpolation",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

# Save the figure
ggsave(
  filename = "../pubfig/supplementary_resolution_comparison.png",
  plot = combined_figure,
  width = 12,
  height = 16,  # Increased height for 3 rows
  dpi = 300
)


# Print summary statistics for the manuscript response
cat("\n=== Summary Statistics for Reviewer Response ===\n")
cat("\nPhysical Subduction - 5 degree resolution:\n")
cat("  Number of grid cells:", nrow(merged_counts_5deg_subd), "\n")
cat("  Mean proportion:", round(mean(merged_counts_5deg_subd$proportion, na.rm = TRUE), 4), "\n")
cat("  Median proportion:", round(median(merged_counts_5deg_subd$proportion, na.rm = TRUE), 4), "\n")
cat("  Range:", round(min(merged_counts_5deg_subd$proportion, na.rm = TRUE), 4), "-",
    round(max(merged_counts_5deg_subd$proportion, na.rm = TRUE), 4), "\n")

cat("\nCarbon Subduction - 5 degree resolution:\n")
cat("  Number of grid cells:", nrow(merged_counts_5deg_carb), "\n")
cat("  Mean proportion:", round(mean(merged_counts_5deg_carb$proportion, na.rm = TRUE), 4), "\n")
cat("  Median proportion:", round(median(merged_counts_5deg_carb$proportion, na.rm = TRUE), 4), "\n")
cat("  Range:", round(min(merged_counts_5deg_carb$proportion, na.rm = TRUE), 4), "-",
    round(max(merged_counts_5deg_carb$proportion, na.rm = TRUE), 4), "\n")

cat("\nPhysical Subduction - 1 degree resolution:\n")
cat("  Number of grid cells:", nrow(merged_counts_1deg_subd), "\n")
cat("  Mean proportion:", round(mean(merged_counts_1deg_subd$proportion, na.rm = TRUE), 4), "\n")
cat("  Median proportion:", round(median(merged_counts_1deg_subd$proportion, na.rm = TRUE), 4), "\n")
cat("  Range:", round(min(merged_counts_1deg_subd$proportion, na.rm = TRUE), 4), "-",
    round(max(merged_counts_1deg_subd$proportion, na.rm = TRUE), 4), "\n")

cat("\nCarbon Subduction - 1 degree resolution:\n")
cat("  Number of grid cells:", nrow(merged_counts_1deg_carb), "\n")
cat("  Mean proportion:", round(mean(merged_counts_1deg_carb$proportion, na.rm = TRUE), 4), "\n")
cat("  Median proportion:", round(median(merged_counts_1deg_carb$proportion, na.rm = TRUE), 4), "\n")
cat("  Range:", round(min(merged_counts_1deg_carb$proportion, na.rm = TRUE), 4), "-",
    round(max(merged_counts_1deg_carb$proportion, na.rm = TRUE), 4), "\n")

cat("\n=== Comparison with GAM Predictions (0.25 degree) ===\n")
cat("\nGAM Physical Subduction (0.25 degree):\n")
cat("  Mean proportion:", round(mean(pred_full_subd$proportion, na.rm = TRUE), 4), "\n")
cat("  Median proportion:", round(median(pred_full_subd$proportion, na.rm = TRUE), 4), "\n")

cat("\nGAM Carbon Subduction (0.25 degree):\n")
cat("  Mean proportion:", round(mean(pred_full_carb$proportion, na.rm = TRUE), 4), "\n")
cat("  Median proportion:", round(median(pred_full_carb$proportion, na.rm = TRUE), 4), "\n")

cat("\nFigure saved to: ../pubfig/supplementary_resolution_comparison.png\n")
