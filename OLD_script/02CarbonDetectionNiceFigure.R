# CARBON DETECTION ALGORITHM - CAT2 Events #
# MK #


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
# Resolve function conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")


# Source the functions required for the algorithm 
source(file = "/data/GLOBARGO/src/02AnnexDetFunBBP700.R")


# Function to downscale data without outlier detection
downscale_data_fun_wo_out <- function(df, bin_width = 40) {
  data <- df %>%
    select(PRES_ADJUSTED, AOU, ABS_SAL,BBP700_ADJUSTED,SPIC, CYCLE_NUMBER, LONGITUDE, LATITUDE, TIME)
  
  pressure_range <- range(data$PRES_ADJUSTED, na.rm = TRUE)
  bins <- seq(
    from = floor(pressure_range[1] / bin_width) * bin_width,
    to = ceiling(pressure_range[2] / bin_width) * bin_width,
    by = bin_width
  )
  
  data$bin <- cut(data$PRES_ADJUSTED, breaks = bins, include.lowest = TRUE, labels = FALSE)
  
  downscaled_data <- data %>%
    group_by(bin) %>%
    mutate(across(
      .cols = -c(PRES_ADJUSTED, LATITUDE, LONGITUDE, CYCLE_NUMBER, TIME),
      .fns = mean,
      na.rm = TRUE
    ))
  
  downscaled_data$PRES_ADJUSTED <- (bins[downscaled_data$bin] + bins[downscaled_data$bin + 1]) / 2
  
  downscaled_data <- unique(downscaled_data)
  return(downscaled_data)
}

# Hyperparameters :
cutoff <- 1.96
resolution <- 40
window <- 60 
#
#5906635


# Fill WMO 
# Example float : WMO = 5904179, CYCLE_NUMBER = 54
wmo <- 5904179
cycle_number <- 54
list_plots <- list()

# False positive in SPIC but not in SAL 
wmo <- 1902455
cycle_number <- 37
# Load data for the current float ('load_float_data' is an internal function to the OneArgo R toolbox)
float_data <- load_float_data(
  float_ids = wmo,
  variables = c(
    "BBP700","BBP700_ADJUSTED",
    "BBP700_ADJUSTED_QC", "LATITUDE", "LONGITUDE",
    "PRES",
    "PRES_ADJUSTED", "CHLA_ADJUSTED",
    "DOXY_ADJUSTED", "PSAL",
    "PSAL_ADJUSTED", "TEMP",
    "TEMP_ADJUSTED","DOXY","DOXY_ADJUSTED"
  ),
  format = "dataframe"
)

float_data_copy <- float_data %>% select(LATITUDE,LONGITUDE,TIME,CYCLE_NUMBER) %>% unique()
float_data_copy$SUBD <- ifelse(float_data_copy$CYCLE_NUMBER == cycle_number,T,F)

float_data_copy

library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)


# 1. Convert float data to sf (WGS84)
float_data_sf <- st_as_sf(
  float_data_copy,
  coords = c("LONGITUDE", "LATITUDE"),
  crs = 4326,   # WGS84
  remove = FALSE
)

# 2. Filter latitudes south of 30°S
float_data_sf_30S <- float_data_sf %>% 
  filter(LATITUDE < -30)

# 3. Get and crop the world map for lat < -30 (this gives the continents)
world <- ne_countries(scale = "medium", returnclass = "sf")
world_30S <- st_crop(
  world,
  xmin = -180, xmax = 180,  # full longitude range
  ymin = -90,  ymax = -30   # latitudes from -90° to -30°
)

# 4. Plot in southern polar stereographic projection
ggplot() +
  # A) Draw continents
  geom_sf(data = world, fill = "lightgrey", color = "lightgrey", inherit.aes = FALSE)+  
  # B) Plot black points for SUBD = FALSE (first layer)
  geom_sf(
    data = filter(float_data_sf_30S, SUBD == FALSE),
    color = "black",
    size = 2,
    alpha = 0.8
  ) +
  
  # C) Plot red points for SUBD = TRUE (second layer, on top)
  geom_sf(
    data = filter(float_data_sf_30S, SUBD == TRUE),
    color = "red",
    size = 3,
    alpha = 1
  ) +
  
  # D) Polar stereographic projection, centered on the South Pole
  coord_sf(
    crs = st_crs("+proj=stere +lat_0=-90 +lon_0=0 +lat_ts=-70 +datum=WGS84 +units=m +no_defs"),
    xlim = c(-6.5e6, 6.5e6),
    ylim = c(-6.5e6, 6.5e6),
    expand = TRUE
  ) +
  
  labs(
    title = "Argo Float  Location (WMO 5904179)",
    subtitle = "Red = subduction event"
  ) +
  theme_minimal()

# Select only the cycle number which are anomalous for efficiency

float_data <- float_data %>% filter(float_data$CYCLE_NUMBER %in% cycle_number) 


# Preprocess data: calculate additional oceanographic parameters
float_data <- float_data %>%
  filter(!is.na(DOXY)) %>%
  group_by(CYCLE_NUMBER) %>%
  mutate(
    ABS_SAL = gsw::gsw_SA_from_SP(
      SP = PSAL_ADJUSTED,
      p = PRES_ADJUSTED,
      longitude = first(LONGITUDE),
      latitude = first(LATITUDE)
    ),
    CONS_TEMP = gsw::gsw_CT_from_t(
      SA = ABS_SAL,
      t = TEMP_ADJUSTED,
      p = PRES_ADJUSTED
    ),
    SAT_DOXY = gsw_O2sol(
      SA = ABS_SAL,
      CT = CONS_TEMP,
      p = PRES_ADJUSTED,
      longitude = first(LONGITUDE),
      latitude = first(LATITUDE)
    ),
    AOU = SAT_DOXY - DOXY_ADJUSTED,
    SIGMA0 = gsw::gsw_sigma0(SA = ABS_SAL, CT = CONS_TEMP),
    SPIC = swSpice(salinity = PSAL_ADJUSTED, temperature = TEMP_ADJUSTED,
                   latitude = first(LATITUDE), longitude = first(LONGITUDE), eos = "unesco")
  )



# Downscale data without outlier detection
cycles_list <- float_data %>%
  group_by(CYCLE_NUMBER) %>%
  group_split()


# For Carbon Detection, we first apply a moving median (k = 3) filter to filter small anomalies
cycles_list <- lapply(cycles_list, function(df) {
  df$BBP700_ADJUSTEDORIG <- df$BBP700_ADJUSTED
  df$BBP700_ADJUSTED <- rollmedian(df$BBP700_ADJUSTED, k = 3, fill = NA)
  return(df)
})    




# Downscale Vertical resolution to 20 meters

downscaled_data_list_20m <- lapply(cycles_list, downscale_data_fun_wo_out,bin_width = 20)
downscaled_data_20m <- bind_rows(downscaled_data_list_20m)
residuals_data_20m <- downscaled_data_20m %>%
  group_by(CYCLE_NUMBER) %>%
  group_modify(~ .x %>%
                 select(PRES_ADJUSTED, AOU,BBP700_ADJUSTED, ABS_SAL,SPIC, LATITUDE, LONGITUDE, TIME) %>%
                 pivot_longer(
                   cols = c(AOU, ABS_SAL,BBP700_ADJUSTED,SPIC),
                   names_to = "VAR",
                   values_to = "VALUE"
                 ) %>%
                 group_by(VAR) %>%
                 mutate(
                   MA_3 = rollmean(VALUE, 3, fill = NA),
                   TM_9 = rollapply(VALUE, 9,
                                    function(x) {
                                      x_subset <- x[
                                        x >= quantile(x, 0.2, na.rm = TRUE) &
                                          x <= quantile(x, 0.8, na.rm = TRUE)
                                      ]
                                      if (length(x_subset) > 0) {
                                        mean(x_subset, na.rm = TRUE)
                                      } else {
                                        NA
                                      }
                                    },
                                    fill = NA
                   ),
                   ROB_RES = MA_3 - TM_9,
                   ROB_RES_RAW = VALUE - TM_9
                 ) %>%
                 mutate(
                   IQRN = IQR(ROB_RES_RAW, na.rm = TRUE) / 1.349,
                   MEDIAN_RES = median(ROB_RES_RAW[ROB_RES_RAW != 0], na.rm = TRUE)
                 ) %>%
                 mutate(
                   SCALE_RES_ROB = ifelse(
                     ROB_RES_RAW == 0,
                     0,
                     (ROB_RES_RAW - MEDIAN_RES) / IQRN
                   )
                 )
  ) %>%
  ungroup()

residuals_data_wf_20m <- residuals_data_20m %>%
  group_by(CYCLE_NUMBER) %>%
  select(PRES_ADJUSTED, SCALE_RES_ROB,VALUE, VAR, CYCLE_NUMBER, LONGITUDE,
         LATITUDE, TIME) %>%
  pivot_wider(names_from = VAR, values_from = c(SCALE_RES_ROB,VALUE))


# outlier calculation is on 40 m
downscaled_data_list <- lapply(cycles_list, downscale_data_fun_wo_out,bin_width = 40)
downscaled_data <- bind_rows(downscaled_data_list)

# Prepare data for residuals calculation
residuals_data <- downscaled_data %>%
  group_by(CYCLE_NUMBER) %>%
  group_modify(~ .x %>%
                 select(PRES_ADJUSTED, AOU,BBP700_ADJUSTED, ABS_SAL, LATITUDE, LONGITUDE, TIME) %>%
                 pivot_longer(
                   cols = c(AOU, ABS_SAL,BBP700_ADJUSTED),
                   names_to = "VAR",
                   values_to = "VALUE"
                 ) %>%
                 group_by(VAR) %>%
                 mutate(
                   MA_3 = rollmean(VALUE, 3, fill = NA),
                   TM_9 = rollapply(VALUE, 9,
                                    function(x) {
                                      x_subset <- x[
                                        x >= quantile(x, 0.2, na.rm = TRUE) &
                                          x <= quantile(x, 0.8, na.rm = TRUE)
                                      ]
                                      if (length(x_subset) > 0) {
                                        mean(x_subset, na.rm = TRUE)
                                      } else {
                                        NA
                                      }
                                    },
                                    fill = NA
                   ),
                   ROB_RES = MA_3 - TM_9,
                   ROB_RES_RAW = VALUE - TM_9
                 ) %>%
                 mutate(
                   IQRN = IQR(ROB_RES_RAW, na.rm = TRUE) / 1.349,
                   MEDIAN_RES = median(ROB_RES_RAW[ROB_RES_RAW != 0], na.rm = TRUE)
                 ) %>%
                 mutate(
                   SCALE_RES_ROB = ifelse(
                     ROB_RES_RAW == 0,
                     0,
                     (ROB_RES_RAW - MEDIAN_RES) / IQRN
                   )
                 )
  ) %>%
  ungroup()



# Detect outliers in downscaled data at 40 meters vertical res and pivot to wide format
residuals_data_wf <- residuals_data %>%
  group_by(CYCLE_NUMBER) %>%
  select(PRES_ADJUSTED, SCALE_RES_ROB,VALUE, VAR, CYCLE_NUMBER, LONGITUDE,
         LATITUDE, TIME) %>%
  pivot_wider(names_from = VAR, values_from = c(SCALE_RES_ROB,VALUE)) %>%
  mutate(OUT.T = ifelse(SCALE_RES_ROB_BBP700_ADJUSTED > 0 & abs(SCALE_RES_ROB_BBP700_ADJUSTED) > 1 &
                          abs(SCALE_RES_ROB_AOU) > cutoff & abs(SCALE_RES_ROB_ABS_SAL) > cutoff & SCALE_RES_ROB_AOU < 0,1,0))


# Identify potential eddy events
potential_eddy_events <- residuals_data_wf %>%
  filter(OUT.T == 1) %>%
  select(CYCLE_NUMBER, PRES_ADJUSTED, LATITUDE,
         LONGITUDE, TIME,SCALE_RES_ROB_ABS_SAL,
         SCALE_RES_ROB_AOU,SCALE_RES_ROB_BBP700_ADJUSTED,
         VALUE_ABS_SAL,VALUE_AOU,VALUE_BBP700_ADJUSTED) %>%
  unique() %>%
  filter(PRES_ADJUSTED >= 200, PRES_ADJUSTED <= 1000)


# Initialize columns to store the check results
potential_eddy_events$AOU_gradient_sign_change <- NA
potential_eddy_events$AOU_second_derivative_check <- NA
potential_eddy_events$ABS_SAL_gradient_sign_change <- NA
potential_eddy_events$ABS_SAL_second_derivative_check <- NA

# Loop over each event
for (i in seq_len(nrow(potential_eddy_events))) {
  cycle_num <- potential_eddy_events$CYCLE_NUMBER[i]
  pres_level <- potential_eddy_events$PRES_ADJUSTED[i]
  
  # Get the profile data at 20m resolution for the current cycle
  profile_data <- residuals_data_wf_20m %>% filter(CYCLE_NUMBER == cycle_num)
  
  # Perform checks for AOU
  checks_aou <- perform_checks(profile_data, target_level = pres_level,
                               variable_name = "VALUE_AOU",
                               second_deriv = 0.001,window=60)
  potential_eddy_events$AOU_gradient_sign_change[i] <- checks_aou$gradient_sign_change
  potential_eddy_events$AOU_second_derivative_check[i] <- checks_aou$second_derivative_check
  
  # Perform checks for ABS_SAL
  checks_abs_sal <- perform_checks(profile_data, target_level = pres_level,
                                   variable_name = "VALUE_ABS_SAL",second_deriv = 0.001,
                                   window=60)
  potential_eddy_events$ABS_SAL_gradient_sign_change[i] <- checks_abs_sal$gradient_sign_change
  potential_eddy_events$ABS_SAL_second_derivative_check[i] <- checks_abs_sal$second_derivative_check
}


filtered_events <- potential_eddy_events %>%
  filter(
    AOU_gradient_sign_change == TRUE,
    ABS_SAL_gradient_sign_change == TRUE)

# Initialize empty lists for plots
prof_plot <- list()
res_plot <- list()
list_plots <- list()

# Plotting profiles and residuals
for (i in seq_along(filtered_events$CYCLE_NUMBER)) {
  current_cycle <- filtered_events$CYCLE_NUMBER[i]
  current_data_20m <- residuals_data_20m %>% filter(CYCLE_NUMBER == current_cycle)
  current_data_40m <- residuals_data %>% filter(CYCLE_NUMBER == current_cycle)
  current_eddy <- filtered_events %>% filter(CYCLE_NUMBER == current_cycle)
  pres_level <- filtered_events$PRES_ADJUSTED[i]
  
  var_labels <- c(
    "AOU" = "Apparent Oxygen Utilization (µmol/kg)",
    "ABS_SAL" = "Absolute Salinity (g/kg)",
    "BBP700_ADJUSTED" = "Backscatter at 700nm (m⁻¹ sr⁻¹)",
    "SPIC" = "Spiciness (kg m⁻³)"
  )
  current_data_20m$VAR <- factor(current_data_20m$VAR, levels = c("ABS_SAL", "SPIC", "AOU", "BBP700_ADJUSTED"))
  # Use consistent colors across all plots
  custom_colors <- c("Observed Values" = "red", "Trimmed Mean (k=9)" = "blue", 
                     "-2 sigma" = "blue", "+2 sigma" = "blue")
  
  # Plotting profiles
  prof_plot[[i]] <- current_data_20m %>%
    ggplot(aes(x = PRES_ADJUSTED, y = VALUE)) +
    facet_wrap(. ~ VAR, scales = "free", labeller = labeller(VAR = var_labels)) +
    coord_flip() +
    scale_x_reverse(limits = c(600, 0), breaks = seq(0, 600, by = 200)) +
    geom_line(aes(y = VALUE, color = "Observations"), size = 1.2) +
    geom_line(aes(y = TM_9, color = "Trimmed Mean (k=9)"), size = 1.2) +
    geom_point(aes(y = VALUE, color = "Observations"), size = 1.2) +
    geom_point(aes(y = TM_9, color = "Trimmed Mean (k=9)"), size = 1.2) +
    theme_bw() +
    labs(x = "Adjusted Pressure (dbar)", y = "") +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),                 # Remove legend title
      text = element_text(size = 16),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      strip.text = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.spacing = unit(3, "lines"),
    ) +
    geom_vline(xintercept = 340, color = "darkgreen", alpha = 0.3, size = 1)  
  

  false_pos_in_spic_not_in_sal <- current_data_20m %>%
    ggplot(aes(x = PRES_ADJUSTED, y = VALUE)) +
    facet_wrap(. ~ VAR, scales = "free", labeller = labeller(VAR = var_labels)) +
    coord_flip() +
    scale_x_reverse(limits = c(600, 0), breaks = seq(0, 600, by = 200)) +
    geom_line(aes(y = VALUE, color = "Observations"), size = 1.2) +
    geom_point(aes(y = VALUE, color = "Observations"), size = 1.2) +
    theme_bw() +
    labs(x = "Adjusted Pressure (dbar)", y = "") +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),                 # Remove legend title
      text = element_text(size = 22),
      axis.text = element_text(size = 25),
      axis.title = element_text(size = 25),
      strip.text = element_text(size = 25),
      legend.text = element_text(size = 25),
      legend.spacing = unit(3, "lines"),
    ) +
    geom_vline(xintercept = 340, color = "darkgreen", alpha = 0.3, size = 1)  
  
  
  
  # Prepare data for residuals plot
  
  
  hline_data <- data.frame(
    VAR = c("AOU", "AOU", "ABS_SAL", "ABS_SAL", "BBP700_ADJUSTED", "BBP700_ADJUSTED"),
    hline = c(-2, 2, -2, 2, -1, 1),
    label = c("-2 sigma", "+2 sigma", "-2 sigma", "+2 sigma", "-1 sigma", "+1 sigma")
  )
  
  
#  # Plotting residuals
#  res_plot[[i]] <- current_data_40m %>%
#    ggplot(aes(x = PRES_ADJUSTED, y = SCALE_RES_ROB)) +
#    scale_x_reverse(limits = c(1200, 0), breaks = seq(0, 1200, by = 40)) +
#    facet_grid(. ~ VAR, scales = "free") +
#    coord_flip() +
#    theme_bw() +
#    geom_point() +
#    labs(x = "Adjusted Pressure (dbar)", y = "") +
#    geom_hline(data = hline_data, aes(yintercept = hline, color = label), inherit.aes = TRUE) +
#    scale_color_manual(name = "Threshold", values = c("-2 sigma" = "blue", "+2 sigma" = "blue")) +
#    theme(legend.position = "bottom") +
#    geom_vline(xintercept = current_eddy$PRES_ADJUSTED, color = "darkgreen", alpha = 0.3, size = 1)
  
  
  # Annotation text
  annotation_text <- paste(
    "Cycle Number:", current_cycle,
    "\nFloat ID (WMO):", wmo,
    "\nLongitude:", current_eddy$LONGITUDE,
    "\nLatitude:", current_eddy$LATITUDE,
    "\nTime:", format(as.POSIXct(current_eddy$TIME, origin = "1970-01-01"), "%Y-%m-%d")
  )
  
  # Combine plots with annotation
  # If you want residual plots and profile plots, uncomment : 
  # combined_plot <- ggarrange(
  #   prof_plot[[i]],
  #   res_plot[[i]],
  #   common.legend = FALSE,
  #   legend = "bottom",
  #   nrow = 2
  # )
  combined_plot <- ggarrange(
      prof_plot[[i]],
      common.legend = FALSE,
      legend = "bottom",
      nrow = 2
    )
  combined_plot <- annotate_figure(
    combined_plot,
    top = text_grob(annotation_text, face = "bold", size = 10)
  )
  ggsave(plot = false_pos_in_spic_not_in_sal,
         filename = "/data/GLOBARGO/src/figures/false_pos_in_spic.png",width = 18,height = 15)
  
  list_plots[[i]] <- combined_plot
}

# Save plots to files
if (length(list_plots) > 0) {
  dir <- paste0("/data/GLOBARGO/figures/NicePlots/", wmo)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  
  for (k in seq_along(list_plots)) {
    cycle_number <- filtered_events$CYCLE_NUMBER[k]
    file_name <- paste0(dir, "/", wmo, "_plot_cycl
# Function to downscale data without outlier detection
downscale_data_fun_wo_out <- function(df, bin_width = 40) {
  data <- df %>%
    select(PRES_ADJUSTED, AOU, ABS_SAL,BBP700_ADJUSTED, CYCLE_NUMBER, LONGITUDE, LATITUDE, TIME)
  
  pressure_range <- range(data$PRES_ADJUSTED, na.rm = TRUE)
  bins <- seq(
    from = floor(pressure_range[1] / bin_width) * bin_width,
    to = ceiling(pressure_range[2] / bin_width) * bin_width,
    by = bin_width
  )
  
  data$bin <- cut(data$PRES_ADJUSTED, breaks = bins, include.lowest = TRUE, labels = FALSE)
  
  downscaled_data <- data %>%
    group_by(bin) %>%
    mutate(across(
      .cols = -c(PRES_ADJUSTED, LATITUDE, LONGITUDE, CYCLE_NUMBER, TIME),
      .fns = mean,
      na.rm = TRUE
    ))
  
  downscaled_data$PRES_ADJUSTED <- (bins[downscaled_data$bin] + bins[downscaled_data$bin + 1]) / 2
  
  downscaled_data <- unique(downscaled_data)
  return(downscaled_data)
}

# Function to compute mean ABS_SAL in the upper 50 meters
compute_mean_ABS_SAL_50m <- function(ABS_SAL_df) {
  mean_ABS_SAL <- ABS_SAL_df %>%
    filter(PRES_ADJUSTED <= 50) %>%
    summarize(mean_ABS_SAL = mean(ABS_SAL, na.rm = TRUE)) %>%
    pull(mean_ABS_SAL)
  
  return(mean_ABS_SAL)
}

# Function to compute mean ABS_SAL at min and max pressure levels around a target pressure
mean_ABS_SAL_at_min_max_levels <- function(data, target_pressure) {
  filtered_data <- data %>%
    filter(!is.na(ABS_SAL)) %>%
    filter(PRES_ADJUSTED >= (target_pressure - 100), PRES_ADJUSTED <= (target_pressure + 100))
  
  min_pressure <- min(filtered_data$PRES_ADJUSTED, na.rm = TRUE)
  max_pressure <- max(filtered_data$PRES_ADJUSTED, na.rm = TRUE)
  
  ABS_SAL_values <- filtered_data %>%
    filter(PRES_ADJUSTED %in% c(min_pressure, max_pressure)) %>%
    pull(ABS_SAL)
  
  mean_ABS_SAL <- mean(ABS_SAL_values, na.rm = TRUE)
  
  return(mean_ABS_SAL)
}
e_", cycle_number, ".png")
    ggsave(file_name, list_plots[[k]], width = 13, height = 15)
  }
}


# Write output 
detected.events.df.carbon <- detected_events_list_carbon %>% bind_rows()

write_csv(detected.events.df.carbon,"/data/GLOBARGO/data/detected_events_unique_with_carbon_cat2.csv")



detected.events.df