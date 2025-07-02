# ============================================
# R Script for Fitting Multiple GAM Models (Optimized for HPC)
# ============================================

# 1. Library Loading and Setup
# ----------------------------

# Load necessary libraries
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

library(tidyverse)    # Data manipulation and visualization
library(mgcv)         # GAM modeling
library(parallel)     # Parallel computing

# 2. Data Reading and Preparation
# -------------------------------

# Read data
df_full <- read_csv("dataframe_full_mld_and_n2.csv")

# Standardize predictors and create Season factor
df_full <- df_full %>%
  mutate(
    log_cleaned_mld_std = as.vector(scale(log(cleaned_mld))),
    log_cleaned_N2_std = as.vector(scale(log(cleaned_N2))),
    Season = cut(
      AdjustedDayOfYear,
      breaks = c(0, 90, 180, 270, 365),
      labels = c("Winter", "Spring", "Summer", "Fall"),
      include.lowest = TRUE
    ),
    Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall"))
  )

# Define bin size (0.25 degrees ~25 km)
bin_size <- 0.25

# Define longitude and latitude bins
longitude_bins <- seq(floor(min(df_full$LONGITUDE, na.rm = TRUE)),
                      ceiling(max(df_full$LONGITUDE, na.rm = TRUE)),
                      by = bin_size)
latitude_bins <- seq(floor(min(df_full$LATITUDE, na.rm = TRUE)),
                     ceiling(max(df_full$LATITUDE, na.rm = TRUE)),
                     by = bin_size)

# Compute bin centers for labeling
lon_centers <- longitude_bins[-length(longitude_bins)] + bin_size / 2
lat_centers <- latitude_bins[-length(latitude_bins)] + bin_size / 2

# Assign bins
df_full <- df_full %>%
  mutate(
    lon_bin = cut(LONGITUDE, breaks = longitude_bins, include.lowest = TRUE, labels = lon_centers),
    lat_bin = cut(LATITUDE, breaks = latitude_bins, include.lowest = TRUE, labels = lat_centers),
    lon_bin = as.numeric(as.character(lon_bin)),
    lat_bin = as.numeric(as.character(lat_bin))
  )

df_full$TIME
# Aggregate data
df_agg <- df_full %>%
  group_by(lat_bin, lon_bin, TIME,Month) %>%
  summarise(
    Anomaly = sum(Anomaly),
    Total = n(),
    cleaned_mld = mean(cleaned_mld, na.rm = TRUE),
    cleaned_N2 = mean(cleaned_N2, na.rm = TRUE),
    LATITUDE = mean(LATITUDE, na.rm = TRUE),
    LONGITUDE = mean(LONGITUDE, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    Season = case_when(
      Month %in% c(12, 1, 2) ~ "DJF",
      Month %in% c(3, 4, 5)  ~ "MAM",
      Month %in% c(6, 7, 8)  ~ "JJA",
      Month %in% c(9, 10, 11) ~ "SON"
    ),
    Season = factor(Season, levels = c("DJF", "MAM", "JJA", "SON")),
    Hemisphere = ifelse(lat_bin > 0,"Northern","Southern")
  )


# Remove rows with zero Total to avoid issues in modeling
df_agg <- df_agg %>%
  filter(Total > 0)

# Standardize aggregated predictors
df_agg <- df_agg %>%
  mutate(
    log_cleaned_mld_std = as.vector(scale(log(cleaned_mld))),
    log_cleaned_N2_std = as.vector(scale(log(cleaned_N2)))
  )

# 3. Parallel Processing Setup
# ----------------------------

# Set the number of cores to use
num_cores <- 32  # Match the number of cores requested in your Slurm script

# 4. Model Specification and Fitting with Saving
# ----------------------------------------------

# Directory to save models
model_dir <- "models"
if (!dir.exists(model_dir)) {
  dir.create(model_dir)
}

# List to store models
models <- list()

# Function to fit, save, and delete model from memory
fit_and_save_model <- function(model_name, model_expression) {
  # Try to fit the model
  models[[model_name]] <<- tryCatch(
    {
      eval(parse(text = model_expression))
    },
    error = function(e) {
      message("Error occurred while fitting ", model_name, ": ", e$message)
      return(NULL)
    }
  )
  
  # If model is successfully fitted, save it and remove it from memory
  if (!is.null(models[[model_name]])) {
    saveRDS(models[[model_name]], file = file.path(model_dir, paste0(model_name, ".rds")))
    message("Model '", model_name, "' has been successfully saved.")
    
    # Remove model from memory
    models[[model_name]] <<- NULL
    gc() # Force garbage collection
  } else {
    message("Model '", model_name, "' could not be fitted and was not saved.")
  }
}

# Adjusted k values to manageable sizes (e.g., k = 200)

# 4.1 Basic GAM with 'sos' basis and optimizations
fit_and_save_model(
  "basic_gam_sos",
  "
  bam(
    cbind(Anomaly, Total - Anomaly) ~ Hemisphere+
      s(log_cleaned_N2_std, bs = 'cr', k = 20) +
      s(log_cleaned_mld_std, bs = 'cr', k = 20) +
      s(LATITUDE, LONGITUDE, bs = 'sos', k = 3000) +
      Season,
    family = binomial,
    weights = Total,
    data = df_agg,
    method = 'fREML',
    discrete = TRUE,
    cluster = num_cores
  )
  "
)

# 4.2 GAM with Spatial Smooth Varying by Season
fit_and_save_model(
  "spatial_by_season",
  "
  bam(
    cbind(Anomaly, Total - Anomaly) ~ Hemisphere +
      s(log_cleaned_N2_std, bs = 'cr', k = 20) +
      s(log_cleaned_mld_std, bs = 'cr', k = 20) +
      s(LATITUDE, LONGITUDE, bs = 'sos', by = Season, k = 3000) +
      Season,
    family = binomial,
    weights = Total,
    data = df_agg,
    method = 'fREML',
    discrete = TRUE,
    cluster = num_cores
  )
  "
)

# 4.3 GAM with Temporal Smooth
fit_and_save_model(
  "temporal_smooth",
  "
  bam(
    cbind(Anomaly, Total - Anomaly) ~ Hemisphere +
      s(log_cleaned_N2_std, bs = 'cr', k = 20) +
      s(log_cleaned_mld_std, bs = 'cr', k = 20) +
      s(LATITUDE, LONGITUDE, bs = 'sos', k = 3000,by= Season) +
      s(Month, bs = 'cc', k = 12) +
      Season,
    family = binomial,
    weights = Total,
    data = df_agg,
    method = 'fREML',
    cluster = num_cores
  )
  "
)

# 4.4 GAM with Interaction between Predictors and Season
fit_and_save_model(
  "predictor_by_season",
  "
  bam(
    cbind(Anomaly, Total - Anomaly) ~ Hemisphere +
      s(log_cleaned_N2_std, by = Season, bs = 'cr', k = 20) +
      s(log_cleaned_mld_std, by = Season, bs = 'cr', k = 20) +
      s(LATITUDE, LONGITUDE, bs = 'sos', k = 3000,by = Season) +
      Season,
    family = binomial,
    weights = Total,
    data = df_agg,
    method = 'fREML',
    cluster = num_cores
  )
  "
)

# 4.5 GAM with Spatial-Temporal Interaction
fit_and_save_model(
  "space_time_interaction",
  "
  bam(
    cbind(Anomaly, Total - Anomaly) ~ Hemisphere +
      s(log_cleaned_N2_std, bs = 'cr', k = 20) +
      s(log_cleaned_mld_std, bs = 'cr', k = 20) +
      te(LATITUDE, LONGITUDE, Month, bs = c('tp', 'tp', 'cc'), k = c(50, 50, 12)) +
      Season,
    family = binomial,
    weights = Total,
    data = df_agg,
    method = 'fREML',
    discrete = TRUE,
    cluster = num_cores
  )
  "
)

# 4.6 GAM with Random Effects
fit_and_save_model(
  "random_effects",
  "
  bam(
    cbind(Anomaly, Total - Anomaly) ~ Hemisphere + 
      s(log_cleaned_N2_std, bs = 'cr', k = 20) +
      s(log_cleaned_mld_std, bs = 'cr', k = 20) +
      s(LATITUDE, LONGITUDE, bs = 'sos', k = 3000) +
      s(lat_bin, lon_bin, bs = 're') +  # Random effect for grid cells
      s(Month, bs = 're') +             # Random effect for months
      Season,
    family = binomial,
    weights = Total,
    data = df_agg,
    method = 'fREML',
    cluster = num_cores
  )
  "
)

# 4.7 GAM with Quasibinomial Family (Overdispersion)
fit_and_save_model(
  "quasibinomial",
  "
  bam(
    cbind(Anomaly, Total - Anomaly) ~ Hemisphere +
      s(log_cleaned_N2_std, bs = 'cr', k = 20) +
      s(log_cleaned_mld_std, bs = 'cr', k = 20) +
      s(LATITUDE, LONGITUDE, bs = 'sos', k = 1000) +
      Season,
    family = quasibinomial(),
    weights = Total,
    data = df_agg,
    method = 'fREML',
    cluster = num_cores
  )
  "
)

# ============================================
# End of Script
# ============================================
