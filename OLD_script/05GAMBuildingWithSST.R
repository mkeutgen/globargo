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
df_wo_sst <- read_csv("~/Documents/GLOBARGO/data/df_N2_MLD_binned.csv")
df_wo_sst
df <- read_csv("~/Documents/GLOBARGO/data/subduction_events_with_sst_gradient.csv")

# By assumption the elementary statistical observation is the occurence of subduction
# in the 1/4 degree gridcell (approx 25 km x 25 km)


# -----------------------------
# 3. Binning Data into 0.25-Degree Grid Cells
# -----------------------------

# Create binning variables
df <- df %>%
  mutate(
    LAT_BIN = floor(LATITUDE * 2) / 2 + 0.25,    # Center of 0.5-degree bin
    LON_BIN = floor(LONGITUDE * 2) / 2 + 0.25
  )

# Aggregate data by grid cells
df_binned <- df_full %>%
  group_by(LAT_BIN, LON_BIN, AdjustedDayOfYear, Hemisphere) %>%
  summarize(
    Anomaly = mean(Anomaly),
    cleaned_mld = mean(cleaned_mld),
    log_cleaned_mld = mean(log_cleaned_mld),
    LATITUDE = mean(LATITUDE),
    LONGITUDE = mean(LONGITUDE),
    .groups = 'drop'
  )



