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
# Resolve function conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Load detected events :
df_abs_sal <- read_csv("/data/GLOBARGO/data/detected_events_abs_sal_var_v5.csv")
df_spic <- read_csv("/data/GLOBARGO/data/detected_events_sens_and_spec_incr.csv")

# Classified datasets in spiciness
df_spic_class <- read_csv("/data/GLOBARGO/data/classification_results_v4.csv")
df_spic_class$WMO <- gsub("_plot", "", df_spic_class$WMO)
df_spic_class$CYCLE_NUMBER <- df_spic_class$Cycle

# Count events
df_abs_sal %>% select(WMO,CYCLE_NUMBER) %>% unique() #7,658 unique profiles when ABS_SAL is used
df_spic %>% select(WMO,CYCLE_NUMBER) %>% unique() # 8351 unique profiles when SPIC is used

df_not_in_spic <- df_abs_sal %>%
  anti_join(df_spic, by = c("WMO", "CYCLE_NUMBER"))


df_not_in_sal <- df_spic %>%
  anti_join(df_abs_sal, by = c("WMO", "CYCLE_NUMBER"))




df_not_in_spic %>% select(WMO,CYCLE_NUMBER) %>% unique() # 3,378 profiles in ABS_SAL But not in SPIC


# Step 2: Define the source and destination folder paths
source_folder <- "/data/GLOBARGO/figures/EddySubductionFiguresSalinityVarV4"
destination_folder <- "/data/GLOBARGO/figures/EddySubductionFiguresInSalinityButNotSpic"

# Create the destination folder if it doesn't exist
if (!dir_exists(destination_folder)) {
  dir_create(destination_folder)
}

# Step 3: Loop through the rows in df_not_in_spic and copy files
# Step 3: Loop through the rows in df_not_in_spic and copy files preserving the folder structure
df_not_in_spic %>%
  rowwise() %>%
  mutate(
    # Construct the source subfolder path
    source_subfolder = file.path(source_folder, as.character(WMO)),
    # Construct the source file path
    source_file = file.path(
      source_subfolder, 
      paste0(WMO, "_plot_cycle_", CYCLE_NUMBER, ".png")
    ),
    # Construct the destination subfolder path
    destination_subfolder = file.path(destination_folder, as.character(WMO)),
    # Construct the destination file path
    destination_file = file.path(
      destination_subfolder, 
      paste0(WMO, "_plot_cycle_", CYCLE_NUMBER, ".png")
    )
  ) %>%
  # Copy files and create subfolders if they don't exist
  do({
    if (file_exists(.$source_file)) {
      if (!dir_exists(.$destination_subfolder)) {
        dir_create(.$destination_subfolder)
      }
      file_copy(.$source_file, .$destination_file, overwrite = TRUE)
    }
    NULL
  })


# See if those classified as cat1 (strong anomalies) and cat2 (weak anomalies)
# in SPIC are also anomalies in ABS_SAL 

df_spic_class$CYCLE_NUMBER <- df_spic_class$Cycle

# Step 1: Filter anomalies IN SALINITY categorized as 1 or 2 in df_spic_class
df_spic_anomalies <- df_spic_class %>%
  filter(Category %in% c(1, 2))


# Step 2: Define the source and destination folder paths
source_folder <- "/data/GLOBARGO/figures/EddySubductionFiguresSalinityVarV4"
destination_folder <- "/data/GLOBARGO/figures/EddySubductionFiguresSpicClass1and2"

# Create the destination folder if it doesn't exist
if (!dir_exists(destination_folder)) {
  dir_create(destination_folder)
}


# Step 3: Loop through the rows in df_spic_anomalies and copy files preserving the folder structure
df_spic_anomalies %>%
  rowwise() %>%
  mutate(
    # Construct the source subfolder path based on WMO
    source_subfolder = file.path(source_folder, as.character(WMO)),
    # Construct the source file path based on WMO and Cycle (CYCLE_NUMBER)
    source_file = file.path(
      source_subfolder, 
      paste0(WMO, "_plot_cycle_", CYCLE_NUMBER, ".png")
    ),
    # Construct the destination subfolder path
    destination_subfolder = file.path(destination_folder, as.character(WMO)),
    # Construct the destination file path
    destination_file = file.path(
      destination_subfolder, 
      paste0(WMO, "_plot_cycle_", CYCLE_NUMBER, ".png")
    )
  ) %>%
  # Copy files and create subfolders if they don't exist
  do({
    if (file_exists(.$source_file)) {
      if (!dir_exists(.$destination_subfolder)) {
        dir_create(.$destination_subfolder)
      }
      file_copy(.$source_file, .$destination_file, overwrite = TRUE)
    }
    NULL
  })

# Extract subset of df_spic_anomalies with corresponding files
df_spic_anomalies_with_files <- df_spic_anomalies %>%
  rowwise() %>%
  mutate(
    # Construct the source subfolder path based on WMO
    source_subfolder = file.path(source_folder, as.character(WMO)),
    # Construct the source file path based on WMO and Cycle (CYCLE_NUMBER)
    source_file = file.path(
      source_subfolder, 
      paste0(WMO, "_plot_cycle_", CYCLE_NUMBER, ".png")
    )
  ) %>%
  # Filter rows where the source file exists
  filter(file_exists(source_file)) %>%
  ungroup()

# Now df_spic_anomalies_with_files contains only the rows with existing files
# Check that they are all contained in the df_abs_sal df : 
df_spic_anomalies_with_files %>% select(WMO,CYCLE_NUMBER)
df_abs_sal %>% select(WMO,CYCLE_NUMBER)


####################
# DEEP ANOMALIES ###
##

df_deep <- df_not_in_spic %>% filter(PRES_ADJUSTED >= 800)

# Step 2: Define the source and destination folder paths
source_folder <- "/data/GLOBARGO/figures/EddySubductionFiguresSalinityVarV5"
destination_folder <- "/data/GLOBARGO/figures/DeepEddySubductionFiguresInSalinity"

# Create the destination folder if it doesn't exist
if (!dir_exists(destination_folder)) {
  dir_create(destination_folder)
}

# Step 3: Loop through the rows in df_not_in_spic and copy files
# Step 3: Loop through the rows in df_not_in_spic and copy files preserving the folder structure
df_deep %>%
  rowwise() %>%
  mutate(
    # Construct the source subfolder path
    source_subfolder = file.path(source_folder, as.character(WMO)),
    # Construct the source file path
    source_file = file.path(
      source_subfolder, 
      paste0(WMO, "_plot_cycle_", CYCLE_NUMBER, ".png")
    ),
    # Construct the destination subfolder path
    destination_subfolder = file.path(destination_folder, as.character(WMO)),
    # Construct the destination file path
    destination_file = file.path(
      destination_subfolder, 
      paste0(WMO, "_plot_cycle_", CYCLE_NUMBER, ".png")
    )
  ) %>%
  # Copy files and create subfolders if they don't exist
  do({
    if (file_exists(.$source_file)) {
      if (!dir_exists(.$destination_subfolder)) {
        dir_create(.$destination_subfolder)
      }
      file_copy(.$source_file, .$destination_file, overwrite = TRUE)
    }
    NULL
  })


