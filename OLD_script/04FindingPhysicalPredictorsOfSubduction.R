# Load necessary libraries
library(gsw)
library(oce)
library(readr)
library(robustbase)
library(zoo)
library(lubridate)
library(tidyverse)
# Read the list of WMO IDs
wmolist <- readRDS("~/Documents/ARGO/BGC_Argo_WMO_PSAL_BBP_DOXY_TEMP.rds")

# Initialize a data frame to store results
mld_results <- data.frame(WMO = integer(), Latitude = numeric(), Longitude = numeric(), Time = as.POSIXct(character()), MLD = numeric())

# Main processing loop over each WMO ID
for (j in seq_along(wmolist)) {
  try({
    # Get current WMO ID
    wmo <- wmolist[j]
    
    # Load data for the current float (assuming 'load_float_data' is an internal function)
    float_data <- load_float_data(
      float_ids = wmo,
      variables = c(
        "DATA_TYPE", "PLATFORM_NUMBER", "BBP700", "BBP700_dPRES",
        "BBP700_ADJUSTED_QC", "LATITUDE", "LONGITUDE", "PROFILE_TEMP_QC",
        "PROFILE_DOXY_QC", "PROFILE_BBP700_QC", "PRES_QC", "PRES",
        "PRES_ADJUSTED", "PROFILE_PSAL_QC", "CHLA_QC", "CHLA_ADJUSTED",
        "CHLA_ADJUSTED_ERROR", "DOXY", "DOXY_QC", "DOXY_ADJUSTED",
        "DOXY_ADJUSTED_QC", "DOXY_ADJUSTED_ERROR", "PSAL", "PSAL_dPRES",
        "PSAL_ADJUSTED", "PSAL_ADJUSTED_QC", "TEMP", "TEMP_QC", "TEMP_dPRES",
        "TEMP_ADJUSTED", "TEMP_ADJUSTED_QC", "TEMP_ADJUSTED_ERROR"
      ),
      format = "dataframe"
    )
    
    # Preprocess data: calculate additional oceanographic parameters
    float_data <- float_data %>%
      filter(!is.na(PRES_ADJUSTED) & !is.na(PSAL_ADJUSTED) & !is.na(TEMP_ADJUSTED)) %>%
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
        SIGMA0 = gsw::gsw_sigma0(ABS_SAL, CONS_TEMP)  # Calculate density anomaly
      )
    
    # Define density threshold for MLD calculation
    density_threshold <- 0.03  # kg/m³ increase from the surface density
    
    # Loop over each profile (cycle) to compute MLD
    for (cycle_data in split(float_data, float_data$CYCLE_NUMBER)) {
      # Surface density
      surface_density <- cycle_data$SIGMA0[1]
      
      # Find the depth where density exceeds the surface density by the threshold
      mld_index <- which(cycle_data$SIGMA0 >= surface_density + density_threshold)[1]
      
      # Determine MLD based on the threshold, or set to NA if no MLD found
      mld <- if (!is.na(mld_index)) cycle_data$PRES_ADJUSTED[mld_index] else NA
      
      # Store WMO, Latitude, Longitude, Time, and MLD for the profile
      mld_results <- rbind(mld_results, data.frame(
        WMO = wmo,
        Latitude = first(cycle_data$LATITUDE),
        Longitude = first(cycle_data$LONGITUDE),
        Time = first(cycle_data$TIME),
        MLD = mld,
        CYCLE_NUMBER = first(cycle_data$CYCLE_NUMBER)
      ))
    }
  }, silent = TRUE)  # Continue to the next float if there's an error
}


mld_results <- read_csv("/data/GLOBARGO/data/mld_results.csv")

# remove ties in the data
mld_results <- mld_results %>%
  group_by(WMO, Time) %>%
  summarise(
    MLD = mean(MLD, na.rm = TRUE),
    CYCLE_NUMBER = first(CYCLE_NUMBER)  # Keep the first CYCLE_NUMBER for reference
  ) %>%
  ungroup()


mld_results <- mld_results %>% filter(MLD > 0)

mld_results %>% group_by(WMO) %>% select(WMO,Latitude,Longitude,Time,MLD,CYCLE_NUMBER)

mld_results %>% filter(WMO == 1901329)

# Calculate the lower and upper percentiles
lower_bound <- quantile(mld_results$MLD, 0.0025, na.rm = TRUE)
upper_bound <- quantile(mld_results$MLD, 0.9975, na.rm = TRUE)



# Create a new column 'cleaned_mld' that sets outliers to NA
mld_results <- mld_results %>%
  mutate(cleaned_mld = ifelse(MLD >= lower_bound & MLD <= upper_bound, MLD, NA))




# Plot histogram using ggplot2
mld_results %>% 
  filter(WMO == 5904677) %>%
  ggplot(aes(x = Time)) +
  geom_histogram(binwidth = 5, fill = "blue", color = "black") +
  scale_x_date(date_breaks = "1 week", date_labels = "%Y %m") +
  labs(title = "Distribution of Time for WMO 5904677", 
       x = "Date", 
       y = "Frequency") +
  theme_minimal()


# Ensure data is properly sorted by WMO, Time, and CYCLE_NUMBER
mld_results <- mld_results %>%
  arrange(WMO, Time, CYCLE_NUMBER)

# Compute 3-bin moving median for each WMO, using a grouped approach
mld_results <- mld_results %>%
  group_by(WMO) %>%
  mutate(
    MLD_MM3 = zoo::rollapply(
      cleaned_mld,
      width = 3,
      FUN = median,
      align = "center",
      partial = TRUE
    )
  ) %>%
  ungroup()

# For rows where the moving median couldn't be computed (i.e., at the edges), use the original cleaned_mld
mld_results$MLD_MM3[is.na(mld_results$MLD_MM3)] <- mld_results$cleaned_mld[is.na(mld_results$MLD_MM3)]

# Display the updated dataframe with the new column
print(head(mld_results, 10))


# Define bins for 'binned_mld' column following  https://doi.org/10.1029/2004JC002378
breaks <- c(0,10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 400, 500, Inf)
labels <- c("0-10","10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", 
            "90-100", "100-125", "125-150", "150-200", "200-300", "300-400", "400-500", "500+")
# Create the 'binned_mld' column using cut
mld_results <- mld_results %>%
  mutate(binned_mld = cut(MLD, breaks = breaks, labels = labels, right = FALSE))



# What about change in MLD? Assuming Argo are Lagrangian, we can compute the difference 
# between the 3 bin centered median MLD 

# Step 1: Sort the data by WMO and Time
# Step 1: Sort the data by WMO and Time
mld_results <- mld_results %>%
  arrange(WMO, Time)

# Step 2: Compute the centered finite difference approximation of the MLD derivative
mld_results <- mld_results %>%
  group_by(WMO) %>%
  mutate(
    # Convert time to numeric days since the first observation (to get consistent units)
    time_numeric = as.numeric(difftime(Time, min(Time), units = "days")),
    
    # Centered difference approximation for the time derivative using dplyr::lag and dplyr::lead
    MLD_rate_of_change = ifelse(
      row_number() > 1 & row_number() < n(),
      (dplyr::lead(MLD_MM3) - dplyr::lag(MLD_MM3)) / (dplyr::lead(time_numeric) - dplyr::lag(time_numeric)),
      NA
    )
  ) %>%
  ungroup()

mld_results <- mld_results %>% mutate(
  MLD_sign_rate_of_change = ifelse(MLD_rate_of_change > 10, 1, ifelse(MLD_rate_of_change < 10, 0,-1))
)


# Save to file
write.csv(mld_results, "/data/GLOBARGO/data/mld_results.csv", row.names = FALSE)




















# If the user is suspicious about quality of the MLD data and the corrections here proposed, they may
# choose to plot the MLD for a suspicious float. Example, very shallow MLD
### MLD VISUAL CHECK ######
wmo <- 6902733 # cycle number 21, mld is estimated at 3.7 m

# Select a few profiles to visualize (e.g., first 5 profiles with valid MLD)
sample_profiles <- mld_results %>% filter(WMO == wmo) %>%
  filter(!is.na(MLD)) 

# Function to plot density profiles with MLD marked
plot_density_profile <- function(wmo, cycle_number, float_data, mld) {
  # Filter data for the given WMO and cycle
  
  profile_data <- float_data %>%
    filter(CYCLE_NUMBER == cycle_number)
  
  # Generate the plot
  ggplot(profile_data, aes(x = SIGMA0, y = PRES_ADJUSTED)) +
    geom_line(color = "blue") +
    geom_hline(yintercept = mld, color = "red", linetype = "dashed") +
    scale_y_reverse() +  # Reverse y-axis so deeper depths are at the bottom
    labs(
      title = paste("Density Profile for WMO", wmo, "Cycle", cycle_number),
      x = "Density Anomaly (kg/m³)",
      y = "Pressure (dbar)"
    ) +
    theme_minimal() +
    annotate("text", x = max(profile_data$SIGMA0, na.rm = TRUE), 
             y = mld, label = paste("MLD =", round(mld, 2)), hjust = 1, vjust = -1)
}

# Loop through the sample profiles and generate plots
for (i in 1:nrow(sample_profiles)) {
  wmo <- sample_profiles$WMO[i]
  cycle_number <- sample_profiles$CYCLE_NUMBER[i]
  mld <- sample_profiles$MLD[i]
  
  # Load data for the specific float and cycle
  float_data <- load_float_data(
    float_ids = wmo,
    variables = c(
      "DATA_TYPE", "PLATFORM_NUMBER", "BBP700", "BBP700_dPRES",
      "BBP700_ADJUSTED_QC", "LATITUDE", "LONGITUDE", "PROFILE_TEMP_QC",
      "PROFILE_DOXY_QC", "PROFILE_BBP700_QC", "PRES_QC", "PRES",
      "PRES_ADJUSTED", "PROFILE_PSAL_QC", "CHLA_QC", "CHLA_ADJUSTED",
      "CHLA_ADJUSTED_ERROR", "DOXY", "DOXY_QC", "DOXY_ADJUSTED",
      "DOXY_ADJUSTED_QC", "DOXY_ADJUSTED_ERROR", "PSAL", "PSAL_dPRES",
      "PSAL_ADJUSTED", "PSAL_ADJUSTED_QC", "TEMP", "TEMP_QC", "TEMP_dPRES",
      "TEMP_ADJUSTED", "TEMP_ADJUSTED_QC", "TEMP_ADJUSTED_ERROR"
    ),
    format = "dataframe"
  )
  
  # Calculate density for the float data
  float_data <- float_data %>%
    mutate(
      ABS_SAL = gsw::gsw_SA_from_SP(PSAL_ADJUSTED, PRES_ADJUSTED, first(LONGITUDE), first(LATITUDE)),
      CONS_TEMP = gsw::gsw_CT_from_t(ABS_SAL, TEMP_ADJUSTED, PRES_ADJUSTED),
      SIGMA0 = gsw::gsw_sigma0(ABS_SAL, CONS_TEMP)
    )
  
  # Plot the density profile
  plot <- plot_density_profile(wmo, cycle_number, float_data, mld)
  
  # Save the plot
  ggsave(filename = paste0("/data/GLOBARGO/figures/MLD_fig/density_profile_wmo_", wmo, "_cycle_", cycle_number, ".png"),
         plot = plot, width = 8, height = 6)
  
  print(plot)
}








# Define the WMO ID and example cycle for visual check
wmo <- 6902733  # Example WMO ID with suspiciously shallow MLD
wmo <- 1902303 # Exemple WMO with MLD at 48 and 30 for cycle number 1 and 3 but only 6 m for cycle number 2
sample_profiles <- mld_results %>%
  filter(WMO == wmo) %>%
  filter(!is.na(MLD)) 

# Function to plot density profiles with MLD marked, including metadata annotation
plot_density_profile <- function(wmo, cycle_number, float_data, mld, latitude, longitude, time) {
  # Filter data for the given cycle
  profile_data <- float_data %>%
    filter(CYCLE_NUMBER == cycle_number)
  
  # Generate the plot
  ggplot(profile_data, aes(x = SIGMA0, y = PRES_ADJUSTED)) +
    geom_line(color = "blue") +
    geom_hline(yintercept = mld, color = "red", linetype = "dashed") +
    scale_y_reverse() +  # Reverse y-axis for depth
    labs(
      title = paste("Density Profile for WMO", wmo, "Cycle", cycle_number),
      subtitle = paste("Lat:", round(latitude, 2), "Lon:", round(longitude, 2), "Time:", format(time, "%Y-%m-%d %H:%M:%S")),
      x = "Density Anomaly (kg/m³)",
      y = "Pressure (dbar)"
    ) +
    theme_minimal() +
    annotate("text", x = max(profile_data$SIGMA0, na.rm = TRUE), y = mld, 
             label = paste("MLD =", round(mld, 2)), hjust = 1, vjust = -1, color = "red") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
}

# Loop through sample profiles to generate plots
for (i in 1:nrow(sample_profiles)) {
  wmo <- sample_profiles$WMO[i]
  cycle_number <- sample_profiles$CYCLE_NUMBER[i]
  mld <- sample_profiles$MLD[i]
  latitude <- sample_profiles$Latitude[i]
  longitude <- sample_profiles$Longitude[i]
  time <- sample_profiles$Time[i]
  
  # Load data for the specific float and cycle
  float_data <- load_float_data(
    float_ids = wmo,
    variables = c(
      "DATA_TYPE", "PLATFORM_NUMBER", "BBP700", "BBP700_dPRES",
      "BBP700_ADJUSTED_QC", "LATITUDE", "LONGITUDE", "PROFILE_TEMP_QC",
      "PROFILE_DOXY_QC", "PROFILE_BBP700_QC", "PRES_QC", "PRES",
      "PRES_ADJUSTED", "PROFILE_PSAL_QC", "CHLA_QC", "CHLA_ADJUSTED",
      "CHLA_ADJUSTED_ERROR", "DOXY", "DOXY_QC", "DOXY_ADJUSTED",
      "DOXY_ADJUSTED_QC", "DOXY_ADJUSTED_ERROR", "PSAL", "PSAL_dPRES",
      "PSAL_ADJUSTED", "PSAL_ADJUSTED_QC", "TEMP", "TEMP_QC", "TEMP_dPRES",
      "TEMP_ADJUSTED", "TEMP_ADJUSTED_QC", "TEMP_ADJUSTED_ERROR"
    ),
    format = "dataframe"
  )
  
  # Calculate density for the float data
  float_data <- float_data %>%
    mutate(
      ABS_SAL = gsw::gsw_SA_from_SP(PSAL_ADJUSTED, PRES_ADJUSTED, first(LONGITUDE), first(LATITUDE)),
      CONS_TEMP = gsw::gsw_CT_from_t(ABS_SAL, TEMP_ADJUSTED, PRES_ADJUSTED),
      SIGMA0 = gsw::gsw_sigma0(ABS_SAL, CONS_TEMP)
    )
  
  # Plot the density profile with metadata annotation
  plot <- plot_density_profile(wmo, cycle_number, float_data, mld, latitude, longitude, time)
  
  # Save the plot
  ggsave(filename = paste0("/data/GLOBARGO/figures/MLD_fig/density_profile_wmo_", wmo, "_cycle_", cycle_number, ".png"),
         plot = plot, width = 8, height = 6)
  
  # Print the plot for quick viewing
  print(plot)
}


###############################
#### N^2 (stratification) #####
###############################


# Initialize a data frame to store results
N2_results <- data.frame(WMO = integer(), Latitude = numeric(), Longitude = numeric(), Time = as.POSIXct(character()), N2 = numeric())
# Main processing loop over each WMO ID
for (j in seq_along(wmolist)) {
  try({
    # Get current WMO ID
    wmo <- wmolist[j]
    
    
    # Load data for the current float (assuming 'load_float_data' is an internal function)
    float_data <- load_float_data(
      float_ids = wmo,
      variables = c(
        "DATA_TYPE", "PLATFORM_NUMBER", "BBP700", "BBP700_dPRES",
        "BBP700_ADJUSTED_QC", "LATITUDE", "LONGITUDE", "PROFILE_TEMP_QC",
        "PROFILE_DOXY_QC", "PROFILE_BBP700_QC", "PRES_QC", "PRES",
        "PRES_ADJUSTED", "PROFILE_PSAL_QC", "CHLA_QC", "CHLA_ADJUSTED",
        "CHLA_ADJUSTED_ERROR", "DOXY", "DOXY_QC", "DOXY_ADJUSTED",
        "DOXY_ADJUSTED_QC", "DOXY_ADJUSTED_ERROR", "PSAL", "PSAL_dPRES",
        "PSAL_ADJUSTED", "PSAL_ADJUSTED_QC", "TEMP", "TEMP_QC", "TEMP_dPRES",
        "TEMP_ADJUSTED", "TEMP_ADJUSTED_QC", "TEMP_ADJUSTED_ERROR"
      ),
      format = "dataframe"
    )
    
    # Preprocess data: calculate additional oceanographic parameters
    float_data <- float_data %>%
      filter(!is.na(PRES_ADJUSTED) & !is.na(PSAL_ADJUSTED) & !is.na(TEMP_ADJUSTED)) %>%
      group_by(CYCLE_NUMBER) %>%
      group_modify(~ {
        df <- .x
        ABS_SAL <- gsw::gsw_SA_from_SP(
          SP = df$PSAL_ADJUSTED,
          p = df$PRES_ADJUSTED,
          longitude = first(df$LONGITUDE),
          latitude = first(df$LATITUDE)
        )
        CONS_TEMP <- gsw::gsw_CT_from_t(
          SA = ABS_SAL,
          t = df$TEMP_ADJUSTED,
          p = df$PRES_ADJUSTED
        )
        N_SQUARED <- gsw::gsw_Nsquared(
          SA = ABS_SAL,
          CT = CONS_TEMP,
          p = df$PRES_ADJUSTED
        )[[1]]  # Extract only the Nsquared component
        
        # Align N_SQUARED with the original profile
        df <- df %>%
          mutate(
            ABS_SAL = ABS_SAL,
            CONS_TEMP = CONS_TEMP,
            N_SQUARED = c(N_SQUARED, NA)  # Add NA for the last row to match the length
          )
        return(df)
      }) %>%
      ungroup()
    
    
    # Extract maximum N^2 for each cycle and store in the results data frame
    max_N2_per_cycle <- float_data %>%
      group_by(CYCLE_NUMBER) %>%
      summarize(
        Max_N2 = max(N_SQUARED, na.rm = TRUE),  # Calculate max N^2 while ignoring NA
        Latitude = first(LATITUDE),            # Take the first latitude for reference
        Longitude = first(LONGITUDE),          # Take the first longitude for reference
        Time = first(TIME) %>% as_date()                     # Take the first timestamp for reference
      ) %>%
      ungroup()
    
    # Append the results to the N2_results data frame
    N2_results <- bind_rows(
      N2_results,
      max_N2_per_cycle %>%
        mutate(WMO = as.integer(wmo))  # Add the WMO column to the results
    )
  })}

# Data Cleaning

N2_results <- N2_results %>%
  filter(!is.na(Max_N2) & Max_N2 > 0 &  is.finite(Max_N2))  # Remove NAs and non-physical values

# Calculate the lower and upper percentiles to identify outliers
lower_bound_N2 <- quantile(N2_results$Max_N2, 0.0025, na.rm = TRUE)
upper_bound_N2 <- quantile(N2_results$Max_N2, 0.9975, na.rm = TRUE)

# Create a new column 'cleaned_N2' that sets outliers to NA
N2_results <- N2_results %>%
  mutate(cleaned_max_N2 = ifelse(Max_N2 >= lower_bound_N2 & Max_N2 <= upper_bound_N2, Max_N2, NA)) %>%
  filter(!is.na(cleaned_max_N2))

N2_results$cleaned_max_N2 %>% summary()

# Print or save results

# Save to file
write.csv(N2_results, "/data/GLOBARGO/src/data/N2_results.csv", row.names = FALSE)


