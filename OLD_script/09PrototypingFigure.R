# 9 PROTOTYPING FIGURE 
library(tidyverse)
library(robustbase)
library(gsw)
library(zoo)
library(oce)


# Let WMO : 
wmo = 5904105


# Find right times : 
# Cycle 20 : 2013-06-09 
# Cycle 23 : "2013-06-24"  
# Cycle 57 : "2013-12-16"
# Cycle 83 : "2014-06-15"
# Cycle 99 : "2014-10-08"
# Cycle 108 : "2014-12-12"
# Cycle 142 : "2015-08-13"

# We need spiciness, AOU and bbp700

data_df <- load_float_data(float_ids = wmo,
                           variables = c("DATA_TYPE", "PLATFORM_NUMBER", "BBP700", "BBP700_dPRES",
                                         "BBP700_ADJUSTED_QC", "LATITUDE", "LONGITUDE", "PROFILE_TEMP_QC",
                                         "PROFILE_DOXY_QC", "PROFILE_BBP700_QC", "PRES_QC", "PRES",
                                         "PRES_ADJUSTED", "PROFILE_PSAL_QC", "CHLA_QC", "CHLA_ADJUSTED",
                                         "CHLA_ADJUSTED_ERROR", "DOXY", "DOXY_QC", "DOXY_ADJUSTED",
                                         "DOXY_ADJUSTED_QC", "DOXY_ADJUSTED_ERROR", "PSAL", "PSAL_dPRES",
                                         "PSAL_ADJUSTED", "PSAL_ADJUSTED_QC", "TEMP", "TEMP_QC", "TEMP_dPRES",
                                         "TEMP_ADJUSTED", "TEMP_ADJUSTED_QC", "TEMP_ADJUSTED_ERROR"),
                           format = "dataframe")


data_df <- data_df %>% filter(!is.na(DOXY)) %>% group_by(CYCLE_NUMBER) %>%
  mutate(SPIC = swSpice(salinity = PSAL_ADJUSTED, temperature = TEMP_ADJUSTED,
                        latitude = first(LATITUDE), longitude = first(LONGITUDE), eos = "unesco"),
         ABS_SAL = gsw::gsw_SA_from_SP(SP = PSAL_ADJUSTED, p = PRES_ADJUSTED,
                                       longitude = first(LONGITUDE), latitude = first(LATITUDE)),
         CONS_TEMP = gsw::gsw_CT_from_t(SA = ABS_SAL, t = TEMP_ADJUSTED, p = PRES_ADJUSTED),
         SAT_DOXY = gsw_O2sol(SA = ABS_SAL, CT = CONS_TEMP, p = PRES_ADJUSTED,
                              longitude = first(LONGITUDE), latitude = first(LATITUDE)),
         AOU = SAT_DOXY - DOXY_ADJUSTED,
         SIGMA0 = gsw::gsw_sigma0(SA = ABS_SAL, CT = CONS_TEMP))

# Moving to python
write_csv(x = data_df,"/data/GLOBARGO/data/wmo_5904105_df.csv")

# We need to make nice section plots 

data <- data_df %>% ungroup %>% select(PRES_ADJUSTED,TIME,DOXY) 

data <- data %>%
  mutate(TIME = as.Date(TIME))

# Create the time-depth plot using ggplot2
ggplot(filtered_data, aes(x = TIME, y = PRES_ADJUSTED, fill = DOXY)) +
  geom_tile(interpolate = TRUE) +  # Use geom_raster to plot, interpolate makes it smoother
  scale_y_reverse() +  # Reverse y-axis to have the surface at the top
  scale_fill_viridis_c(option = "C", direction = -1) +  # Use a viridis color scale similar to the one in the image
  labs(x = "Time", y = "Depth [m]", fill = "DOXY [µmol/kg]", title = "Time-Depth Plot of Dissolved Oxygen (DOXY)") +
  theme_minimal(base_size = 15) +  # Elegant theme with larger base font size
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),  # Remove grid for cleaner look
    legend.position = "right"
  )

filtered_data <- data %>%
  filter(TIME >= as.Date("2013-04-01") & TIME <= as.Date("2013-10-16"))

filtered_data$TIME %>% unique()

# Create the time-depth plot using ggplot2
ggplot(filtered_data, aes(x = TIME, y = PRES_ADJUSTED, color = DOXY)) +
  geom_point(size = 4, alpha = 0.6) +  # Use geom_point to handle irregular data
  scale_y_reverse() +  # Reverse y-axis to have the surface at the top
  scale_color_viridis_c(option = "C", direction = -1) +  # Use a viridis color scale similar to the one in the image
  labs(x = "Time", y = "Depth [m]", color = "DOXY [µmol/kg]", title = "Time-Depth Plot of Dissolved Oxygen (DOXY)") +
  theme_minimal(base_size = 15) +  # Elegant theme with larger base font size
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),  # Remove grid for cleaner look
    legend.position = "right"
  )
filtered_data


library(ggplot2)
library(dplyr)
library(akima)

# Convert TIME to Date type
data <- data %>%
  mutate(TIME = as.Date(TIME))

data$TIME %>% unique()

# Create the time-depth plot using stat_summary_2d
ggplot(data, aes(x = TIME, y = PRES_ADJUSTED, z = DOXY)) +
  stat_summary_2d(fun = mean, bins = 100) +  # Compute the mean DOXY in 2D bins
  scale_y_reverse() +  # Reverse y-axis to have the surface at the top
  scale_fill_viridis_c(option = "C", direction = -1) +  # Use a viridis color scale similar to the one in the image
  labs(x = "Time", y = "Depth [m]", fill = "Mean DOXY [µmol/kg]", title = "Time-Depth Plot of Dissolved Oxygen (DOXY)") +
  theme_minimal(base_size = 15) +  # Elegant theme with larger base font size
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),  # Remove grid for cleaner look
    legend.position = "right"
  )
