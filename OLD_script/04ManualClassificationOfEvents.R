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

# Resolve function conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Load argo dataset
df_argo <- read_csv("/data/GLOBARGO/data/argo_profiles_df.csv")


# Load datasets of manually verified detected subduction events in salinity & spiciness :
df_abs_sal <- read_csv("/data/GLOBARGO/data/detected_events_abs_sal_var_v5.csv")
df_spic <- read_csv("/data/GLOBARGO/data/detected_events_sens_and_spec_incr.csv")

df_partial_salinity_class %>% filter(WMO == "1902304_plot")
# Classified subduction (WITH or WITHOUT carbon) datasets
# Spiciness 
df_spic_class <- read_csv("/data/GLOBARGO/data/anom_in_spic_and_sal_cat1_and2.csv")
# Salinity
df_partial_salinity_class <- read_csv("/data/GLOBARGO/data/classification_results_salinity_fig_not_in_spic.csv")
# Deep anomalies
df_deep <- read_csv("/data/GLOBARGO/data/deep_classification.csv")

df_deep$WMO <- gsub("_plot", "", df_deep$WMO)

df_deep$WMO <- df_deep$WMO %>% as.numeric()

df_deep$CYCLE_NUMBER <- df_deep$Cycle


df_partial_salinity_class <- df_partial_salinity_class %>%
  filter(Category %in% c(1, 2)) %>%
  mutate(WMO = gsub("_plot", "", WMO))  # Remove the 'plot' suffix

df_partial_salinity_class_wcat3 <- df_partial_salinity_class %>%
  filter(Category %in% c(1, 2,3)) %>%
  mutate(WMO = gsub("_plot", "", WMO))  # Remove the 'plot' suffix


df_partial_salinity_class$CYCLE_NUMBER <- df_partial_salinity_class$Cycle
df_partial_salinity_class_wcat3$CYCLE_NUMBER <- df_partial_salinity_class_wcat3$Cycle
df_partial_salinity_class_wcat3$CYCLE_NUMBER <- as.numeric(df_partial_salinity_class_wcat3$CYCLE_NUMBER)

df_partial_salinity_class$WMO <- as.numeric(df_partial_salinity_class$WMO)
df_partial_salinity_class_wcat3$WMO <- as.numeric(df_partial_salinity_class_wcat3$WMO)

df_spic_class$WMO             <- as.numeric(df_spic_class$WMO)    



# Perform the join, 
df_combined <- bind_rows(df_partial_salinity_class, df_spic_class) %>% select(WMO,CYCLE_NUMBER,Category)

df_combined_wcat3 <- bind_rows(df_partial_salinity_class_wcat3, df_spic_class) %>% select(WMO,CYCLE_NUMBER,Category)


# Join df_combined with df_deep, assuming that if a row is common to spic and partial_salinity in terms of their unique
# combination of WMO and Cycle, df_deep should be prefered :

# Remove rows from df_combined that have the same WMO and CYCLE_NUMBER as in df_deep
df_combined_filtered <- df_combined %>%
  anti_join(df_deep, by = c("WMO", "CYCLE_NUMBER"))

df_combined_filtered_wcat3 <- df_combined_wcat3 %>%
  anti_join(df_deep, by = c("WMO", "CYCLE_NUMBER"))


# Combine the filtered df_combined with df_deep
df_final <- bind_rows(df_deep, df_combined_filtered)

df_final_cat3 <- bind_rows(df_deep, df_combined_filtered_wcat3)

# Next we wish to complete this dataframe by retriving
# LON/LAT/TIME info for each of those individual subduction events, df_abs_sal has the required info


# Check for duplicates in df_abs_sal based on WMO and CYCLE_NUMBER
df_abs_sal_unique <- df_abs_sal %>%
  distinct(WMO, CYCLE_NUMBER, .keep_all = TRUE)


# Perform the join again with the unique rows
df_complete <- df_final %>%
  left_join(df_abs_sal_unique %>% 
              select(PRES_ADJUSTED,WMO, CYCLE_NUMBER, LATITUDE, LONGITUDE, TIME), 
            by = c("WMO", "CYCLE_NUMBER"))

# Clean subduction events data
df_complete_clean <- df_complete %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE) & is.finite(LONGITUDE) & is.finite(LATITUDE))

df_argo_clean <- df_argo %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE) & is.finite(LONGITUDE) & is.finite(LATITUDE))


write_csv(df_complete_clean,file = "/data/GLOBARGO/src/data/df_eddy_subduction_anom.csv")
write_csv(df_argo_clean, file = "/data/GLOBARGO/src/data/df_argo_loc.csv")

df_complete_clean <- read_csv(file = "/data/GLOBARGO/src/data/df_eddy_subduction_anom.csv")



# Load datasets of manually verified detected subduction events with carbon  
# Load datasets for carbon classification
df_carbon_cat1 <- read_csv("data/carbon_cat1_class.csv")
df_carbon_cat2 <- read_csv("data/carbon_cat2_class.csv")

# Remove "_plot" suffix and ensure WMO is numeric
df_carbon_cat1 <- df_carbon_cat1 %>%
  mutate(WMO = as.numeric(gsub("_plot", "", WMO)),
         CYCLE_NUMBER = Cycle)

df_carbon_cat2 <- df_carbon_cat2 %>%
  mutate(WMO = as.numeric(gsub("_plot", "", WMO)),
         CYCLE_NUMBER = Cycle)

# Combine carbon datasets
df_carbon_final <- bind_rows(df_carbon_cat1, df_carbon_cat2) %>%
  select(WMO, CYCLE_NUMBER, Category)


# Add location and time data from df_abs_sal
df_carbon_complete <- df_carbon_final %>%
  left_join(df_abs_sal_unique %>% select(WMO,PRES_ADJUSTED, CYCLE_NUMBER, LATITUDE, LONGITUDE, TIME), 
            by = c("WMO", "CYCLE_NUMBER"))

df_carbon_complete_clean <- df_carbon_complete %>%
  filter(!is.na(LONGITUDE) & !is.na(LATITUDE) & is.finite(LONGITUDE) & is.finite(LATITUDE))

# Save the cleaned carbon data
write_csv(df_carbon_complete_clean, file = "data/df_carbon_subduction_anom.csv")
