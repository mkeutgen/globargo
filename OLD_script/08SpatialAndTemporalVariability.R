#CYCLE_NUMBER = 77
library(readr)
library(tidyverse)
library(lubridate)
classification_results_v4 <- read_csv("/data/GLOBARGO/data/classification_results_v4.csv")
classification_results_v4 <- classification_results_v4 %>%  mutate(WMO = str_replace(WMO, "_plot$", ""))
classification_results_v4$CYCLE_NUMBER <- classification_results_v4$Cycle
classification_results_v4$WMO <- as.numeric(classification_results_v4$WMO)
argo_profiles <- read_csv("/data/GLOBARGO/data/argo_profiles_df.csv")

detected.events.df <- read_csv("/data/GLOBARGO/data/detected_events_sens_and_spec_incr.csv")
detected.events.df %>% select(LATITUDE,LONGITUDE,WMO,CYCLE_NUMBER) %>% unique()

merged_df <- detected.events.df %>%
  left_join(classification_results_v4, by = c("WMO", "CYCLE_NUMBER"))

# Remove spurrious Cycle column
merged_df <- merged_df %>% dplyr::select(-Cycle)

# Create a dataframe for the mapping of events, pick lowest value for PRES_ADJUSTED :
unique_sub_loc_df <- merged_df  %>%
  group_by(WMO, CYCLE_NUMBER) %>%
  slice_min(PRES_ADJUSTED) %>%
  ungroup() %>%
  select(CYCLE_NUMBER, LATITUDE, LONGITUDE, WMO, Category, TIME, PRES_ADJUSTED,CONSISTENT_ANOM)


unique_sub_loc_df %>% filter(Category %in% c(1,2))

unique_sub_loc_df %>% filter(Category %in% c(1,2)) %>% filter(CONSISTENT_ANOM == 1)

unique_sub_loc_df %>% filter(Category %in% c(1,2)) %>% filter(CONSISTENT_ANOM == 0)


write_csv(x = unique_sub_loc_df,file="/data/GLOBARGO/src/data/unique_sub_loc.csv")

unique_sub_loc_df %>% filter(WMO==5904677)

merged_df %>% filter(Category == 0) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_histogram(binwidth=38)+theme_bw()
merged_df %>% filter(Category == 1) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_histogram(binwidth=38)+theme_bw()
merged_df %>% filter(Category == 2) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_histogram(binwidth=38)+theme_bw()
merged_df %>% filter(Category == 3) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_histogram(binwidth=38)+theme_bw()


merged_df %>% filter(Category == 0) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_density()+theme_bw()
merged_df %>% filter(Category == 1) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_density()+theme_bw()
merged_df %>% filter(Category == 2) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_density()+theme_bw()
merged_df %>% filter(Category == 3) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_density()+theme_bw()


unique_sub_loc_df %>% filter(Category == 0) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_density()+theme_bw()
unique_sub_loc_df %>% filter(Category == 1) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_density()+theme_bw()
unique_sub_loc_df %>% filter(Category == 2) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_density()+theme_bw()
unique_sub_loc_df %>% filter(Category == 3) %>% ggplot(aes(x=PRES_ADJUSTED)) + geom_density()+theme_bw()


# SEASONALITY, GLM Model, are there more events in winter than during other seasons ? 


cat1_2_df <- unique_sub_loc_df %>% filter(Category %in% c(1,2))


cat1_2_df <- cat1_2_df %>%
  mutate(Key = paste(WMO, CYCLE_NUMBER, sep = "_"))

argo_profiles <- argo_profiles %>%
  mutate(Key = paste(WMO, CYCLE_NUMBER, sep = "_"))

# Step 2: Use left_join to match and create the Event category
argo_profiles <- argo_profiles %>%
  left_join(cat1_2_df %>% select(Key) %>% mutate(Event = 1), by = "Key") %>%
  mutate(Event = ifelse(is.na(Event), 0, Event)) %>%
  select(-Key)  # Optionally remove the Key column if no longer needed

# Step 3: View the result
print(argo_profiles)

event_data <- argo_profiles %>%
  mutate(
    Region = case_when(
      LATITUDE > 30 ~ "NH_ET",
      LATITUDE > 0 & LATITUDE <= 30 ~ "NH_Tropics",
      LATITUDE >= -30 & LATITUDE <= 0 ~ "SH_Tropics",
      LATITUDE < -30 ~ "SH_ET"
    )
  )

# Step 2: Classify Seasons
event_data <- event_data %>%
  mutate(
    Month = month(TIME),
    Season = case_when(
      Month %in% c(12, 1, 2) ~ "DJF",
      Month %in% c(3, 4, 5) ~ "MAM",
      Month %in% c(6, 7, 8) ~ "JJA",
      Month %in% c(9, 10, 11) ~ "SON"
    )
  )

event_data <- event_data %>%
  mutate(
    Hemisphere = ifelse(LATITUDE >= 0, "NH", "SH"),
    Real_Season = case_when(
      Hemisphere == "NH" & Month %in% c(12, 1, 2) ~ "Winter",
      Hemisphere == "NH" & Month %in% c(3, 4, 5) ~ "Spring",
      Hemisphere == "NH" & Month %in% c(6, 7, 8) ~ "Summer",
      Hemisphere == "NH" & Month %in% c(9, 10, 11) ~ "Fall",
      Hemisphere == "SH" & Month %in% c(12, 1, 2) ~ "Summer",
      Hemisphere == "SH" & Month %in% c(3, 4, 5) ~ "Fall",
      Hemisphere == "SH" & Month %in% c(6, 7, 8) ~ "Winter",
      Hemisphere == "SH" & Month %in% c(9, 10, 11) ~ "Spring"
    )
  )

# Interactionist model
glm_model <- glm(Event ~ Real_Season + Region, data = event_data, family = binomial)
summary(glm_model) #AIC: 27317
coefficients <- coef(glm_model)
odds_ratios <- exp(coefficients)

odds_ratio_times_baseline <- odds_ratios * 0.03070336 
odds_ratio_times_baseline/(1+odds_ratio_times_baseline)

