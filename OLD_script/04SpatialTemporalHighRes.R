


# Load necessary libraries
library(conflicted)
# Resolve function conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

library(tidyverse)
library(robustbase)
library(gsw)
library(zoo)
library(oce)
library(ggpubr)
library(segmented)
library(dplyr)
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

# Step 1: Investigate the structure of the data
df_complete_clean <- read_csv(file = "/data/GLOBARGO/src/data/df_eddy_subduction_anom.csv")
df_complete_clean %>% head()
df_argo_clean <- read_csv(file = "/data/GLOBARGO/src/data/df_argo_loc.csv")
df_argo_clean %>% head()

df_mld <- read_csv(file = "/data/GLOBARGO/src/data/mld_results.csv")
# var renaming
df_mld$TIME <- df_mld$Time 

# Step 2: Define months for the four distinct seasons: DJF, MAM, JJA, SON
djf_months <- c(12, 1, 2)   # December, January, February
mam_months <- c(3, 4, 5)    # March, April, May
jja_months <- c(6, 7, 8)    # June, July, August
son_months <- c(9, 10, 11)  # September, October, November

# Step 3: Filter data for each season based on the month of the year
df_argo_djf <- df_argo_clean %>% filter(month(TIME) %in% djf_months)
df_argo_mam <- df_argo_clean %>% filter(month(TIME) %in% mam_months)
df_argo_jja <- df_argo_clean %>% filter(month(TIME) %in% jja_months)
df_argo_son <- df_argo_clean %>% filter(month(TIME) %in% son_months)

df_complete_djf <- df_complete_clean %>% filter(month(TIME) %in% djf_months)
df_complete_mam <- df_complete_clean %>% filter(month(TIME) %in% mam_months)
df_complete_jja <- df_complete_clean %>% filter(month(TIME) %in% jja_months)
df_complete_son <- df_complete_clean %>% filter(month(TIME) %in% son_months)

df_mld_djf <- df_mld %>% filter(month(TIME) %in% djf_months)
df_mld_mam <- df_mld %>% filter(month(TIME) %in% mam_months)
df_mld_jja <- df_mld %>% filter(month(TIME) %in% jja_months)
df_mld_son <- df_mld %>% filter(month(TIME) %in% son_months)


# After discussing with Leon, maybe we don't need to bin the data and it would be a classical case of logistic regression

# Winter df without anomaly 
df_argo_djf$Anomaly <- 0
# Winter df with anomaly
df_complete_djf$Anomaly <- 1
# Bind the two datasets

df_full_djf <- bind_rows(df_argo_djf,df_complete_djf,) %>%
  select(TIME,LATITUDE,LONGITUDE,CYCLE_NUMBER,WMO,Anomaly)

df_full_djf <- df_full_djf %>%
  left_join(
    df_mld_djf %>% select(WMO, CYCLE_NUMBER, TIME, cleaned_mld, binned_mld,MLD_rate_of_change, MLD_sign_rate_of_change),
    by = c("WMO", "CYCLE_NUMBER", "TIME")
  )

# Fit a binomial GAM using latitude and longitude as continuous predictors
gam_model_djf <- gam(
  Anomaly ~ s(LATITUDE, LONGITUDE, bs = "sos", k = 180),  # Spline on sphere basis 
  family = binomial(link = "logit"),
  data = df_full_djf,
  method = "REML"
)

df_full_djf %>% summary()

gam_full_model_djf <- gam(
  Anomaly ~ s(LATITUDE, LONGITUDE, bs = "sos", k = 180) +
    s(as.numeric(TIME), k = 30) +  # Smooth spline for TIME
    MLD_sign_rate_of_change +  MLD_rate_of_change + cleaned_mld + binned_mld,  # Tensor interaction for space-time
  family = binomial(link = "logit"),
  data = df_full_djf,
  method = "REML"
)




# Display the model summary

summary(gam_model_djf)

# Generate predictions on a grid of latitude and longitude
prediction_grid_djf <- expand.grid(
  LATITUDE = seq(min(df_full_djf$LATITUDE), max(df_full_djf$LATITUDE), length = 2*180),
  LONGITUDE = seq(min(df_full_djf$LONGITUDE), max(df_full_djf$LONGITUDE), length = 2*360)
)


# Predict probabilities using the fitted GAM model
prediction_grid_djf$Anomaly_Prob <- predict(gam_model_djf, newdata = prediction_grid_djf, type = "response")

# Define a binned color scale for visualization
binned_color_scale_gam <- scale_fill_viridis_b(
  name = "Anomaly Probability",
  breaks = seq(0, 1, by = 0.05),
  limits = c(0, 0.5),
  oob = scales::squish
)

# Plot the predicted probabilities on a map
ggplot() +
  geom_tile(data = prediction_grid_djf, aes(x = LONGITUDE, y = LATITUDE, fill = Anomaly_Prob)) +
  geom_contour(data = prediction_grid_djf, aes(x = LONGITUDE, y = LATITUDE, z = Anomaly_Prob), color = "white", alpha = 0.3) +
  geom_sf(data = world, fill = "gray80", color = "gray80") +
  coord_sf(xlim = range(prediction_grid_djf$LONGITUDE), ylim = range(prediction_grid_djf$LATITUDE), expand = FALSE) +
  binned_color_scale_gam +
  labs(title = "Estimated Probability of Subduction Events (DJF)",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(legend.position = "bottom")


# AUC-ROC, precision-recall curves

