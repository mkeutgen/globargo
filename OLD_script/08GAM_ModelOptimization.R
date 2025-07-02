# =====================================================
# Script for Modeling Subduction Events in Argo Profiles
# =====================================================

# -----------------------------
# 1. Load Necessary Libraries
# -----------------------------
library(tidyverse)
library(lubridate)
library(mgcv)
library(spdep)
library(sf)
library(sp)
library(spatialreg)
library(ggplot2)
library(ggpubr)
library(viridis)
library(caret)
library(gstat)
library(spaMM)

library(rnaturalearth)
library(rnaturalearthdata)

# -----------------------------
# 2. Data Preparation
# -----------------------------

# Load the data (replace the file path with your actual data path)
df_full_djf <- read_csv("dataframe_full_djf.csv")

# View the structure of the data
glimpse(df_full_djf)

# Since the dataset is large (~33,000 rows), we'll work with a random subset to make computations manageable
set.seed(123)  # For reproducibility
sample_fraction <- 1  # Adjust the fraction as needed
df_sample <- df_full_djf %>% sample_frac(sample_fraction)

# Verify the sample
nrow(df_sample)

# -----------------------------
# 3. Exploratory Data Analysis
# -----------------------------

# 3.1 Check the distribution of the response variable
table(df_sample$Anomaly)

# 3.2 Plot the spatial distribution of anomalies
world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data = df_sample) +
  geom_sf(data = world, fill = "gray90", color = "gray50") +
  geom_point(aes(x = LONGITUDE, y = LATITUDE, color = factor(Anomaly)), alpha = 0.1) +
  scale_color_manual(values = c("blue", "red"), labels = c("No Anomaly", "Anomaly")) +
  labs(title = "Spatial Distribution of Anomalies",
       color = "Anomaly") +
  theme_minimal()

# Observation:
# - Hotspots are visible in certain regions.
# - There are areas with no anomalies, possibly due to low sampling effort.

# -----------------------------
# 4. Assessing Spatial Autocorrelation
# -----------------------------

# 4.1 Create Spatial Weights Matrix for Moran's I Test
coordinates <- df_sample %>% select(LONGITUDE, LATITUDE)
coordinates_sf <- st_as_sf(coordinates, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)

# Create neighbors list using k-nearest neighbors
coords_matrix <- st_coordinates(coordinates_sf)
knn_neighbors <- knearneigh(coords_matrix, k = 8)
nb <- knn2nb(knn_neighbors)
listw <- nb2listw(nb, style = "W")

# 4.2 Moran's I Test on Anomaly Variable
anomaly_factor <- as.numeric(df_sample$Anomaly)
moran_test <- moran.test(anomaly_factor, listw)
print(moran_test)

# Interpretation:
# - Check the p-value to determine if spatial autocorrelation is significant.

# 4.3 Moran's I Test on Model Residuals (after fitting a basic model)
# Fit a basic logistic regression model
basic_model <- glm(Anomaly ~ 1, family = binomial(link = "logit"), data = df_sample)

# Extract residuals
residuals_basic <- residuals(basic_model, type = "pearson")

# Moran's I Test on residuals
moran_residuals <- moran.test(residuals_basic, listw)
print(moran_residuals)

# Interpretation:
# - If significant, consider models that account for spatial autocorrelation.

# -----------------------------
# 5. Model Fitting
# -----------------------------

# 5.1 Logistic Regression Model
logistic_model <- glm(Anomaly ~ cleaned_mld, family = binomial(link = "logit"), data = df_sample)
summary(logistic_model)

# 5.2 Generalized Additive Model (GAM) without Spatial Correlation
gam_model <- gam(
  Anomaly ~ s(LATITUDE, LONGITUDE, bs = "sos", k = 50),
  family = binomial(link = "logit"),
  data = df_sample,
  method = "REML"
)
summary(gam_model)

# Check GAM residuals for spatial autocorrelation
residuals_gam <- residuals(gam_model, type = "pearson")
moran_residuals_gam <- moran.test(residuals_gam, listw)
print(moran_residuals_gam)

# 5.3 GAM with Additional Predictors
gam_model_extended <- gam(
  Anomaly ~ s(LATITUDE, LONGITUDE, bs = "sos", k = 200) +
    s(cleaned_mld) +
    s(MLD_rate_of_change),
  family = binomial(link = "logit"),
  data = df_sample,
  method = "REML"
)
summary(gam_model_extended)
# R-sq.(adj) =  0.0887   Deviance explained = 20.9%

gam_model_extended_smallk <- gam(
  Anomaly ~ s(LATITUDE, LONGITUDE, bs = "sos", k = 50) +
    s(cleaned_mld) +
    s(MLD_rate_of_change),
  family = binomial(link = "logit"),
  data = df_sample,
  method = "REML"
)
summary(gam_model_extended_largek)


# 5.4 GAM using bam() for Large Data
bam_model <- bam(
  Anomaly ~ s(LATITUDE, LONGITUDE, bs = "sos", k = 200),
  family = binomial(link = "logit"),
  data = df_full_djf,
  method = "fREML",
  discrete = TRUE
)
summary(bam_model)

# 5.5 GAM with Spatial Autocorrelation using gamm()
# Note: This may still be computationally intensive
# Try reducing k or use a smaller sample if necessary
gamm_model <- gamm(
  Anomaly ~ s(LATITUDE, LONGITUDE, bs = "sos", k = 10),
  family = binomial(link = "logit"),
  data = df_sample,
  correlation = corExp(form = ~ LONGITUDE + LATITUDE),
  method = "REML"
)
summary(gamm_model$gam)

# 5.6 Spatial Generalized Linear Mixed Model using spaMM
# This approach can handle spatial random effects more efficiently
spaMM_model <- fitme(
  Anomaly ~ 1 + Matern(1 | LONGITUDE + LATITUDE),
  family = binomial(),
  data = df_sample
)
summary(spaMM_model)

# -----------------------------
# 6. Model Comparison and Validation
# -----------------------------

# 6.1 Compare AIC Values
AIC(logistic_model, gam_model, gam_model_extended, gamm_model$gam)

# 6.2 Cross-Validation using caret
# Define training control
train_control <- trainControl(method = "cv", number = 5)

# Train GAM model using cross-validation
set.seed(123)
cv_model <- train(
  Anomaly ~ s(LATITUDE, LONGITUDE, bs = "tp", k = 50),
  data = df_sample,
  method = "gam",
  family = binomial(link = "logit"),
  trControl = train_control
)
print(cv_model)

# 6.3 Residual Diagnostics
# For GAM model
gam.check(gam_model)

# QQ-plot for residuals
qqnorm(residuals(gam_model))
qqline(residuals(gam_model))

# 6.4 ROC Curve and AUC
library(pROC)

# Predict probabilities
df_sample$predicted_prob <- predict(gam_model, type = "response")

# ROC curve
roc_obj <- roc(df_sample$Anomaly, df_sample$predicted_prob)
plot(roc_obj, print.auc = TRUE)

# -----------------------------
# 7. Handling Uneven Sampling and Hotspots
# -----------------------------

# 7.1 Weighting by Sampling Effort
# Create a variable representing sampling effort (e.g., number of observations in a grid cell)

# Define grid size for spatial binning (e.g., 2-degree grid cells)
grid_size <- 0.1

# Bin data into grid cells
df_sample <- df_sample %>%
  mutate(
    LAT_BIN = floor(LATITUDE / grid_size) * grid_size,
    LON_BIN = floor(LONGITUDE / grid_size) * grid_size
  )

# Calculate sampling effort per grid cell
sampling_effort <- df_sample %>%
  group_by(LAT_BIN, LON_BIN) %>%
  summarize(Effort = n())

# Merge sampling effort back into the main data
df_sample <- df_sample %>%
  left_join(sampling_effort, by = c("LAT_BIN", "LON_BIN"))

# 7.2 Include Weights in the Model
gam_model_weighted <- gam(
  Anomaly ~ s(LATITUDE, LONGITUDE, bs = "tp", k = 50),
  family = binomial(link = "logit"),
  data = df_sample,
  weights = Effort,
  method = "REML"
)
summary(gam_model_weighted)

# 7.3 Addressing Zero-Inflation
# If the data has many zeros, consider zero-inflated models

# Install package for zero-inflated models
# install.packages("pscl")
library(pscl)

# Fit a zero-inflated model
# Note: Zero-inflated models are more common for count data, but can be adapted if necessary

# -----------------------------
# 8. Final Model Selection
# -----------------------------

# Based on AIC, cross-validation, and residual diagnostics, choose the best model
# Suppose gam_model_extended performed best

# Refit the chosen model on the full dataset if computationally feasible
best_model <- gam_model_extended  # Replace with the chosen model

# -----------------------------
# 9. Predictions and Visualization
# -----------------------------

# 9.1 Create Prediction Grid
lat_range <- seq(min(df_full_djf$LATITUDE), max(df_full_djf$LATITUDE), length.out = 100)
lon_range <- seq(min(df_full_djf$LONGITUDE), max(df_full_djf$LONGITUDE), length.out = 100)
prediction_grid <- expand.grid(LATITUDE = lat_range, LONGITUDE = lon_range)

# Predict using the best model
prediction_grid$predicted_anomaly <- predict(best_model, newdata = prediction_grid, type = "response")

# 9.2 Plot Predicted Probability Surface
ggplot() +
  geom_tile(data = prediction_grid, aes(x = LONGITUDE, y = LATITUDE, fill = predicted_anomaly)) +
  geom_sf(data = world, fill = "gray90", color = "gray50") +
  scale_fill_viridis(name = "Predicted Probability", option = "viridis") +
  labs(title = "Predicted Probability of Subduction Anomalies") +
  theme_minimal()

# -----------------------------
# 10. Saving the Model and Results
# -----------------------------

# Save the model object
saveRDS(best_model, "best_model.rds")

# Save the predictions
write_csv(prediction_grid, "prediction_grid.csv")

# -----------------------------
# 11. Documentation and Reporting
# -----------------------------

# Document each step, assumptions made, and interpretations of results in a report or notebook.

# -----------------------------
# End of Script
# -----------------------------
