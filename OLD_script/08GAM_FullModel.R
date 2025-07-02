# =====================================================
# Script for Systematic Model Optimization with Error Handling
# =====================================================

# -----------------------------
# 1. Load Necessary Libraries
# -----------------------------
library(tidyverse)
library(lubridate)
library(mgcv)
library(sf)
library(ggplot2)
library(viridis)
library(caret)
library(parallel)

# For world map
library(rnaturalearth)
library(rnaturalearthdata)

library(readr)
model_results_temp <- read_csv("model_results_temp.csv")

model_results_temp %>% View()

## 1 Output Storage
# Check if temporary file exists
if (file.exists("model_results_temp.csv")) {
  model_results <- read_csv("model_results_temp.csv")
  model_counter <- max(model_results$Model_ID) + 1  # Resume from the last model ID
  message("Resuming from previously saved progress. Starting from Model_ID = ", model_counter)
} else {
  # Initialize an empty data frame if no progress exists
  model_results <- data.frame(
    Model_ID = integer(),
    k_space = integer(),
    k_time = integer(),
    k_latitude = integer(),
    interaction = logical(),
    linear_mld = logical(),
    binned = logical(),
    AIC = numeric(),
    stringsAsFactors = FALSE
  )
  model_counter <- 1
}



# -----------------------------
# 2. Data Preparation
# -----------------------------

# Load the data (replace the file path with your actual data path)
df_full <- read_csv("dataframe_full.csv")

# View the structure of the data
glimpse(df_full)

# Remove rows with missing values in key variables
df_full <- df_full %>%
  filter(
    !is.na(LATITUDE),
    !is.na(LONGITUDE),
    !is.na(Anomaly),
    !is.na(cleaned_mld)
  )

# Create Hemisphere variable
df_full$Hemisphere <- factor(ifelse(df_full$LATITUDE >= 0, "Northern", "Southern"))

# Extract date components
df_full$Year <- year(df_full$TIME)
df_full$Month <- month(df_full$TIME)
df_full$DayOfYear <- yday(df_full$TIME)

# Adjust Day of Year for Southern Hemisphere to align seasons
df_full$AdjustedDayOfYear <- ifelse(
  df_full$Hemisphere == "Southern",
  (df_full$DayOfYear + 182.5) %% 365,
  df_full$DayOfYear
)
# Ensure that day 0 is 365
df_full$AdjustedDayOfYear[df_full$AdjustedDayOfYear == 0] <- 365

# Create a log-transformed MLD variable
df_full$log_cleaned_mld <- log(df_full$cleaned_mld + 1)  # Add 1 to avoid log(0)

# -----------------------------
# 3. Binning Data into 0.5-Degree Grid Cells
# -----------------------------

# Create binning variables
df_full <- df_full %>%
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

# -----------------------------
# 4. Hyperparameter Optimization
# -----------------------------

# Define hyperparameters to test
k_space_values <- c(200, 300, 400)       # For spatial smooth
k_time_values <- c(8, 10, 12)            # For temporal smooth
k_latitude_values <- c(40, 50, 60)       # For interaction term
include_interaction <- c(FALSE, TRUE)
include_linear_mld <- c(FALSE, TRUE)     # Include or exclude log_cleaned_mld as linear term
use_binning <- c(FALSE, TRUE)            # Whether to use binned data


# For efficient computation with large data
options(mc.cores = parallel::detectCores())

# Ensure the models directory exists
if (!dir.exists("models")) {
  dir.create("models")
}
# -----------------------------
# 4. Hyperparameter Optimization
# -----------------------------

# Define hyperparameters to test
k_space_values <- c(200, 300, 400)       # For spatial smooth
k_time_values <- c(8, 10, 12)            # For temporal smooth
k_latitude_values <- c(40, 50, 60)       # For interaction term
include_linear_mld <- c(FALSE, TRUE)     # Include or exclude log_cleaned_mld as linear term
use_binning <- c(FALSE, TRUE)            # Whether to use binned data

# Set interaction to always be TRUE
interaction <- TRUE

# For efficient computation with large data
options(mc.cores = parallel::detectCores())

# Ensure the models directory exists
if (!dir.exists("models")) {
  dir.create("models")
}

# Loop over hyperparameters
for (binning in use_binning) {
  # Choose the dataset
  data_used <- if (binning) df_binned else df_full
  data_label <- if (binning) "Binned" else "Raw"
  
  for (k_space in k_space_values) {
    for (k_time in k_time_values) {
      for (linear_mld in include_linear_mld) {
        for (k_latitude in k_latitude_values) {
          # Construct model formula
          terms <- c(
            paste0("s(LATITUDE, LONGITUDE, bs = 'sos', k = ", k_space, ")"),
            paste0("s(AdjustedDayOfYear, bs = 'cc', k = ", k_time, ")"),
            "Hemisphere"
          )
          
          if (linear_mld) {
            terms <- c(terms, "log_cleaned_mld")
          } else {
            terms <- c(terms, "s(log_cleaned_mld)")
          }
          
          # Add interaction term
          terms <- c(
            terms,
            paste0(
              "ti(AdjustedDayOfYear, LATITUDE, bs = c('cc', 'tp'), k = c(",
              k_time, ", ", k_latitude, "))"
            )
          )
          
          formula_str <- paste("Anomaly ~", paste(terms, collapse = " + "))
          formula <- as.formula(formula_str)
          
          # Decide whether to use bam() or gam()
          use_bam <- TRUE
          model_fitted <- FALSE
          
          model_name <- paste0("model_", model_counter, "_", data_label)
          
          # Attempt to fit with bam() using discrete = TRUE
          try({
            model <- bam(
              formula = formula,
              family = binomial(link = "logit"),
              data = data_used,
              method = "fREML",
              discrete = TRUE
            )
            model_fitted <- TRUE
          }, silent = TRUE)
          
          # If bam() with discrete = TRUE fails, try without discrete
          if (!model_fitted) {
            try({
              model <- bam(
                formula = formula,
                family = binomial(link = "logit"),
                data = data_used,
                method = "fREML"
              )
              model_fitted <- TRUE
              message("Fitted model ", model_name, " using bam() without discrete = TRUE")
            }, silent = TRUE)
          }
          
          # If bam() fails, try gam()
          if (!model_fitted) {
            try({
              model <- gam(
                formula = formula,
                family = binomial(link = "logit"),
                data = data_used,
                method = "REML"
              )
              model_fitted <- TRUE
              message("Fitted model ", model_name, " using gam()")
            }, silent = TRUE)
          }
          
          # Check if model fitting was successful
          if (!model_fitted) {
            message("Failed to fit model ", model_name)
            next  # Skip to the next iteration
          }
          
          # Save the model
          saveRDS(model, file = paste0("models/", model_name, ".rds"))
          
          # Append the current model's results to the data frame
          model_results <- rbind(
            model_results,
            data.frame(
              Model_ID = model_counter,
              k_space = k_space,
              k_time = k_time,
              k_latitude = k_latitude,
              interaction = interaction,
              linear_mld = linear_mld,
              binned = binning,
              AIC = AIC(model),
              stringsAsFactors = FALSE
            )
          )
          
          # Save the updated model_results incrementally
          write_csv(model_results, "model_results_temp.csv")
          
          # Remove temporary variables
          rm(model, model_fitted)
          
          # Update counter
          model_counter <- model_counter + 1
          
          # Print progress
          message("Fitted model ", model_name, " using ", data_label, " data")
        }
      }
    }
  }
}

# -----------------------------
# 5. Model Selection
# -----------------------------

# Find the model with the lowest AIC
best_model_info <- model_results %>% arrange(AIC) %>% slice(1)
best_model_id <- best_model_info$Model_ID
best_model_binned <- best_model_info$binned
best_model_label <- ifelse(best_model_binned, "Binned", "Raw")
best_model_file <- paste0("models/model_", best_model_id, "_", best_model_label, ".rds")
best_model <- readRDS(best_model_file)

# Print best model parameters
print("Best Model Parameters:")
print(best_model_info)

# -----------------------------
# 6. Model Evaluation
# -----------------------------

# Summary of the best model
summary(best_model)

# Plot the smooth terms
par(mfrow = c(2, 2))
plot(best_model, pages = 1, shade = TRUE)

# Check residuals
gam.check(best_model)

# -----------------------------
# 7. Cross-Validation
# -----------------------------

# Due to large dataset, we'll perform cross-validation on a sample
set.seed(123)
if (best_model_binned) {
  cv_sample <- df_binned %>% sample_frac(0.1)  # 10% sample
} else {
  cv_sample <- df_full %>% sample_frac(0.1)  # 10% sample
}

# Define cross-validation folds
train_control <- trainControl(method = "cv", number = 5)

# Create a model formula from the best model parameters
terms <- c(
  paste0("s(LATITUDE, LONGITUDE, bs = 'sos', k = ", best_model_info$k_space, ")"),
  paste0("s(AdjustedDayOfYear, bs = 'cc', k = ", best_model_info$k_time, ")"),
  "Hemisphere"
)

if (best_model_info$linear_mld) {
  terms <- c(terms, "log_cleaned_mld")
} else {
  terms <- c(terms, "s(log_cleaned_mld)")
}

if (best_model_info$interaction) {
  terms <- c(
    terms,
    paste0("ti(AdjustedDayOfYear, LATITUDE, bs = c('cc', 'tp'), k = c(", best_model_info$k_time, ", ", best_model_info$k_latitude, "))")
  )
}

cv_formula_str <- paste("Anomaly ~", paste(terms, collapse = " + "))
cv_formula <- as.formula(cv_formula_str)

# Fit the model using caret
cv_model <- train(
  formula = cv_formula,
  data = cv_sample,
  method = "gam",
  family = binomial(link = "logit"),
  trControl = train_control
)

print(cv_model)

# -----------------------------
# 8. Save Final Model and Results
# -----------------------------

# Save the best model
saveRDS(best_model, file = "best_model.rds")

# Save the model results dataframe
write_csv(model_results, "model_results.csv")

# -----------------------------
# 9. Interpretation and Visualization
# -----------------------------

# Predict on a grid for visualization
lat_seq <- seq(min(df_full$LATITUDE), max(df_full$LATITUDE), length.out = 100)
lon_seq <- seq(min(df_full$LONGITUDE), max(df_full$LONGITUDE), length.out = 100)
time_seq <- seq(1, 365, length.out = 12)  # Monthly intervals

prediction_grid <- expand.grid(
  LATITUDE = lat_seq,
  LONGITUDE = lon_seq,
  AdjustedDayOfYear = time_seq,
  Hemisphere = c("Northern", "Southern")
)

# Use mean of log_cleaned_mld
prediction_grid$log_cleaned_mld <- mean(df_full$log_cleaned_mld, na.rm = TRUE)

# Predict using the best model
prediction_grid$PredictedAnomaly <- predict(
  best_model,
  newdata = prediction_grid,
  type = "response"
)

# Plot the spatial distribution of predicted anomalies for a given time
# For example, plot for AdjustedDayOfYear = 183 (approx mid-year)
plot_data <- prediction_grid %>%
  filter(AdjustedDayOfYear == 183)

ggplot() +
  geom_tile(data = plot_data, aes(x = LONGITUDE, y = LATITUDE, fill = PredictedAnomaly)) +
  geom_sf(data = world, fill = "gray80", color = "white") +
  scale_fill_viridis(name = "Predicted Probability", option = "viridis") +
  labs(title = "Predicted Probability of Anomaly (Mid-Year)") +
  theme_minimal()

# -----------------------------
# End of Script
# -----------------------------
