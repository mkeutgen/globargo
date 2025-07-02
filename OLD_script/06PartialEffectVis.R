
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
library(gratia)
library(rnaturalearth)
library(rnaturalearthdata)

# Read  the gam
m <- readRDS("models/spatial_by_season.rds")
m_gam <- readRDS("models/MLD_LAT_LON_Season.Rds")
m_gam %>% summary()

# Extract smooth effects for log(cleaned_mld)
smooth_mld <- smooth_estimates(m_gam, smooth = "s(log(cleaned_mld))")
# Extract smooth effects for log(cleaned_N2)
smooth_n2 <- smooth_estimates(m_gam, smooth = "s(log(cleaned_N2))")

# Back transform
smooth_mld$mld <- exp(smooth_mld$`log(cleaned_mld)`)  # If cleaned_mld was log-transformed earlier
smooth_n2$N2 <- exp(smooth_n2$`log(cleaned_N2)`)     # If cleaned_N2 was log-transformed earlier

# Transform the smooth effects for cleaned_mld to probabilities

smooth_mld <- smooth_mld %>%
  mutate(
    prob = 1 / (1 + exp(-.estimate)),        # Transform to probability
    prob_lower = 1 / (1 + exp(-(.estimate - .se))),  # Lower bound of confidence interval
    prob_upper = 1 / (1 + exp(-(.estimate + .se)))   # Upper bound of confidence interval
  )

# Transform the smooth effects for N2 to probabilities
smooth_n2 <- smooth_n2 %>%
  mutate(
    prob = 1 / (1 + exp(-.estimate)),        # Transform to probability
    prob_lower = 1 / (1 + exp(-(.estimate - .se))),  # Lower bound of confidence interval
    prob_upper = 1 / (1 + exp(-(.estimate + .se)))   # Upper bound of confidence interval
  )


mld_effect <- ggplot(smooth_mld, aes(x = mld, y = prob)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), alpha = 0.2, fill = "blue") +
  labs(
    title = "Effect of MLD on Subduction Probability",
    x = "Mixed Layer Depth (MLD)",
    y = "Predicted Probability"
  ) +
  theme_minimal()+xlim(0,500)


n2_effect <- ggplot(smooth_n2, aes(x = N2, y = prob)) +
  geom_line(color = "red", size = 1) +
  geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), alpha = 0.2, fill = "red") +
  labs(
    title = "Effect of Stratification (N²) on Subduction Probability over Likely Range of N2",
    x = "Brunt-Väisälä Frequency (N^2) log10 scale",
    y = "Predicted Probability"
  ) +
  theme_minimal()+scale_x_log10()+xlim(0,0.01)

# Combine the two plots side by side
combined_plot <- mld_effect + n2_effect + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Extract smooth effects for latitude and longitude 
# Extract smooth effects for log(cleaned_mld)
smooth <- smooth_estimates(m_gam)
write_csv(x = smooth,file="smooth_data.csv")

# Extract unique smooth terms to confirm structure
unique_smooth_terms <- smooth$.smooth %>% unique()
print(unique_smooth_terms)

# Filter smooth effects by each seasonal term
seasons <- c("DJF", "MAM", "JJA", "SON")
seasonal_maps <- list()

for (season in seasons) {
  # Filter for the specific season
  smooth_season <- smooth %>%
    filter(.smooth == paste0("s(LATITUDE,LONGITUDE):Season", season))
  
  # Create the spatial map
  p <- ggplot(smooth_season, aes(x = LONGITUDE, y = LATITUDE, fill = .estimate)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Effect Estimate", option = "viridis") +
    labs(
      title = paste("Spatial Effect for", season),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal()
  
  # Store the plot in a list
  seasonal_maps[[season]] <- p
}

# Combine all plots into one layout
combined_plot <- wrap_plots(seasonal_maps, ncol = 2)
