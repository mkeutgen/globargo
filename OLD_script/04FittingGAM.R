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
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggspatial)
library(viridis)
library(sp)
library(spdep)
library(mgcv)
library(patchwork)
library(gratia)

# Read data to fit the gam
df_full <- read_csv("data/dataframe_full_mld_and_n2.csv") 

# Standardize the log-transformed predictors
df_full <- df_full %>%
  mutate(
    log_cleaned_mld_std = scale(log(cleaned_mld)),
    log_cleaned_N2_std = scale(log(cleaned_N2))
  )

df_full$Anomaly <- df_full$Anomaly %>% as_factor()
df_full$Anomaly_text <- ifelse(df_full$Anomaly == 1,"Yes","No")


ks.test(df_full$cleaned_mld[df_full$Anomaly ==1],
        df_full$cleaned_mld[df_full$Anomaly ==0],alternative="less")

ks.test(df_full$Max_N2[df_full$Anomaly ==1],
        df_full$Max_N2[df_full$Anomaly ==0],alternative="greater")

set.seed(124)       # for reproducibility
B <- 10000          # number of bootstrap iterations
n <- nrow(df_full) # total number of observations

# Storage vectors for the KS statistics in each bootstrap iteration
ks_mld <- numeric(B)
ks_n2  <- numeric(B)

# Original (observed) KS stats for reference:
ks_mld_obs <- ks.test(
  df_full$cleaned_mld[df_full$Anomaly == 1],
  df_full$cleaned_mld[df_full$Anomaly == 0],
  alternative = "less"
)$statistic

ks_n2_obs <- ks.test(
  df_full$Max_N2[df_full$Anomaly == 1],
  df_full$Max_N2[df_full$Anomaly == 0],
  alternative = "greater"
)$statistic

# Bootstrap loop
for (b in seq_len(B)) {
  # 1. Resample row indices (with replacement)
  idx_boot <- sample(seq_len(n), size = n, replace = TRUE)
  df_boot  <- df_full[idx_boot, ]
  
  # 2. Subset MLD in sub vs. no-sub
  sub_mld    <- df_boot$cleaned_mld[df_boot$Anomaly == 1]
  no_sub_mld <- df_boot$cleaned_mld[df_boot$Anomaly == 0]
  
  # 3. Subset N2 in sub vs. no-sub
  sub_n2    <- df_boot$Max_N2[df_boot$Anomaly == 1]
  no_sub_n2 <- df_boot$Max_N2[df_boot$Anomaly == 0]
  
  # 4. Compute the KS statistics in this bootstrap sample
  #    (We'll just use the default two-sided test for the *statistic*,
  #     because we mainly care about the *value* of the KS distance.)
  ks_mld[b] <- suppressWarnings(
    ks.test(sub_mld, no_sub_mld)$statistic
  )
  ks_n2[b] <- suppressWarnings(
    ks.test(sub_n2, no_sub_n2)$statistic
  )
}

# -------------------------------------------
# 5. Compare the distributions of ks_mld vs. ks_n2
# -------------------------------------------

bootstrapped <- tibble(ks_mld, ks_n2) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Parameter",
    values_to = "KS_Distance"
  ) %>%
  mutate(
    Parameter = recode(
      Parameter,
      "ks_mld" = "Mixed-Layer Depth (MLD)",
      "ks_n2"  = "Brunt-Väisälä Frequency (N²)"
    )
  ) %>%
  ggplot() +
  geom_density(aes(x = KS_Distance, fill = Parameter), alpha = 0.7) +
  labs(
    title = "Estimated Density of Simulated KS Distances",
    x = "Bootstrapped Kolmogorov–Smirnov Distances",
    y = "Density",
    fill = NULL
  ) + scale_fill_viridis_d(begin=0.5)+
  theme_bw(base_size = 25) +
  theme(
    panel.border   = element_rect(color = "black", fill = NA, size = 0.5),
    axis.ticks     = element_line(color = "black"),
    plot.title     = element_text(face = "bold", size = 25, hjust = 0.5),
    legend.position = "bottom",
    legend.text    = element_text(size = 25)
  )# For instance, compute the fraction of times that MLD's KS < N2's KS:

ggsave("figures/bootstrapped_ks.png",bootstrapped,width = 15, height = 10, dpi = 300)

mean(ks_mld < ks_n2)

# Or look at the distribution of the difference (N2 minus MLD):
ks_diff <- ks_n2 - ks_mld
hist(ks_diff)
# p-value that the difference is > 0 (i.e., N2 has a larger KS on average):
p_value <- 1-mean(ks_diff > 0)


df_full$log
# Density plot for Mixed Layer Depth (MLD)
density_mld <- ggplot(df_full, aes(x = log_cleaned_mld, fill = factor(Anomaly_text))) +
  geom_density(alpha = 0.7, color = "black", size = 1) +
  labs(
    title = "Estimated PDF of MLD",
    x = "Log(Mixed Layer Depth)",
    y = "Density",
    fill = "Profile contains a subduction anomaly"
  ) +
  theme_bw(base_size = 25) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(face = "bold", size = 25, hjust = 0.5),
    plot.subtitle = element_text(size = 15, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 25),
    legend.text = element_text(size = 25)
  )+ scale_fill_viridis_d(begin=0.5)

# Density plot for Brunt-Väisälä Frequency (N²)
density_N2 <- ggplot(df_full, aes(x = log(cleaned_N2), fill = factor(Anomaly_text))) +
  geom_density(alpha = 0.7, color = "black", size = 1) +
  labs(
    title = "Estimated PDF of N²",
    x = "Log(Brunt–Väisälä Frequency)",
    y = "Density",
    fill = "Profile contains a subduction anomaly"
  ) +
  theme_bw(base_size = 25) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(face = "bold", size = 25, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 25),
    legend.text = element_text(size = 25),
  )+ scale_fill_viridis_d(begin=0.5)

t.test(df_full$cleaned_mld[df_full$Anomaly == 0],df_full$cleaned_mld[df_full$Anomaly == 1])
# t = -15.139, df = 4537, p-value < 2.2e-16
t.test(df_full$log_cleaned_N2[df_full$Anomaly == 0] ,df_full$log_cleaned_N2[df_full$Anomaly == 1])
# t = 1.509, df = 4652.9, p-value = 0.1314 mean not diff

# Combine the two plots side by side
combined_density <- density_mld + density_N2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("figures/density_mld_N2.png", combined_density, width = 20, height = 10, dpi = 300)

#  N2 has a long right tail

# MLD is significant (positive)
m_mld <- glm(formula = Anomaly ~ log(cleaned_mld),
         family = binomial(link = "logit"),
         data = df_full) 

m_mld %>% summary() # 0.25844 log odds for log of MLD
# Stratification is significant (negative)
m_N2 <- glm(formula = Anomaly ~ log(cleaned_N2) ,
         family = binomial(link = "logit"),
         data = df_full) 

m_N2 %>% summary() 

# Compare the 2 models
# Compare AIC values
AIC_mld <- AIC(m_mld)
AIC_N2 <- AIC(m_N2)

print(paste("AIC for MLD model:", AIC_mld))
print(paste("AIC for N2 model:", AIC_N2))
# -0.30904 for log odds of N2

# Stratification is not significant if combined with mld
m_mld_N2 <- glm(formula = Anomaly ~ log(cleaned_N2) + log(cleaned_mld) ,
         family = binomial(link = "logit"),
         data = df_full) 

m_mld_N2 %>% summary() 
# This indicates multicollinearity
m_gam <- gam(formula = Anomaly ~ s(log(cleaned_N2)) + s(log(cleaned_mld)),
             family = binomial(link = "logit"),
             data = df_full)


# Extract smooth effects for log(cleaned_mld)
smooth_mld <- smooth_estimates(m_gam, smooth = "s(log(cleaned_mld))")
# Extract smooth effects for log(cleaned_N2)
smooth_n2 <- smooth_estimates(m_gam, smooth = "s(log(cleaned_N2))")

# # Calculate the range of the smooth effect
mld_range <- max(smooth_mld$.estimate) - min(smooth_mld$.estimate)
N2_range <- max(smooth_n2$.estimate) - min(smooth_n2$.estimate)

# Print the ranges
cat("Effect size for log(cleaned_mld):", mld_range, "\n")
cat("Effect size for log(cleaned_N2):", N2_range, "\n")

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

library(ggplot2)

mld_effect <- ggplot(smooth_mld, aes(x = mld, y = prob)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), alpha = 0.2, fill = "blue") +
  labs(
    title = "Effect of MLD on Subduction Probability",
    x = "Mixed Layer Depth (m)",
    y = "Predicted Probability"
  ) +
  theme_minimal()+scale_x_log10()+ theme_bw(base_size = 25) + theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(face = "bold", size = 25, hjust = 0.5),
    plot.subtitle = element_text(size = 15, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 15)
  )+geom_hline(yintercept = 0.5,show.legend = T)+scale_y_continuous(breaks = seq(0, 1, by = 0.1))


n2_effect <- ggplot(smooth_n2, aes(x = N2, y = prob)) +
  geom_line(color = "red", size = 1) +
  geom_ribbon(aes(ymin = prob_lower, ymax = prob_upper), alpha = 0.2, fill = "red") +
  labs(
    title = "Effect of N² on Subduction Probability",
    x = "Brunt-Väisälä Frequency (N²) log10 scale (s⁻¹)",
    y = "Predicted Probability"
  ) +
  theme_minimal()+scale_x_log10(limits=c(1e-4,1e-2))+
  # Switch to a theme that shows axes and ticks
  theme_bw(base_size = 25) +scale_y_continuous(breaks = seq(0, 1, by = 0.1))+geom_hline(yintercept = 0.5)+
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(face = "bold", size = 25, hjust = 0.5),
    plot.subtitle = element_text(size = 15, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 15)
  )

# Combine the two plots side by side
combined_plot <- mld_effect + n2_effect + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Save the combined plot
ggsave("figures/GAM_combined_effects_plot.png", combined_plot, width = 20, height = 10, dpi = 300)

m_gam %>% summary
# R-sq.(adj) =  0.0103   Deviance explained = 2.53%, so clearly, not much if just using MLD and N2


# Add spatial effect
gam_model_spatial <- bam(
  Anomaly ~  s(log(cleaned_N2)) + s(log(cleaned_mld)) + s(LATITUDE, LONGITUDE, bs = "sos", k = 5),
  family = binomial(link = "logit"),
  data = df_full,
  method = "REML"
)

gam_model_spatial %>% summary()
# Deviance explained = 19.4%

gam_model_spatial_temporal <- gam(
  Anomaly ~  s(log(cleaned_N2),k=15) + s(log(cleaned_mld),k=15) + s(LATITUDE, LONGITUDE, bs = "sos", k = 600)+
    s(AdjustedDayOfYear, bs = 'cc', k = 20),
  family = binomial(link = "logit"),
  data = df_full,
  method = "REML"
)
gam_model_spatial_temporal %>% summary()

# Deviance explained = 19.5%

gam_model_spatial_temporal_int <- gam(
  Anomaly ~  s(log(cleaned_N2)) + s(log(cleaned_mld)) + s(LATITUDE, LONGITUDE, bs = "sos", k = 600)+
    s(AdjustedDayOfYear, bs = 'cc', k = 20) + ti(AdjustedDayOfYear, LATITUDE, bs = c('cc', 'tp'), k = c(20,20)),
  family = binomial(link = "logit"),
  data = df_full,
  method = "REML"
)


# Add seasonal effect :
# Create Season factor
df_full$Season <- cut(
  df_full$AdjustedDayOfYear,
  breaks = c(0, 90, 180, 270, 365),
  labels = c("Winter", "Spring", "Summer", "Fall"),
  include.lowest = TRUE
)

gam_model_spatial_temporal_int_hemi_season <- gam(
  Anomaly ~ Hemisphere + Season + s(log(cleaned_N2)) + s(log(cleaned_mld)) + s(LATITUDE, LONGITUDE, bs = "sos", k = 600)+
    s(AdjustedDayOfYear, bs = 'cc', k = 20) + ti(AdjustedDayOfYear, LATITUDE, bs = c('cc', 'tp'), k = c(20,20)),
  family = binomial(link = "logit"),
  data = df_full,
  method = "REML"
)




# Define grid size to approximately 25 km
# Define bin size
bin_size <- 0.25

# Define longitude and latitude bins
longitude_bins <- seq(floor(min(df_full$LONGITUDE, na.rm = TRUE)), ceiling(max(df_full$LONGITUDE, na.rm = TRUE)), by = bin_size)
latitude_bins <- seq(floor(min(df_full$LATITUDE, na.rm = TRUE)), ceiling(max(df_full$LATITUDE, na.rm = TRUE)), by = bin_size)

# Compute bin centers for labeling
lon_centers <- longitude_bins[-length(longitude_bins)] + bin_size / 2
lat_centers <- latitude_bins[-length(latitude_bins)] + bin_size / 2

# Assign bins to Argo profiles and anomalies
assign_bins <- function(df) {
  df %>%
    mutate(
      lon_bin = cut(LONGITUDE, breaks = longitude_bins, include.lowest = TRUE, labels = lon_centers),
      lat_bin = cut(LATITUDE, breaks = latitude_bins, include.lowest = TRUE, labels = lat_centers)
    ) %>%
    mutate(
      lon_bin = as.numeric(as.character(lon_bin)),
      lat_bin = as.numeric(as.character(lat_bin))
    )
}

df_full_binned <- assign_bins(df_full)



# Optionally, average predictors within grid cells and months
df_agg <- df_full_binned %>%
  group_by(lat_bin, lon_bin, Month) %>%
  summarise(
    Anomaly = sum(Anomaly),
    Total =  n(),
    cleaned_mld = mean(cleaned_mld),
    cleaned_N2 = mean(cleaned_N2),
    LATITUDE = mean(LATITUDE),
    LONGITUDE = mean(LONGITUDE),
  ) %>%
  ungroup()

df_agg <- df_agg %>%
  mutate(
    Season = case_when(
      Month %in% c(12, 1, 2) ~ "DJF",
      Month %in% c(3, 4, 5) ~ "MAM",
      Month %in% c(6, 7, 8) ~ "JJA",
      Month %in% c(9, 10, 11) ~ "SON"
    ),
    Season = factor(Season, levels = c("DJF", "MAM", "JJA", "SON"))
  )

gam_model <- gam(
  cbind(Anomaly, Total - Anomaly) ~
    s(log(cleaned_N2)) +
    s(log(cleaned_mld)) +
    s(LATITUDE, LONGITUDE, by = Season, bs = "sos", k = 500) +
    Season,
  family = binomial,
  data = df_agg
)
# R-sq.(adj) =  0.271   Deviance explained = 42.3%
gam_model %>% summary()

write_rds(gam_model,file = "models/MLD_LAT_LON_Season.Rds")

gam_model <- read_rds(file = "models/MLD_LAT_LON_Season.Rds")
gam_model %>% gam.check()

gam_model_weights_high_k <- gam(
  cbind(Anomaly, Total - Anomaly) ~
    s(log(cleaned_N2)) +
    s(log(cleaned_mld)) +
    s(LATITUDE, LONGITUDE, by = Season, bs = "sos", k = 200) +
    Season,
  family = binomial,
  weights = Total,
  data = df_agg
)
write_rds(gam_model_weights_high_k,file = "gam_model_weights_high_k.Rds")


gam_model_weights_high_k <- bam(
  cbind(Anomaly, Total - Anomaly) ~
    s(log(cleaned_N2)) +
    s(log(cleaned_mld)) +
    s(LATITUDE, LONGITUDE, by = Season, bs = "sos", k = 1000) +
    Season,
  family = binomial,
  weights = Total,
  data = df_agg
)


write_rds(gam_model_weights_high_k,file = "gam_model_weights_high_k.Rds")



gam_model_weights_high_k <- NULL 

gam_model_weights_high_k_int <- bam(
  cbind(Anomaly, Total - Anomaly) ~
    s(log(cleaned_N2)) +
    s(log(cleaned_mld)) +
    ti(AdjustedDayOfYear, LATITUDE, bs = c('cc', 'tp'), k = c(20,20))+
    s(LATITUDE, LONGITUDE, by = Season, bs = "sos", k = 1000) +
    Season,
  family = binomial,
  weights = Total,
  data = df_agg
)

write_rds(gam_model_weights_high_k_int,file = "high_int_mod.Rds")



# Compare AIC values
AIC(gam_model_no_bin, gam_model_bin)

# Plot for non-binned model

plot(gam_model_no_bin, scheme = 2)

# Plot for binned model
plot(gam_model_bin, scheme = 2)



# 2 Spatial Autocorrelation for model without binning

library(spdep)

# Create spatial weights matrix
coordinates <- cbind(df_full_djf$LONGITUDE, df_full_djf$LATITUDE)
neighbors <- knearneigh(coordinates, k = 8,longlat = TRUE)
weights <- nb2listw(knn2nb(neighbors))

# Extract residuals
residuals_no_bin <- residuals(gam_model_no_bin, type = "pearson")

# Moran's I test
moran_test <- moran.test(residuals_no_bin, weights)
print(moran_test)

# Fit a GAMM with spatial correlation
gamm_model <- gamm(
  Anomaly ~ s(LATITUDE, LONGITUDE, bs = "sos", k = 180),
  family = binomial(link = "logit"),
  data = df_full_djf,
  correlation = corExp(form = ~ LONGITUDE + LATITUDE),
  method = "REML"
)
# 2 Spatial Autocorrelation for model WITH binning

library(spdep)

# Create spatial weights matrix
coordinates <- cbind(df_binned$LON_BIN, df_binned$LAT_BIN)
neighbors <- knearneigh(coordinates, k = 8,longlat = TRUE)
weights <- nb2listw(knn2nb(neighbors))

# Extract residuals
residuals_no_bin <- residuals(gam_model_no_bin, type = "pearson")

# Moran's I test
moran_test <- moran.test(residuals_no_bin, weights)
print(moran_test)

# Fit a GAMM with spatial correlation
gamm_model <- gamm(
  Anomaly ~ s(LAT_BIN, LON_BIN, bs = "sos", k = 180),
  family = binomial(link = "logit"),
  data = df_binned,
  correlation = corExp(form = ~ LON_BIN + LAT_BIN),
  method = "REML"
)
