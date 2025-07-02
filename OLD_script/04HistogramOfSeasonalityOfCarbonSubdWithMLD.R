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
library(mgcv)
library(patchwork)

library(fitdistrplus)
library(poweRlaw)
# Resolve function conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")

df_argo_clean <- read_csv("data/df_argo_loc.csv")

df_argo_clean <- df_argo_clean %>% distinct(WMO,TIME,.keep_all = TRUE)

df_complete_clean <- read_csv("data/df_eddy_subduction_anom.csv")
df_carbon_clean <- read_csv("data/df_carbon_subduction_anom.csv")

df_argo_full <- assign_bins(df_argo_clean)
df_argo_djf <- assign_bins(df_argo_djf)
df_argo_mam <- assign_bins(df_argo_mam)
df_argo_jja <- assign_bins(df_argo_jja)
df_argo_son <- assign_bins(df_argo_son)

df_complete_full <- assign_bins(df_complete_clean)
df_complete_djf <- assign_bins(df_complete_djf)
df_complete_mam <- assign_bins(df_complete_mam)
df_complete_jja <- assign_bins(df_complete_jja)
df_complete_son <- assign_bins(df_complete_son)

df_carbon_full <- assign_bins(df_carbon_clean)
df_carbon_djf <- assign_bins(df_carbon_djf)
df_carbon_mam <- assign_bins(df_carbon_mam)
df_carbon_jja <- assign_bins(df_carbon_jja)
df_carbon_son <- assign_bins(df_carbon_son)



# Compute counts of total profiles and anomalies for each season
compute_counts <- function(df_argo, df_complete) {
  total_counts <- df_argo %>%
    group_by(lon_bin, lat_bin) %>%
    summarize(count_total = n(), .groups = 'drop')
  
  anomaly_counts <- df_complete %>%
    group_by(lon_bin, lat_bin) %>%
    summarize(count_anomaly = n(), .groups = 'drop')
  
  # Merge counts and compute proportions
  merged_counts <- full_join(total_counts, anomaly_counts, by = c("lon_bin", "lat_bin")) %>%
    mutate(
      count_anomaly = ifelse(is.na(count_anomaly), 0, count_anomaly),
      proportion = count_anomaly / count_total
    ) %>%
    filter(count_total > 0) %>%
    filter(proportion < 0.8) %>%  # Filter out high proportions (anomalies > 80%)
    filter(!is.na(lat_bin))
  
  return(merged_counts)
}
# Compute merged_counts of SUBDUCTION for whole year and each season

merged_counts_full <- compute_counts(df_argo_full,df_complete_full)
merged_counts_djf <- compute_counts(df_argo_djf, df_complete_djf)
merged_counts_mam <- compute_counts(df_argo_mam, df_complete_mam)
merged_counts_jja <- compute_counts(df_argo_jja, df_complete_jja)
merged_counts_son <- compute_counts(df_argo_son, df_complete_son)

# Compute merged_counts of CARBON SUBDUCTION for each season
merged_carbon_counts_full <- compute_counts(df_argo_full,df_carbon_full)
merged_carbon_counts_djf <- compute_counts(df_argo_djf, df_carbon_djf)
merged_carbon_counts_mam <- compute_counts(df_argo_mam, df_carbon_mam)
merged_carbon_counts_jja <- compute_counts(df_argo_jja, df_carbon_jja)
merged_carbon_counts_son <- compute_counts(df_argo_son, df_carbon_son)



# Now let's do histograms faceted by region showing the probability of (carbon) subduction
# by months :

###############################################################################
# Seasonal Cycle of Subduction (Carbon and Non-Carbon) by Region and Month
###############################################################################

# Add Month and Region Information to Data
# Splitting Northern Tropics and Southern Tropics as a sanity check (should not be seasonal there)
add_month_region <- function(df) {
  df %>%
    mutate(
      month = month(TIME, label = TRUE, abbr = TRUE),  # Extract month as a factor
      region = case_when(
        LATITUDE > 30 & LONGITUDE >= -100 & LONGITUDE <= 20 ~ "North Atlantic",  # Above 30°N
        LATITUDE > 30 & (LONGITUDE < -100 | LONGITUDE > 120) ~ "North Pacific",  # Above 30°N
        LATITUDE >= 0 & LATITUDE <= 30 ~ "Northern Tropics",                     # 0° to 30°N
        LATITUDE < 0 & LATITUDE >= -30 ~ "Southern Tropics",                     # 0° to 30°S
        LATITUDE < -30 ~ "Southern Ocean",                                       # Below 30°S
        TRUE ~ NA_character_  # Exclude undefined regions
      )
    )
}

# Not Splitting Northern Tropics and Southern Tropics as a sanity check (should not be seasonal there)
# add_month_region <- function(df) {
#   df %>%
#     mutate(
#       month = month(TIME, label = TRUE, abbr = TRUE),  # Extract month as a factor
#       region = case_when(
#         LATITUDE > 30 & LONGITUDE >= -100 & LONGITUDE <= 20 ~ "North Atlantic",  # Above 30°N
#         LATITUDE > 30 & (LONGITUDE < -100 | LONGITUDE > 120) ~ "North Pacific",  # Above 30°N
#         LATITUDE <= 30 & LATITUDE >= -30 ~ "Tropics",                     # 0° to 30°S
#         LATITUDE < -30 ~ "Southern Ocean",                                       # Below 30°S
#         TRUE ~ NA_character_  # Exclude undefined regions
#       )
#     )
# }
#
#
df_argo_clean <- add_month_region(df_argo_clean)
df_complete_clean <- add_month_region(df_complete_clean)
df_carbon_clean <- add_month_region(df_carbon_clean)
#
# # Step 2: Compute Probabilities by Region and Month
compute_monthly_probabilities <- function(df_argo, df_events) {
  # Total Argo profiles by region and month
  total_counts <- df_argo %>%
    group_by(region, month) %>%
    summarize(count_total = n(), .groups = "drop")
  
  # Subduction events by region and month
  event_counts <- df_events %>%
    group_by(region, month) %>%
    summarize(count_event = n(), .groups = "drop")
  
  # Merge and compute proportions
  merged <- full_join(total_counts, event_counts, by = c("region", "month")) %>%
    mutate(
      count_event = replace_na(count_event, 0),
      proportion = ifelse(count_total > 0, count_event / count_total, NA)
    )
  return(merged)
}

# Compute probabilities for subduction and carbon subduction
monthly_probs_subduction <- compute_monthly_probabilities(df_argo_clean, df_complete_clean)
monthly_probs_carbon <- compute_monthly_probabilities(df_argo_clean, df_carbon_clean)

monthly_probs_subduction$carbon <- "Subduction with/without carbon"
monthly_probs_carbon$carbon <- "Subduction with carbon"

monthly_probs_binded <- bind_rows(monthly_probs_carbon,monthly_probs_subduction) %>%
  filter(!is.na(region), !is.na(proportion), !is.na(month))  # Remove NA regions, proportions, and months %>%


# Step 3: Plot Monthly Histograms Faceted by Region
plot_monthly_histograms <- function(monthly_probs, title) {
  monthly_probs <- monthly_probs %>%
    filter(!is.na(region), !is.na(proportion), !is.na(month))  # Remove NA regions, proportions, and months
  
  ggplot(monthly_probs, aes(x = month, y = proportion, fill = region)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.8) +
    facet_wrap(~region, ncol = 2,axes="all") +
    scale_fill_viridis_d() +
    labs(
      title = title,
      x = "Month",
      y = "Probability of Subduction",
      fill = "Region"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 15,angle = 45, hjust = 1),
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 16),
      title = element_text(size=15)
    )+  
    guides(fill = "none")
}


# Plot subduction probabilities
plot_subduction <- plot_monthly_histograms(monthly_probs_subduction, "Probability that an Argo Float captures a subduction event")

# Plot carbon subduction probabilities
plot_carbon_subduction <- plot_monthly_histograms(monthly_probs_carbon, "Probability that an Argo Float captures a subduction event with carbon")



####################################################
##### SEASONALITY OF CARBON SUBDUCTION WITH CHLORO #
####################################################
df_chloro_summary <- read_csv(file = "/data/GLOBARGO/data/df_chloro_monthly_climatology_seawifs.csv")
df_mld_summary <- read_csv(file = "data/mld_results.csv")
df_N2_summary <- read_csv(file = "data/N2_results.csv")

df_mld_summary$TIME <- df_mld_summary$Time
df_N2_summary$TIME <- df_N2_summary$Time
df_mld_summary
df_argo_clean

df_mld_merged <- left_join(
  df_mld_summary,
  df_argo_clean,
  by = c("WMO", "Time" = "TIME")
) %>% select(WMO,TIME,LATITUDE,LONGITUDE,MLD)



df_N2_merged <- left_join(
  df_N2_summary,
  df_argo_clean,
  by = c("WMO", "Time" = "TIME")
) %>% select(WMO,TIME,LATITUDE,LONGITUDE,cleaned_N2)


N2_summary <- add_month_region(df_N2_merged) %>%
  group_by(region,month) %>% summarise(N2=mean(cleaned_N2))

N2_df <- add_month_region(df_N2_merged) %>%
  group_by(region,month) 

mld_summary <- add_month_region(df_mld_merged) %>%
  group_by(region,month) %>% summarise(MLD=mean(MLD))

N2_summary$MONTH <- N2_summary$month %>% as.numeric()
mld_summary$MONTH <- mld_summary$month %>% as.numeric()


df_chloro_summary <- df_chloro_summary %>% na.omit()

##
# ----- STEP 1: Merge the datasets -----
monthly_probs_binded %>% str()
# Convert monthly_probs_carbon$month from abbreviated factor to numeric month
# (Assuming monthly_probs_carbon$month is a factor like "Jan", "Feb", etc.)
monthly_probs_binded <- monthly_probs_binded %>%
  mutate(MONTH = as.numeric(factor(month,
                                   levels = c("Jan","Feb","Mar","Apr","May","Jun",
                                              "Jul","Aug","Sep","Oct","Nov","Dec"))))



df_chloro_summary$carbon <- "Subduction with carbon"
mld_summary$carbon <- "Subduction with carbon"
N2_summary$carbon <- "Subduction with carbon"

# Merge by region and month (using an inner join so that only matching rows are kept)
merged_data <- monthly_probs_binded %>%
  left_join(df_chloro_summary, by = c("region", "MONTH","carbon")) %>%
  select(region,MONTH,CHLA_mean,proportion,carbon)

merged_data_with_mld <- merged_data %>% 
  left_join(mld_summary, by = c("region", "MONTH","carbon"))

merged_data_full <- merged_data_with_mld %>% 
  left_join(N2_summary, by = c("region", "MONTH","carbon")) %>%
  select(region,proportion,MONTH,CHLA_mean,MLD,N2,carbon)

write_csv(merged_data_full,"/data/GLOBARGO/src/data/histogram_data_full_for_python.csv") 
merged_data_full$carbon %>% unique()
# Assume your dataset is named 'merged_data_full'
# and contains the columns: region, MONTH, proportion, CHLA_mean, and MLD.
df <- merged_data_full %>% 
  filter(region %in% c("North Atlantic", "Southern Ocean")) %>%
  select(region, MONTH, proportion, CHLA_mean, MLD)


### Plot 1: MLD & Proportion
# Rescale MLD so that its typical values (~100) become similar to the proportion (~0.01)
df_mld <- df %>% 
  mutate(MLD_scaled = MLD * 0.0001)  # e.g. 100 becomes 0.01

### Plot 1: MLD & Proportion
df_mld <- df %>% 
  mutate(MLD_scaled = MLD * 0.0001)  # Rescale MLD (~100 becomes ~0.01)

p_mld <- ggplot(df_mld, aes(x = MONTH)) +
  # Bar plot for proportion using fill mapped to region (viridis colors)
  geom_bar(aes(y = proportion, fill = region), 
           stat = "identity", 
           position = position_dodge(width = 0.9), 
           alpha = 0.6,
           color = "black") +
  # Line plot for MLD with a fixed line color (different from bar colors)
  geom_line(aes(y = MLD_scaled, group = region), 
            color = "darkblue", size = 2) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_fill_viridis_d() +
  scale_y_continuous(
    name = "Proportion",
    limits = c(0, 0.05),
    sec.axis = sec_axis(~ . / 0.0001, name = "Mixed-layer Depth (m)")
  ) +
  facet_wrap(~ region, ncol = 1) +
  labs(title = "", x = "Month") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    # Color secondary y-axis to match the line color
    axis.text.y.right = element_text(color = "darkblue"),
    axis.title.y.right = element_text(color = "darkblue")
  )

### Plot 2: CHLA & Proportion
df_chla <- df %>% 
  mutate(CHLA_scaled = CHLA_mean * 0.1)  # Rescale CHLA (~0.1 becomes ~0.01)

p_chla <- ggplot(df_chla, aes(x = MONTH)) +
  # Bar plot for proportion using fill mapped to region
  geom_bar(aes(y = proportion, fill = region), 
           stat = "identity", 
           position = position_dodge(width = 0.9), 
           alpha = 0.6,
           color = "black") +
  # Line plot for CHLA with a fixed line color distinct from the bar colors
  geom_line(aes(y = CHLA_scaled, group = region), 
            color = "#009E73", size = 2) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_fill_viridis_d() +
  scale_y_continuous(
    name = "Proportion",
    limits = c(0, 0.05),
    sec.axis = sec_axis(~ . / 0.1, name = "Chlorophyll-a (mg m⁻³)")
  ) +
  facet_wrap(~ region, ncol = 1) +
  labs(title = "", x = "Month") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    # Color secondary y-axis to match the line color
    axis.text.y.right = element_text(color = "#009E73"),
    axis.title.y.right = element_text(color = "#009E73")
  )

### Combine the Two Plots Side by Side with a Common Title


combined_plot <- ggarrange(p_mld, p_chla,legend = "none")
annotated_plot <- annotate_figure(
  combined_plot,
  top = text_grob("Carbon subduction peaks at intermediary periods, between the peak in MLD and the peak in chlorophyll", 
                   size = 20)
)
ggsave("figures/carbon_subd_mld_chloro.png",annotated_plot,width = 15,height = 10)



mld_plot <- merged_data_full %>%
  group_by(region) %>% 
  filter(region %in% c("North Atlantic", "Southern Ocean")) %>%
  mutate(
    prop_prop = proportion / sum(proportion),
    prop_MLD = MLD / sum(MLD),
    prop_CHLA = CHLA_mean / sum(CHLA_mean)
  ) %>%
  select(region, MONTH, prop_prop, prop_MLD, prop_CHLA) %>%
  pivot_longer(cols = c("prop_prop", "prop_MLD", "prop_CHLA")) %>%
  ggplot(aes(x = MONTH, y = value, color = name)) +
  facet_wrap(. ~ region)  +
  geom_bar(
    data = ~ filter(.x, name == "prop_prop"),
    stat = "identity",
    position = "dodge",
    alpha = 0.6,
    color = NA
  ) +
  geom_line(
    data = ~ filter(.x, name != "prop_prop"),
    size = 2
  ) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_fill_manual(values = c("prop_prop" = "grey70", "prop_MLD" = NA, "prop_CHLA" = NA), guide = "none") +
  scale_color_discrete(
    breaks = c("prop_prop", "prop_MLD", "prop_CHLA"),
    labels = c("Proportion of Events", "Mixed-layer Depth", "Chlorophyll-a")
  ) + scale_color_viridis_d()+
  theme_minimal(base_size = 25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "bottom") +
  labs(x = "Month", y = "Density", color = "")

n2_plot <- merged_data_full %>%
  group_by(region) %>% 
  filter(region %in% c("North Atlantic", "Southern Ocean")) %>%
  mutate(
    prop_prop = proportion / sum(proportion),
    prop_N2 = N2 / sum(N2),
    prop_CHLA = CHLA_mean / sum(CHLA_mean)
  ) %>%
  select(region, MONTH, prop_prop, prop_N2, prop_CHLA) %>%
  pivot_longer(cols = c("prop_prop", "prop_N2", "prop_CHLA")) %>%
  ggplot(aes(x = MONTH, y = value, color = name)) +
  facet_wrap(. ~ region)  +
  geom_bar(
    data = ~ filter(.x, name == "prop_prop"),
    stat = "identity",
    position = "dodge",
    alpha = 0.6,
    color = NA
  ) +
  geom_line(
    data = ~ filter(.x, name != "prop_prop"),
    size = 2
  ) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_fill_manual(values = c("prop_prop" = "grey70", "prop_N2" = NA, "prop_CHLA" = NA), guide = "none") +
  scale_color_discrete(
    breaks = c("prop_prop", "prop_N2", "prop_CHLA"),
    labels = c("Proportion of Events", "Brunt-Väisäla Frequency (N2)", "Chlorophyll-a")
  ) + scale_color_viridis_d()+
  theme_minimal(base_size = 25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "bottom") +
  labs(x = "Month", y = "Density", color = "")


ggsave(filename = "figures/n2plot_chla_prob_subd.png",n2_plot,width = 10,height = 10)

ggsave(filename = "figures/mldplot_chla_prob_subd.png",mld_plot,width = 10,height = 10)


df <- merged_data_full %>% group_by(region) %>%  mutate(MLD_norm = MLD/sum(MLD),
                                                  N2_norm = N2/sum(N2),
                                                  CHLA_norm = CHLA_mean/sum(CHLA_mean),
                                                  prop_norm = proportion/sum(proportion),
                                                  MLDxCHLA = MLD_norm * CHLA_norm,
                                                  MLDxCHLA = MLDxCHLA/sum(MLDxCHLA))


ggplot(df,aes(x=MONTH,y=prop_norm))+
  facet_wrap(.~ region )+geom_line(aes(y=prop_norm))+geom_line(aes(y=MLDxCHLA),color='pink')
merged_data_full$Prod <- merged_data_full$MLD * merged_data_full$CHLA_mean

merged_data_full

# Make sure region is a factor (so that interactions are interpreted correctly)
merged_data <- merged_data_full %>%
  mutate(region = as.factor(region))

# Define month abbreviations for the x-axis labels.
month_labels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# Assume merged_data_full has columns: region, MONTH (numeric 1-12), CHLA_mean, MLD, N2, and proportion.
# Normalize the observed values per region and create an inverted N2 variable.
df <- merged_data_full %>% 
  group_by(region) %>%
  mutate(
    prop_norm  = proportion / sum(proportion, na.rm = TRUE),
    CHLA_norm  = CHLA_mean / sum(CHLA_mean, na.rm = TRUE),
    MLD_norm   = MLD / sum(MLD, na.rm = TRUE),
    N2_norm    = N2  / sum(N2, na.rm = TRUE),
    # Create N2_inv so that low N2 (i.e. low stratification) becomes high.
    N2_inv = 1/N2,
    N2_norm_inverse = N2_inv / sum(N2_inv, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(Month_label = factor(MONTH, levels = 1:12, labels = month_labels))

# Helper function: for a given data frame compute the product of kernel densities.
# This function estimates densities for var1 and var2 and multiplies the density values (interpolated at each observation).
predict_subduction <- function(data, var1, var2, new_var_name) {
  dens1 <- density(data[[var1]], n = 512, na.rm = TRUE)
  dens2 <- density(data[[var2]], n = 512, na.rm = TRUE)
  
  data <- data %>%
    mutate(
      dens1_val = approx(dens1$x, dens1$y, xout = .data[[var1]])$y,
      dens2_val = approx(dens2$x, dens2$y, xout = .data[[var2]])$y,
      !!new_var_name := dens1_val * dens2_val
    )
  return(data)
}

# For each region, compute the predicted signal based on:
# (a) CHLA_mean and MLD,
# (b) CHLA_mean and N2_inv (for N2, low values become high).
df_predicted <- df %>%
  group_by(region) %>%
  do({
    region_data <- .
    region_data <- predict_subduction(region_data, "CHLA_mean", "MLD", "predicted_MLD")
    region_data <- predict_subduction(region_data, "CHLA_mean", "N2",  "predicted_N2")
    region_data
  }) %>%
  ungroup()

# Normalize the predicted signals within each region (so they sum to one) for comparison.
df_predicted <- df_predicted %>%
  group_by(region) %>%
  mutate(
    predicted_MLD_norm = predicted_MLD / sum(predicted_MLD, na.rm = TRUE),
    predicted_N2_norm  = predicted_N2  / sum(predicted_N2, na.rm = TRUE)
  ) %>%
  ungroup()

# Plot: Observed vs. Predicted using CHLA x MLD.
p_mld <- ggplot(df_predicted, aes(x = MONTH)) +
  facet_wrap(~ region, ncol = 2,scales="free") +
  geom_area(aes(y = prop_norm), fill = "grey70", alpha = 0.5) +
  geom_line(aes(y = prop_norm, group = 1), color = "black", size = 1.2, alpha = 0.8) +
  geom_point(aes(y = prop_norm), color = "black", size = 3, alpha = 0.8) +
  geom_area(aes(y = predicted_MLD_norm), fill = "red", alpha = 0.4) +
  geom_line(aes(y = predicted_MLD_norm, group = 1), color = "red", size = 1.2, linetype = "dashed", alpha = 0.8) +
  geom_point(aes(y = predicted_MLD_norm), color = "red", size = 3, alpha = 0.8) +
  labs(title = "Observed vs Predicted Subduction (CHLA x MLD)",
       x = "Month", y = "Normalized Density") +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.title.y.right = element_text(color = "black", size = 16),
    axis.text.y.right  = element_text(color = "black", size = 12),
    title = element_text(size = 15)
  ) +
  guides(fill = "none")

# Plot: Observed vs. Predicted using CHLA x N2 (with N2 inverted).
p_n2 <- ggplot(df_predicted, aes(x = MONTH)) +
  facet_wrap(~ region, ncol = 2,scales="free") +
  geom_area(aes(y = prop_norm), fill = "grey70", alpha = 0.5) +
  geom_line(aes(y = prop_norm, group = 1), color = "black", size = 1.2, alpha = 0.8) +
  geom_point(aes(y = prop_norm), color = "black", size = 3, alpha = 0.8) +
  geom_area(aes(y = predicted_N2_norm), fill = "blue", alpha = 0.4) +
  geom_line(aes(y = predicted_N2_norm, group = 1), color = "blue", size = 1.2, linetype = "dashed", alpha = 0.8) +
  geom_point(aes(y = predicted_N2_norm), color = "blue", size = 3, alpha = 0.8) +
  labs(title = "Observed vs Predicted Subduction (CHLA x N2)",
       x = "Month", y = "Normalized Density") +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.title.y.right = element_text(color = "black", size = 16),
    axis.text.y.right  = element_text(color = "black", size = 12),
    title = element_text(size = 15)
  ) +
  guides(fill = "none")

# Display the plots.
print(p_mld)
print(p_n2)

# Compute a distance metric (e.g. RMSE) to compare observed vs. predicted densities for each region.
rmse_df <- df_predicted %>%
  group_by(region) %>%
  summarise(
    RMSE_MLD = sqrt(mean((prop_norm - predicted_MLD_norm)^2, na.rm = TRUE)),
    RMSE_N2  = sqrt(mean((prop_norm - predicted_N2_norm)^2, na.rm = TRUE))
  )

print(rmse_df)


# ----- STEP 2: Create a time-adjusted (lag) chlorophyll predictor -----
# We assume that if chlorophyll peaks 1-2 months after the subduction probability,
# then the chlorophyll measured one month later (lag by 1) might be a better predictor.
# Create a new variable (CHLA_lag1) by grouping by region and shifting CHLA_med by 1 month.

# ----- STEP 3: Fit GLM models -----
# Define the candidate predictors
# Assuming merged_data already contains:
#   - count_event, count_total, region, MONTH, CHLA_med, CHLA_mean
# Compute additional lag predictors (lags 1 to 11) by grouping by region

# Assuming merged_data already contains:
#   - count_event, count_total, region, MONTH, CHLA_med, CHLA_mean
# Compute additional lag predictors (lags 1 to 11) by grouping by region

cyclic_lead <- function(x, n) {
  len <- length(x)
  sapply(seq_along(x), function(i) x[((i - 1 + n) %% len) + 1])
}

merged_data <- merged_data %>%
  group_by(region) %>%
  arrange(MONTH) %>%
  mutate(
    CHLA_lead1_mean  = cyclic_lead(CHLA_mean,  1),
    CHLA_lead2_mean  = cyclic_lead(CHLA_mean,  2),
    CHLA_lead3_mean  = cyclic_lead(CHLA_mean,  3),
    CHLA_lead4_mean  = cyclic_lead(CHLA_mean,  4),
    CHLA_lead5_mean  = cyclic_lead(CHLA_mean,  5),
    CHLA_lead6_mean  = cyclic_lead(CHLA_mean,  6),
    CHLA_lead7_mean  = cyclic_lead(CHLA_mean,  7),
    CHLA_lead8_mean  = cyclic_lead(CHLA_mean,  8),
    CHLA_lead9_mean  = cyclic_lead(CHLA_mean,  9),
    CHLA_lead10_mean = cyclic_lead(CHLA_mean, 10),
    CHLA_lead11_mean = cyclic_lead(CHLA_mean, 11)
  ) %>%
  ungroup()









merged_data
plot_monthly_histograms <- function(monthly_probs, chloro_summary, title) {
  # Convert month to numeric if needed
  if (is.factor(monthly_probs$month) || is.ordered(monthly_probs$month)) {
    monthly_probs <- monthly_probs %>%
      mutate(MONTH = as.numeric(factor(month,
                                       levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))))
  } else {
    monthly_probs <- monthly_probs %>% rename(MONTH = month)
  }
  
  # Compute scaling factor for MLD: Scale so that the maximum MLD (from merged_data)
  # corresponds to the maximum probability value in monthly_probs.
  max_mld <- max(merged_data$MLD, na.rm = TRUE)
  scale_factor <- max(monthly_probs$proportion, na.rm = TRUE) / max_mld
  
  # Filter out missing values
  monthly_probs <- monthly_probs %>%
    filter(!is.na(region), !is.na(proportion), !is.na(MONTH))
  
  ggplot(monthly_probs, aes(x = MONTH)) +
    # Bar plot for probability of carbon subduction
    geom_bar(aes(y = proportion, fill = region),
             stat = "identity",
             position = "dodge",
             color = "darkblue",
             alpha = 0.8) +
    # Overlay MLD as a green line and points (scaled to match the probability scale)
    geom_line(data = merged_data,
              aes(x = MONTH, y = MLD * scale_factor, group = region),
              color = "black",
              size = 1.2) +
    geom_point(data = merged_data,
               aes(x = MONTH, y = MLD * scale_factor, group = region),
               color = "black",
               size = 2) +
    facet_wrap(~ region, ncol = 2) +
    scale_fill_viridis_d() +
    labs(
      title = title,
      x = "Month",
      y = "Probability of Subduction"
    ) +
    # Secondary axis for MLD values (restoring original MLD units)
    scale_y_continuous(
      sec.axis = sec_axis(~ . / scale_factor, name = "MLD")
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
      axis.title.y.right = element_text(color = "black", size = 16),
      axis.text.y.right  = element_text(color = "black", size = 12),
      title = element_text(size = 15)
    ) +
    guides(fill = "none")
}

# Example usage:
plot <- plot_monthly_histograms(monthly_probs_carbon, df_chloro_summary,
                                "Monthly Probability of Carbon Subduction \nwith MLD Seasonal Variability")
ggsave("figures/seasonal_carbon_subduction_by_region_with_mld.png", plot = plot, width = 10, height = 8)
