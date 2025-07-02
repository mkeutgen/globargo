# --------------------------------------
# 1. Load Libraries & Define Utility Functions
# --------------------------------------

library(tidyverse)
library(lubridate)
library(ggplot2)
library(zoo)      # for rollmean
library(gsw)      # for oceanographic calculations


df_carbon_clean <- read_csv("/data/GLOBARGO/src/data/df_carbon_subduction_anom.csv")


# Function to compute finite centered differences (if needed)
finiteCenteredDiff <- function(vec) {
  diffVec <- rep(NA, length(vec))
  for(i in 2:(length(vec)-1)) {
    diffVec[i] <- (vec[i+1] - vec[i-1]) / 2
  }
  return(diffVec)
}

# Linear prediction function (used for both AOU and BBP700)
predict_AOU <- function(pressure, m, b) {
  m * pressure + b
}

# Function to replace values with NA outside specified pressure bounds
na_replace_fun <- function(values, pres_level, upper.bound, lower.bound) {
  ifelse(pres_level > upper.bound, NA,
         ifelse(pres_level < lower.bound, NA, values))
}

# Function to group and average pressures that are within a threshold (e.g., 50 m)
average_close_pressures <- function(pressures, threshold = 50) {
  df <- tibble(PRES_ADJUSTED = pressures) %>%
    arrange(PRES_ADJUSTED) %>%
    mutate(group = cumsum(c(1, diff(PRES_ADJUSTED) > threshold))) %>%
    group_by(group) %>%
    summarize(PRES_ADJUSTED = mean(PRES_ADJUSTED), .groups = "drop") %>%
    ungroup() %>%
    select(PRES_ADJUSTED)
  return(df$PRES_ADJUSTED)
}

# Function to compute a single critical pressure from a vector of PRES_ADJUSTED values.
# If several values are close (within threshold), the group with the most points is chosen.
compute_critical_pres <- function(pres_values, threshold = 50) {
  if(length(pres_values) == 0) return(NA)
  df <- tibble(PRES = pres_values) %>%
    arrange(PRES) %>%
    mutate(group = cumsum(c(1, diff(PRES) > threshold)))
  
  groups <- df %>%
    group_by(group) %>%
    summarize(avg_pres = mean(PRES), count = n(), .groups = "drop")
  
  best_group <- groups %>% filter(count == max(count)) %>% slice(1)
  return(best_group$avg_pres)
}

# Function to find maximum AOU (peak) below and above a given pressure
find_max_AOU_above_below_peak <- function(df, peak_pres, range = 150) {
  df_below <- df %>% filter(PRES_ADJUSTED < peak_pres, PRES_ADJUSTED >= (peak_pres - range))
  max_AOU_row_below <- df_below %>%
    filter(AOU == max(AOU)) %>%
    slice(1)
  
  df_above <- df %>% filter(PRES_ADJUSTED > peak_pres, PRES_ADJUSTED <= (peak_pres + range))
  max_AOU_row_above <- df_above %>%
    filter(AOU == max(AOU)) %>%
    slice(1)
  
  return(list(
    max_below = max_AOU_row_below %>% select(PRES_ADJUSTED, AOU, BBP700_ADJUSTED) %>% ungroup(),
    max_above = max_AOU_row_above %>% select(PRES_ADJUSTED, AOU, BBP700_ADJUSTED) %>% ungroup()
  ))
}

# Function to find minimum BBP700 (anomaly) below and above a given pressure
find_min_BBP700_above_below_peak <- function(df, peak_pres, range = 150) {
  df_below <- df %>% filter(PRES_ADJUSTED < peak_pres, PRES_ADJUSTED >= (peak_pres - range))
  min_BBP700_row_below <- df_below %>%
    filter(BBP700_ADJUSTED == min(BBP700_ADJUSTED)) %>%
    slice(1)
  
  df_above <- df %>% filter(PRES_ADJUSTED > peak_pres, PRES_ADJUSTED <= (peak_pres + range))
  min_BBP700_row_above <- df_above %>%
    filter(BBP700_ADJUSTED == min(BBP700_ADJUSTED)) %>%
    slice(1)
  
  return(list(
    min_below = min_BBP700_row_below %>% select(PRES_ADJUSTED, BBP700_ADJUSTED) %>% ungroup(),
    min_above = min_BBP700_row_above %>% select(PRES_ADJUSTED, BBP700_ADJUSTED) %>% ungroup()
  ))
}


add_region <- function(df) {
  df %>%
    mutate(
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


# Function to compute POC from integrated BBP700 anomaly
poc_from_bbp_700_southern_ocean <- function(bbp_700) {
  9.776e4 * (bbp_700 ^ 1.166)
}

poc_from_bbp_700_northatlantic <- function(bbp_700) {
  41550* (bbp_700) 
}

poc_from_bbp_700_subtropical <- function(bbp_700){
  5.8968e4 * bbp_700 + 2.75
}


# --------------------------------------
# 2. Define the Function to Process One Float–Cycle Combination
# --------------------------------------
# Note: The function now takes an extra argument 'anomalies_data' (e.g., df_carbon_clean)
# to compute the critical pressure from the PRES_ADJUSTED column for that float and cycle.
process_cycle <- function(wmo, cycle_number, anomalies_data) {
  message(paste("Processing float", wmo, "cycle", cycle_number))
  
  # Compute critical pressure from anomalies dataset for the given float and cycle
  anom_subset <- anomalies_data %>% filter(WMO == wmo, CYCLE_NUMBER == cycle_number)
  if(nrow(anom_subset) == 0) {
    message(paste("No anomaly records for float", wmo, "cycle", cycle_number))
    return(NULL)
  }
  # Compute a single critical pressure from the PRES_ADJUSTED values.
  critical_pres <- compute_critical_pres(anom_subset$PRES_ADJUSTED, threshold = 50)
  message(paste("Computed critical pressure:", critical_pres))
  
  # Download the float data using OneArgoR toolbox (adjust parameters as needed)
  data_df <- load_float_data(
    float_ids = wmo,
    variables = c(
      "DATA_TYPE", "PLATFORM_NUMBER", "BBP700", "BBP700_dPRES", "BBP700_ADJUSTED_QC",
      "LATITUDE", "LONGITUDE", "PROFILE_TEMP_QC", "PROFILE_DOXY_QC", "PROFILE_BBP700_QC",
      "PRES_QC", "PRES", "PRES_ADJUSTED", "PROFILE_PSAL_QC", "CHLA_QC", "CHLA_ADJUSTED",
      "CHLA_ADJUSTED_ERROR", "DOXY", "DOXY_QC", "DOXY_ADJUSTED", "DOXY_ADJUSTED_QC",
      "DOXY_ADJUSTED_ERROR", "PSAL", "PSAL_dPRES", "PSAL_ADJUSTED", "PSAL_ADJUSTED_QC",
      "TEMP", "TEMP_QC", "TEMP_dPRES", "TEMP_ADJUSTED", "TEMP_ADJUSTED_QC", "TEMP_ADJUSTED_ERROR"
    ),
    format = "dataframe"
  )
  
  # Filter for the desired cycle
  data_df <- data_df %>% filter(CYCLE_NUMBER == cycle_number) %>% ungroup()
  if(nrow(data_df) == 0) {
    message(paste("No float data found for float", wmo, "cycle", cycle_number))
    return(NULL)
  }
  
  # Calculate ancillary variables (only if DOXY is available)
  data_df <- data_df %>% 
    filter(!is.na(DOXY)) %>% 
    group_by(CYCLE_NUMBER) %>%
    mutate(
      ABS_SAL = gsw::gsw_SA_from_SP(SP = PSAL_ADJUSTED, p = PRES_ADJUSTED,
                                    longitude = first(LATITUDE),
                                    latitude = first(LATITUDE)),
      CONS_TEMP = gsw::gsw_CT_from_t(SA = ABS_SAL, t = TEMP_ADJUSTED, p = PRES_ADJUSTED),
      SAT_DOXY = gsw_O2sol(SA = ABS_SAL, CT = CONS_TEMP, p = PRES_ADJUSTED,
                           longitude = first(LONGITUDE), latitude = first(LATITUDE)),
      AOU = SAT_DOXY - DOXY_ADJUSTED
    ) %>%
    ungroup()
  
  # Smooth the BBP700 time series (adjust smoothing window as needed)
  data_df$BBP700_ADJUSTED <- rollmedian(data_df$BBP700, k = 3, na.pad = TRUE)
  
  # Select the relevant columns for further processing
  df <- data_df %>% 
    select(LATITUDE, LONGITUDE, TIME, CYCLE_NUMBER, PRES_ADJUSTED, AOU, BBP700_ADJUSTED)
  
  # Create a long format for plotting multiple anomaly profiles
  df_long <- df %>% pivot_longer(cols = c(AOU, BBP700_ADJUSTED))
  
  # --- Figure 1: Raw Anomaly Profiles ---
  p1 <- ggplot(df_long, aes(x = PRES_ADJUSTED, y = value)) +
    facet_grid(. ~ name, scales = "free") +
    geom_line() +
    coord_flip() +
    scale_x_reverse(limits = c(max(df$PRES_ADJUSTED, na.rm = TRUE), 0),
                    breaks = seq(0, max(df$PRES_ADJUSTED, na.rm = TRUE), by = 40)) +
    geom_vline(xintercept = critical_pres, linetype = "dashed") +
    theme_bw() +
    ggtitle(paste("Anomaly Profiles for float", wmo, "cycle", cycle_number))
  
  anomaly_profile_file <- paste0("anomaly_profiles_", wmo, "_cycle", cycle_number, ".png")
  ggsave(anomaly_profile_file, plot = p1, width = 8, height = 6)
  
  # --- Determine the Peak for Integration Using the Computed Critical Pressure ---
  peak_aou <- find_max_AOU_above_below_peak(df, peak_pres = critical_pres)
  peak_bbp700 <- find_min_BBP700_above_below_peak(df, peak_pres = critical_pres)
  
  pres_bl <- peak_aou$max_below$PRES_ADJUSTED
  aou_bl <- peak_aou$max_below$AOU
  pres_up <- peak_aou$max_above$PRES_ADJUSTED
  aou_up <- peak_aou$max_above$AOU
  
  pres_bbp_bl <- peak_bbp700$min_below$PRES_ADJUSTED
  bbp_bl <- peak_bbp700$min_below$BBP700_ADJUSTED
  pres_bbp_up <- peak_bbp700$min_above$PRES_ADJUSTED
  bbp_up <- peak_bbp700$min_above$BBP700_ADJUSTED
  
  # --- Linear Interpolation for AOU and BBP700 ---
  m_aou <- (aou_up - aou_bl) / (pres_up - pres_bl)
  b_aou <- aou_bl - m_aou * pres_bl
  
  m_bbp <- (bbp_up - bbp_bl) / (pres_bbp_up - pres_bbp_bl)
  b_bbp <- bbp_bl - m_bbp * pres_bbp_bl
  
  df <- df %>%
    mutate(
      predicted_AOU = predict_AOU(PRES_ADJUSTED, m = m_aou, b = b_aou),
      predicted_BBP700_ADJUSTED = predict_AOU(PRES_ADJUSTED, m = m_bbp, b = b_bbp)
    ) %>%
    mutate(
      predicted_AOU = na_replace_fun(predicted_AOU, PRES_ADJUSTED, pres_up, pres_bl),
      predicted_BBP700_ADJUSTED = na_replace_fun(predicted_BBP700_ADJUSTED, PRES_ADJUSTED, pres_bbp_up, pres_bbp_bl)
    )
  
  # --- Figure 2: AOU Observations vs. Prediction ---
  p2 <- ggplot(df, aes(x = PRES_ADJUSTED, y = AOU, color = "Observed AOU")) +
    geom_line(size=2) +
    coord_flip() +
    scale_x_reverse() +
    geom_line(aes(y = predicted_AOU, color = "Predicted AOU (Chen's method)"),size=1.5) +
    theme_bw() +
    theme(legend.position = "bottom") +labs(y="Apparent Oxygen Utilization (µmol/kg)",x="Adjusted Pressure (dbar)")+
    ggtitle(paste("AOU Anomaly for float", wmo, "cycle", cycle_number))+scale_color_viridis_d(begin = 0.5)
  
  aou_plot_file <- paste0("/data/GLOBARGO/figures/IntegratedAnomalies/AOU_anomaly_", wmo, "_cycle", cycle_number, ".png")
  ggsave(aou_plot_file, plot = p2, width = 8, height = 6)
  # First, add the region label based on the LATITUDE and LONGITUDE
  df <- df %>% add_region()
  
  # Compute POC based on region:
  df <- df %>%
    mutate(
      # For observed POC: use the appropriate mapping function
      POC_OBS = case_when(
        region == "Southern Ocean" ~ poc_from_bbp_700_southern_ocean(BBP700_ADJUSTED),
        region %in% c("North Atlantic","North Pacific") ~ poc_from_bbp_700_northatlantic(BBP700_ADJUSTED),
        # For all other regions, apply the subtropical mapping
        region %in% c("Northern Tropics", "Southern Tropics") ~ poc_from_bbp_700_subtropical(BBP700_ADJUSTED),
        TRUE ~ NA_real_
      ),
      # For predicted POC: use the corresponding predicted values
      POC_PRED = case_when(
        region == "Southern Ocean" ~ poc_from_bbp_700_southern_ocean(predicted_BBP700_ADJUSTED),
        region %in% c("North Atlantic","North Pacific") ~ poc_from_bbp_700_northatlantic(predicted_BBP700_ADJUSTED),
        region %in% c("Northern Tropics", "Southern Tropics") ~ poc_from_bbp_700_subtropical(predicted_BBP700_ADJUSTED),
        TRUE ~ NA_real_
      )
    )
  
  # --- Figure 3: BBP700 Observations vs. Prediction ---
  p3 <- ggplot(df, aes(x = PRES_ADJUSTED, y = POC_OBS, color = "Observed POC")) +
    geom_line(size = 2) +
    coord_flip() +
    scale_x_reverse() +
    geom_line(aes(y = POC_PRED, color = "Predicted POC (Chen's method)"), size = 1.5) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(y = "Particulate Organic Carbon (mg/m³)", x = "Adjusted Pressure (dbar)") +
    ggtitle(paste("POC Anomaly for float", wmo, "cycle", cycle_number)) +
    scale_color_viridis_d(begin = 0.5)
  
  poc_plot_file <- paste0("/data/GLOBARGO/figures/IntegratedAnomalies/POC_anomaly_", wmo, "_cycle", cycle_number, ".png")
  ggsave(poc_plot_file, plot = p3, width = 8, height = 6)
  
  # --- Integration of Anomalies ---
  distances <- diff(df$PRES_ADJUSTED)
  differences_aou <- abs(ifelse(df$predicted_AOU - df$AOU > 0,
                                df$predicted_AOU - df$AOU, NA))
  differences_bbp <- abs(ifelse(df$predicted_BBP700_ADJUSTED - df$BBP700_ADJUSTED < 0,
                                df$predicted_BBP700_ADJUSTED - df$BBP700_ADJUSTED, 0))
  
  # Correct pointwise conversion BEFORE integration – 
  # optionally, you might also adjust this conversion based on region if needed.
  # Adjust the pointwise conversion based on region:
  poc_anomaly <- with(df, case_when(
    region == "Southern Ocean" ~ poc_from_bbp_700_southern_ocean(differences_bbp),
    region %in% c("North Atlantic","North Pacific") ~ poc_from_bbp_700_northatlantic(differences_bbp),
    region %in% c("Northern Tropics", "Southern Tropics") ~ poc_from_bbp_700_subtropical(differences_bbp),
    TRUE ~ NA_real_
  ))  
  avg_differences_aou <- (differences_aou[-length(differences_aou)] + differences_aou[-1]) / 2
  trapezoid_areas_aou <- avg_differences_aou * distances
  
  avg_poc_anomaly <- (poc_anomaly[-length(poc_anomaly)] + poc_anomaly[-1]) / 2
  trapezoid_area_poc <- avg_poc_anomaly * distances
  total_integrated_poc <- sum(trapezoid_area_poc, na.rm = TRUE)
  
  # Note: Ensure that the variable 'trapezoid_areas_bbp' is defined if you plan to use it.
  total_area_aou <- sum(trapezoid_areas_aou, na.rm = TRUE)
  #total_area_bbp <- sum(trapezoid_areas_bbp, na.rm = TRUE)
  
  # Return summary results for this float–cycle
  tibble(
    WMO = wmo,
    CYCLE_NUMBER = cycle_number,
    critical_pres = critical_pres,
    total_area_aou = total_area_aou,
    # total_area_bbp = total_area_bbp,
    integrated_poc = total_integrated_poc
  )
}

# --------------------------------------
# 3. Loop Over All Float / Cycle Pairs
# --------------------------------------
# Here we assume your anomalies dataset is stored in 'df_carbon_clean'
# which includes at least the columns: WMO, CYCLE_NUMBER, and PRES_ADJUSTED.
float_cycle_list <- df_carbon_clean %>% distinct(WMO, CYCLE_NUMBER)

# Initialize an empty list to store results
results_list <- list()

# Loop over each combination, process and store the output
for(i in seq_len(nrow(float_cycle_list))) {
  wmo_val <- float_cycle_list$WMO[i]
  cycle_val <- float_cycle_list$CYCLE_NUMBER[i]
  
  res <- tryCatch(
    process_cycle(wmo = wmo_val, cycle_number = cycle_val, anomalies_data = df_carbon_clean),
    error = function(e) {
      message(paste("Error processing float", wmo_val, "cycle", cycle_val, ":", e))
      NULL
    }
  )
  
  if(!is.null(res)) {
    results_list[[length(results_list) + 1]] <- res
  }
}

# Combine the results into a summary tibble and write to CSV
summary_results <- bind_rows(results_list)
summary_results %>% View()


write_csv(summary_results, "/data/GLOBARGO/src/data/integrated_anomalies_summary.csv")
summary_results <- read_csv("/data/GLOBARGO/src/data/integrated_anomalies_summary.csv")

# Filter out extreme values (e.g., values above 5000)
robust_data <- summary_results %>% filter(integrated_poc <= 5000)

ggplot(robust_data, aes(x = integrated_poc)) +
  geom_histogram(bins = 60,position = "dodge") +
  labs(x = "Integrated POC",
       y = "Frequency",
       title = "Histogram of Integrated POC (Values <= 5000, removed 7 outliers)") +
  theme_minimal()+scale_fill_viridis()

integ_poc_fig <- ggplot(robust_data, aes(x = integrated_poc)) +
  geom_histogram(bins = 60, fill = viridis::viridis(1)) +
  labs(x = "Integrated POC",
       y = "Frequency",
       title = "Histogram of Integrated POC (Values <= 5000, removed 7 outliers)") +
  theme_minimal(base_size = 25)

ggsave("figures/integ_poc_fig.png",plot = integ_poc_fig)

robust_data$integrated_poc %>% summary()

View(summary_results)

robust_data <- robust_data %>%
  left_join(
    df_carbon_clean %>% select(WMO, CYCLE_NUMBER, LATITUDE, LONGITUDE,TIME),
    by = c("WMO", "CYCLE_NUMBER")
  )

write_csv(robust_data,file="/data/GLOBARGO/src/data/df_carbon_subduction_anom_with_poc.csv")
