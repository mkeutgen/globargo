# --------------------------------------
# 1. Load Libraries & Define Utility Functions
# --------------------------------------

library(tidyverse)
library(lubridate)
library(ggplot2)
library(zoo)      # for rollmean
library(gsw)      # for oceanographic calculations

df_carbon_clean <- read_csv("/data/GLOBARGO/src/data/df_carbon_subduction_anom.csv")

# Finite centered differences function (if needed)
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
        # Example Mediterranean criteria (adjust as needed)
        LATITUDE > 30 & LATITUDE < 45 & LONGITUDE >= -5 & LONGITUDE <= 40 ~ "Mediterranean",
        LATITUDE > 30 & LONGITUDE >= -100 & LONGITUDE <= 20 ~ "North Atlantic",  
        LATITUDE > 30 & (LONGITUDE < -100 | LONGITUDE > 120) ~ "North Pacific",  
        LATITUDE >= 0 & LATITUDE <= 30 ~ "Northern Tropics",                     
        LATITUDE < 0 & LATITUDE >= -30 ~ "Southern Tropics",                     
        LATITUDE < -30 ~ "Southern Ocean",                                       
        TRUE ~ NA_character_
      )
    )
}
# --- Previous constant conversion functions (kept for reference) ---
poc_from_bbp_700_southern_ocean <- function(bbp_700) {
  9.776e4 * (bbp_700 ^ 1.166)
}

poc_from_bbp_700_northatlantic <- function(bbp_700) {
  41550 * bbp_700 
}

poc_from_bbp_700_subtropical <- function(bbp_700){
  5.8968e4 * bbp_700 + 2.75
}

# --------------------------------------
# NEW FUNCTIONS: Depth-Dependent Conversion & MLD Estimation
# --------------------------------------

# Estimate the Mixed Layer Depth (MLD) using a simple temperature threshold method.
estimate_MLD <- function(PRES, TEMP, delta = 0.2) {
  # Assumes PRES and TEMP are ordered from surface to depth.
  surface_temp <- TEMP[1]
  index <- which(TEMP < (surface_temp - delta))
  if(length(index) == 0) {
    return(max(PRES, na.rm = TRUE))
  } else {
    return(PRES[min(index)])
  }
}

# Function to compute a depth-dependent conversion factor profile and convert bbp700 to POC.
poc_from_bbp700_depth_dependent <- function(bbp700, zvec, mld, region) {
  # Fixed parameters from Gula's method
  c_value <- 12011       # asymptotic conversion factor at depth (mg C/m2)
  b <- -6.57             # exponential slope
  
  # Updated region-dependent surface layer depth (zsurf, in m)
  zsurf_reg <- list(
    NASPG = 15,    # North Atlantic Subpolar Gyre (applies to NA & NP)
    MED   = 14,    # Mediterranean (from Loisel et al., 2001)
    STG   = 21,    # Global Subtropical Gyres
    SO    = 41     # Southern Ocean
  )
  
  # Updated region-dependent surface conversion factor (a_plus_c)
  a_plus_c <- list(
    NASPG = 41550,  # North Atlantic Subpolar Gyre (from Cetinic et al., 2012)
    MED   = 41305,  # Mediterranean (from Loisel et al., 2001)
    STG   = 58968,  # Subtropical Gyres (from Stramski et al., 2008)
    SO    = 31200   # Southern Ocean (from Johnson et al., 2017)
  )
  
  if (! region %in% names(zsurf_reg)) {
    stop("Region not recognized. Please use one of: 'NASPG', 'MED', 'STG', or 'SO'.")
  }
  zsurf <- zsurf_reg[[region]]
  
  a <- a_plus_c[[region]] - c_value
  
  if(mld > max(zvec, na.rm = TRUE)) {
    cfvec <- rep(c_value, length(zvec))
  } else {
    cfvec <- rep(c_value + a, length(zvec))
    cfvec[zvec > zsurf] <- c_value + a * exp(b * (zvec[zvec > zsurf] - zsurf) * 0.001)
    imld <- which.min(abs(zvec - mld))
    if(!is.na(imld)) {
      cfvec[zvec < mld] <- cfvec[imld]
    }
  }
  
  POC_profile <- bbp700 * cfvec
  return(list(conversion_factor = cfvec, POC_profile = POC_profile))
}

# --------------------------------------
# 2. Define the Function to Process One Float–Cycle Combination
# --------------------------------------

process_cycle <- function(wmo, cycle_number, anomalies_data) {
  message(paste("Processing float", wmo, "cycle", cycle_number))
  
  # Subset the anomalies dataset for the given float and cycle.
  anom_subset <- anomalies_data %>% filter(WMO == wmo, CYCLE_NUMBER == cycle_number)
  if(nrow(anom_subset) == 0) {
    message(paste("No anomaly records for float", wmo, "cycle", cycle_number))
    return(NULL)
  }
  
  # Compute a single critical pressure from the PRES_ADJUSTED values.
  critical_pres <- compute_critical_pres(anom_subset$PRES_ADJUSTED, threshold = 50)
  message(paste("Computed critical pressure:", critical_pres))
  
  # Load the float data using OneArgoR toolbox (adjust parameters as needed)
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
  
  # Filter for the desired cycle and ensure TEMP_ADJUSTED is included for MLD estimation.
  data_df <- data_df %>% filter(CYCLE_NUMBER == cycle_number) %>% ungroup()
  if(nrow(data_df) == 0) {
    message(paste("No float data found for float", wmo, "cycle", cycle_number))
    return(NULL)
  }
  
  # Compute ancillary variables (if DOXY is available)
  data_df <- data_df %>% 
    filter(!is.na(DOXY)) %>% 
    group_by(CYCLE_NUMBER) %>%
    mutate(
      ABS_SAL = gsw::gsw_SA_from_SP(SP = PSAL_ADJUSTED, p = PRES_ADJUSTED,
                                    longitude = first(LONGITUDE),
                                    latitude = first(LATITUDE)),
      CONS_TEMP = gsw::gsw_CT_from_t(SA = ABS_SAL, t = TEMP_ADJUSTED, p = PRES_ADJUSTED),
      SAT_DOXY = gsw_O2sol(SA = ABS_SAL, CT = CONS_TEMP, p = PRES_ADJUSTED,
                           longitude = first(LONGITUDE), latitude = first(LATITUDE)),
      AOU = SAT_DOXY - DOXY_ADJUSTED
    ) %>%
    ungroup()
  
  # Smooth the BBP700 time series (adjust smoothing window as needed)
  data_df$BBP700_ADJUSTED <- rollmedian(data_df$BBP700, k = 3, na.pad = TRUE)
  
  # Select the relevant columns (include TEMP_ADJUSTED for MLD estimation)
  df <- data_df %>% 
    select(LATITUDE, LONGITUDE, TIME, CYCLE_NUMBER, PRES_ADJUSTED, AOU, BBP700_ADJUSTED, TEMP_ADJUSTED)
  
  # Create long format for plotting multiple anomaly profiles (if desired)
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
  
  anomaly_profile_file <- paste0("/data/GLOBARGO/figures/IntegratedAnomalies/anomaly_profiles_", wmo, "_cycle", cycle_number, ".png")
  ggsave(anomaly_profile_file, plot = p1, width = 8, height = 6)
  
  # --- Determine the Peak for Integration ---
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
    geom_line(size = 2) +
    coord_flip() +
    scale_x_reverse() +
    geom_line(aes(y = predicted_AOU, color = "Predicted AOU (Chen's method)"), size = 1.5) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(y = "Apparent Oxygen Utilization (µmol/kg)", x = "Adjusted Pressure (dbar)") +
    ggtitle(paste("AOU Anomaly for float", wmo, "cycle", cycle_number)) +
    scale_color_viridis_d(begin = 0.5)
  
  aou_plot_file <- paste0("/data/GLOBARGO/figures/IntegratedAnomalies/AOU_anomaly_", wmo, "_cycle", cycle_number, ".png")
  ggsave(aou_plot_file, plot = p2, width = 8, height = 6)
  
  # --- NEW: Compute POC Using gali's (Gali et al, 2022) Depth-Dependent Conversion ---
  # First, add region information:
  df <- df %>% add_region()
  
  # Map region into method identifiers:
  df <- df %>% add_region() %>%
    mutate(
      region_method = case_when(
        region %in% c("North Atlantic", "North Pacific") ~ "NASPG",
        region == "Southern Ocean" ~ "SO",
        region %in% c("Northern Tropics", "Southern Tropics") ~ "STG",
        region == "Mediterranean" ~ "MED",
        TRUE ~ NA_character_
      )
    )
  
  # Take the first value (assuming spatial location is roughly constant)
  region_cast <- unique(df$region_method)[1]
  
  # Estimate the Mixed Layer Depth (MLD) using TEMP_ADJUSTED.
  mld <- estimate_MLD(df$PRES_ADJUSTED, df$TEMP_ADJUSTED, delta = 0.2)
  message(paste("Estimated MLD:", mld))
  
  # Compute POC profiles using the depth-dependent conversion for both observed and predicted BBP700.
  gali_obs <- poc_from_bbp700_depth_dependent(bbp700 = df$BBP700_ADJUSTED,
                                              zvec = df$PRES_ADJUSTED,
                                              mld = mld,
                                              region = region_cast)
  df$POC_OBS <- gali_obs$POC_profile
  
  gali_pred <- poc_from_bbp700_depth_dependent(bbp700 = df$predicted_BBP700_ADJUSTED,
                                               zvec = df$PRES_ADJUSTED,
                                               mld = mld,
                                               region = region_cast)
  df$POC_PRED <- gali_pred$POC_profile
  
  # --- Figure 3: POC Observations vs. Prediction ---
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
  
  # Adjust the pointwise conversion for POC anomaly using gali's method
  poc_anomaly <- abs(ifelse(df$POC_PRED-df$POC_OBS < 0,df$POC_PRED-df$POC_OBS,0))
  
  
  avg_differences_aou <- (differences_aou[-length(differences_aou)] + differences_aou[-1]) / 2
  trapezoid_areas_aou <- avg_differences_aou * distances
  
  avg_poc_anomaly <- (poc_anomaly[-length(poc_anomaly)] + poc_anomaly[-1]) / 2
  trapezoid_area_poc <- avg_poc_anomaly * distances
  total_integrated_poc <- sum(trapezoid_area_poc, na.rm = TRUE)
  
  total_area_aou <- sum(trapezoid_areas_aou, na.rm = TRUE)
  
  # Return summary results for this float–cycle
  tibble(
    WMO = wmo,
    CYCLE_NUMBER = cycle_number,
    critical_pres = critical_pres,
    total_area_aou = total_area_aou,
    integrated_poc = total_integrated_poc
  )
}

# --------------------------------------
# 3. Loop Over All Float / Cycle Pairs
# --------------------------------------
float_cycle_list <- df_carbon_clean %>% distinct(WMO, CYCLE_NUMBER)

results_list <- list()

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

summary_results <- bind_rows(results_list)
summary_results %>% View()

write_csv(summary_results, "/data/GLOBARGO/src/data/integrated_anomalies_summary_fromgali.csv")
summary_results <- read_csv("/data/GLOBARGO/src/data/integrated_anomalies_summary_fromgali.csv")

# Filter out extreme values (e.g., values above 5000)
robust_data <- summary_results %>% filter(integrated_poc <= 5000)

ggplot(robust_data, aes(x = integrated_poc)) +
  geom_histogram(bins = 60, position = "dodge") +
  labs(x = "Integrated POC",
       y = "Frequency",
       title = "Histogram of Integrated POC (Values <= 5000, removed 7 outliers)") +
  theme_minimal() + scale_fill_viridis()

integ_poc_fig <- ggplot(robust_data, aes(x = integrated_poc)) +
  geom_histogram(bins = 60, fill = viridis::viridis(1)) +
  labs(x = "Integrated POC",
       y = "Frequency",
       title = "Histogram of Integrated POC (Values <= 5000, removed 7 outliers)") +
  theme_minimal(base_size = 25)

ggsave("figures/integ_poc_fig.png", plot = integ_poc_fig)

robust_data$integrated_poc %>% summary()

robust_data <- robust_data %>%
  left_join(
    df_carbon_clean %>% select(WMO, CYCLE_NUMBER, LATITUDE, LONGITUDE, TIME),
    by = c("WMO", "CYCLE_NUMBER")
  )

write_csv(robust_data, file = "/data/GLOBARGO/src/data/df_carbon_subduction_anom_with_poc_fromgali.csv")
