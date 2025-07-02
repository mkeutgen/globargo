
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

theme_map <- function(base_size = 18){
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid      = element_blank(),
      panel.border    = element_rect(colour = "black", fill = NA, size = 0.3),
      axis.ticks      = element_line(colour = "black", size = 0.3),
      axis.text       = element_text(size = base_size * 0.8),
      plot.title      = element_text(hjust = 0.5, face = "bold",
                                     size  = base_size * 1.05,
                                     margin = margin(b = 6)),
      legend.position = "bottom",
      legend.text     = element_text(size = base_size * 0.8),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width  = unit(2.2 ,  "cm"),
      plot.margin       = margin(5, 5, 5, 5)
    )
}


# Linear prediction function (used for both AOU and BBP700)
predict <- function(pressure, m, b) {
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

# Function to find minimum BBP700 (anomaly) below and above a given pressure
find_min_POC_above_below_peak <- function(df, peak_pres, range = 150) {
  df_below <- df %>% filter(PRES_ADJUSTED < peak_pres, PRES_ADJUSTED >= (peak_pres - range))
  min_POC_OBS_row_below <- df_below %>%
    filter(POC_OBS == min(POC_OBS)) %>%
    slice(1)
  
  df_above <- df %>% filter(PRES_ADJUSTED > peak_pres, PRES_ADJUSTED <= (peak_pres + range))
  min_POC_OBS_row_above <- df_above %>%
    filter(POC_OBS == min(POC_OBS)) %>%
    slice(1)
  
  return(list(
    min_below = min_POC_OBS_row_below %>% select(PRES_ADJUSTED, POC_OBS) %>% ungroup(),
    min_above = min_POC_OBS_row_above %>% select(PRES_ADJUSTED, POC_OBS) %>% ungroup()
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


anomalies_data = df_carbon_clean 

wmo <- 5906317
cycle_number <- 25

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
      select(LATITUDE, LONGITUDE, TIME, CYCLE_NUMBER, PRES_ADJUSTED,
             AOU, BBP700_ADJUSTED, TEMP_ADJUSTED,ABS_SAL)
    
    
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
    
    
    
    # Create long format for plotting multiple anomaly profiles (if desired)
    df_long <- df %>% pivot_longer(cols = c(AOU, BBP700_ADJUSTED,ABS_SAL,POC_OBS))
    
    peak_POC_OBS <- find_min_POC_above_below_peak(df, peak_pres = critical_pres)
    
    pres_poc_bl <- peak_POC_OBS$min_below$PRES_ADJUSTED
    poc_bl <- peak_POC_OBS$min_below$POC_OBS
    pres_poc_up <- peak_POC_OBS$min_above$PRES_ADJUSTED
    poc_up <- peak_POC_OBS$min_above$POC_OBS
    
    # --- Linear Interpolation for BBP700 ---
    
    m_poc <- (poc_up - poc_bl) / (pres_poc_up - pres_poc_bl)
    b_poc <- poc_bl - m_poc * pres_poc_bl
    
    df <- df %>%
      mutate(
        predicted_POC = predict(PRES_ADJUSTED, m = m_poc, b = b_poc)
      ) %>%
      mutate(
        predicted_POC = na_replace_fun(predicted_POC, PRES_ADJUSTED, pres_poc_up, pres_poc_bl)
      )
    
    
    saveRDS(object = df,file = "/data/GLOBARGO/src/data/nice_profile_data.Rds") 
    
    # ────────────────────────────────────────────────────────────────
    # 0.  Common palette & theme  (same sizes as your maps)
    # ────────────────────────────────────────────────────────────────
    pal_prof <- viridis::viridis(2, option = "D", end = 0.8)
    names(pal_prof) <- c("Observed", "Interpolated")
    
    theme_profile <- theme_map(base_size = 18) %+replace%
      theme(
        legend.position   = "bottom",
        legend.direction  = "horizontal",
        legend.key.width  = unit(2.2, "cm"),
        legend.key.height = unit(0.5, "cm"),
        axis.text.x       = element_text(angle = 45, hjust = 1,margin = margin(b = 15)),
      )
    
    # helper for the x‑axis (pressure) so all four plots are identical
    scale_pressure <- list(
      coord_flip(),
      scale_x_reverse(limits = c(1000, 0), expand = c(0, 0)),
      labs(x = "Depth (m)")
    )
    
    # ────────────────────────────────────────────────────────────────
    # 1.  Build each panel
    # ────────────────────────────────────────────────────────────────
    POC_plot <- ggplot(df, aes(PRES_ADJUSTED)) +  
      geom_rect(
        aes(xmin = 0, xmax = 200, ymin = -Inf, ymax = Inf),
        fill = "grey90", alpha = 0.1
      )+
      geom_vline(color="#FDE725FF",xintercept = critical_pres,alpha=.2,size=20)+
      geom_line(aes(y = POC_OBS,        colour = "Observed"),     size = 1.4) +
      geom_point(aes(y = POC_OBS, colour = "Observed"), size = 3,alpha=.5) +
      geom_line(aes(y = predicted_POC,  colour = "Interpolated"), size = 1.0) +
      scale_colour_manual(values = pal_prof, name = NULL) +
      scale_pressure +
      labs(y = "Particulate organic carbon (mg m⁻³)",
           title = paste0("c • POC")) +
      theme_profile
    
    #bbp_plot <- ggplot(df, aes(PRES_ADJUSTED)) +
    #  geom_vline(color="#FDE725FF",xintercept = 940,alpha=.1,size=20)+
    #  geom_line(aes(y = BBP700_ADJUSTED, colour = "Observed"), size = 1.4) +
    #  scale_colour_manual(values = pal_prof["Observed"], drop = FALSE, name = NULL) +
    #  scale_pressure +
    #  labs(y = expression("bbp(700 nm) (m"^{-1}*" sr"^{-1}*")"),
    #       title = paste0("c • Backscatter — float ", wmo, " cycle ", cycle_number)) +
    #  theme_profile
    
    
    ABS_SAL_plot <- ggplot(df, aes(PRES_ADJUSTED)) +
      geom_rect(
        aes(xmin = 0, xmax = 200, ymin = -Inf, ymax = Inf),
        fill = "grey90", alpha = 0.1
      )+
      geom_vline(color="#FDE725FF",xintercept = critical_pres,alpha=.2,size=20)+
      geom_line(aes(y = ABS_SAL, colour = "Observed"), size = 1.4) +
      geom_point(aes(y = ABS_SAL, colour = "Observed"), size = 3,alpha=.5) +
      scale_colour_manual(values = pal_prof["Observed"], drop = FALSE, name = NULL) +
      scale_pressure +
      labs(y = "Absolute salinity (g kg⁻¹)",
           title = paste0("a • Salinity"))+ theme_profile
    
    AOU_plot <- ggplot(df, aes(PRES_ADJUSTED)) +
      geom_rect(
        aes(xmin = 0, xmax = 200, ymin = -Inf, ymax = Inf),
        fill = "grey90", alpha = 0.1
      )+geom_vline(color="#FDE725FF",xintercept = critical_pres,alpha=.2,size=20)+
      geom_line(aes(y = AOU, colour = "Observed"), size = 1.4) +
      geom_point(aes(y = AOU, colour = "Observed"), size = 3,alpha=.5) +
      scale_colour_manual(values = pal_prof["Observed"], drop = FALSE, name = NULL) +
      scale_pressure +
      labs(y = "Apparent oxygen utilisation (µmol kg⁻¹)",
           title = paste0("b • AOU")) +
      theme_profile
    
    # ────────────────────────────────────────────────────────────────
    # 2.  Combine & export
    # ────────────────────────────────────────────────────────────────
    profile_grid <- (ABS_SAL_plot | AOU_plot |  POC_plot) +
      plot_layout(guides = "collect") &   # single legend row
      theme(legend.position = "none")

    combined_plot_file <- paste0("/data/GLOBARGO/figures/CombinedPOCFigures/POC_anomaly_", wmo, "_cycle", cycle_number, ".png")
    ggsave(combined_plot_file, plot = profile_grid, width = 16, height = 10)
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

