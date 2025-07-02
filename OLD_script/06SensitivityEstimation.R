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
# Resolve function conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# STILL TO DO : make sure first derivative works
# Implement a test based on 2nd derivative if it's larger than 0.01 for AOU. 


# Read the list of WMO IDs
wmolist <- readRDS("~/Documents/ARGO/BGC_Argo_WMO_PSAL_BBP_DOXY_TEMP.rds")

# Initialize an empty list to store detected eddy events
detected_events_list <- list()

# Hyperparameters :
cutoff <- 1.96
resolution <- 40
window <- 60 
# different wmos tried
wmo <- 5904677 
wmo <- 2902177
# Main processing loop over each WMO ID
for (j in seq_along(wmolist)) {
  try({
    #
    #5906635
    
    
    # Get current WMO ID
    wmo <- wmolist[j]
    
    # Load data for the current float (assuming 'load_float_data' is an internal function)
    float_data <- load_float_data(
      float_ids = wmo,
      variables = c(
        "DATA_TYPE", "PLATFORM_NUMBER", "BBP700", "BBP700_dPRES",
        "BBP700_ADJUSTED_QC", "LATITUDE", "LONGITUDE", "PROFILE_TEMP_QC",
        "PROFILE_DOXY_QC", "PROFILE_BBP700_QC", "PRES_QC", "PRES",
        "PRES_ADJUSTED", "PROFILE_PSAL_QC", "CHLA_QC", "CHLA_ADJUSTED",
        "CHLA_ADJUSTED_ERROR", "DOXY", "DOXY_QC", "DOXY_ADJUSTED",
        "DOXY_ADJUSTED_QC", "DOXY_ADJUSTED_ERROR", "PSAL", "PSAL_dPRES",
        "PSAL_ADJUSTED", "PSAL_ADJUSTED_QC", "TEMP", "TEMP_QC", "TEMP_dPRES",
        "TEMP_ADJUSTED", "TEMP_ADJUSTED_QC", "TEMP_ADJUSTED_ERROR"
      ),
      format = "dataframe"
    )
    
    # Preprocess data: calculate additional oceanographic parameters
    float_data <- float_data %>%
      filter(!is.na(DOXY)) %>%
      group_by(CYCLE_NUMBER) %>%
      mutate(
        ABS_SAL = gsw::gsw_SA_from_SP(
          SP = PSAL_ADJUSTED,
          p = PRES_ADJUSTED,
          longitude = first(LONGITUDE),
          latitude = first(LATITUDE)
        ),
        CONS_TEMP = gsw::gsw_CT_from_t(
          SA = ABS_SAL,
          t = TEMP_ADJUSTED,
          p = PRES_ADJUSTED
        ),
        SAT_DOXY = gsw_O2sol(
          SA = ABS_SAL,
          CT = CONS_TEMP,
          p = PRES_ADJUSTED,
          longitude = first(LONGITUDE),
          latitude = first(LATITUDE)
        ),
        AOU = SAT_DOXY - DOXY_ADJUSTED,
        SIGMA0 = gsw::gsw_sigma0(SA = ABS_SAL, CT = CONS_TEMP)
      )
    
    # Downscale data without outlier detection
    cycles_list <- float_data %>%
      group_by(CYCLE_NUMBER) %>%
      group_split()
    
    
    # Downscale Vertical resolution to 20
    
    downscaled_data_list_20m <- lapply(cycles_list, downscale_data_fun_wo_out,bin_width = 20)
    downscaled_data_20m <- bind_rows(downscaled_data_list_20m)
    residuals_data_20m <- downscaled_data_20m %>%
      group_by(CYCLE_NUMBER) %>%
      group_modify(~ .x %>%
                     select(PRES_ADJUSTED, AOU, ABS_SAL, LATITUDE, LONGITUDE, TIME) %>%
                     pivot_longer(
                       cols = c(AOU, ABS_SAL),
                       names_to = "VAR",
                       values_to = "VALUE"
                     ) %>%
                     group_by(VAR) %>%
                     mutate(
                       MA_3 = rollmean(VALUE, 3, fill = NA),
                       TM_9 = rollapply(VALUE, 9,
                                        function(x) {
                                          x_subset <- x[
                                            x >= quantile(x, 0.2, na.rm = TRUE) &
                                              x <= quantile(x, 0.8, na.rm = TRUE)
                                          ]
                                          if (length(x_subset) > 0) {
                                            mean(x_subset, na.rm = TRUE)
                                          } else {
                                            NA
                                          }
                                        },
                                        fill = NA
                       ),
                       ROB_RES = MA_3 - TM_9,
                       ROB_RES_RAW = VALUE - TM_9
                     ) %>%
                     mutate(
                       IQRN = IQR(ROB_RES_RAW, na.rm = TRUE) / 1.349,
                       MEDIAN_RES = median(ROB_RES_RAW[ROB_RES_RAW != 0], na.rm = TRUE)
                     ) %>%
                     mutate(
                       SCALE_RES_ROB = ifelse(
                         ROB_RES_RAW == 0,
                         0,
                         (ROB_RES_RAW - MEDIAN_RES) / IQRN
                       )
                     )
      ) %>%
      ungroup()
    
    residuals_data_wf_20m <- residuals_data_20m %>%
      group_by(CYCLE_NUMBER) %>%
      select(PRES_ADJUSTED, SCALE_RES_ROB,VALUE, VAR, CYCLE_NUMBER, LONGITUDE,
             LATITUDE, TIME) %>%
      pivot_wider(names_from = VAR, values_from = c(SCALE_RES_ROB,VALUE))
    
    
    # outlier calculation is on 40 m
    downscaled_data_list <- lapply(cycles_list, downscale_data_fun_wo_out,bin_width = 40)
    downscaled_data <- bind_rows(downscaled_data_list)
    
    # Prepare data for residuals calculation
    residuals_data <- downscaled_data %>%
      group_by(CYCLE_NUMBER) %>%
      group_modify(~ .x %>%
                     select(PRES_ADJUSTED, AOU, ABS_SAL, LATITUDE, LONGITUDE, TIME) %>%
                     pivot_longer(
                       cols = c(AOU, ABS_SAL),
                       names_to = "VAR",
                       values_to = "VALUE"
                     ) %>%
                     group_by(VAR) %>%
                     mutate(
                       MA_3 = rollmean(VALUE, 3, fill = NA),
                       TM_9 = rollapply(VALUE, 9,
                                        function(x) {
                                          x_subset <- x[
                                            x >= quantile(x, 0.2, na.rm = TRUE) &
                                              x <= quantile(x, 0.8, na.rm = TRUE)
                                          ]
                                          if (length(x_subset) > 0) {
                                            mean(x_subset, na.rm = TRUE)
                                          } else {
                                            NA
                                          }
                                        },
                                        fill = NA
                       ),
                       ROB_RES = MA_3 - TM_9,
                       ROB_RES_RAW = VALUE - TM_9
                     ) %>%
                     mutate(
                       IQRN = IQR(ROB_RES_RAW, na.rm = TRUE) / 1.349,
                       MEDIAN_RES = median(ROB_RES_RAW[ROB_RES_RAW != 0], na.rm = TRUE)
                     ) %>%
                     mutate(
                       SCALE_RES_ROB = ifelse(
                         ROB_RES_RAW == 0,
                         0,
                         (ROB_RES_RAW - MEDIAN_RES) / IQRN
                       )
                     )
      ) %>%
      ungroup()
    
    
    
    # Detect outliers in downscaled data and pivot to wide format
    residuals_data_wf <- residuals_data %>%
      group_by(CYCLE_NUMBER) %>%
      select(PRES_ADJUSTED, SCALE_RES_ROB,VALUE, VAR, CYCLE_NUMBER, LONGITUDE,
             LATITUDE, TIME) %>%
      pivot_wider(names_from = VAR, values_from = c(SCALE_RES_ROB,VALUE)) %>%
      mutate(OUT_S = ifelse(abs(SCALE_RES_ROB_AOU) > cutoff & abs(SCALE_RES_ROB_ABS_SAL) > cutoff
                            & SCALE_RES_ROB_AOU < 0, 1, 0))
    
    
    # Identify potential eddy events
    potential_eddy_events <- residuals_data_wf %>%
      filter(OUT_S == 1) %>%
      select(CYCLE_NUMBER, PRES_ADJUSTED, LATITUDE,
             LONGITUDE, TIME,SCALE_RES_ROB_ABS_SAL,
             SCALE_RES_ROB_AOU,VALUE_ABS_SAL,VALUE_AOU) %>%
      unique() %>%
      filter(PRES_ADJUSTED >= 200, PRES_ADJUSTED <= 1000)
    
    
    # Initialize columns to store the check results
    potential_eddy_events$AOU_gradient_sign_change <- NA
    potential_eddy_events$AOU_second_derivative_check <- NA
    potential_eddy_events$ABS_SAL_gradient_sign_change <- NA
    potential_eddy_events$ABS_SAL_second_derivative_check <- NA
    
    # Loop over each event
    for (i in seq_len(nrow(potential_eddy_events))) {
      cycle_num <- potential_eddy_events$CYCLE_NUMBER[i]
      pres_level <- potential_eddy_events$PRES_ADJUSTED[i]
      
      # Get the profile data at 20m resolution for the current cycle
      profile_data <- residuals_data_wf_20m %>% filter(CYCLE_NUMBER == cycle_num)
      
      # Perform checks for AOU
      checks_aou <- perform_checks(profile_data, target_level = pres_level,
                                   variable_name = "VALUE_AOU",
                                   second_deriv = 0.001,window=60)
      potential_eddy_events$AOU_gradient_sign_change[i] <- checks_aou$gradient_sign_change
      potential_eddy_events$AOU_second_derivative_check[i] <- checks_aou$second_derivative_check
      
      # Perform checks for ABS_SAL
      checks_abs_sal <- perform_checks(profile_data, target_level = pres_level,
                                       variable_name = "VALUE_ABS_SAL",second_deriv = 0.001,
                                       window=60)
      potential_eddy_events$ABS_SAL_gradient_sign_change[i] <- checks_abs_sal$gradient_sign_change
      potential_eddy_events$ABS_SAL_second_derivative_check[i] <- checks_abs_sal$second_derivative_check
    }
    
    
    filtered_events <- potential_eddy_events %>%
      filter(
        AOU_gradient_sign_change == TRUE,
        ABS_SAL_gradient_sign_change == TRUE)
    

    
   ##########################
    # PLOTTING ROUTINE #####
    ########################
    # Get all unique cycle numbers from your data
    all_cycles <- unique(residuals_data_wf_20m$CYCLE_NUMBER)
    
    # Number of non-anomalous profiles to sample per float
    num_samples <- 40
    
    # Identify cycles without detected anomalies
    non_anomalous_cycles <- setdiff(all_cycles, filtered_events$CYCLE_NUMBER)
    
    # Randomly sample non-anomalous cycles
    set.seed(123)  # For reproducibility
    sampled_non_anomalous_cycles <- sample(non_anomalous_cycles, min(num_samples, length(non_anomalous_cycles)))
    
    # Combine cycles with anomalies and sampled non-anomalous cycles
    #cycles_to_plot <- c(filtered_events$CYCLE_NUMBER, sampled_non_anomalous_cycles)
    
    # If you only want to plot non detected cycles :
    
    cycles_to_plot <- non_anomalous_cycles
    
    # Initialize empty lists for plots
    prof_plot <- list()
    res_plot <- list()
    list_plots <- list()
    
    # Plotting profiles and residuals
    for (i in seq_along(cycles_to_plot)) {
      current_cycle <- cycles_to_plot[i]
      current_data_20m <- residuals_data_20m %>% filter(CYCLE_NUMBER == current_cycle)
      current_data_40m <- residuals_data %>% filter(CYCLE_NUMBER == current_cycle)

      
      # Check if the current cycle has a detected anomaly
      anomaly_detected <- current_cycle %in% filtered_events$CYCLE_NUMBER
      
      # If an anomaly is detected, get the corresponding pressure level
      if (anomaly_detected) {
        current_eddy <- filtered_events %>% filter(CYCLE_NUMBER == current_cycle)
        pres_level <- current_eddy$PRES_ADJUSTED
      } else {
        pres_level <- NA
      }
      
      # Plotting profiles
      prof_plot[[i]] <- current_data_20m %>%
        ggplot(aes(x = PRES_ADJUSTED, y = VALUE)) +
        facet_grid(. ~ VAR, scales = "free") +
        coord_flip() +
        scale_x_reverse(limits = c(900, 0), breaks = seq(0, 900, by = 40)) +
        geom_line(aes(y = VALUE, color = "Observed Values")) +
        geom_line(aes(y = TM_9, color = "Trimmed Mean (k=9)")) +
        theme_bw() +
        labs(x = "Adjusted Pressure (dbar)", y = "") +
        theme(legend.position = "bottom")
      
      # Add vertical line if an anomaly is detected
      if (anomaly_detected) {
        prof_plot[[i]] <- prof_plot[[i]] +
          geom_vline(xintercept = pres_level, color = "darkgreen", alpha = 0.3, size = 1)
      }
      
      # Prepare data for residuals plot
      hline_data <- data.frame(
        VAR = c("AOU", "AOU", "ABS_SAL", "ABS_SAL"),
        hline = c(-2, 2, -2, 2),
        label = c("-2 sigma", "+2 sigma", "-2 sigma", "+2 sigma")
      )
      
      # Plotting residuals
      res_plot[[i]] <- current_data_40m %>%
        ggplot(aes(x = PRES_ADJUSTED, y = SCALE_RES_ROB)) +
        scale_x_reverse(limits = c(1200, 0), breaks = seq(0, 1200, by = 40)) +
        facet_grid(. ~ VAR, scales = "free") +
        coord_flip() +
        theme_bw() +
        geom_point() +
        labs(x = "Adjusted Pressure (dbar)", y = "") +
        geom_hline(data = hline_data, aes(yintercept = hline, color = label), inherit.aes = TRUE) +
        scale_color_manual(name = "Threshold", values = c("-2 sigma" = "blue", "+2 sigma" = "blue")) +
        theme(legend.position = "bottom")
      
      # Add vertical line if an anomaly is detected
      if (anomaly_detected) {
        res_plot[[i]] <- res_plot[[i]] +
          geom_vline(xintercept = pres_level, color = "darkgreen", alpha = 0.3, size = 1)
      }
      
      # Annotation text
      annotation_text <- paste(
        "Cycle Number:", current_cycle,
        "\nFloat ID (WMO):", wmo,
        "\nLongitude:", current_data_20m$LONGITUDE[1],
        "\nLatitude:", current_data_20m$LATITUDE[1],
        "\nTime:", format(as.POSIXct(current_data_20m$TIME[1], origin = "1970-01-01"), "%Y-%m-%d"),
        "\nAnomaly Detected:", ifelse(anomaly_detected, "Yes", "No")
      )
      
      # Combine plots with annotation
      combined_plot <- ggarrange(
        prof_plot[[i]],
        res_plot[[i]],
        common.legend = FALSE,
        legend = "bottom",
        nrow = 2
      )
      combined_plot <- annotate_figure(
        combined_plot,
        top = text_grob(annotation_text, face = "bold", size = 10)
      )
      
      list_plots[[i]] <- combined_plot
    }
      
      
    # Save plots to files
    if (length(list_plots) > 0) {
      dir <- paste0("/data/GLOBARGO/figures/EddySubductionFiguresSensitivityEstimation/", wmo)
      if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
      }
      
      for (k in seq_along(list_plots)) {
        cycle_number <- cycles_to_plot[k]
        file_name <- paste0(dir, "/", wmo, "_plot_cycle_", cycle_number, ".png")
        ggsave(file_name, list_plots[[k]], width = 10, height = 10)
      }
    }
    filtered_events$WMO <- wmo
    # Store the eddy events in the list
    detected_events_list[[j]] <- filtered_events
    
    
  }, silent = TRUE)
}    
    