## Run a single-pass detection using the mock `load_float_data()` implementation.
# This script makes minimal changes to avoid touching the original detection scripts.

# Load libraries used by the detection scripts
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
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Source mock loader and annex functions
source("mock_load_float_data.R")
source("detection_algorithm_physical_subduction_annex_fun.R")

# Choose a single mock WMO
wmolist <- c("MOCK001")

# Minimal loop adapted from detection_algorithm_physical_subduction.R
detected_events_list <- list()
cutoff <- 1.96
resolution <- 40
window <- 60

for (j in seq_along(wmolist)) {
  wmo <- wmolist[j]
  try({
    float_data <- load_float_data(
      float_ids = wmo,
      variables = c("PRES", "PRES_ADJUSTED", "PSAL_ADJUSTED", "TEMP_ADJUSTED", "DOXY_ADJUSTED", "BBP700_ADJUSTED", "LATITUDE", "LONGITUDE"),
      format = "dataframe"
    )

    float_data <- float_data %>%
      filter(!is.na(DOXY_ADJUSTED)) %>%
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

    cycles_list <- float_data %>% group_by(CYCLE_NUMBER) %>% group_split()
    downscaled_data_list <- lapply(cycles_list, downscale_data_fun_wo_out, bin_width = 40)
    downscaled_data <- bind_rows(downscaled_data_list)

    residuals_data <- downscaled_data %>%
      group_by(CYCLE_NUMBER) %>%
      group_modify(~ .x %>%
                     select(PRES_ADJUSTED, AOU, ABS_SAL, LATITUDE, LONGITUDE, TIME) %>%
                     pivot_longer(cols = c(AOU, ABS_SAL), names_to = "VAR", values_to = "VALUE") %>%
                     group_by(VAR) %>%
                     mutate(
                       MA_3 = rollmean(VALUE, 3, fill = NA),
                       TM_9 = rollapply(VALUE, 9, function(x) { x_subset <- x[x >= quantile(x, 0.2, na.rm = TRUE) & x <= quantile(x, 0.8, na.rm = TRUE)]; if (length(x_subset) > 0) mean(x_subset, na.rm = TRUE) else NA }, fill = NA),
                       ROB_RES = MA_3 - TM_9,
                       ROB_RES_RAW = VALUE - TM_9
                     ) %>%
                     mutate(
                       IQRN = IQR(ROB_RES_RAW, na.rm = TRUE) / 1.349,
                       MEDIAN_RES = median(ROB_RES_RAW[ROB_RES_RAW != 0], na.rm = TRUE)
                     ) %>%
                     mutate(
                       SCALE_RES_ROB = ifelse(ROB_RES_RAW == 0, 0, (ROB_RES_RAW - MEDIAN_RES) / IQRN)
                     )
      ) %>% ungroup()

    residuals_data_wf <- residuals_data %>% group_by(CYCLE_NUMBER) %>% select(PRES_ADJUSTED, SCALE_RES_ROB, VALUE, VAR, CYCLE_NUMBER, LONGITUDE, LATITUDE, TIME) %>% pivot_wider(names_from = VAR, values_from = c(SCALE_RES_ROB, VALUE))

    residuals_data_wf <- residuals_data_wf %>% mutate(OUT_S = ifelse(abs(SCALE_RES_ROB_AOU) > cutoff & abs(SCALE_RES_ROB_ABS_SAL) > cutoff & SCALE_RES_ROB_AOU < 0, 1, 0))

    potential_eddy_events <- residuals_data_wf %>% filter(OUT_S == 1) %>% select(CYCLE_NUMBER, PRES_ADJUSTED, LATITUDE, LONGITUDE, TIME, SCALE_RES_ROB_ABS_SAL, SCALE_RES_ROB_AOU, VALUE_ABS_SAL, VALUE_AOU) %>% unique() %>% filter(PRES_ADJUSTED >= 200, PRES_ADJUSTED <= 1000)

    potential_eddy_events$AOU_gradient_sign_change <- NA
    potential_eddy_events$AOU_second_derivative_check <- NA
    potential_eddy_events$ABS_SAL_gradient_sign_change <- NA
    potential_eddy_events$ABS_SAL_second_derivative_check <- NA

    for (i in seq_len(nrow(potential_eddy_events))) {
      cycle_num <- potential_eddy_events$CYCLE_NUMBER[i]
      pres_level <- potential_eddy_events$PRES_ADJUSTED[i]
      profile_data <- residuals_data %>% filter(CYCLE_NUMBER == cycle_num)
      checks_aou <- perform_checks(profile_data, target_level = pres_level, variable_name = "VALUE_AOU", second_deriv = 0.001, window = 60)
      potential_eddy_events$AOU_gradient_sign_change[i] <- checks_aou$gradient_sign_change
      potential_eddy_events$AOU_second_derivative_check[i] <- checks_aou$second_derivative_check
      checks_abs_sal <- perform_checks(profile_data, target_level = pres_level, variable_name = "VALUE_ABS_SAL", second_deriv = 0.001, window = 60)
      potential_eddy_events$ABS_SAL_gradient_sign_change[i] <- checks_abs_sal$gradient_sign_change
      potential_eddy_events$ABS_SAL_second_derivative_check[i] <- checks_abs_sal$second_derivative_check
    }

    filtered_events <- potential_eddy_events %>% filter(AOU_gradient_sign_change == TRUE, ABS_SAL_gradient_sign_change == TRUE)
    filtered_events$WMO <- wmo
    detected_events_list[[j]] <- filtered_events
  }, silent = TRUE)
}

detected.events.df <- detected_events_list %>% bind_rows()
if (nrow(detected.events.df) == 0) {
  message("No events detected in mock run (this is expected with synthetic data).")
} else {
  write_csv(detected.events.df, "../data/detected_physical_subd_events_mock.csv")
  message("Wrote mock results to ../data/detected_physical_subd_events_mock.csv")
}
