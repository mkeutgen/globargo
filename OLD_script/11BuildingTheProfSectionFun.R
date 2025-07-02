# Testing with Llort Float
# TWO FUNCTIONS, ONE TO BUILD THE DATAFRAME WITH SPIC AOU AND BBP
# ONE TO PLOT
library(akima)
library(ggplot2)
library(ggpubr)


WMO <- 5904677
#   plot_mld           : plot mixed layer depth, using either a 
#                        temperature criterion (mld=1) or a density
#                        criterion (mld=2); default: 0=off




unique_sub_loc <- read_csv("/data/GLOBARGO/data/unique_sub_loc.csv")
det_event_cat1 <- read_csv("/data/GLOBARGO/data/detected_events_unique_with_carbon_cat1.csv")

det_event_cat1 %>% filter(WMO == 5904677)

anomalies_det <- det_event_cat1 %>% filter(WMO == 5904677) %>% filter(CYCLE_NUMBER %in% c(25,27)) %>% select(PRES_ADJUSTED,TIME)
anomalies_det$TIME <- as.POSIXct(anomalies_det$TIME)

variables <- c("DOXY","BBP_700","PSAL","NITRATE")
float_profs <- NULL
plot_isopyc=1
plot_mld=1
raw="no"
obs="off"
qc_flags=0:9

# Set the max depth to 1000 meters

max_depth <- 1000

# Downscaling function
downscale_data_fun_wo_out <- function(df,b=20) {
  # Select and pivot data
  data <- df %>%
    select(PRES_ADJUSTED, AOU,BBP700_ADJUSTED, SPIC,CYCLE_NUMBER,LONGITUDE,LATITUDE,TIME,BBP700_ADJUSTED)
  
  # Determine bin edges
  bin_width <- b
  pressure_range <- range(data$PRES_ADJUSTED, na.rm = TRUE)
  bins <- seq(from = floor(pressure_range[1]/bin_width)*bin_width,
              to = ceiling(pressure_range[2]/bin_width)*bin_width,
              by = bin_width)
  
  # Assign data points to bins
  data$bin <- cut(data$PRES_ADJUSTED, breaks = bins, include.lowest = TRUE, labels = FALSE)
  
  # Calculate mean for each bin
  downscaled_data <- data %>%
    group_by(bin) %>%
    mutate(across(-c(PRES_ADJUSTED,LATITUDE,LONGITUDE,CYCLE_NUMBER,TIME), mean, na.rm = TRUE))
  
  # Add bin's central value
  downscaled_data$PRES_ADJUSTED <- (bins[downscaled_data$bin] + bins[downscaled_data$bin + 1]) / 2
  
  downscaled_data <- unique(downscaled_data)
  # Return downscaled data
  return(downscaled_data)
}


gen_df_fun <- function(WMO){
  
  
  # make sure Setting is initialized
  if (exists("Setting")==F) {
    initialize_argo()
  }
  
  # download Sprof files if necessary
  good_float_ids = download_multi_floats(WMO)
  
  loaded = load_float_data(good_float_ids,variables = NULL)
  Data <- loaded$Data
  # Compute MLD :
  Data <- calc_mld(Data,calc_mld_dens = 1)

  Datai = Data[[1]]
  Mdata = loaded$Mdata
  
  
  # vertical interpolation to depths with regular intervals
  #
  #  prs_res = 2
  #  
  #  Datai = depth_interp(Data[[floats[f]]], qc_flags, calc_dens=calc_dens, 
  #                       calc_mld_temp=(plot_mld==1), calc_mld_dens=(plot_mld==2),
  #                       prs_res=prs_res)
  
  # computation of spiciness (SPIC) and apparent oxygen utilization (AOU)
  # based on adjusted variables if available 
  
  SPIC <- matrix(NA, nrow = nrow(Datai$PRES), ncol = ncol(Datai$PRES))
  
  # Loop over each profile (column in the matrix)
  for (i in 1:ncol(Datai$PRES_ADJUSTED)) {
    
    # Extract salinity, temperature, latitude, and longitude for the current profile
    salinity <- if (all(is.na(Datai$PSAL_ADJUSTED[, i]) | is.nan(Datai$PSAL_ADJUSTED[, i])) ) {
      Datai$PSAL[, i]
    } else {
      Datai$PSAL_ADJUSTED[, i]
    }
    
    temperature <- if (all(is.na(Datai$TEMP_ADJUSTED[, i]) | is.nan(Datai$TEMP_ADJUSTED[, i]) ) ) {
      Datai$TEMP[, i]
    } else {
      Datai$TEMP_ADJUSTED[, i]
    }
    
    latitude <- Datai$LATITUDE[1, i]  # Latitude is constant for a given profile
    longitude <- Datai$LONGITUDE[1, i]  # Longitude is constant for a given profile
    
    
    # Compute spiciness for each depth level in the profile
    for (j in 1:nrow(Datai$PRES_ADJUSTED)) {
      SPIC[j, i] <- swSpice(salinity = salinity[j], temperature = temperature[j],
                            latitude = latitude, longitude = longitude, eos = "unesco")
    }
  }
  
  # Add SPIC to Datai list
  Datai$SPIC <- SPIC
  
  
  # Initialize an empty matrix to store AOU values
  AOU <- matrix(NA, nrow = nrow(Datai$PRES_ADJUSTED), ncol = ncol(Datai$PRES_ADJUSTED))
  
  # Loop over each profile (column in the matrix)
  for (i in 1:ncol(Datai$PRES_ADJUSTED)) {
    
    # Conditional handling of adjusted vs non-adjusted variables
    salinity <- if (all(is.na(Datai$PSAL_ADJUSTED[, i]) | is.nan(Datai$PSAL_ADJUSTED[, i])) ) {
      Datai$PSAL[, i]
    } else {
      Datai$PSAL_ADJUSTED[, i]
    }
    
    temperature <- if (all(is.na(Datai$TEMP_ADJUSTED[, i]) | is.nan(Datai$TEMP_ADJUSTED[, i]))) {
      Datai$TEMP[, i]
    } else {
      Datai$TEMP_ADJUSTED[, i]
    }
    
    pressure <- if (all(is.na(Datai$PRES_ADJUSTED[, i]) | is.nan(Datai$PRES_ADJUSTED[, i]))) {
      Datai$PRES[, i]
    } else {
      Datai$PRES_ADJUSTED[, i]
    }
    
    doxy <- if (all(is.na(Datai$DOXY_ADJUSTED[, i]) | is.nan(Datai$DOXY_ADJUSTED[, i]))) {
      Datai$DOXY[, i]
    } else {
      Datai$DOXY_ADJUSTED[, i]
    }
    
    latitude <- Datai$LATITUDE[1, i]  # Latitude is constant for a given profile
    longitude <- Datai$LONGITUDE[1, i]  # Longitude is constant for a given profile
    
    
    # Compute Absolute Salinity (ABS_SAL) and Conservative Temperature (CONS_TEMP) for each depth level
    ABS_SAL <- gsw_SA_from_SP(SP = salinity, p = pressure, longitude = longitude, latitude = latitude)
    CONS_TEMP <- gsw_CT_from_t(SA = ABS_SAL, t = temperature, p = pressure)
    
    # Compute Saturation Dissolved Oxygen (SAT_DOXY)
    SAT_DOXY <- gsw_O2sol(SA = ABS_SAL, CT = CONS_TEMP, p = pressure, longitude = longitude, latitude = latitude)
    
    # Compute AOU for each depth level
    AOU[, i] <- SAT_DOXY - doxy
  }
  
  # Add AOU to Datai list
  Datai$AOU <- AOU
  
  
  # Convert time to universal time coordinates UTC
  
  df = NULL
  for (name in names(Datai)) {
    if (name == "TIME") {
      df[[name]] = as.POSIXct(Datai[[name]], tz="UTC")
    } else {
      df[[name]] = as.vector(Datai[[name]])
    }
  }
  
  df = data.frame(df)
  
  # Downscaling AOU SPIC and BBP700 to uniform 20 meters vertical res 
  list_of_tibbles <- df  %>%
    group_by(CYCLE_NUMBER) %>%
    group_split()
  #  
  #  list_of_tibbles <- lapply(list_of_tibbles, function(df) {
  #    df$BBP700_ADJUSTEDORIG <- df$BBP700_ADJUSTED
  #    df$BBP700_ADJUSTED <- rollmedian(df$BBP700_ADJUSTED, k = 3, fill = NA)
  #    return(df)
  #  })    
  # Downscale data to a uniform resolution of 20 meters then compute residuals at 20 meters
  
  downscaled_ds_list <- lapply(list_of_tibbles,downscale_data_fun_wo_out)
  
  df <- downscaled_ds_list %>% bind_rows()
  
  
  # So 233 is unique number of dates in that float
  # Where does 851 comes from ? format="dataframe"
  list_of_tibbles <- df  %>%
    group_by(TIME) %>%
    group_split()
  list.df <- list()
  
  for (i in seq_along(list_of_tibbles)) {
    tibble_i <- list_of_tibbles[[i]]
    #tibble_i$BBP700_ADJUSTED <- rollmedian(tibble_i$BBP700_ADJUSTED, k = 3, fill = NA)
    # Assign xmin and xmax to each tibble
    tibble_i$xmin <- rep(full_xmin[i], nrow(tibble_i))
    tibble_i$xmax <- rep(full_xmax[i], nrow(tibble_i))
    
    # Replace the tibble in the list with the modified one
    list.df[[i]] <- tibble_i
  }
  
  
  
  df <- list.df %>% bind_rows()
  
  ### END OF DATAFRAME GEN
  
  return(df)
  
}

df <- gen_df_fun(WMO=5904105)


# Function to process data, interpolate, and plot results

generate_plots <- function(df, wmo, start_date, end_date, cycle_number, logbbp = FALSE) {
  anomalies_det <- det_event_cat1 %>% 
    filter(WMO == wmo) %>% 
    filter(CYCLE_NUMBER %in% cycle_number) %>% 
    select(PRES_ADJUSTED, TIME)
  anomalies_det$TIME <- as.POSIXct(anomalies_det$TIME)
  
  df_filtered <- df %>%
    filter(TIME >= as.POSIXct(start_date) & TIME <= as.POSIXct(end_date))
  
  df_filtered_clean <- df_filtered %>%
    filter(!is.na(PRES_ADJUSTED) & !is.na(AOU) & !is.na(SPIC) & !is.na(TIME))
  
  reference_time <- min(df_filtered_clean$TIME)
  df_filtered_clean$TIME_numeric <- as.numeric(difftime(df_filtered_clean$TIME, reference_time, units = "days"))
  
  interp_result_AOU <- with(df_filtered_clean, interp(
    x = TIME_numeric,
    y = PRES_ADJUSTED,
    z = AOU,
    xo = seq(min(TIME_numeric), max(TIME_numeric), length = 200),
    yo = seq(min(PRES_ADJUSTED), max(PRES_ADJUSTED), length = 200)
  ))
  
  interp_result_SPIC <- with(df_filtered_clean, interp(
    x = TIME_numeric,
    y = PRES_ADJUSTED,
    z = SPIC,
    xo = seq(min(TIME_numeric), max(TIME_numeric), length = 200),
    yo = seq(min(PRES_ADJUSTED), max(PRES_ADJUSTED), length = 200)
  ))
  
  interp_result_BBP700 <- with(df_filtered_clean, interp(
    x = TIME_numeric,
    y = PRES_ADJUSTED,
    z = BBP700_ADJUSTED,
    xo = seq(min(TIME_numeric), max(TIME_numeric), length = 200),
    yo = seq(min(PRES_ADJUSTED), max(PRES_ADJUSTED), length = 200)
  ))
  
  interp_df <- data.frame(
    TIME = as.POSIXct(reference_time + interp_result_AOU$x * 24 * 60 * 60),
    PRES_ADJUSTED = rep(interp_result_AOU$y, each = length(interp_result_AOU$x)),
    AOU = as.vector(interp_result_AOU$z),
    SPIC = as.vector(interp_result_SPIC$z),
    BBP700_ADJUSTED = as.vector(interp_result_BBP700$z)
  )
  
  # Extract unique time breaks
  time_breaks <- df_filtered_clean$TIME %>% unique()
  
  plot_AOU <- ggplot(interp_df, aes(x = TIME, y = PRES_ADJUSTED, fill = AOU)) +
    geom_tile() +
    scale_fill_viridis_c(name = "AOU (µmol/kg)") +
    theme_bw() +
    ggtitle("Apparent Oxygen Utilization (Linearly Interpolated)") + 
    scale_y_reverse(limits = c(max(interp_df$PRES_ADJUSTED), 0)) +
    scale_x_datetime(breaks = time_breaks, date_labels = "%d-%b-%Y") +
    geom_point(data = anomalies_det, aes(x = TIME, y = PRES_ADJUSTED), color = "red", size = 3, alpha = 0.5, inherit.aes = FALSE)
  
  plot_SPIC <- ggplot(interp_df, aes(x = TIME, y = PRES_ADJUSTED, fill = SPIC)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Spiciness (kg/m³)") +
    theme_bw() +
    ggtitle("Spiciness (Linearly Interpolated)") + 
    scale_y_reverse(limits = c(max(df_filtered$PRES_ADJUSTED), 0)) +
    scale_x_datetime(breaks = time_breaks, date_labels = "%d-%b-%Y") +
    geom_point(data = anomalies_det, aes(x = TIME, y = PRES_ADJUSTED), color = "red", size = 3, alpha = 0.5, inherit.aes = FALSE)
  
  if (logbbp == FALSE) {
    plot_BBP700 <- ggplot(interp_df, aes(x = TIME, y = PRES_ADJUSTED, fill = BBP700_ADJUSTED)) +
      geom_tile() +
      scale_fill_viridis_c(name = "BBP700 (m⁻¹)") +
      theme_bw() +
      ggtitle("BBP700_ADJUSTED (Linearly Interpolated)") + 
      scale_y_reverse(limits = c(max(df_filtered$PRES_ADJUSTED), 50)) +
      scale_x_datetime(breaks = time_breaks, date_labels = "%d-%b-%Y") +
      geom_point(data = anomalies_det, aes(x = TIME, y = PRES_ADJUSTED), color = "red", size = 3, alpha = 0.5, inherit.aes = FALSE)
  } else {
    plot_BBP700 <- ggplot(interp_df, aes(x = TIME, y = PRES_ADJUSTED, fill = log10(BBP700_ADJUSTED + 1))) +
      geom_tile() +
      scale_fill_viridis_c(name = "log10(BBP700 + 1)") +
      theme_bw() +
      ggtitle("BBP700_ADJUSTED (Log-Transformed, Interpolated)") + 
      scale_y_reverse(limits = c(max(df_filtered$PRES_ADJUSTED), 0)) +
      scale_x_datetime(breaks = time_breaks, date_labels = "%d-%b-%Y") +
      geom_point(data = anomalies_det, aes(x = TIME, y = PRES_ADJUSTED), color = "red", size = 3, alpha = 0.5, inherit.aes = FALSE)
  }
  
  combined_plot <- ggarrange(
    plot_AOU, plot_SPIC, plot_BBP700, 
    ncol = 1,
    align = "v"
  )
  
  return(combined_plot)
}

# Example usage:
df_5904677 <- gen_df_fun(5904677)

combined_plot_5904677 <- generate_plots(df_5904677,wmo = 5904677, start_date = "2016-09-01", end_date = "2016-12-31", cycle_number = c(25,27) )


# Exploring visually noteworthy anomalies : 
det_event_cat1 %>% filter(WMO == 5906312)


df_5906312 <- gen_df_fun(5906312)

combined_plot_5906312 <- generate_plots(df_5906312,wmo = 5906312, start_date = "2023-08-28", end_date = "2023-10-28", cycle_number = c(32))

#6901767
det_event_cat1 %>% filter(WMO == 6901767)

df_6901767 <- gen_df_fun(6901767)

combined_plot_6901767 <- generate_plots(df_5904677,wmo = 6901767, start_date = "2017-12-18", end_date = "2018-02-08", cycle_number = c(179) )

ggsave(combined_plot_6901767,filename = "/data/GLOBARGO/figures/FiguresForPublication/section_plot_6901767.png",width = 12,height = 14)
