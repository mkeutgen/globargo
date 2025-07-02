# Function to calculate the derivative of a variable with respect to pressure

compute_derivatives_vectorized <- function(data) {
  data <- data %>% arrange(PRES_ADJUSTED)
  
  data <- data %>%
    mutate(
      x_next = dplyr::lead(PRES_ADJUSTED),   # Next pressure level
      y_next = dplyr::lead(VALUE),           # Next value
      
      # Compute forward differences for the first derivative
      dVALUE = case_when(
        is.na(x_next) ~ NA_real_,            # No forward difference for the last point
        TRUE ~ (y_next - VALUE) / (x_next - PRES_ADJUSTED)   # Forward difference for other points
      ),
      
      # Auxiliary variables for second derivative computation
      x_prev = dplyr::lag(PRES_ADJUSTED),
      y_prev = dplyr::lag(VALUE),
      h_forward = x_next - PRES_ADJUSTED,
      h_backward = PRES_ADJUSTED - x_prev,
      term1 = (y_next - VALUE) / h_forward,
      term2 = (VALUE - y_prev) / h_backward,
      
      # Compute second derivative using centered differences
      d2VALUE = case_when(
        is.na(x_prev) | is.na(x_next) ~ NA_real_,           # No second derivative at endpoints
        TRUE ~ 2 * (term1 - term2) / (x_next - x_prev)      # Centered difference for second derivative
      )
    ) %>%
    select(-x_prev, -x_next, -y_prev, -y_next, -h_forward, -h_backward, -term1, -term2)
  
  return(data)
}


perform_checks <- function(profile_data, target_level, variable_name,second_deriv,window=100) {
  # Select the variable of interest
  data <- profile_data %>%
    select(CYCLE_NUMBER, PRES_ADJUSTED, VALUE = all_of(variable_name))
  
  # Remove rows with NA in VALUE
  data <- data %>% filter(!is.na(VALUE))
  
  # If there are not enough data points, return FALSE for both checks
  if (nrow(data) < 3) {
    return(list(
      gradient_sign_change = FALSE,
      derivative_close_to_zero = FALSE,
      second_derivative_check = FALSE
    ))
  }
  
  # Compute derivatives
  data <- compute_derivatives_vectorized(data)
  
  # Initialize flags
  gradient_sign_change <- FALSE
  derivative_close_to_zero <- FALSE
  second_derivative_check <- FALSE
  
  # Check for sign change in gradient within ± x= window meters
  subset_data <- data %>%
    filter(PRES_ADJUSTED >= (target_level - window) & PRES_ADJUSTED <= (target_level + window))
  
  # Check for sign change or close to zero in first derivative
  if (nrow(subset_data) >= 2) {
    # Get the signs of the first derivative
    signs <- sign(subset_data$dVALUE)
    # Check if there's a sign change
    sign_changes <- diff(signs)
    gradient_sign_change <- any(sign_changes != 0, na.rm = TRUE)
    
    # Check if any derivative is close to zero (e.g., |dVALUE| < 0.1)
    derivative_close_to_zero <- any(abs(subset_data$dVALUE) < 0.1, na.rm = TRUE)
  }
  
  # Check if second derivative at target level is larger than 0.01
  # Find the index of the point closest to the target level
  idx <- which.min(abs(data$PRES_ADJUSTED - target_level))
  second_derivative_value <- data$d2VALUE[idx]
  second_derivative_check <- !is.na(second_derivative_value) && (second_derivative_value > second_deriv)
  
  # Return the results as a list
  return(list(
    gradient_sign_change = gradient_sign_change,
    derivative_close_to_zero = derivative_close_to_zero,
    second_derivative_check = second_derivative_check
  ))
}





# Function to downscale data and detect outliers
downscale_data_fun <- function(df, bin_width = 20, cutoff = 1.96) {
  data <- df %>%
    select(PRES_ADJUSTED, SCALE_RES_ROB, VAR, CYCLE_NUMBER, LONGITUDE, LATITUDE, TIME) %>%
    pivot_wider(names_from = VAR, values_from = SCALE_RES_ROB)
  
  pressure_range <- range(data$PRES_ADJUSTED, na.rm = TRUE)
  bins <- seq(
    from = floor(pressure_range[1] / bin_width) * bin_width,
    to = ceiling(pressure_range[2] / bin_width) * bin_width,
    by = bin_width
  )
  
  data$bin <- cut(data$PRES_ADJUSTED, breaks = bins, include.lowest = TRUE, labels = FALSE)
  
  downscaled_data <- data %>%
    group_by(bin) %>%
    mutate(across(
      .cols = -c(PRES_ADJUSTED, LATITUDE, LONGITUDE, CYCLE_NUMBER, TIME),
      .fns = ~ mean(., na.rm = TRUE)
    ))
  
  downscaled_data$PRES_ADJUSTED <- (bins[downscaled_data$bin] + bins[downscaled_data$bin + 1]) / 2
  
  downscaled_data <- downscaled_data %>%
    mutate(OUT_S = ifelse(abs(AOU) > cutoff & abs(ABS_SAL) > cutoff & AOU < 0, 1, 0))
  
  return(downscaled_data)
}

# Function to downscale data without outlier detection
downscale_data_fun_wo_out <- function(df, bin_width = 40) {
  data <- df %>%
    select(PRES_ADJUSTED, AOU, ABS_SAL, CYCLE_NUMBER, LONGITUDE, LATITUDE, TIME)
  
  pressure_range <- range(data$PRES_ADJUSTED, na.rm = TRUE)
  bins <- seq(
    from = floor(pressure_range[1] / bin_width) * bin_width,
    to = ceiling(pressure_range[2] / bin_width) * bin_width,
    by = bin_width
  )
  
  data$bin <- cut(data$PRES_ADJUSTED, breaks = bins, include.lowest = TRUE, labels = FALSE)
  
  downscaled_data <- data %>%
    group_by(bin) %>%
    mutate(across(
      .cols = -c(PRES_ADJUSTED, LATITUDE, LONGITUDE, CYCLE_NUMBER, TIME),
      .fns = mean,
      na.rm = TRUE
    ))
  
  downscaled_data$PRES_ADJUSTED <- (bins[downscaled_data$bin] + bins[downscaled_data$bin + 1]) / 2
  
  downscaled_data <- unique(downscaled_data)
  return(downscaled_data)
}

# Function to compute mean ABS_SAL in the upper 50 meters
compute_mean_ABS_SAL_50m <- function(ABS_SAL_df) {
  mean_ABS_SAL <- ABS_SAL_df %>%
    filter(PRES_ADJUSTED <= 50) %>%
    summarize(mean_ABS_SAL = mean(ABS_SAL, na.rm = TRUE)) %>%
    pull(mean_ABS_SAL)
  
  return(mean_ABS_SAL)
}

# Function to compute mean ABS_SAL at min and max pressure levels around a target pressure
mean_ABS_SAL_at_min_max_levels <- function(data, target_pressure) {
  filtered_data <- data %>%
    filter(!is.na(ABS_SAL)) %>%
    filter(PRES_ADJUSTED >= (target_pressure - 100), PRES_ADJUSTED <= (target_pressure + 100))
  
  min_pressure <- min(filtered_data$PRES_ADJUSTED, na.rm = TRUE)
  max_pressure <- max(filtered_data$PRES_ADJUSTED, na.rm = TRUE)
  
  ABS_SAL_values <- filtered_data %>%
    filter(PRES_ADJUSTED %in% c(min_pressure, max_pressure)) %>%
    pull(ABS_SAL)
  
  mean_ABS_SAL <- mean(ABS_SAL_values, na.rm = TRUE)
  
  return(mean_ABS_SAL)
}
