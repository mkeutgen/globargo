# From 

setwd("~/Documents/GLOBARGO")
library(tidyverse)
library(robustbase)
library(gsw)
library(zoo)
library(oce)
library(ggpubr)
wmo <- 5904105
data_df = load_float_data( float_ids= wmo, # specify WMO number
                           variables=c(
                             "DATA_TYPE"
                             ,"PLATFORM_NUMBER"
                             ,"BBP700"
                             ,"BBP700_dPRES"
                             ,"BBP700_ADJUSTED_QC"
                             ,"LATITUDE"
                             ,"LONGITUDE"
                             ,"PROFILE_TEMP_QC"                              
                             ,"PROFILE_DOXY_QC"                                
                             ,"PROFILE_BBP700_QC"                                                
                             ,"PRES_QC"            
                             ,"PRES"  
                             ,"PRES_ADJUSTED"      
                             ,"PROFILE_PSAL_QC"
                             ,"CHLA_QC"  
                             ,"CHLA_ADJUSTED"
                             ,"CHLA_ADJUSTED_ERROR"
                             ,"DOXY"
                             ,"DOXY_QC"
                             ,"DOXY_ADJUSTED"
                             ,"DOXY_ADJUSTED_QC"
                             ,"DOXY_ADJUSTED_ERROR"
                             ,"PSAL"
                             ,"PSAL_dPRES"
                             ,"PSAL_ADJUSTED"
                             ,"PSAL_ADJUSTED_QC"
                             ,"TEMP"                            
                             ,"TEMP_QC"                        
                             ,"TEMP_dPRES"
                             ,"TEMP_ADJUSTED"                  
                             ,"TEMP_ADJUSTED_QC"
                             ,"TEMP_ADJUSTED_ERROR"), # specify variables,
                           format="dataframe" # specify format;  
)



# Calculations of residuals
data_df <- data_df %>% filter(!is.na(DOXY)) %>% group_by(CYCLE_NUMBER) %>%
  mutate(SPIC = swSpice(salinity = PSAL_ADJUSTED,
                        temperature = TEMP_ADJUSTED,
                        latitude = first(LATITUDE),
                        longitude = first(LONGITUDE),
                        eos = "unesco"),
         ABS_SAL = gsw::gsw_SA_from_SP(SP = PSAL_ADJUSTED, p = PRES_ADJUSTED,
                                       longitude = first(LONGITUDE),
                                       latitude = first(LATITUDE)),
         CONS_TEMP =  gsw::gsw_CT_from_t(SA = ABS_SAL, t = TEMP_ADJUSTED,
                                         p = PRES_ADJUSTED),
         SAT_DOXY = gsw_O2sol(SA= ABS_SAL,CT= CONS_TEMP,p= PRES_ADJUSTED,
                              longitude= first(LONGITUDE),latitude= first(LATITUDE)),
         AOU =  SAT_DOXY - DOXY_ADJUSTED,
         SIGMA0 = gsw::gsw_sigma0(SA = ABS_SAL,CT=CONS_TEMP)
  )



# Downscale residuals  
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

# Split
list_of_tibbles <- data_df  %>%
  group_by(CYCLE_NUMBER) %>%
  group_split()

# For Carbon Detection, we first apply a moving median (k = 3) filter to filter small anomalies
list_of_tibbles <- lapply(list_of_tibbles, function(df) {
  df$BBP700_ADJUSTEDORIG <- df$BBP700_ADJUSTED
  df$BBP700_ADJUSTED <- rollmedian(df$BBP700_ADJUSTED, k = 3, fill = NA)
  return(df)
})    
# Downscale data to a uniform resolution of 20 meters then compute residuals at 20 meters

downscaled_ds_list <- lapply(list_of_tibbles,downscale_data_fun_wo_out)

df <- downscaled_ds_list %>% bind_rows()
# matrix of residuals
A <- df %>%
  # Filter out missing DOXY values (not needed in this implementation since we take downscaled data as input)
  # filter(!is.na(DOXY)) %>%
  # Group by CYCLE_NUMBER
  group_by(CYCLE_NUMBER) %>%
  # Apply the computations for each unique CYCLE_NUMBER
  group_modify(~ .x %>%
                 dplyr::select(PRES_ADJUSTED, AOU, SPIC,BBP700_ADJUSTED, LATITUDE, LONGITUDE, TIME) %>%
                 pivot_longer(!c(PRES_ADJUSTED, LATITUDE, LONGITUDE, TIME), names_to = "VAR", values_to = "VALUE") %>%
                 group_by(VAR) %>%
                 mutate(
                   MA_3 = rollmean(VALUE, 3, fill = NA),
                   TM_11 = rollapply(VALUE, 11,
                                     function(x) {
                                       x_subset <- x[x >= quantile(x, 0.2, na.rm = TRUE) & x <= quantile(x, 0.8, na.rm = TRUE)]
                                       if (length(x_subset) > 0) {
                                         mean(x_subset, na.rm = TRUE)
                                       } else {
                                         NA
                                       }
                                     },
                                     fill = NA),
                   TM_9 = rollapply(VALUE, 9,
                                    function(x) {
                                      x_subset <- x[x >= quantile(x, 0.2, na.rm = TRUE) & x <= quantile(x, 0.8, na.rm = TRUE)]
                                      if (length(x_subset) > 0) {
                                        mean(x_subset, na.rm = TRUE)
                                      } else {
                                        NA
                                      }
                                    },
                                    fill = NA),
                   MM_11 = rollmedian(VALUE, 11, fill = NA),
                   ROB.RES = MA_3 - TM_9,
                   ROB.RES.RAW = VALUE - TM_9
                 ) %>%
                 mutate(
                   IQRN = IQR(ROB.RES.RAW, na.rm = TRUE) / 1.349,
                   MEDIAN_RES = median(ROB.RES.RAW[ROB.RES.RAW != 0], na.rm = TRUE)
                 ) %>%
                 mutate(
                   SCALE.RES.ROB = ifelse(ROB.RES.RAW == 0, 0,
                                          (ROB.RES.RAW - MEDIAN_RES) / IQRN)
                 )
  ) %>%
  ungroup()
# Check the result




A$TIME %>% unique()


# Interesting time points 2013-07-15" "2013-07-20" "2013-07-25


current_data <- A %>% filter(TIME == "2013-07-15")
current_eddy <- carb_eddy.id %>% filter(CYCLE_NUMBER == current_cycle)

current_data %>%
  ggplot(aes(x = PRES_ADJUSTED, y = VALUE)) +
  facet_grid(. ~ VAR, scales = "free") +
  coord_flip() +
  scale_x_reverse(limits = c(900, 0), breaks = seq(0, 900, by = 40)) +
  geom_line(aes(y = VALUE, color = "observed values")) +
  geom_point(aes(y = VALUE, color = "observed values"), size = .3) +
  geom_point(aes(y = TM_11, color = "Trimmed Mean (k=11)"), size = .3) +
  geom_line(aes(y = TM_11, color = "Trimmed Mean (k=11)")) +
  theme_bw() +
  labs(x = "Adjusted pressure (dbar)", y = "") +
  theme(legend.position = "bottom")

df <- current_data
df.ds <- downscale_data_fun(df, b = 40)
df.ds <- df.ds %>% ungroup() %>% select(AOU, SPIC, BBP700_ADJUSTED, PRES_ADJUSTED) %>%
  pivot_longer(cols = !PRES_ADJUSTED, names_to = "VAR", values_to = "VALUE")

hline_data <- data.frame(
  VAR = c("AOU", "AOU", "SPIC", "SPIC", "BBP700_ADJUSTED", "BBP700_ADJUSTED"),
  hline = c(-3, 3, -3, 3, -1, 1),
  label = c("-3 sigma", "+3 sigma", "-3 sigma", "+3 sigma", "-1 sigma", "+1 sigma")
)

res.plot[[i]] <- df %>%
  ggplot(aes(x = PRES_ADJUSTED, y = SCALE.RES.ROB)) +
  scale_x_reverse(limits = c(900, 0), breaks = seq(0, 900, by = 40)) +
  facet_grid(. ~ VAR, scales = "free") +
  coord_flip() +
  theme_bw() +
  geom_point(data = df.ds, aes(x = PRES_ADJUSTED, y = VALUE, color = "downscaled residuals to 40 meters")) +
  labs(x = "Adjusted pressure (dbar)", y = "") +
  geom_hline(data = hline_data, aes(yintercept = hline, color = label), inherit.aes = TRUE) +
  scale_color_manual(name = "Threshold", values = c(
    "-3 sigma" = "red",
    "+3 sigma" = "red",
    "-1 sigma" = "blue",
    "+1 sigma" = "blue"
  )) + theme(legend.position = "bottom") +
  geom_vline(
    xintercept = current_eddy$PRES_ADJUSTED[current_eddy$OUT.S == 1],
    color = "darkgreen",
    alpha = .3,
    size = 1
  )
