# Carbon Quantification R Code #

# Code to quantify the magnitude of AOU and SPIC and BBP700 anomalies on Float 5904677
# Input : a profile from an argo float where we know that there is a correct anomaly at given pressure level
# Output : a tibble giving the integrated BBP700 anomaly and AOU anomaly
# (mumol/kg) along with space/time variables.
# Strategy : develop a prototype on Llort Float 5904677 profile 27
# Innovation : integrated anomaly
# Prereq : run Tutorial.R file from OneArgoR toolbox to be able to download float data
# !! File is wrong, integral of BBP_700 is done then mapped to POC but it is not mathematically correct!! 
# !!Check Carbon Quanti on All floats for correction
# Prototype 

wmo <- 5904677
cycle_number <- 27
critical_pres <- c(300,340)

# Functions
finiteCenteredDiff <- function(vec) {
  # Initialize the output vector of the same length as input, filled with NA
  diffVec <- rep(NA, length(vec))
  
  # Compute the finite centered difference for elements that have neighbors
  for(i in 2:(length(vec)-1)) {
    diffVec[i] <- (vec[i+1] - vec[i-1]) / 2
  }
  
  return(diffVec)
}

predict_AOU <- function(pressure,m,b) {
  m * pressure + b
}


na_replace_fun <- function(values,pres_level,upper.bound,lower.bound){
  output <- ifelse(pres_level > upper.bound,NA,
                   ifelse(pres_level < lower.bound,NA,values)
  )
  return(output)
}



# Function to group and average pressures within 100 meters
average_close_pressures <- function(pressures, threshold = 50) {
  # Create a data frame from the vector
  df <- tibble(PRES_ADJUSTED = pressures)
  
  # Group and average values within the threshold
  df_grouped <- df %>%
    arrange(PRES_ADJUSTED) %>%
    mutate(group = cumsum(c(1, diff(PRES_ADJUSTED) > threshold))) %>%
    group_by(group) %>%
    summarize(PRES_ADJUSTED = mean(PRES_ADJUSTED)) %>%
    ungroup() %>%
    select(PRES_ADJUSTED)
  
  return(df_grouped$PRES_ADJUSTED)
}



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
# Select only the cycle number which are anomalous
data_df <- data_df %>% filter(data_df$CYCLE_NUMBER %in% cycle_number) %>% ungroup()

# We know that 
#  | bin| PRES_ADJUSTED| CYCLE_NUMBER| OUT.T|
#  |---:|-------------:|------------:|-----:|
#  |   8|           300|           27|     1|
#  |   9|           340|           27|     1|

# Calculation of ancillary variables SPIC and AOU
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
  )



data_df$BBP700_ADJUSTED <- rollmean(data_df$BBP700_ADJUSTED,k = 7,na.pad = TRUE)

data.df.lf <- data_df %>%
  select(LATITUDE,LONGITUDE,TIME,CYCLE_NUMBER,PRES_ADJUSTED,AOU,SPIC,BBP700_ADJUSTED) %>%
  pivot_longer(cols = !c(LATITUDE,LONGITUDE,TIME,CYCLE_NUMBER,PRES_ADJUSTED))


# First safeguard, if points are close to one another (within 100 meters), choose one of them.
adjusted_critical_pres <- average_close_pressures(critical_pres)

ggplot(data.df.lf,aes(x=PRES_ADJUSTED,y=value))+
  facet_grid(.~name,scales = "free")+geom_line()+coord_flip()+
  scale_x_reverse(limits = c(900, 0), breaks = seq(0, 900, by = 40))+
  geom_vline(xintercept = adjusted_critical_pres)+theme_bw()

# Quantification of the anomalies using Chen's method,


# Then find the peak, in absolute value, in AOU, SPIC and BBP between -100 and 100 meters

df <- data_df %>% filter(!is.na(DOXY)) 

df <- df %>% select(SPIC,AOU,PRES_ADJUSTED,BBP700_ADJUSTED)


# Find max in 100 pres_level above :


find_max_AOU_above_below_peak <- function(df, peak_pres, range = 150) {
  # Filter the dataframe for the range below the peak
  df_below <- df %>%
    filter(PRES_ADJUSTED < peak_pres & PRES_ADJUSTED >= (peak_pres - range))
  
  # Find the row with the maximum AOU in this range below the peak
  max_AOU_row_below <- df_below %>%
    filter(AOU == max(AOU)) %>%
    slice(1) # In case there are multiple, just take the first
  
  # Filter the dataframe for the range above the peak
  df_above <- df %>%
    filter(PRES_ADJUSTED > peak_pres & PRES_ADJUSTED <= (peak_pres + range))
  
  # Find the row with the maximum AOU in this range above the peak
  max_AOU_row_above <- df_above %>%
    filter(AOU == max(AOU)) %>%
    slice(1) # In case there are multiple, just take the first
  
  # Return both pressure levels
  return(list(
    max_below = max_AOU_row_below %>% select(PRES_ADJUSTED,AOU,BBP700_ADJUSTED) %>% ungroup(),
    max_above = max_AOU_row_above %>% select(PRES_ADJUSTED,AOU,BBP700_ADJUSTED) %>% ungroup()
  ))
}
# BBP700 anomaly
find_min_BBP700_above_below_peak <- function(df, peak_pres, range = 150) {
  # Filter the dataframe for the range below the peak
  df_below <- df %>%
    filter(PRES_ADJUSTED < peak_pres & PRES_ADJUSTED >= (peak_pres - range))
  
  # Find the row with the min BBP700 in this range below the peak
  min_BBP700_ADJUSTED_row_below <- df_below %>%
    filter(BBP700_ADJUSTED == min(BBP700_ADJUSTED)) %>%
    slice(1) # In case there are multiple, just take the first
  
  # Filter the dataframe for the range above the peak
  df_above <- df %>%
    filter(PRES_ADJUSTED > peak_pres & PRES_ADJUSTED <= (peak_pres + range))
  
  # Find the row with the maximum AOU in this range above the peak
  min_BBP700_ADJUSTED_row_above <- df_above %>%
    filter(BBP700_ADJUSTED == min(BBP700_ADJUSTED)) %>%
    slice(1) # In case there are multiple, just take the first
  
  # Return both pressure levels
  return(list(
    min_below = min_BBP700_ADJUSTED_row_below %>% select(PRES_ADJUSTED,BBP700_ADJUSTED) %>% ungroup(),
    min_above = min_BBP700_ADJUSTED_row_above %>% select(PRES_ADJUSTED,BBP700_ADJUSTED) %>% ungroup()
  ))
}


peak.up_and_low_aou <- find_max_AOU_above_below_peak(df,peak_pres = adjusted_critical_pres)
peak.up_and_low_bbp700 <- find_min_BBP700_above_below_peak(df,adjusted_critical_pres)



# Extract coordinates of 2 maxima (below and above the detected minimum in AOU) 
pres_bl <- peak.up_and_low_aou$max_below$PRES_ADJUSTED
aou_bl <- peak.up_and_low_aou$max_below$AOU
pres_bbp_bl <- peak.up_and_low_bbp700$min_below$PRES_ADJUSTED
bbp_bl      <- peak.up_and_low_bbp700$min_below$BBP700_ADJUSTED

pres_up <- peak.up_and_low_aou$max_above$PRES_ADJUSTED
aou_up <- peak.up_and_low_aou$max_above$AOU

pres_bbp_up <- peak.up_and_low_bbp700$min_above$PRES_ADJUSTED
bbp_up      <- peak.up_and_low_bbp700$min_above$BBP700_ADJUSTED


# Calculate slope (m) and intercept (b) of the line
m_aou <- (aou_up - aou_bl) /
  (pres_up - pres_bl)
b_aou <- aou_bl - m_aou * pres_bl

m_bbp <- (bbp_up - bbp_bl) /
  (pres_bbp_up - pres_bbp_bl)

b_bbp <- bbp_bl - m_bbp * pres_bbp_bl



# Function to compute predicted AOU at a given pressure level


# Assuming df is your original data frame and you want to plot for this range
df$predicted_AOU <- sapply(df$PRES_ADJUSTED, predict_AOU,m=m_aou,b=b_aou)
df$predicted_BBP700_ADJUSTED <- sapply(df$PRES_ADJUSTED, predict_AOU,m=m_bbp,b=b_bbp)



df$predicted_AOU <- na_replace_fun(df$predicted_AOU,df$PRES_ADJUSTED,pres_up,pres_bl)
df$predicted_BBP700_ADJUSTED <- na_replace_fun(df$predicted_BBP700_ADJUSTED,df$PRES_ADJUSTED,pres_up,pres_bl)

ggplot(df,aes(x=PRES_ADJUSTED,y=AOU,color="AOU observed values"))+geom_line()+geom_point()+coord_flip()+
  scale_x_reverse()+
  geom_line(aes(y=predicted_AOU,color="linear interpolation of AOU following Chen's method (Chen et al, 2021)"))+
  theme_bw()+theme(legend.position="bottom")

ggplot(df,aes(x=PRES_ADJUSTED,y=BBP700_ADJUSTED,color="BBP700_ADJUSTED observed values"))+
  geom_line()+geom_point()+coord_flip()+
  scale_x_reverse()+
  geom_point(aes(y=predicted_BBP700_ADJUSTED,color="linear interpolation of BBP700 following Chen's method (Chen et al, 2021)"))+
  geom_line(aes(y=predicted_BBP700_ADJUSTED,color="linear interpolation of BBP700 following Chen's method (Chen et al, 2021)"))+
  theme_bw()+theme(legend.position="bottom")



# Compute the distances between adjacent pressure levels
distances <- diff(df$PRES_ADJUSTED)

# Calculate absolute differences between predicted and actual AOU for each segment
differences_aou <- abs(df$predicted_AOU - df$AOU)
differences_bbp <- abs(df$predicted_BBP700_ADJUSTED - df$BBP700_ADJUSTED)


avg_differences_aou <- (differences_aou[-length(differences_aou)] + differences_aou[-1]) / 2
avg_differences_bbp <- (differences_bbp[-length(differences_bbp)] + differences_bbp[-1]) / 2

trapezoid_areas_aou <- avg_differences_aou * distances
trapezoid_areas_bbp <- avg_differences_bbp * distances


# Sum up the areas of the trapezoids this gives a depth integrated bbp700/aou anomaly in (m'-1 sr^-1) \times m
total_area_aou <- sum(trapezoid_areas_aou, na.rm = TRUE)
total_area_bbp <- sum(trapezoid_areas_bbp, na.rm = TRUE)

# Relate it to carbon
# POC= 9.776 \times 10^4 \times (bbp_700)^1.166

# Define a function to compute POC from bbp_700
poc_from_bbp_700 <- function(bbp_700) {
  # Apply the given formula:
  9.776e4 * (bbp_700 ^ 1.166)
}
 
poc_from_bbp_700(total_area_bbp) # 2210.865 mg/m^2 day/hour without bbp_700 smoothing
# 1988.211 mg/m^2 
