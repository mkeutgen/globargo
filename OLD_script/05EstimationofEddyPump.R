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
library(R.matlab)

library(fitdistrplus)
library(poweRlaw)
# Resolve function conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")
setwd("/data/GLOBARGO/src/")
# Read the data
df_carbon_with_poc <- read_csv("data/df_carbon_subduction_anom_with_poc_fromgali.csv") %>% filter(integrated_poc > 0.01) 
df_carbon_with_poc
df_argo_clean <- read_csv("data/df_argo_loc.csv")

# Implementing Siegel's correction : 

fseq_df <- read_csv("co2_sequestration_50years.csv")
fseq_df$fseq %>% summary()

fseq_df <- fseq_df %>% na.omit()


# 1. Prepare your data frames
df_carbon <- df_carbon_with_poc %>%
  mutate(depth_meters = critical_pres) %>%  # Approximate pressure to depth
  select(LATITUDE, LONGITUDE, depth_meters, integrated_poc,TIME)

# 2. Create 3D coordinate matrices
coords_carbon <- df_carbon %>% 
  select(LATITUDE, LONGITUDE, depth_meters) %>% 
  as.matrix()

coords_fseq <- fseq_df %>% 
  select(latitude, longitude, depth_meters) %>% 
  as.matrix()

# 3. Find nearest neighbors (takes <1 second even for large datasets)
nn <- nn2(coords_fseq, coords_carbon, k = 1)

# 4. Add the matched fseq values
df_carbon$fseq_50yr <- fseq_df$fseq[nn$nn.idx]

df_carbon_with_poc$fseq_50yr <- df_carbon$fseq_50yr

# clock function
df_carbon_with_poc$MONTH <- month(df_carbon_with_poc$TIME)
pseudo_month <- function(date, lat) {
  ## -------- sanity checks ----------
  if (length(date) != length(lat))
    stop("`date` and `lat` must be the same length.")
  if (!inherits(date, "Date") && !inherits(date, "POSIXt"))
    stop("`date` must be Date/POSIXt; convert first with as.Date() if needed.")
  
  ## -------- shift by 6 months south of the equator ----------
  shifted <- date                    # initialise
  south   <- which(lat < 0)          # indices in S. Hemisphere
  if (length(south))                # only touch if any south
    shifted[south] <- shifted[south] %m+% months(6)
  
  lubridate::month(shifted)          # return integer 1–12
}

subduction_clock <- function(pmonth) {
  frozen <- pmonth %in% c(12, 1, 2, 3)      # Dec–March (clock at 0)
  data.frame(
    age_min = ifelse(frozen, 0L, 0L), # May=0 … Nov=
    age_max = ifelse(frozen, 0L, 3L)       # May=5 … Nov=11
  )
}
df_carbon_with_poc$PSEUDO_MONTH <- pseudo_month(df_carbon_with_poc$TIME,df_carbon_with_poc$LATITUDE)


df_carbon_with_poc <- df_carbon_with_poc %>% 
  mutate(subduction_clock(PSEUDO_MONTH))                         

df_carbon_with_poc %>% View()

df_carbon_with_poc <- df_carbon_with_poc %>% mutate(
#  integrated_poc_bl_1month = integrated_poc*(exp(-age_min/1))^(-1),
#  integrated_poc_up_1month = integrated_poc*(exp(-age_max/1))^(-1),
  integrated_poc_bl_3months = integrated_poc*(exp(-age_min/3))^(-1),
  integrated_poc_up_3months = integrated_poc*(exp(-age_max/3))^(-1),
#  poc_remaining_50yr_bl_1month = integrated_poc_bl_1month * fseq_50yr,
#  poc_remaining_50yr_up_1month = integrated_poc_up_1month * fseq_50yr,
  poc_remaining_50yr_bl_3months = integrated_poc_bl_3months * fseq_50yr,
  poc_remaining_50yr_up_3months = integrated_poc_up_3months * fseq_50yr,
)

df_carbon_with_poc %>% View()
write_rds(df_carbon_with_poc,file = "data/df_carbon_with_POC.Rds")



# --- 1. Split the Data by Season --------------------------------------------
df_carbon_poc_djf <- df_carbon_with_poc %>% filter(month(TIME) %in% c(12, 1, 2))
df_carbon_poc_mam <- df_carbon_with_poc %>% filter(month(TIME) %in% c(3, 4, 5))
df_carbon_poc_jja <- df_carbon_with_poc %>% filter(month(TIME) %in% c(6, 7, 8))
df_carbon_poc_son <- df_carbon_with_poc %>% filter(month(TIME) %in% c(9, 10, 11))





df_carbon_with_poc <- df_carbon_with_poc %>%
  mutate(
    month_val = month(TIME),
    season = case_when(
      month_val %in% c(12, 1, 2)  ~ "DJF",
      month_val %in% c(3, 4, 5)   ~ "MAM",
      month_val %in% c(6, 7, 8)   ~ "JJA",
      month_val %in% c(9, 10, 11) ~ "SON"
    ),
    # Make it a factor in a logical order
    season = factor(season, levels = c("DJF", "MAM", "JJA", "SON"))
  )

m_full_st_sos600_bl <- gam(
  log(integrated_poc_bl_3months) ~ 
    s(LATITUDE, LONGITUDE, by = season, bs = "sos", k = 600) + 
    s(season, bs = "re"),
  data = df_carbon_with_poc,
  method = "REML"
)

m_full_st_sos600_up <- gam(
  log(integrated_poc_up_3months) ~ 
    s(LATITUDE, LONGITUDE, by = season, bs = "sos", k = 600) + 
    s(season, bs = "re"),
  data = df_carbon_with_poc,
  method = "REML"
)


m_full_st_sos600_poc_remaining_bl <- gam(
  log(poc_remaining_50yr_bl_3months) ~ 
    s(LATITUDE, LONGITUDE, by = season, bs = "sos", k = 600) + 
    s(season, bs = "re"),
  data = df_carbon_with_poc,
  method = "REML"
)


m_full_st_sos600_poc_remaining_up <- gam(
  log(poc_remaining_50yr_up_3months) ~ 
    s(LATITUDE, LONGITUDE, by = season, bs = "sos", k = 600) + 
    s(season, bs = "re"),
  data = df_carbon_with_poc,
  method = "REML"
)


# TP k = 200 R-sq.(adj) =  0.178   Deviance explained = 21.9%
# TP k = 600 
# SOS k = 200
# SOS k = 600


# Create base lat-lon grid at 1° resolution
lat_seq <- seq(-90, 89.5, by = 0.25) + 0.5
lon_seq <- seq(-180, 179.5, by = 0.25) + 0.5

# Expand for all seasons
pred_grid <- expand.grid(
  LATITUDE  = lat_seq,
  LONGITUDE = lon_seq,
  season    = c("DJF","MAM","JJA","SON")
) %>%
  mutate(season = factor(season, levels = c("DJF","MAM","JJA","SON")))

pred_grid$log_poc_bl <- predict(m_full_st_sos600_bl, newdata = pred_grid, type = "response")
pred_grid$poc_bl     <- exp(pred_grid$log_poc_bl)

pred_grid$log_poc_up <- predict(m_full_st_sos600_up, newdata = pred_grid, type = "response")
pred_grid$poc_up     <- exp(pred_grid$log_poc_up)

pred_grid$log_poc_50yrs_bl <- predict(m_full_st_sos600_poc_remaining_bl, newdata = pred_grid, type = "response")
pred_grid$poc_50yrs_bl     <- exp(pred_grid$log_poc_50yrs_bl)

pred_grid$log_poc_50yrs_up <- predict(m_full_st_sos600_poc_remaining_up, newdata = pred_grid, type = "response")
pred_grid$poc_50yrs_up     <- exp(pred_grid$log_poc_50yrs_up)



saveRDS(pred_grid,"data/pred_grid_poc_season.Rds")

library(sf)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")

plot_by_season <- function(season_code) {
  # Filter to one season
  df_plot <- pred_grid %>% filter(season == season_code)
  
  ggplot() +
    geom_tile(data = df_plot, aes(x = LONGITUDE, y = LATITUDE, fill = poc_bl)) +
    geom_contour(data = df_plot, aes(x = LONGITUDE, y = LATITUDE, z = poc_bl),
                 color = "white", alpha = 0.3) +
    geom_sf(data = world, fill = "white", color = "white", inherit.aes = FALSE) +
    coord_sf(xlim=c(-180,180), ylim=c(-90,90), expand=FALSE, crs=st_crs(4326)) +
    labs(title = paste("Season:", season_code),
         x="Longitude", y="Latitude") +
    theme_minimal(base_size=20) + scale_fill_viridis(name="Average integratefd POC by event [mg C/m²]")+
    theme(legend.position='right', panel.grid=element_blank())
}

p_djf <- plot_by_season("DJF")
p_mam <- plot_by_season("MAM")
p_jja <- plot_by_season("JJA")
p_son <- plot_by_season("SON")

# Combine with patchwork or ggarrange
library(patchwork)
combined <- p_djf + p_mam + p_jja + p_son + 
  plot_layout(ncol = 2, nrow = 2) +
  plot_annotation(title = "Average POC content of a carbon subduction event per region and season (factor-smooth model)")

ggsave("figures/TimeSpaceVar/4SEASONS/gam_poc_per_area.png",
       combined, width = 18, height = 10)

print(combined)


# Siegel Correction

add_region <- function(df) {
  df %>% 
    mutate(
      region = case_when(
        LATITUDE >  30 & between(LONGITUDE, -100,   20)             ~ "North Atlantic",
        LATITUDE >  30 & (LONGITUDE < -100 | LONGITUDE > 120)       ~ "North Pacific",
        LATITUDE >= 0 & LATITUDE <=  30                              ~ "Northern Tropics",
        LATITUDE <   0 & LATITUDE >= -30                             ~ "Southern Tropics",
        LATITUDE < -30                                               ~ "Southern Ocean",
        TRUE                                                         ~ NA_character_
      )
    )
}



# Estimating POC export :

# Carb prob gives probability that Argo float captures a subduction event per 0.25 degree gridcell and season
# For instance, if p = 0.20, every 5 T, with T a typical timescale, a carbon subduction event happens.
# To get a flux in (mg C m^2 / day ) f
# How to get f from carb_prob_df ? If T is 1 week, then 0.20 * 1/7
carb_prob_df <- readRDS("data/pred_grid_carb_prob.Rds")
pred_grid <- readRDS("data/pred_grid_poc_season.Rds")

carb_prob_df$LONGITUDE <- carb_prob_df$lon_bin
carb_prob_df$LATITUDE <- carb_prob_df$lat_bin


merged_df <- merge(
  carb_prob_df,
  pred_grid,
  by = c("LONGITUDE", "LATITUDE", "season")
)

# Estimating POC export :
# Read in the carb probability data (if not already loaded)
carb_prob_df <- readRDS("data/pred_grid_carb_prob.Rds")
carb_prob_df$LONGITUDE <- carb_prob_df$lon_bin
carb_prob_df$LATITUDE <- carb_prob_df$lat_bin

# Assuming pred_grid is already loaded and has matching LATITUDE and LONGITUDE columns
merged_df <- merge(
  carb_prob_df,
  pred_grid,
  by = c("LONGITUDE", "LATITUDE", "season")
)

# Earth radius in meters (for cell area calculation)
R <- 6371000

# Define season lengths (in days)
season_days <- data.frame(
  season = c("DJF", "MAM", "JJA", "SON"),
  days = c(90, 92, 92, 91)
)

# -----------------------------------------
# 2. Sensitivity Analysis for Different T
# -----------------------------------------
# Define w values (in m/day) for sensitivity analysis
merged_df <- add_region(merged_df)

w_values <- c(20,50,100,200,300,400,500)
# max depth of mixed layer (subd events are detected below this depth)
H <- 200
# w = L / T => T = L/W
# 1/T = W/L
# Create a results data frame to store total export for each T
region_levels <- c("North Atlantic","Southern Ocean",
                   "North Pacific","Northern Tropics","Southern Tropics")

sensitivity_region <- expand.grid(
  w_days = w_values,
  region = region_levels
) %>% 
  mutate(
    export_PgC_bl        = NA_real_,
    export_PgC_up        = NA_real_,
    export50_PgC_bl      = NA_real_,
    export50_PgC_up      = NA_real_
  )

annual_export_region_list <- list()
# Loop over each T value
for (i in seq_along(w_values)) {
  w_val <- w_values[i]
  
  # For each w compute the daily POC flux: poc_flux = poc (mg C /m^-2) * proportion * (1/T) (s^(-1))
  sensitivity_df <- merged_df %>%
    mutate(poc_flux_bl = poc_bl * proportion * w_val/H, # mg C m^2 day^1
           poc_flux_up = poc_up * proportion * w_val/H,
           poc_flux_50years_bl = poc_50yrs_bl * proportion * w_val/H,
           poc_flux_50years_up = poc_50yrs_up * proportion * w_val/H) %>%
    # Merge the season lengths (in days)
    left_join(season_days, by = "season") %>%
    # Calculate the cell area for a 0.25° x 0.25° grid cell (m²)
    mutate(
      delta_rad = 0.25 * pi / 180,
      phi = LATITUDE * pi / 180,
      cell_area = R^2 * delta_rad^2 * cos(phi),      
      export_season_mg_m2_bl = poc_flux_bl * days, #Multiply the daily flux by the number of days in the season to get mg C/m² per season
      export_season_mg_m2_up = poc_flux_up * days,
      export_season_mg_m2_50years_bl = poc_flux_50years_bl * days,
      export_season_mg_m2_50years_up = poc_flux_50years_up * days,
      # Multiply by the grid cell area to get mg C per grid cell per season
      export_season_mg_bl = export_season_mg_m2_bl * cell_area,
      export_season_mg_up = export_season_mg_m2_up * cell_area,
      export_season_mg_50years_bl = export_season_mg_m2_50years_bl * cell_area,
      export_season_mg_50years_up = export_season_mg_m2_50years_up * cell_area
    )
  

  # Sum over seasons to obtain annual export per grid cell
  annual_export_df <- sensitivity_df %>%
    group_by(LONGITUDE, LATITUDE) %>%
    summarise(export_annual_mg_bl = sum(export_season_mg_bl, na.rm = TRUE),
              export_annual_mg_up = sum(export_season_mg_up, na.rm = TRUE),
              export_annual_mg_50years_bl = sum(export_season_mg_50years_bl, na.rm = TRUE),
              export_annual_mg_50years_up = sum(export_season_mg_50years_up, na.rm = TRUE)) %>%
    ungroup()
  # total export per seasons
  annual_export_df <- add_region(annual_export_df)  
  annual_export_region_list[[i]] <- annual_export_df %>% group_by(region) %>% summarize(
    total_export_mg_bl = sum(export_annual_mg_bl,na.rm = T)*1e-18,
    total_export_mg_up = sum(export_annual_mg_up,na.rm = T)*1e-18
  )
  
  
  # Total global annual export in mg C
  total_export_mg_bl <- sum(annual_export_df$export_annual_mg_bl, na.rm = TRUE)
  total_export_mg_up <- sum(annual_export_df$export_annual_mg_up, na.rm = TRUE)
  total_export_mg_50years_bl <- sum(annual_export_df$export_annual_mg_50years_bl, na.rm = TRUE)
  total_export_mg_50years_up <- sum(annual_export_df$export_annual_mg_50years_up, na.rm = TRUE)
  
  # Convert mg C to Pg C (1 mg = 1e-18 Pg)
  total_export_PgC_bl <- total_export_mg_bl * 1e-18
  total_export_PgC_up <- total_export_mg_up * 1e-18

  total_export_PgC_50yrs_bl <- total_export_mg_50years_bl * 1e-18
  total_export_PgC_50yrs_up <- total_export_mg_50years_up * 1e-18
  
  # Store the result for this T value
  sensitivity_results$total_export_PgC_bl[i] <- total_export_PgC_bl
  sensitivity_results$total_export_PgC_up[i] <- total_export_PgC_up
  
  sensitivity_results$total_export_PgC_50yrs_bl[i] <- total_export_PgC_50yrs_bl
  sensitivity_results$total_export_PgC_50yrs_up[i] <- total_export_PgC_50yrs_up
  
}

names(annual_export_region_list) <- w_values <- c(20,50,100,200,300,400,500)

export_region <- annual_export_region_list %>% bind_rows(.id = "w")  %>% group_by(region) %>%
  summarize(min_export = min(total_export_mg_bl,total_export_mg_up),
            max_export = max(total_export_mg_bl,total_export_mg_up)
            )  
saveRDS(object = export_region,file = "data/sensitivity_region.Rds")
saveRDS(object = sensitivity_results,file = "data/sensitivity_res.Rds")
# Print the sensitivity results
print(sensitivity_results)

library(ggplot2)
library(viridis)
sensitivity_results <- readRDS(file = "data/sensitivity_res.Rds")

df <- sensitivity_results
df$total_export_PgC_mean <- (df$total_export_PgC_bl + df$total_export_PgC_up)/2
df$total_export_PgC_50yrs_mean <- (df$total_export_PgC_50yrs_bl + df$total_export_PgC_50yrs_up)/2
# 1) Prepare the tidy data
df_plot <- bind_rows(
  df %>% 
    transmute(
      w_days,
      metric = "Total export",
      min    = total_export_PgC_bl,
      mean   = total_export_PgC_mean,
      max    = total_export_PgC_up
    ),
  df %>% 
    transmute(
      w_days,
      metric = "Export >50 yr",
      min    = total_export_PgC_50yrs_bl,
      mean   = total_export_PgC_50yrs_mean,
      max    = total_export_PgC_50yrs_up
    )
)

# 2) Single‐panel plot
plot_esp <- ggplot(df_plot, aes(
  x = w_days, 
  y = mean, 
  ymin = min, 
  ymax = max, 
  color = metric
)) +
  geom_pointrange(
    size = 2,
    fatten = 2
  ) +
  scale_color_viridis_d(end = 0.8, option = "D") +
  scale_x_continuous(breaks = df$w_days) +
  labs(
    x     = "Characteristic vertical velocity W (m/day)",
    y     = expression("Export rate (PgC/ year)"),
    color = NULL,
    title = "Eddy Subduction Pump Export Estimates for an e-folding timescale of 3 months"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(plot_esp,filename = "figures/sensitivity_ESP_magnitude.png",width=10,height=8)


saveRDS(object = merged_df,"data/merged_dataset_poc_estim.Rds")

merged_df <- readRDS("data/merged_dataset_poc_estim.Rds")

# POC flux assuming T = 4 days (mean speed of 50 m/s)
# POC flux assuming w =  200 m / day (T = L/W)
merged_df <- merged_df %>% mutate(poc_flux50w_bl = poc_bl * proportion * 1/4,
                                  poc_flux50w_50yrs_bl = poc_50yrs_bl * proportion *1/4,
                                  poc_flux50w_up = poc_up * proportion * 1/4,
                                  poc_flux50w_50yrs_up = poc_50yrs_up * proportion * 1/4)


merged_df <- merged_df %>% mutate(poc_flux200w_bl = poc_bl * proportion * 200/200)


annual_mean <- merged_df %>%
  group_by(LONGITUDE, LATITUDE) %>%
  summarize(poc_flux_annual = mean(poc_flux200w_bl, na.rm = TRUE))

stipple_resolution <- 5

argo_bins <- df_argo_clean %>%
  mutate(
    lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
    lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
  ) %>%
  group_by(lon_bin_stipple, lat_bin_stipple) %>%
  summarize(count = n(), .groups = "drop")

full_grid <- expand.grid(
  lon_bin_stipple = seq(-180, 175, by = stipple_resolution),
  lat_bin_stipple = seq(-90, 90, by = stipple_resolution)
) %>%
  as_tibble()

# Fill in any missing bins with count = 0
argo_bins_full <- full_grid %>%
  left_join(argo_bins, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
  mutate(count = ifelse(is.na(count), 0, count))

# 2. Identify undersampled areas (e.g., fewer than 1 profile)
undersampled <- argo_bins_full %>% filter(count < 1)

# 3. Generate corner points for each undersampled grid cell (4 corners per cell)
undersampled_corners <- undersampled %>%
  rowwise() %>%
  mutate(corners = list(
    data.frame(
      LON = c(lon_bin_stipple,
              lon_bin_stipple + stipple_resolution,
              lon_bin_stipple,
              lon_bin_stipple + stipple_resolution),
      LAT = c(lat_bin_stipple,
              lat_bin_stipple,
              lat_bin_stipple + stipple_resolution,
              lat_bin_stipple + stipple_resolution)
    )
  )) %>%
  ungroup() %>%
  unnest(corners)

# 4. Filter the prediction grid to the desired season
df_plot <- annual_mean
world <- ne_countries(scale = "medium", returnclass = "sf")

# 5. Build the plot using the provided discrete scale (common_scale)
ggplot() +
  geom_tile(data = df_plot, aes(x = LONGITUDE, y = LATITUDE, fill = poc_flux_annual)) +
  geom_contour(data = df_plot, aes(x = LONGITUDE, y = LATITUDE, z = poc_flux_annual),
               color = "white", alpha = 0.3,binwidth = 1) +
  geom_sf(data = world, fill = "white", color = "white", inherit.aes = FALSE) +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE, crs = st_crs(4326)) +   # Apply the discrete scale passed as argument
  geom_point(data = undersampled_corners,
             aes(x = LON, y = LAT),
             color = "white", alpha = 1, size = 1, shape = 20) +
  labs(title = paste("Average POC flux",sep = "",x = "Longitude", y = "Latitude")) +
  scale_fill_viridis_b()+
  theme_minimal(base_size = 25) +  # Increased text size for readability
  theme(legend.position = 'top', 
        legend.title = element_text(size = 25, face = "bold", margin = margin(r = 100)),  # Improve legend title
        panel.grid = element_blank(),
        legend.text = element_text(size = 20),  # Improve legend labels
        legend.direction = "horizontal", 
        legend.box = "horizontal",
        legend.key.height = unit(1,"lines"),
        legend.key.width = unit(8,"lines"))+ 
  guides(color = guide_legend(title.position = "top", 
                              # hjust = 0.5 centres the title horizontally
                              title.hjust = 0.5,
                              label.position = "bottom")) +theme_map()




make_discrete_scale <- function(prob_min, prob_max, binwidth = 0.1) {
  # E.g. if prob_max is 0.6 and binwidth is 0.1 => breaks = seq(0, 0.6, 0.1)
  breaks_vec <- seq(prob_min, prob_max, by = binwidth)
  
  scale_fill_viridis_b(
    name = "POC flux (mg C m⁻² day⁻¹)",
    breaks = breaks_vec,
    limits = c(prob_min, prob_max),
    oob = scales::squish
  )
}

poc_scale <- make_discrete_scale(0, 4, binwidth = 0.5)

pred_grid$poc_50yrs


stipple_resolution <- 5
plot_by_season_pocflux <- function(pred_grid, world_data, season_label, common_scale, argo_months) {
  # 1. Aggregate Argo float locations into 5° bins for stippling
  argo_bins <- df_argo_clean %>%
    filter(month(TIME) %in% argo_months) %>%
    mutate(
      lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
      lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
    ) %>%
    group_by(lon_bin_stipple, lat_bin_stipple) %>%
    summarize(count = n(), .groups = "drop")
  
  full_grid <- expand.grid(
    lon_bin_stipple = seq(-180, 175, by = stipple_resolution),
    lat_bin_stipple = seq(-90, 90, by = stipple_resolution)
  ) %>%
    as_tibble()
  
  # Fill in any missing bins with count = 0
  argo_bins_full <- full_grid %>%
    left_join(argo_bins, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
    mutate(count = ifelse(is.na(count), 0, count))
  
  # 2. Identify undersampled areas (e.g., fewer than 1 profile)
  undersampled <- argo_bins_full %>% filter(count < 1)
  
  # 3. Generate corner points for each undersampled grid cell (4 corners per cell)
  undersampled_corners <- undersampled %>%
    rowwise() %>%
    mutate(corners = list(
      data.frame(
        LON = c(lon_bin_stipple,
                lon_bin_stipple + stipple_resolution,
                lon_bin_stipple,
                lon_bin_stipple + stipple_resolution),
        LAT = c(lat_bin_stipple,
                lat_bin_stipple,
                lat_bin_stipple + stipple_resolution,
                lat_bin_stipple + stipple_resolution)
      )
    )) %>%
    ungroup() %>%
    unnest(corners)
  
  # 4. Filter the prediction grid to the desired season
  df_plot <- pred_grid %>% filter(season == season_label)
  
  # 5. Build the plot using the provided discrete scale (common_scale)
  yearly_plot <-  ggplot() +
    geom_tile(data = df_plot, aes(x = LONGITUDE, y = LATITUDE, fill = poc_flux_annual)) +
    geom_contour(data = df_plot, aes(x = LONGITUDE, y = LATITUDE, z = poc_flux_annual),
                 color = "white", alpha = 0.3,binwidth = 0.5) +
    geom_sf(data = world, fill = "white", color = "white", inherit.aes = FALSE) +
    coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE, crs = st_crs(4326)) +
    poc_scale +   # Apply the discrete scale passed as argument
    geom_point(data = undersampled_corners,
               aes(x = LON, y = LAT),
               color = "white", alpha = 1, size = 1, shape = 20) +
    labs(title = paste("Average Yearly POC flux"),
         x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 25) +  # Increased text size for readability
    theme(legend.position = 'bottom', 
          legend.title = element_text(size = 25, face = "bold", margin = margin(r = 100)),  # Improve legend title
          panel.grid = element_blank(),
          legend.text = element_text(size = 20),  # Improve legend labels
          legend.direction = "horizontal", 
          legend.box = "horizontal",
          legend.key.height = unit(1,"lines"),
          legend.key.width = unit(8,"lines"))+ 
    guides(color = guide_legend(title.position = "top", 
                                # hjust = 0.5 centres the title horizontally
                                title.hjust = 0.5,
                                label.position = "bottom")) }
  
  
  ggsave(plot = yearly_plot,"figures/TimeSpaceVar/4SEASONS/poc_map_annual.png",width=18,height = 10)
  
  p_djf <- plot_by_season_pocflux(
    pred_grid = merged_df, 
    world_data = world, 
    season_label = "DJF", 
    common_scale = poc_scale,  # for example, created via make_discrete_scale(...)
    argo_months = c(12, 1, 2)
  )
  
  p_mam <- plot_by_season_pocflux(
    pred_grid = merged_df, 
    world_data = world, 
    season_label = "MAM", 
    common_scale = poc_scale, 
    argo_months = c(3, 4, 5)
  )
  
  p_jja <- plot_by_season_pocflux(
    pred_grid = merged_df, 
    world_data = world, 
    season_label = "JJA", 
    common_scale = poc_scale,  
    argo_months = c(6, 7, 8)
  )
  
  p_son <- plot_by_season_pocflux(
    pred_grid = merged_df, 
    world_data = world, 
    season_label = "SON", 
    common_scale = poc_scale,  
    argo_months = c(9, 10, 11)
  )
  
  
  combined_poc <- ggarrange(p_djf,p_mam,p_jja,p_son,ncol = 2,nrow=2,
                            common.legend = T,legend="bottom")
  ggsave("figures/TimeSpaceVar/4SEASONS/poc_map.png",
         combined_poc, width = 16, height = 12)
  
  make_discrete_scale <- function(prob_min, prob_max, binwidth = 0.1) {
    # E.g. if prob_max is 0.6 and binwidth is 0.1 => breaks = seq(0, 0.6, 0.1)
    breaks_vec <- seq(prob_min, prob_max, by = binwidth)
    
    scale_fill_viridis_b(
      name = "POC flux sequestered for > 50 years (mg C m⁻² day⁻¹)",
      breaks = breaks_vec,
      limits = c(prob_min, prob_max),
      oob = scales::squish
    )
  }
  poc_scale <- make_discrete_scale(0, 15, binwidth = 1)
  
  plot_by_season_pocflux_50years <- function(pred_grid, world_data, season_label, common_scale, argo_months) {
    # 1. Aggregate Argo float locations into 5° bins for stippling
    argo_bins <- df_argo_clean %>%
      filter(month(TIME) %in% argo_months) %>%
      mutate(
        lon_bin_stipple = floor(LONGITUDE / stipple_resolution) * stipple_resolution,
        lat_bin_stipple = floor(LATITUDE / stipple_resolution) * stipple_resolution
      ) %>%
      group_by(lon_bin_stipple, lat_bin_stipple) %>%
      summarize(count = n(), .groups = "drop")
    
    full_grid <- expand.grid(
      lon_bin_stipple = seq(-180, 175, by = stipple_resolution),
      lat_bin_stipple = seq(-90, 90, by = stipple_resolution)
    ) %>%
      as_tibble()
    
    # Fill in any missing bins with count = 0
    argo_bins_full <- full_grid %>%
      left_join(argo_bins, by = c("lon_bin_stipple", "lat_bin_stipple")) %>%
      mutate(count = ifelse(is.na(count), 0, count))
    
    # 2. Identify undersampled areas (e.g., fewer than 1 profile)
    undersampled <- argo_bins_full %>% filter(count < 1)
    
    # 3. Generate corner points for each undersampled grid cell (4 corners per cell)
    undersampled_corners <- undersampled %>%
      rowwise() %>%
      mutate(corners = list(
        data.frame(
          LON = c(lon_bin_stipple,
                  lon_bin_stipple + stipple_resolution,
                  lon_bin_stipple,
                  lon_bin_stipple + stipple_resolution),
          LAT = c(lat_bin_stipple,
                  lat_bin_stipple,
                  lat_bin_stipple + stipple_resolution,
                  lat_bin_stipple + stipple_resolution)
        )
      )) %>%
      ungroup() %>%
      unnest(corners)
    
    # 4. Filter the prediction grid to the desired season
    df_plot <- pred_grid %>% filter(season == season_label)
    
    # 5. Build the plot using the provided discrete scale (common_scale)
    ggplot() +
      geom_tile(data = df_plot, aes(x = LONGITUDE, y = LATITUDE, fill = poc_flux3days50yrs)) +
      geom_contour(data = df_plot, aes(x = LONGITUDE, y = LATITUDE, z = poc_flux3days50yrs),
                   color = "white", alpha = 0.3,binwidth = 1) +
      geom_sf(data = world_data, fill = "white", color = "white", inherit.aes = FALSE) +
      coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE, crs = st_crs(4326)) +
      common_scale +   # Apply the discrete scale passed as argument
      geom_point(data = undersampled_corners,
                 aes(x = LON, y = LAT),
                 color = "white", alpha = 1, size = 1, shape = 20) +
      labs(title = paste("Average POC flux (",season_label,")",sep = ""),
           x = "Longitude", y = "Latitude") +
      theme_minimal(base_size = 25) +  # Increased text size for readability
      theme(legend.position = 'top', 
            legend.title = element_text(size = 20, face = "bold", margin = margin(r = 150)),  # Improve legend title
            panel.grid = element_blank(),
            legend.text = element_text(size = 18),  # Improve legend labels
            legend.direction = "horizontal", 
            legend.box = "horizontal",
            legend.key.height = unit(1,"lines"),
            legend.key.width = unit(5,"lines"))+ 
      guides(color = guide_legend(title.position = "top", 
                                  # hjust = 0.5 centres the title horizontally
                                  title.hjust = 0.5,
                                  label.position = "bottom")) 
  }
  
  p_djf <- plot_by_season_pocflux_50years(
    pred_grid = merged_df, 
    world_data = world, 
    season_label = "DJF", 
    common_scale = poc_scale,  # for example, created via make_discrete_scale(...)
    argo_months = c(12, 1, 2)
  )
  
  p_mam <- plot_by_season_pocflux_50years(
    pred_grid = merged_df, 
    world_data = world, 
    season_label = "MAM", 
    common_scale = poc_scale,  # for example, created via make_discrete_scale(...)
    argo_months = c(3, 4, 5)
  )
  
  p_jja <- plot_by_season_pocflux_50years(
    pred_grid = merged_df, 
    world_data = world, 
    season_label = "JJA", 
    common_scale = poc_scale,  # for example, created via make_discrete_scale(...)
    argo_months = c(6, 7, 8)
  )
  
  p_son <- plot_by_season_pocflux_50years(
    pred_grid = merged_df, 
    world_data = world, 
    season_label = "SON", 
    common_scale = poc_scale,  # for example, created via make_discrete_scale(...)
    argo_months = c(9, 10, 11)
  )
  
  
  combined_poc_50yrs <- ggarrange(p_djf,p_mam,p_jja,p_son,ncol = 2,nrow=2,
                                  common.legend = T,legend="bottom")
  ggsave("figures/TimeSpaceVar/4SEASONS/poc_map_50yrs.png",
         combined_poc_50yrs, width = 16, height = 12)
  
  
