
library(tidync)   # install.packages("tidync")
src <- tidync("data/ArgoTS_eddydiffusivity_20052015_1deg.nc")
print(src)        # lists every grid, variable and dimension
hyper_vars(src)   # compact table of variables in the active grid


library(readr)
kappa_surface <- read_csv("data/eddy_diffusivity_top200m.csv")
kappa_surface <- read_csv("data/eddy_diffusivity_top500m.csv")
kappa_surface_interior <- read_csv("data/eddy_diffusivity_interior.csv")

kappa_surface$lat <- kappa_surface$lat -63.5
kappa_surface_interior$lat <- kappa_surface_interior$lat -63.5

kappa_surface <- kappa_surface_interior

df_carbon_with_poc <- read_csv("data/df_carbon_subduction_anom_with_poc_fromgali.csv") %>% filter(integrated_poc > 0.01) 
df_carbon <- df_carbon_with_poc 

df_carbon$LATITUDE %>% summary()
df_carbon$LONGITUDE %>% summary()


ggplot(kappa_surface, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = kh200)) +                        # colourful pixels
  geom_contour(aes(z = kh200), colour = "white", 
               alpha = 0.15, bins = 8) +               # faint isolines
  geom_sf(data = world, fill = "grey80", colour = "black",
          linewidth = 0.2, inherit.aes = FALSE) +      # continents
  coord_sf(expand = TRUE) +
  scale_fill_viridis_b(trans = "log10",                # κ spans orders‑of‑mag
                       option = "turbo",
                       name = expression(K[h]~(m^2~s^{-1}))) +
  labs(title = "Argo‑derived horizontal diffusivity (3° bins)",
       x = "Longitude", y = "Latitude") +
  theme_minimal() 

ggplot(kappa_surface_interior, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = kh200)) +                        # colourful pixels
  geom_contour(aes(z = kh200), colour = "white", 
               alpha = 0.15, bins = 8) +               # faint isolines
  geom_sf(data = world, fill = "grey80", colour = "black",
          linewidth = 0.2, inherit.aes = FALSE) +      # continents
  coord_sf(expand = TRUE) +
  scale_fill_viridis(trans = "log10",                # κ spans orders‑of‑mag
                       option = "turbo",
                       name = expression(K[h]~(m^2~s^{-1}))) +
  labs(title = "Argo‑derived horizontal diffusivity (3° bins)",
       x = "Longitude", y = "Latitude") +
  theme_minimal() 


# matrices must have same column order → (lat, lon)
coords_carbon <- as.matrix(df_carbon      %>% select(LATITUDE,  LONGITUDE))
coords_kappa  <- as.matrix(kappa_surface  %>% select(lat,      lon))

idx <- nn2(coords_kappa, coords_carbon, k = 1)$nn.idx 



df_carbon$kappa_hsurf <- kappa_surface$kh200[idx]     


L_half <- 25000      # metres  (1 km half‑width filament)


df_carbon <- df_carbon %>% 
  mutate(tau_days = ((L_half^2) / (pi^2 * kappa_hsurf)) / 86400)   # convert s → d


# ── 1. prep a 1°×1° grid (optional but helps geom_tile) ───────────────────────
df_grid <- df_carbon %>% 
  mutate(lon_bin = floor(LONGITUDE),             # 1‑degree bins
         lat_bin = floor(LATITUDE)) %>% 
  group_by(lon_bin, lat_bin) %>% 
  summarise(kapp_surf = mean(kappa_hsurf,na.rm=TRUE),
            tau_days  = mean(tau_days,na.rm=TRUE))

# ── 2. world polygons for context ─────────────────────────────────────────────
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# ── 3. quick map ──────────────────────────────────────────────────────────────
ggplot(df_grid, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = kapp_surf)) +                        # colourful pixels
  geom_contour(aes(z = kapp_surf), colour = "white", 
               alpha = 0.15, bins = 8) +               # faint isolines
  geom_sf(data = world, fill = "grey80", colour = "black",
          linewidth = 0.2, inherit.aes = FALSE) +      # continents
  coord_sf(expand = FALSE) +
  scale_fill_viridis(trans = "log10",                # κ spans orders‑of‑mag
                       option = "turbo",
                       name = expression(K[h]~(m^2~s^{-1}))) +
  labs(title = "Argo‑derived horizontal diffusivity (3° bins)",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

df_grid$tau_days_esp <- 10*df_grid$tau_days

ggplot(df_grid, aes(x = lon_bin, y = lat_bin)) +
  geom_tile(aes(fill = tau_days_esp)) +                        # colourful pixels
  geom_contour(aes(z = tau_days_esp), colour = "white", 
               alpha = 0.15, bins = 8) +               # faint isolines
  geom_sf(data = world, fill = "grey80", colour = "black",
          linewidth = 0.2, inherit.aes = FALSE) +      # continents
  coord_sf(expand = FALSE) +
  scale_fill_viridis(trans="log10",option = "turbo",
                       name = "e-folding timescale (days)") +
  labs(title = "Argo‑derived e-folding timescale (3° bins) average eddy diffusivity over 500m",
       x = "Longitude", y = "Latitude") +
  theme_minimal()


library(dplyr)      # tidy tools


df_carbon$tau_days_esp <- 10*df_carbon$tau_days
 
df_carbon$original_poc_min <- df_carbon$integrated_poc *exp(1/df_carbon$tau_days_esp)    
