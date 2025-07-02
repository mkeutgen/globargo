# Partie2BIS detection downscaled
setwd("~/Documents/GLOBARGO")

library(tidyverse)
library(robustbase)
library(gsw)
library(zoo)
library(oce)
library(ggpubr)
set.seed(10)
wmolist <- readRDS("~/Documents/ARGO/BGC_Argo_WMO_PSAL_BBP_DOXY_TEMP.rds")
# Read the manually classified dataframe
classified_df <- read.csv("/data/GLOBARGO/data/classification_results_v3.csv")
classified_df$WMO <- gsub("_plot", "", classified_df$WMO)

#Create two distinct df with high and low confidence subduction events
cat1_events <- classified_df %>% filter(Category==1)
cat2_events <- classified_df %>% filter(Category==2) 

wmo_cat1 <- cat1_events$WMO %>% unique()
wmo_cat2 <- cat2_events$WMO %>% unique()

detected.events.list <- list()
# j, indicator for the WMO
for (j in seq_along(wmo_cat1) ) {  #seq_along(wmo_cat1)
  try({
    # Start with first WMO of the random sample
    wmo <- wmo_cat1[j]
    #wmo <- wmo_cat2[j]
    #wmo <- 5904470
    #wmo <-  5905073  # with soooo much zero residuals, for instance prof 65 has sooo tiny little bumps and it'd
    # be nice if there was a way for the method not to detect it... However 29 and 40 are legit AF
    #wmo <- 5905073
    # Cleaned Llort WMO
    # wmo <- 5904677 # prof 27
    #wmo <- 7901099
    #wmo <-5906581
    #wmo <- 5906635
    #wmo <- 6903024 # Problem of very small residuals, blowing up residuals
    #wmo <- 5904179
    #wmo <- 1902332
    cycle_number <- cat1_events %>% filter(WMO== wmo) %>% select(Cycle) %>% unique() %>% as_vector()
    #cycle_number <- cat2_events %>% filter(WMO== wmo) %>% select(Cycle) %>%  unique() %>% as_vector()
    
    
    # LOT OF ANOMALIES IFF
    # Download the whole float data for this WMO
    
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
    data_df <- data_df %>% filter(data_df$CYCLE_NUMBER %in% cycle_number)
    
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
    
    
    
    list_of_tibbles <- A  %>%
      group_by(CYCLE_NUMBER) %>%
      group_split()
    
    
    
    downscale_data_fun <- function(df,b=20,cutoff=2) {
      # Select and pivot data
      data <- df %>%
        select(PRES_ADJUSTED, SCALE.RES.ROB, VAR,CYCLE_NUMBER,LONGITUDE,LATITUDE,TIME) %>%
        pivot_wider(names_from = VAR, values_from = SCALE.RES.ROB)
      
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
      
      # Additional mutations
      downscaled_data <- downscaled_data %>%
        mutate(OUT.S =  ifelse(abs(AOU) > cutoff & abs(SPIC) > cutoff & AOU < 0, 1, 0),
               OUT.T = ifelse(BBP700_ADJUSTED> 0 & abs(BBP700_ADJUSTED) > 1 & abs(AOU) > cutoff & abs(SPIC) > cutoff & AOU < 0, 1, 0))
      
      # Return downscaled data
      return(downscaled_data)
    }
    
    output <- lapply(list_of_tibbles,downscale_data_fun,b=40)
    output <- output %>% bind_rows()
    # Matrix of downscaled residuals
    B <- output %>% bind_rows()
    
    prop_zero <- B %>% select(BBP700_ADJUSTED,AOU,SPIC,CYCLE_NUMBER,PRES_ADJUSTED,TIME) %>%
      unique() %>%
      pivot_longer(cols = c("AOU","SPIC","BBP700_ADJUSTED"),names_to = "VAR",values_to = "VALUE") %>%
      group_by(CYCLE_NUMBER,VAR) %>%  summarize(zero_proportion = mean(VALUE == 0, na.rm = TRUE),
                                                .groups = 'drop')
    
    pivot_proportion <- prop_zero %>%
      pivot_wider(names_from = VAR, values_from = zero_proportion)
    
    # Now filter for cycles where both AOU and SPIC have <= 50% zero proportions
    selected_cycles <- pivot_proportion %>%
      filter(AOU <= 0.5, SPIC <= 0.5) %>%
      pull(CYCLE_NUMBER) # Extracting the cycle numbers as a vector
    
    
    # It seems 55 % of prop zero would be acceptable cutoff, if more than half of the residuals are zero
    # then profile is seen as too monotonic
    
    
    # Find time where outlying and associated pressure level :
     # !!!! in this case OUT.T == 1 so only CARBON outlying profile are selected
    carb_eddy.id <- B %>% filter(OUT.T==1) %>% unique()
    carb_eddy.id <- carb_eddy.id %>% filter(PRES_ADJUSTED <= 700) %>% filter(PRES_ADJUSTED >= 200)
    
    carb_eddy.id$WMO <- wmo
    detected.events.list[[j]] <- carb_eddy.id
    
    
    # Plotting :
    prof.plot <- list()
    res.plot  <- list()
    list.plots <- list()
    
    
    for (i in seq_along(carb_eddy.id$CYCLE_NUMBER)) {
      
      current_cycle <- carb_eddy.id$CYCLE_NUMBER[i]
      current_data <- A %>% filter(CYCLE_NUMBER == current_cycle)
      current_eddy <- carb_eddy.id %>% filter(CYCLE_NUMBER == current_cycle)
      
      prof.plot[[i]] <- current_data %>%
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
        theme(legend.position = "bottom") +
        geom_vline(
          xintercept = current_eddy$PRES_ADJUSTED[current_eddy$OUT.S == 1],
          color = "darkgreen",
          alpha = 0.3,
          size = 1
        ) +
        geom_vline(
          xintercept = current_eddy$PRES_ADJUSTED[current_eddy$OUT.T == 1],
          color = "black",
          alpha = 0.3,
          size = 1
        )
      
      df <- current_data
      df.ds <- downscale_data_fun(df, b = 40)
      df.ds <- df.ds %>% ungroup() %>% select(AOU, SPIC, BBP700_ADJUSTED, PRES_ADJUSTED) %>%
        pivot_longer(cols = !PRES_ADJUSTED, names_to = "VAR", values_to = "VALUE")
      
      hline_data <- data.frame(
        VAR = c("AOU", "AOU", "SPIC", "SPIC", "BBP700_ADJUSTED", "BBP700_ADJUSTED"),
        hline = c(-2, 2, -2, 2, -1, 1),
        label = c("-2 sigma", "+2 sigma", "-2 sigma", "+2 sigma", "-1 sigma", "+1 sigma")
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
          "-2 sigma" = "red",
          "+2 sigma" = "red",
          "-1 sigma" = "blue",
          "+1 sigma" = "blue"
        )) + theme(legend.position = "bottom") +
        geom_vline(
          xintercept = current_eddy$PRES_ADJUSTED[current_eddy$OUT.S == 1],
          color = "darkgreen",
          alpha = .3,
          size = 1
        )
      
      if (any(current_eddy$OUT.T == 1)) {
        res.plot[[i]] <- res.plot[[i]] +
          geom_vline(
            xintercept = current_eddy$PRES_ADJUSTED[current_eddy$OUT.T == 1],
            color = "black",
            alpha = 1,
            size = 1
          )
      }
      
      combined_plot <- ggarrange(prof.plot[[i]], res.plot[[i]], common.legend = F, legend = "bottom", nrow = 2)
      
      annotation_text <- paste("Cycle Number: ", current_cycle,
                               "\nFloat ID (WMO): ", wmo,
                               "\nLongitude: ", current_eddy$LONGITUDE,
                               "\nLatitude: ", current_eddy$LATITUDE,
                               "\nTime: ", format(as.POSIXct(current_eddy$TIME, origin = "1970-01-01"), "%Y-%m-%d"))
      
      combined_plot <- annotate_figure(combined_plot, 
                                       top = text_grob(annotation_text, face = "bold", size = 10))
      
      list.plots[[i]] <- combined_plot
    }
    
    
    
    
    # Check if list.plots is not empty
    if (length(list.plots) > 0) {
      # Create directory if it doesn't exist
      dir <- paste0("/data/GLOBARGO/figures/CarbonFiguresCat1/", wmo)
      if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
      }
      
      # Save each plot in the list to the directory
      for (k in seq_along(list.plots)) {
        cycle_number <- carb_eddy.id$CYCLE_NUMBER[k]
        file_name <- paste0(dir, "/", wmo, "_plot_cycle_",cycle_number, ".png")  # or .pdf, .jpeg as per your requirement
        ggsave(file_name, list.plots[[k]], width = 10, height = 10)  # Adjust size as needed
      }
    }
    
    
    
    
    
  },silent=TRUE)
}

# Convert the 'WMO' column to a consistent data type (e.g., character) for each element of the list
detected.events.list <- lapply(detected.events.list, function(x) {
  x$WMO <- as.character(x$WMO)
  return(x)
})

# Combine the list into a single data frame using bind_rows() 948
detected.events.df <- detected.events.list %>% bind_rows()
detected.events.df 

write_csv(detected.events.df,"/data/GLOBARGO/data/detected_events_unique_with_carbon_cat1.csv")

#write_csv(detected.events.df,"/data/GLOBARGO/data/detected_events_unique_with_carbon_cat2.csv")

#cat1_df <- read_csv("~/Documents/ARGO/data/detected_events_unique_with_carbon_cat1.csv")
cat1_df$OUT.T %>% sum()
detected.events.df$OUT.T %>% sum()
detected.events.df %>% nrow()

# Carbon Detection
carbon_cat1 <- read_csv("/data/GLOBARGO/data/detected_events_unique_with_carbon_cat1.csv")
carbon_cat2 <- read_csv("/data/GLOBARGO/data/detected_events_unique_with_carbon_cat2.csv")

