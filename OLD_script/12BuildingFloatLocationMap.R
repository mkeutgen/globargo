float_ids <- 6901767
float_ids <- 5906312

subduction_events

det_event_cat1 <- read_csv("/data/GLOBARGO/data/detected_events_unique_with_carbon_cat1.csv")

subduction_events <- det_event_cat1 %>% filter(WMO == 6901767)
subduction_events <- det_event_cat1 %>% filter(WMO == 5906312)

highlight_value=179
highlight_value=32

highlight_type='cycle_number'
highlight_trajectory <- function(float_ids, 
                                 color='multiple',
                                 float_profs=NULL,
                                 position=NULL,
                                 title="Float Trajectories",
                                 highlight_type='cycle_number',  # Either 'time' or 'cycle_number'
                                 highlight_value=NULL,  # Time or Cycle number to highlight
                                 return_ggplot=TRUE) {
  
  # Check if settings are initialized
  if (!exists("Setting")) {
    initialize_argo()
  }
  
  subduction_events <- det_event_cat1 %>% filter(WMO == float_ids)
  
  
  # Download Sprof files if needed
  good_float_ids = download_multi_floats(float_ids)
  
  if (length(good_float_ids) == 0) {
    warning('No valid floats found.')
    return(1)
  }
  
  # Load float data
  loaded = load_float_data(float_ids=good_float_ids)
  Data = loaded$Data
  
  # If position argument is set to 'first' or 'last', adjust longitude/latitude accordingly
  if (!is.null(position)) {
    nfloats = length(Data)
    if (position == "first") {
      for (f in 1:nfloats) {
        Data[[f]]$LONGITUDE <- Data[[f]]$LONGITUDE[,1]
        Data[[f]]$LATITUDE <- Data[[f]]$LATITUDE[,1]
      }
    } else if (position == "last") {
      for (f in 1:nfloats) {
        Data[[f]]$LONGITUDE <- Data[[f]]$LONGITUDE[,ncol(Data[[f]]$LONGITUDE)]
        Data[[f]]$LATITUDE <- Data[[f]]$LATITUDE[,ncol(Data[[f]]$LATITUDE)]
      }
    }
  }
  
  # Prepare for trajectory plotting
  g1 = plot_trajectories(Data=Data, color=color, title=title)+theme(legend.position = "bottom")
  
  # Highlight specific times or cycle numbers for subduction events
  if (!is.null(subduction_events) && !is.null(highlight_value)) {
    
    # Filter subduction events based on highlight_type
    if (highlight_type == 'time') {
      sub_event_filtered <- subduction_events %>% 
        filter(TIME == highlight_value)
    } else if (highlight_type == 'cycle_number') {
      sub_event_filtered <- subduction_events %>% 
        filter(CYCLE_NUMBER == highlight_value)
    }
    
    if (nrow(sub_event_filtered) > 0) {
      # Add highlighted points for subduction events
      g1 = g1 +
        geom_point(data=sub_event_filtered, 
                   aes(x=LONGITUDE, y=LATITUDE), 
                   color="green", 
                   size=3, 
                   alpha=0.8)
    } else {
      message("No subduction events found for the given value.")
    }
  }
  
  # Return the ggplot object or display it
  if (return_ggplot) {
    return(g1)
  } else {
    x11()
    plot(g1)
  }
}

traj_5906312 <- highlight_trajectory(float_ids = 5906312,color = "multiple",
                     title = "Float Trajectories",highlight_type = "cycle_number",
                      highlight_value = 32)

ggsave(traj_5906312,
       filename = "/data/GLOBARGO/figures/FiguresForPublication/location_float_5906312.png",
       width = 12,height = 12)

ggsave(g1,filename = "/data/GLOBARGO/figures/FiguresForPublication/location_float_6901767.png",width = 12,height = 12)
