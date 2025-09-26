## Lightweight mock of load_float_data() for local testing when OneArgo is unavailable.
# The mock returns a data.frame with requested variable names and plausible shapes.
load_float_data <- function(float_ids, variables = NULL, format = "dataframe") {
  # float_ids can be single string/number or vector; we'll use the first one for metadata
  id <- if (length(float_ids) >= 1) as.character(float_ids[1]) else "MOCK001"
  set.seed(42)

  # create a small set of profiles (5 cycles) and depth levels (0:1200 every 20m)
  cycles <- seq(1, 5)
  pres <- seq(0, 1200, by = 20)

  rows <- length(cycles) * length(pres)
  df_list <- list()
  t0 <- as.numeric(Sys.time())
  for (i in seq_along(cycles)) {
    n <- length(pres)
    df_list[[i]] <- data.frame(
      PLATFORM_NUMBER = rep(id, n),
      WMO = rep(id, n),
      CYCLE_NUMBER = rep(cycles[i], n),
      PRES = pres,
      PRES_ADJUSTED = pres,
      LATITUDE = rep(0 + rnorm(1, 0, 0.1), n),
      LONGITUDE = rep(0 + rnorm(1, 0, 0.1), n),
      TIME = rep(t0 + (i * 86400), n)
    )
  }
  df <- do.call(rbind, df_list)

  # add requested variables with simple vertical structure/noise
  if (!is.null(variables)) {
    for (var in variables) {
      if (!var %in% names(df)) {
        if (grepl("PSAL", var, ignore.case = TRUE)) {
          df[[var]] <- 35 + rnorm(nrow(df), 0, 0.01) - df$PRES_ADJUSTED * 0.00008
        } else if (grepl("TEMP", var, ignore.case = TRUE)) {
          df[[var]] <- 15 - df$PRES_ADJUSTED * 0.01 + rnorm(nrow(df), 0, 0.1)
        } else if (grepl("DOXY", var, ignore.case = TRUE)) {
          df[[var]] <- 200 - df$PRES_ADJUSTED * 0.02 + rnorm(nrow(df), 0, 2)
        } else if (grepl("BBP700", var, ignore.case = TRUE)) {
          df[[var]] <- pmax(0, 0.001 * exp(-df$PRES_ADJUSTED / 200) + rnorm(nrow(df), 0, 0.0002))
        } else if (grepl("LATITUDE|LONGITUDE|TIME|PRES", var, ignore.case = TRUE)) {
          # already exists in df
        } else {
          # generic fallback for QC/other fields
          df[[var]] <- rnorm(nrow(df), 0, 1)
        }
      }
    }
  }

  if (format == "dataframe") return(df)
  return(df)
}
