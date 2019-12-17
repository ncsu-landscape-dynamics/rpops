#' Validates the accuracy of the calibrated reproductive rate and dispersal scales of the pops model.
#' 
#' This function uses the quantity, allocation, and configuration disagreement to validate the model
#' across the landscape using the parameters from the calibrate function. Ideally the model is calibrated
#' with 2 or more years of data and validated for the last year or if you have 6 or more years of data 
#' then the model can be validated for the final 2 years.
#'
#' @inheritParams pops
#' @param infected_years_file years of initial infection/infestation as individual locations of a pest or pathogen in raster format
#' @param num_iterations how many iterations do you want to run to allow the calibration to converge at least 10 
#' @param number_of_cores enter how many cores you want to use (default = NA). If not set uses the # of CPU cores - 1. must be an integer >= 1
#' @param success_metric Choose which success metric to use for calibration. Choices are "quantity", "quantity and configuration", "residual error" and "odds ratio". Default is "quantity"
#' @param mask Raster file used to provide a mask to remove 0's that are not true negatives from comparisons (e.g. mask out lakes and oceans from statics if modeling terrestrial species).
#'
#' @importFrom raster raster values as.matrix xres yres stack reclassify cellStats nlayers
#' @importFrom stats runif rnorm
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar%
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom iterators icount
#' @importFrom lubridate interval time_length
#' 
#' @return a dataframe of the variables saved and their success metrics for each run
#' @export
#'
#' @examples
#' \dontrun{
#' infected_years_file <- system.file("extdata", "SODexample", "initial_infection2005.tif", 
#' package = "PoPS")
#' num_iterations <- 100
#' number_of_cores <- NA
#' 
#' infected_file <- system.file("extdata", "SODexample", "initial_infection2004.tif", 
#' package = "PoPS")
#' host_file <- system.file("extdata", "SODexample", "host.tif", package = "PoPS")
#' total_plants_file <- system.file("extdata", "SODexample", "all_plants.tif", package = "PoPS")
#' temperature_coefficient_file <- system.file("extdata", "SODexample", "weather.tif", 
#' package = "PoPS")
#' treatments_file <- system.file("extdata", "SODexample", "management2005.tif", package = "PoPS")
#' 
#' params <- validate(infected_years_file, num_iterations, number_of_cores,
#' infected_file, host_file, total_plants_file, reproductive_rate = 1.0,
#' use_lethal_temperature = FALSE, temp = TRUE, precip = FALSE, management = TRUE, 
#' mortality_on = TRUE, temperature_file = "", temperature_coefficient_file, 
#' precipitation_coefficient_file ="", treatments_file,
#' season_month_start = 1, season_month_end = 12, time_step = "month",
#' start_date = '2005-01-01', end_date = 2005-12-31', treatment_dates = c(2005),
#' dispersal_kern = "cauchy", percent_short_distance_dispersal = 1.0,
#' short_distance_scale = 20.57, long_distance_scale = 0.0,
#' lethal_temperature = -12.87, lethal_temperature_month = 1,
#' mortality_rate = 0.05, mortality_time_lag = 2,
#' wind_dir = "NONE", kappa = 0)
#' }
#' 
validate <- function(infected_years_file, num_iterations, number_of_cores = NA,
                     infected_file, host_file, total_plants_file, 
                     temp = FALSE, temperature_coefficient_file = "", 
                     precip = FALSE, precipitation_coefficient_file = "", 
                     time_step = "month", reproductive_rate = 3.0,
                     season_month_start = 1, season_month_end = 12, 
                     start_date = '2008-01-01', end_date = '2008-12-31',  
                     use_lethal_temperature = FALSE, temperature_file = "",
                     lethal_temperature = -12.87, lethal_temperature_month = 1,
                     mortality_on = FALSE, mortality_rate = 0, mortality_time_lag = 0, 
                     management = FALSE, treatment_dates = c(0), treatments_file = "",
                     treatment_method = "ratio",
                     percent_natural_dispersal = 1.0,
                     natural_kernel_type = "cauchy", anthropogenic_kernel_type = "cauchy",
                     natural_distance_scale = 21, anthropogenic_distance_scale = 0.0,
                     natural_dir = "NONE", natural_kappa = 0, 
                     anthropogenic_dir = "NONE", anthropogenic_kappa = 0, 
                     pesticide_duration = 0, pesticide_efficacy = 1.0,
                     mask = NULL, success_metric = "quantity", output_frequency = "year"
                     ){ 
  
  if (success_metric == "quantity") {
    configuration = FALSE
  } else if (success_metric == "quantity and configuration") {
    configuration = TRUE
  } else if (success_metric == "odds_ratio") {
    configuration = FALSE
  } else {
    return("Success metric must be one of 'quantity', 'quantity and configuration', or 'odds_ratio'")
  }
  
  if (!treatment_method %in% c("ratio", "all infected")) {
    return("treatment method is not one of the valid treatment options")
  }
  
  if (!file.exists(infected_file)) {
    return("Infected file does not exist") 
  }
  
  if (!(raster::extension(infected_file) %in% c(".grd", ".tif", ".img"))) {
    return("Infected file is not one of '.grd', '.tif', '.img'")
  }
  
  if (!file.exists(host_file)) {
    return("Host file does not exist") 
  }
  
  if (!(raster::extension(host_file) %in% c(".grd", ".tif", ".img"))) {
    return("Host file is not one of '.grd', '.tif', '.img'")
  }
  
  if (!file.exists(total_plants_file)) {
    return("Total plants file does not exist") 
  }
  
  if (!(raster::extension(total_plants_file) %in% c(".grd", ".tif", ".img"))) {
    return("Total plants file is not one of '.grd', '.tif', '.img'")
  }
  
  if (!(time_step %in% list("week", "month", "day"))) {
    return("Time step must be one of 'week', 'month' or 'day'")
  }
  
  if (class(end_date) != "character" || class(start_date) != "character" || class(as.Date(end_date, format="%Y-%m-%d")) != "Date" || class(as.Date(start_date, format="%Y-%m-%d")) != "Date" || is.na(as.Date(end_date, format="%Y-%m-%d")) || is.na(as.Date(start_date, format="%Y-%m-%d"))){
    return("End time and/or start time not of type numeric and/or in format YYYY")
  }
  
  if (!(output_frequency %in% list("week", "month", "day", "year", "time_step"))) {
    return("Time step must be one of 'week', 'month' or 'day'")
  }
  
  if (output_frequency == "day") {
    if (time_step == "week" || time_step == "month") {
      return("Output frequency is more frequent than time_step. The minimum output_frequency you can use is the time_step of your simulation. You can set the output_frequency to 'time_step' to default to most frequent output possible")
    }
  }
  
  if (output_frequency == "week") {
    if (time_step == "month") {
      return("Output frequency is more frequent than time_step. The minimum output_frequency you can use is the time_step of your simulation. You can set the output_frequency to 'time_step' to default to most frequent output possible")
    }
  }
  
  duration <- lubridate::interval(start_date, end_date)
  
  if (time_step == "week") {
    number_of_time_steps <- ceiling(time_length(duration, "week"))
  } else if (time_step == "month") {
    number_of_time_steps <- ceiling(time_length(duration, "month"))
  } else if (time_step == "day") {
    number_of_time_steps <- ceiling(time_length(duration, "day"))
  }
  
  number_of_years <- ceiling(time_length(duration, "year"))
  
  if (output_frequency == "week") {
    number_of_outputs <- ceiling(lubridate::time_length(duration, "week"))
  } else if (output_frequency == "month") {
    number_of_outputs <- ceiling(lubridate::time_length(duration, "month"))
  } else if (output_frequency == "day") {
    number_of_outputs <- ceiling(lubridate::time_length(duration, "day"))
  } else if (output_frequency == "year") {
    number_of_outputs <- ceiling(lubridate::time_length(duration, "year"))
  } else if (output_frequency == "time_step") {
    number_of_outputs <- number_of_time_steps
  }
  
  infected <- raster::raster(infected_file)
  infected <- raster::reclassify(infected, matrix(c(NA,0), ncol = 2, byrow = TRUE), right = NA)
  host <- raster::raster(host_file)
  host <- raster::reclassify(host, matrix(c(NA,0), ncol = 2, byrow = TRUE), right = NA)
  total_plants <- raster::raster(total_plants_file)
  total_plants <- raster::reclassify(total_plants, matrix(c(NA, 0), ncol = 2, byrow = TRUE), right = NA)
  
  if (!(raster::extent(infected) == raster::extent(host) && raster::extent(infected) == raster::extent(total_plants))) {
    return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
  }
  
  if (!(raster::xres(infected) == raster::xres(host) && raster::xres(infected) == raster::xres(total_plants) && raster::yres(infected) == raster::yres(host) && raster::yres(infected) == raster::yres(total_plants))) {
    return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  }
  
  if (!(raster::compareCRS(host,infected) && raster::compareCRS(host, total_plants))) {
    return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
  }
  
  susceptible <- host - infected
  susceptible <- raster::reclassify(susceptible, matrix(c(NA,0), ncol = 2, byrow = TRUE), right = NA)
  susceptible[susceptible < 0] <- 0
  
  if (use_lethal_temperature == TRUE  && !file.exists(temperature_file)) {
    return("Temperature file does not exist")
  }
  
  if (use_lethal_temperature == TRUE  && !(raster::extension(temperature_file) %in% c(".grd", ".tif", ".img"))) {
    return("Temperature file is not one of '.grd', '.tif', '.img'")
  }
  
  if (use_lethal_temperature == TRUE) {
    temperature_stack <- raster::stack(temperature_file)
    temperature_stack <- raster::reclassify(temperature_stack, matrix(c(NA,0), ncol = 2, byrow = TRUE), right = NA)
    
    if (!(raster::extent(infected) == raster::extent(temperature_stack))) {
      return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
    }
    
    if (!(raster::xres(infected) == raster::xres(temperature_stack) && raster::yres(infected) == raster::yres(temperature_stack))) {
      return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
    }
    
    if (!(raster::compareCRS(infected, temperature_stack))) {
      return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
    }
    
    temperature <- list(raster::as.matrix(temperature_stack[[1]]))
    for(i in 2:number_of_years) {
      temperature[[i]] <- raster::as.matrix(temperature_stack[[i]])
    }
  } else {
    temperature <- host
    raster::values(temperature) <- 1
    temperature <- list(raster::as.matrix(temperature))
  }
  
  if (temp == TRUE  && !file.exists(temperature_coefficient_file)) {
    return("Temperature coefficient file does not exist")
  }
  
  if (temp == TRUE  && !(raster::extension(temperature_coefficient_file) %in% c(".grd", ".tif", ".img"))) {
    return("Temperature coefficient file is not one of '.grd', '.tif', '.img'")
  }
  
  if (precip == TRUE  && !file.exists(precipitation_coefficient_file)) {
    return("Precipitation coefficient file does not exist")
  }
  
  if (precip == TRUE  && !(raster::extension(precipitation_coefficient_file) %in% c(".grd", ".tif", ".img"))) {
    return("Precipitation coefficient file is not one of '.grd', '.tif', '.img'")
  }
  
  weather <- FALSE
  if (temp == TRUE) {
    temperature_coefficient <- raster::stack(temperature_coefficient_file)
    
    if (!(raster::extent(infected) == raster::extent(temperature_coefficient))) {
      return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
    }
    
    if (!(raster::xres(infected) == raster::xres(temperature_coefficient) && raster::yres(infected) == raster::yres(temperature_coefficient))) {
      return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
    }
    
    if (!(raster::compareCRS(infected, temperature_coefficient))) {
      return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
    }
    
    weather <- TRUE
    weather_coefficient_stack <- temperature_coefficient
    if (precip ==TRUE){
      precipitation_coefficient <- raster::stack(precipitation_coefficient_file)
      
      if (!(raster::extent(infected) == raster::extent(precipitation_coefficient))) {
        return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
      }
      
      if (!(raster::xres(infected) == raster::xres(precipitation_coefficient) && raster::yres(infected) == raster::yres(precipitation_coefficient))) {
        return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
      }
      
      if (!(raster::compareCRS(infected, precipitation_coefficient))) {
        return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
      }
      
      weather_coefficient_stack <- weather_coefficient_stack * precipitation_coefficient
    }
  } else if(precip == TRUE){
    precipitation_coefficient <- raster::stack(precipitation_coefficient_file)
    
    if (!(raster::extent(infected) == raster::extent(precipitation_coefficient))) {
      return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
    }
    
    if (!(raster::xres(infected) == raster::xres(precipitation_coefficient) && raster::yres(infected) == raster::yres(precipitation_coefficient))) {
      return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
    }
    
    if (!(raster::compareCRS(infected, precipitation_coefficient))) {
      return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
    }
    
    weather <- TRUE
    weather_coefficient_stack <- precipitation_coefficient
  }
  
  if (weather == TRUE){
    weather_coefficient_stack <- raster::reclassify(weather_coefficient_stack, matrix(c(NA,0), ncol = 2, byrow = TRUE), right = NA)
    weather_coefficient <- list(raster::as.matrix(weather_coefficient_stack[[1]]))
    for(i in 2:number_of_time_steps) {
      weather_coefficient[[i]] <- raster::as.matrix(weather_coefficient_stack[[i]])
    }
  } else {
    weather_coefficient <- host
    raster::values(weather_coefficient) <- 1
    weather_coefficient <- list(raster::as.matrix(weather_coefficient))
  }
  
  if (management == TRUE  && !file.exists(treatments_file)) {
    return("Treatments file does not exist")
  }
  
  if (management == TRUE  && !(raster::extension(treatments_file) %in% c(".grd", ".tif", ".img"))) {
    return("Treatments file is not one of '.grd', '.tif', '.img'")
  }
  
  
  if (management == TRUE) {
    treatment_stack <- raster::stack(treatments_file)
    treatment_stack <- raster::reclassify(treatment_stack, matrix(c(NA,0), ncol = 2, byrow = TRUE), right = NA)
    
    if (!(raster::extent(infected) == raster::extent(treatment_stack))) {
      return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
    }
    
    if (!(raster::xres(infected) == raster::xres(treatment_stack) && raster::yres(infected) == raster::yres(treatment_stack))) {
      return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
    }
    
    if (!(raster::compareCRS(infected, treatment_stack))) {
      return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
    }
    
    if (length(treatments_file) != length(treatment_dates)) {
      return("Length of list for treatment dates and treatments_file must be equal")
    }
    
    if (length(pesticide_duration) != length(treatment_dates)) {
      return("Length of list for treatment dates and pesticide_duration must be equal")
    }
    
    if (pesticide_duration[1] > 0) {
      treatment_maps <- list(raster::as.matrix(treatment_stack[[1]] * pesticide_efficacy))
    } else {
      treatment_maps <- list(raster::as.matrix(treatment_stack[[1]]))
    }
    if (raster::nlayers(treatment_stack) >= 2) {
      for(i in 2:raster::nlayers(treatment_stack)) {
        if (pesticide_duration[i] > 0) {
          treatment_maps[[i]] <- raster::as.matrix(treatment_stack[[i]] * pesticide_efficacy)
        } else {
          treatment_maps[[i]] <- raster::as.matrix(treatment_stack[[i]])
          
        }
      }
    }
  } else {
    treatment_map <- host
    raster::values(treatment_map) <- 0
    treatment_maps <- list(raster::as.matrix(treatment_map))
  }
  
  if(percent_natural_dispersal == 1.0) {
    use_anthropogenic_kernel = FALSE
  } else if (percent_natural_dispersal < 1.0  && percent_natural_dispersal >= 0.0) {
    use_anthropogenic_kernel = TRUE
  } else {
    return("Percent natural dispersal must be between 0.0 and 1.0")
  }
  
  ew_res <- raster::xres(susceptible)
  ns_res <- raster::yres(susceptible)
  num_cols <- raster::ncol(susceptible)
  num_rows <- raster::nrow(susceptible)
  
  mortality_tracker <- infected
  raster::values(mortality_tracker) <- 0
  
  infected <- raster::as.matrix(infected)
  susceptible <- raster::as.matrix(susceptible)
  total_plants <- raster::as.matrix(total_plants)
  mortality_tracker <- raster::as.matrix(mortality_tracker)
  mortality <- mortality_tracker
  resistant <- mortality_tracker
  
  # reference <- raster(infected_years_file)
  ## Load observed data on occurence
  infection_years <- stack(infected_years_file)
  ## set up reclassification matrix for binary reclassification
  rcl <- c(1, Inf, 1, 0, 0.99, NA)
  rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
  ## reclassify to binary values
  infection_years <- reclassify(infection_years, rclmat)
  ## Get rid of NA values to make comparisons
  infection_years <- raster::reclassify(infection_years, matrix(c(NA,0), ncol = 2, byrow = TRUE), right = NA)
  
  num_layers_infected_years <- raster::nlayers(infection_years)
  if (num_layers_infected_years < number_of_outputs) {
    return(paste("The infection years file must have enough layers to match the number of outputs from the model. The number of layers of your infected year file is", num_layers_infected_years, "and the number of outputs is", number_of_time_steps))
  }
  
  if (is.na(number_of_cores) || number_of_cores > parallel::detectCores()) {
    core_count <- parallel::detectCores() - 1
  } else {
    core_count <- number_of_cores
  }
  cl <- makeCluster(core_count)
  registerDoParallel(cl)
  

  qa <- foreach::foreach (icount(num_iterations), .combine = rbind, .packages = c("raster", "PoPS", "foreach")) %dopar% {
    random_seed <- round(stats::runif(1, 1, 1000000))
    data <- PoPS::pops_model(random_seed = random_seed, 
                             use_lethal_temperature = use_lethal_temperature, 
                             lethal_temperature = lethal_temperature, lethal_temperature_month = lethal_temperature_month,
                             infected = infected,
                             susceptible = susceptible,
                             total_plants = total_plants,
                             mortality_on = mortality_on,
                             mortality_tracker = mortality_tracker,
                             mortality = mortality,
                             treatment_maps = treatment_maps,
                             treatment_dates = treatment_dates,
                             pesticide_duration = pesticide_duration,
                             resistant = resistant,
                             weather = weather,
                             temperature = temperature,
                             weather_coefficient = weather_coefficient,
                             ew_res = ew_res, ns_res = ns_res, num_rows = num_rows, num_cols = num_cols,
                             time_step = time_step, reproductive_rate = reproductive_rate,
                             mortality_rate = mortality_rate, mortality_time_lag = mortality_time_lag,
                             season_month_start = season_month_start, season_month_end = season_month_end,
                             start_date = start_date, end_date = end_date,
                             treatment_method = treatment_method,
                             natural_kernel_type = natural_kernel_type, anthropogenic_kernel_type = anthropogenic_kernel_type, 
                             use_anthropogenic_kernel = use_anthropogenic_kernel, percent_natural_dispersal = percent_natural_dispersal,
                             natural_distance_scale = natural_distance_scale, anthropogenic_distance_scale = anthropogenic_distance_scale, 
                             natural_dir = natural_dir, natural_kappa = natural_kappa,
                             anthropogenic_dir = anthropogenic_dir, anthropogenic_kappa = anthropogenic_kappa,
                             output_frequency = output_frequency
    )
    
    comp_year <- raster(infected_file)
    all_disagreement <- foreach(q = 1:length(data$infected), .combine = rbind, .packages =c("raster", "PoPS", "foreach"), .final = colSums) %dopar% {
      comp_year[] <- data$infected[[q]]
      comp_year <- reclassify(comp_year, rclmat)
      to.all_disagreement <- quantity_allocation_disagreement(infection_years[[q]], comp_year, configuration, mask)
    }
    
    # comparison <- reference
    # raster::values(comparison) <- data$infected[[1]]
    # 
    # comparison <- raster::reclassify(comparison, rclmat)
    
    # to.qa <- PoPS::quantity_allocation_disagreement(reference, comparison)
    to.qa <- data.frame(t(all_disagreement))
  }
  
  parallel::stopCluster(cl)
  
  return(qa)
}
