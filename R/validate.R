#' Validates the accuracy of the calibrated reproductive rate and dispersal scales of the pops model.
#' 
#' This function uses the quantity, allocation, and configuration disagreement to validate the model
#' across the landscape using the parameters from the calibrate function. Ideally the model is calibrated
#' with 2 or more years of data and validated for the last year or if you have 6 or more years of data 
#' then the model can be validated for the final 2 years.
#'
#' @inheritParams pops
#' @param infected_years_file years of initial infection/infestation as individual locations of a pest or pathogen in raster format
#' @param num_iterations how many iterations do you want to run to allow the calibration to converge
#' @param number_cores enter how many cores you want to use (default = NA). If not set uses the # of CPU cores - 1.
#'
#' @importFrom raster raster values as.matrix xres yres stack reclassify cellStats nlayers
#' @importFrom stats runif rnorm
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar%
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom iterators icount
#' @return a dataframe of the variables saved and their success metrics for each run
#' @export
#'
#' @examples
#' infected_years_file <- system.file("extdata", "SODexample", "initial_infection2005.tif", 
#' package = "PoPS")
#' num_iterations <- 100
#' number_cores <- NA
#' 
#' infected_file <- system.file("extdata", "SODexample", "initial_infection2004.tif", 
#' package = "PoPS")
#' host_file <- system.file("extdata", "SODexample", "host.tif", package = "PoPS")
#' total_plants_file <- system.file("extdata", "SODexample", "all_plants.tif", package = "PoPS")
#' temperature_coefficient_file <- system.file("extdata", "SODexample", "weather.tif", 
#' package = "PoPS")
#' treatments_file <- system.file("extdata", "SODexample", "management2005.tif", package = "PoPS")
#' 
#' params <- validate(infected_years_file, num_iterations, number_cores,
#' infected_file, host_file, total_plants_file, reproductive_rate = 1.0,
#' use_lethal_temperature = FALSE, temp = TRUE, precip = FALSE, management = TRUE, 
#' mortality_on = TRUE, temperature_file = "", temperature_coefficient_file, 
#' precipitation_coefficient_file ="", treatments_file,
#' season_month_start = 1, season_month_end = 12, time_step = "month",
#' start_time = 2005, end_time = 2005, treatment_years = c(2005),
#' dispersal_kern = "cauchy", percent_short_distance_dispersal = 1.0,
#' short_distance_scale = 20.57, long_distance_scale = 0.0,
#' lethal_temperature = -12.87, lethal_temperature_month = 1,
#' mortality_rate = 0.05, mortality_time_lag = 2,
#' wind_dir = "NONE", kappa = 0)
#' 
validate <- function(infected_years_file, num_iterations, number_cores = NA,
                      infected_file, host_file, total_plants_file, reproductive_rate = 3.0,
                      use_lethal_temperature = FALSE, temp = FALSE, precip = FALSE, management = FALSE, mortality_on = FALSE,
                      temperature_file = "", temperature_coefficient_file = "", 
                      precipitation_coefficient_file ="", treatments_file = "",
                      season_month_start = 1, season_month_end = 12, time_step = "month",
                      start_time = 2018, end_time = 2020, treatment_years = c(0),
                      dispersal_kern = "cauchy", percent_short_distance_dispersal = 1.0,
                      short_distance_scale = 59, long_distance_scale = 0.0,
                      lethal_temperature = -12.87, lethal_temperature_month = 1,
                      mortality_rate = 0, mortality_time_lag = 0, treatment_method = "ratio",
                      wind_dir = "NONE", kappa = 0){ 
  
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
  
  if (class(end_time) != "numeric" || nchar(end_time) != 4 || class(start_time) != "numeric" || nchar(start_time) != 4){
    return("End time and/or start time not of type numeric and/or in format YYYY")
  }
  
  if (time_step == "week") {
    number_of_time_steps <- (end_time-start_time+1)*52
  } else if (time_step == "month") {
    number_of_time_steps <- (end_time-start_time+1)*12
  } else if (time_step == "day") {
    number_of_time_steps <- (end_time-start_time+1)*365
  }
  
  number_of_years <- end_time-start_time+1
  
  infected <- raster::raster(infected_file)
  infected[is.na(infected)] <- 0
  host <- raster::raster(host_file)
  host[is.na(host)] <- 0
  total_plants <- raster::raster(total_plants_file)
  total_plants[is.na(total_plants)] <- 0
  
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
  susceptible[is.na(susceptible)] <- 0
  susceptible[susceptible < 0] <- 0
  
  if (use_lethal_temperature == TRUE  && !file.exists(temperature_file)) {
    return("Temperature file does not exist")
  }
  
  if (use_lethal_temperature == TRUE  && !(raster::extension(temperature_file) %in% c(".grd", ".tif", ".img"))) {
    return("Temperature file is not one of '.grd', '.tif', '.img'")
  }
  
  if (use_lethal_temperature == TRUE) {
    temperature_stack <- raster::stack(temperature_file)
    temperature_stack[is.na(temperature_stack)] <- 0
    
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
    weather_coefficient_stack[is.na(weather_coefficient_stack)] <- 0
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
    treatment_stack[is.na(treatment_stack)] <- 0
    
    if (!(raster::extent(infected) == raster::extent(treatment_stack))) {
      return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
    }
    
    if (!(raster::xres(infected) == raster::xres(treatment_stack) && raster::yres(infected) == raster::yres(treatment_stack))) {
      return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
    }
    
    if (!(raster::compareCRS(infected, treatment_stack))) {
      return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
    }
    
    treatment_maps <- list(raster::as.matrix(treatment_stack[[1]]))
    if (raster::nlayers(treatment_stack) >= 2) {
      for(i in 2:raster::nlayers(treatment_stack)) {
        treatment_maps[[i]] <- raster::as.matrix(treatment_stack[[i]])
      }
    }
    treatment_years = treatment_years
  } else {
    treatment_map <- host
    raster::values(treatment_map) <- 0
    treatment_maps = list(raster::as.matrix(treatment_map))
  }
  
  ew_res <- raster::xres(susceptible)
  ns_res <- raster::yres(susceptible)
  
  mortality_tracker <- infected
  raster::values(mortality_tracker) <- 0
  
  infected <- raster::as.matrix(infected)
  susceptible <- raster::as.matrix(susceptible)
  total_plants <- raster::as.matrix(total_plants)
  mortality_tracker <- raster::as.matrix(mortality_tracker)
  mortality <- mortality_tracker
  
  reference <- raster(infected_years_file)
  
  core_count <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(core_count)
  doParallel::registerDoParallel(cl)

  qa <- foreach::foreach (icount(num_iterations), .combine = rbind, .packages = c("raster", "PoPS")) %dopar% {
    random_seed <- round(stats::runif(1, 1, 1000000))
    data <- PoPS::pops_model(random_seed = random_seed, 
                       lethal_temperature = lethal_temperature, use_lethal_temperature = use_lethal_temperature, lethal_temperature_month = lethal_temperature_month,
                       reproductive_rate = reproductive_rate, 
                       weather = weather, mortality_on = mortality_on,
                       short_distance_scale = short_distance_scale, infected = infected,
                       susceptible = susceptible, mortality_tracker = mortality_tracker, mortality = mortality,
                       total_plants = total_plants, 
                       treatment_maps = treatment_maps, treatment_years = treatment_years,
                       temperature = temperature,
                       weather_coefficient = weather_coefficient, 
                       ew_res = ew_res, ns_res = ns_res,
                       time_step = time_step, mortality_rate = mortality_rate, mortality_time_lag = mortality_time_lag,
                       season_month_start = season_month_start, season_month_end = season_month_end,
                       start_time = start_time, end_time = end_time,
                       dispersal_kern = dispersal_kern, percent_short_distance_dispersal = percent_short_distance_dispersal,
                       long_distance_scale = long_distance_scale, treatment_method = treatment_method,
                       wind_dir = wind_dir, kappa = kappa)
    
    comparison <- reference
    raster::values(comparison) <- data$infected[[1]]
    
    rcl <- c(1, Inf, 1, 0, 0.99, NA)
    rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
    ## reclassify to binary values
    reference <- raster::reclassify(reference, rclmat)
    comparison <- raster::reclassify(comparison, rclmat)
    
    to.qa <- PoPS::quantity_allocation_disagreement(reference, comparison)
  }
  
  parallel::stopCluster(cl)
  
  return(qa)
}
