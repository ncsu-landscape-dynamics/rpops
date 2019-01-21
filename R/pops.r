#' PoPS (Pest or Pathogen Spread) model
#' 
#' A dynamic species distribution model for pest or pathogen spread in forest or agricultural ecosystems. The model is process based
#' meaning that it uses understanding of the effect of weather on reproduction and survival of the pest/pathogen in order to simulate
#' spread of the pest/pathogen into the future. 
#'
#' @param infected_file path to raster file with initial infections
#' @param host_file path to raster file with number of hosts
#' @param total_plants_file path to raster file with number of total plants
#' @param reproductive_rate number of spores or pest units produced by a single host under optimal weather conditions
#' @param use_lethal_temperature does your pest or pathogen have a temperature at which it cannot survive 
#' @param temp allows the use of temperature coefficients to modify spread 
#' @param precip allows the use of precipitation coefficients to modify spread
#' @param temperature_file path to raster file with temperature data for minimum temperature
#' @param temperature_coefficient_file path to raster file with 
#' @param precipitation_coefficient_file path to raster file with 
#' @param season_month_start when does spread first start occurring in the year for your pest or pathogen
#' @param season_month_end when does spread end during the year for your pest or pathogen
#' @param time_step how often should spread occur
#' @param start_time first year to start the simulation
#' @param end_time last year of the simulation
#' @param dispersal_kern what type of dispersal kernel should be used
#' @param percent_short_distance_dispersal  what percentage of dispersal is short range versus long range
#' @param short_distance_scale distance scale parameter for short range dispersal kernel
#' @param long_distance_scale distance scale parameter for long range dispersal kernel
#' @param lethal_temperature the temperature at which mortality occurs for your pest or pathogen
#' @param lethal_temperature_month the month in which mortality occurs
#' @param wind_dir sets the wind direction 
#' @param kappa sets the strength of the wind direction in the von-mises distribution
#' @param random_seed sets the random seed for the simulation used for reproducibility
#' @param management boolean to allow use of managemnet
#' @param mortality_on  boolean to turn host mortality on and off
#' @param treatments_file path to raster file with treatment data by years
#' @param treatment_years years in which to apply treatment
#' @param mortality_rate rate at which mortality occurs
#' @param mortality_time_lag time lag from infection until mortality can occur in years
#' 
#' @useDynLib PoPS, .registration = TRUE
#' @importFrom raster raster values as.matrix xres yres stack
#' @importFrom Rcpp sourceCpp evalCpp
#' @importFrom  stats runif
#' @return list of infected and susceptible per year
#' @export
#'
#' @examples This example 
#' infected_file <-  system.file("extdata", "SODexample", "initial_infections.tif", package = "PoPS")
#' host_file <- system.file("extdata", "SODexample", "host.tif", package = "PoPS")
#' total_plants_file <- system.file("extdata", "SODexample", "all_plants.tif", package = "PoPS")
#' temperature_coefficient_file <- system.file("extdata", "SODexample", "weather.tif", package = "PoPS")
#' treatments_file <- system.file("extdata", "SODexample", "management.tif", package = "PoPS")
#' 
pops <- function(infected_file, host_file, total_plants_file, reproductive_rate = 3.0,
                 use_lethal_temperature = FALSE, temp = FALSE, precip = FALSE, management = FALSE, mortality_on = FALSE,
                 temperature_file = "", temperature_coefficient_file = "", 
                 precipitation_coefficient_file ="", treatments_file = "",
                 season_month_start = 1, season_month_end = 12, time_step = "month",
                 start_time = 2018, end_time = 2020, treatment_years = c(0),
                 dispersal_kern = "cauchy", percent_short_distance_dispersal = 1.0,
                 short_distance_scale = 59, long_distance_scale = 0.0,
                 lethal_temperature = -12.87, lethal_temperature_month = 1,
                 mortality_rate = 0, mortality_time_lag = 0,
                 wind_dir = "NONE", kappa = 0, random_seed = NULL){ 
  

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
  
  if (is.null(random_seed)) {
    random_seed = round(stats::runif(1, 1, 1000000))
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

    treatment_stack <- stack(treatments_file)
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
  
  data <- pops_model(random_seed = random_seed, 
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
           long_distance_scale = long_distance_scale,
           wind_dir = wind_dir, kappa = kappa)
  
  return(data)
  
}

