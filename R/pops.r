#' PoPS (Pest or Pathogen Spread) model
#'
#' @param infected_file path to raster file with initial infections
#' @param host_file path to raster file with number of hosts
#' @param total_plants_file path to raster file with number of total plants
#' @param reproductive_rate 
#' @param use_lethal_temperature 
#' @param temp 
#' @param precip 
#' @param temperature_file path to raster file with temperature data for minimum temperature
#' @param temperature_coefficient_file path to raster file with 
#' @param precipitation_coefficient_file path to raster file with 
#' @param season_month_start 
#' @param season_month_end 
#' @param time_step 
#' @param start_time 
#' @param end_time 
#' @param dispersal_kern 
#' @param percent_short_distance_dispersal 
#' @param short_distance_scale 
#' @param long_distance_scale 
#' @param lethal_temperature 
#' @param lethal_temperature_month 
#' @param wind_dir sets the wind direction 
#' @param kappa sets the strength of the wind direction in the von-mises distribution
#' @param random_seed sets the random seed for the simulation used for reproducibility
#'
#' @importFrom raster raster values as.matrix xres yres stack
#' @importFrom Rcpp sourceCpp
#' @return 
#' @export
#'
#' @examples None
#' 
pops <- function(infected_file, host_file, total_plants_file, reproductive_rate = 3.0,
                 use_lethal_temperature = FALSE, temp = FALSE, precip = FALSE,
                 temperature_file = "", temperature_coefficient_file = "", precipitation_coefficient_file ="",
                 season_month_start = 1, season_month_end = 12, time_step = "month",
                 start_time = 2018, end_time = 2022,
                 dispersal_kern = "cauchy", percent_short_distance_dispersal = 1.0,
                 short_distance_scale = 59, long_distance_scale = 0.0,
                 lethal_temperature = -12.87, lethal_temperature_month = 1,
                 wind_dir = "NONE", kappa = 0, random_seed = 42){ 
  
  Sys.setenv("PKF_CXXFLAGS"="-std=c++11")
  Rcpp::sourceCpp("C:/Users/Chris/Desktop/rpops/pops.cpp")
  
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
    number_of_time_steps <- (end_time-start_time+1)*52 +1
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
  susceptible <- host - infected
  susceptible[is.na(susceptible)] <- 0
  total_plants <- raster::raster(total_plants_file)
  total_plants[is.na(total_plants)] <- 0
  
  if (use_lethal_temperature == TRUE  && !file.exists(temperature_file)) {
    return("Temperature file does not exist")
  }
  
  if (use_lethal_temperature == TRUE  && !(raster::extension(temperature_file) %in% c(".grd", ".tif", ".img"))) {
    return("Temperature file is not one of '.grd', '.tif', '.img'")
  }
  
  if (use_lethal_temperature == TRUE) {
    temperature_stack <- raster::stack(temperature_file)
    temperature_stack[is.na(temperature_stack)] <- 0
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
    weather <- TRUE
    weather_coefficient_stack <- temperature_coefficient
    if (precip ==TRUE){
      precipitation_coefficient <- raster::stack(precipitation_coefficient_file)
      weather_coefficient_stack <- weather_coefficient_stack * precipitation_coefficient
    }
  } else if(precip == TRUE){
    precipitation_coefficient <- raster::stack(precipitation_coefficient_file)
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
  
  ew_res <- raster::xres(susceptible)
  ns_res <- raster::yres(susceptible)
  
  mortality_tracker <- infected
  raster::values(mortality_tracker) <- 0
  
  infected <- raster::as.matrix(infected)
  susceptible <- raster::as.matrix(susceptible)
  total_plants <- raster::as.matrix(total_plants)
  mortality_tracker <- raster::as.matrix(mortality_tracker)
  
  data <- pops_model(random_seed = random_seed, 
           lethal_temperature = lethal_temperature, use_lethal_temperature = use_lethal_temperature, lethal_temperature_month = lethal_temperature_month,
           reproductive_rate = reproductive_rate, 
           weather = weather, short_distance_scale = short_distance_scale, infected = infected,
           susceptible = susceptible, mortality_tracker = mortality_tracker,
           total_plants = total_plants, temperature = temperature,
           weather_coefficient = weather_coefficient, 
           ew_res = ew_res, ns_res = ns_res,
           time_step = time_step,
           season_month_start = season_month_start, season_month_end = season_month_end,
           start_time = start_time, end_time = end_time,
           dispersal_kern = dispersal_kern, percent_short_distance_dispersal = percent_short_distance_dispersal,
           long_distance_scale = long_distance_scale,
           wind_dir = wind_dir, kappa = kappa)
  
  return(data)
  
}

