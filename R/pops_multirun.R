#' PoPS (Pest or Pathogen Spread) model Multiple Runs
#' 
#' A dynamic species distribution model for pest or pathogen spread in forest or agricultural ecosystems. The model is process based
#' meaning that it uses understanding of the effect of weather on reproduction and survival of the pest/pathogen in order to forecast
#' spread of the pest/pathogen into the future. 
#'
#' @inheritParams pops
#' @param num_iterations how many iterations do you want to run to allow the calibration to converge at least 10 
#' @param number_cores enter how many cores you want to use (default = NA). If not set uses the # of CPU cores - 1. must be an integer >= 1
#'
#' @importFrom raster raster values as.matrix xres yres stack reclassify cellStats nlayers
#' @importFrom stats runif rnorm
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar%
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom iterators icount
#' @return list of infected and susceptible per year
#' @export
#'
#' @examples
#' \dontrun{
#' infected_file <-  system.file("extdata", "SODexample", "initial_infection2001.tif", 
#' package = "PoPS")
#' host_file <- system.file("extdata", "SODexample", "host.tif", package = "PoPS")
#' total_plants_file <- system.file("extdata", "SODexample", "all_plants.tif", package = "PoPS")
#' temperature_coefficient_file <- system.file("extdata", "SODexample", "weather.tif", package = "
#' PoPS")
#' treatments_file <- system.file("extdata", "SODexample", "management.tif", package = "PoPS")
#' 
#' data <- pops(infected_file, host_file, total_plants_file, reproductive_rate = 1.0,
#' use_lethal_temperature = FALSE, temp = TRUE, precip = FALSE, management = TRUE, 
#' mortality_on = TRUE, temperature_file = "", temperature_coefficient_file, 
#' precipitation_coefficient_file ="", treatments_file,
#' season_month_start = 1, season_month_end = 12, time_step = "week",
#' start_time = 2001, end_time = 2005, treatment_years = c(2001,2002,2003,2004,2005),
#' natural_kernel_type = "cauchy", percent_natural_dispersal = 1.0,
#' natural_distance_scale = 20.57, anthropogenic_distance_scale = 0.0,
#' lethal_temperature = -12.87, lethal_temperature_month = 1,
#' mortality_rate = 0.05, mortality_time_lag = 2,
#' treatment_date = 12, natural_dir = "NONE", kappa = 0, random_seed = NULL)
#' }
#' 
#' 
pops_multirun <- function(infected_file, host_file, total_plants_file, 
                 temp = FALSE, temperature_coefficient_file = "", 
                 precip = FALSE, precipitation_coefficient_file = "", 
                 time_step = "month", reproductive_rate = 3.0,
                 season_month_start = 1, season_month_end = 12, 
                 start_time = 2018, end_time = 2020, 
                 use_lethal_temperature = FALSE, temperature_file = "",
                 lethal_temperature = -12.87, lethal_temperature_month = 1,
                 mortality_on = FALSE, mortality_rate = 0, mortality_time_lag = 0, 
                 management = FALSE, treatment_years = c(0), treatments_file = "",
                 treatment_method = "ratio", treatment_month = 12,
                 percent_natural_dispersal = 1.0,
                 natural_kernel_type = "cauchy", anthropogenic_kernel_type = "cauchy",
                 natural_distance_scale = 21, anthropogenic_distance_scale = 0.0,
                 natural_dir = "NONE", natural_kappa = 0, 
                 anthropogenic_dir = "NONE", anthropogenic_kappa = 0,
                 num_iterations = 100, number_cores = NA,
                 random_seed = NULL){ 
  
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
    treatment_method = treatment_method
  } else {
    treatment_map <- host
    raster::values(treatment_map) <- 0
    treatment_maps = list(raster::as.matrix(treatment_map))
    treatment_method = treatment_method
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
  
  if (is.na(number_cores) || number_cores > detectCores()) {
    core_count <- detectCores() - 1
  } else {
    core_count <- number_cores
  }
  cl <- makeCluster(core_count)
  registerDoParallel(cl)
  years <- seq(start_time, end_time, 1)
  rcl <- c(1, Inf, 1, 0, 0.99, NA)
  rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
  
  infected_stack <- foreach::foreach(i = 1:num_iterations, .combine = c, .packages = c("raster", "PoPS"), .export = ls(globalenv())) %dopar% {
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
                             treatment_years = treatment_years,
                             weather = weather,
                             temperature = temperature,
                             weather_coefficient = weather_coefficient,
                             ew_res = ew_res, ns_res = ns_res, num_rows = num_rows, num_cols = num_cols,
                             time_step = time_step, reproductive_rate = reproductive_rate,
                             mortality_rate = mortality_rate, mortality_time_lag = mortality_time_lag,
                             season_month_start = season_month_start, season_month_end = season_month_end,
                             start_time = start_time, end_time = end_time,
                             treatment_month = treatment_month, treatment_method = treatment_method,
                             natural_kernel_type = natural_kernel_type, anthropogenic_kernel_type = anthropogenic_kernel_type, 
                             use_anthropogenic_kernel = use_anthropogenic_kernel, percent_natural_dispersal = percent_natural_dispersal,
                             natural_distance_scale = natural_distance_scale, anthropogenic_distance_scale = anthropogenic_distance_scale, 
                             natural_dir = natural_dir, natural_kappa = natural_kappa,
                             anthropogenic_dir = anthropogenic_dir, anthropogenic_kappa = anthropogenic_kappa)
    
    comp_years <- raster::stack(lapply(1:length(data$infected_before_treatment), function(i) host))
    susceptible_runs <- raster::stack(lapply(1:length(data$infected_before_treatment), function(i) host))
    
    for (q in 1:raster::nlayers(comp_years)) {
      comp_years[[q]] <- data$infected[[q]]
      susceptible_runs[[q]] <- data$susceptible[[q]]
    }
    
    number_infected <- data$number_infected
    spread_rate <- data$rates
    infected_area <- data$area_infected
    single_run <- comp_years
    comp_years <- raster::reclassify(comp_years, rclmat)
    comp_years[is.na(comp_years)] <- 0
    infected_stack <- comp_years
    data <- list(single_run, infected_stack, number_infected, susceptible_runs, infected_area, spread_rate)
  }
  
  stopCluster(cl)
  single_runs <- infected_stack[seq(1,length(infected_stack),6)]
  probability_runs <- infected_stack[seq(2,length(infected_stack),6)]
  number_infected_runs <- infected_stack[seq(3,length(infected_stack),6)]
  susceptible_runs <- infected_stack[seq(4,length(infected_stack),6)]
  area_infected_runs <- infected_stack[seq(5,length(infected_stack),6)]
  spread_rate_runs <- infected_stack[seq(6,length(infected_stack),6)]
  
  prediction <- probability_runs[[1]]
  prediction[prediction > 0] <- 0
  infected_area <- data.frame(t(years))
  infected_number <- data.frame(t(years))
  west_rates <- data.frame(t(years))
  east_rates <- data.frame(t(years))
  south_rates <- data.frame(t(years))
  north_rates <- data.frame(t(years))
  max_values <- data.frame(t(years))
  
  for (i in 1:length(probability_runs)) {
    prediction <- prediction + probability_runs[[i]]
    infected_number[i,] <- number_infected_runs[[i]]
    infected_area[i,] <- area_infected_runs[[i]]
    rates <- do.call(rbind, spread_rate_runs[[i]])
    west_rates[i,] <- rates[,4]
    east_rates[i,] <- rates[,3]
    south_rates[i,] <- rates[,2]
    north_rates[i,] <- rates[,1]
    max_values[i,] <- raster::maxValue(single_runs[[i]])
  }
  
  
  probability <- (prediction/(length(probability_runs))) * 100
  
  infected_areas <- round(sapply(infected_area, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  number_infecteds <- round(sapply(infected_number, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  west_rate <- round(sapply(west_rates, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  east_rate <- round(sapply(east_rates, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  south_rate <- round(sapply(south_rates, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  north_rate <- round(sapply(north_rates, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  
  which_median <- function(x) raster::which.min(abs(x - median(x)))
  
  median_run_index <- which_median(infected_number[[1]])
  
  single_run <- single_runs[[median_run_index]]
  susceptible_run <- susceptible_runs[[median_run_index]]
  
  single_run_out <- single_run[[1]]
  susceptible_run_out <- susceptible_run[[1]]
  
  data <- list(probability, single_run_out, number_infecteds, infected_areas, west_rate, east_rate, south_rate, north_rate)
  
  return(data)
  
}

