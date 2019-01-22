#' Calibrates the reproductive rate and dispersal scales of the pops model.
#' 
#' Markov Chain Monte Carlo approximation is used to estimate the reproductive rate and the short distance scale parameters. Model 
#' accuracy is gauged using a custom quantity allocation disagreement function to assess accuracy of spatial configuration. The 
#' calibration uses this metric to determine if an MCMC run is kept either because it improves the results or randomly gets kept 
#' despite being worse. We recommend running calibration for at least 10,000 iterations but even more will provide a better result. 
#' If the model converges and doesn't improve for awhile it will exist calibration prior to reaching the total number of iterations specified.
#'
#' @inheritParams pops
#' @param infected_years_file years of initial infection/infestation as individual locations of a pest or pathogen in raster format
#' @param num_interations how many iterations do you want to run to allow the calibration to converge
#' @param start_reproductive_rate starting reproductive rate for MCMC calibration 
#' @param start_short_distance_scale starting short distance scale parameter for MCMC calibration
#' @param sd_reproductive_rate starting standard deviation for reproductive rate for MCMC calibration
#' @param sd_short_distance_scale starting standard deviation for short distance scale for MCMC calibration
#'
#' @importFrom raster raster values as.matrix xres yres stack reclassify cellStats nlayers
#' @importFrom  stats runif rnorm
#' @return a dataframe of the variables saved and their success metrics for each run
#' @export
#'
#' @examples
#' infected_years_file <- system.file("extdata", "SODexample", "initial_infections.tif", package = "PoPS")
#' num_iterations <- 100
#' start_reproductive_rate <- 0.5
#' start_short_distance_scale <- 20
#' sd_reproductive_rate <- 0.2
#' sd_short_distance_scale <- 1
#' infected_file <- system.file("extdata", "SODexample", "initial_infections.tif", package = "PoPS")
#' host_file <- system.file("extdata", "SODexample", "host.tif", package = "PoPS")
#' total_plants_file <- system.file("extdata", "SODexample", "all_plants.tif", package = "PoPS")
#' temperature_coefficient_file <- system.file("extdata", "SODexample", "weather.tif", package = "PoPS")
#' treatments_file <- system.file("extdata", "SODexample", "management.tif", package = "PoPS")
#' 
#' params <- calibrate(infected_years_file, num_interations, start_reproductive_rate, 
#' start_short_distance_scale, sd_reproductive_rate, sd_short_distance_scale,
#' infected_file, host_file, total_plants_file, reproductive_rate = 1.0,
#' use_lethal_temperature = FALSE, temp = TRUE, precip = FALSE, management = TRUE, mortality_on = TRUE,
#' temperature_file = "", temperature_coefficient_file, 
#' precipitation_coefficient_file ="", treatments_file,
#' season_month_start = 1, season_month_end = 12, time_step = "month",
#' start_time = 2001, end_time = 2005, treatment_years = c(2001,2002,2003,2004,2005),
#' dispersal_kern = "cauchy", percent_short_distance_dispersal = 1.0,
#' short_distance_scale = 20.57, long_distance_scale = 0.0,
#' lethal_temperature = -12.87, lethal_temperature_month = 1,
#' mortality_rate = 0.05, mortality_time_lag = 2,
#' wind_dir = "NONE", kappa = 0)
#' 
calibrate <- function(infected_years_file, num_interations, start_reproductive_rate, 
                      start_short_distance_scale, sd_reproductive_rate, sd_short_distance_scale,
                      infected_file, host_file, total_plants_file, reproductive_rate = 3.0,
                      use_lethal_temperature = FALSE, temp = FALSE, precip = FALSE, management = FALSE, mortality_on = FALSE,
                      temperature_file = "", temperature_coefficient_file = "", 
                      precipitation_coefficient_file ="", treatments_file = "",
                      season_month_start = 1, season_month_end = 12, time_step = "month",
                      start_time = 2018, end_time = 2020, treatment_years = c(0),
                      dispersal_kern = "cauchy", percent_short_distance_dispersal = 1.0,
                      short_distance_scale = 59, long_distance_scale = 0.0,
                      lethal_temperature = -12.87, lethal_temperature_month = 1,
                      mortality_rate = 0, mortality_time_lag = 0,
                      wind_dir = "NONE", kappa = 0){ 
  
  
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
    number_of_time_steps <- (end_time-start_time+1)*52+2
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
  
  ## set the parameter function to only need the parameters that chanage
  param_func <- function(reproductive_rate, short_distance_scale) {
    
    random_seed = round(stats::runif(1, 1, 1000000))
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
  
  ## Load observed data on occurence
  infection_years <- stack(infected_years_file)
  ## calculate total infections per year
  total_infections <- raster::cellStats(infection_years, 'sum')
  if (length(total_infections) > number_of_years){
    total_infections <- total_infections[1:number_of_years]
  }
  ## set up reclassification matrix for binary reclassification
  rcl <- c(1, Inf, 1, 0, 0.99, NA)
  rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
  ## reclassify to binary values
  infection_years <- raster::reclassify(infection_years, rclmat)
  ## Get rid of NA values to make comparisons
  infection_years[is.na(infection_years)] <- 0
  
  ## Create function for MCMC runs
  MCMC = function(num_iterations, start_reproductive_rate, start_short_distance_scale, sd_reproductive_rate, sd_short_distance_scale){
    params = data.frame(reproductive_rate = rep(0,num_iterations+1), short_distance_scale = rep(0,num_iterations+1), total_disagreement = rep(0,num_iterations+1), number_of_infected_difference = rep(0,num_iterations+1), directional_disagreement = rep(0,num_iterations+1))
    params$reproductive_rate[1] = start_reproductive_rate
    params$short_distance_scale[1] = start_short_distance_scale
    data <- param_func(start_reproductive_rate, start_short_distance_scale)
    reject_count = 0
    
    ## set up comparison
    comp_years <- stack(lapply(1:length(data$infected_before_treatment), function(i) infected_file))
    for (p in 1:raster::nlayers(comp_years)) {
      comp_years[[p]] <- data$infected_before_treatment[[p]]
    }

    comp_total_infections = raster::cellStats(comp_years, 'sum')
    if (length(comp_total_infections) > min(length(comp_total_infections), length(total_infections)) || length(total_infections) > min(length(comp_total_infections), length(total_infections))) {
      comp_total_infections = comp_total_infections[1:min(length(comp_total_infections), length(total_infections))]
      total_infections = total_infections[1:min(length(comp_total_infections), length(total_infections))]
    }

    comp_years <- raster::reclassify(comp_years, rclmat)
    comp_years[is.na(comp_years)] <- 0

    total_disagreement = 0
    for (j in 1:min(raster::nlayers(comp_years), raster::nlayers(infection_years))) {
      total_disagreement[j] = quantity_allocation_disagreement(infection_years[[j]], comp_years[[j]])$total_disagreement
    }

    params$total_disagreement[1] = sum(total_disagreement)
    params$number_of_infected_difference[1] = sum(abs(total_infections - comp_total_infections))
    
    current_reproductive_rate = start_reproductive_rate
    proposed_reproductive_rate =  0.0
    while (proposed_reproductive_rate <= 0) {
      proposed_reproductive_rate = round(rnorm(1,mean=current_reproductive_rate,sd= sd_reproductive_rate), digits = 1)
    }
    
    current_short_distance_scale = start_short_distance_scale
    proposed_short_distance_scale = 0.0
    while (proposed_short_distance_scale <= 0.0) {
      proposed_short_distance_scale = round(abs(rnorm(1, mean=current_short_distance_scale, sd=sd_short_distance_scale)), digits = 1)
    }
    
    params$reproductive_rate[2] <- proposed_reproductive_rate
    params$short_distance_scale <- proposed_short_distance_scale
    i = 2
    while(i <= num_iterations){

      
      data <- param_func(proposed_reproductive_rate, proposed_short_distance_scale)
      
      ## set up comparison
      comp_years <- stack(lapply(1:length(data$infected_before_treatment), function(i) infected_file))
      for (p in 1:raster::nlayers(comp_years)) {
        comp_years[[p]] <- data$infected_before_treatment[[p]]
      }
      
      comp_total_infections = raster::cellStats(comp_years, 'sum')
      if (length(comp_total_infections) > min(length(comp_total_infections), length(total_infections)) || length(total_infections) > min(length(comp_total_infections), length(total_infections))) {
        comp_total_infections = comp_total_infections[1:min(length(comp_total_infections), length(total_infections))]
        total_infections = total_infections[1:min(length(comp_total_infections), length(total_infections))]
      }
      
      comp_years <- raster::reclassify(comp_years, rclmat)
      comp_years[is.na(comp_years)] <- 0
      
      total_disagreement <- 0
      directional_disagreement <- 0
      all_disagreement = data.frame(quantity_disagreement = 0, allocation_disagreement = 0, total_disagreement = 0 , omission = 0, commission = 0 ,number_of_infected_comp =0, directional_disagreement = 0)
      for (p in 1:min(raster::nlayers(comp_years), raster::nlayers(infection_years))) {
        all_disagreement[p,] = quantity_allocation_disagreement(infection_years[[p]], comp_years[[p]])
        total_disagreement[p] = all_disagreement$total_disagreement[p]
        directional_disagreement[p] <- all_disagreement$directional_disagreement[p]
      }
      
      params$total_disagreement[i] = sum(total_disagreement)
      params$directional_disagreement[i] = sum(directional_disagreement)
      params$number_of_infected_difference[i] = sum(abs(total_infections - comp_total_infections))
      
      accept = FALSE
      
      if(params$total_disagreement[i] <= params$total_disagreement[i-1]){ # accept change if model improves or doesn't change
        current_short_distance_scale = proposed_short_distance_scale
        proposed_short_distance_scale = 0
        while (proposed_short_distance_scale <= 0) {
          if (params$directional_disagreement[i] <= 0) {
            proposed_short_distance_scale = round(current_short_distance_scale - abs(current_short_distance_scale - rnorm(1, mean = current_short_distance_scale, sd = sd_short_distance_scale)), digits = 1)
          } else if (params$directional_disagreement[i] > 0) {
            proposed_short_distance_scale = round(current_short_distance_scale + abs(current_short_distance_scale - rnorm(1, mean = current_short_distance_scale, sd = sd_short_distance_scale)), digits = 1)
          }
        }
        accept = TRUE
      } else if ((1 - ((params$total_disagreement[i] - params$total_disagreement[i-1])/(params$total_disagreement[i] + params$total_disagreement[i-1]))) < runif(1)) {
        # accept change randomly if model is worse than previous run
        current_short_distance_scale = proposed_short_distance_scale
        proposed_short_distance_scale = 0
        while (proposed_short_distance_scale <= 0) {
          if (params$directional_disagreement[i] <= 0) {
            proposed_short_distance_scale = round(current_short_distance_scale - abs(current_short_distance_scale - rnorm(1, mean = current_short_distance_scale, sd = sd_short_distance_scale)), digits = 1)
          } else if (params$directional_disagreement[i] > 0) {
            proposed_short_distance_scale = round(current_short_distance_scale + abs(current_short_distance_scale - rnorm(1, mean = current_short_distance_scale, sd = sd_short_distance_scale)), digits = 1)
          }
        }
        accept = TRUE
      } else {
        # otherwise "reject" move, and stay where we are
        proposed_short_distance_scale = current_short_distance_scale
      }
      
      if(params$number_of_infected_difference[i] <= params$number_of_infected_difference[i-1]){
        # accept change if model improves or doesn't change
        current_reproductive_rate = proposed_reproductive_rate
        proposed_reproductive_rate =  0
        while (proposed_reproductive_rate <= 0) {
          if (sum(total_infections - comp_total_infections) <= 0) {
            proposed_reproductive_rate = round(current_reproductive_rate - abs(current_reproductive_rate - rnorm(1, mean = current_reproductive_rate, sd = sd_reproductive_rate)), digits = 1)
          } else if (sum(total_infections - comp_total_infections) > 0) {
            proposed_reproductive_rate = round(current_reproductive_rate + abs(current_reproductive_rate - rnorm(1, mean = current_reproductive_rate, sd = sd_reproductive_rate)), digits = 1)
          }
        }
        accept = TRUE
      } else if ((1 - ((params$number_of_infected_difference[i] - params$number_of_infected_difference[i-1])/(params$number_of_infected_difference[i] + params$number_of_infected_difference[i-1]))) < runif(1)) {
        # accept change randomly if model is worse than previous run
        current_reproductive_rate = proposed_reproductive_rate
        proposed_reproductive_rate =  0
        while (proposed_reproductive_rate <= 0) {
          if (sum(total_infections - comp_total_infections) <= 0) {
            proposed_reproductive_rate = round(current_reproductive_rate - abs(current_reproductive_rate - rnorm(1, mean = current_reproductive_rate, sd = sd_reproductive_rate)), digits = 1)
          } else if (sum(total_infections - comp_total_infections) > 0) {
            proposed_reproductive_rate = round(current_reproductive_rate + abs(current_reproductive_rate - rnorm(1, mean = current_reproductive_rate, sd = sd_reproductive_rate)), digits = 1)
          }
        }
        accept = TRUE
      } else {
        # otherwise "reject" move, and stay where we are
        proposed_reproductive_rate = current_reproductive_rate
      }
      
      if (accept == TRUE) {
        i = i +1
        params$short_distance_scale[i] = proposed_short_distance_scale
        params$reproductive_rate[i] = proposed_reproductive_rate 
        reject_count = 0
      } else {
        i = i
        reject_count = reject_count + 1
      }
      
      if (reject_count >= 40) {
        break
      }
      print(i)
    }
    params <- params[1:i,]
    return(params)
  }
  
  
  params <- MCMC(num_iterations, start_reproductive_rate, start_short_distance_scale, sd_reproductive_rate, sd_short_distance_scale)
  
  return(params)
}

