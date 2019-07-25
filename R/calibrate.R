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
#' @param num_iterations how many iterations do you want to run to allow the calibration to converge
#' @param start_reproductive_rate starting reproductive rate for MCMC calibration 
#' @param start_short_distance_scale starting short distance scale parameter for MCMC calibration
#' @param number_of_cores number of cores to use for calibration
#' @param sd_reproductive_rate starting standard deviation for reproductive rate for MCMC calibration
#' @param sd_short_distance_scale starting standard deviation for short distance scale for MCMC calibration
#' @param success_metric Choose which success metric to use for calibration. Choices are "quantity", "quantity and configuration", and "odds_ratio". Default = "quantity"
#' @param mask Used to provide a mask to remove 0's that are not true negatives from comparisons. 
#'
#' @importFrom raster raster values as.matrix xres yres stack reclassify cellStats nlayers extent extension compareCRS getValues
#' @importFrom stats runif rnorm
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar% %do% %:% foreach
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom iterators icount
#' 
#' @return a dataframe of the variables saved and their success metrics for each run
#' 
#' @export 
#'


calibrate <- function(infected_years_file, num_iterations, start_reproductive_rate, number_of_cores,
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
                      mortality_rate = 0, mortality_time_lag = 0, treatment_method = "ratio",
                      treatment_month = 12, wind_dir = "NONE", kappa = 0, 
                      mask = NULL, success_metric = "quantity"){ 
  
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
  
  if (!(extension(infected_file) %in% c(".grd", ".tif", ".img"))) {
    return("Infected file is not one of '.grd', '.tif', '.img'")
  }
  
  if (!file.exists(host_file)) {
    return("Host file does not exist") 
  }
  
  if (!(extension(host_file) %in% c(".grd", ".tif", ".img"))) {
    return("Host file is not one of '.grd', '.tif', '.img'")
  }
  
  if (!file.exists(total_plants_file)) {
    return("Total plants file does not exist") 
  }
  
  if (!(extension(total_plants_file) %in% c(".grd", ".tif", ".img"))) {
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
  
  infected <- raster(infected_file)
  infected[is.na(infected)] <- 0
  host <- raster(host_file)
  host[is.na(host)] <- 0
  total_plants <- raster(total_plants_file)
  total_plants[is.na(total_plants)] <- 0
  
  if (!(extent(infected) == extent(host) && extent(infected) == extent(total_plants))) {
    return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
  }
  
  if (!(xres(infected) == xres(host) && xres(infected) == xres(total_plants) && yres(infected) == yres(host) && yres(infected) == yres(total_plants))) {
    return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  }
  
  if (!(compareCRS(host,infected) && compareCRS(host, total_plants))) {
    return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
  }
  
  susceptible <- host - infected
  susceptible[is.na(susceptible)] <- 0
  susceptible[susceptible < 0] <- 0
  
  if (use_lethal_temperature == TRUE  && !file.exists(temperature_file)) {
    return("Temperature file does not exist")
  }
  
  if (use_lethal_temperature == TRUE  && !(extension(temperature_file) %in% c(".grd", ".tif", ".img"))) {
    return("Temperature file is not one of '.grd', '.tif', '.img'")
  }
  
  if (use_lethal_temperature == TRUE) {
    temperature_stack <- stack(temperature_file)
    temperature_stack[is.na(temperature_stack)] <- 0
    
    if (!(extent(infected) == extent(temperature_stack))) {
      return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
    }
    
    if (!(xres(infected) == xres(temperature_stack) && yres(infected) == yres(temperature_stack))) {
      return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
    }
    
    if (!(compareCRS(infected, temperature_stack))) {
      return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
    }
    
    temperature <- list(as.matrix(temperature_stack[[1]]))
    for(i in 2:number_of_years) {
      temperature[[i]] <- as.matrix(temperature_stack[[i]])
    }
  } else {
    temperature <- host
    temperature[] <- 1
    temperature <- list(as.matrix(temperature))
  }
  
  if (temp == TRUE  && !file.exists(temperature_coefficient_file)) {
    return("Temperature coefficient file does not exist")
  }
  
  if (temp == TRUE  && !(extension(temperature_coefficient_file) %in% c(".grd", ".tif", ".img"))) {
    return("Temperature coefficient file is not one of '.grd', '.tif', '.img'")
  }
  
  if (precip == TRUE  && !file.exists(precipitation_coefficient_file)) {
    return("Precipitation coefficient file does not exist")
  }
  
  if (precip == TRUE  && !(extension(precipitation_coefficient_file) %in% c(".grd", ".tif", ".img"))) {
    return("Precipitation coefficient file is not one of '.grd', '.tif', '.img'")
  }
  
  weather <- FALSE
  if (temp == TRUE) {
    temperature_coefficient <- stack(temperature_coefficient_file)
    
    if (!(extent(infected) == extent(temperature_coefficient))) {
      return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
    }
    
    if (!(xres(infected) == xres(temperature_coefficient) && yres(infected) == yres(temperature_coefficient))) {
      return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
    }
    
    if (!(compareCRS(infected, temperature_coefficient))) {
      return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
    }
    
    weather <- TRUE
    weather_coefficient_stack <- temperature_coefficient
    if (precip ==TRUE){
      precipitation_coefficient <- stack(precipitation_coefficient_file)
      
      if (!(extent(infected) == extent(precipitation_coefficient))) {
        return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
      }
      
      if (!(xres(infected) == xres(precipitation_coefficient) && yres(infected) == yres(precipitation_coefficient))) {
        return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
      }
      
      if (!(compareCRS(infected, precipitation_coefficient))) {
        return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
      }
      
      weather_coefficient_stack <- weather_coefficient_stack * precipitation_coefficient
    }
  } else if(precip == TRUE){
    precipitation_coefficient <- stack(precipitation_coefficient_file)
    
    if (!(extent(infected) == extent(precipitation_coefficient))) {
      return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
    }
    
    if (!(xres(infected) == xres(precipitation_coefficient) && yres(infected) == yres(precipitation_coefficient))) {
      return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
    }
    
    if (!(compareCRS(infected, precipitation_coefficient))) {
      return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
    }
    
    weather <- TRUE
    weather_coefficient_stack <- precipitation_coefficient
  }
  
  if (weather == TRUE){
    weather_coefficient_stack[is.na(weather_coefficient_stack)] <- 0
    weather_coefficient <- list(as.matrix(weather_coefficient_stack[[1]]))
    for(i in 2:number_of_time_steps) {
      weather_coefficient[[i]] <- as.matrix(weather_coefficient_stack[[i]])
    }
  } else {
    weather_coefficient <- host
    weather_coefficient[] <- 1
    weather_coefficient <- list(as.matrix(weather_coefficient))
  }
  
  if (management == TRUE  && !file.exists(treatments_file)) {
    return("Treatments file does not exist")
  }
  
  if (management == TRUE  && !(extension(treatments_file) %in% c(".grd", ".tif", ".img"))) {
    return("Treatments file is not one of '.grd', '.tif', '.img'")
  }
  
  if (management == TRUE) {
    
    treatment_stack <- stack(treatments_file)
    treatment_stack[is.na(treatment_stack)] <- 0
    
    if (!(extent(infected) == extent(treatment_stack))) {
      return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
    }
    
    if (!(xres(infected) == xres(treatment_stack) && yres(infected) == yres(treatment_stack))) {
      return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
    }
    
    if (!(compareCRS(infected, treatment_stack))) {
      return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
    }
    
    treatment_maps <- list(as.matrix(treatment_stack[[1]]))
    if (nlayers(treatment_stack) >= 2) {
      for(i in 2:nlayers(treatment_stack)) {
        treatment_maps[[i]] <- as.matrix(treatment_stack[[i]])
      }
    }
    treatment_years <- treatment_years
    treatment_method <- treatment_method
  } else {
    treatment_map <- host
    treatment_map[] <- 0
    treatment_maps <- list(as.matrix(treatment_map))
    treatment_method <- treatment_method
  }
  
  ew_res <- xres(susceptible)
  ns_res <- yres(susceptible)
  
  mortality_tracker <- infected
  mortality_tracker[] <- 0
  
  infected <- as.matrix(infected)
  susceptible <- as.matrix(susceptible)
  total_plants <- as.matrix(total_plants)
  mortality_tracker <- as.matrix(mortality_tracker)
  mortality <- mortality_tracker
  
  ## Load observed data on occurence
  infection_years <- stack(infected_years_file)
  ## calculate total infections per year in the landscape
  total_infections <- cellStats(infection_years, 'sum')
  if (length(total_infections) > number_of_years){
    total_infections <- total_infections[1:number_of_years]
  }
  ## set up reclassification matrix for binary reclassification
  rcl <- c(1, Inf, 1, 0, 0.99, NA)
  rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
  ## reclassify to binary values
  infection_years <- reclassify(infection_years, rclmat)
  ## Get rid of NA values to make comparisons
  infection_years[is.na(infection_years)] <- 0
  
  ## set the parameter function to only need the parameters that chanage
  param_func <- function(reproductive_rate, short_distance_scale) {
    random_seed <- round(runif(1, 1, 1000000))
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
                       long_distance_scale = long_distance_scale, treatment_method = treatment_method,
                       treatment_month = treatment_month, wind_dir = wind_dir, kappa = kappa)
    return(data)
  }
  
  data <- param_func(start_reproductive_rate, start_short_distance_scale)
  
  ## set up comparison
  
  comp_year <- raster(infected_file)
  all_disagreement <- foreach(q = 1:length(data$infected_before_treatment), .combine = rbind, .packages =c("raster", "PoPS"), .final = colSums) %do% {
    comp_year[] <- data$infected_before_treatment[[q]]
    comp_year <- reclassify(comp_year, rclmat)
    to.all_disagreement <- quantity_allocation_disagreement(infection_years[[q]], comp_year, configuration, mask)
  }
  
  ## save current state of the system
  current <- best <- data.frame(t(all_disagreement), reproductive_rate = start_reproductive_rate, short_distance_scale = start_short_distance_scale)
  
  ## create parallel environment
  if (is.na(number_of_cores)) {
    core_count <- detectCores() - 1
  } else {
    core_count <- number_of_cores
  }
  cl <- makeCluster(core_count)
  registerDoParallel(cl)
  
  proposed_reproductive_rate <-  0
  while (proposed_reproductive_rate <= 0) {
    proposed_reproductive_rate <- round(rnorm(1, mean = best$reproductive_rate, sd = sd_reproductive_rate), digits = 1)
  }
  
  proposed_short_distance_scale <- 0
  while (proposed_short_distance_scale <= 0) {
    proposed_short_distance_scale <- round(rnorm(1, mean = best$short_distance_scale, sd = sd_short_distance_scale), digits = 0)
  }
  
  params <- foreach(icount(num_iterations), .combine = rbind, .packages = c("raster", "PoPS", "foreach", "iterators"), .inorder = TRUE) %do% {
    average_disagreements_odds_ratio <- foreach(p = 1:10, .combine = rbind, .packages = c("raster", "PoPS", "foreach"), .final = colMeans) %dopar% {
      disagreements_odds_ratio <- data.frame(reproductive_rate = 0, short_distance_scale = 0, total_disagreement = 0, quantity_disagreement = 0, allocation_disagreement = 0, odds_ratio = 0)
      
      data <- param_func(proposed_reproductive_rate, proposed_short_distance_scale)
      
      # set up comparison
      all_disagreement <- foreach(q = 1:length(data$infected_before_treatment), .combine = rbind, .packages =c("raster", "PoPS"), .final = colSums) %dopar% {
        comp_year[] <- data$infected_before_treatment[[q]]
        comp_year <- reclassify(comp_year, rclmat)
        to.all_disagreement <- quantity_allocation_disagreement(infection_years[[q]], comp_year, configuration, mask)
      }
      
      to.average_disagreements_odds_ratio <- all_disagreement
    }
    
    proposed <- data.frame(t(average_disagreements_odds_ratio), reproductive_rate = proposed_reproductive_rate, short_distance_scale = proposed_short_distance_scale)
    
    if (proposed$allocation_disagreement == 0) {proposed$allocation_disagreement <- 1}
    if (proposed$quantity_disagreement == 0) {proposed$quantity_disagreement <- 1}
    if (proposed$total_disagreement == 0) {proposed$total_disagreement <- 1}
    if (proposed$configuration_disagreement == 0) {proposed$configuration_disagreement <- 0.01}
    allocation_test <- min(1,  current$allocation_disagreement / proposed$allocation_disagreement)
    quantity_test <- min(1, current$quantity_disagreement / proposed$quantity_disagreement)
    total_disagreement_test <- min(1,  current$total_disagreement / proposed$total_disagreement)
    configuration_test <- min(1, current$configuration_disagreement / proposed$configuration_disagreement)
    oddsratio_test <- min(1, proposed$odds_ratio / current$odds_ratio)
    
    quantity_pass <- runif(1) <= quantity_test
    allocation_pass <- runif(1) <= allocation_test
    total_pass <- runif(1) <= total_disagreement_test
    configuration_pass <- runif(1) <= configuration_test
    oddsratio_pass <- runif(1) <= oddsratio_test 
    
    if (success_metric == "quantity") {
      if ( quantity_pass == TRUE) {
        current <- proposed
        if (current$quantity_disagreement <= best$quantity_disagreement) {
          best <- current
        }
        param <- current
        proposed_reproductive_rate <-  0
        while (proposed_reproductive_rate <= 0) {
          proposed_reproductive_rate <- round(rnorm(1, mean = best$reproductive_rate, sd = sd_reproductive_rate), digits = 1)
        }
        proposed_short_distance_scale <- 0
        while (proposed_short_distance_scale <= 0) {
          proposed_short_distance_scale <- round(rnorm(1, mean = best$short_distance_scale, sd = sd_short_distance_scale), digits = 0)
        }
        to.params <- param
      }
    } else if (success_metric == "quantity and configuration") {
      if ( quantity_pass == TRUE && configuration_pass == TRUE) {
        current <- proposed
        if (current$quantity_disagreement <= best$quantity_disagreement) {
          best <- current
        }
        param <- current
        proposed_reproductive_rate <- 0
        while (proposed_reproductive_rate <= 0) {
          proposed_reproductive_rate <- round(rnorm(1, mean = best$reproductive_rate, sd_reproductive_rate), digits = 1)
        }
        proposed_short_distance_scale <- 0.0
        while (proposed_short_distance_scale <= 0) {
          proposed_short_distance_scale <- round(rnorm(1, mean = best$short_distance_scale, sd_short_distance_scale), digits = 0)
        }
        to.params <- param
      } else if (configuration_pass == TRUE && quantity_pass == FALSE) {
        current <- proposed
        if (current$quantity_disagreement <= best$quantity_disagreement) {
          best <- current
        }
        param <- current
        proposed_short_distance_scale <- 0
        while (proposed_short_distance_scale <= 0) {
          proposed_short_distance_scale <- round(rnorm(1, mean = best$short_distance_scale, sd_short_distance_scale), digits = 0)
        }
        to.params <- param
      } else if (quantity_pass == TRUE && configuration_pass == FALSE) {
        current <- proposed
        if (current$quantity_disagreement <= best$quantity_disagreement) {
          best <- current
        }
        param <- current
        proposed_reproductive_rate <-  0
        while (proposed_reproductive_rate <= 0) {
          proposed_reproductive_rate <- round(rnorm(1, mean = best$reproductive_rate, sd_reproductive_rate), digits = 1)
        }
        to.params <- param
      }
    } else if (success_metric == "odds_ratio") {
      if ( oddsratio_pass == TRUE  ) { # accept change if model improves or doesn't change
        current <- proposed
        if(current$odds_ratio >= best$odds_ratio) {
          best <- current
        }
        param <- current
        proposed_reproductive_rate <-  0
        while (proposed_reproductive_rate <= 0) {
          proposed_reproductive_rate <- round(rnorm(1, mean = best$reproductive_rate, sd = sd_reproductive_rate), digits = 1)
        }
        
        proposed_short_distance_scale <- 0
        while (proposed_short_distance_scale <= 0) {
          proposed_short_distance_scale <- round(rnorm(1, mean = best$short_distance_scale, sd = sd_short_distance_scale), digits = 0)
        }
        to.params <- param
      } 
    } else {
      return("Success metric must be one of 'quantity', 'quantity and configuration', or 'odds_ratio'")
    }
      
  }
  stopCluster(cl)

  return(params)
}