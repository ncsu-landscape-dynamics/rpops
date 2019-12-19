#' Calibrates the reproductive rate and dispersal scales of the pops model.
#' 
#' Markov Chain Monte Carlo approximation is used to estimate the reproductive rate and the short distance scale parameters. Model 
#' accuracy is gauged using a custom quantity allocation disagreement function to assess accuracy of spatial configuration. The 
#' calibration uses this metric to determine if an MCMC run is kept either because it improves the results or randomly gets kept 
#' despite being worse. We recommend running calibration for at least 10,000 iterations but even more will provide a better result. 
#' If the model converges and doesn't improve for awhile it will exist calibration prior to reaching the total number of iterations specified.
#'
#' @inheritParams pops
#' @param infected_years_file Raster file with years of initial infection/infestation as individual locations of a pest or pathogen
#' @param num_iterations how many iterations do you want to run to allow the calibration to converge (recommend a minimum of at least 100,000 but preferably 1 million).
#' @param number_of_observations the number of observations used for this calibartion 
#' @param prior_number_of_observations the number of total observations from previous calibrations used to weight the posterior distributions (if this is a new calibration this value takes the form of a prior weight (0 - 1))
#' @param prior_reproductive_rate the prior reproductive rate for MCMC calibration as a list with mean and standard deviation ( e.g. c(mean, sd)) with mean and sd being numeric values or as a 2-column matrix with value and probability as columns probabilites must sum to 1
#' @param prior_natural_distance_scale the prior natural distance scale for MCMC calibration as a list with mean and standard deviation ( e.g. c(mean, sd)) with mean and sd being numeric values or as a 2-column matrix with value and probability as columns probabilites must sum to 1
#' @param prior_percent_natural_dispersal the prior percent natural distance for MCMC calibration as a list with mean and standard deviation ( e.g. c(mean, sd)) with mean and sd being numeric values or as a 2-column matrix with value and probability as columns probabilites must sum to 1
#' @param prior_anthropogenic_distance_scale the prior anthropogenic distance scale for MCMC calibration as a list with mean and standard deviation ( e.g. c(mean, sd)) with mean and sd being numeric values or as a 2-column matrix with value and probability as columns probabilites must sum to 1
#' @param number_of_cores number of cores to use for calibration (defaults to the number of cores available on the machine) integer value >= 1
#' @param success_metric Choose which success metric to use for calibration. Choices are "quantity", "quantity and configuration", "residual error" and "odds ratio". Default is "quantity"
#' @param mask Raster file used to provide a mask to remove 0's that are not true negatives from comparisons (e.g. mask out lakes and oceans from statics if modeling terrestrial species). 
#'
#' @importFrom raster raster values as.matrix xres yres stack reclassify cellStats nlayers extent extension compareCRS getValues
#' @importFrom stats runif rnorm
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar% %do% %:% foreach
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom iterators icount
#' @importFrom lubridate interval time_length
#' 
#' @return a dataframe of the variables saved and their success metrics for each run
#' 
#' @export 
#'
#' @examples
#' \dontrun{
#' }

calibrate <- function(infected_years_file, num_iterations,  number_of_cores = NA,
                      number_of_observations, prior_number_of_observations,
                      prior_reproductive_rate,
                      prior_natural_distance_scale,
                      prior_percent_natural_dispersal = c(1.0,0), 
                      prior_anthropogenic_distance_scale = c(1000,0),
                      infected_file, host_file, total_plants_file, 
                      temp = FALSE, temperature_coefficient_file = "", 
                      precip = FALSE, precipitation_coefficient_file = "", 
                      time_step = "month", 
                      season_month_start = 1, season_month_end = 12, 
                      start_date = '2008-01-01', end_date = '2008-12-31', 
                      use_lethal_temperature = FALSE, temperature_file = "",
                      lethal_temperature = -12.87, lethal_temperature_month = 1,
                      mortality_on = FALSE, mortality_rate = 0, mortality_time_lag = 0, 
                      management = FALSE, treatment_dates = c(0), treatments_file = "",
                      treatment_method = "ratio",
                      natural_kernel_type = "cauchy", anthropogenic_kernel_type = "cauchy",
                      natural_dir = "NONE", natural_kappa = 0, 
                      anthropogenic_dir = "NONE", anthropogenic_kappa = 0,
                      pesticide_duration = c(0), pesticide_efficacy = 1.0,
                      mask = NULL, success_metric = "quantity", output_frequency = "year") { 
  
  metric_check <- metric_checks(success_metric)
  if (metric_check$checks_passed){
    configuration <- metric_check$configuration
  } else {
    return(metric_check$failed_check)
  }
  
  treatment_metric_check <- treatment_metric_checks(treatment_method)
  if (!treatment_metric_check$checks_passed) {
    return(treatment_metric_check$failed_check)
  }
  
  time_check <- time_checks(end_date, start_date, time_step, output_frequency)
  if(time_check$checks_passed) {
    number_of_time_steps <- time_check$number_of_time_steps
    number_of_years <- time_check$number_of_years
    number_of_outputs <- time_check$number_of_outputs
  } else {
    return(time_check$failed_check)
  }
  
  ## Setup for reproductive rate to be passed in as either mean and sd or a 2 column data.frame with value and probability as columns 1 and 2
  if (class(prior_reproductive_rate) == "numeric" && length(prior_reproductive_rate) == 2) {
    prior_reproductive_rate <- matrix(prior_reproductive_rate, ncol = 2)
  } 
  
  if (class(prior_reproductive_rate) %in% c("matrix", "data.frame") && ncol(prior_reproductive_rate) == 2) {
    if (class(prior_reproductive_rate) == "matrix" && nrow(prior_reproductive_rate) == 1) {
      start_reproductive_rate <- prior_reproductive_rate[1]
      sd_reproductive_rate <- prior_reproductive_rate[2]
    } else if (class(prior_reproductive_rate) == "data.frame" && nrow(prior_reproductive_rate) == 1) {
      start_reproductive_rate <- prior_reproductive_rate[[1]]
      sd_reproductive_rate <- 0
    } else if (class(prior_reproductive_rate) %in% c("matrix", "data.frame") && nrow(prior_reproductive_rate) > 1) {
      names(prior_reproductive_rate) <- c('var', 'prob')
      start_reproductive_rate <- prior_reproductive_rate$var[prior_reproductive_rate$prob == max(prior_reproductive_rate$prob)]
      if(length(start_reproductive_rate) > 1) {
        start_reproductive_rate <- mean(start_reproductive_rate)
      }
      sd_reproductive_rate <- sd(prior_reproductive_rate$var)
    }
  } else {
    return("Incorrect format for prior reproductive rate")
  }
  
  ## Setup for natural dispersal scale to be passed in as either mean and sd or a 2 column data.frame with value and probability as columns 1 and 2
  if (class(prior_natural_distance_scale) == "numeric" && length(prior_natural_distance_scale) == 2) {
    prior_natural_distance_scale <- matrix(prior_natural_distance_scale, ncol = 2)
  } 
  
  if (class(prior_natural_distance_scale) %in% c("matrix", "data.frame") && ncol(prior_natural_distance_scale) == 2) {
    if (class(prior_natural_distance_scale) == "matrix" && nrow(prior_natural_distance_scale) == 1) {
      start_natural_distance_scale <- prior_natural_distance_scale[1]
      sd_natural_distance_scale <- prior_natural_distance_scale[2]
    } else if (class(prior_natural_distance_scale) == "data.frame" && nrow(prior_natural_distance_scale) == 1) {
      start_natural_distance_scale <- prior_natural_distance_scale[[1]]
      sd_natural_distance_scale <- 0
    } else if (class(prior_natural_distance_scale) %in% c("matrix", "data.frame") && nrow(prior_natural_distance_scale) > 1) {
      names(prior_natural_distance_scale) <- c('var', 'prob')
      start_natural_distance_scale <- prior_natural_distance_scale$var[prior_natural_distance_scale$prob == max(prior_natural_distance_scale$prob)]
      if(length(start_natural_distance_scale) > 1) {
        start_natural_distance_scale <- mean(start_natural_distance_scale)
      }
      sd_natural_distance_scale <- sd(prior_natural_distance_scale$var)
    }
  } else {
    return("Incorrect format for prior natural distance scale")
  }
  
  ## Setup for anthropogenic dispersal scale to be passed in as either mean and sd or a 2 column data.frame with value and probability as columns 1 and 2
  if (class(prior_anthropogenic_distance_scale) == "numeric" && length(prior_anthropogenic_distance_scale) == 2) {
    prior_anthropogenic_distance_scale <- matrix(prior_anthropogenic_distance_scale, ncol = 2)
  } 
  
  if (class(prior_anthropogenic_distance_scale) %in% c("matrix", "data.frame") && ncol(prior_anthropogenic_distance_scale) == 2) {
    if (class(prior_anthropogenic_distance_scale) == "matrix" && nrow(prior_anthropogenic_distance_scale) == 1) {
      start_anthropogenic_distance_scale <- prior_anthropogenic_distance_scale[1]
      sd_anthropogenic_distance_scale <- prior_anthropogenic_distance_scale[2]
    } else if (class(prior_anthropogenic_distance_scale) == "data.frame" && nrow(prior_anthropogenic_distance_scale) == 1) {
      start_anthropogenic_distance_scale <- prior_anthropogenic_distance_scale[[1]]
      sd_anthropogenic_distance_scale <- 0
    } else if (class(prior_anthropogenic_distance_scale) %in% c("matrix", "data.frame") && nrow(prior_anthropogenic_distance_scale) > 1) {
      names(prior_anthropogenic_distance_scale) <- c('var', 'prob')
      start_anthropogenic_distance_scale <- prior_anthropogenic_distance_scale$var[prior_anthropogenic_distance_scale$prob == max(prior_anthropogenic_distance_scale$prob)]
      if(length(start_anthropogenic_distance_scale) > 1) {
        start_anthropogenic_distance_scale <- mean(start_anthropogenic_distance_scale)
      }
      sd_anthropogenic_distance_scale <- sd(prior_anthropogenic_distance_scale$var)
    }
  } else {
    return("Incorrect format for prior athropogenic distance scale")
  }
    
  ## Setup for percent natural distance to be passed in as either mean and sd or a 2 column data.frame with value and probability as columns 1 and 2
  if (class(prior_percent_natural_dispersal) == "numeric" && length(prior_percent_natural_dispersal) == 2) {
    prior_percent_natural_dispersal <- matrix(prior_percent_natural_dispersal, ncol = 2)
  } 
    
  if (class(prior_percent_natural_dispersal) %in% c("matrix", "data.frame") && ncol(prior_percent_natural_dispersal) == 2) {
    if (class(prior_percent_natural_dispersal) == "matrix" && nrow(prior_percent_natural_dispersal) == 1) {
      start_percent_natural_dispersal <- prior_percent_natural_dispersal[1]
      sd_percent_natural_dispersal <- prior_percent_natural_dispersal[2]
    } else if (class(prior_percent_natural_dispersal) == "data.frame" && nrow(prior_percent_natural_dispersal) == 1) {
      start_percent_natural_dispersal <- prior_percent_natural_dispersal[[1]]
      sd_percent_natural_dispersal <- 0
    } else if (class(prior_percent_natural_dispersal) %in% c("matrix", "data.frame") && nrow(prior_percent_natural_dispersal) > 1) {
      names(prior_percent_natural_dispersal) <- c('var', 'prob')
      start_percent_natural_dispersal <- prior_percent_natural_dispersal$var[prior_percent_natural_dispersal$prob == max(prior_percent_natural_dispersal$prob)]
      if(length(start_percent_natural_dispersal) > 1) {
        start_percent_natural_dispersal <- mean(start_percent_natural_dispersal)
      }
      sd_percent_natural_dispersal <- sd(prior_percent_natural_dispersal$var)
    }
  } else {
    return("Incorrect format for prior percent natural distance")
  }
  
  infected_check <- initial_raster_checks(infected_file)
  if (infected_check$checks_passed) {
    infected <- infected_check$raster
    if (raster::nlayers(infected) > 1) {
      infected <- output_from_raster_mean_and_sd(infected)
    }
  } else {
    return(infected_check$failed_check)
  }
  
  host_check <- secondary_raster_checks(host_file, infected)
  if (host_check$checks_passed) {
    host <- host_check$raster
    if (raster::nlayers(host) > 1) {
      host <- output_from_raster_mean_and_sd(host)
    }
  } else {
    return(host_check$failed_check)
  }
  
  total_plants_check <- secondary_raster_checks(total_plants_file, infected)
  if (total_plants_check$checks_passed) {
    total_plants <- total_plants_check$raster
    if (raster::nlayers(total_plants) > 1) {
      total_plants <- output_from_raster_mean_and_sd(total_plants)
    }
  } else {
    return(total_plants_check$failed_check)
  }
  
  susceptible <- host - infected
  susceptible[susceptible < 0] <- 0
  
  if (use_lethal_temperature == TRUE) {
    temperature_check <- secondary_raster_checks(temperature_file, infected)
    if (temperature_check$checks_passed) {
      temperature_stack <- temperature_check$raster
    } else {
      return(temperature_check$failed_check)
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
  
  weather <- FALSE
  if (temp == TRUE) {
    temperature_coefficient_check <- secondary_raster_checks(temperature_coefficient_file, infected)
    if (temperature_coefficient_check$checks_passed) {
      temperature_coefficient <- temperature_coefficient_check$raster
    } else {
      return(temperature_coefficient_check$failed_check)
    }
    
    weather <- TRUE
    weather_coefficient_stack <- temperature_coefficient
    if (precip ==TRUE){
      precipitation_coefficient_check <- secondary_raster_checks(precipitation_coefficient_file, infected)
      if (precipitation_coefficient_check$checks_passed) {
        precipitation_coefficient <- precipitation_coefficient_check$raster
      } else {
        return(precipitation_coefficient_check$failed_check)
      }
      
      weather_coefficient_stack <- weather_coefficient_stack * precipitation_coefficient
    }
  } else if(precip == TRUE){
    precipitation_coefficient_check <- secondary_raster_checks(precipitation_coefficient_file, infected)
    if (precipitation_coefficient_check$checks_passed) {
      precipitation_coefficient <- precipitation_coefficient_check$raster
    } else {
      return(precipitation_coefficient_check$failed_check)
    }
    
    weather <- TRUE
    weather_coefficient_stack <- precipitation_coefficient
  }
  
  if (weather == TRUE){
    # weather_coefficient_stack <- raster::reclassify(weather_coefficient_stack, matrix(c(NA,0), ncol = 2, byrow = TRUE), right = NA)
    weather_coefficient <- list(raster::as.matrix(weather_coefficient_stack[[1]]))
    for(i in 2:number_of_time_steps) {
      weather_coefficient[[i]] <- raster::as.matrix(weather_coefficient_stack[[i]])
    }
  } else {
    weather_coefficient <- host
    raster::values(weather_coefficient) <- 1
    weather_coefficient <- list(raster::as.matrix(weather_coefficient))
  }
  
  if (management == TRUE) {
    treatments_check <- secondary_raster_checks(treatments_file, infected)
    if (treatments_check$checks_passed) {
      treatment_stack <- treatments_check$raster
    } else {
      return(treatments_check$failed_check)
    }
    
    treatment_check <- treatment_checks(treatment_stack, treatments_file, pesticide_duration, treatment_dates)
    if (treatment_check$checks_passed) {
      treatment_maps <- treatment_check$treatment_maps
    } else {
      return(treatment_check$failed_check)
    }
  } else {
    treatment_map <- host
    raster::values(treatment_map) <- 0
    treatment_maps <- list(raster::as.matrix(treatment_map))
  }
  
  percent_check <- percent_checks(start_percent_natural_dispersal)
  if (percent_check$checks_passed){
    use_anthropogenic_kernel <- percent_check$use_anthropogenic_kernel
  } else {
    return(percent_check$failed_check)
  }
  
  ew_res <- xres(susceptible)
  ns_res <- yres(susceptible)
  num_cols <- raster::ncol(susceptible)
  num_rows <- raster::nrow(susceptible)
  
  mortality_tracker <- infected
  mortality_tracker[] <- 0
  
  infected <- as.matrix(infected)
  susceptible <- as.matrix(susceptible)
  total_plants <- as.matrix(total_plants)
  mortality_tracker <- as.matrix(mortality_tracker)
  mortality <- mortality_tracker
  resistant <- mortality_tracker
  
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
  infection_years <- raster::reclassify(infection_years, matrix(c(NA,0), ncol = 2, byrow = TRUE), right = NA)
  num_layers_infected_years <- raster::nlayers(infection_years)
  
  if (num_layers_infected_years < number_of_outputs) {
    return(paste("The infection years file must have enough layers to match the number of outputs from the model. The number of layers of your infected year file is", num_layers_infected_years, "and the number of outputs is", number_of_time_steps))
  }
  
  ## set the parameter function to only need the parameters that chanage
  param_func <- function(reproductive_rate, natural_distance_scale, anthropogenic_distance_scale, percent_natural_dispersal) {
    random_seed <- round(runif(1, 1, 1000000))
    data <- pops_model(random_seed = random_seed, 
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
    return(data)
  }
  
  data <- param_func(start_reproductive_rate, start_natural_distance_scale, start_anthropogenic_distance_scale, start_percent_natural_dispersal)
  
  ## set up comparison
  
  comp_year <- raster(infected_file)
  all_disagreement <- foreach(q = 1:length(data$infected), .combine = rbind, .packages =c("raster", "PoPS"), .final = colSums) %do% {
    comp_year[] <- data$infected[[q]]
    comp_year <- reclassify(comp_year, rclmat)
    to.all_disagreement <- quantity_allocation_disagreement(infection_years[[q]], comp_year, configuration, mask)
  }
  
  ## save current state of the system
  current <- best <- data.frame(t(all_disagreement), reproductive_rate = start_reproductive_rate, natural_distance_scale = start_natural_distance_scale, anthropogenic_distance_scale = start_anthropogenic_distance_scale, percent_natural_dispersal = start_percent_natural_dispersal)
  
  ## create parallel environment
  if (is.na(number_of_cores) || number_of_cores > parallel::detectCores()) {
    core_count <- parallel::detectCores() - 1
  } else {
    core_count <- number_of_cores
  }
  cl <- makeCluster(core_count)
  registerDoParallel(cl)
  
  proposed_reproductive_rate <-  0
  while (proposed_reproductive_rate <= 0) {
    proposed_reproductive_rate <- round(rnorm(1, mean = best$reproductive_rate, sd = sd_reproductive_rate), digits = 1)
  }
  
  proposed_natural_distance_scale <- 0
  while (proposed_natural_distance_scale <= 0) {
    proposed_natural_distance_scale <- round(rnorm(1, mean = best$natural_distance_scale, sd = sd_natural_distance_scale), digits = 0)
  }
  
  proposed_anthropogenic_distance_scale <- 0
  while (proposed_anthropogenic_distance_scale <= 0) {
    proposed_anthropogenic_distance_scale <- round(rnorm(1, mean = best$anthropogenic_distance_scale, sd = sd_anthropogenic_distance_scale)/10, digits = 0) * 10
  }
  
  proposed_percent_natural_dispersal <- 0
  while (proposed_percent_natural_dispersal <= 0 || proposed_percent_natural_dispersal > 1.000) {
    proposed_percent_natural_dispersal <- round(rnorm(1, mean = best$percent_natural_dispersal, sd = sd_percent_natural_dispersal), digits = 3)
  }
  
  params <- foreach(icount(num_iterations), .combine = rbind, .packages = c("raster", "PoPS", "foreach", "iterators"), .inorder = TRUE) %do% {
    average_disagreements_odds_ratio <- foreach(p = 1:10, .combine = rbind, .packages = c("raster", "PoPS", "foreach"), .final = colMeans) %dopar% {
      # disagreements_odds_ratio <- data.frame(reproductive_rate = 0, natural_distance_scale = 0, total_disagreement = 0, quantity_disagreement = 0, allocation_disagreement = 0, odds_ratio = 0)
      
      data <- param_func(proposed_reproductive_rate, proposed_natural_distance_scale, proposed_anthropogenic_distance_scale, proposed_percent_natural_dispersal)
      
      # set up comparison
      all_disagreement <- foreach(q = 1:length(data$infected), .combine = rbind, .packages =c("raster", "PoPS"), .final = colSums) %dopar% {
        comp_year[] <- data$infected[[q]]
        comp_year <- reclassify(comp_year, rclmat)
        to.all_disagreement <- quantity_allocation_disagreement(infection_years[[q]], comp_year, configuration, mask)
      }
      
      to.average_disagreements_odds_ratio <- all_disagreement
    }
    
    proposed <- data.frame(t(average_disagreements_odds_ratio), reproductive_rate = proposed_reproductive_rate, natural_distance_scale = proposed_natural_distance_scale, anthropogenic_distance_scale = proposed_anthropogenic_distance_scale, percent_natural_dispersal = proposed_percent_natural_dispersal)
    
    # make sure no proposed statistics are 0 or the calculation fails instead set them all to the lowest possible non-zero value
    if (proposed$allocation_disagreement == 0) {proposed$allocation_disagreement <- 1}
    if (proposed$quantity_disagreement == 0) {proposed$quantity_disagreement <- 1}
    if (proposed$total_disagreement == 0) {proposed$total_disagreement <- 1}
    if (proposed$configuration_disagreement == 0) {proposed$configuration_disagreement <- 0.01}
    if (proposed$residual_error == 0) {proposed$residual_error <- 1}
    # Set up tests for to see if new variable is an improvement in performance metrics
    allocation_test <- min(1,  current$allocation_disagreement / proposed$allocation_disagreement)
    quantity_test <- min(1, current$quantity_disagreement / proposed$quantity_disagreement)
    total_disagreement_test <- min(1,  current$total_disagreement / proposed$total_disagreement)
    configuration_test <- min(1, current$configuration_disagreement / proposed$configuration_disagreement)
    oddsratio_test <- min(1, proposed$odds_ratio / current$odds_ratio) # odds ratio is treated differently than all the other metrics as it is the only one where higher numbers means better model performance
    residual_error_test <- min(1, current$residual_error / proposed$residual_error)
    
    quantity_pass <- runif(1) <= quantity_test
    allocation_pass <- runif(1) <= allocation_test
    total_pass <- runif(1) <= total_disagreement_test
    configuration_pass <- runif(1) <= configuration_test
    oddsratio_pass <- runif(1) <= oddsratio_test 
    residual_error_pass <- runif(1) <= residual_error_test
    
    if (success_metric == "quantity") {
      if ( quantity_pass ) {
        current <- proposed
        if (current$quantity_disagreement <= best$quantity_disagreement) {
          best <- current
        }
        param <- current
        proposed_reproductive_rate <-  0
        while (proposed_reproductive_rate <= 0) {
          proposed_reproductive_rate <- round(rnorm(1, mean = best$reproductive_rate, sd = sd_reproductive_rate), digits = 1)
        }
        
        proposed_natural_distance_scale <- 0
        while (proposed_natural_distance_scale <= 0) {
          proposed_natural_distance_scale <- round(rnorm(1, mean = best$natural_distance_scale, sd = sd_natural_distance_scale), digits = 0)
        }
        
        proposed_anthropogenic_distance_scale <- 0
        while (proposed_anthropogenic_distance_scale <= 0) {
          proposed_anthropogenic_distance_scale <- round(rnorm(1, mean = best$anthropogenic_distance_scale, sd = sd_anthropogenic_distance_scale)/10, digits = 0) * 10
        }
        
        proposed_percent_natural_dispersal <- 0
        while (proposed_percent_natural_dispersal <= 0 || proposed_percent_natural_dispersal > 1.000) {
          proposed_percent_natural_dispersal <- round(rnorm(1, mean = best$percent_natural_dispersal, sd = sd_percent_natural_dispersal), digits = 3)
        }
        
        to.params <- param
      }
    } else if (success_metric == "quantity and configuration") {
      if ( quantity_pass && configuration_pass ) {
        current <- proposed
        if (current$quantity_disagreement <= best$quantity_disagreement) {
          best <- current
        }
        param <- current
        proposed_reproductive_rate <- 0
        while (proposed_reproductive_rate <= 0) {
          proposed_reproductive_rate <- round(rnorm(1, mean = best$reproductive_rate, sd_reproductive_rate), digits = 1)
        }
        
        proposed_natural_distance_scale <- 0.0
        while (proposed_natural_distance_scale <= 0) {
          proposed_natural_distance_scale <- round(rnorm(1, mean = best$natural_distance_scale, sd_natural_distance_scale), digits = 0)
        }
        proposed_anthropogenic_distance_scale <- 0
        while (proposed_anthropogenic_distance_scale <= 0) {
          proposed_anthropogenic_distance_scale <- round(rnorm(1, mean = best$anthropogenic_distance_scale, sd = sd_anthropogenic_distance_scale)/10, digits = 0) * 10
        }
        
        proposed_percent_natural_dispersal <- 0
        while (proposed_percent_natural_dispersal <= 0 || proposed_percent_natural_dispersal > 1.000) {
          proposed_percent_natural_dispersal <- round(rnorm(1, mean = best$percent_natural_dispersal, sd = sd_percent_natural_dispersal), digits = 3)
        }
        
        to.params <- param
      } else if (configuration_pass && quantity_pass == FALSE) {
        current <- proposed
        if (current$quantity_disagreement <= best$quantity_disagreement) {
          best <- current
        }
        param <- current
        proposed_natural_distance_scale <- 0
        while (proposed_natural_distance_scale <= 0) {
          proposed_natural_distance_scale <- round(rnorm(1, mean = best$natural_distance_scale, sd_natural_distance_scale), digits = 0)
        }
        
        proposed_anthropogenic_distance_scale <- 0
        while (proposed_anthropogenic_distance_scale <= 0) {
          proposed_anthropogenic_distance_scale <- round(rnorm(1, mean = best$anthropogenic_distance_scale, sd = sd_anthropogenic_distance_scale)/10, digits = 0) * 10
        }
        
        proposed_percent_natural_dispersal <- 0
        while (proposed_percent_natural_dispersal <= 0 || proposed_percent_natural_dispersal > 1.000) {
          proposed_percent_natural_dispersal <- round(rnorm(1, mean = best$percent_natural_dispersal, sd = sd_percent_natural_dispersal), digits = 3)
        }
        to.params <- param
      } else if (quantity_pass && configuration_pass == FALSE) {
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
    } else if (success_metric == "odds ratio") {
      if ( oddsratio_pass ) { # accept change if model improves or doesn't change
        current <- proposed
        if(current$odds_ratio >= best$odds_ratio) {
          best <- current
        }
        param <- current
        proposed_reproductive_rate <-  0
        while (proposed_reproductive_rate <= 0) {
          proposed_reproductive_rate <- round(rnorm(1, mean = best$reproductive_rate, sd = sd_reproductive_rate), digits = 1)
        }
        
        proposed_natural_distance_scale <- 0
        while (proposed_natural_distance_scale <= 0) {
          proposed_natural_distance_scale <- round(rnorm(1, mean = best$natural_distance_scale, sd = sd_natural_distance_scale), digits = 0)
        }
        proposed_anthropogenic_distance_scale <- 0
        while (proposed_anthropogenic_distance_scale <= 0) {
          proposed_anthropogenic_distance_scale <- round(rnorm(1, mean = best$anthropogenic_distance_scale, sd = sd_anthropogenic_distance_scale), digits = 3)
        }
        
        proposed_percent_natural_dispersal <- 0
        while (proposed_percent_natural_dispersal <= 0 || proposed_percent_natural_dispersal > 1.000) {
          proposed_percent_natural_dispersal <- round(rnorm(1, mean = best$percent_natural_dispersal, sd = sd_percent_natural_dispersal), digits = 3)
        }
        
        to.params <- param
      } 
    } else if (success_metric == "residual error") {
      if ( residual_error_pass ) { # accept change if model improves or doesn't change
        current <- proposed
        if(current$residual_error >= best$residual_error) {
          best <- current
        }
        param <- current
        proposed_reproductive_rate <-  0
        while (proposed_reproductive_rate <= 0) {
          proposed_reproductive_rate <- round(rnorm(1, mean = best$reproductive_rate, sd = sd_reproductive_rate), digits = 1)
        }
        
        proposed_natural_distance_scale <- 0
        while (proposed_natural_distance_scale <= 0) {
          proposed_natural_distance_scale <- round(rnorm(1, mean = best$natural_distance_scale, sd = sd_natural_distance_scale), digits = 0)
        }
        proposed_anthropogenic_distance_scale <- 0
        while (proposed_anthropogenic_distance_scale <= 0) {
          proposed_anthropogenic_distance_scale <- round(rnorm(1, mean = best$anthropogenic_distance_scale, sd = sd_anthropogenic_distance_scale)/10, digits = 0) * 10
        }
        
        proposed_percent_natural_dispersal <- 0
        while (proposed_percent_natural_dispersal <= 0 || proposed_percent_natural_dispersal > 1.000) {
          proposed_percent_natural_dispersal <- round(rnorm(1, mean = best$percent_natural_dispersal, sd = sd_percent_natural_dispersal), digits = 3)
        }
        
        to.params <- param
      } 
    } else {
      return("Success metric must be one of 'quantity', 'quantity and configuration', 'residual error', or 'odds ratio'")
    }
      
  }
  stopCluster(cl)

  if (prior_number_of_observations < 1) {
    prior_weight <- prior_number_of_observations
    total_number_of_observations <- number_of_observations + round(number_of_observations * prior_number_of_observations)
    weight <- 1 - prior_weight
  } else if (prior_number_of_observations >= 1) {
    total_number_of_observations <- prior_number_of_observations + number_of_observations
    prior_weight <- prior_number_of_observations/total_number_of_observations
    weight <- 1 - prior_weight
  }
  
  ## Use prior and calibrated parameters to update to posteriors
  count <- 10000000
  # Reproductive Rate
  if (class(prior_reproductive_rate) == "matrix" && nrow(prior_reproductive_rate) == 1) {
    prior_reproductive_rates <- round(rnorm(count, start_reproductive_rate, sd_reproductive_rate), digits = 1)
    prior_reproductive_rates <- as.data.frame(table(prior_reproductive_rates))
    prior_reproductive_rates$prior_reproductive_rates <- as.numeric(as.character(prior_reproductive_rates$prior_reproductive_rates))
    prior_reproductive_rates$prob <- round(prior_reproductive_rates$Freq/count, digits = 3)
    prior_reproductive_rates <- prior_reproductive_rates[prior_reproductive_rates$prob > 0.00,]
  } else if (class(prior_reproductive_rate) %in% c("matrix", "data.frame") && nrow(prior_reproductive_rate) > 1) {
    prior_reproductive_rates <- prior_reproductive_rate[ , 1:2]
    names(prior_reproductive_rates) <- c('prior_reproductive_rates', 'prob')
    
  }

  calibration_count <- nrow(params)
  calibrated_reproductive_rates <- as.data.frame(table(params$reproductive_rate))
  calibrated_reproductive_rates$Var1 <- as.numeric(as.character(calibrated_reproductive_rates$Var1))
  calibrated_reproductive_rates$prob <- round(calibrated_reproductive_rates$Freq/calibration_count, digits = 3)
  
  min_reproductive_rate <- min(min(prior_reproductive_rates$prior_reproductive_rates), min(calibrated_reproductive_rates$Var1))
  max_reproductive_rate <- max(max(prior_reproductive_rates$prior_reproductive_rates), max(calibrated_reproductive_rates$Var1))
  
  reproductive_rates <- data.frame(reproductive_rate = round(seq(min_reproductive_rate, max_reproductive_rate, 0.1), digits = 1), 
                                   prior_probability = rep(0, length(seq(min_reproductive_rate, max_reproductive_rate, 0.1))),
                                   calibrated_probability = rep(0, length(seq(min_reproductive_rate, max_reproductive_rate, 0.1))),
                                   posterior_probability = rep(0, length(seq(min_reproductive_rate, max_reproductive_rate, 0.1))))
  
  for (i in 1:nrow(reproductive_rates)) {
    if (length(prior_reproductive_rates$prob[prior_reproductive_rates$prior_reproductive_rates == reproductive_rates$reproductive_rate[i]]) > 0) {
      reproductive_rates$prior_probability[i] <- prior_reproductive_rates$prob[prior_reproductive_rates$prior_reproductive_rates == reproductive_rates$reproductive_rate[i]]
    }
    if (length(calibrated_reproductive_rates$prob[calibrated_reproductive_rates$Var1 == reproductive_rates$reproductive_rate[i]]) > 0) {
      reproductive_rates$calibrated_probability[i] <- calibrated_reproductive_rates$prob[calibrated_reproductive_rates$Var1 == reproductive_rates$reproductive_rate[i]]
    }
  }
  reproductive_rates$posterior_probability <- round(reproductive_rates$prior_probability*prior_weight + reproductive_rates$calibrated_probability * weight, digits = 3)
  colSums(reproductive_rates)
  posterior_reproductive_rates <- reproductive_rates[,c(1,4)]
  
  # Natural Distance Scale
  if (class(prior_natural_distance_scale) == "matrix" && nrow(prior_natural_distance_scale) == 1) {
    prior_natural_distance_scales <- round(rnorm(count, start_natural_distance_scale, sd_natural_distance_scale), digits = 0)
    prior_natural_distance_scales <- as.data.frame(table(prior_natural_distance_scales))
    prior_natural_distance_scales$prior_natural_distance_scales <- as.numeric(as.character(prior_natural_distance_scales$prior_natural_distance_scales))
    prior_natural_distance_scales$prob <- round(prior_natural_distance_scales$Freq/count, digits = 3)
    prior_natural_distance_scales <- prior_natural_distance_scales[prior_natural_distance_scales$prob > 0.000,]
  } else if (class(prior_natural_distance_scale) %in% c("matrix", "data.frame") && nrow(prior_natural_distance_scale) > 1) {
    prior_natural_distance_scales <- prior_natural_distance_scale[ , 1:2]
    names(prior_natural_distance_scales) <- c('prior_natural_distance_scales', 'prob')
  }
  
  calibration_count <- nrow(params)
  calibrated_natural_distance_scales <- as.data.frame(table(params$natural_distance_scale))
  calibrated_natural_distance_scales$Var1 <- as.numeric(as.character(calibrated_natural_distance_scales$Var1))
  calibrated_natural_distance_scales$prob <- round(calibrated_natural_distance_scales$Freq/calibration_count, digits = 3)
  
  min_natural_distance_scale <- min(min(prior_natural_distance_scales$prior_natural_distance_scales), min(calibrated_natural_distance_scales$Var1))
  max_natural_distance_scale <- max(max(prior_natural_distance_scales$prior_natural_distance_scales), max(calibrated_natural_distance_scales$Var1))
  
  natural_distance_scales <- data.frame(natural_distance_scale = round(seq(min_natural_distance_scale, max_natural_distance_scale, 1), digits = 1), 
                                   prior_probability = rep(0, length(seq(min_natural_distance_scale, max_natural_distance_scale, 1))),
                                   calibrated_probability = rep(0, length(seq(min_natural_distance_scale, max_natural_distance_scale, 1))),
                                   posterior_probability = rep(0, length(seq(min_natural_distance_scale, max_natural_distance_scale, 1))))
  
  for (i in 1:nrow(natural_distance_scales)) {
    if (length(prior_natural_distance_scales$prob[prior_natural_distance_scales$prior_natural_distance_scales == natural_distance_scales$natural_distance_scale[i]]) > 0) {
      natural_distance_scales$prior_probability[i] <- prior_natural_distance_scales$prob[prior_natural_distance_scales$prior_natural_distance_scales == natural_distance_scales$natural_distance_scale[i]]
    }
    if (length(calibrated_natural_distance_scales$prob[calibrated_natural_distance_scales$Var1 == natural_distance_scales$natural_distance_scale[i]]) > 0) {
      natural_distance_scales$calibrated_probability[i] <- calibrated_natural_distance_scales$prob[calibrated_natural_distance_scales$Var1 == natural_distance_scales$natural_distance_scale[i]]
    }
  }
  natural_distance_scales$posterior_probability <- round(natural_distance_scales$prior_probability*prior_weight + natural_distance_scales$calibrated_probability * weight, digits = 3)
  colSums(natural_distance_scales)
  posterior_natural_distance_scales <- natural_distance_scales[,c(1,4)]
  
  # Anthropogenic Distance Scale
  if (class(prior_anthropogenic_distance_scale) %in% c("matrix", "data.frame") && nrow(prior_anthropogenic_distance_scale) == 1) {
    prior_anthropogenic_distance_scales <- round(rnorm(count, start_anthropogenic_distance_scale, sd_anthropogenic_distance_scale)/10, digits = 0)*10
    prior_anthropogenic_distance_scales <- as.data.frame(table(prior_anthropogenic_distance_scales))
    prior_anthropogenic_distance_scales$prior_anthropogenic_distance_scales <- as.numeric(as.character(prior_anthropogenic_distance_scales$prior_anthropogenic_distance_scales))
    prior_anthropogenic_distance_scales$prob <- round(prior_anthropogenic_distance_scales$Freq/count, digits = 3)
    prior_anthropogenic_distance_scales <- prior_anthropogenic_distance_scales[prior_anthropogenic_distance_scales$prob > 0.000,]
  } else if (class(prior_anthropogenic_distance_scale) %in% c("matrix", "data.frame") && nrow(prior_anthropogenic_distance_scale) > 1) {
    prior_anthropogenic_distance_scales <- prior_anthropogenic_distance_scale[ , 1:2]
    names(prior_anthropogenic_distance_scales) <- c('prior_anthropogenic_distance_scales', 'prob')
  }
  
  calibration_count <- nrow(params)
  calibrated_anthropogenic_distance_scales <- as.data.frame(table(params$anthropogenic_distance_scale))
  calibrated_anthropogenic_distance_scales$Var1 <- as.numeric(as.character(calibrated_anthropogenic_distance_scales$Var1))
  calibrated_anthropogenic_distance_scales$prob <- round(calibrated_anthropogenic_distance_scales$Freq/calibration_count, digits = 3)
  
  min_anthropogenic_distance_scale <- min(min(prior_anthropogenic_distance_scales$prior_anthropogenic_distance_scales), min(calibrated_anthropogenic_distance_scales$Var1))
  max_anthropogenic_distance_scale <- max(max(prior_anthropogenic_distance_scales$prior_anthropogenic_distance_scales), max(calibrated_anthropogenic_distance_scales$Var1))
  
  anthropogenic_distance_scales <- data.frame(anthropogenic_distance_scale = round(seq(min_anthropogenic_distance_scale, max_anthropogenic_distance_scale, 10), digits = 1), 
                                        prior_probability = rep(0, length(seq(min_anthropogenic_distance_scale, max_anthropogenic_distance_scale, 10))),
                                        calibrated_probability = rep(0, length(seq(min_anthropogenic_distance_scale, max_anthropogenic_distance_scale, 10))),
                                        posterior_probability = rep(0, length(seq(min_anthropogenic_distance_scale, max_anthropogenic_distance_scale, 10))))
  
  for (i in 1:nrow(anthropogenic_distance_scales)) {
    if (length(prior_anthropogenic_distance_scales$prob[prior_anthropogenic_distance_scales$prior_anthropogenic_distance_scales == anthropogenic_distance_scales$anthropogenic_distance_scale[i]]) > 0) {
      anthropogenic_distance_scales$prior_probability[i] <- prior_anthropogenic_distance_scales$prob[prior_anthropogenic_distance_scales$prior_anthropogenic_distance_scales == anthropogenic_distance_scales$anthropogenic_distance_scale[i]]
    }
    if (length(calibrated_anthropogenic_distance_scales$prob[calibrated_anthropogenic_distance_scales$Var1 == anthropogenic_distance_scales$anthropogenic_distance_scale[i]]) > 0) {
      anthropogenic_distance_scales$calibrated_probability[i] <- calibrated_anthropogenic_distance_scales$prob[calibrated_anthropogenic_distance_scales$Var1 == anthropogenic_distance_scales$anthropogenic_distance_scale[i]]
    }
  }
  anthropogenic_distance_scales$posterior_probability <- round(anthropogenic_distance_scales$prior_probability*prior_weight + anthropogenic_distance_scales$calibrated_probability * weight, digits = 3)
  colSums(anthropogenic_distance_scales)
  posterior_anthropogenic_distance_scales <- anthropogenic_distance_scales[,c(1,4)]
  
  # Percent Natural Dispersal
  if (class(prior_percent_natural_dispersal) %in% c("matrix", "data.frame") && nrow(prior_percent_natural_dispersal) == 1) {
    prior_percent_natural_dispersals <- round(rnorm(count, start_percent_natural_dispersal, sd_percent_natural_dispersal), digits = 3)
    prior_percent_natural_dispersals <- as.data.frame(table(prior_percent_natural_dispersals))
    prior_percent_natural_dispersals$prior_percent_natural_dispersals <- as.numeric(as.character(prior_percent_natural_dispersals$prior_percent_natural_dispersals))
    prior_percent_natural_dispersals$prob <- round(prior_percent_natural_dispersals$Freq/count, digits = 3)
    prior_percent_natural_dispersals <- prior_percent_natural_dispersals[prior_percent_natural_dispersals$prob > 0.000,]
    prior_percent_natural_dispersals$prob[prior_percent_natural_dispersals$prior_percent_natural_dispersals == 1.000] <- sum(prior_percent_natural_dispersals$prob[prior_percent_natural_dispersals$prior_percent_natural_dispersals >= 1.0])
    prior_percent_natural_dispersals <- prior_percent_natural_dispersals[prior_percent_natural_dispersals$prior_percent_natural_dispersals <= 1.0,]
  } else if (class(prior_percent_natural_dispersal) %in% c("matrix", "data.frame") && nrow(prior_percent_natural_dispersal) > 1) {
    prior_percent_natural_dispersals <- prior_percent_natural_dispersal[ , 1:2]
    names(prior_percent_natural_dispersals) <- c('prior_percent_natural_dispersals', 'prob')
  }
  
  calibration_count <- nrow(params)
  calibrated_percent_natural_dispersals <- as.data.frame(table(params$percent_natural_dispersal))
  calibrated_percent_natural_dispersals$Var1 <- as.numeric(as.character(calibrated_percent_natural_dispersals$Var1))
  calibrated_percent_natural_dispersals$prob <- round(calibrated_percent_natural_dispersals$Freq/calibration_count, digits = 3)
  
  min_percent_natural_dispersal <- min(min(prior_percent_natural_dispersals$prior_percent_natural_dispersals), min(calibrated_percent_natural_dispersals$Var1))
  max_percent_natural_dispersal <- max(max(prior_percent_natural_dispersals$prior_percent_natural_dispersals), max(calibrated_percent_natural_dispersals$Var1))
  
  percent_natural_dispersals <- data.frame(percent_natural_dispersal = round(seq(min_percent_natural_dispersal, max_percent_natural_dispersal, 0.001), digits = 3), 
                                        prior_probability = rep(0, length(seq(min_percent_natural_dispersal, max_percent_natural_dispersal, 0.001))),
                                        calibrated_probability = rep(0, length(seq(min_percent_natural_dispersal, max_percent_natural_dispersal, 0.001))),
                                        posterior_probability = rep(0, length(seq(min_percent_natural_dispersal, max_percent_natural_dispersal, 0.001))))
  
  for (i in 1:nrow(percent_natural_dispersals)) {
    if (length(prior_percent_natural_dispersals$prob[prior_percent_natural_dispersals$prior_percent_natural_dispersals == percent_natural_dispersals$percent_natural_dispersal[i]]) > 0) {
      percent_natural_dispersals$prior_probability[i] <- prior_percent_natural_dispersals$prob[prior_percent_natural_dispersals$prior_percent_natural_dispersals == percent_natural_dispersals$percent_natural_dispersal[i]]
    }
    if (length(calibrated_percent_natural_dispersals$prob[calibrated_percent_natural_dispersals$Var1 == percent_natural_dispersals$percent_natural_dispersal[i]]) > 0) {
      percent_natural_dispersals$calibrated_probability[i] <- calibrated_percent_natural_dispersals$prob[calibrated_percent_natural_dispersals$Var1 == percent_natural_dispersals$percent_natural_dispersal[i]]
    }
  }
  percent_natural_dispersals$posterior_probability <- round(percent_natural_dispersals$prior_probability*prior_weight + percent_natural_dispersals$calibrated_probability * weight, digits = 3)
  colSums(percent_natural_dispersals)
  posterior_percent_natural_dispersals <- percent_natural_dispersals[,c(1,4)]
  
  outputs <- list(posterior_reproductive_rates, posterior_natural_distance_scales, posterior_anthropogenic_distance_scales, posterior_percent_natural_dispersals, total_number_of_observations, reproductive_rates, natural_distance_scales, anthropogenic_distance_scales, percent_natural_dispersals, params)
  names(outputs) <- c('posterior_reproductive_rates', 'posterior_natural_distance_scales', 'posterior_anthropogenic_distance_scales', 'posterior_percent_natural_dispersals', 'total_number_of_observations', 'reproductive_rates', 'natural_distance_scales', 'anthropogenic_distance_scales', 'percent_natural_dispersals', 'raw_calibration_data')
  
  return(outputs)
}
