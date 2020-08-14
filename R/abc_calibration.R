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
#' @param number_of_observations the number of observations used for this calibartion 
#' @param number_of_generations the number of generations to use to decrease the uncertainty in the parameter estimation (too many and it will take a long time, too few and your parameter sets will be too wide)
#' @param generation_size how many accepted parameter sets should occur in each generation
#' @param prior_number_of_observations the number of total observations from previous calibrations used to weight the posterior distributions (if this is a new calibration this value takes the form of a prior weight (0 - 1))
#' @param params_to_estimate A list of booleans specificing which parameters to estimate ordered from (reproductive_rate, natural_dispersal_distance, percent_natural_dispersal, anthropogenic_dispersal_distance, natural kappa, and anthropogenic kappa)
#' @param success_metric Choose which success metric to use for calibration. Choices are "number of locations", "number of locations and total distance", "number of locations, number of infections, and total distance", or "residual error". Default is "number of locations and total distance"
#' @param prior_means A vector of the means of your parameters you are estimating in order from (reproductive_rate, natural_dispersal_distance, percent_natural_dispersal, anthropogenic_dispersal_distance, natural kappa, and anthropogenic kappa)
#' @param prior_cov_matrix A covariance matrix from the previous years posterior parameter estimation ordered from (reproductive_rate, natural_dispersal_distance, percent_natural_dispersal, anthropogenic_dispersal_distance, natural kappa, and anthropogenic kappa)
#' @param mask Raster file used to provide a mask to remove 0's that are not true negatives from comparisons (e.g. mask out lakes and oceans from statics if modeling terrestrial species). 
#' @param checks A list of the 4 starting check values in order of # of locations, total min distance, residual error, and # infected. default is (500,500000, 100000, 1000). Starting check values can play a role in speed of calibration and in success of calibration.
#' @param natural_kappa sets the strength of the natural direction in the von-mises distribution numeric value between 0.01 and 12
#' @param anthropogenic_kappa sets the strength of the anthropogenic direction in the von-mises distribution numeric value between 0.01 and 12
#' 
#' @importFrom raster raster values as.matrix xres yres stack reclassify cellStats nlayers extent extension compareCRS getValues calc extract rasterToPoints pointDistance
#' @importFrom stats runif rnorm cov
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar% %do% %:% foreach
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom iterators icount
#' @importFrom lubridate interval time_length
#' @importFrom MASS mvrnorm
#' 
#' @return a dataframe of the variables saved and their success metrics for each run
#' 
#' @export 
#'
#' @examples
#' \dontrun{
#' }

abc_calibration <- function(infected_years_file, 
                      number_of_observations,
                      prior_number_of_observations,
                      prior_means, prior_cov_matrix, 
                      params_to_estimate = c(T, T, T, T, F, F),
                      number_of_generations = 7,
                      generation_size = 1000,
                      checks = c(500,500000, 100000, 1000),
                      infected_file,
                      host_file, 
                      total_plants_file, 
                      temp = FALSE, 
                      temperature_coefficient_file = "", 
                      precip = FALSE, 
                      precipitation_coefficient_file = "", 
                      model_type = "SI", 
                      latency_period = 0,
                      time_step = "month", 
                      season_month_start = 1,
                      season_month_end = 12, 
                      start_date = '2008-01-01',
                      end_date = '2008-12-31', 
                      use_lethal_temperature = FALSE, 
                      temperature_file = "",
                      lethal_temperature = -12.87, 
                      lethal_temperature_month = 1,
                      mortality_on = FALSE, 
                      mortality_rate = 0, 
                      mortality_time_lag = 0, 
                      management = FALSE, 
                      treatment_dates = c("2008-01-01"), 
                      treatments_file = "",
                      treatment_method = "ratio",
                      natural_kernel_type = "cauchy", 
                      anthropogenic_kernel_type = "cauchy",
                      natural_dir = "NONE", 
                      natural_kappa = 0, 
                      anthropogenic_dir = "NONE", 
                      anthropogenic_kappa = 0,
                      pesticide_duration = c(0), 
                      pesticide_efficacy = 1.0,
                      mask = NULL, 
                      success_metric = "number of locations and total distance", 
                      output_frequency = "year",
                      movements_file = "", 
                      use_movements = FALSE,
                      start_exposed = FALSE,
                      generate_stochasticity = TRUE,
                      establishment_stochasticity = TRUE,
                      movement_stochasticity = TRUE,
                      deterministic = FALSE,
                      establishment_probability = 0.5,
                      dispersal_percentage = 0.99) {
  
  if (model_type == "SEI" && latency_period <= 0) {
    return("Model type is set to SEI but the latency period is less than 1")
  } else if (model_type == "SI" && latency_period > 0) {
    latency_period <- 0
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
  
  if (use_movements) {
    movements_check <- movement_checks(movements_file, infected, start_date, end_date)
    if (movements_check$checks_passed) {
      movements <- movements_check$movements
      movements_dates <- movements_check$movements_dates
      movements_r <- movements_check$movements_r
    } else {
      return(movements_check$failed_check)
    }
  } else {
    movements <- list(0,0,0,0,0)
    movements_dates <- start_date
  }
  
  if (use_lethal_temperature == TRUE) {
    temperature_check <- secondary_raster_checks(temperature_file, infected)
    if (temperature_check$checks_passed) {
      temperature_stack <- temperature_check$raster
    } else {
      return(temperature_check$failed_check)
    }
    
    temperature <- list(as.matrix(temperature_stack[[1]]))
    if (number_of_years > 1) {
      for(i in 2:number_of_years) {
        temperature[[i]] <- as.matrix(temperature_stack[[i]])
      }
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
    
    treatment_check <- treatment_checks(treatment_stack, treatments_file, pesticide_duration, treatment_dates, pesticide_efficacy)
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
  exposed <- list(mortality_tracker)
  
  if (latency_period > 1){
    for (ex in 2:(latency_period + 1)) {
      exposed[[ex]] <- mortality_tracker
    }
  }
  
  if (model_type == "SEI" & start_exposed) {
    exposed[[latency_period + 1]] <- infected
    infected <- mortality_tracker
  }
  
  ## Load observed data on occurence
  infection_years <- stack(infected_years_file)
  infection_years[] <- as.integer(infection_years[])
  ## Get rid of NA values to make comparisons
  infection_years <- raster::reclassify(infection_years, matrix(c(NA,0), ncol = 2, byrow = TRUE), right = NA)
  num_layers_infected_years <- raster::nlayers(infection_years)
  
  if (num_layers_infected_years < number_of_outputs) {
    return(paste("The infection years file must have enough layers to match the number of outputs from the model. The number of layers of your infected year file is", num_layers_infected_years, "and the number of outputs is", number_of_time_steps))
  }
  
  use_anthropogenic_kernel <- TRUE
  ## set the parameter function to only need the parameters that chanage
  param_func <- function(reproductive_rate, natural_distance_scale, anthropogenic_distance_scale, percent_natural_dispersal, natural_kappa, anthropogenic_kappa) {
    random_seed <- round(runif(1, 1, 1000000))
    data <- PoPS::pops_model(random_seed = random_seed, 
                             use_lethal_temperature = use_lethal_temperature, 
                             lethal_temperature = lethal_temperature, 
                             lethal_temperature_month = lethal_temperature_month,
                             infected = infected,
                             exposed = exposed,
                             susceptible = susceptible,
                             total_plants = total_plants,
                             mortality_on = mortality_on,
                             mortality_tracker = mortality_tracker,
                             mortality = mortality,
                             treatment_maps = treatment_maps,
                             treatment_dates = treatment_dates,
                             pesticide_duration = pesticide_duration,
                             resistant = resistant,
                             use_movements = use_movements, 
                             movements = movements,
                             movements_dates = movements_dates,
                             weather = weather,
                             temperature = temperature,
                             weather_coefficient = weather_coefficient,
                             ew_res = ew_res, 
                             ns_res = ns_res, 
                             num_rows = num_rows, 
                             num_cols = num_cols,
                             time_step = time_step, 
                             reproductive_rate = reproductive_rate,
                             mortality_rate = mortality_rate, 
                             mortality_time_lag = mortality_time_lag,
                             season_month_start = season_month_start, 
                             season_month_end = season_month_end,
                             start_date = start_date, 
                             end_date = end_date,
                             treatment_method = treatment_method,
                             natural_kernel_type = natural_kernel_type, 
                             anthropogenic_kernel_type = anthropogenic_kernel_type, 
                             use_anthropogenic_kernel = use_anthropogenic_kernel, 
                             percent_natural_dispersal = percent_natural_dispersal,
                             natural_distance_scale = natural_distance_scale, 
                             anthropogenic_distance_scale = anthropogenic_distance_scale, 
                             natural_dir = natural_dir, 
                             natural_kappa = natural_kappa,
                             anthropogenic_dir = anthropogenic_dir, 
                             anthropogenic_kappa = anthropogenic_kappa,
                             output_frequency = output_frequency,
                             model_type_ = model_type,
                             latency_period = latency_period,
                             generate_stochasticity = generate_stochasticity,
                             establishment_stochasticity = establishment_stochasticity,
                             movement_stochasticity = movement_stochasticity,
                             deterministic = deterministic,
                             establishment_probability = establishment_probability,
                             dispersal_percentage = dispersal_percentage
    )
    return(data)
  }
  
  num_particles <- number_of_generations * generation_size
  
  total_particles <- 1
  current_particles <- 1
  proposed_particles <- 1
  current_bin <- 1
  
  parameters_kept <- matrix(ncol = 10, nrow = num_particles)
  acc_rate <- 1
  acc_rates <- matrix(ncol = 1, nrow = number_of_generations)
  infected_checks <- matrix(ncol = 1, nrow = number_of_generations)
  locs_checks <- matrix(ncol = 1, nrow = number_of_generations)
  dist_checks <- matrix(ncol = 1, nrow = number_of_generations)
  res_error_checks <- matrix(ncol = 1, nrow = number_of_generations)
  
  ## calculate comparison metrics for input data (still need to add in configuration metrics to this)
  num_infected_data <-  c()
  num_infected_data <- length(number_of_outputs)
  num_locs_data <- c()
  num_locs_data <- length(number_of_outputs)
  infected_data_points <- vector(mode = "list", length = number_of_outputs)

  for (y in 1:nlayers(infection_years)) {
    inf_year <- infection_years[[y]]
    num_infected_data[[y]] <- sum(inf_year[inf_year > 0])
    num_locs_data[[y]] <- sum(inf_year[inf_year > 0] > 0)
    if (success_metric %in% c("number of locations and total distance", "number of locations, number of infections, and total distance")) {
      infected_data_points[[y]] <- rasterToPoints(inf_year, fun=function(x){x>0}, spatial = TRUE)
    }
  }

  ## calculate total infections per output in the landscape
  total_infections <- cellStats(infection_years, 'sum')
  if (length(total_infections) > number_of_outputs){
    total_infections <- total_infections[1:number_of_outputs]
  }
  
  locs_check <- checks[1]
  dist_check <- checks[2]
  res_error_check <- checks[3]
  inf_check <- checks[4]
  infected_sims <- infection_years
  infected_sim <- infection_years[[1]]

  while (current_bin <= number_of_generations) {
    
    while (current_particles <= generation_size){

      if (current_bin == 1){
        proposed_reproductive_rate <- round(runif(1, 0.055, 8), digits = 2)
        proposed_natural_distance_scale <- round(runif(1, 0.5, 100), digits = 1)
        if (params_to_estimate[3]) {
          proposed_percent_natural_dispersal <- round(runif(1, 0.93, 1), digits = 3)
        } else {
          proposed_percent_natural_dispersal <- 1.0
        }
        if (params_to_estimate[4]) {
          proposed_anthropogenic_distance_scale <- round(runif(1, 30, 80), digits = 0)*100
        } else {
          proposed_anthropogenic_distance_scale <- 0
        }
        if (params_to_estimate[5]) {
          proposed_natural_kappa <- round(runif(1, 0, 5), digits = 1)
        } else {
          proposed_natural_kappa <- natural_kappa
        }
        if (params_to_estimate[6]) {
          proposed_anthropogenic_kappa <- round(runif(1, 0, 5), digits = 1)
        } else {
          proposed_anthropogenic_kappa <- anthropogenic_kappa
        }
      } else {
        proposed_parameters <- mvrnorm(1, parameter_means, parameter_cov_matrix)
        while(proposed_parameters[1] <= 0 || proposed_parameters[2] <= 0) {
          proposed_parameters <- mvrnorm(1, parameter_means, parameter_cov_matrix)
        }
        proposed_reproductive_rate <- proposed_parameters[1]
        proposed_natural_distance_scale <- proposed_parameters[2]
        if (params_to_estimate[3]) {
          proposed_percent_natural_dispersal <- proposed_parameters[3]
          if (proposed_percent_natural_dispersal > 1.000) {proposed_percent_natural_dispersal <- 1.000} 
        } else {
          proposed_percent_natural_dispersal <- 1.0
        }
        if (params_to_estimate[4]) {
          proposed_anthropogenic_distance_scale <- proposed_parameters[4]
        } else {
          proposed_anthropogenic_distance_scale <- 0
        }
        if (params_to_estimate[5]) {
          proposed_natural_kappa <- proposed_parameters[5]
          if (proposed_natural_kappa <= 0.000) {proposed_natural_kappa <- 0}
        } else {
          proposed_natural_kappa <- natural_kappa
        }
        if (params_to_estimate[6]) {
          proposed_anthropogenic_kappa <- proposed_parameters[6]
          if (proposed_anthropogenic_kappa <= 0.000) {proposed_anthropogenic_kappa <- 0}
        } else {
          proposed_anthropogenic_kappa <- anthropogenic_kappa
        }
      }
      
      data <- param_func(proposed_reproductive_rate, proposed_natural_distance_scale, proposed_anthropogenic_distance_scale, proposed_percent_natural_dispersal, proposed_natural_kappa, proposed_anthropogenic_kappa)
      
      ## calculate comparison metrics for simulation data (still need to add in configuration metrics to this)
      num_infected_simulated <-  c()
      num_infected_simulated <- length(number_of_outputs)
      num_locs_simulated <- c()
      num_locs_simulated <- length(number_of_outputs)
      infected_sim_points <- vector(mode = "list", length = number_of_outputs)
      dist <- vector(mode = "list", length = number_of_outputs)
      dist_diffs <- vector(mode = "list", length = number_of_outputs)
      residual_diffs <- c()
      residual_diffs <- length(number_of_outputs)
      
      for (y in 1:nlayers(infection_years)) {
        if (nlayers(infection_years) > 1) {
          infected_sims[[y]][] <- data$infected[[y]]
          infected_sim[] <- data$infected[[y]]
        } else {
          infected_sim[] <- data$infected[[y]]
        }
        
        if (!is.null(mask)){
          infected_sim[is.na(mask)] <- 0
        }
        
        diff_raster <- infection_years[[y]] - infected_sim
        residual_diffs[[y]] <- abs(sum(diff_raster[diff_raster != 0]))
        
        num_infected_simulated[[y]] <- sum(infected_sim[infected_sim > 0])
        num_locs_simulated[[y]] <- sum(infected_sim[infected_sim > 0] > 0)
        if (success_metric %in% c("number of locations and total distance", "number of locations, number of infections, and total distance")) {
          infected_sim_points[[y]] <- rasterToPoints(infected_sim, fun=function(x){x>0}, spatial = TRUE)
          dist[[y]] <- pointDistance(infected_sim_points[[y]], infected_data_points[[y]], lonlat = FALSE)
          if (class(dist) == "matrix") {
            dist_diffs[[y]] <- apply(dist[[y]], 2, min)
          } else {
            dist_diffs[[y]] <- dist[[y]]
          }
        }
      }

      if (success_metric %in% c("number of locations and total distance", "number of locations, number of infections, and total distance")) {
        all_distances <- function(dist_diffs){
          dist_diffs <- round(sqrt(sum(dist_diffs^2)), digits = 0)
          return(dist_diffs)
        }
        dist_diffs <- lapply(dist_diffs, all_distances)
        dist_diffs <- unlist(dist_diffs, recursive = TRUE, use.names = TRUE)
      } else {
        dist_diffs <- 0
      }
      
      
      num_differences <- sqrt((num_infected_data - num_infected_simulated)^2)
      locs_diffs <- sqrt((num_locs_data - num_locs_simulated)^2)
      
      num_difference <- sum(num_differences)
      locs_diff <- sum(locs_diffs)
      residual_diff <- sum(residual_diffs)
      dist_diff <- sum(dist_diffs)
      
      diff_checks <- FALSE
      if (success_metric == "number of locations and total distance") {
        if (locs_diff <= locs_check && dist_diff <= dist_check) {
          diff_checks <- TRUE
        }
      } else if (success_metric == "number of locations"){
        if (locs_diff <= locs_check) {
          diff_checks <- TRUE
        }
      } else if (success_metric == "residual error") {
        if (residual_diff <= res_error_check) {
          diff_checks <- TRUE
        }
      } else if (success_metric == "number of locations, number of infections, and total distance") {
        if (locs_diff <= locs_check && dist_diff <= dist_check && num_difference <= inf_check) {
          diff_checks <- TRUE
        }
      } else {
        return("success metric must be one of 'number of locations and total distance', 'number of locations', and 'residual error'")
      }
      
      if (diff_checks) {
        parameters_kept[total_particles, ] <- c(proposed_reproductive_rate, proposed_natural_distance_scale, proposed_percent_natural_dispersal, proposed_anthropogenic_distance_scale, proposed_natural_kappa, proposed_anthropogenic_kappa, num_difference, locs_diff, dist_diff, residual_diff)
        current_particles <- current_particles + 1
        total_particles <- total_particles + 1
        proposed_particles <- proposed_particles + 1
      } else {
        proposed_particles <- proposed_particles + 1
      }
      acc_rate <- current_particles/proposed_particles
      print(paste("The current generation is ", current_bin, " and the current particle is ", current_particles, " and the current acceptance rate is ", acc_rate, sep = ""))
    }
    
    start_index <- current_bin * generation_size - generation_size + 1
    end_index <- current_bin * generation_size
    parameter_means <- colMeans(parameters_kept[start_index:end_index,1:6])
    parameter_cov_matrix <- cov(parameters_kept[start_index:end_index,1:6])
    reproductive_rate_generation <- as.data.frame(table(parameters_kept[start_index:end_index, 1]))
    natural_distance_scale_generation <- as.data.frame(table(parameters_kept[start_index:end_index, 2]))
    percent_natural_dispersal_generation <- as.data.frame(table(parameters_kept[start_index:end_index, 3]))
    anthro_dis_scale <- parameters_kept[start_index:end_index, 3:4]
    anthro_dis_scale <- anthro_dis_scale[anthro_dis_scale[,1] < 1.000,]
    anthropogenic_distance_scale_generation <- as.data.frame(table(anthro_dis_scale[,2]))
    
    reproductive_rate_generation$Freq <- reproductive_rate_generation$Freq/generation_size
    natural_distance_scale_generation$Freq <- natural_distance_scale_generation$Freq/generation_size
    percent_natural_dispersal_generation$Freq <- percent_natural_dispersal_generation$Freq/generation_size
    anthropogenic_distance_scale_generation$Freq <- anthropogenic_distance_scale_generation$Freq/nrow(anthro_dis_scale[anthro_dis_scale[,1] < 1.000,])
    
    reproductive_rate_generation$Var1 <- as.numeric(as.character(reproductive_rate_generation$Var1))
    natural_distance_scale_generation$Var1 <- as.numeric(as.character(natural_distance_scale_generation$Var1))
    percent_natural_dispersal_generation$Var1 <- as.numeric(as.character(percent_natural_dispersal_generation$Var1))
    anthropogenic_distance_scale_generation$Var1 <- as.numeric(as.character(anthropogenic_distance_scale_generation$Var1))
    
    current_particles <- 1
    proposed_particles <- 1
    acc_rates[current_bin] <- acc_rate
    infected_checks[current_bin] <- inf_check
    locs_checks[current_bin] <- locs_check
    dist_checks[current_bin] <- dist_check
    inf_check <- median(parameters_kept[start_index:end_index, 7])
    locs_check <- median(parameters_kept[start_index:end_index, 8])
    dist_check <- median(parameters_kept[start_index:end_index, 9])
    res_error_check <- median(parameters_kept[start_index:end_index, 10])
    current_bin <- current_bin + 1
  }

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
  calibrated_means <- colMeans(parameters_kept[start_index:end_index,1:6])
  calibrated_cov_matrix <- cov(parameters_kept[start_index:end_index,1:6])
  posterior_check <- bayesian_MNN_checks(prior_means, prior_cov_matrix, calibrated_means, calibrated_cov_matrix, prior_weight, weight)
  
  if(posterior_check$checks_passed) {
    posterior_means <- posterior_check$posterior_means
    posterior_cov_matrix <- posterior_check$posterior_cov_matrix
  } else {
    return(posterior_check$failed_check)
  }

  outputs <- list(posterior_means, posterior_cov_matrix, total_number_of_observations, parameters_kept)
  names(outputs) <- c('posterior_means', 'posterior_cov_matrix', 'total_number_of_observations', 'raw_calibration_data')

  return(outputs)
}


