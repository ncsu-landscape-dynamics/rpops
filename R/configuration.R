# This function is designed to be a single input and output for parsing
# and outputting all checks for all functions
Sys.setenv("R_TESTS" = "")

configuration <- function(config) {

  # ensures correct model type
  if (config$model_type %in%
      c("SEI", "susceptible-exposed-infected", "susceptible_exposed_infected",
        "Susceptible-Exposed-Infected", "Susceptible_Exposed_Infected")) {
    config$model_type <- "SEI"
  } else if (config$model_type %in%
             c("SI", "susceptible-infected", "susceptible_infected",
               "Susceptible-Infected", "Susceptible_Infected")) {
    config$model_type <- "SI"
  } else {
    config$failure <-
      "Model type is not a valid type options are 'SI' or 'SEI'"
    return(config)
  }

  # ensures latent period is correct for type of model selected
  if (config$model_type == "SEI" && config$latency_period <= 0) {
    config$failure <-
      "Model type is set to SEI but the latency period is less than 1"
    return(config)
  } else if (config$model_type == "SI" && config$latency_period > 0) {
    config$latency_period <- 0
  }

  # ensure correct treatment method
  if (!config$treatment_method %in% c("ratio", "all infected")) {
    config$failure <-
      "treatment method is not one of the valid treatment options"
    return(config)
  }

  if (is.null(config$random_seed)) {
    config$random_seed <- round(stats::runif(1, 1, 1000000))
  }

  ## check output and timestep are correct.
  time_check <- time_checks(config$end_date, config$start_date,
                            config$time_step, config$output_frequency)
  if (time_check$checks_passed) {
    config$number_of_time_steps <- time_check$number_of_time_steps
    config$number_of_years <- time_check$number_of_years
    config$number_of_outputs <- time_check$number_of_outputs
    config$output_frequency <- time_check$output_frequency
    config$quarantine_frequency <- config$output_frequency
    config$quarantine_frequency_n <- config$output_frequency_n
    config$spreadrate_frequency <- config$output_frequency
    config$spreadrate_frequency_n <- config$output_frequency_n
  } else {
    config$failure <- time_check$failed_check
  }

  infected_check <- initial_raster_checks(config$infected_file)
  if (infected_check$checks_passed) {
    infected <- infected_check$raster
    if (raster::nlayers(infected) > 1) {
      infected <- output_from_raster_mean_and_sd(infected)
    }
  } else {
    config$failure <- infected_check$failed_check
    return(config)
  }

  host_check <- secondary_raster_checks(config$host_file, infected)
  if (host_check$checks_passed) {
    host <- host_check$raster
    if (raster::nlayers(host) > 1) {
      host <- output_from_raster_mean_and_sd(host)
    }
  } else {
    config$failure <- host_check$failed_check
    return(config)
  }

  total_populations_check <- secondary_raster_checks(
    config$total_populations_file, infected)
  if (total_populations_check$checks_passed) {
    total_populations <- total_populations_check$raster
    if (raster::nlayers(total_populations) > 1) {
      total_populations <- output_from_raster_mean_and_sd(total_populations)
    }
  } else {
    config$failure <- total_populations_check$failed_check
    return(config)
  }

  susceptible <- host - infected
  susceptible[susceptible < 0] <- 0

  if (config$use_lethal_temperature == TRUE) {
    temperature_check <- secondary_raster_checks(
      config$temperature_file, infected)
    if (temperature_check$checks_passed) {
      temperature_stack <- temperature_check$raster
    } else {
      config$failure <- temperature_check$failed_check
      return(config)
    }

    temperature <- list(raster::as.matrix(temperature_stack[[1]]))
    for (i in 2:config$number_of_years) {
      temperature[[i]] <- raster::as.matrix(temperature_stack[[i]])
    }
  } else {
    temperature <- host
    raster::values(temperature) <- 1
    temperature <- list(raster::as.matrix(temperature))
  }

  config$temperature <- temperature

  config$weather <- FALSE
  if (config$temp == TRUE) {
    temperature_coefficient_check <- secondary_raster_checks(
      config$temperature_coefficient_file, infected)
    if (temperature_coefficient_check$checks_passed) {
      temperature_coefficient <- temperature_coefficient_check$raster
    } else {
      config$failure <- temperature_coefficient_check$failed_check
      return(config)
    }

    config$weather <- TRUE
    weather_coefficient_stack <- temperature_coefficient
    if (config$precip == TRUE) {
      precipitation_coefficient_check <- secondary_raster_checks(
        config$precipitation_coefficient_file, infected)
      if (precipitation_coefficient_check$checks_passed) {
        precipitation_coefficient <- precipitation_coefficient_check$raster
      } else {
        config$failure <- precipitation_coefficient_check$failed_check
        return(config)
      }

      weather_coefficient_stack <- weather_coefficient_stack *
        precipitation_coefficient
    }
  } else if (config$precip == TRUE) {
    precipitation_coefficient_check <- secondary_raster_checks(
      config$precipitation_coefficient_file, infected)
    if (precipitation_coefficient_check$checks_passed) {
      precipitation_coefficient <- precipitation_coefficient_check$raster
    } else {
      config$failure <- precipitation_coefficient_check$failed_check
      return(config)
    }

    config$weather <- TRUE
    weather_coefficient_stack <- precipitation_coefficient
  }

  if (config$weather == TRUE) {
    weather_coefficient <- list(raster::as.matrix(
      weather_coefficient_stack[[1]]))
    for (i in 2:config$number_of_time_steps) {
      weather_coefficient[[i]] <- raster::as.matrix(
        weather_coefficient_stack[[i]])
    }
  } else {
    weather_coefficient <- host
    raster::values(weather_coefficient) <- 1
    weather_coefficient <- list(raster::as.matrix(weather_coefficient))
  }

  config$weather_coefficient <- weather_coefficient

  if (config$management == TRUE) {
    treatments_check <- secondary_raster_checks(config$treatments_file,
                                                infected)
    if (treatments_check$checks_passed) {
      treatment_stack <- treatments_check$raster
    } else {
      config$failure <- treatments_check$failed_check
      return(config)
    }

    treatment_check <- treatment_checks(treatment_stack, config$treatments_file,
                                        config$pesticide_duration,
                                        config$treatment_dates,
                                        config$pesticide_efficacy)
    if (treatment_check$checks_passed) {
      config$treatment_maps <- treatment_check$treatment_maps
    } else {
      config$failure <- treatment_check$failed_check
      return(config)
    }
  } else {
    treatment_map <- host
    raster::values(treatment_map) <- 0
    config$treatment_maps <- list(raster::as.matrix(treatment_map))
    config$treatment_dates <- c(config$start_date)
  }

  config$ew_res <- raster::xres(susceptible)
  config$ns_res <- raster::yres(susceptible)
  config$num_cols <- raster::ncol(susceptible)
  config$num_rows <- raster::nrow(susceptible)

  if (config$use_movements) {
    movements_check <- movement_checks(config$movements_file, infected,
                                       config$start_date, config$end_date)
    if (movements_check$checks_passed) {
      config$movements <- movements_check$movements
      config$movements_dates <- movements_check$movements_dates
    } else {
      config$failure <- movements_check$failed_check
      return(config)
    }
  } else {
    config$movements <- list(0, 0, 0, 0, 0)
    config$movements_dates <- config$start_date
  }

  mortality_tracker <- infected
  raster::values(mortality_tracker) <- 0

  infected <- raster::as.matrix(infected)
  config$susceptible <- raster::as.matrix(susceptible)
  config$total_populations <- raster::as.matrix(total_populations)
  mortality_tracker <- raster::as.matrix(mortality_tracker)
  config$mortality <- mortality_tracker
  config$resistant <- mortality_tracker
  exposed <- list(mortality_tracker)

  if (config$model_type == "SEI" & config$latency_period > 1) {
    for (ex in 2:(config$latency_period + 1)) {
      exposed[[ex]] <- mortality_tracker
    }
  }

  if (config$model_type == "SEI" & config$start_exposed) {
    exposed[[config$latency_period + 1]] <- infected
    infected <- mortality_tracker
  }

  if (config$use_quarantine) {
    quarantine_check <- secondary_raster_checks(
      config$quarantine_areas_file, host)
    if (quarantine_check$checks_passed) {
      quarantine_areas <- quarantine_check$raster
      config$quarantine_areas <- raster::as.matrix(quarantine_areas)
    } else {
      config$failure <- quarantine_check$failed_check
      return(config)
    }
  } else {
    # set quarantine areas to all zeros. meaning no quarantine areas are
    # considered
    config$quarantine_areas <- mortality_tracker
  }

  config$mortality_tracker <- mortality_tracker
  config$exposed <- exposed
  config$infected <- infected

  config$rcl <- c(1, Inf, 1, 0, 0.99, NA)
  config$rclmat <- matrix(config$rcl, ncol = 3, byrow = TRUE)

  if (config$function_name %in% c("validate", "multirun")) {
    if (is.na(config$number_of_cores) ||
        config$number_of_cores > parallel::detectCores()) {
      core_count <- parallel::detectCores() - 1
    } else {
      core_count <- config$number_of_cores
    }
    config$core_count <- core_count
  }

  if (config$function_name %in% c("validate", "pops", "multirun")) {

    if (nrow(config$parameter_cov_matrix) != 6 |
        ncol(config$parameter_cov_matrix) != 6) {
      config$failure <- "parameter covariance matrix is not 6 x 6"
      return(config)
    }

    if (length(config$parameter_means) != 6) {
      config$failure <- "parameter means is not a vector of length 6"
      return(config)
    }

    parameters <- MASS::mvrnorm(config$number_of_iterations,
                                config$parameter_means,
                                config$parameter_cov_matrix)
    while (any(parameters[, 1] < 0 |
               parameters[, 2] < 0 |
               parameters[, 3] > 1 |
               parameters[, 3] < 0 |
               parameters[, 4] < 0 |
               parameters[, 5] < 0 |
               parameters[, 6] < 0)) {
      config$number_of_draws <- nrow(parameters[parameters[, 1] < 0 |
                                                parameters[, 2] < 0 |
                                                parameters[, 3] > 1 |
                                                parameters[, 3] < 0 |
                                                parameters[, 4] < 0 |
                                                parameters[, 5] < 0 |
                                                parameters[, 6] < 0, ])

      parameters[parameters[, 1] < 0 |
                   parameters[, 2] < 0 |
                   parameters[, 3] > 1 |
                   parameters[, 3] < 0 |
                   parameters[, 4] < 0 |
                   parameters[, 5] < 0 |
                   parameters[, 6] < 0, ] <-
        MASS::mvrnorm(config$number_of_draws,
                      config$parameter_means,
                      config$parameter_cov_matrix)
    }
    config$reproductive_rate <- parameters[, 1]
    config$natural_distance_scale <- parameters[, 2]
    config$percent_natural_dispersal <- parameters[, 3]
    config$anthropogenic_distance_scale <- parameters[, 4]
    config$natural_kappa <- parameters[, 5]
    config$anthropogenic_kappa <- parameters[, 6]

    if (any(config$percent_natural_dispersal < 1.0)) {
      config$use_anthropogenic_kernel <- TRUE
    } else {
      config$use_anthropogenic_kernel <- FALSE
    }
  }

  if (config$function_name %in% c("validate", "calibrate")) {
    config$use_anthropogenic_kernel <- TRUE

    # Load observed data on occurence
    infection_years <- raster::stack(config$infected_years_file)
    infection_years[] <- as.integer(infection_years[])
    # Get rid of NA values to make comparisons
    infection_years <-
      raster::reclassify(infection_years,
                         matrix(c(NA, 0), ncol = 2, byrow = TRUE), right = NA)
    config$num_layers_infected_years <- raster::nlayers(infection_years)

    if (config$num_layers_infected_years < config$number_of_outputs) {
      config$failure <-
      paste("The infection years file must have enough layers to match the
            number of outputs from the model. The number of layers of your
            infected year file is", config$num_layers_infected_years, "and the
            number of outputs is", config$number_of_time_steps)
      return(config)
    }
    config$infection_years <- infection_years
  }

  if (config$function_name %in% c("validate") |
      (config$function_name %in% c("calibrate") &&
       config$calibration_method == "MCMC")) {
    metric_check <- metric_checks(config$success_metric)
    if (metric_check$checks_passed) {
      config$configuration <- metric_check$configuration
    } else {
      config$failure <- metric_check$failed_check
      return(config)
    }
  }

  if (config$function_name %in% c("calibrate") &&
      config$calibration_method == "ABC") {
    config$num_particles <-
      config$number_of_generations * config$generation_size

    config$total_particles <- 1
    config$current_particles <- 1
    config$proposed_particles <- 1
    config$current_bin <- 1
  }

  return(config)
  }
