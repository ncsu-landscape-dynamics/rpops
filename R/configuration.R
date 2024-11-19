#' PoPS (configuration
#'
#' Function with a single input and output list for parsing, transforming,
#' and performing all checks for all functions to run the pops c++ model
#'
#' @param config list of all data necessary used to set up c++ model
#'
#' @importFrom terra app rast xres yres classify extract ext as.points ncol nrow
#' nlyr rowFromCell colFromCell values as.matrix rowFromCell colFromCell crs
#' rowColFromCell global vect
#' @importFrom stats runif rnorm median sd
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar% %do%
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom lubridate interval time_length mdy %within%
#' @importFrom aws.s3 head_object save_object
#'
#' @return config list with all data ready for pops C++ or error message
#'
#' @export

configuration <- function(config) {

  config$rcl <- c(1, Inf, 1, 0, 0.99, NA)
  config$rclmat <- matrix(config$rcl, ncol = 3, byrow = TRUE)

  if (is.null(config$random_seed)) {
    config$random_seed <- sample(1:999999999999, config$number_of_iterations, replace = FALSE)
  }

  if (config$multiple_random_seeds) {
    if (!is.null(config$file_random_seeds)) {
      ## check random seed file
      random_seeds_file_check <- random_seeds_file_checks(config$file_random_seeds)
      if (!random_seeds_file_check$checks_passed) {
        config$failure <- random_seeds_file_check$failed_check
        return(config)
      } else {
        config$random_seeds <- random_seeds_file_check$random_seeds
      }
    } else {
      config$random_seeds <- create_random_seeds(config$number_of_iterations)
    }
  } else {
    config$random_seeds <- create_random_seeds(1)
  }

  if (config$write_outputs %notin% output_list) {
    config$failure <- write_outputs_error
  }

  if (config$write_outputs %in% output_write_list) {
    if (!base::dir.exists(config$output_folder_path)) {
      config$failure <- output_path_error
    }
  }

  if (config$natural_kernel_type %notin% kernel_list) {
    config$failure <- natural_kernel_error
    return(config)
  }

  if (config$anthropogenic_kernel_type %notin% kernel_list) {
    config$failure <- anthropogenic_kernel_error
    return(config)
  }

  # ensures correct model type
  if (config$model_type %in% si_list) {
    config$model_type <- "SEI"
  } else if (config$model_type %in% sei_list) {
    config$model_type <- "SI"
  } else {
    config$failure <- model_type_error
    return(config)
  }

  seasons <- seq(1, 12, 1)
  if (config$season_month_start %in% seasons &&
      config$season_month_end %in% seasons) {
    season_month_start_end <- c()
    season_month_start_end$start_month <- as.integer(config$season_month_start)
    season_month_start_end$end_month <- as.integer(config$season_month_end)
    config$season_month_start_end <- season_month_start_end
  } else {
    config$failure <- season_month_error
  }

  # ensures latent period is correct for type of model selected
  if (config$model_type == "SEI" && config$latency_period <= 0) {
    config$failure <- latency_period_error
    return(config)
  } else if (config$model_type == "SI" && config$latency_period > 0) {
    config$latency_period <- 0
  }

  # ensure correct treatment method
  if (!config$treatment_method %in% treatment_list) {
    config$failure <- treatment_option_error
    return(config)
  }

  # check output and timestep are correct.
  time_check <- time_checks(
    config$end_date, config$start_date,
    config$time_step, config$output_frequency, config$output_frequency_n
  )
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
    return(config)
  }

  # check that network movement is one of the correct options
  if (config$network_movement %notin% network_movement_options) {
    config$failure <- network_movement_error
    return(config)
  }

  # check that weather_type is correct
  if (config$weather_type %notin% weather_type_list) {
    config$failure <- weather_type_error
    return(config)
  }

  # check that initial raster file exists
  if (config$function_name %in% aws_bucket_list) {
    infected_check <- initial_raster_checks(config$infected_file, config$use_s3, config$bucket)
  } else {
    infected_check <- initial_raster_checks(config$infected_file)
  }
  if (infected_check$checks_passed) {
    infected <- infected_check$raster
    infected <- terra::classify(infected, matrix(c(NA, 0), ncol = 2, byrow = TRUE), right = NA)
  } else {
    config$failure <- infected_check$failed_check
    if (config$failure == file_exists_error) {
      config$failure <- detailed_file_exists_error(config$infected_file)
    }
    return(config)
  }

  zero_matrix <- matrix(0, nrow = nrow(infected[[1]]), ncol = ncol(infected[[1]]))
  one_matrix <- matrix(1, nrow = nrow(infected[[1]]), ncol = ncol(infected[[1]]))

  # check that host raster has the same crs, resolution, and extent
  if (config$function_name %in% aws_bucket_list) {
    host_check <- secondary_raster_checks(config$host_file, infected, config$use_s3, config$bucket)
  } else {
    host_check <- secondary_raster_checks(config$host_file, infected)
  }
  if (host_check$checks_passed) {
    host <- host_check$raster
    config$host <- host
  } else {
    config$failure <- host_check$failed_check
    if (config$failure == file_exists_error) {
      config$failure <- detailed_file_exists_error(config$host_file)
    }
    return(config)
  }

  # check that total populations raster has the same crs, resolution, and extent
  if (config$function_name %in% aws_bucket_list) {
    total_populations_check <-
      secondary_raster_checks(config$total_populations_file, infected, config$use_s3, config$bucket)
  } else {
    total_populations_check <- secondary_raster_checks(config$total_populations_file, infected)
  }
  if (total_populations_check$checks_passed) {
    total_populations <- total_populations_check$raster
  } else {
    config$failure <- total_populations_check$failed_check
    if (config$failure == file_exists_error) {
      config$failure <- detailed_file_exists_error(config$total_populations_file)
    }
    return(config)
  }

  # check that soils raster has the same crs, resolutin, and extent.
  if (config$use_soils) {
    config$soil_survival_steps <- ceiling(1 / config$dispersers_to_soils_percentage)
    soil_reservoirs <- list(zero_matrix)
    for (sr in 2:(config$soil_survival_steps)) {
      soil_reservoirs[[sr]] <- zero_matrix
    }
    if (config$start_with_soil_populations) {
      if (config$function_name %in% aws_bucket_list) {
        soils_check <-
          secondary_raster_checks(
            config$soil_starting_pest_file, infected, config$use_s3, config$bucket)
      } else {
        soils_check <- secondary_raster_checks(config$soil_starting_pest_file, infected)
      }
      if (soils_check$checks_passed) {
        soil_pests <- soils_check$raster
        soil_reservoirs[[config$soil_survival_steps]] <- terra::as.matrix(soil_pests, wide = TRUE)
      } else {
        config$failure <- soils_check$failed_check
        if (config$failure == file_exists_error) {
          config$failure <- detailed_file_exists_error(config$soil_starting_pest_file)
        }
        return(config)
      }
    }
    config$soil_reservoirs <- soil_reservoirs
  } else {
    config$soil_reservoirs <- list(matrix(0, 2, 2))
  }
  
  # check that survival_rates raster has the same crs, resolution, and extent
  if (config$use_survival_rates == TRUE) {
    if (config$function_name %in% aws_bucket_list) {
      survival_rate_check <-
        secondary_raster_checks(config$survival_rates_file, infected, config$use_s3, config$bucket)
    } else {
      survival_rate_check <- secondary_raster_checks(config$survival_rates_file, infected)
    }
    if (survival_rate_check$checks_passed) {
      survival_rates_stack <- survival_rate_check$raster
    } else {
      config$failure <- survival_rate_check$failed_check
      if (config$failure == file_exists_error) {
        config$failure <- detailed_file_exists_error(config$survival_rates_file)
      }
      return(config)
    }

    config$survival_rates <- list(terra::as.matrix(survival_rates_stack[[1]], wide = TRUE))
    if (terra::nlyr(survival_rates_stack) > 1) {
      for (i in 2:config$number_of_years) {
        config$survival_rates[[i]] <- terra::as.matrix(survival_rates_stack[[i]], wide = TRUE)
      }
    }
  } else {
    config$survival_rates <- list(matrix(1, 2, 2))
  }

  # check that temperature raster has the same crs, resolution, and extent
  if (config$use_lethal_temperature == TRUE) {
    if (config$function_name %in% aws_bucket_list) {
      temperature_check <-
        secondary_raster_checks(config$temperature_file, infected, config$use_s3, config$bucket)
    } else {
      temperature_check <- secondary_raster_checks(config$temperature_file, infected)
    }
    if (temperature_check$checks_passed) {
      temperature_stack <- temperature_check$raster
    } else {
      config$failure <- temperature_check$failed_check
      if (config$failure == file_exists_error) {
        config$failure <- detailed_file_exists_error(config$temperature_file)
      }
      return(config)
    }

    config$temperature <- list(terra::as.matrix(temperature_stack[[1]], wide = TRUE))
    if (terra::nlyr(temperature_stack) > 1) {
      for (i in 2:config$number_of_years) {
        config$temperature[[i]] <- terra::as.matrix(temperature_stack[[i]], wide = TRUE)
      }
    }
  } else {
    config$temperature <- list(matrix(1, 2, 2))
  }
  
  # check that temp and precip rasters have the same crs, resolution, and extent
  config$weather <- FALSE
  if (config$temp == TRUE) {
    if (config$function_name %in% aws_bucket_list) {
      temperature_coefficient_check <-
        secondary_raster_checks(config$temperature_coefficient_file, infected,
                                config$use_s3, config$bucket)
      if (config$weather_type == "probabilistic") {
        temperature_coefficient_sd_check <-
          secondary_raster_checks(config$temperature_coefficient_sd_file, infected,
                                  config$use_s3, config$bucket)
      }
    } else {
      temperature_coefficient_check <-
        secondary_raster_checks(config$temperature_coefficient_file, infected)
      if (config$weather_type == "probabilistic") {
        temperature_coefficient_sd_check <-
          secondary_raster_checks(config$temperature_coefficient_sd_file, infected)
      }
    }

    if (temperature_coefficient_check$checks_passed) {
      temperature_coefficient <- temperature_coefficient_check$raster
    } else {
      config$failure <- temperature_coefficient_check$failed_check
      if (config$failure == file_exists_error) {
        config$failure <- detailed_file_exists_error(config$temperature_coefficient_file)
      }
      return(config)
    }

    if (config$weather_type == "probabilistic") {
      if (temperature_coefficient_sd_check$checks_passed) {
        temperature_coefficient_sd <- temperature_coefficient_sd_check$raster
      } else {
        config$failure <- temperature_coefficient_sd_check$failed_check
        if (config$failure == file_exists_error) {
          config$failure <- detailed_file_exists_error(config$temperature_coefficient_sd_file)
        }
        return(config)
      }
    }

    config$weather <- TRUE
    weather_coefficient_stack <- temperature_coefficient
    if (config$weather_type == "probabilistic") {
      weather_coefficient_sd_stack <- temperature_coefficient_sd
    }

    if (config$precip == TRUE) {
      if (config$function_name %in% aws_bucket_list) {
        precipitation_coefficient_check <-
          secondary_raster_checks(config$precipitation_coefficient_file, infected,
                                  config$use_s3, config$bucket)
        if (config$weather_type == "probabilistic") {
          precipitation_coefficient_sd_check <-
            secondary_raster_checks(config$precipitation_coefficient_sd_file, infected,
                                    config$use_s3, config$bucket)
        }

      } else {
        precipitation_coefficient_check <-
          secondary_raster_checks(config$precipitation_coefficient_file, infected)
        if (config$weather_type == "probabilistic") {
          precipitation_coefficient_sd_check <-
            secondary_raster_checks(config$precipitation_coefficient_sd_file, infected)
        }
      }

      if (precipitation_coefficient_check$checks_passed) {
        precipitation_coefficient <- precipitation_coefficient_check$raster
      } else {
        config$failure <- precipitation_coefficient_check$failed_check
        if (config$failure == file_exists_error) {
          config$failure <- detailed_file_exists_error(config$precipitation_coefficient_file)
        }
        return(config)
      }

      if (config$weather_type == "probabilistic") {
        if (precipitation_coefficient_sd_check$checks_passed) {
          precipitation_coefficient_sd <- precipitation_coefficient_sd_check$raster
        } else {
          config$failure <- precipitation_coefficient_sd_check$failed_check
          if (config$failure == file_exists_error) {
            config$failure <- detailed_file_exists_error(config$precipitation_coefficient_sd_file)
          }
          return(config)
        }
      }

      weather_coefficient_stack <- weather_coefficient_stack * precipitation_coefficient
      if (config$weather_type == "probabilistic") {
        # compute sd from combined sd of the two rasters hard coded 10 years as our current
        weather_coefficient_sd_stack <-
          combined_sd(temperature_coefficient_sd, precipitation_coefficient_sd,
                      temperature_coefficient, precipitation_coefficient, 10, 10)
      }
    }
  } else if (config$precip == TRUE) {
    if (config$function_name %in% aws_bucket_list) {
      precipitation_coefficient_check <-
        secondary_raster_checks(config$precipitation_coefficient_file, infected,
                                config$use_s3, config$bucket)
      if (config$weather_type == "probabilistic") {
        precipitation_coefficient_sd_check <-
          secondary_raster_checks(config$precipitation_coefficient_sd_file, infected,
                                  config$use_s3, config$bucket)
      }
    } else {
      precipitation_coefficient_check <-
        secondary_raster_checks(config$precipitation_coefficient_file, infected)
      if (config$weather_type == "probabilistic") {
        precipitation_coefficient_sd_check <-
          secondary_raster_checks(config$precipitation_coefficient_sd_file, infected)
      }
    }

    if (precipitation_coefficient_check$checks_passed) {
      precipitation_coefficient <- precipitation_coefficient_check$raster
    } else {
      config$failure <- precipitation_coefficient_check$failed_check
      if (config$failure == file_exists_error) {
        config$failure <- detailed_file_exists_error(config$precipitation_coefficient_file)
      }
      return(config)
    }

    if (config$weather_type == "probabilistic") {
      if (precipitation_coefficient_sd_check$checks_passed) {
        precipitation_coefficient_sd <- precipitation_coefficient_sd_check$raster
      } else {
        config$failure <- precipitation_coefficient_sd_check$failed_check
        if (config$failure == file_exists_error) {
          config$failure <- detailed_file_exists_error(config$precipitation_coefficient_sd_file)
        }
        return(config)
      }
    }

    config$weather <- TRUE

    weather_coefficient_stack <- precipitation_coefficient
    if (config$weather_type == "probabilistic") {
      weather_coefficient_sd_stack <- precipitation_coefficient_sd
    }
  }

  if (config$weather == TRUE) {
    config$weather_size <- terra::nlyr(weather_coefficient_stack)
    if (config$weather_type == "deterministic") {
      if (config$number_of_time_steps > config$weather_size) {
        config$failure <- weather_size_deterministic_error
        return(config)
      }
    }
    
    config$weather_coefficient <- vector("list", config$weather_size)
    for(i in seq_along(config$weather_coefficient)) {
      current_month <- i %% 12
      current_month <- ifelse(current_month == 0, 12, current_month)
      
      if(current_month >= config$season_month_start && current_month
          <= config$season_month_end) {
        config$weather_coefficient[[i]] <- terra::as.matrix(weather_coefficient_stack[[i]],
                                                     wide = TRUE)
        } else {
          config$weather_coefficient[[i]] <- matrix(0, 2, 2)
        }
      }

    if (config$weather_type == "probabilistic") {
      if (config$number_of_time_steps > config$weather_size) {
        config$failure <- weather_size_probabilitic_error
        return(config)
      }

      if (config$weather_size != terra::nlyr(weather_coefficient_sd_stack)) {
        config$failure <- weather_sd_layer_error
        return(config)
      }
      
      config$weather_coefficient_sd <- vector("list", config$weather_size)
      for(i in seq_along(config$weather_coefficient_sd)) {
        current_month <- i %% 12
        current_month <- ifelse(current_month == 0, 12, current_month)
        
        if(current_month >= config$season_month_start && current_month
           <= config$season_month_end) {
          config$weather_coefficient_sd[[i]] <- terra::as.matrix(weather_coefficient_sd_stack[[i]],
                                                       wide = TRUE)
          } else {
            config$weather_coefficient_sd[[i]] <- matrix(0, 2, 2)
          }
        }
      } else {
        config$weather_coefficient_sd <- list(zero_matrix)
        }
    } else {
      config$weather_size <- 1
      config$weather_type <- "None"
      config$weather_coefficient <- list(one_matrix)
      config$weather_coefficient_sd <- list(zero_matrix)
      config$weather_type <- "none"
    }
  
  rm(one_matrix)

  if (config$management == TRUE) {
    if (config$function_name %in% aws_bucket_list) {
      treatments_check <-
        secondary_raster_checks(config$treatments_file, infected, config$use_s3, config$bucket)
    } else {
      treatments_check <- secondary_raster_checks(config$treatments_file, infected)
    }

    if (treatments_check$checks_passed) {
      treatment_stack <- treatments_check$raster
    } else {
      config$failure <- treatments_check$failed_check
      if (config$failure == file_exists_error) {
        config$failure <- detailed_file_exists_error(config$treatments_file)
      }
      return(config)
    }

    treatment_check <- treatment_checks(
      treatment_stack, config$treatments_file,
      config$pesticide_duration,
      config$treatment_dates,
      config$pesticide_efficacy
    )
    if (treatment_check$checks_passed) {
      config$treatment_maps <- treatment_check$treatment_maps
    } else {
      config$failure <- treatment_check$failed_check
      return(config)
    }
  } else {
    config$treatment_maps <- list(zero_matrix)
    config$treatment_dates <- c(config$start_date)
  }

  # setup up movements to be used in the model converts from lat/long to i/j
  if (config$use_movements) {
    movements_check <-
      movement_checks(config$movements_file, infected, config$start_date, config$end_date)
    if (movements_check$checks_passed) {
      config$movements <- movements_check$movements
      config$movements_dates <- movements_check$movements_dates
    } else {
      config$failure <- movements_check$failed_check
      if (config$failure == file_exists_error) {
        config$failure <- detailed_file_exists_error(config$movements_file)
      }
      return(config)
    }
  } else {
    config$movements <- list(0, 0, 0, 0, 0)
    config$movements_dates <- config$start_date
  }

  exposed <- list(zero_matrix)
  config$total_exposed <- zero_matrix

  if (config$model_type == "SEI" && config$latency_period > 1) {
    for (ex in 2:(config$latency_period + 1)) {
      exposed[[ex]] <- zero_matrix
    }
  }

  if (config$model_type == "SEI" && config$start_exposed) {
    if (config$function_name %in% aws_bucket_list) {
      exposed_check <-
        secondary_raster_checks(config$exposed_file, infected, config$use_s3, config$bucket)
    } else {
      exposed_check <- secondary_raster_checks(config$exposed_file, infected)
    }
    if (exposed_check$checks_passed) {
      exposed2 <- exposed_check$raster
      if (config$use_initial_condition_uncertainty) {
        if (terra::nlyr(exposed2) == 2) {
          config$exposed_mean <- terra::as.matrix(exposed2[[1]], wide = TRUE)
          config$exposed_sd <- terra::as.matrix(exposed2[[2]], wide = TRUE)
        } else {
          config$failure <- initial_cond_uncert_error
          return(config)
        }
      } else {
        config$exposed_mean <- terra::as.matrix(exposed2[[1]], wide = TRUE)
        config$exposed_sd <- matrix(0, 2, 2)
      }
    } else {
      config$failure <- exposed_check$failed_check
      if (config$failure == file_exists_error) {
        config$failure <- detailed_file_exists_error(config$exposed_file)
      }
      return(config)
    }
  } else {
    config$exposed_mean <- zero_matrix
    config$exposed_sd <- matrix(0, 2, 2)
  }

  # create spatial indices for computational speed up.
  suitable <- host[[1]] + infected[[1]]
  if (config$use_host_uncertainty && terra::nlyr(host) > 1) {
    suitable <- suitable + host[[2]]
  }
  if (config$use_initial_condition_uncertainty && terra::nlyr(infected) > 1) {
    suitable <- suitable + infected[[2]]
  }
  if (config$model_type == "SEI" && config$start_exposed) {
    suitable <- suitable + exposed2[[1]]
    if (config$use_initial_condition_uncertainty && terra::nlyr(exposed2) > 1) {
      suitable <- suitable + exposed2[[2]]
    }
  }
  
  suitable <- terra::classify(suitable, matrix(c(-Inf, 0, NA), ncol = 3, byrow = TRUE),
                              right = FALSE)
  suitable_cells <- which(!is.na(values(suitable)))
  suitable_row <- as.integer(terra::rowFromCell(suitable, suitable_cells) - 1L)
  suitable_col <- as.integer(terra::colFromCell(suitable, suitable_cells) - 1L)
  config$spatial_indices <- lapply(seq_along(suitable_row),
                                   function(i) c(suitable_row[i],
                                                 suitable_col[i]))
  names(config$spatial_indices) <- NULL
  rm(suitable, suitable_cells, suitable_row, suitable_col)

  res <- c()
  res$ew_res <- terra::xres(infected)
  res$ns_res <- terra::yres(infected)
  config$res <- res
  rows_cols <- c()
  rows_cols$num_rows <- terra::nrow(infected)
  rows_cols$num_cols <- terra::ncol(infected)
  config$rows_cols <- rows_cols

  if (config$use_host_uncertainty) {
    if (terra::nlyr(host) == 2) {
      config$host_mean <- terra::as.matrix(host[[1]], wide = TRUE)
      config$host_sd <- terra::as.matrix(host[[2]], wide = TRUE)
    } else {
      config$failure <- host_uncert_error
      return(config)
    }
  } else {
    config$host_mean <- terra::as.matrix(host[[1]], wide = TRUE)
    config$host_sd <- matrix(0, 2, 2)
  }

  if (!is.null(config$mask)) {
    if (config$function_name %in% aws_bucket_list) {
      mask_check <- secondary_raster_checks(config$mask, infected, config$use_s3, config$bucket)
    } else {
      mask_check <- secondary_raster_checks(config$mask, infected)
    }
    if (mask_check$checks_passed) {
      mask <- mask_check$raster
      mask <- terra::classify(mask, config$rclmat)
      host_mask <- terra::classify(host[[1]], config$rclmat)
      if (config$use_host_uncertainty && terra::nlyr(host) > 1) {
        host_mask2 <- terra::classify(host[[2]], config$rclmat)
        host_mask <- terra::mask(host_mask2, host_mask, maskvalues = NA, updatevalue = NA)
      }
      mask <- terra::mask(host_mask, mask, maskvalues = NA, updatevalue = NA)
      config$mask <- mask
      config$mask_matrix <- terra::as.matrix(mask, wide = TRUE)
    } else {
      config$failure <- mask_check$failed_check
      if (config$failure == file_exists_error) {
        config$failure <- detailed_file_exists_error(config$mask)
      }
      return(config)
    }
  } else {
    mask <- terra::classify(host[[1]], config$rclmat)
    config$mask <- mask
    config$mask_matrix <- terra::as.matrix(mask, wide = TRUE)
  }

  if (config$use_initial_condition_uncertainty) {
    if (terra::nlyr(infected) == 2) {
      config$infected_mean <- terra::as.matrix(infected[[1]], wide = TRUE)
      config$infected_sd <- terra::as.matrix(infected[[2]], wide = TRUE)
    } else {
      config$failure <- initial_cond_uncert_error
      return(config)
    }
  } else {
    config$infected_mean <- terra::as.matrix(infected[[1]], wide = TRUE)
    config$infected_sd <- matrix(0, 2, 2)
  }

  exposed[[config$latency_period + 1]] <- config$exposed_mean
  config$total_exposed <- config$exposed_mean
  config$exposed <- exposed

  susceptible_mean <- config$host_mean - config$infected_mean - config$exposed_mean
  susceptible_mean[susceptible_mean < 0] <- 0
  config$susceptible_mean <- terra::as.matrix(susceptible_mean, wide = TRUE)

  config$total_populations <- terra::as.matrix(total_populations, wide = TRUE)
  config$mortality <- zero_matrix
  config$resistant <- zero_matrix

  # check that quarantine raster has the same crs, resolution, and extent
  if (config$use_quarantine) {
    if (config$function_name %in% aws_bucket_list) {
      quarantine_check <-
        secondary_raster_checks(config$quarantine_areas_file, host, config$use_s3, config$bucket)
    } else {
      quarantine_check <- secondary_raster_checks(config$quarantine_areas_file, host)
    }

    if (quarantine_check$checks_passed) {
      quarantine_areas <- quarantine_check$raster
      config$quarantine_areas <- terra::as.matrix(quarantine_areas, wide = TRUE)
    } else {
      config$failure <- quarantine_check$failed_check
      if (config$failure == file_exists_error) {
        config$failure <- detailed_file_exists_error(config$quarantine_areas_file)
      }
      return(config)
    }
  } else {
    # set quarantine areas to all zeros. meaning no quarantine areas are considered
    config$quarantine_areas <- matrix(0, 2, 2)
  }

  mortality_tracker <- zero_matrix
  rm(zero_matrix)
  mortality_tracker2 <- list(mortality_tracker)
  if (config$mortality_on) {
    mortality_length <- 1 / config$mortality_rate + config$mortality_time_lag
    for (mt in 2:(mortality_length)) {
      mortality_tracker2[[mt]] <- mortality_tracker
    }
  }
  # add currently infected cells to last element of the mortality tracker so
  # that mortality occurs at the appropriate interval
  if (config$mortality_on) {
    mortality_tracker2[[length(mortality_tracker2)]] <- config$infected_mean
  }


  config$mortality_tracker <- mortality_tracker2

  if (config$function_name %in% parallel_function_list) {
    if (is.na(config$number_of_cores) ||
      config$number_of_cores > parallel::detectCores()) {
      core_count <- parallel::detectCores() - 1
    } else {
      core_count <- config$number_of_cores
    }
    config$core_count <- core_count
  }

  if (config$function_name %in% parameter_draw_list) {

    if (nrow(config$parameter_cov_matrix) != 8 ||
      ncol(config$parameter_cov_matrix) != 8) {
      config$failure <- covariance_mat_error
      return(config)
    }

    if (length(config$parameter_means) != 8) {
      config$failure <- paramter_means_error
      return(config)
    }

    if (config$anthropogenic_kernel_type != "network") {
      config$parameter_means[7] <- config$res$ew_res / 2
    }

    if (config$anthropogenic_kernel_type != "network") {
      config$parameter_means[8] <-
        min(config$rows_cols$num_cols, config$rows_cols$num_rows) * config$res$ew_res
    }

    if (config$parameter_means[7] < config$res$ew_res / 2) {
      config$failure <- network_min_distance_small_error
      return(config)
    }

    if (config$parameter_means[7] > config$parameter_means[8]) {
      config$failure <- network_min_distance_large_error
      return(config)
    }

    if (config$parameter_means[8] >
        (min(config$rows_cols$num_cols, config$rows_cols$num_rows) * config$res$ew_res)) {
      config$failure <- network_max_distance_large_error
      return(config)
    }

    config <- draw_parameters(config)

    if (any(config$percent_natural_dispersal < 1.0)) {
      config$use_anthropogenic_kernel <- TRUE
    } else {
      config$use_anthropogenic_kernel <- FALSE
    }
  }

  if (config$function_name %in% val_cal_list) {
    config$use_anthropogenic_kernel <- TRUE
    # Load observed data on occurrence
    config$infection_years <- terra::rast(config$infected_years_file)
    config$infection_years[] <- as.integer(config$infection_years[])
    config$num_layers_infected_years <- terra::nlyr(config$infection_years)

    if (config$num_layers_infected_years < config$number_of_outputs) {
      config$failure <-
        infection_years_length_error(config$num_layers_infected_years, config$number_of_time_steps)
      return(config)
    }

    config$infection_years2 <- list(terra::as.matrix(config$infection_years[[1]], wide = TRUE))
    if (terra::nlyr(config$infection_years) > 1) {
      for (i in 2:terra::nlyr(config$infection_years)) {
        config$infection_years2[[i]] <- terra::as.matrix(config$infection_years[[i]], wide = TRUE)
      }
    }
  }

  if (config$function_name %in% c("calibrate") &&
    config$calibration_method == "ABC") {
    config$num_particles <- config$number_of_generations * config$generation_size
    config$total_particles <- 1
    config$current_particles <- 1
    config$proposed_particles <- 1
    config$current_bin <- 1
  }

  if (config$function_name == "auto-manage") {
    ## management module information
    config$num_cells <-
      round((config$budget / config$cost_per_meter_sq) / (config$ew_res * config$ns_res))
    config$buffer_cells <- config$buffer / config$ew_res
    config$years_simulated <- length(config$years)
  }

  # add / to output folder path if not provided by user.
  if (substr(config$output_folder_path, nchar(config$output_folder_path),
             nchar(config$output_folder_path)) == "/") {
    config$output_folder_path <- config$output_folder_path
  } else {
    config$output_folder_path <- paste(config$output_folder_path, "/", sep = "")
  }

  config$crs <- terra::crs(config$host)
  config$xmax <- terra::xmax(config$host)
  config$xmin <- terra::xmin(config$host)
  config$ymax <- terra::ymax(config$host)
  config$ymin <- terra::ymin(config$host)
  bounding_box <- c()
  bounding_box$north <- config$ymax
  bounding_box$south <- config$ymin
  bounding_box$west <- config$xmin
  bounding_box$east <- config$xmax
  config$bounding_box <- bounding_box
  gc()
  return(config)
}
