#' PoPS (configuration
#'
#' Function with a single input and output list for parsing, transforming,
#' and performing all checks for all functions to run the pops c++ model
#'
#' @param config list of all data necessary used to set up c++ model
#'
#' @importFrom raster
#' cellStats  calc extract
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
  "%notin%" <- Negate("%in%")

  # Check that all data has same length if using multiple species currently
  # only implemented for auto manage
  if (config$function_name == "auto-manage") {
    multispecies_check <-
      multispecies_checks(
        config$species,
        config$infected_files,
        config$parameter_means,
        config$parameter_cov_matrix,
        config$natural_kernel_type,
        config$anthropogenic_kernel_type,
        config$natural_dir,
        config$anthropogenic_dir,
        config$model_type,
        config$host_file,
        config$total_populations_file,
        config$temp,
        config$temperature_coefficient_file,
        config$precip,
        config$precipitation_coefficient_file,
        config$latency_period,
        config$time_step,
        config$season_month_start,
        config$season_month_end,
        config$use_lethal_temperature,
        config$temperature_file,
        config$lethal_temperature,
        config$lethal_temperature_month,
        config$mortality_on,
        config$mortality_rate,
        config$mortality_time_lag,
        config$movements_file,
        config$use_movements,
        config$start_exposed,
        config$quarantine_areas_file,
        config$use_quarantine,
        config$use_spreadrates
      )
    if (!multispecies_check$checks_passed) {
      config$failure <- multispecies_check$failed_check
      return(config)
    }
  }

  if (config$function_name == "sensitivity") {
    config$sensitivity_rcl <- c(0.10, Inf, 1, 0, 0.10, 0)
    config$sensitivity_rclmat <- matrix(config$rcl, ncol = 3, byrow = TRUE)
  }

  config$rcl <- c(1, Inf, 1, 0, 0.99, NA)
  config$rclmat <- matrix(config$rcl, ncol = 3, byrow = TRUE)

  config$output_list <- c("all_simulations", "summary_outputs", "None")
  config$output_write_list <- c("all_simulations", "summary_outputs")

  if (config$write_outputs %notin% config$output_list) {
    config$failure <- write_outputs_error
  }

  if (config$write_outputs %in% config$output_write_list) {
    if (!base::dir.exists(config$output_folder_path)) {
      config$failure <- output_path_error
    }
  }

  # ensures correct model type
  if (config$model_type %in%
    c(
      "SEI",
      "susceptible-exposed-infected",
      "susceptible_exposed_infected",
      "Susceptible-Exposed-Infected",
      "Susceptible_Exposed_Infected"
    )) {
    config$model_type <- "SEI"
  } else if (config$model_type %in%
    c(
      "SI",
      "susceptible-infected",
      "susceptible_infected",
      "Susceptible-Infected",
      "Susceptible_Infected"
    )) {
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
  if (!config$treatment_method %in% c("ratio", "all infected")) {
    config$failure <- treatment_option_error
    return(config)
  }

  if (is.null(config$random_seed)) {
    config$random_seed <- round(stats::runif(1, 1, 1000000))
  }

  ## check output and timestep are correct.
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

  # check that initial raster file exists
  if (config$function_name %in% c("casestudy_creation", "model_api")) {
    infected_check <-
      initial_raster_checks(config$infected_file, config$use_s3, config$bucket)
  } else {
    infected_check <- initial_raster_checks(config$infected_file)
  }
  if (infected_check$checks_passed) {
    infected <- infected_check$raster
    if (terra::nlyr(infected) > 1) {
      infected <- output_from_raster_mean_and_sd(infected)
    }
  } else {
    config$failure <- infected_check$failed_check
    return(config)
  }

  # check that host raster has the same crs, resolution, and extent
  if (config$function_name %in% c("casestudy_creation", "model_api")) {
    host_check <-
      secondary_raster_checks(
        config$host_file, infected,
        config$use_s3, config$bucket
      )
  } else {
    host_check <- secondary_raster_checks(config$host_file, infected)
  }
  if (host_check$checks_passed) {
    host <- host_check$raster
    if (terra::nlyr(host) > 1) {
      host <- output_from_raster_mean_and_sd(host)
    }
    config$host <- host
  } else {
    config$failure <- host_check$failed_check
    return(config)
  }

  if (!is.null(config$mask)) {
    if (config$function_name %in% c("casestudy_creation", "model_api")) {
      mask_check <-
        secondary_raster_checks(
          config$mask, infected,
          config$use_s3, config$bucket
        )
    } else {
      mask_check <- secondary_raster_checks(config$mask, infected)
    }
    if (mask_check$checks_passed) {
      mask <- mask_check$raster
      mask <- terra::classify(mask, config$rclmat)
      host_mask <- terra::classify(host, config$rclmat)
      mask <- terra::mask(host_mask, mask, maskvalues = NA, updatevalue = NA)
      config$mask <- mask
      config$mask_matrix <- terra::as.matrix(mask, wide = TRUE)
    } else {
      config$failure <- mask_check$failed_check
      return(config)
    }
  } else {
    mask <- terra::classify(host, config$rclmat)
    config$mask <- mask
    config$mask_matrix <- terra::as.matrix(mask, wide = TRUE)
  }

  suitable <- host + infected
  suitable_points <- terra::as.points(suitable)
  names(suitable_points) <- "data"
  suitable_points <- suitable_points[suitable_points$data > 0]
  suitable_cells <-
    terra::extract(suitable, suitable_points, cells = TRUE)$cell
  suitable_row <- terra::rowFromCell(suitable, suitable_cells)
  suitable_row <- suitable_row - 1
  suitable_row <- as.integer(suitable_row)
  suitable_col <- terra::colFromCell(suitable, suitable_cells)
  suitable_col <- suitable_col - 1
  suitable_col <- as.integer(suitable_col)
  spatial_indices2 <- data.frame(row = suitable_row, col = suitable_col)
  spatial_indices2 <- unname(spatial_indices2)
  spatial_indices2 <- as.matrix(spatial_indices2)
  spatial_indices <- list()
  # movements_date
  for (i in seq_len(terra::nrow(spatial_indices2))) {
    spatial_indices[[i]] <- spatial_indices2[i, 1:2]
  }

  spatial_indices <- unname(spatial_indices)
  config$spatial_indices <- spatial_indices

  # check that total populations raster has the same crs, resolution, and extent
  if (config$function_name %in% c("casestudy_creation", "model_api")) {
    total_populations_check <-
      secondary_raster_checks(
        config$total_populations_file, infected,
        config$use_s3, config$bucket
      )
  } else {
    total_populations_check <-
      secondary_raster_checks(config$total_populations_file, infected)
  }
  if (total_populations_check$checks_passed) {
    total_populations <- total_populations_check$raster
    if (terra::nlyr(total_populations) > 1) {
      total_populations <- output_from_raster_mean_and_sd(total_populations)
    }
  } else {
    config$failure <- total_populations_check$failed_check
    return(config)
  }

  susceptible <- host - infected
  susceptible[susceptible < 0] <- 0

  # check that temperature raster has the same crs, resolution, and extent
  if (config$use_lethal_temperature == TRUE) {
    if (config$function_name %in% c("casestudy_creation", "model_api")) {
      temperature_check <-
        secondary_raster_checks(
          config$temperature_file, infected,
          config$use_s3, config$bucket
        )
    } else {
      temperature_check <-
        secondary_raster_checks(config$temperature_file, infected)
    }
    if (temperature_check$checks_passed) {
      temperature_stack <- temperature_check$raster
    } else {
      config$failure <- temperature_check$failed_check
      return(config)
    }

    temperature <- list(terra::as.matrix(temperature_stack[[1]],
      wide = TRUE
    ))
    if (terra::nlyr(temperature_stack) > 1) {
      for (i in 2:config$number_of_years) {
        temperature[[i]] <- terra::as.matrix(temperature_stack[[i]],
          wide = TRUE
        )
      }
    }
  } else {
    temperature <- host
    terra::values(temperature) <- 1
    temperature <- list(terra::as.matrix(temperature,
      wide = TRUE
    ))
  }

  config$temperature <- temperature

  # check that temp and precip rasters have the same crs, resolution, and extent
  config$weather <- FALSE
  if (config$temp == TRUE) {
    if (config$function_name %in% c("casestudy_creation", "model_api")) {
      temperature_coefficient_check <-
        secondary_raster_checks(
          config$temperature_coefficient_file, infected,
          config$use_s3, config$bucket
        )
    } else {
      temperature_coefficient_check <-
        secondary_raster_checks(config$temperature_coefficient_file, infected)
    }
    if (temperature_coefficient_check$checks_passed) {
      temperature_coefficient <- temperature_coefficient_check$raster
    } else {
      config$failure <- temperature_coefficient_check$failed_check
      return(config)
    }

    config$weather <- TRUE
    weather_coefficient_stack <- temperature_coefficient
    if (config$precip == TRUE) {
      if (config$function_name %in% c("casestudy_creation", "model_api")) {
        precipitation_coefficient_check <-
          secondary_raster_checks(
            config$precipitation_coefficient_file,
            infected, config$use_s3, config$bucket
          )
      } else {
        precipitation_coefficient_check <-
          secondary_raster_checks(
            config$precipitation_coefficient_file,
            infected
          )
      }
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
    if (config$function_name %in% c("casestudy_creation", "model_api")) {
      precipitation_coefficient_check <-
        secondary_raster_checks(
          config$precipitation_coefficient_file, infected,
          config$use_s3, config$bucket
        )
    } else {
      precipitation_coefficient_check <-
        secondary_raster_checks(config$precipitation_coefficient_file, infected)
    }
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
    weather_coefficient <- list(terra::as.matrix(
      weather_coefficient_stack[[1]],
      wide = TRUE
    ))
    for (i in 2:terra::nlyr(weather_coefficient_stack)) {
      weather_coefficient[[i]] <- terra::as.matrix(
        weather_coefficient_stack[[i]],
        wide = TRUE
      )
    }
  } else {
    weather_coefficient <- host
    terra::values(weather_coefficient) <- 1
    weather_coefficient <- list(terra::as.matrix(weather_coefficient,
      wide = TRUE
    ))
  }

  config$weather_coefficient <- weather_coefficient

  if (config$management == TRUE) {
    if (config$function_name %in% c("casestudy_creation", "model_api")) {
      treatments_check <-
        secondary_raster_checks(
          config$treatments_file, infected,
          config$use_s3, config$bucket
        )
    } else {
      treatments_check <-
        secondary_raster_checks(config$treatments_file, infected)
    }

    if (treatments_check$checks_passed) {
      treatment_stack <- treatments_check$raster
    } else {
      config$failure <- treatments_check$failed_check
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
    treatment_map <- host
    treatment_map[] <- 0
    config$treatment_maps <- list(terra::as.matrix(treatment_map,
      wide = TRUE
    ))
    config$treatment_dates <- c(config$start_date)
  }

  res <- c()
  res$ew_res <- terra::xres(susceptible)
  res$ns_res <- terra::yres(susceptible)
  config$res <- res
  rows_cols <- c()
  rows_cols$num_rows <- terra::nrow(susceptible)
  rows_cols$num_cols <- terra::ncol(susceptible)
  config$rows_cols <- rows_cols
  # setup up movements to be used in the model converts from lat/long to i/j
  if (config$use_movements) {
    movements_check <- movement_checks(
      config$movements_file, infected,
      config$start_date, config$end_date
    )
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
  terra::values(mortality_tracker) <- 0
  mortality_tracker <- terra::as.matrix(mortality_tracker,
    wide = TRUE
  )
  exposed <- list(mortality_tracker)
  config$total_exposed <- mortality_tracker

  if (config$model_type == "SEI" & config$latency_period > 1) {
    for (ex in 2:(config$latency_period + 1)) {
      exposed[[ex]] <- mortality_tracker
    }
  }

  if (config$model_type == "SEI" & config$start_exposed) {
    if (config$function_name %in% c("casestudy_creation", "model_api")) {
      exposed_check <-
        secondary_raster_checks(
          config$exposed_file, infected,
          config$use_s3, config$bucket
        )
    } else {
      exposed_check <- secondary_raster_checks(config$exposed_file, infected)
    }
    if (exposed_check$checks_passed) {
      exposed2 <- exposed_check$raster
      susceptible <- susceptible - exposed2
      susceptible[susceptible < 0] <- 0
      exposed[[config$latency_period + 1]] <-
        terra::as.matrix(exposed2, wide = TRUE)
      config$total_exposed <- terra::as.matrix(exposed2, wide = TRUE)
    } else {
      config$failure <- exposed_check$failed_check
      return(config)
    }
  }

  kernel_list <- c(
    "cauchy",
    "Cauchy",
    "exponential",
    "Exponential",
    "uniform",
    "Uniform",
    "deterministic neighbor",
    "deterministic-neighbor",
    "power law",
    "power-law",
    "Power-law",
    "Power-Law",
    "Power Law",
    "Power law",
    "hyperbolic secant",
    "hyperbolic-secant",
    "Hyperbolic-secant",
    "Hyperbolic-Secant",
    "Hyperbolic secant",
    "Hyperbolic Secant",
    "gamma",
    "Gamma",
    # "exponential power",
    # "exponential-power",
    # "Exponential-power",
    # "Exponential-Power",
    # "Exponential power",
    "weibull",
    "Weibull",
    "normal",
    # "log normal",
    # "log-normal",
    # "Log-normal",
    # "Log-Normal",
    # "Log normal",
    # "Log Normal",
    "logistic",
    "Logistic"
  )

  if (config$natural_kernel_type %notin% kernel_list) {
    config$failure <- natural_kernel_error
    return(config)
  }

  if (config$anthropogenic_kernel_type %notin% kernel_list) {
    config$failure <- anthropogenic_kernel_error
    return(config)
  }

  infected <- terra::as.matrix(infected,
    wide = TRUE
  )
  config$susceptible <- terra::as.matrix(susceptible,
    wide = TRUE
  )
  config$total_populations <- terra::as.matrix(total_populations,
    wide = TRUE
  )
  config$total_hosts <- terra::as.matrix(host, wide = TRUE)

  config$mortality <- mortality_tracker
  config$resistant <- mortality_tracker

  # check that quarantine raster has the same crs, resolution, and extent
  if (config$use_quarantine) {
    if (config$function_name %in% c("casestudy_creation", "model_api")) {
      quarantine_check <-
        secondary_raster_checks(
          config$quarantine_areas_file, host,
          config$use_s3, config$bucket
        )
    } else {
      quarantine_check <-
        secondary_raster_checks(config$quarantine_areas_file, host)
    }

    if (quarantine_check$checks_passed) {
      quarantine_areas <- quarantine_check$raster
      config$quarantine_areas <- terra::as.matrix(quarantine_areas, wide = TRUE)
    } else {
      config$failure <- quarantine_check$failed_check
      return(config)
    }
  } else {
    # set quarantine areas to all zeros. meaning no quarantine areas are
    # considered
    config$quarantine_areas <- mortality_tracker
  }

  mortality_tracker2 <- list(mortality_tracker)
  if (config$mortality_on) {
    mortality_length <- 1 / config$mortality_rate + config$mortality_time_lag

    for (mt in 2:(mortality_length)) {
      mortality_tracker2[[mt]] <- mortality_tracker
    }
  }
  ## add currently infected cells to last element of the mortality tracker so
  ## that mortality occurs at the appropriate interval
  mortality_tracker2[[length(mortality_tracker2)]] <- infected

  config$mortality_tracker <- mortality_tracker2
  config$exposed <- exposed
  config$infected <- infected

  if (config$function_name %in% c("validate", "multirun", "sensitivity")) {
    if (is.na(config$number_of_cores) ||
      config$number_of_cores > parallel::detectCores()) {
      core_count <- parallel::detectCores() - 1
    } else {
      core_count <- config$number_of_cores
    }
    config$core_count <- core_count
  }

  if (config$function_name %in%
    c("validate", "pops", "multirun", "sensitivity", "casestudy_creation")) {
    if (nrow(config$parameter_cov_matrix) != 6 |
      ncol(config$parameter_cov_matrix) != 6) {
      config$failure <- covariance_mat_error
      return(config)
    }

    if (length(config$parameter_means) != 6) {
      config$failure <- paramter_means_error
      return(config)
    }

    parameters <- MASS::mvrnorm(
      config$number_of_iterations,
      config$parameter_means,
      config$parameter_cov_matrix
    )
    while (any(parameters[, 1] < 0 |
      parameters[, 2] <= 0 |
      parameters[, 3] > 1 |
      parameters[, 3] <= 0 |
      parameters[, 4] <= 0 |
      parameters[, 5] < 0 |
      parameters[, 6] < 0)) {
      config$number_of_draws <- nrow(parameters[parameters[, 1] < 0 |
        parameters[, 2] <= 0 |
        parameters[, 3] > 1 |
        parameters[, 3] <= 0 |
        parameters[, 4] <= 0 |
        parameters[, 5] < 0 |
        parameters[, 6] < 0, ])
      if (is.null(config$number_of_draws)) {
        config$number_of_draws <- 1
      }

      parameters[parameters[, 1] < 0 |
        parameters[, 2] <= 0 |
        parameters[, 3] > 1 |
        parameters[, 3] <= 0 |
        parameters[, 4] <= 0 |
        parameters[, 5] < 0 |
        parameters[, 6] < 0, ] <-
        MASS::mvrnorm(
          config$number_of_draws,
          config$parameter_means,
          config$parameter_cov_matrix
        )
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
    infection_years <- terra::rast(config$infected_years_file)
    infection_years[] <- as.integer(infection_years[])
    config$num_layers_infected_years <- terra::nlyr(infection_years)

    if (config$num_layers_infected_years < config$number_of_outputs) {
      config$failure <-
        paste("The infection years file must have enough layers to match the
            number of outputs from the model. The number of layers of your
            infected year file is", config$num_layers_infected_years, "and the
            number of outputs is", config$number_of_time_steps, sep = " ")
      return(config)
    }

    # add / to output folder path if not provided by user.
    if (substr(config$output_folder_path, nchar(config$output_folder_path),
               nchar(config$output_folder_path)) == "/") {
      config$output_folder_path <- config$output_folder_path
    } else {
      config$output_folder_path <- paste(config$output_folder_path, "/", sep = "")
    }

    infection_years2 <- list(terra::as.matrix(infection_years[[1]], wide = TRUE))
    if (terra::nlyr(infection_years) > 1) {
      for (i in 2:terra::nlyr(infection_years)) {
        infection_years2[[i]] <- terra::as.matrix(infection_years[[i]], wide = TRUE)
      }
    }
    config$infection_years <- infection_years
    config$infection_years2 <- infection_years2
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

  config$crs <- terra::crs(config$host)
  config$xmax <- terra::xmax(config$host)
  config$xmin <- terra::xmin(config$host)
  config$ymax <- terra::ymax(config$host)
  config$ymin <- terra::ymin(config$host)

  return(config)
}
