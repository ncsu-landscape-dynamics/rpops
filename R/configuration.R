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

  # Check for correct kernel options
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
    "Logistic",
    "network",
    "Network"
  )

  if (config$natural_kernel_type %notin% kernel_list) {
    config$failure <- natural_kernel_error
    return(config)
  }

  if (config$anthropogenic_kernel_type %notin% kernel_list) {
    config$failure <- anthropogenic_kernel_error
    return(config)
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
  network_movement_options <- c("walk", "jump", "teleport")
  if (config$network_movement %notin% network_movement_options) {
    config$failure <- network_movement_error
  }

  # check that initial raster file exists
  if (config$function_name %in% c("casestudy_creation", "model_api")) {
    infected_check <- initial_raster_checks(config$infected_file, config$use_s3, config$bucket)
  } else {
    infected_check <- initial_raster_checks(config$infected_file)
  }
  if (infected_check$checks_passed) {
    infected <- infected_check$raster
    infected <- terra::classify(infected, matrix(c(NA, 0), ncol = 2, byrow = TRUE), right = NA)
  } else {
    config$failure <- infected_check$failed_check
    return(config)
  }

  zero_matrix <- infected
  terra::values(zero_matrix) <- 0
  zero_matrix <- terra::as.matrix(zero_matrix, wide = TRUE)

  one_matrix <- infected
  terra::values(one_matrix) <- 0
  one_matrix <- terra::as.matrix(one_matrix, wide = TRUE)

  # check that host raster has the same crs, resolution, and extent
  if (config$function_name %in% c("casestudy_creation", "model_api")) {
    host_check <- secondary_raster_checks(config$host_file, infected, config$use_s3, config$bucket)
  } else {
    host_check <- secondary_raster_checks(config$host_file, infected)
  }
  if (host_check$checks_passed) {
    host <- host_check$raster
    config$host <- host
  } else {
    config$failure <- host_check$failed_check
    return(config)
  }

  # check that total populations raster has the same crs, resolution, and extent
  if (config$function_name %in% c("casestudy_creation", "model_api")) {
    total_populations_check <-
      secondary_raster_checks(config$total_populations_file, infected, config$use_s3, config$bucket)
  } else {
    total_populations_check <- secondary_raster_checks(config$total_populations_file, infected)
  }
  if (total_populations_check$checks_passed) {
    total_populations <- total_populations_check$raster
  } else {
    config$failure <- total_populations_check$failed_check
    return(config)
  }

  # check that survival_rates raster has the same crs, resolution, and extent
  if (config$use_survival_rates == TRUE) {
    if (config$function_name %in% c("casestudy_creation", "model_api")) {
      survival_rate_check <-
        secondary_raster_checks(config$survival_rates_file, infected, config$use_s3, config$bucket)
    } else {
      survival_rate_check <- secondary_raster_checks(config$survival_rates_file, infected)
    }
    if (survival_rate_check$checks_passed) {
      survival_rates_stack <- survival_rate_check$raster
    } else {
      config$failure <- survival_rate_check$failed_check
      return(config)
    }

    survival_rates <- list(terra::as.matrix(survival_rates_stack[[1]], wide = TRUE))
    if (terra::nlyr(survival_rates_stack) > 1) {
      for (i in 2:config$number_of_years) {
        survival_rates[[i]] <- terra::as.matrix(survival_rates_stack[[i]], wide = TRUE)
      }
    }
  } else {
    survival_rates <- list(one_matrix)
  }

  config$survival_rates <- survival_rates

  # check that temperature raster has the same crs, resolution, and extent
  if (config$use_lethal_temperature == TRUE) {
    if (config$function_name %in% c("casestudy_creation", "model_api")) {
      temperature_check <-
        secondary_raster_checks(config$temperature_file, infected, config$use_s3, config$bucket)
    } else {
      temperature_check <- secondary_raster_checks(config$temperature_file, infected)
    }
    if (temperature_check$checks_passed) {
      temperature_stack <- temperature_check$raster
    } else {
      config$failure <- temperature_check$failed_check
      return(config)
    }

    temperature <- list(terra::as.matrix(temperature_stack[[1]], wide = TRUE))
    if (terra::nlyr(temperature_stack) > 1) {
      for (i in 2:config$number_of_years) {
        temperature[[i]] <- terra::as.matrix(temperature_stack[[i]], wide = TRUE)
      }
    }
  } else {
    temperature <- list(one_matrix)
  }

  config$temperature <- temperature

  # check that temp and precip rasters have the same crs, resolution, and extent
  config$weather <- FALSE
  if (config$temp == TRUE) {
    if (config$function_name %in% c("casestudy_creation", "model_api")) {
      temperature_coefficient_check <-
        secondary_raster_checks(config$temperature_coefficient_file, infected,
                                config$use_s3, config$bucket)
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
          secondary_raster_checks(config$precipitation_coefficient_file, infected,
                                  config$use_s3, config$bucket)
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

      weather_coefficient_stack <- weather_coefficient_stack * precipitation_coefficient
    }
  } else if (config$precip == TRUE) {
    if (config$function_name %in% c("casestudy_creation", "model_api")) {
      precipitation_coefficient_check <-
        secondary_raster_checks(config$precipitation_coefficient_file, infected,
                                config$use_s3, config$bucket)
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
    weather_coefficient <- list(terra::as.matrix(weather_coefficient_stack[[1]], wide = TRUE))
    for (i in 2:terra::nlyr(weather_coefficient_stack)) {
      weather_coefficient[[i]] <- terra::as.matrix(weather_coefficient_stack[[i]], wide = TRUE)
    }
  } else {
    weather_coefficient <- list(one_matrix)
  }

  config$weather_coefficient <- weather_coefficient

  if (config$management == TRUE) {
    if (config$function_name %in% c("casestudy_creation", "model_api")) {
      treatments_check <-
        secondary_raster_checks(config$treatments_file, infected, config$use_s3, config$bucket)
    } else {
      treatments_check <- secondary_raster_checks(config$treatments_file, infected)
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
      return(config)
    }
  } else {
    config$movements <- list(0, 0, 0, 0, 0)
    config$movements_dates <- config$start_date
  }

  exposed <- list(zero_matrix)
  config$total_exposed <- zero_matrix

  if (config$model_type == "SEI" & config$latency_period > 1) {
    for (ex in 2:(config$latency_period + 1)) {
      exposed[[ex]] <- zero_matrix
    }
  }

  if (config$model_type == "SEI" & config$start_exposed) {
    if (config$function_name %in% c("casestudy_creation", "model_api")) {
      exposed_check <-
        secondary_raster_checks(config$exposed_file, infected, config$use_s3, config$bucket)
    } else {
      exposed_check <- secondary_raster_checks(config$exposed_file, infected)
    }
    if (exposed_check$checks_passed) {
      exposed2 <- exposed_check$raster
      if (config$use_initial_condition_uncertainty) {
        if (terra::nlyr(exposed2) == 2) {
          exposed_mean <- terra::as.matrix(exposed2[[1]], wide = TRUE)
          exposed_sd <- terra::as.matrix(exposed2[[2]], wide = TRUE)
        } else {
          config$failure <- initial_cond_uncert_error
          return(config)
        }
      } else {
        exposed_mean <- terra::as.matrix(exposed2, wide = TRUE)
        exposed_sd <- zero_matrix
      }
    } else {
      config$failure <- exposed_check$failed_check
      return(config)
    }
  } else {
    exposed_mean <- zero_matrix
    exposed_sd <- zero_matrix
  }

  config$exposed_mean <- exposed_mean
  config$exposed_sd <- exposed_sd

  # create spatial indices for computational speed up.
  suitable <- host[[1]] + infected[[1]]
  if (config$use_host_uncertainty && terra::nlyr(host) > 1) {
    suitable <- suitable + host[[2]]
  }
  if (config$use_initial_condition_uncertainty && terra::nlyr(infected) > 1) {
    suitable <- suitable + infected[[2]]
  }
  if (config$model_type == "SEI" & config$start_exposed) {
    suitable <- suitable + exposed2[[1]]
    if (config$use_initial_condition_uncertainty && terra::nlyr(exposed2) > 1) {
      suitable <- suitable + exposed2[[2]]
    }
  }
  suitable_points <- terra::as.points(suitable)
  names(suitable_points) <- "data"
  suitable_points <- suitable_points[suitable_points$data > 0]
  suitable_cells <- terra::extract(suitable, suitable_points, cells = TRUE)$cell
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
  for (i in seq_len(terra::nrow(spatial_indices2))) {
    spatial_indices[[i]] <- spatial_indices2[i, 1:2]
  }
  spatial_indices <- unname(spatial_indices)
  config$spatial_indices <- spatial_indices

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
      host_mean <- terra::as.matrix(host[[1]], wide = TRUE)
      host_sd <- terra::as.matrix(host[[2]], wide = TRUE)
    } else {
      config$failure <- host_uncert_error
      return(config)
    }
  } else {
    host_mean <- terra::as.matrix(host, wide = TRUE)
    host_sd <- zero_matrix
  }
  config$host_mean <- host_mean
  config$host_sd <- host_sd

  if (!is.null(config$mask)) {
    if (config$function_name %in% c("casestudy_creation", "model_api")) {
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
      return(config)
    }
  } else {
    mask <- terra::classify(host, config$rclmat)
    config$mask <- mask
    config$mask_matrix <- terra::as.matrix(mask, wide = TRUE)
  }

  if (config$use_initial_condition_uncertainty) {
    if (terra::nlyr(infected) == 2) {
      infected_mean <- terra::as.matrix(infected[[1]], wide = TRUE)
      infected_sd <- terra::as.matrix(infected[[2]], wide = TRUE)
    } else {
      config$failure <- initial_cond_uncert_error
      return(config)
    }
  } else {
    infected_mean <- terra::as.matrix(infected, wide = TRUE)
    infected_sd <- zero_matrix
  }

  config$infected_mean <- infected_mean
  config$infected_sd <- infected_sd

  exposed[[config$latency_period + 1]] <- exposed_mean
  config$total_exposed <- exposed_mean
  config$exposed <- exposed

  susceptible_mean <- host_mean - infected_mean - exposed_mean
  susceptible_mean[susceptible_mean < 0] <- 0
  config$susceptible_mean <- terra::as.matrix(susceptible_mean, wide = TRUE)

  config$total_populations <- terra::as.matrix(total_populations, wide = TRUE)
  # config$total_hosts_mean <- terra::as.matrix(host_mean, wide = TRUE)
  config$mortality <- zero_matrix
  config$resistant <- zero_matrix

  # check that quarantine raster has the same crs, resolution, and extent
  if (config$use_quarantine) {
    if (config$function_name %in% c("casestudy_creation", "model_api")) {
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
      return(config)
    }
  } else {
    # set quarantine areas to all zeros. meaning no quarantine areas are considered
    config$quarantine_areas <- zero_matrix
  }

  mortality_tracker <- zero_matrix
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
    mortality_tracker2[[length(mortality_tracker2)]] <- infected_mean
  }


  config$mortality_tracker <- mortality_tracker2

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

    if (nrow(config$parameter_cov_matrix) != 8 |
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

    if (config$parameter_means[8] > (min(config$rows_cols$num_cols, config$rows_cols$num_rows) * config$res$ew_res)) {
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

  if (config$function_name %in% c("validate", "calibrate")) {
    config$use_anthropogenic_kernel <- TRUE
    # Load observed data on occurrence
    infection_years <- terra::rast(config$infected_years_file)
    infection_years[] <- as.integer(infection_years[])
    config$num_layers_infected_years <- terra::nlyr(infection_years)

    if (config$num_layers_infected_years < config$number_of_outputs) {
      config$failure <-
        infection_years_length_error(config$num_layers_infected_years, config$number_of_time_steps)
      return(config)
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

  return(config)
}
