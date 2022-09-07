#' PoPS configuration
#'
#' Function with a single input and output list for parsing, transforming,
#' and performing all checks for all functions to run the pops c++ model
#'
#' @param config the .xlsx or .csv file that contains all of the values needed to run the model.
#' This file should be in the same folder with all of the files needed to run the simulation.
#'
#' @importFrom raster cellStats calc extract
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

configuration <- function(config_file) {

  config_file <- "C:/Users/Chris/Desktop/pops_test/pops_model_inputs.xlsx"
  config_data <- readxl::read_xlsx(path = config_file , sheet = 1)

  config <- c()
  config$failure <- NULL ## maybe remove
  config$input_folder <- paste0(dirname(config_file), "/")
  config$output_folder_path <- paste(config$input_folder, "outputs/", sep = "")
  config$calibrate_folder_path <- paste(config$input_folder, "calibrate/", sep = "")
  config$validate_folder_path <- paste(config$input_folder, "validate/", sep = "")
  config$forecast_folder_path <- paste(config$input_folder, "forecast/", sep = "")
  if (!base::dir.exists(config$output_folder_path)) {
    stop(output_path_error)
  }
  if (!base::dir.exists(config$calibrate_folder_path)) {
    stop(calibrate_path_error)
  }
  if (!base::dir.exists(config$validate_folder_path)) {
    stop(validate_path_error)
  }
  if (!base::dir.exists(config$forecast_folder_path)) {
    stop(forecast_path_error)
  }

  # Pest/pathogen information
  config$species_name <- config_data$value[config_data$variable_name == 'species_name']
  config$latency_period <-
    as.integer(config_data$value[config_data$variable_name == 'latency_period'])
  config$start_exposed <-
    as.logical(config_data$value[config_data$variable_name == 'start_exposed'])

  # Dates and timings
  config$time_step <- config_data$value[config_data$variable_name == 'time_step']
  config$output_frequency <- config_data$value[config_data$variable_name == 'output_frequency']
  config$output_frequency_n <-
    as.integer(config_data$value[config_data$variable_name == 'output_frequency_n'])

  # Host data
  config$host_species_name <- config_data$value[config_data$variable_name == 'host_species_name']
  config$mortality_on <- as.logical(config_data$value[config_data$variable_name == 'mortality_on'])
  config$mortality_rate <-
    as.numeric(config_data$value[config_data$variable_name == 'mortality_rate'])
  config$mortality_time_lag <-
    as.integer(config_data$value[config_data$variable_name == 'mortality_time_lag'])
  config$mortality_frequency <-
    config_data$value[config_data$variable_name == 'mortality_frequency']
  config$mortality_frequency_n <-
    as.integer(config_data$value[config_data$variable_name == 'mortality_frequency_n'])

  # Weather values
  config$season_month_start <-
    as.integer(config_data$value[config_data$variable_name == 'season_month_start'])
  config$season_month_end <-
    as.integer(config_data$value[config_data$variable_name == 'season_month_end'])

  if (config$season_month_start %in% seasons && config$season_month_end %in% seasons) {
    season_month_start_end <- c()
    season_month_start_end$start_month <- as.integer(config$season_month_start)
    season_month_start_end$end_month <- as.integer(config$season_month_end)
    config$season_month_start_end <- season_month_start_end
  } else {
    stop(season_month_error)
  }

  config$use_temperature <-
    as.logical(config_data$value[config_data$variable_name == 'use_temperature'])
  config$use_precipitation <-
    as.logical(config_data$value[config_data$variable_name == 'use_precipitation'])
  config$use_lethal_temperature <-
    as.logical(config_data$value[config_data$variable_name == 'use_lethal_temperature'])
  config$lethal_temperature <-
    as.numeric(config_data$value[config_data$variable_name == 'lethal_temperature'])
  config$lethal_temperature_month <-
    as.integer(config_data$value[config_data$variable_name == 'lethal_temperature_month'])

  # Survival Rates
  config$use_survival_rates <-
    as.logical(config_data$value[config_data$variable_name == 'use_survival_rates'])
  config$survival_rate_month <-
    as.integer(config_data$value[config_data$variable_name == 'survival_rate_month'])
  config$survival_rate_day <-
    as.integer(config_data$value[config_data$variable_name == 'survival_rate_day'])

  # Dispersal information (doesn't change across stages)
  config$natural_kernel_type <-
    config_data$value[config_data$variable_name == 'natural_kernel_type']
  config$anthropogenic_kernel_type <-
    config_data$value[config_data$variable_name == 'anthropogenic_kernel_type']

  if (config$natural_kernel_type %notin% kernel_list) {
    stop(natural_kernel_error)
  }
  if (config$anthropogenic_kernel_type %notin% kernel_list) {
    stop(anthropogenic_kernel_error)
  }

  config$natural_dir <-
    config_data$value[config_data$variable_name == 'natural_direction']
  config$anthropogenic_dir <-
    config_data$value[config_data$variable_name == 'anthropogenic_direction']

  if (config$natural_dir %notin% directions) {
    stop(natural_direction_error)
  }

  if (config$anthropogenic_dir %notin% directions) {
    stop(anthropogenic_direction_error)
  }

  config$network_file <-
    paste0(config$input_folder, config_data$value[config_data$variable_name == 'network_file'])
  config$network_movement <-
    config_data$value[config_data$variable_name == 'network_movement']
  # check that network movement is one of the correct options
  if (config$network_movement %notin% network_movement_options) {
    config$failure <- network_movement_error
  }

  # Overpopulation Movements
  config$use_overpopulation_movements <-
    as.logical(config_data$value[config_data$variable_name == 'use_overpopulation_movements'])
  config$overpopulation_percentage <-
    as.numeric(config_data$value[config_data$variable_name == 'overpopulation_percentage'])
  config$leaving_percentage <-
    as.numeric(config_data$value[config_data$variable_name == 'leaving_percentage'])
  config$leaving_scale_coefficient <-
    as.numeric(config_data$value[config_data$variable_name == 'leaving_scale_coefficient'])

  # Movement information (the movementfile changes across stages)
  config$use_movements <-
    as.logical(config_data$value[config_data$variable_name == 'use_movements'])

  # Model information (parameters change across stages)
  config$model_type <- config_data$value[config_data$variable_name == 'model_type']
  if (config$model_type %in% sei_model_names) {
    config$model_type <- "SEI"
  } else if (config$model_type %in% si_model_names) {
    config$model_type <- "SI"
  } else {
    stop(model_type_error)
  }

  # ensures latent period is correct for type of model selected
  if (config$model_type == "SEI" && config$latency_period <= 0) {
    stop(latency_period_error)
  } else if (config$model_type == "SI" && config$latency_period > 0) {
    config$latency_period <- 0
  }
  config$mask <- config_data$value[config_data$variable_name == 'mask']

  # determine if you should calculate additional statistics
  config$use_spreadrates <-
    as.logical(config_data$value[config_data$variable_name == 'use_spreadrates'])
  config$use_quarantine <-
    as.logical(config_data$value[config_data$variable_name == 'use_quarantine'])
  config$quarantine_areas_file <-
    paste0(config$input_folder,
           config_data$value[config_data$variable_name == 'quarantine_areas_file'])

  # Variables for switching between stochastic and deterministic simulations.
  config$generate_stochasticity <-
    as.logical(config_data$value[config_data$variable_name == 'generate_stochasticity'])
  config$movement_stochasticity <-
    as.logical(config_data$value[config_data$variable_name == 'movement_stochasticity'])
  config$dispersal_stochasticity <-
    as.logical(config_data$value[config_data$variable_name == 'dispersal_stochasticity'])
  config$dispersal_percentage <-
    as.numeric(config_data$value[config_data$variable_name == 'dispersal_percentage'])
  config$establishment_stochasticity <-
    as.logical(config_data$value[config_data$variable_name == 'establishment_stochasticity'])
  config$establishment_probability <-
    as.numeric(config_data$value[config_data$variable_name == 'establishment_probability'])
  if (is.na(config_data$value[config_data$variable_name == 'random_seed'])) {
    config$random_seed <- round(stats::runif(1, 1, 1000000))
  } else {
    config$random_seed <- as.integer(config_data$value[config_data$variable_name == 'random_seed'])
  }

  # S3 bucket
  config$use_s3 <- as.logical(config_data$value[config_data$variable_name == 'use_s3'])

  # Control Measures (change dates and files change across stages)
  config$use_management <-
    as.logical(config_data$value[config_data$variable_name == 'use_management'])
  # config$treatment_dates <- config_data$value[config_data$variable_name == 'treatment_dates']
  config$treatments_file <-
    paste0(config$input_folder, config_data$value[config_data$variable_name == 'treatments_file'])
  config$treatment_method <- config_data$value[config_data$variable_name == 'treatment_method']
  # config$pesticide_duration <-
  #   as.integer(config_data$value[config_data$variable_name == 'pesticide_duration'])
  # config$pesticide_efficacy <-
  #   as.numeric(config_data$value[config_data$variable_name == 'pesticide_efficacy'])

  # ensure correct treatment method
  if (!config$treatment_method %in% c("ratio", "all infected")) {
    stop(treatment_option_error)
  }

  # reclass
  config$rcl <- c(1, Inf, 1, 0, 0.99, NA)
  config$rclmat <- matrix(config$rcl, ncol = 3, byrow = TRUE)

  # check that initial raster file exists
  forecast_infected_file <- paste0(config$forecast_folder_path,
                      config_data$value[config_data$variable_name == 'forecast_infection_data'])

  if (config$use_s3) {
    infected_check <- initial_raster_checks(forecast_infected_file, config$use_s3, config$bucket)
  } else {
    infected_check <- initial_raster_checks(forecast_infected_file)
  }
  if (infected_check$checks_passed) {
    infected <- infected_check$raster
    if (terra::nlyr(infected) > 1) {
      infected <- output_from_raster_mean_and_sd(infected)
    }
    infected <- terra::classify(infected, matrix(c(NA, 0), ncol = 2, byrow = TRUE), right = NA)
  } else {
    stop(infected_check$failed_check)
  }

  config$crs <- terra::crs(infected)
  config$xmax <- terra::xmax(infected)
  config$xmin <- terra::xmin(infected)
  config$ymax <- terra::ymax(infected)
  config$ymin <- terra::ymin(infected)
  bounding_box <- c()
  bounding_box$north <- config$ymax
  bounding_box$south <- config$ymin
  bounding_box$west <- config$xmin
  bounding_box$east <- config$xmax
  config$bounding_box <- bounding_box



  # set up mask matrix
  if (!file.exists(config$mask)) {
    if (config$use_s3) {
      mask_check <- secondary_raster_checks(config$mask, infected, config$use_s3, config$bucket)
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
    stop(time_check$failed_check)
  }

  # calibrate data
  config$start_calibration_date <-
    config_data$value[config_data$variable_name == 'start_calibration_date']
  config$end_calibration_date <-
    config_data$value[config_data$variable_name == 'end_calibration_date']

  config$host_file <-
    paste0(config$input_folder, config_data$value[config_data$variable_name == 'host_file'])
  config$total_populations_file <-
    paste0(config$input_folder,
           config_data$value[config_data$variable_name == 'total_populations_file'])
  config$temperature_coefficient_file <-
    paste0(config$input_folder,
           config_data$value[config_data$variable_name == 'temperature_coefficient_file'])
  config$precipitation_coefficient_file <-
    paste0(config$input_folder,
           config_data$value[config_data$variable_name == 'precipitation_coefficient_file'])
  config$temperature_file <-
    paste0(config$input_folder, config_data$value[config_data$variable_name == 'temperature_file'])

  config$movements_file <- config_data$value[config_data$variable_name == 'movements_file']

  # calibration outputs
  config$parameter_means <- as.numeric(read.csv(
    paste0(config$input_folder, config_data$value[config_data$variable_name == 'parameter_means']),
    header = FALSE))
  config$parameter_cov_matrix <- as.matrix(read.csv(
    paste0(config$input_folder,config_data$value[config_data$variable_name == 'parameter_cov_matrix']),
    header = FALSE))

  # validate data

  # forecast data
  config$start_forecast_date <-
    config_data$value[config_data$variable_name == 'start_forecast_date']
  config$end_forecast_date <-
    config_data$value[config_data$variable_name == 'end_forecast_date']
  config$forecast_infected_file <-
    paste0(config$input_folder,
           config_data$value[config_data$variable_name == 'forecast_infection_data'])


  if (is.na(config$number_of_cores) || config$number_of_cores > parallel::detectCores()) {
    config$core_count <- parallel::detectCores() - 1
  } else {
    config$core_count <- config$number_of_cores
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

  # check that host raster has the same crs, resolution, and extent
  if (config$function_name %in% c("casestudy_creation", "model_api")) {
    host_check <- secondary_raster_checks(config$host_file, infected, config$use_s3, config$bucket)
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

  suitable <- host + infected
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

  # check that total populations raster has the same crs, resolution, and extent
  if (config$function_name %in% c("casestudy_creation", "model_api")) {
    total_populations_check <-
      secondary_raster_checks(config$total_populations_file, infected, config$use_s3, config$bucket)
  } else {
    total_populations_check <- secondary_raster_checks(config$total_populations_file, infected)
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
    survival_rates <- host
    terra::values(survival_rates) <- 1
    survival_rates <- list(terra::as.matrix(survival_rates, wide = TRUE))
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
    temperature <- host
    terra::values(temperature) <- 1
    temperature <- list(terra::as.matrix(temperature, wide = TRUE))
  }

  config$temperature <- temperature

  # check that temp and precip rasters have the same crs, resolution, and extent
  config$weather <- FALSE
  if (config$use_temperature == TRUE) {
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
    weather_coefficient <- host
    terra::values(weather_coefficient) <- 1
    weather_coefficient <- list(terra::as.matrix(weather_coefficient, wide = TRUE))
  }

  config$weather_coefficient <- weather_coefficient


  ## Checks for management
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
    treatment_map <- host
    treatment_map[] <- 0
    config$treatment_maps <- list(terra::as.matrix(treatment_map, wide = TRUE))
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


  mortality_tracker <- infected
  terra::values(mortality_tracker) <- 0
  mortality_tracker <- terra::as.matrix(mortality_tracker, wide = TRUE)
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
        secondary_raster_checks(config$exposed_file, infected, config$use_s3, config$bucket)
    } else {
      exposed_check <- secondary_raster_checks(config$exposed_file, infected)
    }
    if (exposed_check$checks_passed) {
      exposed2 <- exposed_check$raster
      susceptible <- susceptible - exposed2
      susceptible[susceptible < 0] <- 0
      exposed[[config$latency_period + 1]] <- terra::as.matrix(exposed2, wide = TRUE)
      config$total_exposed <- terra::as.matrix(exposed2, wide = TRUE)
    } else {
      config$failure <- exposed_check$failed_check
      return(config)
    }
  }

  infected <- terra::as.matrix(infected, wide = TRUE)
  config$susceptible <- terra::as.matrix(susceptible, wide = TRUE)
  config$total_populations <- terra::as.matrix(total_populations, wide = TRUE)
  config$total_hosts <- terra::as.matrix(host, wide = TRUE)
  config$mortality <- mortality_tracker
  config$resistant <- mortality_tracker

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

  if (config$function_name %in% c("calibrate") && config$calibration_method == "ABC") {
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

  return(config)
}
