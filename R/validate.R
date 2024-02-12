#' Validates the accuracy of the calibrated reproductive rate and dispersal
#' scales of the pops model.
#'
#' This function uses the quantity, allocation, and configuration disagreement
#' to validate the model across the landscape using the parameters from the
#' calibrate function. Ideally the model is calibrated with 2 or more years of
#' data and validated for the last year or if you have 6 or more years of data
#' then the model can be validated for the final 2 years.
#'
#' @inheritParams pops
#' @param infected_years_file years of initial infection/infestation as
#' individual locations of a pest or pathogen in raster format
#' @param number_of_iterations how many iterations do you want to run to allow
#' the calibration to converge at least 10
#' @param number_of_cores enter how many cores you want to use (default = NA).
#' If not set uses the # of CPU cores - 1. must be an integer >= 1
#' @param mask Raster file used to provide a mask to remove 0's that are not
#' true negatives from comparisons (e.g. mask out lakes and oceans from statics
#' if modeling terrestrial species).
#' @param parameter_means the parameter means from the abc calibration function
#' (posterior means)
#' @param parameter_cov_matrix the parameter covariance matrix from the ABC
#' calibration function (posterior covariance matrix)
#' @param write_outputs Either c("summary_outputs", "all_simulations", or "None"). If not
#' "None" output folder path must be provided.
#' @param output_folder_path this is the full path with either / or \\ (e.g.,
#' "C:/user_name/desktop/pops_sod_2020_2023/outputs/")
#' @param point_file  file for point comparison if not provided skips calculations
#' @param use_distance Boolean if you want to compare distance between
#' simulations and observations. Default is FALSE.
#' @param use_configuration Boolean if you want to use configuration
#' disagreement for comparing model runs. Default is FALSE.
#'
#' @importFrom terra app rast xres yres classify extract ext as.points ncol nrow project
#' nlyr rowFromCell colFromCell values as.matrix rowFromCell colFromCell crs vect
#' @importFrom stats runif rnorm
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar%
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom lubridate interval time_length mdy %within%
#' @importFrom MASS mvrnorm
#' @importFrom Metrics rmse
#' @importFrom utils write.csv read.table read.csv
#' @importFrom methods is
#'
#' @return a data frame of statistical measures of model performance.
#' @export
#'
validate <- function(infected_years_file,
                     number_of_iterations = 10,
                     number_of_cores = NA,
                     parameter_means,
                     parameter_cov_matrix,
                     pest_host_table,
                     competency_table,
                     infected_file_list,
                     host_file_list,
                     total_populations_file,
                     temp = FALSE,
                     temperature_coefficient_file = "",
                     precip = FALSE,
                     precipitation_coefficient_file = "",
                     model_type = "SI",
                     latency_period = 0,
                     time_step = "month",
                     season_month_start = 1,
                     season_month_end = 12,
                     start_date = "2008-01-01",
                     end_date = "2008-12-31",
                     use_survival_rates = FALSE,
                     survival_rate_month = 3,
                     survival_rate_day = 15,
                     survival_rates_file = "",
                     use_lethal_temperature = FALSE,
                     temperature_file = "",
                     lethal_temperature = -12.87,
                     lethal_temperature_month = 1,
                     mortality_frequency = "year",
                     mortality_frequency_n = 1,
                     management = FALSE,
                     treatment_dates = c(""),
                     treatments_file = "",
                     treatment_method = "ratio",
                     natural_kernel_type = "cauchy",
                     anthropogenic_kernel_type = "cauchy",
                     natural_dir = "NONE",
                     anthropogenic_dir = "NONE",
                     pesticide_duration = 0,
                     pesticide_efficacy = 1.0,
                     mask = NULL,
                     output_frequency = "year",
                     output_frequency_n = 1,
                     movements_file = "",
                     use_movements = FALSE,
                     start_exposed = FALSE,
                     generate_stochasticity = TRUE,
                     establishment_stochasticity = TRUE,
                     movement_stochasticity = TRUE,
                     dispersal_stochasticity = TRUE,
                     establishment_probability = 0.5,
                     dispersal_percentage = 0.99,
                     quarantine_areas_file = "",
                     use_quarantine = FALSE,
                     use_spreadrates = FALSE,
                     use_overpopulation_movements = FALSE,
                     overpopulation_percentage = 0,
                     leaving_percentage = 0,
                     leaving_scale_coefficient = 1,
                     exposed_file_list = "",
                     write_outputs = "None",
                     output_folder_path = "",
                     point_file = "",
                     network_filename = "",
                     network_movement = "walk",
                     use_distance = FALSE,
                     use_configuration = FALSE,
                     use_initial_condition_uncertainty = FALSE,
                     use_host_uncertainty = FALSE,
                     weather_type = "deterministic",
                     temperature_coefficient_sd_file = "",
                     precipitation_coefficient_sd_file = "",
                     dispersers_to_soils_percentage = 0,
                     quarantine_directions = "",
                     multiple_random_seeds = FALSE,
                     file_random_seeds = NULL,
                     use_soils = FALSE,
                     soil_starting_pest_file = "",
                     start_with_soil_populations = FALSE,
                     county_level_infection_data = FALSE) {
  config <- c()
  config$infected_years_file <- infected_years_file
  config$infected_file_list <- infected_file_list
  config$host_file_list <- host_file_list
  config$total_populations_file <- total_populations_file
  config$parameter_means <- parameter_means
  config$parameter_cov_matrix <- parameter_cov_matrix
  config$temp <- temp
  config$temperature_coefficient_file <- temperature_coefficient_file
  config$precip <- precip
  config$precipitation_coefficient_file <- precipitation_coefficient_file
  config$model_type <- model_type
  config$latency_period <- latency_period
  config$time_step <- time_step
  config$season_month_start <- season_month_start
  config$season_month_end <- season_month_end
  config$start_date <- start_date
  config$end_date <- end_date
  config$use_lethal_temperature <- use_lethal_temperature
  config$temperature_file <- temperature_file
  config$lethal_temperature <- lethal_temperature
  config$lethal_temperature_month <- lethal_temperature_month
  config$use_survival_rates <- use_survival_rates
  config$survival_rate_month <- survival_rate_month
  config$survival_rate_day <- survival_rate_day
  config$survival_rates_file <- survival_rates_file
  config$management <- management
  config$treatment_dates <- treatment_dates
  config$treatments_file <- treatments_file
  config$treatment_method <- treatment_method
  config$natural_kernel_type <- natural_kernel_type
  config$anthropogenic_kernel_type <- anthropogenic_kernel_type
  config$natural_dir <- natural_dir
  config$anthropogenic_dir <- anthropogenic_dir
  config$pesticide_duration <- pesticide_duration
  config$pesticide_efficacy <- pesticide_efficacy
  config$mask <- mask
  config$output_frequency <- output_frequency
  config$output_frequency_n <- output_frequency_n
  config$movements_file <- movements_file
  config$use_movements <- use_movements
  config$start_exposed <- start_exposed
  config$generate_stochasticity <- generate_stochasticity
  config$establishment_stochasticity <- establishment_stochasticity
  config$movement_stochasticity <- movement_stochasticity
  config$dispersal_stochasticity <- dispersal_stochasticity
  config$establishment_probability <- establishment_probability
  config$dispersal_percentage <- dispersal_percentage
  config$quarantine_areas_file <- quarantine_areas_file
  config$quarantine_directions <- quarantine_directions
  config$use_quarantine <- use_quarantine
  config$use_spreadrates <- use_spreadrates
  config$use_overpopulation_movements <- use_overpopulation_movements
  config$overpopulation_percentage <- overpopulation_percentage
  config$leaving_percentage <- leaving_percentage
  config$leaving_scale_coefficient <- leaving_scale_coefficient
  config$number_of_iterations <- number_of_iterations
  config$number_of_cores <- number_of_cores
  # add function name for use in configuration function to skip
  # function specific specific configurations namely for validation and
  # calibration.
  config$function_name <- "validate"
  config$failure <- NULL
  config$exposed_file_list <- exposed_file_list
  config$write_outputs <- write_outputs
  config$output_folder_path <- output_folder_path
  config$mortality_frequency <- mortality_frequency
  config$mortality_frequency_n <- mortality_frequency_n
  config$point_file <- point_file
  config$network_filename <- network_filename
  config$network_movement <- network_movement
  config$use_configuration <- use_configuration
  config$use_distance <- use_distance
  config$use_initial_condition_uncertainty <- use_initial_condition_uncertainty
  config$use_host_uncertainty <- use_host_uncertainty
  config$weather_type <- weather_type
  config$temperature_coefficient_sd_file <- temperature_coefficient_sd_file
  config$precipitation_coefficient_sd_file <- precipitation_coefficient_sd_file
  config$dispersers_to_soils_percentage <- dispersers_to_soils_percentage
  config$multiple_random_seeds <- multiple_random_seeds
  config$file_random_seeds <- file_random_seeds
  config$use_soils <- use_soils
  config$soil_starting_pest_file <- soil_starting_pest_file
  config$start_with_soil_populations <- start_with_soil_populations
  config$county_level_infection_data <- county_level_infection_data
  config$pest_host_table <- pest_host_table
  config$competency_table <- competency_table

  config <- configuration(config)

  if (!is.null(config$failure)) {
    stop(config$failure)
  }

  if (config$multiple_random_seeds && is.null(config$file_random_seeds) &&
      dir.exists(config$output_folder_path)) {
    write.csv(config$random_seeds, paste0(config$output_folder_path, "validation_random_seeds.csv"),
              row.names = FALSE)
  }

  i <- NULL

  cl <- makeCluster(config$core_count)
  registerDoParallel(cl)

  qa <-
    foreach::foreach(
      i = 1:number_of_iterations,
      .combine = rbind,
      .packages = c("terra", "PoPS", "foreach")
    ) %dopar% {

      set.seed(config$random_seed[[i]])
      config <- draw_parameters(config) # draws parameter set for the run
      config <- host_pool_setup(config)
      while (any(config$total_hosts > config$total_populations) ||
            any(config$total_exposed > config$total_populations) ||
            any(config$total_infecteds > config$total_populations)) {
        config <- host_pool_setup(config)
      }
      config$competency_table_list <- competency_table_list_creator(config$competency_table)
      config$pest_host_table_list <- pest_host_table_list_creator(config$pest_host_table)

      data <- pops_model(
        random_seed = config$random_seed[i],
        multiple_random_seeds = config$multiple_random_seeds,
        random_seeds = as.matrix(config$random_seeds[i, ])[1, ],
        use_lethal_temperature = config$use_lethal_temperature,
        lethal_temperature = config$lethal_temperature,
        lethal_temperature_month = config$lethal_temperature_month,
        use_survival_rates = config$use_survival_rates,
        survival_rate_month = config$survival_rate_month,
        survival_rate_day = config$survival_rate_day,
        host_pools = config$host_pools,
        total_populations = config$total_populations,
        competency_table = config$competency_table_list,
        pest_host_table = config$pest_host_table_list,
        mortality_on = config$mortality_on,
        quarantine_areas = config$quarantine_areas,
        quarantine_directions = config$quarantine_directions,
        treatment_maps = config$treatment_maps,
        treatment_dates = config$treatment_dates,
        pesticide_duration = config$pesticide_duration,
        use_movements = config$use_movements,
        movements = config$movements,
        movements_dates = config$movements_dates,
        weather = config$weather,
        temperature = config$temperature,
        survival_rates = config$survival_rates,
        weather_coefficient = config$weather_coefficient,
        weather_coefficient_sd = config$weather_coefficient_sd,
        res = config$res,
        rows_cols = config$rows_cols,
        time_step = config$time_step,
        reproductive_rate = config$reproductive_rate,
        spatial_indices = config$spatial_indices,
        season_month_start_end = config$season_month_start_end,
        soil_reservoirs = config$soil_reservoirs,
        start_date = config$start_date,
        end_date = config$end_date,
        treatment_method = config$treatment_method,
        natural_kernel_type = config$natural_kernel_type,
        anthropogenic_kernel_type = config$anthropogenic_kernel_type,
        use_anthropogenic_kernel = config$use_anthropogenic_kernel,
        percent_natural_dispersal = config$percent_natural_dispersal,
        natural_distance_scale = config$natural_distance_scale,
        anthropogenic_distance_scale = config$anthropogenic_distance_scale,
        natural_dir = config$natural_dir,
        natural_kappa = config$natural_kappa,
        anthropogenic_dir = config$anthropogenic_dir,
        anthropogenic_kappa = config$anthropogenic_kappa,
        output_frequency = config$output_frequency,
        output_frequency_n = config$output_frequency_n,
        quarantine_frequency = config$quarantine_frequency,
        quarantine_frequency_n = config$quarantine_frequency_n,
        use_quarantine = config$use_quarantine,
        spreadrate_frequency = config$spreadrate_frequency,
        spreadrate_frequency_n = config$spreadrate_frequency_n,
        mortality_frequency = config$mortality_frequency,
        mortality_frequency_n = config$mortality_frequency_n,
        use_spreadrates = config$use_spreadrates,
        model_type_ = config$model_type,
        latency_period = config$latency_period,
        generate_stochasticity = config$generate_stochasticity,
        establishment_stochasticity = config$establishment_stochasticity,
        movement_stochasticity = config$movement_stochasticity,
        dispersal_stochasticity = config$dispersal_stochasticity,
        establishment_probability = config$establishment_probability,
        dispersal_percentage = config$dispersal_percentage,
        use_overpopulation_movements = config$use_overpopulation_movements,
        overpopulation_percentage = config$overpopulation_percentage,
        leaving_percentage = config$leaving_percentage,
        leaving_scale_coefficient = config$leaving_scale_coefficient,
        bbox = config$bounding_box,
        network_min_distance = config$network_min_distance,
        network_max_distance = config$network_max_distance,
        network_filename = config$network_filename,
        network_movement = config$network_movement,
        weather_size = config$weather_size,
        weather_type = config$weather_type,
        dispersers_to_soils_percentage = config$dispersers_to_soils_percentage,
        use_soils = config$use_soils)

      all_disagreement <- calculate_all_stats(config, data)
      all_disagreement
    }

  parallel::stopCluster(cl)

  output_list <- list()
  for (j in 1:max(qa$output)) {
    output_step <- qa[qa$output == j, ]
    assign(paste("output_step_", j, sep = ""), output_step)
    output_list[[paste0("output_step_", j)]] <- output_step
    if (config$write_outputs %in% config$output_write_list) {
      file_name <- paste(config$output_folder_path, "output_step_", j, ".csv", sep = "")
      write.csv(output_step, file = file_name)
    }
    if (j == 1) {
      cum_output_step <- output_step
      assign(paste("cum_output_step_", j, sep = ""), cum_output_step / j)
      output_list[[paste0("cum_output_step_", j)]] <- cum_output_step
      if (config$write_outputs %in% config$output_write_list) {
        file_name <- paste(config$output_folder_path, "cum_output_step_", j, ".csv", sep = "")
        write.csv(cum_output_step, file = file_name)
      }
    }
    else {
      assign(paste("cum_output_step", sep = ""), (cum_output_step + output_step))
      assign(paste("cum_output_step_", j, sep = ""), cum_output_step / j)
      output_list[[paste0("cum_output_step_", j)]] <- cum_output_step
      if (config$write_outputs %in% config$output_write_list) {
        file_name <- paste(config$output_folder_path, "cum_output_step_", j, ".csv", sep = "")
        write.csv(cum_output_step, file = file_name)
      }
    }
  }

  if (config$write_outputs %in% config$output_write_list) {
    file_name <- paste(config$output_folder_path, "validation_outputs.rdata", sep = "")
    save(output_list, file = file_name)
  }

  return(output_list)
}
