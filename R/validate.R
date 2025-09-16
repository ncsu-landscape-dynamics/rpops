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
validate <- function(config) {

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
      while (any(config$total_hosts > config$total_populations, na.rm = TRUE) ||
            any(config$total_exposed > config$total_populations, na.rm = TRUE) ||
            any(config$total_infecteds > config$total_populations, na.rm = TRUE)) {
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
        network_min_distances = config$network_min_distances,
        network_max_distances = config$network_max_distances,
        network_filenames = config$network_filenames,
        network_movement_types = config$network_movement_types,
        network_weights = config$network_weights,
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
