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
#' @importFrom terra app rast xres yres classify extract ext as.points ncol nrow
#' nlyr rowFromCell colFromCell values as.matrix rowFromCell colFromCell crs vect
#' @importFrom stats runif rnorm
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar%
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom lubridate interval time_length mdy %within%
#' @importFrom MASS mvrnorm
#' @importFrom Metrics rmse
#' @importFrom utils write.csv
#'
#' @return a data frame of statistical measures of model performance.
#' @export
#'
validate <- function(config) {

  i <- NULL

  cl <- makeCluster(config$core_count)
  registerDoParallel(cl)

  qa <-
    foreach::foreach(
      i = 1:number_of_iterations,
      .combine = rbind,
      .packages = c("terra", "PoPS", "foreach")
    ) %dopar% {

      config$random_seed <- round(stats::runif(1, 1, 1000000))
      config <- draw_parameters(config)

      data <- pops_model(
        random_seed = config$random_seed,
        use_lethal_temperature = config$use_lethal_temperature,
        lethal_temperature = config$lethal_temperature,
        lethal_temperature_month = config$lethal_temperature_month,
        use_survival_rates = config$use_survival_rates,
        survival_rate_month = config$survival_rate_month,
        survival_rate_day = config$survival_rate_day,
        infected = config$infected,
        total_exposed = config$total_exposed,
        exposed = config$exposed,
        susceptible = config$susceptible,
        total_populations = config$total_populations,
        total_hosts = config$total_hosts,
        mortality_on = config$mortality_on,
        mortality_tracker = config$mortality_tracker,
        mortality = config$mortality,
        quarantine_areas = config$quarantine_areas,
        treatment_maps = config$treatment_maps,
        treatment_dates = config$treatment_dates,
        pesticide_duration = config$pesticide_duration,
        resistant = config$resistant,
        use_movements = config$use_movements,
        movements = config$movements,
        movements_dates = config$movements_dates,
        weather = config$weather,
        temperature = config$temperature,
        survival_rates = config$survival_rates,
        weather_coefficient = config$weather_coefficient,
        res = config$res,
        rows_cols = config$rows_cols,
        time_step = config$time_step,
        reproductive_rate = config$reproductive_rate,
        spatial_indices = config$spatial_indices,
        season_month_start_end = config$season_month_start_end,
        mortality_rate = config$mortality_rate,
        mortality_time_lag = config$mortality_time_lag,
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
        network_movement = config$network_movement
      )


      all_disagreement <-
        foreach(
          q = seq_len(length(data$infected)), .combine = rbind,
          .packages = c("terra", "PoPS")
        ) %do% {
          # need to assign reference, comparison, and mask in inner loop since
          # terra objects are pointers and pointers using %dopar%
          comparison <- terra::rast(config$infected_file)
          reference <- terra::rast(config$infected_file)
          terra::values(comparison) <- data$infected[[q]]
          terra::values(reference) <- config$infection_years2[[q]]
          mask <- terra::rast(config$infected_file)
          terra::values(mask) <- config$mask_matrix
          ad <-
            quantity_allocation_disagreement(reference,
                                             comparison,
                                             use_configuration = config$use_configuration,
                                             mask = mask,
                                             use_distance = config$use_distance)
          if (file.exists(config$point_file)) {
            obs_data <- terra::vect(config$point_file)
            obs_data <- terra::project(obs_data, comparison)
            s <- extract(comparison, obs_data)
            names(s) <- c("ID", paste("sim_value_output_", q, sep = ""))
            s <- s[2]
            obs_data <- cbind(obs_data, s)
            ## calculate true positive, true negatives, false positives, false
            ## negatives, and other statistics and add them to the data frame
            ## for export
            ad$points_true_positive <-
              nrow(obs_data[obs_data$positive > 0 & obs_data$sim_value_output_1 > 0, ])
            ad$points_false_negative <-
              nrow(obs_data[obs_data$positive > 0 & obs_data$sim_value_output_1 == 0, ])
            ad$points_false_positive <-
              nrow(obs_data[obs_data$positive == 0 & obs_data$sim_value_output_1 > 0, ])
            ad$points_true_negative <-
              nrow(obs_data[obs_data$positive == 0 & obs_data$sim_value_output_1 == 0, ])
            ad$points_total_obs <-
              ad$points_true_negative + ad$points_true_positive +
              ad$points_false_negative + ad$points_false_positive
            ad$points_accuracy <-
              (ad$points_true_negative + ad$points_true_positive) / ad$points_total_obs
            ad$points_precision <-
              ad$points_true_positive / (ad$points_true_positive + ad$points_false_positive)
            ad$points_recall <-
              ad$points_true_positive / (ad$points_true_positive + ad$points_false_negative)
            ad$points_specificiity <-
              ad$points_true_negative / (ad$points_true_negative + ad$points_false_positive)
          }
          ad$ouput <- q
          ad
        }

      data.frame(all_disagreement)
    }

  parallel::stopCluster(cl)

  output_list <- list()
  for (j in 1:max(qa$ouput)) {
    output_step <- qa[qa$ouput == j, ]
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
