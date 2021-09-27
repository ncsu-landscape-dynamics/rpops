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
#' @param success_metric Choose which success metric to use for calibration.
#' Choices are "quantity", "quantity and configuration", "residual error" and
#' "odds ratio". Default is "quantity"
#' @param mask Raster file used to provide a mask to remove 0's that are not
#' true negatives from comparisons (e.g. mask out lakes and oceans from statics
#' if modeling terrestrial species).
#' @param parameter_means the parameter means from the abc calibration function
#' (posterior means)
#' @param parameter_cov_matrix the parameter covariance matrix from the abc
#' calibration function (posterior covairance matrix)
#' @param write_outputs Either c("summary_outputs", or "None"). If not
#' "None" output folder path must be provided.
#' @param output_folder_path this is the full path with either / or \\ (e.g.,
#' "C:/user_name/desktop/pops_sod_2020_2023/outputs/")
#' @param point_file  file for point comparison if not provided skips
#' calculations
#'
#' @importFrom terra app rast xres yres classify extract ext as.points ncol nrow
#' nlyr rowFromCell colFromCell values as.matrix rowFromCell colFromCell crs
#' @importFrom stats runif rnorm
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar%
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom lubridate interval time_length mdy %within%
#' @importFrom MASS mvrnorm
#'
#' @return a dataframe of the variables saved and their success metrics for
#' each run
#' @export
#'
validate <- function(infected_years_file,
                     number_of_iterations = 10,
                     number_of_cores = NA,
                     parameter_means,
                     parameter_cov_matrix,
                     infected_file,
                     host_file,
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
                     use_lethal_temperature = FALSE,
                     temperature_file = "",
                     lethal_temperature = -12.87,
                     lethal_temperature_month = 1,
                     mortality_on = FALSE,
                     mortality_rate = 0,
                     mortality_time_lag = 0,
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
                     success_metric = "quantity",
                     output_frequency = "year",
                     output_frequency_n = 1,
                     movements_file = "",
                     use_movements = FALSE,
                     start_exposed = FALSE,
                     generate_stochasticity = TRUE,
                     establishment_stochasticity = TRUE,
                     movement_stochasticity = TRUE,
                     deterministic = FALSE,
                     establishment_probability = 0.5,
                     dispersal_percentage = 0.99,
                     quarantine_areas_file = "",
                     use_quarantine = FALSE,
                     use_spreadrates = FALSE,
                     use_overpopulation_movements = FALSE,
                     overpopulation_percentage = 0,
                     leaving_percentage = 0,
                     leaving_scale_coefficient = 1,
                     exposed_file = "",
                     write_outputs = "None",
                     output_folder_path = "",
                     point_file = "") {
  config <- c()
  config$infected_years_file <- infected_years_file
  config$infected_file <- infected_file
  config$host_file <- host_file
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
  config$mortality_on <- mortality_on
  config$mortality_rate <- mortality_rate
  config$mortality_time_lag <- mortality_time_lag
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
  config$success_metric <- success_metric
  config$output_frequency <- output_frequency
  config$output_frequency_n <- output_frequency_n
  config$movements_file <- movements_file
  config$use_movements <- use_movements
  config$start_exposed <- start_exposed
  config$generate_stochasticity <- generate_stochasticity
  config$establishment_stochasticity <- establishment_stochasticity
  config$movement_stochasticity <- movement_stochasticity
  config$deterministic <- deterministic
  config$establishment_probability <- establishment_probability
  config$dispersal_percentage <- dispersal_percentage
  config$quarantine_areas_file <- quarantine_areas_file
  config$use_quarantine <- use_quarantine
  config$use_spreadrates <- use_spreadrates
  config$use_overpopulation_movements <- use_overpopulation_movements
  config$overpopulation_percentage <- overpopulation_percentage
  config$leaving_percentage <- leaving_percentage
  config$leaving_scale_coefficient <- leaving_scale_coefficient
  config$number_of_iterations <- number_of_iterations
  config$number_of_cores <- number_of_cores
  # add function name for use in configuration function to skip
  # function specific specifc configurations namely for validation and
  # calibration.
  config$function_name <- "validate"
  config$failure <- NULL
  config$exposed_file <- exposed_file
  config$write_outputs <- write_outputs
  config$output_folder_path <- output_folder_path
  config$mortality_frequency <- mortality_frequency
  config$mortality_frequency_n <- mortality_frequency_n
  config$point_file <- point_file

  config <- configuration(config)

  if (!is.null(config$failure)) {
    stop(config$failure)
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

      config$random_seed <- round(stats::runif(1, 1, 1000000))
      data <- pops_model(
        random_seed = config$random_seed,
        use_lethal_temperature = config$use_lethal_temperature,
        lethal_temperature = config$lethal_temperature,
        lethal_temperature_month =
          config$lethal_temperature_month,
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
        weather_coefficient = config$weather_coefficient,
        res = config$res,
        rows_cols = config$rows_cols,
        time_step = config$time_step,
        reproductive_rate = config$reproductive_rate[i],
        spatial_indices = config$spatial_indices,
        season_month_start_end = config$season_month_start_end,
        mortality_rate = config$mortality_rate,
        mortality_time_lag = config$mortality_time_lag,
        start_date = config$start_date,
        end_date = config$end_date,
        treatment_method = config$treatment_method,
        natural_kernel_type = config$natural_kernel_type,
        anthropogenic_kernel_type =
          config$anthropogenic_kernel_type,
        use_anthropogenic_kernel =
          config$use_anthropogenic_kernel,
        percent_natural_dispersal =
          config$percent_natural_dispersal[i],
        natural_distance_scale =
          config$natural_distance_scale[i],
        anthropogenic_distance_scale =
          config$anthropogenic_distance_scale[i],
        natural_dir = config$natural_dir,
        natural_kappa = config$natural_kappa[i],
        anthropogenic_dir = config$anthropogenic_dir,
        anthropogenic_kappa = config$anthropogenic_kappa[i],
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
        generate_stochasticity =
          config$generate_stochasticity,
        establishment_stochasticity =
          config$establishment_stochasticity,
        movement_stochasticity = config$movement_stochasticity,
        deterministic = config$deterministic,
        establishment_probability =
          config$establishment_probability,
        dispersal_percentage = config$dispersal_percentage,
        use_overpopulation_movements = config$use_overpopulation_movements,
        overpopulation_percentage = config$overpopulation_percentage,
        leaving_percentage = config$leaving_percentage,
        leaving_scale_coefficient = config$leaving_scale_coefficient
      )


      all_disagreement <-
        foreach(
          q = seq_len(length(data$infected)), .combine = rbind,
          .packages = c("terra", "PoPS")
        ) %do% {
          # need to assign reference, comp_year, and mask in inner loop since
          # terra objects are pointers and pointers using %dopar%
          comp_year <- terra::rast(config$infected_file)
          reference <- terra::rast(config$infected_file)
          terra::values(comp_year) <- data$infected[[q]]
          terra::values(reference) <- config$infection_years2[[q]]
          mask <- terra::rast(config$infected_file)
          terra::values(mask) <- config$mask_matrix
          ad <-
            quantity_allocation_disagreement(reference,
                                             comp_year,
                                             config$configuration,
                                             mask)
          if (file.exists(config$point_file)) {
            obs_data <- vect(config$point_file)
            obs_data <- terra::project(obs_data, comp_year)
            s <- extract(comp_year, obs_data)
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
              points_true_negative + points_true_positive + points_false_negative + points_false_positive
            ad$points_accuracy <-
              (points_true_negative + points_true_positive) / points_total_obs
            ad$points_precision <-
              points_true_positive / (points_true_positive + points_false_positive)
            ad$points_recall <-
              points_true_positive / (points_true_positive + points_false_negative)
            ad$points_specificiity <-
              points_true_negative / (points_true_negative + points_false_positive)

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
      write.csv(output_step, ffOut(paste("output_step_", j, ".csv", sep = "")))
    }
    if (j == 1) {
      cum_output_step <- output_step
      assign(paste("cum_output_step_", j, sep = ""), cum_output_step/j)
      output_list[[paste0("cum_output_step_", j)]] <- cum_output_step
      if (config$write_outputs %in% config$output_write_list) {
        write.csv(cum_output_step, ffOut(paste("cum_output_step_", j, ".csv", sep = "")))
      }
    }
    else {
      assign(paste("cum_output_step", sep = ""), (cum_output_step + output_step))
      assign(paste("cum_output_step_", j, sep = ""), cum_output_step/j)
      output_list[[paste0("cum_output_step_", j)]] <- cum_output_step
      if (config$write_outputs %in% config$output_write_list) {
        write.csv(cum_output_step, ffOut(paste("cum_output_step_", j, ".csv", sep = "")))
      }

    }
  }

  if (config$write_outputs %in% config$output_write_list) {
    save(output_list, file = ffOut("validation_outputs.rdata"))
  }

  return(output_list)
}
