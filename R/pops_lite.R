#' PoPS (Pest or Pathogen Spread) Lite Model Multiple Runs
#' 
#' A process-based model designed to forecast the spread of pests or pathogens 
#' in forest or agricultural ecosystems. It incorporates the effects of weather 
#' and environmental factors on reproduction and survival. The model runs 
#' multiple stochastic simulations, propagating uncertainty in parameters, 
#' initial conditions, and drivers, similar to `pops_multirun`. However, this 
#' function exports raw outputs directly from `pops_model`. 
#'
#' Use this function if raw simulation data is needed or if extent and/or 
#' resolution constraints prevent using `pops_multirun`.

#'
#' @inheritParams pops
#' 
#' @param config_file Path to config file produced when calling `configuration`. 
#' The config file includes all data necessary used to set up c++ PoPS model
#' 
#' @importFrom terra app rast xres yres classify extract ext as.points ncol nrow project
#' nlyr rowFromCell colFromCell values as.matrix rowFromCell colFromCell crs vect
#' @importFrom stats runif rnorm median sd
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar% %do%
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom lubridate interval time_length mdy %within%
#' @importFrom utils write.csv read.table read.csv
#' @importFrom methods is
#'
#' @return list of infected and susceptible per year
#' @export

pops_lite <- function(config_file = "") {
  config <- readRDS(config_file)
  
  if (!is.null(config$failure)) {
    stop(config$failure)
  }
  
  uid <- generate_unique_id()
  
  if (!dir.exists(config$output_folder_path)) {
    suppressWarnings(dir.create(config$output_folder_path, recursive = TRUE))
  }
  
  
  if (config$multiple_random_seeds &&
      is.null(config$file_random_seeds) &&
      dir.exists(config$output_folder_path)) {
    write.csv(
      config$random_seeds,
      paste0(
        config$output_folder_path,
        uid,
        "_forecast_random_seeds.csv"
      ),
      row.names = FALSE
    )
  }
  
  config$crs <- terra::crs(terra::rast(config$host_file))
  
  i <- NULL
  cl <- parallel::makeCluster(config$core_count)
  doParallel::registerDoParallel(cl)
  foreach::foreach(
    i = seq_len(config$number_of_iterations),
    .packages = c("PoPS", "terra")
  ) %dopar% {
    config <- update_config(config)
    
    data <- PoPS::pops_model(
      random_seed = config$random_seed[1],
      multiple_random_seeds = config$multiple_random_seeds,
      random_seeds = as.matrix(config$random_seeds[i, ])[1, ],
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
      quarantine_directions = config$quarantine_directions,
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
      weather_coefficient_sd = config$weather_coefficient_sd,
      res = config$res,
      rows_cols = config$rows_cols,
      time_step = config$time_step,
      reproductive_rate = config$reproductive_rate,
      spatial_indices = config$spatial_indices,
      season_month_start_end = config$season_month_start_end,
      soil_reservoirs = config$soil_reservoirs,
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
      network_movement = config$network_movement,
      weather_size = config$weather_size,
      weather_type = config$weather_type,
      dispersers_to_soils_percentage = config$dispersers_to_soils_percentage,
      use_soils = config$use_soils
    )
    
    data[c("spatial_indices",
           "soil_reservoirs",
           "total_populations",
           "total_exposed")] <- NULL
    
    if (config$model_type == "SI") {
      data[c("exposed", "total_exposed", "resistant")] <- NULL
    }
    
    if (!config$mortality_on) {
      data$mortality <- NULL
    }
    
    if (!config$use_quarantine) {
      data[c(
        "quarantine_escape",
        "quarantine_escape_directions",
        "quarantine_escape_distance"
      )] <- NULL
    }
    
    # Remove any null in data
    data <- data[!sapply(data, is.null)]
    gc()
    
    saveRDS(
      data,
      file = paste0(config$output_folder_path, uid, "_pops_output_", i, ".rds"),
      compress = TRUE
    )
    rm(data)
  }
  stopCluster(cl)
  return(cat("Raw PoPS runs outputs saved to output_path"))
}
