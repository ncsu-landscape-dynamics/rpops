#' @title PoPS (Pest or Pathogen Spread) model
#'
#' @description Wrapper for pops_model_cpp only used internally in other functions.
#'
#' @param config a list of all of the necessary data, parameters, files, and paths formatted for
#' the C++ model created from configuration and modified slightly by the other functions for each
#' simulation.
#'
#' @return list of vector matrices of infected and susceptible hosts per
#' simulated year and associated statistics (e.g. spread rate)
#' @export
#'

pops_model <-
  function(config) {

    # List of overpopulation parameters of type double
    overpopulation_config <- c()
    overpopulation_config$overpopulation_percentage <- config$overpopulation_percentage
    overpopulation_config$leaving_percentage <- config$leaving_percentage
    overpopulation_config$leaving_scale_coefficient <- config$leaving_scale_coefficient


    # List of frequency n parameters
    frequencies_n_config <- c()
    frequencies_n_config$output_frequency_n <- config$output_frequency_n
    frequencies_n_config$quarantine_frequency_n <- config$quarantine_frequency_n
    frequencies_n_config$spreadrate_frequency_n <- config$spreadrate_frequency_n
    frequencies_n_config$mortality_frequency_n <- config$mortality_frequency_n

    # Network configuration
    network_config <- NULL
    network_data_config <- NULL
    if (config$anthropogenic_kernel_type == "network") {
      network_config <- c()
      network_config$network_min_distances <- config$network_min_distances
      network_config$network_max_distances <- config$network_max_distances
      network_config$network_movement_types <- config$network_movement_types
      network_config$network_weights <- config$network_weights

      network_data_config <- c()
      network_data_config$network_filenames <- network_files
    }

    # List of frequencies type string
    frequency_config <- c()
    frequency_config$time_step <- config$time_step
    frequency_config$mortality_frequency <- config$mortality_frequency
    frequency_config$spreadrate_frequency <- config$spreadrate_frequency
    frequency_config$quarantine_frequency <- config$quarantine_frequency
    frequency_config$output_frequency <- config$output_frequency

    # List of all booleans
    bool_config <- c()
    bool_config$use_lethal_temperature <- config$use_lethal_temperature
    bool_config$mortality_on <- config$mortality_on
    bool_config$use_movements <- config$use_movements
    bool_config$weather <- config$weather
    bool_config$use_anthropogenic_kernel <- config$use_anthropogenic_kernel
    bool_config$use_quarantine <- config$use_quarantine
    bool_config$use_spreadrates <- config$use_spreadrates
    bool_config$generate_stochasticity <- config$generate_stochasticity
    bool_config$establishment_stochasticity <- config$establishment_stochasticity
    bool_config$movement_stochasticity <- config$movement_stochasticity
    bool_config$dispersal_stochasticity <- config$dispersal_stochasticity
    bool_config$use_overpopulation_movements <- config$use_overpopulation_movements
    bool_config$use_survival_rate <- config$use_overwinter_survival
    bool_config$use_soils <- config$use_soils


    data <-
      suppressWarnings(pops_model_cpp(random_seed = config$random_seed,
                     multiple_random_seeds = config$use_multiple_random_seeds,
                     random_seeds = config$random_seeds,
                     lethal_temperature = config$lethal_temperature,
                     lethal_temperature_month = config$lethal_temperature_month,
                     host_pools = config$host_pools,
                     total_populations = config$total_populations,
                     competency_table = config$competency_table_list,
                     pest_host_table = config$pest_host_table_list,
                     quarantine_areas = config$quarantine_areas,
                     quarantine_directions = config$quarantine_directions,
                     treatment_maps = config$treatment_maps,
                     treatment_dates = config$treatment_dates,
                     pesticide_duration = config$pesticide_durations,
                     movements = config$movements,
                     movements_dates = config$movements_dates,
                     temperature = config$temperature,
                     survival_rates = config$survival_rates,
                     weather_coefficient = config$weather_coefficient,
                     weather_coefficient_sd = config$weather_coefficient_sd,
                     bbox = config$bounding_box,
                     res = config$res,
                     rows_cols = config$rows_cols,
                     soil_reservoirs = config$soil_reservoirs,
                     reproductive_rate = config$reproductive_rate,
                     spatial_indices = config$spatial_indices,
                     season_month_start_end = config$season_month_start_end,
                     frequency_config = frequency_config,
                     bool_config = bool_config,
                     start_date = config$start_date,
                     end_date = config$end_date,
                     treatment_method = config$treatment_method,
                     natural_kernel_type = config$natural_kernel_type,
                     anthropogenic_kernel_type = config$anthropogenic_kernel_type,
                     percent_natural_dispersal = config$percent_natural_dispersal,
                     natural_distance_scale = config$natural_distance_scale,
                     anthropogenic_distance_scale = config$anthropogenic_distance_scale,
                     natural_dir = config$natural_dir,
                     natural_kappa = config$natural_kappa,
                     anthropogenic_dir = config$anthropogenic_dir,
                     anthropogenic_kappa = config$anthropogenic_kappa,
                     frequencies_n_config = frequencies_n_config,
                     model_type_ = config$model_type,
                     latency_period = config$latency_period,
                     establishment_probability = config$establishment_probability,
                     dispersal_percentage = config$dispersal_percentage,
                     survival_rate_month = config$overwinter_survival_rate_month,
                     survival_rate_day = config$overwinter_survival_rate_day,
                     overpopulation_config = overpopulation_config,
                     network_config = network_config,
                     network_data_config = network_data_config,
                     weather_size = config$weather_size,
                     weather_type = config$weather_type,
                     dispersers_to_soils_percentage = config$dispersers_to_soils_percentage
    ))

    return(data)
  }
