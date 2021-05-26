#' @title PoPS (Pest or Pathogen Spread) model
#'
#' @description Wrapper for pops_model_cpp. A dynamic species distribution model for pest or
#' pathogen spread in forest or agricultural ecosystems. The model is process
#' based meaning that it uses understanding of the effect of weather and other
#' environmental factors on reproduction and survival of the pest/pathogen in
#' order to forecast spread of the pest/pathogen into the future.
#'
#' @inheritParams pops
#' @param weather Boolean that is true if weather is used
#' @param infected matrix of infected hosts
#' @param susceptible matrix of susceptible hosts
#' @param mortality_tracker matrix of 0's to track mortality per year
#' @param mortality matrix to track cumulative mortality
#' @param resistant matrix to track resistant population over time
#' @param total_populations  matrix of total populations
#' @param total_hosts matrix of all hosts
#' @param treatment_maps list of matrices where treatment or management has
#' occurred in a given year
#' @param temperature vector of matrices of temperature values used to check
#' against lethal temperature
#' @param weather_coefficient vector of matrices of weather coefficients
#' @param res  vector of east/west resolution and north/south resolution
#' @param num_rows number of rows in the raster files
#' @param num_cols number of columns in the raster files
#' @param use_anthropogenic_kernel A boolean that turns on the use of the
#' anthropogenic or long distance dispersal portion of the 2 scale dispersal
#' kernel function
#' @param movements a matrix with columns lon_from, lat_from, lon_to, lat_to,
#' number of animals, and date.
#' @param use_movements this is a boolean to turn on use of the movement module.
#' @param movements_dates this is a list of dates passed as strings in the
#' format 'YYYY-MM-DD'
#' @param exposed vector of matrices of the exposed class for use with "SEI"
#' model type
#' @param total_exposed sum of all exposed cohorts in exposed class for use with
#' "SEI" model type
#' @param model_type_ What type of model most represents your system. Options
#' are "SEI" (Susceptible - Exposed - Infected/Infested) or "SI"
#' (Susceptible - Infected/Infested). Default value is "SI".
#' @param reproductive_rate number of spores or pest units produced by a single
#'  host under optimal weather conditions
#' @param percent_natural_dispersal  what percentage of dispersal is natural
#' range versus anthropogenic range value between 0 and 1
#' @param natural_distance_scale distance scale parameter for natural range
#' dispersal kernel numeric value > 0
#' @param anthropogenic_distance_scale distance scale parameter for
#' anthropogenic range dispersal kernel numeric value > 0
#' @param natural_kappa sets the strength of the natural direction in the
#' von-mises distribution numeric value between 0.01 and 12
#' @param anthropogenic_kappa sets the strength of the anthropogenic direction
#' in the von-mises distribution numeric value between 0.01 and 12
#' @param quarantine_areas areas that are set as quarantined for computing
#' escape from quarantine statistics.
#' @param quarantine_frequency sets how often the quarantine statistics are
#' calculated either ('year', 'month', 'week', 'day' or 'time step') (default is year)
#' @param quarantine_frequency_n sets number units ('year', 'month', 'week',
#' 'day' or 'time step') in which to calculate and export quarantine statistics.
#' @param spreadrate_frequency sets how often the spread rate statistics are
#' calculated either ('year', 'month', 'week', 'day' or 'time step')
#' (default is year)
#' @param spreadrate_frequency_n sets number units ('year', 'month', 'week',
#' 'day' or 'time step') in which to calculate and export spread rate statistics.
#' @param spatial_indices list of all spatial locations with suitable hosts
#' @return list of vector matrices of infected and susceptible hosts per
#' simulated year and associated statistics (e.g. spread rate)
#' @export
#'

pops_model <-
  function(random_seed,
           use_lethal_temperature,
           lethal_temperature,
           lethal_temperature_month,
           infected,
           total_exposed,
           exposed,
           susceptible,
           total_populations,
           total_hosts,
           mortality_on,
           mortality_tracker,
           mortality,
           quarantine_areas,
           treatment_maps,
           treatment_dates,
           pesticide_duration,
           resistant,
           use_movements,
           movements,
           movements_dates,
           weather,
           temperature,
           weather_coefficient,
           res,
           num_rows,
           num_cols,
           time_step,
           reproductive_rate,
           spatial_indices,
           mortality_rate = 0.0,
           mortality_time_lag = 2,
           season_month_start = 1,
           season_month_end = 12,
           start_date = "2018-01-01",
           end_date = "2018-12-31",
           treatment_method = "ratio",
           natural_kernel_type = "cauchy",
           anthropogenic_kernel_type = "cauchy",
           use_anthropogenic_kernel = FALSE,
           percent_natural_dispersal = 0.0,
           natural_distance_scale = 21,
           anthropogenic_distance_scale = 0.0,
           natural_dir = "NONE",
           natural_kappa = 0,
           anthropogenic_dir = "NONE",
           anthropogenic_kappa = 0,
           output_frequency = "year",
           output_frequency_n = 1,
           quarantine_frequency = "year",
           quarantine_frequency_n = 1,
           use_quarantine = FALSE,
           spreadrate_frequency = "year",
           spreadrate_frequency_n = 1,
           mortality_frequency = "year",
           mortality_frequency_n = 1,
           use_spreadrates = FALSE,
           model_type_ = "SI",
           latency_period = 0,
           generate_stochasticity = TRUE,
           establishment_stochasticity = TRUE,
           movement_stochasticity = TRUE,
           deterministic = FALSE,
           establishment_probability = 0,
           dispersal_percentage = 0.99,
           use_overpopulation_movements = FALSE,
           overpopulation_percentage = 0.0,
           leaving_percentage = 0.0,
           leaving_scale_coefficient = 1.0) {

    # List of overpopulation parameters of type double
    overpopulation_config = c()
    overpopulation_config$overpopulation_percentage <- overpopulation_percentage
    overpopulation_config$leaving_percentage <- leaving_percentage
    overpopulation_config$leaving_scale_coefficient <- leaving_scale_coefficient

    data <-
      pops_model_cpp(random_seed = random_seed,
                     use_lethal_temperature = use_lethal_temperature,
                     lethal_temperature = lethal_temperature,
                     lethal_temperature_month = lethal_temperature_month,
                     infected = infected,
                     total_exposed = total_exposed,
                     exposed = exposed,
                     susceptible = susceptible,
                     total_populations = total_populations,
                     total_hosts = total_hosts,
                     mortality_on = mortality_on,
                     mortality_tracker = mortality_tracker,
                     mortality = mortality,
                     quarantine_areas = quarantine_areas,
                     treatment_maps = treatment_maps,
                     treatment_dates = treatment_dates,
                     pesticide_duration = pesticide_duration,
                     resistant = resistant,
                     use_movements = use_movements,
                     movements = movements,
                     movements_dates = movements_dates,
                     weather = weather,
                     temperature = temperature,
                     weather_coefficient = weather_coefficient,
                     res = res,
                     num_rows = num_rows,
                     num_cols = num_cols,
                     time_step = time_step,
                     reproductive_rate = reproductive_rate,
                     spatial_indices = spatial_indices,
                     mortality_rate = mortality_rate,
                     mortality_time_lag = mortality_time_lag,
                     season_month_start = season_month_start,
                     season_month_end = season_month_end,
                     start_date = start_date,
                     end_date = end_date,
                     treatment_method = treatment_method,
                     natural_kernel_type = natural_kernel_type,
                     anthropogenic_kernel_type =
                       anthropogenic_kernel_type,
                     use_anthropogenic_kernel =
                       use_anthropogenic_kernel,
                     percent_natural_dispersal =
                       percent_natural_dispersal,
                     natural_distance_scale = natural_distance_scale,
                     anthropogenic_distance_scale =
                       anthropogenic_distance_scale,
                     natural_dir = natural_dir,
                     natural_kappa = natural_kappa,
                     anthropogenic_dir = anthropogenic_dir,
                     anthropogenic_kappa = anthropogenic_kappa,
                     output_frequency = output_frequency,
                     output_frequency_n = output_frequency_n,
                     quarantine_frequency = quarantine_frequency,
                     quarantine_frequency_n = quarantine_frequency_n,
                     use_quarantine = use_quarantine,
                     spreadrate_frequency = spreadrate_frequency,
                     spreadrate_frequency_n = spreadrate_frequency_n,
                     mortality_frequency = mortality_frequency,
                     mortality_frequency_n = mortality_frequency_n,
                     use_spreadrates = use_spreadrates,
                     model_type_ = model_type_,
                     latency_period = latency_period,
                     generate_stochasticity = generate_stochasticity,
                     establishment_stochasticity =
                       establishment_stochasticity,
                     movement_stochasticity = movement_stochasticity,
                     deterministic = deterministic,
                     establishment_probability =
                       establishment_probability,
                     dispersal_percentage = dispersal_percentage,
                     use_overpopulation_movements = use_overpopulation_movements,
                     overpopulation_config = overpopulation_config
    )

  }
