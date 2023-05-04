#' @title PoPS (Pest or Pathogen Spread) model
#'
#' @description Wrapper for pops_model_cpp. A dynamic species distribution model
#' for pest or pathogen spread in forest or agricultural ecosystems. The model
#' is process based meaning that it uses understanding of the effect of weather
#' and other environmental factors on reproduction and survival of the
#' pest/pathogen in order to forecast spread of the pest/pathogen into the
#' future.
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
#' @param rows_cols vector of number of rows and columns in the raster files
#' @param season_month_start_end vector of months when spread starts and stops
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
#' calculated either ('year', 'month', 'week', 'day' or 'time step')
#' (default is year)
#' @param quarantine_frequency_n sets number units ('year', 'month', 'week',
#' 'day' or 'time step') in which to calculate and export quarantine statistics.
#' @param spreadrate_frequency sets how often the spread rate statistics are
#' calculated either ('year', 'month', 'week', 'day' or 'time step')
#' (default is year)
#' @param spreadrate_frequency_n sets number units ('year', 'month', 'week',
#' 'day' or 'time step') in which to calculate and export spread rate
#' statistics.
#' @param spatial_indices list of all spatial locations with suitable hosts
#' @param bbox bounding box for network kernel
#' @param network_min_distance minimum time a propagule rides on the network. Used if
#' anthropogenic_kernel_type = 'network'.
#' @param network_max_distance maximum time a propagule rides on the network. Used if
#' anthropogenic_kernel_type = 'network'.
#' @param survival_rates vector of matrices of survival rates used to determine percentage of
#' overwinter population that emerges
#' @param weather_size the number of matrices in a list or layers in a raster object
#' @param weather_type string indicating how the weather data is passed in  either
#' as a mean and standard deviation to represent uncertainty ("probablisticc") or as a time
#' series ("deterministic")
#' @param dispersers_to_soils_percentage range from 0 to 1 representing the percentage
#' of dispersers that fall to the soil and survive.
#'
#' @return list of vector matrices of infected and susceptible hosts per
#' simulated year and associated statistics (e.g. spread rate)
#' @export
#'

pops_model <-
  function(random_seed,
           use_lethal_temperature,
           lethal_temperature,
           lethal_temperature_month,
           use_survival_rates,
           survival_rate_month,
           survival_rate_day,
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
           survival_rates,
           weather_coefficient,
           res,
           rows_cols,
           time_step,
           reproductive_rate,
           spatial_indices,
           season_month_start_end,
           mortality_rate = 0.0,
           mortality_time_lag = 2,
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
           dispersal_stochasticity = TRUE,
           establishment_probability = 0,
           dispersal_percentage = 0.99,
           use_overpopulation_movements = FALSE,
           overpopulation_percentage = 0.0,
           leaving_percentage = 0.0,
           leaving_scale_coefficient = 1.0,
           bbox = NULL,
           network_min_distance = 0,
           network_max_distance = 0,
           network_filename = "",
           network_movement = "walk",
           weather_size = 0,
           weather_type = "deterministic",
           dispersers_to_soils_percentage = 0) {

    # List of overpopulation parameters of type double
    overpopulation_config <- c()
    overpopulation_config$overpopulation_percentage <- overpopulation_percentage
    overpopulation_config$leaving_percentage <- leaving_percentage
    overpopulation_config$leaving_scale_coefficient <- leaving_scale_coefficient


    # List of frequency n parameters
    frequencies_n_config <- c()
    frequencies_n_config$output_frequency_n <- output_frequency_n
    frequencies_n_config$quarantine_frequency_n <- quarantine_frequency_n
    frequencies_n_config$spreadrate_frequency_n <- spreadrate_frequency_n
    frequencies_n_config$mortality_frequency_n <- mortality_frequency_n

    # Network configuration
    network_config <- NULL;
    network_data_config <- NULL;
    if (!(is.na(network_filename) || is.null(network_filename) || network_filename == '')) {
      network_config <- c()
      network_config$network_min_distance <- network_min_distance
      network_config$network_max_distance <- network_max_distance
      network_config$network_movement <- network_movement

      network_data_config <- c()
      network_data_config$network_filename <- network_filename
    }

    # List of frequencies type string
    frequency_config <- c()
    frequency_config$time_step <- time_step
    frequency_config$mortality_frequency <- mortality_frequency
    frequency_config$spreadrate_frequency <- spreadrate_frequency
    frequency_config$quarantine_frequency <- quarantine_frequency
    frequency_config$output_frequency <- output_frequency

    # List of all booleans
    bool_config <- c()
    bool_config$use_lethal_temperature <- use_lethal_temperature
    bool_config$mortality_on <- mortality_on
    bool_config$use_movements <- use_movements
    bool_config$weather <- weather
    bool_config$use_anthropogenic_kernel <- use_anthropogenic_kernel
    bool_config$use_quarantine <- use_quarantine
    bool_config$use_spreadrates <- use_spreadrates
    bool_config$generate_stochasticity <- generate_stochasticity
    bool_config$establishment_stochasticity <- establishment_stochasticity
    bool_config$movement_stochasticity <- movement_stochasticity
    bool_config$dispersal_stochasticity <- dispersal_stochasticity
    bool_config$use_overpopulation_movements <- use_overpopulation_movements
    bool_config$use_survival_rate <- use_survival_rates


    data <-
      pops_model_cpp(random_seed = random_seed,
                     lethal_temperature = lethal_temperature,
                     lethal_temperature_month = lethal_temperature_month,
                     infected = infected,
                     total_exposed = total_exposed,
                     exposed = exposed,
                     susceptible = susceptible,
                     total_populations = total_populations,
                     total_hosts = total_hosts,
                     mortality_tracker = mortality_tracker,
                     mortality = mortality,
                     quarantine_areas = quarantine_areas,
                     treatment_maps = treatment_maps,
                     treatment_dates = treatment_dates,
                     pesticide_duration = pesticide_duration,
                     resistant = resistant,
                     movements = movements,
                     movements_dates = movements_dates,
                     temperature = temperature,
                     survival_rates = survival_rates,
                     weather_coefficient = weather_coefficient,
                     bbox = bbox,
                     res = res,
                     rows_cols = rows_cols,
                     reproductive_rate = reproductive_rate,
                     spatial_indices = spatial_indices,
                     season_month_start_end = season_month_start_end,
                     frequency_config = frequency_config,
                     bool_config = bool_config,
                     mortality_rate = mortality_rate,
                     mortality_time_lag = mortality_time_lag,
                     start_date = start_date,
                     end_date = end_date,
                     treatment_method = treatment_method,
                     natural_kernel_type = natural_kernel_type,
                     anthropogenic_kernel_type = anthropogenic_kernel_type,
                     percent_natural_dispersal = percent_natural_dispersal,
                     natural_distance_scale = natural_distance_scale,
                     anthropogenic_distance_scale = anthropogenic_distance_scale,
                     natural_dir = natural_dir,
                     natural_kappa = natural_kappa,
                     anthropogenic_dir = anthropogenic_dir,
                     anthropogenic_kappa = anthropogenic_kappa,
                     frequencies_n_config = frequencies_n_config,
                     model_type_ = model_type_,
                     latency_period = latency_period,
                     dispersal_percentage = dispersal_percentage,
                     survival_rate_month = survival_rate_month,
                     survival_rate_day = survival_rate_day,
                     establishment_probability = establishment_probability,
                     overpopulation_config = overpopulation_config,
                     network_config = network_config,
                     network_data_config = network_data_config,
                     weather_size = weather_size,
                     weather_type = weather_type,
                     dispersers_to_soils_percentage = dispersers_to_soils_percentag
    )

  }
