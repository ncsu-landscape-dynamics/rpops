#' @title PoPS (Pest or Pathogen Spread) model
#'
#' @description A dynamic species distribution model for pest or pathogen spread
#' in forest or agricultural ecosystems. The model is process based meaning that
#' it uses understanding of the effect of weather and other environmental
#' factors on reproduction and survival of the pest/pathogen in order to
#' forecast spread of the pest/pathogen into the future. This function performs
#' a single stochastic realization of the model and is predominantly used for
#' automated tests of model features.
#'
#' @param infected_file Raster file with initial infections. Units for
#' infections are based on data availability and the way the units used for your
#' host file is created (e.g. percent area, # of hosts per cell, etc.).
#' @param host_file path to raster files with number of hosts and standard
#' deviation on those estimates can be based in 3 formats (a single file with
#' number of hosts, a single file with 2 layers number of hosts and standard
#' deviation, or two files 1 with number of hosts and the other with standard
#' deviation of those estimates). The units for this can be of many formats the
#' two most common that we use are either percent area (0 to 100) or # of hosts
#' in the cell. Usually depends on data available and estimation methods.
#' @param total_populations_file path to raster file with number of total
#' populations of all hosts and non-hosts. This depends on how your host data is
#' set up. If host is percent area then this should be a raster with values
#' that are 100 anywhere with host. If host file is # of hosts in a cell then
#' this should be a raster with values that are the max of the host raster any
#' where the # of hosts is greater than 0.
#' @param temp boolean that allows the use of temperature coefficients to
#' modify spread (TRUE or FALSE)
#' @param temperature_coefficient_file path to raster file with temperature
#' coefficient data for the timestep and and time period specified (e.g. if timestep
#' = week and start_date = 2017_01_01 and end_date = 2019_12_31 this file would
#' have 52 * 3 bands = 156 bands with data being weekly precipitation
#' coefficients). We convert raw precipitation values to coefficients that
#' affect the reproduction and survival of the pest all values in the raster are
#' between 0 and 1.
#' @param precip boolean that allows the use of precipitation coefficients to
#' modify spread (TRUE or FALSE)
#' @param precipitation_coefficient_file Raster file with precipitation
#' coefficient data for the timestep and time period specified (e.g. if timestep
#' = week and start_date = 2017_01_01 and end_date = 2019_12_31 this file would
#' have 52 * 3 bands = 156 bands with data being weekly precipitation
#' coefficients). We convert raw precipitation values to coefficients that
#' affect the reproduction and survival of the pest all values in the raster are
#' between 0 and 1.
#' @param time_step how often should spread occur options: ('day', 'week', 'month').
#' @param season_month_start when does spread first start occurring in the year
#' for your pest or pathogen (integer value between 1 and 12)
#' @param season_month_end when does spread end during the year for your pest
#' or pathogen (integer value between 1 and 12)
#' @param start_date date to start the simulation with format ('YYYY_MM_DD')
#' @param end_date date to end the simulation with format ('YYYY_MM_DD')
#' @param use_lethal_temperature a boolean to answer the question: does your
#' pest or pathogen have a temperature at which it cannot survive?
#' (TRUE or FALSE)
#' @param temperature_file path to raster file with temperature data for
#' minimum temperature
#' @param lethal_temperature the temperature in degrees C at which lethal
#' temperature related mortality occurs for your pest or pathogen (-50 to 60)
#' @param lethal_temperature_month the month in which lethal temperature
#' related mortality occurs for your pest or pathogen integer value between 1
#' and 12
#' @param mortality_on  boolean to turn host mortality on and off
#' (TRUE or FALSE)
#' @param mortality_rate rate at which mortality occurs value between 0 and 1
#' @param mortality_time_lag time lag from infection until mortality can occur
#' in time steps integer >= 1
#' @param mortality_frequency sets the frequency of mortality calculations occur
#' either ('year', 'month', week', 'day', 'time step', or 'every_n_steps')
#' @param mortality_frequency_n sets number of units from mortality_frequency in
#' which to run the mortality calculation if mortality_frequency is
#' 'every_n_steps'. Must be an integer >= 1.
#' @param management boolean to allow use of management (TRUE or FALSE)
#' @param treatment_dates dates in which to apply treatment list with format
#' ('YYYY_MM_DD') (needs to be the same length as treatment_file and
#' pesticide_duration)
#' @param treatments_file path to raster files with treatment data by dates.
#' Needs to be a list of files the same length as treatment_dates and
#' pesticide_duration.
#' @param treatment_method what method to use when applying treatment one of
#' ("ratio" or "all infected"). ratio removes a portion of all infected and
#' susceptibles, all infected removes all infected a portion of susceptibles.
#' @param natural_kernel_type what type of dispersal kernel should be used for
#' natural dispersal. Current dispersal kernel options are
#'  ('cauchy', 'exponential', 'uniform', 'deterministic neighbor','power law',
#'  'hyperbolic secant', 'gamma', 'weibull',
#'  'logistic')
#' @param anthropogenic_kernel_type what type of dispersal kernel should be
#' used for anthropogenic dispersal. Current dispersal kernel options are
#'  ('cauchy', 'exponential', 'uniform', 'deterministic neighbor','power law',
#'  'hyperbolic secant', 'gamma', 'weibull', 'logistic', 'network')
#' @param natural_dir sets the predominate direction of natural dispersal
#' usually due to wind values ('N', 'NW', 'W', 'SW', 'S', 'SE', 'E', 'NE',
#' 'NONE')
#' @param anthropogenic_dir sets the predominate direction of anthropogenic
#' dispersal usually due to human movement typically over long distances (e.g.
#' nursery trade, movement of firewood, etc..) ('N', 'NW', 'W', 'SW', 'S',
#' 'SE', 'E', 'NE', 'NONE')
#' @param pesticide_duration how long does the pestcide (herbicide, vaccine,
#' etc..) last before the host is susceptible again. If value is 0 treatment
#' is a culling (i.e. host removal) not a pesticide treatment. (needs to be the
#'  same length as treatment_dates and treatment_file)
#' @param pesticide_efficacy how effictive is the pesticide at preventing the
#' disease or killing the pest (if this is 0.70 then when applied it
#' successfully treats 70 percent of the plants or animals)
#' @param random_seed sets the random seed for the simulation used for reproducibility
#' @param output_frequency sets when outputs occur either ('year', 'month',
#' 'week', 'day', 'time step', or 'every_n_steps')
#' @param output_frequency_n sets number of units from output_frequency in which
#' to export model results if mortality_frequency is 'every_n_steps'.
#' Must be an integer >= 1.
#' @param movements_file this is a csv file with columns lon_from, lat_from,
#' lon_to, lat_to, number of animals, and date.
#' @param use_movements this is a boolean to turn on use of the movement module.
#' @param latency_period How many times steps does it take to for exposed
#' populations become infected/infested. This is an integer value and must be
#' greater than 0 if model type is SEI.
#' @param model_type What type of model most represents your system. Options
#' are "SEI" (Susceptible - Exposed - Infected/Infested) or "SI"
#' (Susceptible - Infected/Infested). Default value is "SI".
#' @param parameter_means A vector of the means of the model parameters
#' (reproductive_rate, natural_dispersal_distance, percent_natural_dispersal,
#' anthropogenic_dispersal_distance, natural kappa, anthropogenic kappa,
#' network_min_distance, and network_max_distance). 1x8 vector.
#' @param parameter_cov_matrix A covariance matrix from the previous years
#' posterior parameter estimation ordered from (reproductive_rate,
#' natural_dispersal_distance, percent_natural_dispersal,
#' anthropogenic_dispersal_distance, natural kappa, anthropogenic kappa,
#' network_min_distance, and network_max_distance) Should be 8x8 matrix.
#' @param start_exposed Do your initial conditions start as exposed or infected
#' (only used if model_type is "SEI"). Default False. If this is TRUE need to
#' have both an infected_file (this can be a raster of all 0's) and exposed_file
#' @param generate_stochasticity Boolean to indicate whether to use
#' stochasticity in reproductive functions default is TRUE
#' @param establishment_stochasticity Boolean to indicate whether to use
#' stochasticity in establishment functions default is TRUE
#' @param movement_stochasticity Boolean to indicate whether to use
#' stochasticity in movement functions default is TRUE
#' @param dispersal_stochasticity Boolean to indicate whether to use a stochasticity in the
#' dispersal kernel default is TRUE
#' @param establishment_probability Threshold to determine establishment if
#' establishment_stochasticity is FALSE (range 0 to 1, default = 0.5)
#' @param dispersal_percentage  Percentage of dispersal used to calculate the
#' bounding box for deterministic dispersal
#' @param quarantine_areas_file path to raster file with quarantine boundaries
#' used in calculating likelihood of quarantine escape if use_quarantine is
#' TRUE
#' @param use_quarantine boolean to indicate whether or not there is a
#' quarantine area if TRUE must pass in a raster file indicating the
#' quarantine areas (default = FALSE)
#' @param quarantine_directions string with comma separated directions to include
#' in the quarantine direction analysis, e.g., 'N,E'. By default all directions
#' (N, S, E, W) are considered
#' @param use_spreadrates boolean to indicate whether or not to calculate
#' spread rates
#' @param use_overpopulation_movements boolean to indicate whether to use
#' the overpopulation pest movement module (driven by the natural kernel with
#' its scale parameter modified by a coefficient)
#' @param overpopulation_percentage percentage of occupied hosts when the cell
#' is considered to be overpopulated
#' @param leaving_percentage percentage of pests leaving an overpopulated cell
#' @param leaving_scale_coefficient coefficient to multiply scale parameter of
#' the natural kernel (if applicable)
#' @param exposed_file a file with the exposed for the current
#' @param mask Raster file used to provide a mask to remove 0's that are not
#' true negatives from comparisons (e.g. mask out lakes and oceans from statics
#' if modeling terrestrial species). This can also be used to mask out areas
#' that can't be managed in the auto_manage function.
#' @param network_filename entire file path for the network file. Used if
#' anthropogenic_kernel_type = 'network'.
#' @param use_survival_rates boolean to indicate if the model will use survival rates to limit the
#' survival or emergence of overwintering generations.
#' @param survival_rate_month What month do over wintering generations emerge. We suggest using the
#' month before for this parameter as it is when the survival rates raster will be applied.
#' @param survival_rate_day what day should the survival rates be applied
#' @param survival_rates_file Raster file with survival rates from 0 to 1 representing the
#' percentage of emergence for a cell.
#' @param network_movement What movement type do you want to use in the network kernel either
#' "walk", "jump", or "teleport". "walk" allows dispersing units to leave the network at any cell
#' along the edge. "jump" automatically moves to the nearest node when moving through the network.
#' "teleport" moves from node to node most likely used for airport and seaport networks.
#' @param use_initial_condition_uncertainty boolean to indicate whether or not to propagate and
#' partition uncertainty from initial conditions. If TRUE the infected_file needs to have 2 layers
#' one with the mean value and one with the standard deviation. If an SEI model is used the
#' exposed_file needs to have 2 layers one with the mean value and one with the standard
#' deviation
#' @param use_host_uncertainty boolean to indicate whether or not to propagate and partition
#' uncertainty from host data. If TRUE the host_file needs to have 2 layers one with the mean value
#' and one with the standard deviation.
#' @param weather_type string indicating how the weather data is passed in  either
#' as a mean and standard deviation to represent uncertainty ("probabilistic") or as a time
#' series ("deterministic")
#' @param dispersers_to_soils_percentage range from 0 to 1 representing the percentage
#' of dispersers that fall to the soil and survive.
#' @param soil_starting_pest_file path to the raster file with the starting
#' amount of pest or pathogen.
#' @param multiple_random_seeds boolean to indicate if the model should use multiple random seeds
#' (allows for performing uncertainty partitioning) or a single random seed (backwards
#' compatibility option). Default is FALSE.
#' @param random_seeds A file path to the file with the .csv file containing random_seeds table. 
#' Use a file if you are trying to recreate an exact analysis otherwise we suggest leaving the 
#' default. Default is Null which draws the seed numbers for each. 
#'
#' @useDynLib PoPS, .registration = TRUE
#' @importFrom terra app rast xres yres classify extract ext as.points ncol nrow
#' nlyr rowFromCell colFromCell values as.matrix rowFromCell colFromCell crs
#' @importFrom Rcpp sourceCpp evalCpp
#' @importFrom  stats runif
#' @importFrom lubridate interval time_length mdy %within%
#' @importFrom utils read.csv
#' @importFrom sp SpatialPointsDataFrame CRS spTransform
#' @importFrom  methods is
#' @return list of infected and susceptible per year
#' @export
#'

pops <- function(infected_file,
                 host_file,
                 total_populations_file,
                 parameter_means,
                 parameter_cov_matrix,
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
                 pesticide_duration = c(0),
                 pesticide_efficacy = 1.0,
                 random_seed = NULL,
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
                 exposed_file = "",
                 mask = NULL,
                 network_filename = "",
                 network_movement = "walk",
                 use_initial_condition_uncertainty = FALSE,
                 use_host_uncertainty = FALSE,
                 weather_type = "deterministic",
                 dispersers_to_soils_percentage = 0,
                 quarantine_directions = "",
                 multiple_random_seeds = FALSE,
                 random_seeds = NULL) {

  config <- c()
  config$random_seed <- random_seed
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
  config$use_survival_rates <- use_survival_rates
  config$survival_rate_month <- survival_rate_month
  config$survival_rate_day <- survival_rate_day
  config$survival_rates_file <- survival_rates_file
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
  config$use_quarantine <- use_quarantine
  config$quarantine_directions <- quarantine_directions
  config$use_spreadrates <- use_spreadrates
  config$use_overpopulation_movements <- use_overpopulation_movements
  config$overpopulation_percentage <- overpopulation_percentage
  config$leaving_percentage <- leaving_percentage
  config$leaving_scale_coefficient <- leaving_scale_coefficient
  # added number of iterations to config to avoid multiple if else statements
  # in configuration function used to determine number of draws from parameter
  # distribution
  config$number_of_iterations <- 2
  # add function name for use in configuration function to skip
  # function specific specific configurations namely for validation and
  # calibration.
  config$function_name <- "pops"
  config$failure <- NULL
  config$exposed_file <- exposed_file
  config$write_outputs <- "None"
  config$output_folder_path <- ""
  config$mortality_frequency <- mortality_frequency
  config$mortality_frequency_n <- mortality_frequency_n

  config$network_filename <- network_filename
  config$network_movement <- network_movement
  config$use_initial_condition_uncertainty <- use_initial_condition_uncertainty
  config$use_host_uncertainty <- use_host_uncertainty
  config$weather_type <- weather_type
  config$dispersers_to_soils_percentage <- dispersers_to_soils_percentage
  config$multiple_random_seeds <- multiple_random_seeds
  config$random_seeds <-random_seeds

  config <- configuration(config)

  if (!is.null(config$failure)) {
    stop(config$failure)
  }

  config <- draw_parameters(config) # draws parameter set for the run

  if (config$use_initial_condition_uncertainty) {
    config$infected <-  matrix_norm_distribution(config$infected_mean, config$infected_sd)
    exposed2 <- matrix_norm_distribution(config$exposed_mean, config$exposed_sd)
    exposed <- config$exposed
    exposed[[config$latency_period + 1]] <- exposed2
    config$exposed <- exposed
  } else {
    config$infected <- config$infected_mean
    exposed2 <- config$exposed_mean
    exposed <- config$exposed
    exposed[[config$latency_period + 1]] <- exposed2
    config$exposed <- exposed
  }

  if (config$use_host_uncertainty) {
    config$host <- matrix_norm_distribution(config$host_mean, config$host_sd)
  } else {
    config$host <- config$host_mean
  }

  susceptible <- config$host - config$infected - exposed2
  susceptible[susceptible < 0] <- 0

  config$susceptible <- susceptible
  config$total_hosts <- config$host
  config$total_exposed <- exposed2

  if (config$mortality_on) {
    mortality_tracker2 <- config$mortality_tracker
    mortality_tracker2[[length(mortality_tracker2)]] <- config$infected
    config$mortality_tracker <- mortality_tracker2
  }

  data <- pops_model(random_seed = config$random_seed[1],
                     multiple_random_seeds = config$multiple_random_seeds,
                     random_seeds = config$random_seeds,
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
                     res = config$res,
                     rows_cols = config$rows_cols,
                     time_step = config$time_step,
                     reproductive_rate = config$reproductive_rate[1],
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
                     percent_natural_dispersal = config$percent_natural_dispersal[1],
                     natural_distance_scale = config$natural_distance_scale[1],
                     anthropogenic_distance_scale = config$anthropogenic_distance_scale[1],
                     natural_dir = config$natural_dir,
                     natural_kappa = config$natural_kappa[1],
                     anthropogenic_dir = config$anthropogenic_dir,
                     anthropogenic_kappa = config$anthropogenic_kappa[1],
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
                     network_min_distance = config$network_min_distance[1],
                     network_max_distance = config$network_max_distance[1],
                     network_filename = config$network_filename,
                     network_movement = config$network_movement,
                     weather_size = config$weather_size,
                     weather_type = config$weather_type,
                     dispersers_to_soils_percentage = config$dispersers_to_soils_percentage
  )

  return(data)
}
