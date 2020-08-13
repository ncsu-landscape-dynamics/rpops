#' PoPS (Pest or Pathogen Spread) model
#' 
#' A dynamic species distribution model for pest or pathogen spread in forest or agricultural ecosystems. The model is process based
#' meaning that it uses understanding of the effect of weather and other environmental factors on reproduction and survival of the pest/pathogen in order to forecast
#' spread of the pest/pathogen into the future. 
#'
#' @param infected_file path to raster file with initial infections 
#' @param host_file path to raster files with number of hosts and standard deviation on those estimates can be based in 3 formats (a single file with number of hosts, a single file with 2 layers number of hosts and standard deviation, or two files 1 with number of hosts and the other with standard deviation of those estimates)
#' @param total_plants_file path to raster file with number of total plants
#' @param temp boolean that allows the use of temperature coefficients to modify spread (TRUE or FALSE)
#' @param temperature_coefficient_file path to raster file with temperature coefficient data for the timestep and number of years specified
#' @param precip boolean that allows the use of precipitation coefficients to modify spread (TRUE or FALSE)
#' @param precipitation_coefficient_file path to raster file with precipitation coefficient data for the timestep and number of years specified
#' @param time_step how often should spread occur options: ('day', 'week', 'month')
#' @param season_month_start when does spread first start occurring in the year for your pest or pathogen (integer value between 1 and 12)
#' @param season_month_end when does spread end during the year for your pest or pathogen (integer value between 1 and 12)
#' @param start_date date to start the simulation with format ('YYYY_MM_DD')
#' @param end_date date to end the simulation with format ('YYYY_MM_DD')
#' @param use_lethal_temperature a boolean to answer the question: does your pest or pathogen have a temperature at which it cannot survive? (TRUE or FALSE)
#' @param temperature_file path to raster file with temperature data for minimum temperature
#' @param lethal_temperature the temperature in degrees C at which lethal temperature related mortality occurs for your pest or pathogen (-50 and 60)
#' @param lethal_temperature_month the month in which lethal temperature related mortality occurs for your pest or pathogen integer value between 1 and 12
#' @param mortality_on  boolean to turn host mortality on and off (TRUE or FALSE)
#' @param mortality_rate rate at which mortality occurs value between 0 and 1 
#' @param mortality_time_lag time lag from infection until mortality can occur in years integer >= 1
#' @param management boolean to allow use of managemnet (TRUE or FALSE)
#' @param treatment_dates dates in which to apply treatment list with format ('YYYY_MM_DD') (needs to be the same length as treatment_file and pesticide_duration)
#' @param treatments_file path to raster file with treatment data by dates (needs to be the same length as treatment_dates and pesticide_duration)
#' @param treatment_method what method to use when applying treatment one of ("ratio" or "all infected"). ratio removes a portion of all infected and susceptibles, all infected removes all infected a portion of susceptibles.
#' @param natural_kernel_type what type of dispersal kernel should be used for natural dispersal can be ('cauchy', 'exponential')
#' @param anthropogenic_kernel_type what type of dispersal kernel should be used for anthropogenic dispersal can be ('cauchy', 'exponential')
#' @param natural_dir sets the predominate direction of natural dispersal usually due to wind values ('N', 'NW', 'W', 'SW', 'S', 'SE', 'E', 'NE', 'NONE')
#' @param anthropogenic_dir sets the predominate direction of anthropogenic dispersal usually due to human movement typically over long distances (e.g. nursery trade, movement of firewood, etc..) ('N', 'NW', 'W', 'SW', 'S', 'SE', 'E', 'NE', 'NONE')
#' @param pesticide_duration how long does the pestcide (herbicide, vaccine, etc..) last before the host is susceptible again. If value is 0 treatment is a culling (i.e. host removal) not a pesticide treatment. (needs to be the same length as treatment_dates and treatment_file)
#' @param pesticide_efficacy how effictive is the pesticide at preventing the disease or killing the pest (if this is 0.70 then when applied it successfully treats 70 percent of the plants or animals)
#' @param random_seed sets the random seed for the simulation used for reproducibility
#' @param output_frequency sets when outputs occur either ('year', 'month' or 'time step')
#' @param movements_file this is a csv file with columns lon_from, lat_from, lon_to, lat_to, number of animals, and date.
#' @param use_movements this is a boolean to turn on use of the movement module.
#' @param latency_period How many times steps does it take to for exposed populations become infected/infested. This is an integer value and must be greater than 0 if model type is SEI.
#' @param model_type What type of model most represents your sysetm. Options are "SEI" (Susceptible - Exposed - Infected/Infested) or "SI" (Susceptible - Infected/Infested). Default value is "SI".
#' @param parameter_means A vector of the means of the model parameters (reproductive_rate, natural_dispersal_distance, percent_natural_dispersal, anthropogenic_dispersal_distance, natural kappa, and anthropogenic kappa)
#' @param parameter_cov_matrix A covariance matrix from the previous years posterior parameter estimation ordered from (reproductive_rate, natural_dispersal_distance, percent_natural_dispersal, anthropogenic_dispersal_distance, natural kappa, and anthropogenic kappa)
#' @param start_exposed Do your initial conditions start as exposed or infected (only used if model_type is "SEI")
#' @useDynLib PoPS, .registration = TRUE
#' @importFrom raster raster values as.matrix xres yres stack extent calc extract rasterToPoints crs rowColFromCell
#' @importFrom Rcpp sourceCpp evalCpp
#' @importFrom  stats runif
#' @importFrom lubridate interval time_length mdy %within%
#' @importFrom utils read.csv 
#' @importFrom sp SpatialPointsDataFrame CRS spTransform
#' @return list of infected and susceptible per year
#' @export
#'
#' @examples
#' \dontrun{
#' infected_file <-  system.file("extdata", "SODexample", "initial_infection2001.tif", 
#' package = "PoPS")
#' host_file <- system.file("extdata", "SODexample", "host.tif", package = "PoPS")
#' total_plants_file <- system.file("extdata", "SODexample", "all_plants.tif", package = "PoPS")
#' temperature_coefficient_file <- system.file("extdata", "SODexample", "weather.tif", package = "
#' PoPS")
#' treatments_file <- system.file("extdata", "SODexample", "management.tif", package = "PoPS")
#' 
#' data <- pops(
#' infected_file, 
#' host_file, 
#' total_plants_file,
#' use_lethal_temperature = FALSE, 
#' temp = TRUE, precip = FALSE, 
#' management = TRUE, 
#' mortality_on = TRUE, 
#' temperature_file = "", 
#' temperature_coefficient_file, 
#' precipitation_coefficient_file ="", 
#' treatments_file,
#' season_month_start = 1, 
#' season_month_end = 12, 
#' time_step = "week",
#' start_date = '2001-01-01', 
#' end_date = 2005-12-31', 
#' treatment_dates = c('2001-12-24'),
#' natural_kernel_type = "cauchy",
#' lethal_temperature = -12.87, 
#' lethal_temperature_month = 1,
#' mortality_rate = 0.05, 
#' mortality_time_lag = 2,
#' treatment_date = 12, 
#' natural_dir = "NONE", 
#' kappa = 0, 
#' random_seed = NULL, 
#' output_frequency = "yearly")
#' }
#' 
pops <- function(infected_file, 
                 host_file, 
                 total_plants_file, 
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
                 start_date = '2008-01-01',
                 end_date = '2008-12-31', 
                 use_lethal_temperature = FALSE, 
                 temperature_file = "",
                 lethal_temperature = -12.87, 
                 lethal_temperature_month = 1,
                 mortality_on = FALSE, 
                 mortality_rate = 0, 
                 mortality_time_lag = 0, 
                 management = FALSE, 
                 treatment_dates = c('2008-12-24'), 
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
                 movements_file = "", 
                 use_movements = FALSE,
                 start_exposed = FALSE,
                 generate_stochasticity = FALSE,
                 establishment_stochasticity = FALSE,
                 movement_stochasticity = FALSE,
                 deterministic = FALSE,
                 establishment_probability = 0,
                 dispersal_percentage = 0.99){
  
  if (model_type == "SEI" && latency_period <= 0) {
    return("Model type is set to SEI but the latency period is less than 1")
  } else if (model_type == "SI" && latency_period > 0) {
    latency_period <- 0
  } 
  
  if (nrow(parameter_cov_matrix) != 6 | ncol(parameter_cov_matrix) != 6) {
    return("parameter covariance matrix is not 6 x 6")
  }
  
  if (length(parameter_means) != 6) {
    return("parameter means is not a vector of length 6")
  }
  
  parameters <- MASS::mvrnorm(1, parameter_means, parameter_cov_matrix)
  # names(parameters) <- c('reproductive_rate', 'natural_dispersal_distance', 'percent_natural_dispersal', 'anthropogenic_dispersal_distance', 'natural kappa', 'anthropogenic kappa')
  while(parameters[1] < 0 || parameters[2] < 0) {
    parameters <- mvrnorm(1, parameter_means, parameter_cov_matrix)
  }
  reproductive_rate <- parameters[1]
  natural_distance_scale <- parameters[2]
  percent_natural_dispersal <- parameters[3]
  if (percent_natural_dispersal > 1.000) {percent_natural_dispersal <- 1.000} 
  anthropogenic_distance_scale <- parameters[4]
  natural_kappa <- parameters[5]
  if (natural_kappa < 0.000) {natural_kappa <- 0}
  anthropogenic_kappa <- parameters[6]
  if (anthropogenic_kappa < 0.000) {anthropogenic_kappa <- 0}

  if (percent_natural_dispersal < 1.0) {
    use_anthropogenic_kernel <- TRUE
  } else {
    use_anthropogenic_kernel <- FALSE
  }
  
  treatment_metric_check <- treatment_metric_checks(treatment_method)
  if (!treatment_metric_check$checks_passed) {
    return(treatment_metric_check$failed_check)
  }
  
  time_check <- time_checks(end_date, start_date, time_step, output_frequency)
  if(time_check$checks_passed) {
    number_of_time_steps <- time_check$number_of_time_steps
    number_of_years <- time_check$number_of_years
    number_of_outputs <- time_check$number_of_outputs
  } else {
    return(time_check$failed_check)
  }
  
  if (is.null(random_seed)) {
    random_seed = round(stats::runif(1, 1, 1000000))
  }
  
  infected_check <- initial_raster_checks(infected_file)
  if (infected_check$checks_passed) {
    infected <- infected_check$raster
    if (raster::nlayers(infected) > 1) {
      infected <- output_from_raster_mean_and_sd(infected)
    }
  } else {
    return(infected_check$failed_check)
  }
  
  host_check <- secondary_raster_checks(host_file, infected)
  if (host_check$checks_passed) {
    host <- host_check$raster
    if (raster::nlayers(host) > 1) {
      host <- output_from_raster_mean_and_sd(host)
    }
  } else {
    return(host_check$failed_check)
  }
  
  total_plants_check <- secondary_raster_checks(total_plants_file, infected)
  if (total_plants_check$checks_passed) {
    total_plants <- total_plants_check$raster
    if (raster::nlayers(total_plants) > 1) {
      total_plants <- output_from_raster_mean_and_sd(total_plants)
    }
  } else {
    return(total_plants_check$failed_check)
  }
  
  susceptible <- host - infected
  susceptible[susceptible < 0] <- 0
  
  if (use_movements) {
    movements_check <- movement_checks(movements_file, infected, start_date, end_date)
    if (movements_check$checks_passed) {
      movements <- movements_check$movements
      movements_dates <- movements_check$movements_dates
      movements_r <- movements_check$movements_r
    } else {
      return(movements_check$failed_check)
    }
  } else {
    movements <- list(0,0,0,0,0)
    movements_dates <- start_date
  }
  
  if (use_lethal_temperature == TRUE) {
    temperature_check <- secondary_raster_checks(temperature_file, infected)
    if (temperature_check$checks_passed) {
      temperature_stack <- temperature_check$raster
    } else {
      return(temperature_check$failed_check)
    }
    
    temperature <- list(raster::as.matrix(temperature_stack[[1]]))
    for(i in 2:number_of_years) {
      temperature[[i]] <- raster::as.matrix(temperature_stack[[i]])
    }
  } else {
    temperature <- host
    raster::values(temperature) <- 1
    temperature <- list(raster::as.matrix(temperature))
  }
  
  weather <- FALSE
  if (temp == TRUE) {
    temperature_coefficient_check <- secondary_raster_checks(temperature_coefficient_file, infected)
    if (temperature_coefficient_check$checks_passed) {
      temperature_coefficient <- temperature_coefficient_check$raster
    } else {
      return(temperature_coefficient_check$failed_check)
    }
    
    weather <- TRUE
    weather_coefficient_stack <- temperature_coefficient
    if (precip ==TRUE){
      precipitation_coefficient_check <- secondary_raster_checks(precipitation_coefficient_file, infected)
      if (precipitation_coefficient_check$checks_passed) {
        precipitation_coefficient <- precipitation_coefficient_check$raster
      } else {
        return(precipitation_coefficient_check$failed_check)
      }
      
      weather_coefficient_stack <- weather_coefficient_stack * precipitation_coefficient
    }
  } else if(precip == TRUE){
    precipitation_coefficient_check <- secondary_raster_checks(precipitation_coefficient_file, infected)
    if (precipitation_coefficient_check$checks_passed) {
      precipitation_coefficient <- precipitation_coefficient_check$raster
    } else {
      return(precipitation_coefficient_check$failed_check)
    }
    
    weather <- TRUE
    weather_coefficient_stack <- precipitation_coefficient
  }
  
  if (weather == TRUE){
    # weather_coefficient_stack <- raster::reclassify(weather_coefficient_stack, matrix(c(NA,0), ncol = 2, byrow = TRUE), right = NA)
    weather_coefficient <- list(raster::as.matrix(weather_coefficient_stack[[1]]))
    for(i in 2:number_of_time_steps) {
      weather_coefficient[[i]] <- raster::as.matrix(weather_coefficient_stack[[i]])
    }
  } else {
    weather_coefficient <- host
    raster::values(weather_coefficient) <- 1
    weather_coefficient <- list(raster::as.matrix(weather_coefficient))
  }
  
  if (management == TRUE) {
    treatments_check <- secondary_raster_checks(treatments_file, infected)
    if (treatments_check$checks_passed) {
      treatment_stack <- treatments_check$raster
    } else {
      return(treatments_check$failed_check)
    }
    
    treatment_check <- treatment_checks(treatment_stack, treatments_file, pesticide_duration, treatment_dates, pesticide_efficacy)
    if (treatment_check$checks_passed) {
      treatment_maps <- treatment_check$treatment_maps
    } else {
      return(treatment_check$failed_check)
    }
  } else {
    treatment_map <- host
    raster::values(treatment_map) <- 0
    treatment_maps <- list(raster::as.matrix(treatment_map))
  }
  
  ew_res <- raster::xres(susceptible)
  ns_res <- raster::yres(susceptible)
  num_cols <- raster::ncol(susceptible)
  num_rows <- raster::nrow(susceptible)
  
  mortality_tracker <- infected
  raster::values(mortality_tracker) <- 0
  
  infected <- raster::as.matrix(infected)
  susceptible <- raster::as.matrix(susceptible)
  total_plants <- raster::as.matrix(total_plants)
  mortality_tracker <- raster::as.matrix(mortality_tracker)
  mortality <- mortality_tracker
  resistant <- mortality_tracker
  exposed <- list(mortality_tracker)
  
  if (latency_period > 1){
    for (ex in 2:(latency_period + 1)) {
      exposed[[ex]] <- mortality_tracker
    }
  }
  
  if (model_type == "SEI" & start_exposed) {
    exposed[[latency_period + 1]] <- infected
    infected <- mortality_tracker
  }
  
  data <- PoPS::pops_model(random_seed = random_seed, 
                           use_lethal_temperature = use_lethal_temperature, 
                           lethal_temperature = lethal_temperature, 
                           lethal_temperature_month = lethal_temperature_month,
                           infected = infected,
                           exposed = exposed,
                           susceptible = susceptible,
                           total_plants = total_plants,
                           mortality_on = mortality_on,
                           mortality_tracker = mortality_tracker,
                           mortality = mortality,
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
                           ew_res = ew_res, 
                           ns_res = ns_res, 
                           num_rows = num_rows, 
                           num_cols = num_cols,
                           time_step = time_step, 
                           reproductive_rate = reproductive_rate,
                           mortality_rate = mortality_rate, 
                           mortality_time_lag = mortality_time_lag,
                           season_month_start = season_month_start, 
                           season_month_end = season_month_end,
                           start_date = start_date, 
                           end_date = end_date,
                           treatment_method = treatment_method,
                           natural_kernel_type = natural_kernel_type, 
                           anthropogenic_kernel_type = anthropogenic_kernel_type, 
                           use_anthropogenic_kernel = use_anthropogenic_kernel, 
                           percent_natural_dispersal = percent_natural_dispersal,
                           natural_distance_scale = natural_distance_scale, 
                           anthropogenic_distance_scale = anthropogenic_distance_scale, 
                           natural_dir = natural_dir, 
                           natural_kappa = natural_kappa,
                           anthropogenic_dir = anthropogenic_dir, 
                           anthropogenic_kappa = anthropogenic_kappa,
                           output_frequency = output_frequency, 
                           model_type_ = model_type,
                           latency_period = latency_period,
                           generate_stochasticity = generate_stochasticity,
                           establishment_stochasticity = establishment_stochasticity,
                           movement_stochasticity = movement_stochasticity,
                           deterministic = deterministic,
                           establishment_probability = establishment_probability,
                           dispersal_percentage = dispersal_percentage
  )
  
  return(data)
  
}

