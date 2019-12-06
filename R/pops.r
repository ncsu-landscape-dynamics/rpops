#' PoPS (Pest or Pathogen Spread) model
#' 
#' A dynamic species distribution model for pest or pathogen spread in forest or agricultural ecosystems. The model is process based
#' meaning that it uses understanding of the effect of weather on reproduction and survival of the pest/pathogen in order to forecast
#' spread of the pest/pathogen into the future. 
#'
#' @param infected_file path to raster file with initial infections 
#' @param host_file path to raster file with number of hosts
#' @param total_plants_file path to raster file with number of total plants
#' @param temp boolean that allows the use of temperature coefficients to modify spread (TRUE or FALSE)
#' @param temperature_coefficient_file path to raster file with temperature coefficient data for the timestep and number of years specified
#' @param precip boolean that allows the use of precipitation coefficients to modify spread (TRUE or FALSE)
#' @param precipitation_coefficient_file path to raster file with precipitation coefficient data for the timestep and number of years specified
#' @param time_step how often should spread occur options: ('day', 'week', 'month')
#' @param reproductive_rate number of spores or pest units produced by a single host under optimal weather conditions 
#' @param season_month_start when does spread first start occurring in the year for your pest or pathogen (integer value between 1 and 12)
#' @param season_month_end when does spread end during the year for your pest or pathogen (integer value between 1 and 12)
#' @param start_time first year to start the simulation (needs to be a 4 digit year)
#' @param end_time last year of the simulation (needs to be a 4 digit year)
#' @param use_lethal_temperature a boolean to answer the question: does your pest or pathogen have a temperature at which it cannot survive? (TRUE or FALSE)
#' @param temperature_file path to raster file with temperature data for minimum temperature
#' @param lethal_temperature the temperature in degrees C at which lethal temperature related mortality occurs for your pest or pathogen (-50 and 60)
#' @param lethal_temperature_month the month in which lethal temperature related mortality occurs for your pest or pathogen integer value between 1 and 12
#' @param mortality_on  boolean to turn host mortality on and off (TRUE or FALSE)
#' @param mortality_rate rate at which mortality occurs value between 0 and 1 
#' @param mortality_time_lag time lag from infection until mortality can occur in years integer >= 1
#' @param management boolean to allow use of managemnet (TRUE or FALSE)
#' @param treatment_dates dates in which to apply treatment list with format (YYYY_MM_DD)
#' @param treatments_file path to raster file with treatment data by dates
#' @param treatment_method what method to use when applying treatment one of ("ratio" or "all infected"). ratio removes a portion of all infected and susceptibles, all infected removes all infected a portion of susceptibles.
#' @param percent_natural_dispersal  what percentage of dispersal is natural range versus anthropogenic range value between 0 and 1
#' @param natural_kernel_type what type of dispersal kernel should be used for natural dispersal can be ('cauchy', 'exponential')
#' @param anthropogenic_kernel_type what type of dispersal kernel should be used for anthropogenic dispersal can be ('cauchy', 'exponential')
#' @param natural_distance_scale distance scale parameter for natural range dispersal kernel numeric value > 0 
#' @param anthropogenic_distance_scale distance scale parameter for anthropogenic range dispersal kernel numeric value > 0 
#' @param natural_dir sets the predominate direction of natural dispersal usually due to wind values ('N', 'NW', 'W', 'SW', 'S', 'SE', 'E', 'NE', 'NONE')
#' @param natural_kappa sets the strength of the natural direction in the von-mises distribution numeric value between 0.01 and 12
#' @param anthropogenic_dir sets the predominate direction of anthropogenic dispersal usually due to human movement typically over long distances (e.g. nursery trade, movement of firewood, etc..) ('N', 'NW', 'W', 'SW', 'S', 'SE', 'E', 'NE', 'NONE')
#' @param anthropogenic_kappa sets the strength of the anthropogenic direction in the von-mises distribution numeric value between 0.01 and 12
#' @param pesticide_duration how long does the pestcide (herbicide, vaccine, etc..) last before the host is susceptible again. If value is 0 treatment is a culling (i.e. host removal) not a pesticide treatment.
#' @param pesticide_efficacy how effictive is the pesticide at preventing the disease or killing the pest (if this is 0.70 then when applied it successfully treats 70 percent of the plants or animals)
#' @param random_seed sets the random seed for the simulation used for reproducibility
#' 
#' @useDynLib PoPS, .registration = TRUE
#' @importFrom raster raster values as.matrix xres yres stack extent
#' @importFrom Rcpp sourceCpp evalCpp
#' @importFrom  stats runif
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
#' data <- pops(infected_file, host_file, total_plants_file, reproductive_rate = 1.0,
#' use_lethal_temperature = FALSE, temp = TRUE, precip = FALSE, management = TRUE, 
#' mortality_on = TRUE, temperature_file = "", temperature_coefficient_file, 
#' precipitation_coefficient_file ="", treatments_file,
#' season_month_start = 1, season_month_end = 12, time_step = "week",
#' start_time = 2001, end_time = 2005, treatment_dates = c('2001-12-24'),
#' natural_kernel_type = "cauchy", percent_natural_dispersal = 1.0,
#' natural_distance_scale = 20.57, anthropogenic_distance_scale = 0.0,
#' lethal_temperature = -12.87, lethal_temperature_month = 1,
#' mortality_rate = 0.05, mortality_time_lag = 2,
#' treatment_date = 12, natural_dir = "NONE", kappa = 0, random_seed = NULL)
#' }
#' 
pops <- function(infected_file, host_file, total_plants_file, 
                 temp = FALSE, temperature_coefficient_file = "", 
                 precip = FALSE, precipitation_coefficient_file = "", 
                 time_step = "month", reproductive_rate = 3.0,
                 season_month_start = 1, season_month_end = 12, 
                 start_time = 2018, end_time = 2020, 
                 use_lethal_temperature = FALSE, temperature_file = "",
                 lethal_temperature = -12.87, lethal_temperature_month = 1,
                 mortality_on = FALSE, mortality_rate = 0, mortality_time_lag = 0, 
                 management = FALSE, treatment_dates = c('2000-12-24'), treatments_file = "",
                 treatment_method = "ratio",
                 percent_natural_dispersal = 1.0,
                 natural_kernel_type = "cauchy", anthropogenic_kernel_type = "cauchy",
                 natural_distance_scale = 21, anthropogenic_distance_scale = 0.0,
                 natural_dir = "NONE", natural_kappa = 0, 
                 anthropogenic_dir = "NONE", anthropogenic_kappa = 0,
                 pesticide_duration = c(0), pesticide_efficacy = 1.0,
                 random_seed = NULL){ 
  
  if (!treatment_method %in% c("ratio", "all infected")) {
    return("treatment method is not one of the valid treatment options")
  }
  
  if (!file.exists(infected_file)) {
    return("Infected file does not exist") 
  }
  
  if (!(raster::extension(infected_file) %in% c(".grd", ".tif", ".img"))) {
    return("Infected file is not one of '.grd', '.tif', '.img'")
  }
  
  if (!file.exists(host_file)) {
    return("Host file does not exist") 
  }
  
  if (!(raster::extension(host_file) %in% c(".grd", ".tif", ".img"))) {
    return("Host file is not one of '.grd', '.tif', '.img'")
  }
  
  if (!file.exists(total_plants_file)) {
    return("Total plants file does not exist") 
  }
  
  if (!(raster::extension(total_plants_file) %in% c(".grd", ".tif", ".img"))) {
    return("Total plants file is not one of '.grd', '.tif', '.img'")
  }
  
  if (!(time_step %in% list("week", "month", "day"))) {
    return("Time step must be one of 'week', 'month' or 'day'")
  }
  
  if (class(end_time) != "numeric" || nchar(end_time) != 4 || class(start_time) != "numeric" || nchar(start_time) != 4){
    return("End time and/or start time not of type numeric and/or in format YYYY")
  }
  
  if (is.null(random_seed)) {
    random_seed = round(stats::runif(1, 1, 1000000))
  }
  
  if (time_step == "week") {
    number_of_time_steps <- (end_time-start_time+1)*52
  } else if (time_step == "month") {
    number_of_time_steps <- (end_time-start_time+1)*12
  } else if (time_step == "day") {
    number_of_time_steps <- (end_time-start_time+1)*365
  }
  
  number_of_years <- end_time-start_time+1
  
  infected <- raster::raster(infected_file)
  infected[is.na(infected)] <- 0
  host <- raster::raster(host_file)
  host[is.na(host)] <- 0
  total_plants <- raster::raster(total_plants_file)
  total_plants[is.na(total_plants)] <- 0
  
  if (!(raster::extent(infected) == raster::extent(host) && raster::extent(infected) == raster::extent(total_plants))) {
    return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
  }
  
  if (!(raster::xres(infected) == raster::xres(host) && raster::xres(infected) == raster::xres(total_plants) && raster::yres(infected) == raster::yres(host) && raster::yres(infected) == raster::yres(total_plants))) {
    return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  }
  
  if (!(raster::compareCRS(host,infected) && raster::compareCRS(host, total_plants))) {
    return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
  }
  
  susceptible <- host - infected
  susceptible[is.na(susceptible)] <- 0
  susceptible[susceptible < 0] <- 0
  
  if (use_lethal_temperature == TRUE  && !file.exists(temperature_file)) {
    return("Temperature file does not exist")
  }
  
  if (use_lethal_temperature == TRUE  && !(raster::extension(temperature_file) %in% c(".grd", ".tif", ".img"))) {
    return("Temperature file is not one of '.grd', '.tif', '.img'")
  }
  
  if (use_lethal_temperature == TRUE) {
    temperature_stack <- raster::stack(temperature_file)
    temperature_stack[is.na(temperature_stack)] <- 0
    
    if (!(raster::extent(infected) == raster::extent(temperature_stack))) {
      return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
    }
    
    if (!(raster::xres(infected) == raster::xres(temperature_stack) && raster::yres(infected) == raster::yres(temperature_stack))) {
      return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
    }
    
    if (!(raster::compareCRS(infected, temperature_stack))) {
      return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
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
  
  if (temp == TRUE  && !file.exists(temperature_coefficient_file)) {
    return("Temperature coefficient file does not exist")
  }
  
  if (temp == TRUE  && !(raster::extension(temperature_coefficient_file) %in% c(".grd", ".tif", ".img"))) {
    return("Temperature coefficient file is not one of '.grd', '.tif', '.img'")
  }
  
  if (precip == TRUE  && !file.exists(precipitation_coefficient_file)) {
    return("Precipitation coefficient file does not exist")
  }
  
  if (precip == TRUE  && !(raster::extension(precipitation_coefficient_file) %in% c(".grd", ".tif", ".img"))) {
    return("Precipitation coefficient file is not one of '.grd', '.tif', '.img'")
  }
  
  weather <- FALSE
  if (temp == TRUE) {
    temperature_coefficient <- raster::stack(temperature_coefficient_file)
    
    if (!(raster::extent(infected) == raster::extent(temperature_coefficient))) {
      return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
    }
    
    if (!(raster::xres(infected) == raster::xres(temperature_coefficient) && raster::yres(infected) == raster::yres(temperature_coefficient))) {
      return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
    }
    
    if (!(raster::compareCRS(infected, temperature_coefficient))) {
      return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
    }
    
    weather <- TRUE
    weather_coefficient_stack <- temperature_coefficient
    if (precip ==TRUE){
      precipitation_coefficient <- raster::stack(precipitation_coefficient_file)
      
      if (!(raster::extent(infected) == raster::extent(precipitation_coefficient))) {
        return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
      }
      
      if (!(raster::xres(infected) == raster::xres(precipitation_coefficient) && raster::yres(infected) == raster::yres(precipitation_coefficient))) {
        return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
      }
      
      if (!(raster::compareCRS(infected, precipitation_coefficient))) {
        return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
      }
      
      weather_coefficient_stack <- weather_coefficient_stack * precipitation_coefficient
    }
  } else if(precip == TRUE){
    precipitation_coefficient <- raster::stack(precipitation_coefficient_file)
    
    if (!(raster::extent(infected) == raster::extent(precipitation_coefficient))) {
      return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
    }
    
    if (!(raster::xres(infected) == raster::xres(precipitation_coefficient) && raster::yres(infected) == raster::yres(precipitation_coefficient))) {
      return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
    }
    
    if (!(raster::compareCRS(infected, precipitation_coefficient))) {
      return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
    }
    
    weather <- TRUE
    weather_coefficient_stack <- precipitation_coefficient
  }
  
  if (weather == TRUE){
    weather_coefficient_stack[is.na(weather_coefficient_stack)] <- 0
    weather_coefficient <- list(raster::as.matrix(weather_coefficient_stack[[1]]))
    for(i in 2:number_of_time_steps) {
      weather_coefficient[[i]] <- raster::as.matrix(weather_coefficient_stack[[i]])
    }
  } else {
    weather_coefficient <- host
    raster::values(weather_coefficient) <- 1
    weather_coefficient <- list(raster::as.matrix(weather_coefficient))
  }
  
  if (management == TRUE  && !file.exists(treatments_file)) {
    return("Treatments file does not exist")
  }
  
  if (management == TRUE  && !(raster::extension(treatments_file) %in% c(".grd", ".tif", ".img"))) {
    return("Treatments file is not one of '.grd', '.tif', '.img'")
  }
  
  if (management == TRUE) {
    treatment_stack <- raster::stack(treatments_file)
    treatment_stack[is.na(treatment_stack)] <- 0
    
    if (!(raster::extent(infected) == raster::extent(treatment_stack))) {
      return("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
    }
    
    if (!(raster::xres(infected) == raster::xres(treatment_stack) && raster::yres(infected) == raster::yres(treatment_stack))) {
      return("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
    }
    
    if (!(raster::compareCRS(infected, treatment_stack))) {
      return("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
    }
    
    if (length(treatments_file) != length(treatment_dates)) {
      return("Length of list for treatment dates and treatments_file must be equal")
    }
    if (length(pesticide_duration) != length(treatment_dates)) {
      return("Length of list for treatment dates and pesticide_duration must be equal")
    }
    if (pesticide_duration[1] > 0) {
      treatment_stack[[1]] <- treatment_stack[[1]] * pesticide_efficacy
    }
    treatment_maps <- list(raster::as.matrix(treatment_stack[[1]]))
    if (raster::nlayers(treatment_stack) >= 2) {
      for(i in 2:raster::nlayers(treatment_stack)) {
        if (pesticide_duration[i] > 0) {
          treatment_stack[[i]] <- treatment_stack[[i]] * pesticide_efficacy
        }
        treatment_maps[[i]] <- raster::as.matrix(treatment_stack[[i]])
      }
    }
    treatment_dates <- treatment_dates
    pesticide_duration <- pesticide_duration
    treatment_method <- treatment_method
  } else {
    treatment_map <- host
    raster::values(treatment_map) <- 0
    treatment_maps <- list(raster::as.matrix(treatment_map))
    treatment_method <- treatment_method
    pesticide_duration <- pesticide_duration
  }
  
  if(percent_natural_dispersal == 1.0) {
    use_anthropogenic_kernel = FALSE
  } else if (percent_natural_dispersal < 1.0  && percent_natural_dispersal >= 0.0) {
    use_anthropogenic_kernel = TRUE
  } else {
    return("Percent natural dispersal must be between 0.0 and 1.0")
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
  
  data <- PoPS::pops_model(random_seed = random_seed, 
                     use_lethal_temperature = use_lethal_temperature, 
                     lethal_temperature = lethal_temperature, lethal_temperature_month = lethal_temperature_month,
                     infected = infected,
                     susceptible = susceptible,
                     total_plants = total_plants,
                     mortality_on = mortality_on,
                     mortality_tracker = mortality_tracker,
                     mortality = mortality,
                     treatment_maps = treatment_maps,
                     treatment_dates = treatment_dates,
                     pesticide_duration = pesticide_duration,
                     resistant = resistant,
                     weather = weather,
                     temperature = temperature,
                     weather_coefficient = weather_coefficient,
                     ew_res = ew_res, ns_res = ns_res, num_rows = num_rows, num_cols = num_cols,
                     time_step = time_step, reproductive_rate = reproductive_rate,
                     mortality_rate = mortality_rate, mortality_time_lag = mortality_time_lag,
                     season_month_start = season_month_start, season_month_end = season_month_end,
                     start_time = start_time, end_time = end_time,
                     treatment_method = treatment_method,
                     natural_kernel_type = natural_kernel_type, anthropogenic_kernel_type = anthropogenic_kernel_type, 
                     use_anthropogenic_kernel = use_anthropogenic_kernel, percent_natural_dispersal = percent_natural_dispersal,
                     natural_distance_scale = natural_distance_scale, anthropogenic_distance_scale = anthropogenic_distance_scale, 
                     natural_dir = natural_dir, natural_kappa = natural_kappa,
                     anthropogenic_dir = anthropogenic_dir, anthropogenic_kappa = anthropogenic_kappa
  )
  
  return(data)
  
}

