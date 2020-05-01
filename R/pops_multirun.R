#' PoPS (Pest or Pathogen Spread) model Multiple Runs
#' 
#' A dynamic species distribution model for pest or pathogen spread in forest or agricultural ecosystems. The model is process based
#' meaning that it uses understanding of the effect of weather on reproduction and survival of the pest/pathogen in order to forecast
#' spread of the pest/pathogen into the future. 
#'
#' @inheritParams pops
#' @param num_iterations how many iterations do you want to run to allow the calibration to converge at least 10 
#' @param number_of_cores enter how many cores you want to use (default = NA). If not set uses the # of CPU cores - 1. must be an integer >= 1
#'
#' @importFrom raster raster values as.matrix xres yres stack reclassify cellStats nlayers calc extract rasterToPoints
#' @importFrom stats runif rnorm median sd
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar%
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom iterators icount
#' @importFrom lubridate interval time_length year
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
#' start_date = '2001-01-01', end_date = 2005-12-31', treatment_dates = c(2001,2002,2003,2004,2005),
#' natural_kernel_type = "cauchy", percent_natural_dispersal = 1.0,
#' natural_distance_scale = 20.57, anthropogenic_distance_scale = 0.0,
#' lethal_temperature = -12.87, lethal_temperature_month = 1,
#' mortality_rate = 0.05, mortality_time_lag = 2,
#' treatment_date = 12, natural_dir = "NONE", kappa = 0, random_seed = NULL)
#' }
#' 
pops_multirun <- function(infected_file, 
                          host_file, 
                          total_plants_file, 
                          temp = FALSE, 
                          temperature_coefficient_file = "", 
                          precip = FALSE, 
                          precipitation_coefficient_file = "",
                          model_type = "SI", 
                          latency_period = 0,
                          time_step = "month", 
                          reproductive_rate = 3.0,
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
                          treatment_dates = c(0), 
                          treatments_file = "",
                          treatment_method = "ratio",
                          percent_natural_dispersal = 1.0,
                          natural_kernel_type = "cauchy", 
                          anthropogenic_kernel_type = "cauchy",
                          natural_distance_scale = 21, 
                          anthropogenic_distance_scale = 0.0,
                          natural_dir = "NONE", 
                          natural_kappa = 0, 
                          anthropogenic_dir = "NONE", 
                          anthropogenic_kappa = 0,
                          num_iterations = 100, 
                          number_of_cores = NA,
                          pesticide_duration = 0,
                          pesticide_efficacy = 1.0,
                          random_seed = NULL,
                          output_frequency = "year",
                          movements_file = "", 
                          use_movements = FALSE){ 
  
  if (model_type == "SEI" && latency_period <= 0) {
    return("Model type is set to SEI but the latency period is less than 1")
  } else if (model_type == "SI" && latency_period > 0) {
    latency_period <- 0
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
  
  percent_check <- percent_checks(percent_natural_dispersal)
  if (percent_check$checks_passed){
    use_anthropogenic_kernel <- percent_check$use_anthropogenic_kernel
  } else {
    return(percent_check$failed_check)
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
  
  reproductive_rate_check <- uncertainty_check(reproductive_rate, round_to = 1, n = num_iterations)
  if (reproductive_rate_check$checks_passed) {
    reproductive_rate <- reproductive_rate_check$value
  } else {
    return(reproductive_rate_check$failed_check)
  }
  
  natural_distance_scale_check <- uncertainty_check(natural_distance_scale, round_to = 0, n = num_iterations)
  if (natural_distance_scale_check$checks_passed) {
    natural_distance_scale <- natural_distance_scale_check$value
  } else {
    return(natural_distance_scale_check$failed_check)
  }
  
  anthropogenic_distance_scale_check <- uncertainty_check(anthropogenic_distance_scale, round_to = 0, n = num_iterations)
  if (anthropogenic_distance_scale_check$checks_passed) {
    anthropogenic_distance_scale <- anthropogenic_distance_scale_check$value
  } else {
    return(anthropogenic_distance_scale_check$failed_check)
  }
  
  percent_natural_dispersal_check <- uncertainty_check(percent_natural_dispersal, round_to = 3, n = num_iterations)
  if (percent_natural_dispersal_check$checks_passed) {
    percent_natural_dispersal <- percent_natural_dispersal_check$value
  } else {
    return(percent_natural_dispersal_check$failed_check)
  }
  
  years <- seq(year(start_date), year(end_date), 1)
  rcl <- c(1, Inf, 1, 0, 0.99, NA)
  rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
  
  if (is.na(number_of_cores) || number_of_cores > parallel::detectCores()) {
    core_count <- parallel::detectCores() - 1
  } else {
    core_count <- number_of_cores
  }
  cl <- makeCluster(core_count)
  registerDoParallel(cl)
  
  infected_stack <- foreach::foreach(i = 1:num_iterations, .combine = c, .packages = c("raster", "PoPS"), .export = ls(globalenv())) %dopar% {
    random_seed <- round(stats::runif(1, 1, 1000000))
    data <- pops_model(random_seed = random_seed, 
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
                       reproductive_rate = reproductive_rate[i],
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
                       percent_natural_dispersal = percent_natural_dispersal[i],
                       natural_distance_scale = natural_distance_scale[i], 
                       anthropogenic_distance_scale = anthropogenic_distance_scale[i], 
                       natural_dir = natural_dir, 
                       natural_kappa = natural_kappa,
                       anthropogenic_dir = anthropogenic_dir, 
                       anthropogenic_kappa = anthropogenic_kappa,
                       output_frequency = output_frequency,
                       model_type_ = model_type,
                       latency_period = latency_period
    )
    
    comp_years <- raster::stack(lapply(1:length(data$infected), function(i) host))
    susceptible_runs <- raster::stack(lapply(1:length(data$infected), function(i) host))
    
    for (q in 1:raster::nlayers(comp_years)) {
      comp_years[[q]] <- data$infected[[q]]
      susceptible_runs[[q]] <- data$susceptible[[q]]
    }
    
    number_infected <- data$number_infected
    spread_rate <- data$rates
    infected_area <- data$area_infected
    single_run <- comp_years
    comp_years <- raster::reclassify(comp_years, rclmat)
    comp_years <- raster::reclassify(comp_years, matrix(c(NA,0), ncol = 2, byrow = TRUE), right = NA)
    infected_stack <- comp_years
    data <- list(single_run, infected_stack, number_infected, susceptible_runs, infected_area, spread_rate)
  }
  
  stopCluster(cl)
  single_runs <- infected_stack[seq(1,length(infected_stack),6)]
  probability_runs <- infected_stack[seq(2,length(infected_stack),6)]
  number_infected_runs <- infected_stack[seq(3,length(infected_stack),6)]
  susceptible_runs <- infected_stack[seq(4,length(infected_stack),6)]
  area_infected_runs <- infected_stack[seq(5,length(infected_stack),6)]
  spread_rate_runs <- infected_stack[seq(6,length(infected_stack),6)]
  
  prediction <- probability_runs[[1]]
  prediction[prediction > 0] <- 0
  infected_area <- data.frame(t(years))
  infected_number <- data.frame(t(years))
  west_rates <- data.frame(t(years))
  east_rates <- data.frame(t(years))
  south_rates <- data.frame(t(years))
  north_rates <- data.frame(t(years))
  max_values <- data.frame(t(years))
  
  for (i in 1:length(probability_runs)) {
    prediction <- prediction + probability_runs[[i]]
    infected_number[i,] <- number_infected_runs[[i]]
    infected_area[i,] <- area_infected_runs[[i]]
    rates <- do.call(rbind, spread_rate_runs[[i]])
    west_rates[i,] <- rates[,4]
    east_rates[i,] <- rates[,3]
    south_rates[i,] <- rates[,2]
    north_rates[i,] <- rates[,1]
    max_values[i,] <- raster::maxValue(single_runs[[i]])
  }
  
  
  probability <- (prediction/(length(probability_runs))) * 100
  
  infected_areas <- round(sapply(infected_area, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  number_infecteds <- round(sapply(infected_number, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  west_rate <- round(sapply(west_rates, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  east_rate <- round(sapply(east_rates, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  south_rate <- round(sapply(south_rates, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  north_rate <- round(sapply(north_rates, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  
  which_median <- function(x) raster::which.min(abs(x - median(x)))
  
  median_run_index <- which_median(infected_number[[1]])
  
  single_run <- single_runs[[median_run_index]]
  susceptible_run <- susceptible_runs[[median_run_index]]
  
  single_run_out <- single_run
  susceptible_run_out <- susceptible_run
  
  outputs <- list(probability, single_run_out, number_infecteds, infected_areas, west_rate, east_rate, south_rate, north_rate)
  names(outputs) <- c('probability', 'single_run_out', 'number_infecteds', 'infected_areas', 'west_rate', 'east_rate', 'south_rate', 'north_rate')
  
  return(outputs)
  
}

