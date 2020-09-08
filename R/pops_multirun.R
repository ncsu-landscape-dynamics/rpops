#' PoPS (Pest or Pathogen Spread) model Multiple Runs
#' 
#' A dynamic species distribution model for pest or pathogen spread in forest or agricultural ecosystems. The model is process based
#' meaning that it uses understanding of the effect of weather on reproduction and survival of the pest/pathogen in order to forecast
#' spread of the pest/pathogen into the future. 
#'
#' @inheritParams pops
#' @param number_of_iterations how many iterations do you want to run to allow the calibration to converge at least 10 
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
#' total_populations_file <- system.file("extdata", "SODexample", "all_plants.tif", package = "PoPS")
#' temperature_coefficient_file <- system.file("extdata", "SODexample", "weather.tif", package = "
#' PoPS")
#' treatments_file <- system.file("extdata", "SODexample", "management.tif", package = "PoPS")
#' 
#' data <- pops(infected_file, 
#' host_file, 
#' total_populations_file, 
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
#' treatment_dates = c(2001,2002,2003,2004,2005),
#' natural_kernel_type = "cauchy", 
#' lethal_temperature = -12.87, 
#' lethal_temperature_month = 1,
#' mortality_rate = 0.05, 
#' mortality_time_lag = 2,
#' treatment_date = 12, 
#' natural_dir = "NONE", 
#' kappa = 0, 
#' random_seed = NULL)
#' }
#' 
pops_multirun <- function(infected_file, 
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
                          natural_kernel_type = "cauchy", 
                          anthropogenic_kernel_type = "cauchy",
                          natural_dir = "NONE", 
                          anthropogenic_dir = "NONE", 
                          number_of_iterations = 100, 
                          number_of_cores = NA,
                          pesticide_duration = 0,
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
                          deterministic = FALSE,
                          establishment_probability = 0.5,
                          dispersal_percentage = 0.99,
                          quarantine_areas_file = "",
                          use_quarantine = FALSE,
                          use_spreadrates = FALSE){ 
  
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
  
  parameters <- data.frame(MASS::mvrnorm(number_of_iterations, parameter_means, parameter_cov_matrix))
  names(parameters) <- c('reproductive_rate', 'natural_dispersal_distance', 'percent_natural_dispersal', 'anthropogenic_dispersal_distance', 'natural kappa', 'anthropogenic kappa')
  while(any(parameters[,1] < 0) || any(parameters[,2] < 0)) {
    parameters[parameters[,1] < 0 | parameters[,2] <= 0,] <- mvrnorm(nrow(parameters[parameters[,1] < 0 | parameters[,2] < 0,]), parameter_means, parameter_cov_matrix)
  }
  reproductive_rate <- parameters[[1]]
  natural_distance_scale <- parameters[[2]]
  percent_natural_dispersal <- parameters[[3]]
  if (any(percent_natural_dispersal > 1.000)) {percent_natural_dispersal[percent_natural_dispersal > 1.00] <- 1.000} 
  anthropogenic_distance_scale <- parameters[[4]]
  natural_kappa <- parameters[[5]]
  if (any(natural_kappa < 0.000)) {natural_kappa[natural_kappa < 0.000] <- 0}
  anthropogenic_kappa <- parameters[[6]]
  if (any(anthropogenic_kappa < 0.000)) {anthropogenic_kappa[anthropogenic_kappa < 0.000] <- 0}
  
  if (any(percent_natural_dispersal < 1.0)) {
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
    quarantine_frequency <- output_frequency
    quarantine_frequency_n <- output_frequency_n
    spreadrate_frequency <- output_frequency
    spreadrate_frequency_n <- output_frequency_n
    output_frequency <- time_check$output_frequency
  } else {
    return(time_check$failed_check)
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
  
  total_populations_check <- secondary_raster_checks(total_populations_file, infected)
  if (total_populations_check$checks_passed) {
    total_populations <- total_populations_check$raster
    if (raster::nlayers(total_populations) > 1) {
      total_populations <- output_from_raster_mean_and_sd(total_populations)
    }
  } else {
    return(total_populations_check$failed_check)
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
  total_populations <- raster::as.matrix(total_populations)
  mortality_tracker <- raster::as.matrix(mortality_tracker)
  mortality <- mortality_tracker
  resistant <- mortality_tracker
  exposed <- list(mortality_tracker)
  
  if (use_quarantine){
    quarantine_check <- secondary_raster_checks(quarantine_areas_file, host)
    if (quarantine_check$checks_passed) {
      quarantine_areas <- quarantine_check$raster
      quarantine_areas <- raster::as.matrix(quarantine_areas)
    } else {
      return(quarantine_check$failed_check)
    }
  } else {
    # set quarantine areas to all zeros (meaning no quarantine areas are considered)
    quarantine_areas <- mortality_tracker
  }
  
  if (latency_period > 1){
    for (ex in 2:(latency_period + 1)) {
      exposed[[ex]] <- mortality_tracker
    }
  }
  
  if (model_type == "SEI" & start_exposed) {
    exposed[[latency_period + 1]] <- infected
    infected <- mortality_tracker
  }
  
  # years <- seq(year(start_date), year(end_date), 1)
  rcl <- c(1, Inf, 1, 0, 0.99, NA)
  rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
  
  if (is.na(number_of_cores) || number_of_cores > parallel::detectCores()) {
    core_count <- parallel::detectCores() - 1
  } else {
    core_count <- number_of_cores
  }
  cl <- makeCluster(core_count)
  registerDoParallel(cl)
  
  infected_stack <- foreach::foreach(i = 1:number_of_iterations, .combine = c, .packages = c("raster", "PoPS")) %dopar% {
    random_seed <- round(stats::runif(1, 1, 1000000))
    data <- PoPS::pops_model(random_seed = random_seed, 
                             use_lethal_temperature = use_lethal_temperature, 
                             lethal_temperature = lethal_temperature, 
                             lethal_temperature_month = lethal_temperature_month,
                             infected = infected,
                             exposed = exposed,
                             susceptible = susceptible,
                             total_populations = total_populations,
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
                             natural_kappa = natural_kappa[i],
                             anthropogenic_dir = anthropogenic_dir, 
                             anthropogenic_kappa = anthropogenic_kappa[i],
                             output_frequency = output_frequency,
                             output_frequency_n = output_frequency_n,
                             quarantine_frequency = quarantine_frequency,
                             quarantine_frequency_n = quarantine_frequency_n,
                             use_quarantine = use_quarantine,
                             spreadrate_frequency = spreadrate_frequency,
                             spreadrate_frequency_n = spreadrate_frequency_n,
                             use_spreadrates = use_spreadrates,
                             model_type_ = model_type,
                             latency_period = latency_period,
                             generate_stochasticity = generate_stochasticity,
                             establishment_stochasticity = establishment_stochasticity,
                             movement_stochasticity = movement_stochasticity,
                             deterministic = deterministic,
                             establishment_probability = establishment_probability,
                             dispersal_percentage = dispersal_percentage
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
    ## add quarantine here
    quarantine_escape <- data$quarantine_escape
    quarantine_escape_distance <- data$quarantine_escape_distance
    # quarantine_escape_direction <- data$quarantine_escape_direction
    quarantine_escape_direction <- 0
    
    infected_stack <- comp_years
    data <- list(single_run, infected_stack, number_infected, susceptible_runs, infected_area, spread_rate, quarantine_escape, quarantine_escape_distance, quarantine_escape_direction)
  }
  
  stopCluster(cl)
  single_runs <- infected_stack[seq(1,length(infected_stack),9)]
  probability_runs <- infected_stack[seq(2,length(infected_stack),9)]
  number_infected_runs <- infected_stack[seq(3,length(infected_stack),9)]
  susceptible_runs <- infected_stack[seq(4,length(infected_stack),9)]
  area_infected_runs <- infected_stack[seq(5,length(infected_stack),9)]
  spread_rate_runs <- infected_stack[seq(6,length(infected_stack),9)]
  ## add quarantine here
  quarantine_escape_runs <- infected_stack[seq(7,length(infected_stack),9)]
  quarantine_escape_distance_runs <- infected_stack[seq(8,length(infected_stack),9)]
  quarantine_escape_directions_runs <- infected_stack[seq(9,length(infected_stack),9)]
  
  prediction <- probability_runs[[1]]
  prediction[prediction > 0] <- 0
  escape_probability <- data.frame(t(rep(0, nlayers(probability_runs[[1]]))))
  infected_area <- data.frame(t(rep(0, nlayers(probability_runs[[1]]))))
  infected_number <- data.frame(t(rep(0, nlayers(probability_runs[[1]]))))
  west_rates <- data.frame(t(rep(0, nlayers(probability_runs[[1]]))))
  east_rates <- data.frame(t(rep(0, nlayers(probability_runs[[1]]))))
  south_rates <- data.frame(t(rep(0, nlayers(probability_runs[[1]]))))
  north_rates <- data.frame(t(rep(0, nlayers(probability_runs[[1]]))))
  max_values <- data.frame(t(rep(0, nlayers(probability_runs[[1]]))))
  ## add quarantine here
  quarantine_escapes <- data.frame(t(rep(0, nlayers(probability_runs[[1]]))))
  quarantine_escape_distances <- data.frame(t(rep(0, nlayers(probability_runs[[1]]))))
  quarantine_escape_directions <- data.frame(t(rep(0, nlayers(probability_runs[[1]]))))
  
  for (i in 1:length(probability_runs)) {
    prediction <- prediction + probability_runs[[i]]
    infected_number[i,] <- number_infected_runs[[i]]
    infected_area[i,] <- area_infected_runs[[i]]
    rates <- do.call(rbind, spread_rate_runs[[i]])
    if (!is.null(rates)) {
      west_rates[i,] <- rates[,4]
      east_rates[i,] <- rates[,3]
      south_rates[i,] <- rates[,2]
      north_rates[i,] <- rates[,1]
    }
    ## add quarantine here
    if (use_quarantine & length(quarantine_escape_runs[[i]]) == nlayers(probability_runs[[i]])) {
      escape_probability <- escape_probability + quarantine_escape_runs[[i]]
      # quarantine_escapes[i,] <- quarantine_escape_runs[[i]]
      quarantine_escape_distances <- quarantine_escape_distance_runs[[i]]
      quarantine_escape_directions <- quarantine_escape_directions_runs[[i]]
    }
    max_values[i,] <- raster::maxValue(single_runs[[i]])
  }
  
  probability <- (prediction/(length(probability_runs))) * 100
  
  infected_areas <- round(sapply(infected_area, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  number_infecteds <- round(sapply(infected_number, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  west_rate <- round(sapply(west_rates, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  east_rate <- round(sapply(east_rates, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  south_rate <- round(sapply(south_rates, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  north_rate <- round(sapply(north_rates, function(x) c( "Mean"= mean(x,na.rm=TRUE),"Stand dev" = sd(x))), digits = 0)
  ## add quarantine here
  if (use_quarantine) {
    escape_probability <- escape_probability/length(probability_runs) * 100
    north_distance_to_quarantine <- 0
  } else {
    escape_probability <- data.frame(t(rep(NA, nlayers(probability_runs[[1]]))))
    north_distance_to_quarantine <- data.frame(t(rep(NA, nlayers(probability_runs[[1]]))))
    south_distance_to_quarantine <- data.frame(t(rep(NA, nlayers(probability_runs[[1]]))))
    east_distance_to_quarantine <- data.frame(t(rep(NA, nlayers(probability_runs[[1]]))))
    west_distance_to_quarantine <- data.frame(t(rep(NA, nlayers(probability_runs[[1]]))))
  }
  
  
  which_median <- function(x) raster::which.min(abs(x - median(x)))
  
  median_run_index <- which_median(infected_number[[1]])
  
  single_run <- single_runs[[median_run_index]]
  susceptible_run <- susceptible_runs[[median_run_index]]
  
  single_run_out <- single_run
  susceptible_run_out <- susceptible_run
  
  raster_stacks_list <- list()
  simulation_mean_stack <- stack()
  simulation_sd_stack <- stack()
  for (q in 1:nlayers(single_runs[[1]])){
    raster_stacks <- stack()
    for (j in 1:length(single_runs)) {
      raster_stacks <- stack(raster_stacks, single_runs[[j]][[q]])
    }
    simulation_mean <- raster::calc(raster_stacks, mean)
    simulation_sd <- raster::calc(raster_stacks, sd)
    simulation_mean_stack <- stack(simulation_mean_stack, simulation_mean)
    simulation_sd_stack <-stack(simulation_sd_stack, simulation_sd)
  }

  
  # simulation_mean <- raster::calc(raster_stacks, mean)
  # simulation_sd <- raster::calc(raster_stacks, sd)
  
  outputs <- list(probability, simulation_mean_stack, simulation_sd_stack, single_run_out, number_infecteds, infected_areas, west_rate, east_rate, south_rate, north_rate)
  names(outputs) <- c('probability', 'simulation_mean', 'simulation_sd', 'single_run_out', 'number_infecteds', 'infected_areas', 'west_rate', 'east_rate', 'south_rate', 'north_rate')
  
  return(outputs)
  
}

