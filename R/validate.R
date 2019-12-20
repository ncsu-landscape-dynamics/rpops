#' Validates the accuracy of the calibrated reproductive rate and dispersal scales of the pops model.
#' 
#' This function uses the quantity, allocation, and configuration disagreement to validate the model
#' across the landscape using the parameters from the calibrate function. Ideally the model is calibrated
#' with 2 or more years of data and validated for the last year or if you have 6 or more years of data 
#' then the model can be validated for the final 2 years.
#'
#' @inheritParams pops
#' @param infected_years_file years of initial infection/infestation as individual locations of a pest or pathogen in raster format
#' @param num_iterations how many iterations do you want to run to allow the calibration to converge at least 10 
#' @param number_of_cores enter how many cores you want to use (default = NA). If not set uses the # of CPU cores - 1. must be an integer >= 1
#' @param success_metric Choose which success metric to use for calibration. Choices are "quantity", "quantity and configuration", "residual error" and "odds ratio". Default is "quantity"
#' @param mask Raster file used to provide a mask to remove 0's that are not true negatives from comparisons (e.g. mask out lakes and oceans from statics if modeling terrestrial species).
#'
#' @importFrom raster raster values as.matrix xres yres stack reclassify cellStats nlayers calc extract rasterToPoints
#' @importFrom stats runif rnorm
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar%
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom iterators icount
#' @importFrom lubridate interval time_length
#' 
#' @return a dataframe of the variables saved and their success metrics for each run
#' @export
#'
#' @examples
#' \dontrun{
#' infected_years_file <- system.file("extdata", "SODexample", "initial_infection2005.tif", 
#' package = "PoPS")
#' num_iterations <- 100
#' number_of_cores <- NA
#' 
#' infected_file <- system.file("extdata", "SODexample", "initial_infection2004.tif", 
#' package = "PoPS")
#' host_file <- system.file("extdata", "SODexample", "host.tif", package = "PoPS")
#' total_plants_file <- system.file("extdata", "SODexample", "all_plants.tif", package = "PoPS")
#' temperature_coefficient_file <- system.file("extdata", "SODexample", "weather.tif", 
#' package = "PoPS")
#' treatments_file <- system.file("extdata", "SODexample", "management2005.tif", package = "PoPS")
#' 
#' params <- validate(infected_years_file, num_iterations, number_of_cores,
#' infected_file, host_file, total_plants_file, reproductive_rate = 1.0,
#' use_lethal_temperature = FALSE, temp = TRUE, precip = FALSE, management = TRUE, 
#' mortality_on = TRUE, temperature_file = "", temperature_coefficient_file, 
#' precipitation_coefficient_file ="", treatments_file,
#' season_month_start = 1, season_month_end = 12, time_step = "month",
#' start_date = '2005-01-01', end_date = 2005-12-31', treatment_dates = c(2005),
#' dispersal_kern = "cauchy", percent_short_distance_dispersal = 1.0,
#' short_distance_scale = 20.57, long_distance_scale = 0.0,
#' lethal_temperature = -12.87, lethal_temperature_month = 1,
#' mortality_rate = 0.05, mortality_time_lag = 2,
#' wind_dir = "NONE", kappa = 0)
#' }
#' 
validate <- function(infected_years_file, num_iterations, number_of_cores = NA,
                     infected_file, host_file, total_plants_file, 
                     temp = FALSE, temperature_coefficient_file = "", 
                     precip = FALSE, precipitation_coefficient_file = "", 
                     time_step = "month", reproductive_rate = 3.0,
                     season_month_start = 1, season_month_end = 12, 
                     start_date = '2008-01-01', end_date = '2008-12-31',  
                     use_lethal_temperature = FALSE, temperature_file = "",
                     lethal_temperature = -12.87, lethal_temperature_month = 1,
                     mortality_on = FALSE, mortality_rate = 0, mortality_time_lag = 0, 
                     management = FALSE, treatment_dates = c(0), treatments_file = "",
                     treatment_method = "ratio",
                     percent_natural_dispersal = 1.0,
                     natural_kernel_type = "cauchy", anthropogenic_kernel_type = "cauchy",
                     natural_distance_scale = 21, anthropogenic_distance_scale = 0.0,
                     natural_dir = "NONE", natural_kappa = 0, 
                     anthropogenic_dir = "NONE", anthropogenic_kappa = 0, 
                     pesticide_duration = 0, pesticide_efficacy = 1.0,
                     mask = NULL, success_metric = "quantity", output_frequency = "year"
                     ){ 
  
  metric_check <- metric_checks(success_metric)
  if (metric_check$checks_passed){
    configuration <- metric_check$configuration
  } else {
    return(metric_check$failed_check)
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
  
  # reference <- raster(infected_years_file)
  ## Load observed data on occurence
  infection_years <- stack(infected_years_file)
  ## set up reclassification matrix for binary reclassification
  rcl <- c(1, Inf, 1, 0, 0.99, NA)
  rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
  ## reclassify to binary values
  infection_years <- reclassify(infection_years, rclmat)
  ## Get rid of NA values to make comparisons
  infection_years <- raster::reclassify(infection_years, matrix(c(NA,0), ncol = 2, byrow = TRUE), right = NA)
  
  num_layers_infected_years <- raster::nlayers(infection_years)
  if (num_layers_infected_years < number_of_outputs) {
    return(paste("The infection years file must have enough layers to match the number of outputs from the model. The number of layers of your infected year file is", num_layers_infected_years, "and the number of outputs is", number_of_time_steps))
  }
  
  if (is.na(number_of_cores) || number_of_cores > parallel::detectCores()) {
    core_count <- parallel::detectCores() - 1
  } else {
    core_count <- number_of_cores
  }
  cl <- makeCluster(core_count)
  registerDoParallel(cl)
  

  qa <- foreach::foreach (icount(num_iterations), .combine = rbind, .packages = c("raster", "PoPS", "foreach")) %dopar% {
    random_seed <- round(stats::runif(1, 1, 1000000))
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
                             time_step = time_step, reproductive_rate = reproductive_rate[i],
                             mortality_rate = mortality_rate, mortality_time_lag = mortality_time_lag,
                             season_month_start = season_month_start, season_month_end = season_month_end,
                             start_date = start_date, end_date = end_date,
                             treatment_method = treatment_method,
                             natural_kernel_type = natural_kernel_type, anthropogenic_kernel_type = anthropogenic_kernel_type, 
                             use_anthropogenic_kernel = use_anthropogenic_kernel, percent_natural_dispersal = percent_natural_dispersal[i],
                             natural_distance_scale = natural_distance_scale[i], anthropogenic_distance_scale = anthropogenic_distance_scale[i], 
                             natural_dir = natural_dir, natural_kappa = natural_kappa,
                             anthropogenic_dir = anthropogenic_dir, anthropogenic_kappa = anthropogenic_kappa,
                             output_frequency = output_frequency
    )
    
    comp_year <- raster(infected_file)
    all_disagreement <- foreach(q = 1:length(data$infected), .combine = rbind, .packages =c("raster", "PoPS", "foreach"), .final = colSums) %dopar% {
      comp_year[] <- data$infected[[q]]
      comp_year <- reclassify(comp_year, rclmat)
      to.all_disagreement <- quantity_allocation_disagreement(infection_years[[q]], comp_year, configuration, mask)
    }
    
    # comparison <- reference
    # raster::values(comparison) <- data$infected[[1]]
    # 
    # comparison <- raster::reclassify(comparison, rclmat)
    
    # to.qa <- PoPS::quantity_allocation_disagreement(reference, comparison)
    to.qa <- data.frame(t(all_disagreement))
  }
  
  parallel::stopCluster(cl)
  
  return(qa)
}
