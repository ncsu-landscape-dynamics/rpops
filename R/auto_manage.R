#' PoPS (Pest or Pathogen Spread) automated management
#' 
#' A dynamic species distribution model for pest or pathogen spread in forest or agricultural ecosystems. The model is process based
#' meaning that it uses understanding of the effect of weather on reproduction and survival of the pest/pathogen in order to forecast
#' spread of the pest/pathogen into the future. 
#'
#' @inheritParams pops
#' @param infected_files file path to the infected species files for the start of the simulation
#' @param num_iterations how many iterations do you want to run to allow the calibration to converge at least 10 
#' @param number_of_cores enter how many cores you want to use (default = NA). If not set uses the # of CPU cores - 1. must be an integer >= 1
#' @param cost_per_meter_sq the cost of treatment per square meter
#' @param budget the total budget to spend on management
#' @param buffer the size of the buffer to include around managed locations
#' @param treatment_priority how to prioritize which of many species to treat options are 'equal' or 'ranked' where equal selects locations based on guidelines from any of the species and rank selects from the species in ranked order and then 
#' @param treatment_rank binary 0 or 1 for the species that is the most important
#' @param selection_method the method for determining the management strategy to use. must be one of 'Foci', 'Border', or'Points'.
#' @param selection_priority how to prioritize locations for management must be one of "group size", "host", or "infected"
#' @param points used if selection_method is points
#' @param treatment_efficacy The overall efficacy of the treatment
#' @param species a list of the species names for naming outputs files must be the same length and infected_files
#' @param direction_first boolean to indicate where or not direction is the first priortity in sorting (if false first sorting priority goes to the selection_method) 
#' @param anthropogenic_kappa sets the strength of the anthropogenic direction in the von-mises distribution numeric value between 0.01 and 12
#' @param natural_kappa sets the strength of the natural direction in the von-mises distribution numeric value between 0.01 and 12
#' @param reproductive_rate number of spores or pest units produced by a single host under optimal weather conditions 
#' @param percent_natural_dispersal  what percentage of dispersal is natural range versus anthropogenic range value between 0 and 1
#' @param natural_distance_scale distance scale parameter for natural range dispersal kernel numeric value > 0 
#' @param anthropogenic_distance_scale distance scale parameter for anthropogenic range dispersal kernel numeric value > 0
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
#' infected_files <-  system.file("extdata", "SODexample", "initial_infection2001.tif", 
#' package = "PoPS")
#' host_file <- system.file("extdata", "SODexample", "host.tif", package = "PoPS")
#' total_plants_file <- system.file("extdata", "SODexample", "all_plants.tif", package = "PoPS")
#' temperature_coefficient_file <- system.file("extdata", "SODexample", "weather.tif", package = "
#' PoPS")
#' treatments_file <- system.file("extdata", "SODexample", "management.tif", package = "PoPS")
#' 
#' data <- pops(infected_files, host_file, total_plants_file, reproductive_rate = 1.0,
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
auto_manage <- function(infected_files, 
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
                        use_movements = FALSE,
                        cost_per_meter_sq = 1.37, 
                        budget = 1500000, 
                        buffer = 600,
                        treatment_priority = "equal", 
                        treatment_rank = c(1), 
                        selection_method = 'Points', 
                        selection_priority = 'group size',
                        points = points, 
                        treatment_efficacy = 1, 
                        species = c('species1'), 
                        direction_first = TRUE,
                        start_exposed = FALSE,
                        generate_stochasticity = TRUE,
                        establishment_stochasticity = TRUE,
                        movement_stochasticity = TRUE,
                        deterministic = FALSE,
                        establishment_probability = 0.5,
                        dispersal_percentage = 0.99) {
  
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
  
  multispecies_check <- multispecies_checks(species, infected_files, reproductive_rate, percent_natural_dispersal, natural_kernel_type, anthropogenic_kernel_type, 
                                  natural_distance_scale, anthropogenic_distance_scale, natural_dir, natural_kappa, anthropogenic_dir, anthropogenic_kappa)
  if (!multispecies_check$checks_passed){
    return(percent_check$failed_check)
  }
  
  infected_check <- initial_raster_checks(infected_files)
  if (infected_check$checks_passed) {
    infected <- infected_check$raster
    # if (raster::nlayers(infected) > 1) {
    #   infected <- output_from_raster_mean_and_sd(infected)
    # }
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
  
  infected_speci <- infected
  susceptible_speci <- stack()
  for (r in 1:length(infected_files)) {
    infected_name <-  paste("infected", r, sep = "")
    susceptible <- host - infected[[r]]
    susceptible[susceptible < 0] <- 0

    susceptible_speci <- stack(susceptible_speci, susceptible)
  }
  
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
  
  infected_species <- list(as.matrix(infected_speci[[1]]))
  susceptible_species <- list(as.matrix(susceptible_speci[[1]]))
  for(u in 2:nlayers(infected_speci)) {
    infected_species[[u]] <- as.matrix(infected_speci[[u]])
    susceptible_species[[u]] <- as.matrix(susceptible_speci[[u]])
  }
  
  if (use_lethal_temperature == TRUE) {
    temperature_check <- secondary_raster_checks(temperature_file, infected)
    if (temperature_check$checks_passed) {
      temperature_stack <- temperature_check$raster
    } else {
      return(temperature_check$failed_check)
    }
    
    temperature <- list(raster::as.matrix(temperature_stack[[1]]))
    for(o in 2:number_of_years) {
      temperature[[o]] <- raster::as.matrix(temperature_stack[[o]])
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
    for(h in 2:number_of_time_steps) {
      weather_coefficient[[h]] <- raster::as.matrix(weather_coefficient_stack[[h]])
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
  
  mortality_tracker <- infected[[1]]
  raster::values(mortality_tracker) <- 0
  
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
  exposed_list <- list(exposed)
  
  if (length(infected_files) > 1) {
    for (z in 2:length(infected_files)) {
      exposed_list[[z]] <- exposed
    }
  }
  
  if (model_type == "SEI" & start_exposed) {
    exposed[[latency_period + 1]] <- infected
    infected <- mortality_tracker
  }

  years <- seq(year(start_date), year(end_date), 1)
  rcl <- c(1, Inf, 1, 0, 0.99, NA)
  rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
  
  ## management module information
  num_cells <- round((budget/cost_per_meter_sq)/(ew_res*ns_res))
  buffer_cells <- buffer/ew_res
  years_simulated <- length(years)
  
  random_seeds <- round(stats::runif(num_iterations, 1, 1000000))

  treatment_speci <- raster()
  
  if (is.na(number_of_cores) || number_of_cores > parallel::detectCores()) {
    core_count <- parallel::detectCores() - 1
  } else {
    core_count <- number_of_cores
  }
  
  run_years <-   foreach(y = 1:years_simulated, .combine = rbind, .packages = c("raster", "PoPS", "foreach", "lubridate")) %do% {

    if (treatment_priority == "equal") {
      treatment_speci <- raster::stackApply(infected_speci, indices = rep(1, raster::nlayers(infected_speci)), fun = sum)
      print("works")
    } else if (treatment_priority == "ranked") {
      for (m in 1:raster::nlayers(infected_speci)) {
        if (treatment_rank[[m]]) {
          treatment_speci <- infected_speci[[m]]
        }
      }
    }
    print("start")
    treatment <- treatmentAuto(rasts = infected_speci, rasts2 = susceptible_speci, 
                               method = selection_method, priority = selection_priority,
                               number_of_locations = num_cells, points = points, 
                               treatment_efficacy = treatment_efficacy, 
                               buffer_cells = buffer_cells, direction_first = direction_first, 
                               treatment_rank = treatment_rank, treatment_priority = treatment_priority)
    
    treatment_dates <- paste(years[1], "-12", "-01", sep = "")
    treatment_maps <- list(as.matrix(treatment))
    management <- TRUE
    print("end_treatment")
    
    tests <-   foreach(i = 1:length(infected_files), .combine = rbind, .packages = c("raster", "PoPS", "foreach")) %do% {
      
      cl <- makeCluster(core_count)
      registerDoParallel(cl)
      
      infected_stack <- foreach(p = 1:num_iterations, .combine = rbind, .packages = c("raster", "PoPS", "foreach"), .export = ls(globalenv())) %dopar% {
        
        data <- pops_model(random_seed = random_seeds[p], 
                           use_lethal_temperature = use_lethal_temperature, 
                           lethal_temperature = lethal_temperature, 
                           lethal_temperature_month = lethal_temperature_month,
                           infected = infected_species[[i]],
                           exposed = exposed[[i]],
                           susceptible = susceptible_species[[i]],
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
                           reproductive_rate = reproductive_rate[[i]],
                           mortality_rate = mortality_rate, 
                           mortality_time_lag = mortality_time_lag,
                           season_month_start = season_month_start, 
                           season_month_end = season_month_end,
                           start_date = start_date, 
                           end_date = end_date,
                           treatment_method = treatment_method,
                           natural_kernel_type = natural_kernel_type[[i]], 
                           anthropogenic_kernel_type = anthropogenic_kernel_type[[i]], 
                           use_anthropogenic_kernel = use_anthropogenic_kernel, 
                           percent_natural_dispersal = percent_natural_dispersal[[i]],
                           natural_distance_scale = natural_distance_scale[[i]], 
                           anthropogenic_distance_scale = anthropogenic_distance_scale[[i]], 
                           natural_dir = natural_dir[[i]], natural_kappa = natural_kappa[[i]],
                           anthropogenic_dir = anthropogenic_dir[[i]], 
                           anthropogenic_kappa = anthropogenic_kappa[[i]],
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
        
        infected_runs <- raster::stack(lapply(1:length(data$infected), function(x) host))
        susceptible_runs <- raster::stack(lapply(1:length(data$infected), function(x) host))
        
        for (q in 1:raster::nlayers(infected_runs)) {
          infected_runs[[q]] <- data$infected[[q]]
          susceptible_runs[[q]] <- data$susceptible[[q]]
        }
        
        prob_runs <- raster::reclassify(infected_runs, rclmat)
        prob_runs[is.na(prob_runs)] <- 0
        number_infected <- data$number_infected
        spread_rate <- data$rates
        infected_area <- data$area_infected
        to.infected_stack <- list(infected_runs, susceptible_runs, number_infected, infected_area, spread_rate, prob_runs)
      }
      
      stopCluster(cl)
      infected_runs <- infected_stack[1:10]
      susceptible_runs <- infected_stack[11:20]
      number_infected_runs <- infected_stack[21:30]
      area_infected_runs <- infected_stack[31:40]
      spread_rate_runs <- infected_stack[41:50]
      probability_runs <- infected_stack[51:60]
      
      prediction <- probability_runs[[1]]
      prediction[prediction > 0] <- 0
      infected_area <- data.frame(t(years))
      infected_number <- data.frame(t(years))
      west_rates <- data.frame(t(years))
      east_rates <- data.frame(t(years))
      south_rates <- data.frame(t(years))
      north_rates <- data.frame(t(years))
      
      for (k in 1:length(infected_runs)) {
        prediction <- prediction + probability_runs[[k]]
        infected_number[k,] <- number_infected_runs[[k]]
        infected_area[k,] <- area_infected_runs[[k]]
        rates <- do.call(rbind, spread_rate_runs[[k]])
        west_rates[k,] <- rates[,4]
        east_rates[k,] <- rates[,3]
        south_rates[k,] <- rates[,2]
        north_rates[k,] <- rates[,1]
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
      
      infected_run <- infected_runs[[median_run_index]]
      susceptible_run <- susceptible_runs[[median_run_index]]
      
      to.test <- list(infected_run, susceptible_run, probability, number_infecteds, infected_areas, west_rate, east_rate, north_rate, south_rate)
    }
    
    if (y == 1) {
      infected_speci <- stack()
      susceptible_speci <- stack()
      probabilities <- stack()
      
      infections_out <- c()
      susceptibles_out <- c()
      probabilities_out <- c()
      number_infecteds_out <- c()
      infected_areas_out <- c()
      west_rate_out <- c()
      east_rate_out <- c()
      north_rate_out <- c()
      south_rate_out <- c()
    } else if (y > 1) {
      infections_out <- c()
      susceptibles_out <- c()
      probabilities_out <- c()
      number_infecteds_out <- c()
      infected_areas_out <- c()
      west_rate_out <- c()
      east_rate_out <- c()
      north_rate_out <- c()
      south_rate_out <- c()
      
      infected_speci <- stack()
      susceptible_speci <- stack()
      probabilities <- stack()
    }
    
    for (t in 1:length(infected_species)) {
      infected_speci <- stack(infected_speci, tests[[t]][[1]])
      susceptible_speci <- stack(susceptible_speci, tests[[t+length(infected_species)]][[1]])
      probabilities <- stack(probabilities, tests[[t+(2*length(infected_species))]][[1]])
      
      infections_out[[t]] <- c(tests[[t]])
      susceptibles_out[[t]] <- c(tests[[t+length(infected_species)]])
      probabilities_out[[t]] <- c(tests[[t+(2*length(infected_species))]])
      number_infecteds_out[[t]] <- c(tests[[t+(3*length(infected_species))]])
      infected_areas_out[[t]] <- c(tests[[t+(4*length(infected_species))]])
      west_rate_out[[t]] <- c(tests[[t+(5*length(infected_species))]])
      east_rate_out[[t]] <- c(tests[[t+(6*length(infected_species))]])
      north_rate_out[[t]] <- c(tests[[t+(7*length(infected_species))]])
      south_rate_out[[t]] <- c(tests[[t+(8*length(infected_species))]])
      
      infected_species[[t]] <- as.matrix(tests[[t]][[1]])
      susceptible_species[[t]] <- as.matrix(tests[[t+length(infected_species)]][[1]])
    }
    
    start_date <- paste(year(start_date) + 1, "-01", "-01", sep = "")
    end_date <- end_date 
    if (y < years_simulated) {
      years <- seq(year(start_date), year(end_date), 1)
    } else {
      years <- c(year(start_date))
    }
    print("outer loop")
    
    outputs <- list(infections_out, susceptibles_out, probabilities_out, number_infecteds_out, infected_areas_out, west_rate_out, east_rate_out, north_rate_out, south_rate_out, treatment)
    
  }
  
  year_names <- list()
  
  infecteds <- list()
  susceptibles <- list()
  probabilities <- list()
  number_infecteds <- list()
  infected_areas <- list()
  west_rates <- list()
  east_rates <- list()
  north_rates <- list()
  south_rates <- list()
  treatments <- list()
  
  for (sp in 1:length(infected_files)) {
    infecs <- list()
    susceps <- list()
    probs <- list()
    number_infects <- list()
    infected_ars <- list()
    west_rats <- list()
    east_rats <- list()
    north_rats <- list()
    south_rats <- list()
    for (l in 1:years_simulated) {
      year_names[l] <- paste(years - years_simulated - 1 + l)
      infecs[[l]] <- run_years[[l]][[sp]][[1]]
      susceps[[l]] <- run_years[[l + years_simulated]][[sp]][[1]]
      probs[[l]] <- run_years[[l + years_simulated * 2]][[sp]][[1]]
      number_infects[[l]] <- run_years[[l + years_simulated * 3]][[sp]][[1]]
      infected_ars[[l]] <- run_years[[l + years_simulated * 4]][[sp]][[1]]
      west_rats[[l]] <- run_years[[l + years_simulated * 5]][[sp]][[1]]
      east_rats[[l]] <- run_years[[l + years_simulated * 6]][[sp]][[1]]
      north_rats[[l]] <- run_years[[l + years_simulated * 7]][[sp]][[1]]
      south_rats[[l]] <- run_years[[l + years_simulated * 8]][[sp]][[1]]
      treatments[[l]] <- run_years[[l + years_simulated * 9]][[1]]
    }
    names(infecs) <- year_names
    names(susceps) <- year_names
    names(probs) <- year_names
    names(number_infects) <- year_names
    names(infected_ars) <- year_names
    names(west_rats) <- year_names
    names(east_rats) <- year_names
    names(north_rats) <- year_names
    names(south_rats) <- year_names
    names(treatments) <- year_names
    infecteds[[sp]] <- infecs
    susceptibles[[sp]] <- susceps
    probabilities[[sp]] <- probs
    number_infecteds[[sp]] <- number_infects
    infected_areas[[sp]] <- infected_ars
    west_rates[[sp]] <- west_rats
    east_rates[[sp]] <- east_rats
    north_rates[[sp]] <- north_rats
    south_rates[[sp]] <- south_rats
  }
  
  names(infecteds) <- species
  names(susceptibles) <- species
  names(probabilities) <- species
  names(number_infecteds) <- species
  names(infected_areas) <- species
  names(west_rates) <- species
  names(east_rates) <- species
  names(north_rates) <- species
  names(south_rates) <- species
  
  outputs <- list(infecteds, susceptibles, probabilities, number_infecteds, infected_areas, west_rates, east_rates, north_rates, south_rates, treatments)
  names(outputs) <- c('infecteds', 'susceptibles', 'probabilities', 'number_infecteds', 'infected_areas', 'west_rates', 'east_rates', 'north_rates', 'south_rates', 'treatments')
  return(outputs)
}
  