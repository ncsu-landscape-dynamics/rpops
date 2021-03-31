#' PoPS (Pest or Pathogen Spread) automated management
#'
#' This function allows for treating multiple species with the same management
#' strategy. Was designed for multiple lineages of the same species that share
#' common hosts and are in the same study area. This functions loop straight
#' through each stochastic realisation with the same treatment.
#'
#' @inheritParams pops
#' @param infected_files file path to the infected species files for the start
#' of the simulation
#' @param number_of_iterations how many iterations do you want to run to allow the
#' calibration to converge at least 10
#' @param number_of_cores enter how many cores you want to use (default = NA).
#' If not set uses the # of CPU cores - 1. must be an integer >= 1
#' @param cost_per_meter_sq the cost of treatment per square meter
#' @param budget the total budget to spend on management
#' @param buffer the size of the buffer to include around managed locations
#' @param treatment_priority how to prioritize which of many species to treat
#' options are 'equal' or 'ranked' where equal selects locations based on
#' guidelines from any of the species and rank selects from the species in
#' ranked order and then
#' @param treatment_rank binary 0 or 1 for the species that is the most
#' important
#' @param selection_method the method for determining the management strategy
#' to use. must be one of 'Foci', 'Border', or'Points'.
#' @param selection_priority how to prioritize locations for management must be
#' one of "group size", "host", or "infected"
#' @param points used if selection_method is points
#' @param treatment_efficacy The overall efficacy of the treatment
#' @param species a list of the species names for naming outputs files must be
#' the same length and infected_files
#' @param direction_first boolean to indicate if direction is the
#' first priortity in sorting (if false first sorting priority goes to the
#' selection_method)
#'
#' @importFrom terra global rast xres yres classify extract ext as.points ncol
#' nrow nlyr rowFromCell colFromCell values as.matrix rowFromCell colFromCell
#' crs app patches
#' @importFrom stats runif rnorm median sd
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar%
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom iterators icount
#' @importFrom lubridate interval time_length year
#' @return list of infected and susceptible per year
#' @export
#'
auto_manage <- function(infected_files,
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
                        treatment_dates = c(""),
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
                        dispersal_percentage = 0.99,
                        quarantine_areas_file = "",
                        use_quarantine = FALSE,
                        use_spreadrates = FALSE,
                        use_overpopulation_movements = FALSE,
                        overpopulation_percentage = 0,
                        leaving_percentage = 0,
                        leaving_scale_coefficient = 1) {

  config <- c()
  config$random_seed <- random_seed
  config$infected_files <- infected_files
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
  config$deterministic <- deterministic
  config$establishment_probability <- establishment_probability
  config$dispersal_percentage <- dispersal_percentage
  config$quarantine_areas_file <- quarantine_areas_file
  config$use_quarantine <- use_quarantine
  config$use_spreadrates <- use_spreadrates
  config$use_overpopulation_movements <- use_overpopulation_movements
  config$overpopulation_percentage <- overpopulation_percentage
  config$leaving_percentage <- leaving_percentage
  config$leaving_scale_coefficient <- leaving_scale_coefficient
  config$number_of_iterations <- number_of_iterations
  config$number_of_cores <- number_of_cores
  config$species <- species
  # add function name for use in configuration function to skip
  # function specific specific configurations namely for validation and
  # calibration.
  config$function_name <- "auto-manage"
  config$failure <- NULL

  config <- configuration(config)

  if (!is.null(config$failure)) {
    return(config$failure)
  }

  infected_speci <- infected
  susceptible_speci <- stack()
  for (r in 1:length(infected_files)) {
    infected_name <-  paste("infected", r, sep = "")
    susceptible <- host - infected[[r]]
    susceptible[susceptible < 0] <- 0

    susceptible_speci <- stack(susceptible_speci, susceptible)
  }

  infected_species <- list(as.matrix(infected_speci[[1]]))
  susceptible_species <- list(as.matrix(susceptible_speci[[1]]))
  for(u in 2:nlayers(infected_speci)) {
    infected_species[[u]] <- as.matrix(infected_speci[[u]])
    susceptible_species[[u]] <- as.matrix(susceptible_speci[[u]])
  }

  config$random_seeds <-
    round(stats::runif(config$number_of_iterations, 1, 1000000))

  # treatment_speci <- raster()

  i <- p <- y <- NULL

  if (is.na(number_of_cores) || number_of_cores > parallel::detectCores()) {
    core_count <- parallel::detectCores() - 1
  } else {
    core_count <- number_of_cores
  }

  run_years <-   foreach(y = 1:config$number_of_outputs, .combine = rbind, .packages = c("raster", "PoPS", "foreach", "lubridate")) %do% {

    # if (treatment_priority == "equal") {
    #   treatment_speci <- raster::stackApply(infected_speci, indices = rep(1, raster::nlayers(infected_speci)), fun = sum)
    #   print("works")
    # } else if (treatment_priority == "ranked") {
    #   for (m in 1:raster::nlayers(infected_speci)) {
    #     if (treatment_rank[[m]]) {
    #       treatment_speci <- infected_speci[[m]]
    #     }
    #   }
    # }
    print("start")
    treatment <- treatment_auto(rasts = infected_speci, rasts2 = susceptible_speci,
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

      infected_stack <- foreach(p = 1:number_of_iterations, .combine = rbind, .packages = c("raster", "PoPS", "foreach"), .export = ls(globalenv())) %dopar% {

        data <- pops_model(random_seed = random_seeds[p],
                           use_lethal_temperature = use_lethal_temperature,
                           lethal_temperature = lethal_temperature,
                           lethal_temperature_month = lethal_temperature_month,
                           infected = infected_species[[i]],
                           exposed = exposed[[i]],
                           susceptible = susceptible_species[[i]],
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
                           reproductive_rate = reproductive_rate[[i]],
                           spatial_indices = config$spatial_indices,
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
                           dispersal_percentage = dispersal_percentage,
                           use_overpopulation_movements = use_overpopulation_movements,
                           overpopulation_percentage = overpopulation_percentage,
                           leaving_percentage = leaving_percentage,
                           leaving_scale_coefficient = leaving_scale_coefficient
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
    if (y < config$number_of_outputs) {
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
    for (l in 1:config$number_of_outputs) {
      year_names[l] <- paste(years - config$number_of_outputs - 1 + l)
      infecs[[l]] <- run_years[[l]][[sp]][[1]]
      susceps[[l]] <- run_years[[l + config$number_of_outputs]][[sp]][[1]]
      probs[[l]] <- run_years[[l + config$number_of_outputs * 2]][[sp]][[1]]
      number_infects[[l]] <- run_years[[l + config$number_of_outputs * 3]][[sp]][[1]]
      infected_ars[[l]] <- run_years[[l + config$number_of_outputs * 4]][[sp]][[1]]
      west_rats[[l]] <- run_years[[l + config$number_of_outputs * 5]][[sp]][[1]]
      east_rats[[l]] <- run_years[[l + config$number_of_outputs * 6]][[sp]][[1]]
      north_rats[[l]] <- run_years[[l + config$number_of_outputs * 7]][[sp]][[1]]
      south_rats[[l]] <- run_years[[l + config$number_of_outputs * 8]][[sp]][[1]]
      treatments[[l]] <- run_years[[l + config$number_of_outputs * 9]][[1]]
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
