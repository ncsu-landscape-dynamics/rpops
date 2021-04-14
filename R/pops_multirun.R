#' PoPS (Pest or Pathogen Spread) model Multiple Runs
#'
#' A dynamic species distribution model for pest or pathogen spread in forest
#' or agricultural ecosystems. The model is process based meaning that it uses
#' understanding of the effect of weather and other environmental factors on
#' reproduction and survival of the pest/pathogen in order to forecast spread
#' of the pest/pathogen into the future. Run multiple stochasitic simulations,
#' propogating uncertainty in parameters, initial conditions, and drivers.
#' The model is process based meaning that it uses understanding of the effect
#' of weather on reproduction and survival of the pest/pathogen in order to
#' forecast spread of the pest/pathogen into the future.
#'
#' @inheritParams pops
#' @param number_of_iterations how many iterations do you want to run to allow
#' the calibration to converge at least 10
#' @param number_of_cores enter how many cores you want to use (default = NA).
#' If not set uses the # of CPU cores - 1. must be an integer >= 1
#'
#' @importFrom terra app rast xres yres classify extract ext as.points ncol nrow
#' nlyr rowFromCell colFromCell values as.matrix rowFromCell colFromCell crs
#' @importFrom stats runif rnorm median sd
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar% %do%
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom lubridate interval time_length mdy %within%
#' @return list of infected and susceptible per year
#' @export
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
                          start_date = "2008-01-01",
                          end_date = "2008-12-31",
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
                          leaving_scale_coefficient = 1,
                          exposed_file = "",
                          mask = NULL) {

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
  # add function name for use in configuration function to skip
  # function specific specifc configurations namely for validation and
  # calibration.
  config$function_name <- "multirun"
  config$failure <- NULL
  config$exposed_file <- exposed_file
  config$mask <- mask

  config <- configuration(config)

  if (!is.null(config$failure)) {
    return(config$failure)
  }

  config$crs <- terra::crs(config$host)
  i <- NULL

  cl <- makeCluster(config$core_count)
  registerDoParallel(cl)

  infected_stack <-
    foreach::foreach(i = seq_len(config$number_of_iterations),
                     .combine = c,
                     .packages = c("PoPS", "terra")) %dopar% {

    config$random_seed <- round(stats::runif(1, 1, 1000000))
    data <- PoPS::pops_model(random_seed = config$random_seed,
                       use_lethal_temperature = config$use_lethal_temperature,
                       lethal_temperature = config$lethal_temperature,
                       lethal_temperature_month =
                         config$lethal_temperature_month,
                       infected = config$infected,
                       exposed = config$exposed,
                       susceptible = config$susceptible,
                       total_populations = config$total_populations,
                       mortality_on = config$mortality_on,
                       mortality_tracker = config$mortality_tracker,
                       mortality = config$mortality,
                       quarantine_areas = config$quarantine_areas,
                       treatment_maps = config$treatment_maps,
                       treatment_dates = config$treatment_dates,
                       pesticide_duration = config$pesticide_duration,
                       resistant = config$resistant,
                       use_movements = config$use_movements,
                       movements = config$movements,
                       movements_dates = config$movements_dates,
                       weather = config$weather,
                       temperature = config$temperature,
                       weather_coefficient = config$weather_coefficient,
                       ew_res = config$ew_res,
                       ns_res = config$ns_res,
                       num_rows = config$num_rows,
                       num_cols = config$num_cols,
                       time_step = config$time_step,
                       reproductive_rate = config$reproductive_rate[i],
                       spatial_indices = config$spatial_indices,
                       mortality_rate = config$mortality_rate,
                       mortality_time_lag = config$mortality_time_lag,
                       season_month_start = config$season_month_start,
                       season_month_end = config$season_month_end,
                       start_date = config$start_date,
                       end_date = config$end_date,
                       treatment_method = config$treatment_method,
                       natural_kernel_type = config$natural_kernel_type,
                       anthropogenic_kernel_type =
                         config$anthropogenic_kernel_type,
                       use_anthropogenic_kernel =
                         config$use_anthropogenic_kernel,
                       percent_natural_dispersal =
                         config$percent_natural_dispersal[i],
                       natural_distance_scale =
                         config$natural_distance_scale[i],
                       anthropogenic_distance_scale =
                         config$anthropogenic_distance_scale[i],
                       natural_dir = config$natural_dir,
                       natural_kappa = config$natural_kappa[i],
                       anthropogenic_dir = config$anthropogenic_dir,
                       anthropogenic_kappa = config$anthropogenic_kappa[i],
                       output_frequency = config$output_frequency,
                       output_frequency_n = config$output_frequency_n,
                       quarantine_frequency = config$quarantine_frequency,
                       quarantine_frequency_n = config$quarantine_frequency_n,
                       use_quarantine = config$use_quarantine,
                       spreadrate_frequency = config$spreadrate_frequency,
                       spreadrate_frequency_n = config$spreadrate_frequency_n,
                       use_spreadrates = config$use_spreadrates,
                       model_type_ = config$model_type,
                       latency_period = config$latency_period,
                       generate_stochasticity =
                         config$generate_stochasticity,
                       establishment_stochasticity =
                         config$establishment_stochasticity,
                       movement_stochasticity = config$movement_stochasticity,
                       deterministic = config$deterministic,
                       establishment_probability =
                         config$establishment_probability,
                       dispersal_percentage = config$dispersal_percentage,
                       use_overpopulation_movements = config$use_overpopulation_movements,
                       overpopulation_percentage = config$overpopulation_percentage,
                       leaving_percentage = config$leaving_percentage,
                       leaving_scale_coefficient = config$leaving_scale_coefficient
    )

    run <- c()
    run$single_run <- data$infected
    run$comp_years <- data$infected
    run$number_infected <- data$number_infected
    run$susceptible_runs <- data$susceptible
    run$infected_area <- data$area_infected
    run$spread_rate <- data$rates
    run$quarantine_escape <- data$quarantine_escape
    run$quarantine_escape_distance <- data$quarantine_escape_distance
    run$quarantine_escape_direction <- data$quarantine_escape_directions
    run$exposed_runs <- data$exposed

    run

  }

  stopCluster(cl)
  single_runs <- infected_stack[seq(1, length(infected_stack), 10)]
  probability_runs <- infected_stack[seq(2, length(infected_stack), 10)]
  number_infected_runs <- infected_stack[seq(3, length(infected_stack), 10)]
  susceptible_runs <- infected_stack[seq(4, length(infected_stack), 10)]
  area_infected_runs <- infected_stack[seq(5, length(infected_stack), 10)]
  spread_rate_runs <- infected_stack[seq(6, length(infected_stack), 10)]
  quarantine_escape_runs <-
    infected_stack[seq(7, length(infected_stack), 10)]
  quarantine_escape_distance_runs <-
    infected_stack[seq(8, length(infected_stack), 10)]
  quarantine_escape_directions_runs <-
    infected_stack[seq(9, length(infected_stack), 10)]
  exposed_runs <- infected_stack[seq(10, length(infected_stack), 10)]

  prediction <- probability_runs[[1]]
  for (w in seq_len(length(prediction))) {
    prediction[[w]] <- 0
  }
  escape_probability <-
    data.frame(t(rep(0, length(probability_runs[[1]]))))
  infected_area <- data.frame(t(rep(0, length(probability_runs[[1]]))))
  infected_number <- data.frame(t(rep(0, length(probability_runs[[1]]))))
  west_rates <- data.frame(t(rep(0, length(probability_runs[[1]]))))
  east_rates <- data.frame(t(rep(0, length(probability_runs[[1]]))))
  south_rates <- data.frame(t(rep(0, length(probability_runs[[1]]))))
  north_rates <- data.frame(t(rep(0, length(probability_runs[[1]]))))
  max_values <- data.frame(t(rep(0, length(probability_runs[[1]]))))
  quarantine_escapes <-
    data.frame(t(rep(0, length(probability_runs[[1]]))))
  quarantine_escape_distances <-
    data.frame(t(rep(0, length(probability_runs[[1]]))))
  quarantine_escape_directions <-
    data.frame(t(rep(0, length(probability_runs[[1]]))))

  for (p in seq_len(length(probability_runs))) {
    for (w in seq_len(length(prediction))) {
      prob <- probability_runs[[p]][[w]]
      max_values[p, w] <- max(prob)
      prob[prob <= 1] <- 0
      prob[prob > 1] <- 1
      prediction[[w]] <- prediction[[w]] + prob
    }
    infected_number[p, ] <- number_infected_runs[[p]]
    infected_area[p, ] <- area_infected_runs[[p]]
    rates <- do.call(rbind, spread_rate_runs[[p]])
    if (!is.null(rates)) {
      west_rates[p, ] <- rates[, 4]
      east_rates[p, ] <- rates[, 3]
      south_rates[p, ] <- rates[, 2]
      north_rates[p, ] <- rates[, 1]
    }else {
      west_rates[p, ] <- 0
      east_rates[p, ] <- 0
      south_rates[p, ] <- 0
      north_rates[p, ] <- 0
    }

    if (config$use_quarantine &
        length(quarantine_escape_runs[[p]]) ==
        length(probability_runs[[p]])) {
      escape_probability <- escape_probability + quarantine_escape_runs[[p]]
      quarantine_escapes[p, ] <- quarantine_escape_runs[[p]]
      quarantine_escape_distances <- quarantine_escape_distance_runs[[p]]
      quarantine_escape_directions <- quarantine_escape_directions_runs[[p]]
    }
  }

  probability <- prediction
  for (w in seq_len(length(prediction))) {
    probability[[w]] <- (prediction[[w]] / (length(probability_runs))) * 100
  }

  infected_areas <-
    round(sapply(infected_area, function(x) c("Mean" = mean(x, na.rm = TRUE),
                                              "Stand dev" = sd(x))), digits = 0)
  number_infecteds <-
    round(sapply(infected_number, function(x) c("Mean" = mean(x, na.rm = TRUE),
                                                "Stand dev" = sd(x))),
          digits = 0)
  west_rate <-
    round(sapply(west_rates, function(x) c("Mean" = mean(x, na.rm = TRUE),
                                           "Stand dev" = sd(x))), digits = 0)
  east_rate <-
    round(sapply(east_rates, function(x) c("Mean" = mean(x, na.rm = TRUE),
                                           "Stand dev" = sd(x))), digits = 0)
  south_rate <-
    round(sapply(south_rates, function(x) c("Mean" = mean(x, na.rm = TRUE),
                                            "Stand dev" = sd(x))), digits = 0)
  north_rate <-
    round(sapply(north_rates, function(x) c("Mean" = mean(x, na.rm = TRUE),
                                            "Stand dev" = sd(x))), digits = 0)

  west_rate[is.na(west_rate)] <- 0
  east_rate[is.na(east_rate)] <- 0
  south_rate[is.na(south_rate)] <- 0
  north_rate[is.na(north_rate)] <- 0

  ## add quarantine here
  if (use_quarantine) {
    escape_probability <- escape_probability / length(probability_runs) * 100
    if (
      length(quarantine_escape_distances[quarantine_escape_directions == "N"]) >
      0) {
      north_distance_to_quarantine <-
        round(sapply(
          quarantine_escape_distances[quarantine_escape_directions == "N"],
          function(x) c("Mean" = mean(x, na.rm = TRUE),
                        "Stand dev" = sd(x))), digits = 0)
    } else {
      north_distance_to_quarantine <-
        data.frame(t(rep(NA, length(probability_runs[[1]]))))
    }

    if (
      length(quarantine_escape_distances[quarantine_escape_directions == "S"]) >
      0) {
      south_distance_to_quarantine <-
        round(sapply(
          quarantine_escape_distances[quarantine_escape_directions == "S"],
          function(x) c("Mean" = mean(x, na.rm = TRUE),
                        "Stand dev" = sd(x, na.rm = TRUE))), digits = 0)
    } else {
      south_distance_to_quarantine <-
        data.frame(t(rep(NA, length(probability_runs[[1]]))))
    }

    if (
      length(quarantine_escape_distances[quarantine_escape_directions == "E"]) >
      0) {
      east_distance_to_quarantine <-
        round(sapply(
          quarantine_escape_distances[quarantine_escape_directions == "E"],
          function(x) c("Mean" = mean(x, na.rm = TRUE),
                        "Stand dev" = sd(x))), digits = 0)
    } else {
      east_distance_to_quarantine <-
        data.frame(t(rep(NA, length(probability_runs[[1]]))))
    }

    if (
      length(quarantine_escape_distances[quarantine_escape_directions == "W"]) >
      0) {
      west_distance_to_quarantine <-
        round(sapply(
          quarantine_escape_distances[quarantine_escape_directions == "W"],
          function(x) c("Mean" = mean(x, na.rm = TRUE),
                        "Stand dev" = sd(x))), digits = 0)
    } else {
      west_distance_to_quarantine <-
        data.frame(t(rep(NA, length(probability_runs[[1]]))))
    }

  } else {
    escape_probability <-
      data.frame(t(rep(NA, length(probability_runs[[1]]))))
    north_distance_to_quarantine <-
      data.frame(t(rep(NA, length(probability_runs[[1]]))))
    south_distance_to_quarantine <-
      data.frame(t(rep(NA, length(probability_runs[[1]]))))
    east_distance_to_quarantine <-
      data.frame(t(rep(NA, length(probability_runs[[1]]))))
    west_distance_to_quarantine <-
      data.frame(t(rep(NA, length(probability_runs[[1]]))))
  }

  which_median <- function(x) raster::which.min(abs(x - median(x)))

  median_run_index <- which_median(infected_number[[ncol(infected_number)]])
  min_run_index <- which.min(infected_number[[ncol(infected_number)]])
  max_run_index <- which.max(infected_number[[ncol(infected_number)]])

  single_run <- single_runs[[median_run_index]]
  susceptible_run <- susceptible_runs[[median_run_index]]
  exposed_run <- exposed_runs[[median_run_index]]

  min_run <- single_runs[[min_run_index]]
  max_run <- single_runs[[max_run_index]]

  for (q in seq_len(length(single_runs[[1]]))) {
    for (j in seq_len(length(single_runs))) {
      if (j == 1) {
        raster_stacks <- list(single_runs[[j]][[q]])
      } else {
        raster_stacks[[j]] <- single_runs[[j]][[q]]
      }
    }
    raster_stacks2 <- do.call(cbind, raster_stacks)
    raster_stacks2 <-
      array(raster_stacks2, dim = c(dim(raster_stacks[[1]]),
                                  length(raster_stacks)))
    sim_mean <-
      round(apply(raster_stacks2, c(1, 2), mean, na.rm = TRUE), digits = 0)
    sim_sd <- apply(raster_stacks2, c(1, 2), sd, na.rm = TRUE)

    simulation_mean <-
      terra::rast(nrow = config$num_rows, ncol = config$num_cols,
                  xmin = config$xmin, xmax = config$xmax,
                  ymin = config$ymin, ymax = config$ymax, crs = config$crs)
    simulation_sd <- simulation_mean
    simulation_max <- simulation_mean
    simulation_min <- simulation_mean
    simulation_probability <- simulation_mean
    simulation_median <- simulation_mean
    simulation_susceptible <- simulation_mean
    terra::values(simulation_mean) <- sim_mean
    names(simulation_mean) <- "mean"

    terra::values(simulation_sd) <- sim_sd
    names(simulation_sd) <- "sd"

    terra::values(simulation_max) <- max_run[[q]]
    names(simulation_max) <- "max"

    terra::values(simulation_min) <- min_run[[q]]
    names(simulation_min) <- "min"

    terra::values(simulation_median) <- single_run[[q]]
    names(simulation_median) <- "median"

    terra::values(simulation_probability) <- probability[[q]]
    names(simulation_probability) <- "probability"

    terra::values(simulation_susceptible) <- susceptible_run[[q]]
    names(simulation_susceptible) <- "susceptible"

    if (q == 1) {
      simulation_mean_stack <- simulation_mean
      simulation_sd_stack <- simulation_sd
      simulation_min_stack <- simulation_min
      simulation_max_stack <- simulation_max
      simulation_median_stack <- simulation_median
      simulation_probability_stack <- simulation_probability
    } else {
      simulation_mean_stack <- c(simulation_mean_stack, simulation_mean)
      simulation_sd_stack <- c(simulation_sd_stack, simulation_sd)
      simulation_min_stack <- c(simulation_min_stack, simulation_min)
      simulation_max_stack <- c(simulation_max_stack, simulation_max)
      simulation_median_stack <- c(simulation_median_stack, simulation_median)
      simulation_probability_stack <-
        c(simulation_probability_stack, simulation_probability)
    }
  }

  if (!is.null(config$mask)) {
    simulation_probability_stack[is.na(config$mask)] <- 0
    simulation_mean_stack[is.na(config$mask)] <- 0
    simulation_sd_stack[is.na(config$mask)] <- 0
    simulation_min_stack[is.na(config$mask)] <- 0
    simulation_max_stack[is.na(config$mask)] <- 0
    # single_run[is.na(config$mask)] <- 0
    # susceptible_run[is.na(config$mask)] <- 0
    # exposed_run[is.na(config$mask)] <- 0
  }

  outputs <-
    list(simulation_probability_stack,
         simulation_mean_stack,
         simulation_sd_stack,
         simulation_min_stack,
         simulation_max_stack,
         single_run,
         susceptible_run,
         exposed_run,
         number_infecteds,
         infected_areas,
         west_rate,
         east_rate,
         south_rate,
         north_rate,
         escape_probability,
         north_distance_to_quarantine,
         south_distance_to_quarantine,
         east_distance_to_quarantine,
         west_distance_to_quarantine)

  names(outputs) <-
    c("probability",
      "simulation_mean",
      "simulation_sd",
      "simulation_min",
      "simulation_max",
      "single_run",
      "susceptible_run",
      "exposed_run",
      "number_infecteds",
      "infected_areas",
      "west_rate",
      "east_rate",
      "south_rate",
      "north_rate",
      "escape_probability",
      "north_distance_to_quarantine",
      "south_distance_to_quarantine",
      "east_distance_to_quarantine",
      "west_distance_to_quarantine")

  return(outputs)
}
