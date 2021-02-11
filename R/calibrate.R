#' Calibrates the reproductive rate and dispersal scales of the pops model.
#'
#' Either Approximate Bayesian Computation or Markov Chain Monte Carlo
#' Approximation is used to estimate relevant model parameters. Model accuracy
#' is gauged using a custom quantity allocation disagreement function to assess
#' accuracy of spatial configuration. We test number of predictions, number of
#' predicted locations, cumulative distance to nearest infection. The
#' calibration uses these metrics to determine if a run is kept if it is under
#' a threshold. either because it improves the results or randomly gets kept
#' despite being worse. We recommend running calibration for at least 10,000
#' iterations but even more will provide a better result. If the model converges
#' and doesn't improve for awhile it will exist calibration prior to reaching
#' the total number of iterations specified.
#'
#' @inheritParams pops
#' @param infected_years_file Raster file with years of initial
#' infection/infestation as individual locations of a pest or pathogen
#' @param number_of_observations the number of observations used for this
#' calibartion
#' @param number_of_generations the number of generations to use to decrease
#' the uncertainty in the parameter estimation (too many and it will take a
#' long time, too few and your parameter sets will be too wide)
#' @param generation_size how many accepted parameter sets should occur in each
#' generation
#' @param prior_number_of_observations the number of total observations from
#' previous calibrations used to weight the posterior distributions (if this is
#' a new calibration this value takes the form of a prior weight (0 - 1))
#' @param params_to_estimate A list of booleans specificing which parameters to
#' estimate ordered from (reproductive_rate, natural_dispersal_distance,
#' percent_natural_dispersal, anthropogenic_dispersal_distance, natural kappa,
#' and anthropogenic kappa)
#' @param success_metric Choose which success metric to use for calibration.
#' Choices are "number of locations", "number of locations and total distance",
#' "number of locations, number of infections, and total distance", or
#' "residual error" for 'ABC" method. Choices are "quantity", "quantity and
#' configuration", "residual error" and "odds ratio" for the 'MCMC' method.
#' Default is "number of locations and total distance"
#' @param prior_means A vector of the means of your parameters you are
#' estimating in order from (reproductive_rate, natural_dispersal_distance,
#' percent_natural_dispersal, anthropogenic_dispersal_distance, natural kappa,
#' and anthropogenic kappa)
#' @param prior_cov_matrix A covariance matrix from the previous years
#' posterior parameter estimation ordered from (reproductive_rate,
#' natural_dispersal_distance, percent_natural_dispersal,
#' anthropogenic_dispersal_distance, natural kappa, and anthropogenic kappa)
#' @param mask Raster file used to provide a mask to remove 0's that are not
#' true negatives from comparisons (e.g. mask out lakes and oceans from statics
#' if modeling terrestrial species).
#' @param checks A list of the 4 starting check values in order of # of
#' locations, total min distance, residual error, and # infected. default is
#' (500,500000, 100000, 1000). Starting check values can play a role in speed
#' of calibration and in success of calibration.
#' @param natural_kappa sets the strength of the natural direction in the
#' von-mises distribution numeric value between 0.01 and 12
#' @param anthropogenic_kappa sets the strength of the anthropogenic direction
#' in the von-mises distribution numeric value between 0.01 and 12
#' @param calibration_method choose which method of calibration to use either
#' 'ABC' (Approximate Bayesian Computation) or 'MCMC' (Markov Chain Monte Carlo
#' Approximation)
#' @param number_of_iterations how many iterations do you want to run to allow
#' the calibration to converge (recommend a minimum of at least 100,000 but
#' preferably 1 million).
#'
#' @importFrom terra global rast xres yres classify extract ext as.points ncol
#' nrow nlyr rowFromCell colFromCell values as.matrix rowFromCell colFromCell
#' crs app
#' @importFrom stats runif rnorm cov
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar% %do% %:% foreach
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom lubridate interval time_length mdy %within%
#' @importFrom MASS mvrnorm
#'
#' @return a dataframe of the variables saved and their success metrics for
#' each run
#'
#' @export

calibrate <- function(infected_years_file,
                      number_of_observations,
                      prior_number_of_observations,
                      prior_means,
                      prior_cov_matrix,
                      params_to_estimate = c(T, T, T, T, F, F),
                      number_of_generations = 7,
                      generation_size = 1000,
                      checks = c(500, 500000, 100000, 1000),
                      infected_file,
                      host_file,
                      total_populations_file,
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
                      natural_kappa = 0,
                      anthropogenic_dir = "NONE",
                      anthropogenic_kappa = 0,
                      pesticide_duration = c(0),
                      pesticide_efficacy = 1.0,
                      mask = NULL,
                      success_metric =
                        "number of locations and total distance",
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
                      calibration_method = "ABC",
                      number_of_iterations = 100000) {

  # add all data to config list
  config <- c()
  config$infected_years_file <- infected_years_file
  config$number_of_observations <- number_of_observations
  config$prior_number_of_observations <- prior_number_of_observations
  config$prior_means <- prior_means
  config$prior_cov_matrix <- prior_cov_matrix
  config$params_to_estimate <- params_to_estimate
  config$number_of_generations <- number_of_generations
  config$generation_size <- generation_size
  config$checks <- checks
  config$infected_file <- infected_file
  config$host_file <- host_file
  config$total_populations_file <- total_populations_file
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
  config$natural_kappa <- natural_kappa
  config$anthropogenic_dir <- anthropogenic_dir
  config$anthropogenic_kappa <- anthropogenic_kappa
  config$pesticide_duration <- pesticide_duration
  config$pesticide_efficacy <- pesticide_efficacy
  config$mask <- mask
  config$success_metric <- success_metric
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
  config$calibration_method <- calibration_method
  config$number_of_iterations <- number_of_iterations
  # add function name for use in configuration function to skip
  # function specific specifc configurations namely for validation and
  # calibration.
  config$function_name <- "calibrate"
  config$failure <- NULL

  # call configuration function to perform data checks and transform data into
  # format used in pops c++
  config <- configuration(config)

  if (!is.null(config$failure)) {
    return(config$failure)
  }

  # set the parameter function to only need the parameters that chanage so that
  # each call to param func needs to pass in the parameters being calibrated
  param_func <-
    function(reproductive_rate,
             natural_distance_scale,
             anthropogenic_distance_scale,
             percent_natural_dispersal,
             natural_kappa,
             anthropogenic_kappa) {
      config$random_seed <- round(runif(1, 1, 1000000))
      data <- pops_model(
        random_seed = config$random_seed,
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
        reproductive_rate = reproductive_rate,
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
          percent_natural_dispersal,
        natural_distance_scale = natural_distance_scale,
        anthropogenic_distance_scale =
          anthropogenic_distance_scale,
        natural_dir = config$natural_dir,
        natural_kappa = natural_kappa,
        anthropogenic_dir = config$anthropogenic_dir,
        anthropogenic_kappa = anthropogenic_kappa,
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
        dispersal_percentage = config$dispersal_percentage
      )
      return(data)
    }

  # Check which calibration method is being used either Approximate Bayesian
  # Computation or Markov Chain Monte Carlo.
  if (config$calibration_method == "ABC") {
    # set up data structures for storing results
    parameters_kept <- matrix(ncol = 10, nrow = config$num_particles)
    acceptance_rate <- 1
    acceptance_rates <- matrix(ncol = 1, nrow = config$number_of_generations)
    infected_checks <- matrix(ncol = 1, nrow = config$number_of_generations)
    location_checks <- matrix(ncol = 1, nrow = config$number_of_generations)
    distance_checks <- matrix(ncol = 1, nrow = config$number_of_generations)
    residual_error_checks <-
      matrix(ncol = 1, nrow = config$number_of_generations)

    # calculate comparison metrics for input data (still need to add in
    # configuration metrics to this)
    num_infected_data <- c()
    num_infected_data <- length(config$number_of_outputs)
    num_locs_data <- c()
    num_locs_data <- length(config$number_of_outputs)
    infected_data_points <-
      vector(mode = "list", length = config$number_of_outputs)

    for (y in seq_len(terra::nlyr(config$infection_years))) {
      inf_year <- config$infection_years[[y]]
      num_infected_data[[y]] <- sum(inf_year[inf_year > 0])
      num_locs_data[[y]] <- sum(inf_year[inf_year > 0] > 0)
      if (success_metric %in%
          c(
            "number of locations and total distance",
            "number of locations, number of infections, and total distance"
          )) {
        infected_data_points[[y]] <-
          terra::as.points(inf_year, crs = terra::crs(inf_year),
                           fun = function(x) { x > 0},
                           spatial = TRUE)
        }
    }

    # calculate total infections per output in the landscape
    total_infections <- terra::global(config$infection_years, "sum")
    if (length(total_infections) > config$number_of_outputs) {
      total_infections <- total_infections[1:config$number_of_outputs]
    }
    # assign thresholds for summary static values to be compared to the
    # difference between the observed and simulated data
    location_check <- config$checks[1] # number of locations
    distance_check <- config$checks[2] # minimum spatial distance
    residual_error_check <- config$checks[3] # residual error used when
    # observations are very accurate with very few areas not sampled
    infected_check <- config$checks[4] # number of pests found (either infected
    # trees or pests)

    # create raster data structures for storing simulated data for comparison
    infected_sims <- config$infection_years
    infected_sim <- config$infection_years[[1]]

    # loop through until all generations are complete
    while (config$current_bin <= config$number_of_generations) {
      # loop until all # of parameter sets kept equals the generation size
      while (config$current_particles <= config$generation_size) {

        # draw a set of proposed parameters if current generation is 1 draw from
        # a uniform distribution otherwise draw from a multivarite normal
        # distribution with mean and covariance matrix based on the previous
        # generation values
        if (config$current_bin == 1) {
          proposed_reproductive_rate <- round(runif(1, 0.055, 8), digits = 2)
          proposed_natural_distance_scale <-
            round(runif(1, 0.5, 100), digits = 1)
          if (params_to_estimate[3]) {
            proposed_percent_natural_dispersal <-
              round(runif(1, 0.93, 1), digits = 3)
          } else {
            proposed_percent_natural_dispersal <- 1.0
          }
          if (params_to_estimate[4]) {
            proposed_anthropogenic_distance_scale <-
              round(runif(1, 30, 80), digits = 0) * 100
          } else {
            proposed_anthropogenic_distance_scale <- 0.1
          }
          if (params_to_estimate[5]) {
            proposed_natural_kappa <- round(runif(1, 0, 5), digits = 1)
          } else {
            proposed_natural_kappa <- natural_kappa
          }
          if (params_to_estimate[6]) {
            proposed_anthropogenic_kappa <- round(runif(1, 0, 5), digits = 1)
          } else {
            proposed_anthropogenic_kappa <- anthropogenic_kappa
          }
        } else {
          # draw from the multivariate normal distribution and ensure that
          # parameters are within their allowed range
          proposed_parameters <-
            mvrnorm(1, config$parameter_means, config$parameter_cov_matrix)
          while (proposed_parameters[1] <= 0 |
                 proposed_parameters[2] <= 0 |
                 proposed_parameters[3] > 1 |
                 proposed_parameters[3] < 0 |
                 proposed_parameters[4] <= 0 |
                 proposed_parameters[5] < 0 |
                 proposed_parameters[6] < 0) {
            proposed_parameters <-
              mvrnorm(1, config$parameter_means, config$parameter_cov_matrix)
          }
          proposed_reproductive_rate <- proposed_parameters[1]
          proposed_natural_distance_scale <- proposed_parameters[2]
          proposed_percent_natural_dispersal <- proposed_parameters[3]
          proposed_anthropogenic_distance_scale <- proposed_parameters[4]
          proposed_natural_kappa <- proposed_parameters[5]
          proposed_anthropogenic_kappa <- proposed_parameters[6]
        }

        # run the model with the proposed parameter set
        data <-
          param_func(
            proposed_reproductive_rate,
            proposed_natural_distance_scale,
            proposed_anthropogenic_distance_scale,
            proposed_percent_natural_dispersal,
            proposed_natural_kappa,
            proposed_anthropogenic_kappa
          )

        # create holding data structures for running comparison
        num_infected_simulated <- c()
        num_infected_simulated <- length(config$number_of_outputs)
        num_locs_simulated <- c()
        num_locs_simulated <- length(config$number_of_outputs)
        infected_sim_points <-
          vector(mode = "list", length = config$number_of_outputs)
        dist <- vector(mode = "list", length = config$number_of_outputs)
        distance_differences <-
          vector(mode = "list", length = config$number_of_outputs)
        residual_differences <- c()
        residual_differences <- length(config$number_of_outputs)

        # calculate comparison metrics for simulation data for each time step in
        # the simulation
        ## To DO add in configuration metrics to this)
        for (y in seq_len(terra::nlyr(config$infection_years))) {
          if (terra::nlyr(config$infection_years) > 1) {
            terra::values(infected_sims[[y]]) <- data$infected[[y]]
            terra::values(infected_sim) <- data$infected[[y]]
          } else {
            terra::values(infected_sim) <- data$infected[[y]]
          }

          if (!is.null(config$mask)) {
            infected_sim[is.na(config$mask)] <- 0
          }

          # calculate residual error for each time step
          diff_raster <- config$infection_years[[y]] - infected_sim
          residual_differences[[y]] <- sum(diff_raster[diff_raster != 0])

          # calculate number of infection in the simulation
          num_infected_simulated[[y]] <- sum(infected_sim[infected_sim > 0])

          num_locs_simulated[[y]] <- sum(infected_sim[infected_sim > 0] > 0)
          if (success_metric %in%
              c(
                "number of locations and total distance",
                "number of locations, number of infections, and total distance"
              )
          ) {
            infected_sim_points[[y]] <-
              terra::as.points(infected_sim, crs = crs(infected_sim),
                             fun = function(x) {
                               x > 0
                             },
                             spatial = TRUE
              )
            dist[[y]] <-
              terra::distance(infected_sim_points[[y]],
                            infected_data_points[[y]])
            if (class(dist) == "matrix") {
              distance_differences[[y]] <- apply(dist[[y]], 2, min)
            } else {
              distance_differences[[y]] <- dist[[y]]
            }
          }
        }

        if (success_metric %in%
            c(
              "number of locations and total distance",
              "number of locations, number of infections, and total distance"
            )
        ) {
          all_distances <- function(distance_differences) {
            distance_differences <-
              round(sqrt(sum(distance_differences^2)), digits = 0)
            return(distance_differences)
          }
          distance_differences <- lapply(distance_differences, all_distances)
          distance_differences <-
            unlist(distance_differences, recursive = TRUE, use.names = TRUE)
        } else {
          distance_differences <- 0
        }

        number_infected_differences <-
          sqrt((num_infected_data - num_infected_simulated)^2)
        location_differences <- sqrt((num_locs_data - num_locs_simulated)^2)

        number_infected_difference <- sum(number_infected_differences)
        location_difference <- sum(location_differences)
        residual_difference <- sum(residual_differences)
        distance_difference <- sum(distance_differences)

        # Check
        diff_checks <- FALSE
        if (success_metric == "number of locations and total distance") {
          if (location_difference <= location_check &&
              distance_difference <= distance_check) {
            diff_checks <- TRUE
          }
        } else if (success_metric == "number of locations") {
          if (location_difference <= location_check) {
            diff_checks <- TRUE
          }
        } else if (success_metric == "residual error") {
          if (residual_difference <= residual_error_check) {
            diff_checks <- TRUE
          }
        } else if (success_metric ==
                "number of locations, number of infections, and total distance"
        ) {
          if (location_difference <= location_check &&
              distance_difference <= distance_check &&
              number_infected_difference <= infected_check) {
            diff_checks <- TRUE
          }
        } else {
          return("success metric must be one of 'number of locations and total
               distance', 'number of locations', and 'residual error'")
        }

        if (diff_checks) {
          parameters_kept[config$total_particles, ] <-
            c(
              proposed_reproductive_rate,
              proposed_natural_distance_scale,
              proposed_percent_natural_dispersal,
              proposed_anthropogenic_distance_scale,
              proposed_natural_kappa,
              proposed_anthropogenic_kappa,
              number_infected_difference,
              location_difference,
              distance_difference,
              residual_difference
            )
          config$current_particles <- config$current_particles + 1
          config$total_particles <- config$total_particles + 1
          config$proposed_particles <- config$proposed_particles + 1
        } else {
          config$proposed_particles <- config$proposed_particles + 1
        }
        acceptance_rate <- config$current_particles / config$proposed_particles
        acceptance_rate_info <-
          paste("The current generation is ", config$current_bin, " and the
            current particle is ", config$current_particles,
                " and the current acceptance rate is ", acceptance_rate,
                sep = "")
        print(acceptance_rate_info)
      }

      start_index <- config$current_bin * generation_size - generation_size + 1
      end_index <- config$current_bin * generation_size
      config$parameter_means <-
        colMeans(parameters_kept[start_index:end_index, 1:6])
      config$parameter_cov_matrix <-
        cov(parameters_kept[start_index:end_index, 1:6])

      config$current_particles <- 1
      config$proposed_particles <- 1
      acceptance_rates[config$current_bin] <- acceptance_rate
      infected_checks[config$current_bin] <- infected_check
      location_checks[config$current_bin] <- location_check
      distance_checks[config$current_bin] <- distance_check
      residual_error_checks[config$current_bin] <- residual_error_check
      infected_check <- median(parameters_kept[start_index:end_index, 7])
      location_check <- median(parameters_kept[start_index:end_index, 8])
      distance_check <- median(parameters_kept[start_index:end_index, 9])
      residual_error_check <- median(parameters_kept[start_index:end_index, 10])
      config$current_bin <- config$current_bin + 1
    }

    calibrated_means <- colMeans(parameters_kept[start_index:end_index, 1:6])
    calibrated_cov_matrix <- cov(parameters_kept[start_index:end_index, 1:6])

  } else if (config$calibration_method == "MCMC") {

    proposed_reproductive_rate <- round(runif(1, 0.05, 8), digits = 2)
    proposed_natural_distance_scale <- round(runif(1, 0.5, 100), digits = 1)
    if (config$params_to_estimate[3]) {
      proposed_percent_natural_dispersal <-
        round(runif(1, 0.93, 1), digits = 3)
    } else {
      proposed_percent_natural_dispersal <- 1.0
    }
    if (config$params_to_estimate[4]) {
      proposed_anthropogenic_distance_scale <-
        round(runif(1, 30, 80), digits = 0) * 100
    } else {
      proposed_anthropogenic_distance_scale <- 0.1
    }
    if (config$params_to_estimate[5]) {
      proposed_natural_kappa <- round(runif(1, 0, 5), digits = 1)
    } else {
      proposed_natural_kappa <- natural_kappa
    }
    if (config$params_to_estimate[6]) {
      proposed_anthropogenic_kappa <- round(runif(1, 0, 5), digits = 1)
    } else {
      proposed_anthropogenic_kappa <- anthropogenic_kappa
    }

    data <-
      param_func(
        proposed_reproductive_rate,
        proposed_natural_distance_scale,
        proposed_anthropogenic_distance_scale,
        proposed_percent_natural_dispersal,
        proposed_natural_kappa,
        proposed_anthropogenic_kappa
      )

    comp_year <- terra::rast(config$infected_file)
    all_disagreement <-
      foreach::foreach(
        q = seq_len(length(data$infected)),
        .combine = rbind,
        .packages = c("terra", "PoPS"),
        .final = colSums
      ) %do% {
        terra::values(comp_year) <- data$infected[[q]]
        quantity_allocation_disagreement(
          config$infection_years[[q]],
          comp_year,
          config$configuration,
          config$mask
        )
      }

    ## save current state of the system
    current <- best <-
      data.frame(t(all_disagreement),
                 reproductive_rate = proposed_reproductive_rate,
                 natural_distance_scale = proposed_natural_distance_scale,
                 anthropogenic_distance_scale =
                   proposed_anthropogenic_distance_scale,
                 percent_natural_dispersal = proposed_percent_natural_dispersal,
                 natural_kappa = proposed_natural_kappa,
                 anthropogenic_kappa = proposed_anthropogenic_kappa
      )

    params <-
      data.frame(quantity_disagreement = rep(0, config$number_of_iterations),
                 allocation_disagreement = rep(0, config$number_of_iterations),
                 total_disagreement = rep(0, config$number_of_iterations),
                 configuration_disagreement =
                   rep(0, config$number_of_iterations),
                 omission = rep(0, config$number_of_iterations),
                 commission = rep(0, config$number_of_iterations),
                 true_positives = rep(0, config$number_of_iterations),
                 true_negatives = rep(0, config$number_of_iterations),
                 odds_ratio = rep(0, config$number_of_iterations),
                 residual_error = rep(0, config$number_of_iterations),
                 true_infected = rep(0, config$number_of_iterations),
                 simulated_infected = rep(0, config$number_of_iterations),
                 infected_difference = rep(0, config$number_of_iterations),
                 reproductive_rate = rep(0, config$number_of_iterations),
                 natural_distance_scale = rep(0, config$number_of_iterations),
                 anthropogenic_distance_scale =
                   rep(0, config$number_of_iterations),
                 percent_natural_dispersa = rep(0, config$number_of_iterations),
                 natural_kappa = rep(0, config$number_of_iterations),
                 anthropogenic_kappa = rep(0, config$number_of_iterations))

    for (i in seq_len(config$number_of_iterations)) {

      proposed_reproductive_rate <- 0
      while (proposed_reproductive_rate <= 0.1) {
        proposed_reproductive_rate <-
          round(rnorm(1, mean = current$reproductive_rate,
                      sd = current$reproductive_rate / 10), digits = 1)
      }

      proposed_natural_distance_scale <- 0
      while (proposed_natural_distance_scale <= 0.1) {
        proposed_natural_distance_scale <-
          round(rnorm(1, mean = current$natural_distance_scale,
                      sd = current$natural_distance_scale / 10), digits = 0)
      }

      if (config$params_to_estimate[3]) {
        proposed_percent_natural_dispersal <- 0
        while (proposed_percent_natural_dispersal < 0.93 ||
               proposed_percent_natural_dispersal >= 1) {
          proposed_percent_natural_dispersal <-
            round(rnorm(1, mean = current$percent_natural_dispersal,
                        sd = current$percent_natural_dispersal / 20),
                  digits = 3)
        }
      } else {
        proposed_percent_natural_dispersal <- 1.0
      }

      if (config$params_to_estimate[4]) {
        proposed_anthropogenic_distance_scale <- 0
        while (proposed_anthropogenic_distance_scale <= 0.01 |
               proposed_anthropogenic_distance_scale > 100000) {
          proposed_anthropogenic_distance_scale <-
            round(rnorm(1, mean = current$anthropogenic_distance_scale,
                        sd = current$anthropogenic_distance_scale / 20),
                  digits = 0)
        }
      } else {
        proposed_anthropogenic_distance_scale <- 0.1
      }

      if (config$params_to_estimate[5]) {
        proposed_natural_kappa <- 0
        while (proposed_natural_kappa <= 0 ||
               proposed_natural_kappa > 4.000) {
          proposed_natural_kappa <-
            round(rnorm(1, mean = current$natural_kappa,
                        sd = current$natural_kappa / 20), digits = 3)
        }
      } else {
        proposed_natural_kappa <- natural_kappa
      }

      if (config$params_to_estimate[6]) {
        proposed_anthropogenic_kappa <- 0
        while (proposed_anthropogenic_kappa <= 0 ||
               proposed_anthropogenic_kappa > 4.000) {
          proposed_anthropogenic_kappa <-
            round(rnorm(1, mean = current$anthropogenic_kappa,
                        sd = current$anthropogenic_kappa / 20), digits = 3)
        }
      } else {
        proposed_anthropogenic_kappa <- anthropogenic_kappa
      }

      data <-
        param_func(
          proposed_reproductive_rate,
          proposed_natural_distance_scale,
          proposed_anthropogenic_distance_scale,
          proposed_percent_natural_dispersal,
          proposed_natural_kappa,
          proposed_anthropogenic_kappa
        )

      # set up comparison
      comp_year <- terra::rast(config$infected_file)
      all_disagreement <-
        foreach::foreach(
          q = seq_len(length(data$infected)),
          .combine = rbind,
          .packages = c("terra", "PoPS"),
          .final = colSums
        ) %do% {
          terra::values(comp_year) <- data$infected[[q]]
          quantity_allocation_disagreement(
            config$infection_years[[q]],
            comp_year,
            config$configuration,
            config$mask
          )
        }

      proposed <-
        data.frame(t(all_disagreement),
                   reproductive_rate = proposed_reproductive_rate,
                   natural_distance_scale = proposed_natural_distance_scale,
                   anthropogenic_distance_scale =
                     proposed_anthropogenic_distance_scale,
                   percent_natural_dispersal =
                     proposed_percent_natural_dispersal,
                   natural_kappa = proposed_natural_kappa,
                   anthropogenic_kappa = proposed_anthropogenic_kappa
        )

      # make sure no proposed statistics are 0 or the calculation fails
      # instead set them all to the lowest possible non-zero value
      if (proposed$allocation_disagreement == 0) {
        proposed$allocation_disagreement <- 1
      }
      if (proposed$quantity_disagreement == 0) {
        proposed$quantity_disagreement <- 1
      }
      if (proposed$total_disagreement == 0) {
        proposed$total_disagreement <- 1
      }
      if (proposed$configuration_disagreement == 0) {
        proposed$configuration_disagreement <- 0.01
      }
      if (proposed$residual_error == 0) {
        proposed$residual_error <- 1
      }
      # Set up tests for to see if new variable is an improvement in
      # performance metrics
      quantity_test <-
        min(1, current$quantity_disagreement / proposed$quantity_disagreement)
      configuration_test <-
        min(1, current$configuration_disagreement /
              proposed$configuration_disagreement)
      # odds ratio is treated differently than all the other metrics as it is
      # the only one where higher numbers means better model performance
      oddsratio_test <-
        min(1, proposed$odds_ratio / current$odds_ratio)
      residual_error_test <-
        min(1, current$residual_error / proposed$residual_error)

      quantity_pass <- runif(1) <= quantity_test
      configuration_pass <- runif(1) <= configuration_test
      oddsratio_pass <- runif(1) <= oddsratio_test
      residual_error_pass <- runif(1) <= residual_error_test

      proposed_accepted <- FALSE
      if (success_metric == "quantity" & quantity_pass) {
        proposed_accepted <- TRUE
      } else if (success_metric == "quantity and configuration" &
                 quantity_pass &
                 configuration_pass) {
        proposed_accepted <- TRUE
      } else if (success_metric == "odds ratio" & oddsratio_pass) {
        proposed_accepted <- TRUE
      } else if (success_metric == "residual error" & residual_error_pass) {
        proposed_accepted <- TRUE
      } else {
        return("Success metric must be one of 'quantity', 'quantity and
                 configuration', 'residual error', or 'odds ratio'")
      }

      if (proposed_accepted) {
        current <- proposed
        if (current$quantity_disagreement <= best$quantity_disagreement) {
          best <- current
        }
      }

      param <- current

      print(i)
      params[i, ] <- param
    }

    if (config$number_of_iterations > 10000) {
      start_index <- 5000
    } else {
      start_index <- number_of_iterations / 2
    }

    calibrated_means <-
      colMeans(params[start_index:config$number_of_iterations, 14:19])
    calibrated_cov_matrix <-
      cov(params[start_index:config$number_of_iterations, 14:19])

    parameters_kept <- params

  } else {
    return("Calibration method must be one of 'ABC' or 'MCMC'")
  }


  if (prior_number_of_observations < 1) {
    prior_weight <- prior_number_of_observations
    total_number_of_observations <-
      number_of_observations +
      round(number_of_observations * prior_number_of_observations)
    weight <- 1 - prior_weight
  } else if (prior_number_of_observations >= 1) {
    total_number_of_observations <-
      prior_number_of_observations + number_of_observations
    prior_weight <- prior_number_of_observations / total_number_of_observations
    weight <- 1 - prior_weight
  }

  # Use prior and calibrated parameters to update to posteriors
  posterior_check <-
    bayesian_mnn_checks(
      prior_means,
      prior_cov_matrix,
      calibrated_means,
      calibrated_cov_matrix,
      prior_weight,
      weight
    )

  if (posterior_check$checks_passed) {
    posterior_means <- posterior_check$posterior_means
    posterior_cov_matrix <- posterior_check$posterior_cov_matrix
  } else {
    return(posterior_check$failed_check)
  }

  outputs <-
    list(
      posterior_means, posterior_cov_matrix,
      total_number_of_observations, parameters_kept
    )
  names(outputs) <-
    c(
      "posterior_means", "posterior_cov_matrix",
      "total_number_of_observations", "raw_calibration_data"
    )
  return(outputs)
}
