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
#' infection/infestation as individual locations of a pest or pathogen. This is
#' a multiband raster file (e.g. .tif) with each band representing a unique time
#' step (e.g. band 1 = year 1 .... band 6 = year 6 or band 1 = week 1 .... band
#' 6 = week 6). This needs to allign with both the time step selection and start
#' and end dates selection. Units for infections are based on data availability
#' and the way the units used for your host file creation (e.g. percent area, #
#' of hosts per cell, etc.). This doesn't include the start year which passed in
#' in the initial_infected_file (e.g. if we had observation data from 2017,
#' 2018, and 2019 the 2017 raster file would be the initial_infected_file and a
#' dual band raster file would have band 1 = 2018 and band 2 = 2019 observations)
#' @param number_of_observations the number of observations used for this
#' calibration. Useful if using previous calibration. This is used to
#' weight the parameters when updating parameters when new data becomes
#' available. Example if we have 2,000 observations in 2019 and had 1,000
#' observations in 2018 and 1,000 in 2017, we would use 2,000 here and 2,000 for
#' our prior_number_of_observations.
#' @param number_of_generations the number of generations to use to decrease
#' the uncertainty in the parameter estimation (too many and it will take a
#' long time, too few and your parameter sets will be too wide). This is an ABC
#' implementation naming convention but should be set to greater than 7 for
#' robust calibrations. There is a tradeoff between computational time and model
#' accuracy the larger this number gets. Usually 7 to 9 is the ideal range.
#' @param generation_size how many accepted parameter sets should occur in each
#' generation. For example if generation size is 1,000 then the simulation runs
#' until 1,000 model runs are less than the threshold value.
#' We recommend running at least 1,000 but the greater this number the more
#' accurate the model parameters selected will be.
#' @param prior_number_of_observations the number of total observations from
#' previous calibrations used to weight the posterior distributions (if this is
#' a new calibration this value takes the form of a prior weight (0 - 1)). This
#' is used to weight the parameters when updating parameters when new data
#' becomes available. Example if we have 2,000 observations in 2019 and had
#' 1,000 observations in 2018 and 1,000 in 2017, we would use 2,000 here and
#' 2,000 for our number_of_observations.
#' @param params_to_estimate A list of booleans specificing which parameters to
#' estimate ordered from (reproductive_rate, natural_dispersal_distance,
#' percent_natural_dispersal, anthropogenic_dispersal_distance, natural kappa,
#' and anthropogenic kappa)
#' @param prior_means A vector of the means of your parameters you are
#' estimating in order from (reproductive_rate, natural_dispersal_distance,
#' percent_natural_dispersal, anthropogenic_dispersal_distance, natural kappa,
#' and anthropogenic kappa). This is used when updating a parameter set from a
#' previous calibration using the iterative framework.
#' @param prior_cov_matrix A covariance matrix from the previous years
#' posterior parameter estimation ordered from (reproductive_rate,
#' natural_dispersal_distance, percent_natural_dispersal,
#' anthropogenic_dispersal_distance, natural kappa, and anthropogenic kappa).
#' This is used when updating a parameter set from a previous calibration using
#' the iterative framework.
#' @param mask Raster file used to provide a mask to remove 0's that are not
#' true negatives from comparisons (e.g. mask out lakes and oceans from statics
#' if modeling terrestrial species). A numerical value represents the area you
#' want to calculate statistics on and an NA value represents the area to remove
#' from the statistics.
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
#' @param verbose Boolean with true printing current status of calibration,
#' (e.g. the current generation, current particle, and the acceptance rate).
#' Defaults if FALSE.
#' @param write_outputs Either c("summary_outputs", or "None"). If not
#' "None" output folder path must be provided.
#' @param output_folder_path this is the full path with either / or \\ (e.g.,
#' "C:/user_name/desktop/pops_sod_2020_2023/outputs/")
#' @param use_distance Boolean if you want to compare distance between
#' simulations and observations. Default is FALSE.
#' @param use_rmse Boolean if you want to calibrate based on rmse. This is
#' useful if you have very good population level observations. Default is FALSE.
#'
#' @importFrom terra global rast xres yres classify extract ext as.points ncol
#' nrow nlyr rowFromCell colFromCell values as.matrix rowFromCell colFromCell
#' crs app vect
#' @importFrom stats runif rnorm cov
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar% %do% %:% foreach
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom lubridate interval time_length mdy %within%
#' @importFrom MASS mvrnorm
#' @importFrom Metrics rmse
#' @importFrom utils write.csv
#'
#' @return a dataframe of the variables saved and their success metrics for
#' each run
#'
#' @export

calibrate <- function(infected_years_file,
                      number_of_observations = 1,
                      prior_number_of_observations = 0,
                      prior_means = c(0, 0, 0, 0, 0, 0),
                      prior_cov_matrix = matrix(0, 6, 6),
                      params_to_estimate = c(T, T, T, T, F, F),
                      number_of_generations = 7,
                      generation_size = 1000,
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
                      mortality_frequency = "year",
                      mortality_frequency_n = 1,
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
                      calibration_method = "ABC",
                      number_of_iterations = 100000,
                      exposed_file = "",
                      verbose = TRUE,
                      write_outputs = "None",
                      output_folder_path = "",
                      use_distance = FALSE,
                      use_rmse = FALSE) {

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
  config$calibration_method <- calibration_method
  config$number_of_iterations <- number_of_iterations
  config$exposed_file <- exposed_file
  # add function name for use in configuration function to skip
  # function specific specific configurations namely for validation and
  # calibration.
  config$function_name <- "calibrate"
  config$failure <- NULL
  config$write_outputs <- write_outputs
  config$output_folder_path <- output_folder_path
  config$mortality_frequency <- mortality_frequency
  config$mortality_frequency_n <- mortality_frequency_n
  config$use_distance <- use_distance
  config$use_rmse <- use_rmse

  # call configuration function to perform data checks and transform data into
  # format used in pops c++
  config <- configuration(config)

  if (!is.null(config$failure)) {
    return(config$failure)
  }

  # set the parameter function to only need the parameters that change so that
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
        lethal_temperature_month = config$lethal_temperature_month,
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
        res = config$res,
        rows_cols = config$rows_cols,
        time_step = config$time_step,
        reproductive_rate = reproductive_rate,
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
        percent_natural_dispersal = percent_natural_dispersal,
        natural_distance_scale = natural_distance_scale,
        anthropogenic_distance_scale = anthropogenic_distance_scale,
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
        mortality_frequency = config$mortality_frequency,
        mortality_frequency_n = config$mortality_frequency_n,
        use_spreadrates = config$use_spreadrates,
        model_type_ = config$model_type,
        latency_period = config$latency_period,
        generate_stochasticity = config$generate_stochasticity,
        establishment_stochasticity = config$establishment_stochasticity,
        movement_stochasticity = config$movement_stochasticity,
        deterministic = config$deterministic,
        establishment_probability = config$establishment_probability,
        dispersal_percentage = config$dispersal_percentage,
        use_overpopulation_movements = config$use_overpopulation_movements,
        overpopulation_percentage = config$overpopulation_percentage,
        leaving_percentage = config$leaving_percentage,
        leaving_scale_coefficient = config$leaving_scale_coefficient
      )
      return(data)
    }

  # Check which calibration method is being used either Approximate Bayesian
  # Computation or Markov Chain Monte Carlo.
  if (config$calibration_method == "ABC") {
    # set up data structures for storing results
    parameters_kept <- matrix(ncol = 12, nrow = config$num_particles)
    parameters_test <- matrix(ncol = 12, nrow = 200)
    acceptance_rate <- 1
    acceptance_rates <- matrix(ncol = 1, nrow = config$number_of_generations)

    accuracy_thresholds <-
      matrix(ncol = 1, nrow = config$number_of_generations)
    precision_thresholds <-
      matrix(ncol = 1, nrow = config$number_of_generations)
    recall_thresholds <-
      matrix(ncol = 1, nrow = config$number_of_generations)
    specificity_thresholds <-
      matrix(ncol = 1, nrow = config$number_of_generations)
    rmse_thresholds <-
      matrix(ncol = 1, nrow = config$number_of_generations)
    distance_thresholds <-
      matrix(ncol = 1, nrow = config$number_of_generations)

    # assign thresholds for summary static values to be compared to the
    accuracy_threshold <- 0.70 # starting threshold for model accuracy
    precision_threshold <- 0.70 # starting threshold for model precision
    recall_threshold <- 0.70 # starting threshold for model recall
    specificity_threshold <- 0.70 # starting threshold for model
    rmse_threshold <- 7 # starting threshold for RMSE (root mean squared error)
    distance_threshold <- 1000 # starting threshold for distance between simulated
    # and observed data in units
    acceptance_rate_particle_check <- seq(60, 200, 20)

    # loop through until all generations are complete
    while (config$current_bin <= config$number_of_generations) {
      # loop until all # of parameter sets kept equals the generation size
      while (config$current_particles <= config$generation_size) {
        # draw a set of proposed parameters if current generation is 1 draw from
        # a uniform distribution otherwise draw from a multivariate normal
        # distribution with mean and co-variance matrix based on the previous
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
            MASS::mvrnorm(1, config$parameter_means, config$parameter_cov_matrix)
          while (proposed_parameters[1] < 0.1 |
                 proposed_parameters[2] < 0.1 |
                 proposed_parameters[3] > 1.00 |
                 proposed_parameters[3] <= 0.92 |
                 proposed_parameters[4] < 0.1 |
                 proposed_parameters[5] < 0 |
                 proposed_parameters[6] < 0) {
            proposed_parameters <-
              MASS::mvrnorm(1, config$parameter_means, config$parameter_cov_matrix)
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

        # calculate comparison metrics for simulation data for each time step in
        # the simulation
        all_disagreement <-
          foreach::foreach(
            q = seq_len(length(data$infected)),
            .combine = rbind,
            .packages = c("terra", "PoPS"),
            .final = colSums
          ) %do% {
            comparison <- terra::rast(config$infected_file)
            reference <- terra::rast(config$infected_file)
            terra::values(comparison) <- data$infected[[q]]
            terra::values(reference) <- config$infection_years2[[q]]
            mask <- terra::rast(config$infected_file)
            terra::values(mask) <- config$mask_matrix
            quantity_allocation_disagreement(reference,
                                             comparison,
                                             use_configuration = FALSE,
                                             mask = mask,
                                             use_distance = config$use_distance)
          }

        all_disagreement <- as.data.frame(t(all_disagreement))
        all_disagreement <- all_disagreement / length(data$infected)
        accuracy <- all_disagreement$accuracy
        precision <- all_disagreement$precision
        recall <- all_disagreement$recall
        specificity <- all_disagreement$specificity
        rmse <- all_disagreement$rmse
        distance_difference <- all_disagreement$distance_difference

        # Check that statistics are improvements
        model_improved <- FALSE
        if (accuracy >= accuracy_threshold &&
            precision >= precision_threshold &&
            recall >= recall_threshold &&
            specificity >= specificity_threshold) {

          if (config$use_distance) {
            if (distance_difference <= distance_threshold) {
              model_improved <- TRUE
            } else {
              model_improved <- FALSE
            }
          } else {
            model_improved <- TRUE
          }

          if (config$use_rmse) {
            if (rmse <= rmse_threshold) {
              model_improved <- TRUE
            } else {
              model_improved <- FALSE
            }
          } else {
            model_improved <- TRUE
          }
        }

        if (model_improved & config$total_particles <= config$num_particles) {
          parameters_kept[config$total_particles, ] <-
            c(
              proposed_reproductive_rate,
              proposed_natural_distance_scale,
              proposed_percent_natural_dispersal,
              proposed_anthropogenic_distance_scale,
              proposed_natural_kappa,
              proposed_anthropogenic_kappa,
              accuracy,
              precision,
              recall,
              specificity,
              rmse,
              distance_difference
            )
          if (config$current_bin == 1 && config$proposed_particles <= 200) {
            parameters_test[config$proposed_particles, ] <-
              c(
                proposed_reproductive_rate,
                proposed_natural_distance_scale,
                proposed_percent_natural_dispersal,
                proposed_anthropogenic_distance_scale,
                proposed_natural_kappa,
                proposed_anthropogenic_kappa,
                accuracy,
                precision,
                recall,
                specificity,
                rmse,
                distance_difference
              )
          }
          config$current_particles <- config$current_particles + 1
          config$total_particles <- config$total_particles + 1
          config$proposed_particles <- config$proposed_particles + 1
        } else {
          if (config$current_bin == 1 && config$proposed_particles <= 200) {
            parameters_test[config$proposed_particles, ] <-
              c(
                proposed_reproductive_rate,
                proposed_natural_distance_scale,
                proposed_percent_natural_dispersal,
                proposed_anthropogenic_distance_scale,
                proposed_natural_kappa,
                proposed_anthropogenic_kappa,
                accuracy,
                precision,
                recall,
                specificity,
                rmse,
                distance_difference
              )
          }
          config$proposed_particles <- config$proposed_particles + 1
        }

        acceptance_rate <- config$current_particles / config$proposed_particles
        acceptance_rate_info <- paste(
                            "generation:            ",
                            config$current_bin,
                            "\nparticle:              ",
                            config$current_particles,
                            "\nacceptance rate:       ",
                            format(acceptance_rate, digits = 5),
                            "\naccuracy:              ",
                            accuracy,
                            "\naccuracy threshold:    ",
                            accuracy_threshold,
                            "\nprecision:             ",
                            precision,
                            "\nprecision threshold:   ",
                            precision_threshold,
                            "\nrecall:                ",
                            recall,
                            "\nrecall threshold:      ",
                            recall_threshold,
                            "\nspecificity:           ",
                            specificity,
                            "\nspecificity threshold: ",
                            specificity_threshold,
                            "\nrmse:                  ",
                            rmse,
                            "\nrmse threshold:        ",
                            rmse_threshold,
                            "\ndistance difference:   ",
                            distance_difference,
                            "\ndistance threshold:    ",
                            distance_threshold,
                            "\n\n",
                            sep = " ")

        ## Check that acceptance rates are within a range for the first generation
        ## if the acceptance rate is less than 5% or greater than 15% adjust the
        ## thresholds to bring the acceptance rate within that range.
        if (config$proposed_particles %in% acceptance_rate_particle_check &&
            config$current_bin == 1
            ) {
          if (acceptance_rate < 0.05) {
            accuracy_threshold <-
              mean(c(median(parameters_test[, 7], na.rm = TRUE), accuracy_threshold)) - 0.03
            precision_threshold <-
              mean(c(median(parameters_test[, 8], na.rm = TRUE), precision_threshold)) - 0.03
            recall_threshold <-
              mean(c(median(parameters_test[, 9], na.rm = TRUE), recall_threshold)) - 0.03
            specificity_threshold <-
              mean(c(median(parameters_test[, 10], na.rm = TRUE), specificity_threshold)) - 0.03
            rmse_threshold <-
              mean(c(median(parameters_test[, 11], na.rm = TRUE), rmse_threshold)) + 2
            distance_threshold <-
              mean(c(median(parameters_test[, 12], na.rm = TRUE), distance_threshold)) + 10
            ## reset starting point of parameters kept and acceptance rate
            parameters_kept <- matrix(ncol = 12, nrow = config$num_particles)
            parameters_test <- matrix(ncol = 12, nrow = 200)
            config$current_particles <- 1
            config$total_particles <- 1
            config$proposed_particles <- 1
          } else if (acceptance_rate > 0.15) {
            accuracy_threshold <- median(parameters_kept[, 7], na.rm = TRUE)
            precision_threshold <- median(parameters_kept[, 8], na.rm = TRUE)
            recall_threshold <- median(parameters_kept[, 9], na.rm = TRUE)
            specificity_threshold <- median(parameters_test[, 10], na.rm = TRUE)
            rmse_threshold <- median(parameters_test[, 11], na.rm = TRUE)
            distance_threshold <- median(parameters_test[, 12], na.rm = TRUE)
            ## reset starting point of parameters kept and acceptance rate
            parameters_kept <- matrix(ncol = 12, nrow = config$num_particles)
            parameters_test <- matrix(ncol = 12, nrow = 200)
            config$current_particles <- 1
            config$total_particles <- 1
            config$proposed_particles <- 1
          }
        }

        if (verbose) {
          cat(acceptance_rate_info)
        }
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
      accuracy_thresholds[config$current_bin] <- accuracy_threshold
      precision_thresholds[config$current_bin] <- precision_threshold
      recall_thresholds[config$current_bin] <- recall_threshold
      rmse_thresholds[config$current_bin] <- rmse_threshold
      distance_threshold[config$current_bin] <- distance_threshold
      accuracy_threshold <- median(parameters_kept[start_index:end_index, 7])
      precision_threshold <- median(parameters_kept[start_index:end_index, 8])
      recall_threshold <- median(parameters_kept[start_index:end_index, 9])
      specificity_threshold <- median(parameters_kept[start_index:end_index, 10])
      rmse_threshold <- median(parameters_kept[start_index:end_index, 11])
      distance_threshold <- median(parameters_kept[start_index:end_index, 12])
      config$current_bin <- config$current_bin + 1
    }

    calibrated_means <- colMeans(parameters_kept[start_index:end_index, 1:6])
    calibrated_cov_matrix <- cov(parameters_kept[start_index:end_index, 1:6])

  } else if (config$calibration_method == "MCMC") {
    proposed_reproductive_rate <- round(runif(1, 0.05, 8), digits = 2)
    proposed_natural_distance_scale <- round(runif(1, 0.5, 100), digits = 1)
    if (config$params_to_estimate[3]) {
      proposed_percent_natural_dispersal <-
        round(runif(1, 0.93, 1.000), digits = 3)
    } else {
      proposed_percent_natural_dispersal <- 1.0
    }
    if (config$params_to_estimate[4]) {
      proposed_anthropogenic_distance_scale <-
        round(runif(1, 30, 100), digits = 0) * 100
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

    all_disagreement <-
      foreach::foreach(
        q = seq_len(length(data$infected)),
        .combine = rbind,
        .packages = c("terra", "PoPS"),
        .final = colSums
      ) %do% {
        comparison <- terra::rast(config$infected_file)
        reference <- terra::rast(config$infected_file)
        terra::values(comparison) <- data$infected[[q]]
        terra::values(reference) <- config$infection_years2[[q]]
        mask <- terra::rast(config$infected_file)
        terra::values(mask) <- config$mask_matrix
        quantity_allocation_disagreement(reference,
                                         comparison,
                                         use_configuration = FALSE,
                                         mask = mask,
                                         use_distance = config$use_distance)
      }

    all_disagreement <- as.data.frame(t(all_disagreement))
    all_disagreement <- all_disagreement / length(data$infected)
    accuracy <- all_disagreement$accuracy
    precision <- all_disagreement$precision
    recall <- all_disagreement$recall
    specificity <- all_disagreement$specificity
    rmse <- all_disagreement$rmse
    distance_difference <- all_disagreement$distance_difference

    ## save current state of the system
    current <-
      data.frame(all_disagreement[, c("accuracy", "precision", "recall",
                                      "specificity", "rmse",
                                      "distance_difference", "false_negatives",
                                      "false_positives", "true_positives",
                                      "true_negatives", "odds_ratio")],
                 reproductive_rate = proposed_reproductive_rate,
                 natural_distance_scale = proposed_natural_distance_scale,
                 anthropogenic_distance_scale =
                   proposed_anthropogenic_distance_scale,
                 percent_natural_dispersal = proposed_percent_natural_dispersal,
                 natural_kappa = proposed_natural_kappa,
                 anthropogenic_kappa = proposed_anthropogenic_kappa
      )

    params <-
      data.frame(accuracy = rep(0, config$number_of_iterations),
                 precision = rep(0, config$number_of_iterations),
                 recall = rep(0, config$number_of_iterations),
                 specificity = rep(0, config$number_of_iterations),
                 rmse = rep(0, config$number_of_iterations),
                 distance_difference = rep(0, config$number_of_iterations),
                 false_negatives = rep(0, config$number_of_iterations),
                 false_positives = rep(0, config$number_of_iterations),
                 true_positives = rep(0, config$number_of_iterations),
                 true_negatives = rep(0, config$number_of_iterations),
                 odds_ratio = rep(0, config$number_of_iterations),
                 reproductive_rate = rep(0, config$number_of_iterations),
                 natural_distance_scale = rep(0, config$number_of_iterations),
                 anthropogenic_distance_scale =
                   rep(0, config$number_of_iterations),
                 percent_natural_dispersal = rep(0, config$number_of_iterations),
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
      while (proposed_natural_distance_scale <= 1) {
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
        while (proposed_anthropogenic_distance_scale <= 1 |
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
      all_disagreement <-
        foreach::foreach(
          q = seq_len(length(data$infected)),
          .combine = rbind,
          .packages = c("terra", "PoPS"),
          .final = colSums
        ) %do% {
          comparison <- terra::rast(config$infected_file)
          reference <- terra::rast(config$infected_file)
          terra::values(comparison) <- data$infected[[q]]
          terra::values(reference) <- config$infection_years2[[q]]
          mask <- terra::rast(config$infected_file)
          terra::values(mask) <- config$mask_matrix
          quantity_allocation_disagreement(reference,
                                           comparison,
                                           use_configuration = FALSE,
                                           mask = mask,
                                           use_distance = config$use_distance)
        }

      all_disagreement <- as.data.frame(t(all_disagreement))
      all_disagreement <- all_disagreement / length(data$infected)
      proposed <-
        data.frame(all_disagreement[, c("accuracy", "precision", "recall",
                                        "specificity", "rmse",
                                        "distance_difference", "false_negatives",
                                        "false_positives", "true_positives",
                                        "true_negatives", "odds_ratio")],
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
      if (proposed$accuracy == 0) {
        proposed$accuracy <- 0.001
      }
      if (proposed$precision == 0) {
        proposed$precisiont <- 0.001
      }
      if (proposed$recall == 0) {
        proposed$recall <- 0.001
      }
      if (proposed$specificity == 0) {
        proposed$specificity <- 0.001
      }
      if (proposed$rmse == 0) {
        proposed$rmse <- 0.001
      }
      if (proposed$distance_difference == 0) {
        proposed$distance_difference <- 0.001
      }
      # Set up tests for to see if new variable is an improvement in
      # performance metrics for accuracy, precision, recall, and specificity
      # higher values are better so the proposed parameter is in the numerator,
      # for rmse and distance lower values are improvements and the proposed
      # value is in the denominator.
      accurracy_test <- min(1, proposed$accuracy / current$accuracy)
      precision_test <- min(1, proposed$precision / current$precision)
      recall_test <- min(1, proposed$recall / current$recall)
      specificity_test <- min(1, proposed$specificity / current$specificity)
      rmse_test <- min(1, current$rmse / proposed$rmse)
      distance_test <- min(1, current$distance / proposed$distance)

      accurracy_pass <- runif(1) <= accurracy_test
      precision_pass <- runif(1) <= precision_test
      recall_pass <- runif(1) <= recall_test
      specificity_pass <- runif(1) <= specificity_test
      rmse_pass <- runif(1) <= rmse_test
      distance_pass <- runif(1) <= distance_test

      proposed_accepted <- FALSE
      if (accurracy_pass && precision_pass && recall_pass && specificity_pass) {

        if (config$use_distance) {
          if (distance_pass) {
            proposed_accepted <- TRUE
          } else {
            proposed_accepted <- FALSE
          }
        } else {
          proposed_accepted <- TRUE
        }

        if (config$use_rmse) {
          if (rmse_pass) {
            proposed_accepted <- TRUE
          } else {
            proposed_accepted <- FALSE
          }
        } else {
          proposed_accepted <- TRUE
        }
      }

      if (proposed_accepted) {
        current <- proposed
      }

      param <- current
      if (verbose) {
        print(i)
      }
      params[i, ] <- param
    }

    if (config$number_of_iterations > 10000) {
      start_index <- 5000
    } else {
      start_index <- number_of_iterations / 2
    }

    calibrated_means <-
      colMeans(params[start_index:config$number_of_iterations,
                      c("reproductive_rate",
                        "natural_distance_scale",
                        "anthropogenic_distance_scale",
                        "percent_natural_dispersal",
                        "natural_kappa",
                        "anthropogenic_kappa")])
    calibrated_cov_matrix <-
      cov(params[start_index:config$number_of_iterations,
                 c("reproductive_rate",
                   "natural_distance_scale",
                   "anthropogenic_distance_scale",
                   "percent_natural_dispersal",
                   "natural_kappa",
                   "anthropogenic_kappa")])

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

  if (config$write_outputs %in% config$output_write_list) {
    file_name <- paste(config$output_folder_path, "calibration_outputs.rdata", sep = "")
    save(outputs, file = file_name)
    file_name <- paste(config$output_folder_path, "posterior_means.csv", sep = "")
    write.csv(posterior_means, file_name, row.names = FALSE)
    file_name <- paste(config$output_folder_path, "posterior_cov_matrix.csv", sep = "")
    write.csv(posterior_cov_matrix, file_name, row.names = FALSE)
    file_name <- paste(config$output_folder_path, "raw_calibration_data.csv", sep = "")
    write.csv(parameters_kept, file_name, row.names = FALSE)
  }

  return(outputs)
}
