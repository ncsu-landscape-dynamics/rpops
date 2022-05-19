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

calibrate <- function(config) {
  # set the parameter function to only need the parameters that change so that
  # each call to param func needs to pass in the parameters being calibrated
  param_func <-
    function(reproductive_rate,
             natural_distance_scale,
             anthropogenic_distance_scale,
             percent_natural_dispersal,
             natural_kappa,
             anthropogenic_kappa,
             network_min_distance,
             network_max_distance) {
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
        leaving_scale_coefficient = config$leaving_scale_coefficient,
        bbox = config$bounding_box,
        network_min_distance = network_min_distance,
        network_max_distance = network_max_distance,
        network_filename = config$network_filename
      )
      return(data)
    }

  # Check which calibration method is being used either Approximate Bayesian
  # Computation or Markov Chain Monte Carlo.
  if (config$calibration_method == "ABC") {
    # set up data structures for storing results
    parameters_kept <- matrix(ncol = 15, nrow = config$num_particles)
    parameters_test <- matrix(ncol = 15, nrow = 200)
    acceptance_rate <- 1
    acceptance_rates <- matrix(ncol = 1, nrow = config$number_of_generations)

    accuracy_thresholds <- matrix(ncol = 1, nrow = config$number_of_generations)
    precision_thresholds <- matrix(ncol = 1, nrow = config$number_of_generations)
    recall_thresholds <- matrix(ncol = 1, nrow = config$number_of_generations)
    specificity_thresholds <- matrix(ncol = 1, nrow = config$number_of_generations)
    rmse_thresholds <- matrix(ncol = 1, nrow = config$number_of_generations)
    distance_thresholds <- matrix(ncol = 1, nrow = config$number_of_generations)

    # assign thresholds for summary static values to be compared to the
    accuracy_threshold <- 0.70 # starting threshold for model accuracy
    precision_threshold <- 0.70 # starting threshold for model precision
    recall_threshold <- 0.70 # starting threshold for model recall
    specificity_threshold <- 0.70 # starting threshold for model
    rmse_threshold <- 7 # starting threshold for RMSE (root mean squared error)
    distance_threshold <- 1000 # starting threshold for distance between simulated
    # and observed data in units of the CRS
    mcc_threshold <- 0.50 # starting threshold for Mathews Correlation coefficient
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
          proposed_natural_distance_scale <- round(runif(1, 0.5, 100), digits = 1)
          proposed_percent_natural_dispersal <- round(runif(1, 0.93, 1), digits = 3)
          proposed_anthropogenic_distance_scale <- round(runif(1, 30, 80), digits = 0) * 100
          if (config$natural_dir != "NONE" && config$natural_kappa <= 0) {
            proposed_natural_kappa <- round(runif(1, 0, 5), digits = 1)
          } else {
            proposed_natural_kappa <- config$natural_kappa
          }
          if (config$anthropogenic_dir != "NONE" && config$anthropogenic_kappa <= 0) {
            proposed_anthropogenic_kappa <- round(runif(1, 0, 5), digits = 1)
          } else {
            proposed_anthropogenic_kappa <- config$anthropogenic_kappa
          }
          if (config$anthropogenic_kernel_type == "network") {
            proposed_network_min_distance <-
              round(runif(1, config$res$ew_res / 2, config$res$ew_res * 10), digits = 0)
            proposed_network_max_distance <-
              round(runif(1, config$res$ew_res * 10, config$res$ew_res *
                            min(config$rows_cols$num_cols, config$rows_cols$num_rows)), digits = 0)
          } else {
            proposed_network_min_distance <- config$res$ew_res / 2
            proposed_network_max_distance <-
              min(config$rows_cols$num_cols, config$rows_cols$num_rows) * config$res$ew_res
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
                 proposed_parameters[6] < 0 |
                 proposed_parameters[7] < config$res$ew_res / 2 |
                 proposed_parameters[7] > proposed_parameters[8] |
                 proposed_parameters[8] >
                 min(config$rows_cols$num_cols, config$rows_cols$num_rows) * config$res$ew_res) {
            proposed_parameters <-
              MASS::mvrnorm(1, config$parameter_means, config$parameter_cov_matrix)
          }
          proposed_reproductive_rate <- proposed_parameters[1]
          proposed_natural_distance_scale <- proposed_parameters[2]
          proposed_percent_natural_dispersal <- proposed_parameters[3]
          proposed_anthropogenic_distance_scale <- proposed_parameters[4]
          proposed_natural_kappa <- proposed_parameters[5]
          proposed_anthropogenic_kappa <- proposed_parameters[6]
          proposed_network_min_distance <- proposed_parameters[7]
          proposed_network_max_distance <- proposed_parameters[8]
        }

        # run the model with the proposed parameter set
        data <-
          param_func(
            proposed_reproductive_rate,
            proposed_natural_distance_scale,
            proposed_anthropogenic_distance_scale,
            proposed_percent_natural_dispersal,
            proposed_natural_kappa,
            proposed_anthropogenic_kappa,
            proposed_network_min_distance,
            proposed_network_max_distance
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
        mcc <- all_disagreement$mcc

        # Check that statistics are improvements
        model_improved <- FALSE
        if (config$use_mcc) {
          if (mcc >= mcc_threshold) {
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
        } else {
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
              proposed_network_min_distance,
              proposed_network_max_distance,
              accuracy,
              precision,
              recall,
              specificity,
              rmse,
              distance_difference,
              mcc
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
                proposed_network_min_distance,
                proposed_network_max_distance,
                accuracy,
                precision,
                recall,
                specificity,
                rmse,
                distance_difference,
                mcc
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
                proposed_network_min_distance,
                proposed_network_max_distance,
                accuracy,
                precision,
                recall,
                specificity,
                rmse,
                distance_difference,
                mcc
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
                            "\nMCC:                   ",
                            mcc,
                            "\nMCC threshold:         ",
                            mcc_threshold,
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
              mean(c(median(parameters_test[, 9], na.rm = TRUE), accuracy_threshold)) - 0.03
            precision_threshold <-
              mean(c(median(parameters_test[, 10], na.rm = TRUE), precision_threshold)) - 0.03
            recall_threshold <-
              mean(c(median(parameters_test[, 11], na.rm = TRUE), recall_threshold)) - 0.03
            specificity_threshold <-
              mean(c(median(parameters_test[, 12], na.rm = TRUE), specificity_threshold)) - 0.03
            rmse_threshold <-
              mean(c(median(parameters_test[, 13], na.rm = TRUE), rmse_threshold)) + 2
            distance_threshold <-
              mean(c(median(parameters_test[, 14], na.rm = TRUE), distance_threshold)) + 10
            mcc_threshold <-
              mean(c(median(parameters_test[, 15], na.rm = TRUE), mcc_threshold)) - 0.02
            ## reset starting point of parameters kept and acceptance rate
            parameters_kept <- matrix(ncol = 15, nrow = config$num_particles)
            parameters_test <- matrix(ncol = 15, nrow = 200)
            config$current_particles <- 1
            config$total_particles <- 1
            config$proposed_particles <- 1
          } else if (acceptance_rate > 0.15) {
            accuracy_threshold <- median(parameters_kept[, 9], na.rm = TRUE)
            precision_threshold <- median(parameters_kept[, 10], na.rm = TRUE)
            recall_threshold <- median(parameters_kept[, 11], na.rm = TRUE)
            specificity_threshold <- median(parameters_kept[, 12], na.rm = TRUE)
            rmse_threshold <- median(parameters_kept[, 13], na.rm = TRUE)
            distance_threshold <- median(parameters_kept[, 14], na.rm = TRUE)
            mcc_threshold <- median(parameters_kept[, 15], na.rm = TRUE)
            ## reset starting point of parameters kept and acceptance rate
            parameters_kept <- matrix(ncol = 15, nrow = config$num_particles)
            parameters_test <- matrix(ncol = 15, nrow = 200)
            config$current_particles <- 1
            config$total_particles <- 1
            config$proposed_particles <- 1
          }
        }

        if (config$verbose) {
          cat(acceptance_rate_info)
        }
      }

      start_index <- config$current_bin * config$generation_size - config$generation_size + 1
      end_index <- config$current_bin * config$generation_size
      config$parameter_means <- colMeans(parameters_kept[start_index:end_index, 1:8])
      config$parameter_cov_matrix <- cov(parameters_kept[start_index:end_index, 1:8])

      config$current_particles <- 1
      config$proposed_particles <- 1
      acceptance_rates[config$current_bin] <- acceptance_rate
      accuracy_thresholds[config$current_bin] <- accuracy_threshold
      precision_thresholds[config$current_bin] <- precision_threshold
      recall_thresholds[config$current_bin] <- recall_threshold
      rmse_thresholds[config$current_bin] <- rmse_threshold
      distance_threshold[config$current_bin] <- distance_threshold
      accuracy_threshold <- median(parameters_kept[start_index:end_index, 9])
      precision_threshold <- median(parameters_kept[start_index:end_index, 10])
      recall_threshold <- median(parameters_kept[start_index:end_index, 11])
      specificity_threshold <- median(parameters_kept[start_index:end_index, 12])
      rmse_threshold <- median(parameters_kept[start_index:end_index, 13])
      distance_threshold <- median(parameters_kept[start_index:end_index, 14])
      mcc_threshold <- median(parameters_kept[start_index:end_index, 15])
      config$current_bin <- config$current_bin + 1
    }

    calibrated_means <- colMeans(parameters_kept[start_index:end_index, 1:8])
    calibrated_cov_matrix <- cov(parameters_kept[start_index:end_index, 1:8])

  } else if (config$calibration_method == "MCMC") {

    proposed_reproductive_rate <- round(runif(1, 0.055, 8), digits = 2)
    proposed_natural_distance_scale <- round(runif(1, 0.5, 100), digits = 1)
    proposed_percent_natural_dispersal <- round(runif(1, 0.93, 1), digits = 3)
    proposed_anthropogenic_distance_scale <- round(runif(1, 30, 80), digits = 0) * 100
    if (config$natural_dir != "NONE" && config$natural_kappa <= 0) {
      proposed_natural_kappa <- round(runif(1, 0, 5), digits = 1)
    } else {
      proposed_natural_kappa <- config$natural_kappa
    }
    if (config$natural_dir != "NONE" && config$natural_kappa <= 0) {
      proposed_anthropogenic_kappa <- round(runif(1, 0, 5), digits = 1)
    } else {
      proposed_anthropogenic_kappa <- config$anthropogenic_kappa
    }
    if (config$anthropogenic_kernel_type == "network") {
      proposed_network_min_distance <-
        round(runif(1, config$res$ew_res / 2, config$res$ew_res * 10), digits = 0)
      proposed_network_max_distance <-
        round(runif(1, config$res$ew_res * 10, config$res$ew_res *
                      min(config$rows_cols$num_cols, config$rows_cols$num_rows)), digits = 0)
    } else {
      proposed_network_min_distance <- config$res$ew_res / 2
      proposed_network_max_distance <-
        min(config$rows_cols$num_cols, config$rows_cols$num_rows) * config$res$ew_res
    }

    data <-
      param_func(
        proposed_reproductive_rate,
        proposed_natural_distance_scale,
        proposed_anthropogenic_distance_scale,
        proposed_percent_natural_dispersal,
        proposed_natural_kappa,
        proposed_anthropogenic_kappa,
        proposed_network_min_distance,
        proposed_network_max_distance
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
    mcc <- all_disagreement$mcc

    ## save current state of the system
    current <-
      data.frame(all_disagreement[, c("accuracy", "precision", "recall",
                                      "specificity", "rmse",
                                      "distance_difference", "false_negatives",
                                      "false_positives", "true_positives",
                                      "true_negatives", "odds_ratio", "mcc")],
                 reproductive_rate = proposed_reproductive_rate,
                 natural_distance_scale = proposed_natural_distance_scale,
                 anthropogenic_distance_scale = proposed_anthropogenic_distance_scale,
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
                 mcc = rep(0, config$number_of_iterations),
                 reproductive_rate = rep(0, config$number_of_iterations),
                 natural_distance_scale = rep(0, config$number_of_iterations),
                 anthropogenic_distance_scale = rep(0, config$number_of_iterations),
                 percent_natural_dispersal = rep(0, config$number_of_iterations),
                 natural_kappa = rep(0, config$number_of_iterations),
                 anthropogenic_kappa = rep(0, config$number_of_iterations),
                 network_min_distance = rep(0, config$number_of_iterations),
                 network_max_distance = rep(0, config$number_of_iterations))

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

      proposed_percent_natural_dispersal <- 0
      while (proposed_percent_natural_dispersal < 0.93 ||
             proposed_percent_natural_dispersal >= 1) {
        proposed_percent_natural_dispersal <-
          round(rnorm(1, mean = current$percent_natural_dispersal,
                      sd = current$percent_natural_dispersal / 20), digits = 3)
      }

      proposed_anthropogenic_distance_scale <- 0
      while (proposed_anthropogenic_distance_scale <= 1 |
             proposed_anthropogenic_distance_scale > 100000) {
        proposed_anthropogenic_distance_scale <-
          round(rnorm(1, mean = current$anthropogenic_distance_scale,
                      sd = current$anthropogenic_distance_scale / 20), digits = 0)
      }

      proposed_natural_kappa <- 0
      if (config$natural_dir != "NONE" && config$natural_kappa <= 0) {
        proposed_natural_kappa <-
          round(rnorm(1, mean = current$natural_kappa, sd = current$natural_kappa / 20), digits = 3)
      } else {
        proposed_natural_kappa <- config$natural_kappa
      }

      proposed_anthropogenic_kappa <- 0
      if (config$anthropogenic_dir != "NONE" && config$anthropogenic_kappa <= 0) {
        proposed_anthropogenic_kappa <-
          round(rnorm(1, mean = current$anthropogenic_kappa,
                      sd = current$anthropogenic_kappa / 20), digits = 3)
      } else {
        proposed_anthropogenic_kappa <- config$anthropogenic_kappa
      }

      if (config$anthropogenic_kernel_type == "network") {
        proposed_network_min_distance <- 0
        while (proposed_network_min_distance < config$res$ew_res / 2) {
          proposed_network_min_distance <-
            round(rnorm(1, mean = current$network_min_distance,
                        sd = current$network_min_distance / 20), digits = 0)
        }
        proposed_network_max_distance <-
          round(runif(1, config$res$ew_res * 10, config$res$ew_res *
                        min(config$rows_cols$num_cols, config$rows_cols$num_rows)), digits = 0)

        proposed_network_max_distance <- 0
        while (proposed_network_max_distance < proposed_network_min_distance ||
               proposed_network_max_distance >
               (config$res$ew_res * min(config$rows_cols$num_cols, config$rows_cols$num_rows))) {
          proposed_network_max_distance <-
            round(rnorm(1, mean = current$network_max_distance,
                        sd = current$network_max_distance / 20), digits = 0)
        }
      } else {
        proposed_network_min_distance <- config$res$ew_res / 2
        proposed_network_max_distance <-
          min(config$rows_cols$num_cols, config$rows_cols$num_rows) * config$res$ew_res
      }

      data <-
        param_func(
          proposed_reproductive_rate,
          proposed_natural_distance_scale,
          proposed_anthropogenic_distance_scale,
          proposed_percent_natural_dispersal,
          proposed_natural_kappa,
          proposed_anthropogenic_kappa,
          proposed_network_min_distance,
          proposed_network_max_distance
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
                                        "true_negatives", "odds_ratio", "mcc")],
                   reproductive_rate = proposed_reproductive_rate,
                   natural_distance_scale = proposed_natural_distance_scale,
                   anthropogenic_distance_scale = proposed_anthropogenic_distance_scale,
                   percent_natural_dispersal = proposed_percent_natural_dispersal,
                   natural_kappa = proposed_natural_kappa,
                   anthropogenic_kappa = proposed_anthropogenic_kappa,
                   network_min_distance = proposed_network_min_distance,
                   network_max_distance = proposed_network_max_distance
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
      mcc_test <- min(1, current$mcc / proposed$mcc)

      accurracy_pass <- runif(1) <= accurracy_test
      precision_pass <- runif(1) <= precision_test
      recall_pass <- runif(1) <= recall_test
      specificity_pass <- runif(1) <= specificity_test
      rmse_pass <- runif(1) <= rmse_test
      distance_pass <- runif(1) <= distance_test
      mcc_pass <- runif(1) <= mcc_test

      proposed_accepted <- FALSE
      if (config$use_mcc) {
        if (mcc_pass) {

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
      } else {
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
      }


      if (proposed_accepted) {
        current <- proposed
      }

      param <- current
      if (config$verbose) {
        print(i)
      }
      params[i, ] <- param
    }

    if (config$number_of_iterations > 10000) {
      start_index <- 5000
    } else {
      start_index <- config$number_of_iterations / 2
    }

    calibrated_means <-
      colMeans(params[start_index:config$number_of_iterations,
                      c("reproductive_rate",
                        "natural_distance_scale",
                        "percent_natural_dispersal",
                        "anthropogenic_distance_scale",
                        "natural_kappa",
                        "anthropogenic_kappa",
                        "network_min_distance",
                        "network_max_distance")])

    calibrated_cov_matrix <-
      cov(params[start_index:config$number_of_iterations,
                 c("reproductive_rate",
                   "natural_distance_scale",
                   "percent_natural_dispersal",
                   "anthropogenic_distance_scale",
                   "natural_kappa",
                   "anthropogenic_kappa",
                   "network_min_distance",
                   "network_max_distance")])

    parameters_kept <- params

  } else {
    return("Calibration method must be one of 'ABC' or 'MCMC'")
  }

  if (config$prior_number_of_observations < 1) {
    prior_weight <- config$prior_number_of_observations
    total_number_of_observations <- config$number_of_observations +
      round(config$number_of_observations * config$prior_number_of_observations)
    weight <- 1 - prior_weight
  } else if (config$prior_number_of_observations >= 1) {
    total_number_of_observations <-
      config$prior_number_of_observations + config$number_of_observations
    prior_weight <- config$prior_number_of_observations / total_number_of_observations
    weight <- 1 - prior_weight
  }

  # Use prior and calibrated parameters to update to posteriors
  posterior_check <-
    bayesian_mnn_checks(
      config$prior_means,
      config$prior_cov_matrix,
      calibrated_means,
      calibrated_cov_matrix,
      prior_weight,
      weight
    )

  if (posterior_check$checks_passed) {
    posterior_means <- as.numeric(posterior_check$posterior_means)
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
