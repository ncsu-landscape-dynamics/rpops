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
#' @inheritParams pops_simulate
#'
#' @importFrom terra global rast xres yres classify extract ext as.points ncol project
#' nrow nlyr rowFromCell colFromCell values as.matrix rowFromCell colFromCell
#' crs app vect
#' @importFrom stats runif rnorm cov
#' @importFrom lubridate interval time_length mdy %within%
#' @importFrom MASS mvrnorm
#' @importFrom Metrics rmse
#' @importFrom utils write.csv read.table read.csv
#' @importFrom methods is
#'
#' @return a dataframe of the variables saved and their success metrics for each run
#'
#' @export

calibrate <- function(config_rds_file) {
  config <- readRDS(config_rds_file)

  if (!is.null(config$failure)) {
    stop(config$failure)
  }

  if (config$success_metric %notin% success_metric_options) {
    stop(success_metric_error)
  }

  config <- set_success_metrics(config)

  # set the parameter function to only need the parameters that change so that
  # each call to param func needs to pass in the parameters being calibrated
  param_func <-
    function(reproductive_rate,
             natural_distance_scale,
             anthropogenic_distance_scale,
             percent_natural_dispersal,
             natural_kappa,
             anthropogenic_kappa) {

      config$reproductive_rate <- reproductive_rate
      config$natural_distance_scale <- natural_distance_scale
      config$anthropogenic_distance_scale <- anthropogenic_distance_scale
      config$percent_natural_dispersal <- percent_natural_dispersal
      config$natural_kappa <- natural_kappa
      config$anthropogenic_kapp <- anthropogenic_kappa
      config$random_seed <- as.integer(sample.int(1e9, 1, replace = FALSE))
      set.seed(config$random_seed)
      config <- draw_parameters(config) # draws parameter set for the run
      config <- host_pool_setup(config)
      while (any(config$total_hosts > config$total_populations, na.rm = TRUE) ||
             any(config$total_exposed > config$total_populations, na.rm = TRUE) ||
             any(config$total_infecteds > config$total_populations, na.rm = TRUE)) {
        config <- host_pool_setup(config)
      }
      config$competency_table_list <- competency_table_list_creator(config$competency_table)
      config$pest_host_table_list <- pest_host_table_list_creator(config$pest_host_table)
      config$random_seeds <- as.matrix(config$random_seeds_list[1, ])[1, ]

      data <- pops_model(config)
      return(data)
    }

  # Check which calibration method is being used either Approximate Bayesian
  # Computation or Markov Chain Monte Carlo.
  if (config$calibration_method == "ABC") {
    # set up data structures for storing results
    config$parameters_kept <- matrix(ncol = 16, nrow = config$num_particles)
    config$parameters_test <- matrix(ncol = 16, nrow = 200)
    config$acceptance_rate <- 1
    config$acceptance_rates <- matrix(ncol = 1, nrow = config$number_of_generations)

    config$quantity_thresholds <- matrix(ncol = 1, nrow = config$number_of_generations)
    config$allocation_threshold <- matrix(ncol = 1, nrow = config$number_of_generations)
    config$configuration_threshold <- matrix(ncol = 1, nrow = config$number_of_generations)
    config$accuracy_thresholds <- matrix(ncol = 1, nrow = config$number_of_generations)
    config$precision_thresholds <- matrix(ncol = 1, nrow = config$number_of_generations)
    config$recall_thresholds <- matrix(ncol = 1, nrow = config$number_of_generations)
    config$specificity_thresholds <- matrix(ncol = 1, nrow = config$number_of_generations)
    config$rmse_thresholds <- matrix(ncol = 1, nrow = config$number_of_generations)
    config$distance_thresholds <- matrix(ncol = 1, nrow = config$number_of_generations)
    config$mcc_thresholds <- matrix(ncol = 1, nrow = config$number_of_generations)

    # assign thresholds for summary static values to be compared to the
    config$quantity_threshold <- 40 # starting threshold for quantity disagreement
    config$allocation_threshold <- 40 # starting threshold for allocation disagreement
    config$configuration_threshold <- 0.20 # starting threshold for configuration disagreement
    config$accuracy_threshold <- 0.70 # starting threshold for model accuracy
    config$precision_threshold <- 0.70 # starting threshold for model precision
    config$recall_threshold <- 0.70 # starting threshold for model recall
    config$specificity_threshold <- 0.70 # starting threshold for model
    config$rmse_threshold <- 5 # starting threshold for RMSE (root mean squared error)
    config$distance_threshold <- 1000 # starting threshold for distance between simulated
    # and observed data in units
    config$mcc_threshold <- 0.50 # starting threshold for Mathews Correlation Coefficient
    config$acceptance_rate_particle_check <- seq(60, 200, 20)
    config$save_particle_check <- seq(100, 900, 100)
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
          if (config$res$ew_res > 1000 || config$res$ns_res > 1000) {
            proposed_natural_distance_scale <- round(runif(1, 0.5, 500), digits = 1) * 10
            if (config$params_to_estimate[4]) {
              proposed_anthropogenic_distance_scale <- round(runif(1, 30, 800), digits = 0) * 100
            } else {
              proposed_anthropogenic_distance_scale <- 0.1
            }
          } else {
            proposed_natural_distance_scale <- round(runif(1, 0.5, 500), digits = 1)
            if (config$params_to_estimate[4]) {
              proposed_anthropogenic_distance_scale <- round(runif(1, 30, 80), digits = 0) * 100
            } else {
              proposed_anthropogenic_distance_scale <- 0.1
            }
          }
          if (config$params_to_estimate[3]) {
            proposed_percent_natural_dispersal <- round(runif(1, 0.87, 1), digits = 3)
          } else {
            proposed_percent_natural_dispersal <- 1.0
          }
          if (config$params_to_estimate[5]) {
            proposed_natural_kappa <- round(runif(1, 0, 5), digits = 1)
          } else {
            proposed_natural_kappa <- config$natural_kappa
          }
          if (config$params_to_estimate[6]) {
            proposed_anthropogenic_kappa <- round(runif(1, 0, 5), digits = 1)
          } else {
            proposed_anthropogenic_kappa <- config$anthropogenic_kappa
          }
        } else {
          # draw from the multivariate normal distribution and ensure that
          # parameters are within their allowed range
          proposed_parameters <-
            MASS::mvrnorm(1, config$parameter_means, config$parameter_cov_matrix)
          while (proposed_parameters[1] < 0.1 ||
                 proposed_parameters[2] < 0.1 ||
                 proposed_parameters[3] > 1.00 ||
                 proposed_parameters[3] <= 0.85 ||
                 proposed_parameters[4] < 0.1 ||
                 proposed_parameters[5] < 0 ||
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

        data <-
          param_func(
            proposed_reproductive_rate,
            proposed_natural_distance_scale,
            proposed_anthropogenic_distance_scale,
            proposed_percent_natural_dispersal,
            proposed_natural_kappa,
            proposed_anthropogenic_kappa
          )

        # calculate comparison metrics for simulation data for each time step in the simulation
        all_disagreement <- calculate_all_stats(config, data)
        all_disagreement <- colSums(all_disagreement)

        all_disagreement <- as.data.frame(t(all_disagreement))
        all_disagreement <- all_disagreement / length(data$host_pools[[1]]$infected)
        config$quantity <- all_disagreement$quantity_disagreement
        config$allocation <- all_disagreement$allocation_disagreement
        config$configuration_dis <- all_disagreement$configuration_disagreement
        config$accuracy <- all_disagreement$accuracy
        config$precision <- all_disagreement$precision
        config$recall <- all_disagreement$recall
        config$specificity <- all_disagreement$specificity
        config$rmse <- all_disagreement$rmse
        config$distance_difference <- all_disagreement$distance_difference
        config$mcc <- all_disagreement$mcc

        # Check that statistics are improvements
        model_improved <- TRUE
        if (config$use_quantity && model_improved) {
          if (config$quantity <= config$quantity_threshold) {
            model_improved <- TRUE
          } else {
            model_improved <- FALSE
          }
        }

        if (config$use_allocation && model_improved) {
          if (config$allocation <= config$allocation_threshold) {
            model_improved <- TRUE
          } else {
            model_improved <- FALSE
          }
        }

        if (config$use_configuration && model_improved) {
          if (config$configuration_dis <= config$configuration_threshold) {
            model_improved <- TRUE
          } else {
            model_improved <- FALSE
          }
        }

        if (config$use_accuracy && model_improved) {
          if (config$accuracy >= config$accuracy_threshold) {
            model_improved <- TRUE
          } else {
            model_improved <- FALSE
          }
        }

        if (config$use_precision && model_improved) {
          if (config$precision >= config$precision_threshold) {
            model_improved <- TRUE
          } else {
            model_improved <- FALSE
          }
        }

        if (config$use_recall && model_improved) {
          if (config$recall >= config$recall_threshold) {
            model_improved <- TRUE
          } else {
            model_improved <- FALSE
          }
        }

        if (config$use_specificity && model_improved) {
          if (config$specificity >= config$specificity_threshold) {
            model_improved <- TRUE
          } else {
            model_improved <- FALSE
          }
        }

        if (config$use_mcc && model_improved) {
          if (config$mcc >= config$mcc_threshold) {
            model_improved <- TRUE
          } else {
            model_improved <- FALSE
          }
        }

        if (config$use_distance && model_improved) {
          if (config$distance_difference <= config$distance_threshold) {
            model_improved <- TRUE
          } else {
            model_improved <- FALSE
          }
        }

        if (config$use_rmse && model_improved) {
          if (config$rmse <= config$rmse_threshold) {
            model_improved <- TRUE
          } else {
            model_improved <- FALSE
          }
        }

        if (model_improved && config$total_particles <= config$num_particles) {
          config$parameters_kept[config$total_particles, ] <-
            c(
              proposed_reproductive_rate,
              proposed_natural_distance_scale,
              proposed_percent_natural_dispersal,
              proposed_anthropogenic_distance_scale,
              proposed_natural_kappa,
              proposed_anthropogenic_kappa,
              config$accuracy,
              config$precision,
              config$recall,
              config$specificity,
              config$rmse,
              config$distance_difference,
              config$mcc,
              config$quantity,
              config$allocation,
              config$configuration_dis
            )

          if (config$current_bin == 1 && config$proposed_particles <= 200) {
            config$parameters_test[config$proposed_particles, ] <-
              c(
                proposed_reproductive_rate,
                proposed_natural_distance_scale,
                proposed_percent_natural_dispersal,
                proposed_anthropogenic_distance_scale,
                proposed_natural_kappa,
                proposed_anthropogenic_kappa,
                config$accuracy,
                config$precision,
                config$recall,
                config$specificity,
                config$rmse,
                config$distance_difference,
                config$mcc,
                config$quantity,
                config$allocation,
                config$configuration_dis
              )
          }
          config$current_particles <- config$current_particles + 1
          config$total_particles <- config$total_particles + 1
          config$proposed_particles <- config$proposed_particles + 1
        } else {
          if (config$current_bin == 1 && config$proposed_particles <= 200) {
            config$parameters_test[config$proposed_particles, ] <-
              c(
                proposed_reproductive_rate,
                proposed_natural_distance_scale,
                proposed_percent_natural_dispersal,
                proposed_anthropogenic_distance_scale,
                proposed_natural_kappa,
                proposed_anthropogenic_kappa,
                config$accuracy,
                config$precision,
                config$recall,
                config$specificity,
                config$rmse,
                config$distance_difference,
                config$mcc,
                config$quantity,
                config$allocation,
                config$configuration_dis
              )
          }
          config$proposed_particles <- config$proposed_particles + 1
        }

        config$acceptance_rate <- config$current_particles / config$proposed_particles
        config <- create_cal_print(config)

        if (config$current_particle %in% config$save_particle_check) {
          saveRDS(config, file.path(config$output_path, "calibration_not_complete.rds"))
        }
        ## Check that acceptance rates are within a range for the first generation
        ## if the acceptance rate is less than 5% or greater than 15% adjust the
        ## thresholds to bring the acceptance rate within that range.
        if (config$proposed_particles %in% config$acceptance_rate_particle_check &&
            config$current_bin == 1
            ) {
          if (config$acceptance_rate < 0.05) {
            config$accuracy_threshold <-
              mean(c(median(config$parameters_test[, 7], na.rm = TRUE), config$accuracy_threshold)) - 0.03
            config$precision_threshold <-
              mean(c(median(config$parameters_test[, 8], na.rm = TRUE), config$precision_threshold))
              - 0.03
            config$recall_threshold <-
              mean(c(median(config$parameters_test[, 9], na.rm = TRUE), config$recall_threshold)) - 0.03
            config$specificity_threshold <-
              mean(c(median(config$parameters_test[, 10], na.rm = TRUE), config$specificity_threshold))
              - 0.03
            config$rmse_threshold <-
              mean(c(median(config$parameters_test[, 11], na.rm = TRUE), config$rmse_threshold)) + 2
            config$distance_threshold <-
              mean(c(median(config$parameters_test[, 12], na.rm = TRUE), config$distance_threshold)) + 10
            config$mcc_threshold <-
              mean(c(median(config$parameters_test[, 13], na.rm = TRUE), config$mcc_threshold)) - 0.02
            config$quantity_threshold <-
              mean(c(median(config$parameters_test[, 14], na.rm = TRUE), config$quantity_threshold)) + 0.02
            config$allocation_threshold <-
              mean(c(median(config$parameters_test[, 15], na.rm = TRUE),
                     config$allocation_threshold)) + 0.02
            config$configuration_threshold <-
              mean(c(median(config$parameters_test[, 16], na.rm = TRUE),
                     config$configuration_threshold)) + 0.02
            ## reset starting point of parameters kept and acceptance rate
            config$parameters_kept <- matrix(ncol = 16, nrow = config$num_particles)
            config$parameters_test <- matrix(ncol = 16, nrow = 200)
            config$current_particles <- 1
            config$total_particles <- 1
            config$proposed_particles <- 1
          } else if (config$acceptance_rate > 0.15) {
            config$accuracy_threshold <- median(config$parameters_kept[, 7], na.rm = TRUE)
            config$precision_threshold <- median(config$parameters_kept[, 8], na.rm = TRUE)
            config$recall_threshold <- median(config$parameters_kept[, 9], na.rm = TRUE)
            config$specificity_threshold <- median(config$parameters_kept[, 10], na.rm = TRUE)
            config$rmse_threshold <- median(config$parameters_kept[, 11], na.rm = TRUE)
            config$distance_threshold <- median(config$parameters_kept[, 12], na.rm = TRUE)
            config$mcc_threshold <- median(config$parameters_kept[, 13], na.rm = TRUE)
            config$quantity_threshold <- median(config$parameters_kept[, 14], na.rm = TRUE)
            config$allocation_threshold <- median(config$parameters_kept[, 15], na.rm = TRUE)
            config$configuration_threshold <- median(config$parameters_kept[, 16], na.rm = TRUE)
            ## reset starting point of parameters kept and acceptance rate
            config$parameters_kept <- matrix(ncol = 16, nrow = config$num_particles)
            config$parameters_test <- matrix(ncol = 16, nrow = 200)
            config$current_particles <- 1
            config$total_particles <- 1
            config$proposed_particles <- 1
          }
        }

        if (config$verbose) {
          cat(config$acceptance_rate_info)
        }
      }

      start_index <- config$current_bin * config$generation_size - config$generation_size + 1
      end_index <- config$current_bin * config$generation_size
      config$parameter_means <- colMeans(config$parameters_kept[start_index:end_index, 1:6])
      config$parameter_cov_matrix <- cov(config$parameters_kept[start_index:end_index, 1:6])

      config$current_particles <- 1
      config$proposed_particles <- 1
      config$quantity_thresholds <- config$quantity_threshold
      config$allocation_threshold <- config$allocation_threshold
      config$configuration_threshold <- config$configuration_dis_threshold
      config$acceptance_rates[config$current_bin] <- config$acceptance_rate
      config$accuracy_thresholds[config$current_bin] <- config$accuracy_threshold
      config$precision_thresholds[config$current_bin] <- config$precision_threshold
      config$recall_thresholds[config$current_bin] <- config$recall_threshold
      config$rmse_thresholds[config$current_bin] <- config$rmse_threshold
      config$distance_thresholds[config$current_bin] <- config$distance_threshold
      config$specificity_thresholds[config$current_bin] <- config$specificity_threshold
      config$mcc_thresholds[config$current_bin] <- config$mcc_threshold
      config$accuracy_threshold <- median(config$parameters_kept[start_index:end_index, 7])
      config$precision_threshold <- median(config$parameters_kept[start_index:end_index, 8])
      config$recall_threshold <- median(config$parameters_kept[start_index:end_index, 9])
      config$specificity_threshold <- median(config$parameters_kept[start_index:end_index, 10])
      config$rmse_threshold <- median(config$parameters_kept[start_index:end_index, 11])
      config$distance_threshold <- median(config$parameters_kept[start_index:end_index, 12])
      config$mcc_threshold <- median(config$parameters_kept[start_index:end_index, 13])
      config$quantity_threshold <- median(config$parameters_kept[start_index:end_index, 14])
      config$allocation_threshold <- median(config$parameters_kept[start_index:end_index, 15])
      config$configuration_threshold <- median(config$parameters_kept[start_index:end_index, 16])
      config$current_bin <- config$current_bin + 1

      if (config$current_bin < config$number_of_generations) {
        saveRDS(config, file.path(config$output_path, "calibration_not_complete.rds"))
      }
    }

    calibrated_means <- colMeans(config$parameters_kept[start_index:end_index, 1:6])
    calibrated_cov_matrix <- cov(config$parameters_kept[start_index:end_index, 1:6])

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
      proposed_anthropogenic_distance_scale <- round(runif(1, 30, 100), digits = 0) * 100
    } else {
      proposed_anthropogenic_distance_scale <- 0.1
    }
    if (config$params_to_estimate[5]) {
      proposed_natural_kappa <- round(runif(1, 0, 5), digits = 1)
    } else {
      proposed_natural_kappa <- config$natural_kappa
    }
    if (config$params_to_estimate[6]) {
      proposed_anthropogenic_kappa <- round(runif(1, 0, 5), digits = 1)
    } else {
      proposed_anthropogenic_kappa <- config$anthropogenic_kappa
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

    all_disagreement <- calculate_all_stats(config, data)
    all_disagreement <- colSums(all_disagreement)

    all_disagreement <- as.data.frame(t(all_disagreement))
    all_disagreement <- all_disagreement / length(data$host_pools[[1]]$infected)
    config$accuracy <- all_disagreement$accuracy
    config$precision <- all_disagreement$precision
    config$recall <- all_disagreement$recall
    config$specificity <- all_disagreement$specificity
    config$rmse <- all_disagreement$rmse
    config$distance_difference <- all_disagreement$distance_difference
    config$mcc <- all_disagreement$mcc

    ## save current state of the system
    current <-
      data.frame(all_disagreement[, c("quantity_disagreement", "allocation_disagreement",
                                      "configuration_disagreement", "accuracy", "precision",
                                      "recall", "specificity", "rmse", "distance_difference",
                                      "false_negatives", "false_positives", "true_positives",
                                      "true_negatives", "odds_ratio", "mcc")],
                 reproductive_rate = proposed_reproductive_rate,
                 natural_distance_scale = proposed_natural_distance_scale,
                 anthropogenic_distance_scale = proposed_anthropogenic_distance_scale,
                 percent_natural_dispersal = proposed_percent_natural_dispersal,
                 natural_kappa = proposed_natural_kappa,
                 anthropogenic_kappa = proposed_anthropogenic_kappa
      )

    params <-
      data.frame(quantity = rep(0, config$number_of_iterations),
                 allocation = rep(0, config$number_of_iterations),
                 configuration = rep(0, config$number_of_iterations),
                 accuracy = rep(0, config$number_of_iterations),
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
                        sd = current$percent_natural_dispersal / 20), digits = 3)
        }
      } else {
        proposed_percent_natural_dispersal <- 1.0
      }

      if (config$params_to_estimate[4]) {
        proposed_anthropogenic_distance_scale <- 0
        while (proposed_anthropogenic_distance_scale <= 1 ||
               proposed_anthropogenic_distance_scale > 100000) {
          proposed_anthropogenic_distance_scale <-
            round(rnorm(1, mean = current$anthropogenic_distance_scale,
                        sd = current$anthropogenic_distance_scale / 20), digits = 0)
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
        proposed_natural_kappa <- config$natural_kappa
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
        proposed_anthropogenic_kappa <- config$anthropogenic_kappa
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
      all_disagreement <- calculate_all_stats(config, data)
      all_disagreement <- colSums(all_disagreement)

      all_disagreement <- as.data.frame(t(all_disagreement))
      all_disagreement <- all_disagreement / length(data$host_pools[[1]]$infected)
      proposed <-
        data.frame(all_disagreement[, c("quantity_disagreement", "allocation_disagreement",
                                        "configuration_disagreement", "accuracy", "precision",
                                        "recall", "specificity", "rmse", "distance_difference",
                                        "false_negatives", "false_positives", "true_positives",
                                        "true_negatives", "odds_ratio", "mcc")],
                   reproductive_rate = proposed_reproductive_rate,
                   natural_distance_scale = proposed_natural_distance_scale,
                   anthropogenic_distance_scale = proposed_anthropogenic_distance_scale,
                   percent_natural_dispersal = proposed_percent_natural_dispersal,
                   natural_kappa = proposed_natural_kappa,
                   anthropogenic_kappa = proposed_anthropogenic_kappa
        )


      # make sure no proposed statistics are 0 or the calculation fails
      # instead set them all to the lowest possible non-zero value
      if (proposed$quantity_disagreement == 0) {
        proposed$quantity_disagreement <- 0.001
      }
      if (proposed$allocation_disagreement == 0) {
        proposed$allocation_disagreement <- 0.001
      }
      if (proposed$configuration_disagreement == 0) {
        proposed$configuration_disagreement <- 0.001
      }
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
      quantity_test <- min(1, current$quantity_disagreement / proposed$quantity_disagreement)
      allocation_test <- min(1, current$allocation_disagreement / proposed$allocation_disagreement)
      configuration_test <-
        min(1, current$configuration_disagreement / proposed$configuration_disagreement)
      rmse_test <- min(1, current$rmse / proposed$rmse)
      distance_test <- min(1, current$distance / proposed$distance)

      accurracy_test <- min(1, proposed$accuracy / current$accuracy)
      precision_test <- min(1, proposed$precision / current$precision)
      recall_test <- min(1, proposed$recall / current$recall)
      specificity_test <- min(1, proposed$specificity / current$specificity)
      mcc_test <- min(1, proposed$mcc / current$mcc)

      quantity_pass <- runif(1) <= quantity_test
      allocation_pass <- runif(1) <= allocation_test
      configuration_pass <- runif(1) <= configuration_test
      accurracy_pass <- runif(1) <= accurracy_test
      precision_pass <- runif(1) <= precision_test
      recall_pass <- runif(1) <= recall_test
      specificity_pass <- runif(1) <= specificity_test
      rmse_pass <- runif(1) <= rmse_test
      distance_pass <- runif(1) <= distance_test
      mcc_pass <- runif(1) <= mcc_test

      proposed_accepted <- TRUE
      if (config$use_quantity && proposed_accepted) {
        if (quantity_pass) {
          proposed_accepted <- TRUE
        } else {
          proposed_accepted <- FALSE
        }
      }

      if (config$use_allocation && proposed_accepted) {
        if (allocation_pass) {
          proposed_accepted <- TRUE
        } else {
          proposed_accepted <- FALSE
        }
      }

      if (config$use_configuration && proposed_accepted) {
        if (configuration_pass) {
          proposed_accepted <- TRUE
        } else {
          proposed_accepted <- FALSE
        }
      }

      if (config$use_accuracy && proposed_accepted) {
        if (accurracy_pass) {
          proposed_accepted <- TRUE
        } else {
          proposed_accepted <- FALSE
        }
      }

      if (config$use_precision && proposed_accepted) {
        if (precision_pass) {
          proposed_accepted <- TRUE
        } else {
          proposed_accepted <- FALSE
        }
      }

      if (config$use_recall && proposed_accepted) {
        if (recall_pass) {
          proposed_accepted <- TRUE
        } else {
          proposed_accepted <- FALSE
        }
      }

      if (config$use_specificity && proposed_accepted) {
        if (specificity_pass) {
          proposed_accepted <- TRUE
        } else {
          proposed_accepted <- FALSE
        }
      }

      if (config$use_mcc && proposed_accepted) {
        if (mcc_pass) {
          proposed_accepted <- TRUE
        } else {
          proposed_accepted <- FALSE
        }
      }

      if (config$use_distance && proposed_accepted) {
        if (distance_pass) {
          proposed_accepted <- TRUE
        } else {
          proposed_accepted <- FALSE
        }
      }

      if (config$use_rmse && proposed_accepted) {
        if (rmse_pass) {
          proposed_accepted <- TRUE
        } else {
          proposed_accepted <- FALSE
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
                        "anthropogenic_kappa")])

    calibrated_cov_matrix <-
      cov(params[start_index:config$number_of_iterations,
                 c("reproductive_rate",
                   "natural_distance_scale",
                   "percent_natural_dispersal",
                   "anthropogenic_distance_scale",
                   "natural_kappa",
                   "anthropogenic_kappa")])

    config$parameters_kept <- params

  } else {
    return("Calibration method must be one of 'ABC' or 'MCMC'")
  }

  if (config$prior_number_of_observations < 1) {
    config$prior_weight <- config$prior_number_of_observations
    config$total_number_of_observations <- config$number_of_observations +
      round(config$number_of_observations * config$prior_number_of_observations)
    config$weight <- 1 - config$prior_weight
  } else if (config$prior_number_of_observations >= 1) {
    config$total_number_of_observations <-
      config$prior_number_of_observations + config$number_of_observations
    config$prior_weight <- config$prior_number_of_observations / config$total_number_of_observations
    config$weight <- 1 - config$prior_weight
  }

  # Use prior and calibrated parameters to update to posteriors
  posterior_check <-
    bayesian_mnn_checks(
      config$prior_means,
      config$prior_cov_matrix,
      calibrated_means,
      calibrated_cov_matrix,
      config$prior_weight,
      config$weight
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
      config$total_number_of_observations, config$parameters_kept
    )
  names(outputs) <-
    c(
      "posterior_means", "posterior_cov_matrix",
      "total_number_of_observations", "raw_calibration_data"
    )

  config$parameters_kept <- as.data.frame(config$parameters_kept)
  names(config$parameters_kept) <-  c("reproductive_rate",
                               "natural_distance_scale",
                               "percent_natural_dispersal",
                               "anthropogenic_distance_scale",
                               "natural_kappa",
                               "anthropogenic_kappa",
                               "accuracy",
                               "precision",
                               "recall",
                               "specificity",
                               "rmse",
                               "distance_difference",
                               "mcc",
                               "quantity_disagreement",
                               "allocation_disagreement",
                               "configuration_disagreement")

  file_name <- paste(config$output_path, "calibration_outputs.rdata", sep = "")
  save(outputs, file = file_name)
  file_name <- paste(config$output_path, "posterior_means.csv", sep = "")
  write.csv(posterior_means, file_name, row.names = FALSE)
  file_name <- paste(config$output_path, "posterior_cov_matrix.csv", sep = "")
  write.csv(posterior_cov_matrix, file_name, row.names = FALSE)
  file_name <- paste(config$output_path, "raw_calibration_data.csv", sep = "")
  write.csv(config$parameters_kept, file_name, row.names = FALSE)

  return(outputs)
}
