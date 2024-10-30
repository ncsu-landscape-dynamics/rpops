# These functions are designed to reduce code complexity and the need to copy
# and past code across main functions
"%notin%" <- Negate("%in%")

set_success_metrics <- function(config) {
  config$use_quantity <- FALSE
  config$use_allocation <- FALSE
  config$use_configuration <- FALSE
  config$use_accuracy <- FALSE
  config$use_precision <- FALSE
  config$use_recall <- FALSE
  config$use_specificity <- FALSE
  config$use_rmse <- FALSE
  config$use_distance <- FALSE
  config$use_mcc <- FALSE

  if (config$success_metric %in% quantity_list) {
    config$use_quantity <- TRUE
  }

  if (config$success_metric %in% allocation_list) {
    config$use_allocation <- TRUE
  }

  if (config$success_metric %in% configuration_list) {
    config$use_configuration <- TRUE
  }

  if (config$success_metric %in% accurracy_list) {
    config$use_accuracy <- TRUE
  }

  if (config$success_metric %in% precision_list) {
    config$use_precision <- TRUE

  }

  if (config$success_metric %in% recall_list) {
    config$use_recall <- TRUE
  }

  if (config$success_metric %in% specificity_list) {
    config$use_specificity <- TRUE
  }

  if (config$success_metric %in% rmse_list) {
    config$use_rmse <- TRUE
  }

  if (config$success_metric %in% distance_list) {
    config$use_distance <- TRUE
  }

  if (config$success_metric %in% mcc_list) {
    config$use_mcc <- TRUE
  }

  return(config)
}

create_cal_print <- function(config) {
  config$acceptance_rate_info <- paste(
    "generation:                ",
    config$current_bin,
    "\nparticle:                ",
    config$current_particles,
    "\nacceptance rate:         ",
    format(config$acceptance_rate, digits = 5), sep = " ")

  if (config$use_quantity) {
    config$acceptance_rate_info <- paste(config$acceptance_rate_info,
                                         "\nquantity:                ",
                                         config$quantity,
                                         "\nquantity threshold:      ",
                                         config$allocation_threshold,
                                         sep = " ")
  }

  if (config$use_allocation) {
    config$acceptance_rate_info <- paste(config$acceptance_rate_info,
                                         "\nallocation:              ",
                                         config$allocation,
                                         "\nallocation threshold:    ",
                                         config$allocation_threshold,
                                         sep = " ")
  }

  if (config$use_configuration) {
    config$acceptance_rate_info <- paste(config$acceptance_rate_info,
                                         "\nconfiguration:           ",
                                         config$configuration_dis,
                                         "\nconfiguration threshold: ",
                                         config$configuration_threshold,
                                         sep = " ")
  }

  if (config$use_accuracy) {
    config$acceptance_rate_info <- paste(config$acceptance_rate_info,
                                         "\naccuracy:                ",
                                         config$accuracy,
                                         "\naccuracy threshold:      ",
                                         config$accuracy_threshold,
                                         sep = " ")
  }

  if (config$use_precision) {
    config$acceptance_rate_info <- paste(config$acceptance_rate_info,
                                         "\nprecision:               ",
                                         config$precision,
                                         "\nprecision threshold:     ",
                                         config$precision_threshold,
                                         sep = " ")

  }

  if (config$use_recall) {
    config$acceptance_rate_info <- paste(config$acceptance_rate_info,
                                         "\nrecall:                  ",
                                         config$recall,
                                         "\nrecall threshold:        ",
                                         config$recall_threshold,
                                         sep = " ")
  }

  if (config$use_specificity) {
    config$acceptance_rate_info <- paste(config$acceptance_rate_info,
                                         "\nspecificity:             ",
                                         config$specificity,
                                         "\nspecificity threshold:   ",
                                         config$specificity_threshold,
                                         sep = " ")
  }

  if (config$use_rmse) {
    config$acceptance_rate_info <- paste(config$acceptance_rate_info,
                                         "\nrmse:                    ",
                                         config$rmse,
                                         "\nrmse threshold:          ",
                                         config$rmse_threshold,
                                         sep = " ")
  }

  if (config$use_distance) {
    config$acceptance_rate_info <- paste(config$acceptance_rate_info,
                                         "\ndistance difference:     ",
                                         config$distance_difference,
                                         "\ndistance threshold:      ",
                                         config$distance_threshold,
                                         sep = " ")
  }

  if (config$use_mcc) {
    config$acceptance_rate_info <- paste(config$acceptance_rate_info,
                                         "\nMCC:                     ",
                                         config$mcc,
                                         "\nMCC threshold:           ",
                                         config$mcc_threshold,
                                         sep = " ")
  }

  config$acceptance_rate_info <- paste(config$acceptance_rate_info,
                                       "\n\n",
                                       sep = " ")
  return(config)
}

draw_parameters <- function(config) {
  parameters <-
    MASS::mvrnorm(1, config$parameter_means, config$parameter_cov_matrix)
  while (any(parameters[1] < 0 |
             parameters[2] <= 0 |
             parameters[3] > 1 |
             parameters[3] <= 0 |
             parameters[4] <= 0 |
             parameters[5] < 0 |
             parameters[6] < 0 |
             parameters[7] < config$res$ew_res / 2 |
             parameters[7] > parameters[8] |
             parameters[8] >
             min(config$rows_cols$num_cols, config$rows_cols$num_rows) * config$res$ew_res)) {

    config$number_of_draws <-
      nrow(parameters[parameters[1] < 0 |
                        parameters[2] <= 0 |
                        parameters[3] > 1 |
                        parameters[3] <= 0 |
                        parameters[4] <= 0 |
                        parameters[5] < 0 |
                        parameters[6] < 0 |
                        parameters[7] < config$res$ew_res / 2 |
                        parameters[7] > parameters[8] |
                        parameters[8] >
                        min(config$rows_cols$num_cols, config$rows_cols$num_rows) *
                        config$res$ew_res
      ])

    if (is.null(config$number_of_draws)) {
      config$number_of_draws <- 1
    }

    parameters[parameters[1] < 0 |
                 parameters[2] <= 0 |
                 parameters[3] > 1 |
                 parameters[3] <= 0 |
                 parameters[4] <= 0 |
                 parameters[5] < 0 |
                 parameters[6] < 0 |
                 parameters[7] < config$res$ew_res / 2 |
                 parameters[7] > parameters[8] |
                 parameters[8] >
                 (min(config$rows_cols$num_cols, config$rows_cols$num_rows) * config$res$ew_res)] <-
      MASS::mvrnorm(
        config$number_of_draws,
        config$parameter_means,
        config$parameter_cov_matrix
      )
  }
  config$reproductive_rate <- parameters[1]
  config$natural_distance_scale <- parameters[2]
  config$percent_natural_dispersal <- parameters[3]
  config$anthropogenic_distance_scale <- parameters[4]
  config$natural_kappa <- parameters[5]
  config$anthropogenic_kappa <- parameters[6]
  config$network_min_distance <- parameters[7]
  config$network_max_distance <- parameters[8]

  return(config)
}

create_random_seeds <- function(n) {
  random_seeds <-
    data.frame(disperser_generation = sample(1:999999999, n, replace = FALSE),
               natural_dispersal = sample(1:999999999, n, replace = FALSE),
               anthropogenic_dispersal = sample(1:999999999999, n, replace = FALSE),
               establishment = sample(1:999999999, n, replace = FALSE),
               weather = sample(1:999999999, n, replace = FALSE),
               movement = sample(1:999999999, n, replace = FALSE),
               overpopulation = sample(1:999999999, n, replace = FALSE),
               survival_rate = sample(1:999999999, n, replace = FALSE),
               soil = sample(1:999999999, n, replace = FALSE))

  return(random_seeds)
}

# creates a matrix from a matrix of mean values and a matrix of standard deviations. The two
# matrices must be the same size.
matrix_norm_distribution <- function(mean_matrix, sd_matrix) {
  new_matrix <- matrix(
    round(rnorm(length(mean_matrix), mean = mean_matrix[], sd = sd_matrix[])),
    nrow = nrow(mean_matrix),
    ncol = ncol(mean_matrix)
  )
  new_matrix[is.na(new_matrix) | new_matrix < 0] <- 0
  return(new_matrix)
}

# Uncertainty propagation for raster data sets, expects a spatRaster with 2
# layers (mean and standard deviation)
output_from_raster_mean_and_sd <- function(x) {
  x[[1]] <- terra::classify(x[[1]], matrix(c(-Inf, 0, 0), ncol = 3, byrow = TRUE))
  x[[2]] <- terra::classify(x[[2]], matrix(c(-Inf, 0, 0), ncol = 3, byrow = TRUE))
  fun <- function(x) {
    round(rnorm(1, mean = x[1], sd = x[2]), digits = 0)
  }
  x2 <- suppressWarnings(terra::app(x, fun))
  return(x2)
}

# Combine two standard deviation spatRasters
combined_sd <- function(v1, v2, m1, m2, n1, n2) {
  (((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 1)) +
    (((n1 * n2) * (m1 - m2)^2) / ((n1 + n2) * (n1 + n2 - 1)))
}

# Create a unique identifier for the pops_run_lite output filenames
generate_uid <- function(length = 10) {
  uid <- paste0(
    sample(c(letters[1:23], 1:9), length, replace = TRUE),
    collapse = ""
  )
  return(uid)
}

# Add time step for the pops_run_lite output filenames
format_output_date <- function(start_date, output_frequency, time_step, multiplier = 0) {
  # Ensure the start date is in Date format and check for validity
  if (is.na(as.Date(start_date, format = "%Y-%m-%d"))) {
    stop("Invalid start_date provided. Please use a valid date in 'YYYY-MM-DD' format.")
  }
  start_date <- as.Date(start_date)
  
  # Function to increment the date based on a unit and multiplier
  advance_date <- function(date, unit, multiplier) {
    switch(unit,
           "year" = date + lubridate::years(multiplier),
           "month" = date %m+% lubridate::months(multiplier),
           "week" = date + lubridate::weeks(multiplier),
           "day" = date + lubridate::days(multiplier))
  }
  
  # Determine the correct format and advance the date
  formatted_date <- if (output_frequency %in% c("year", "month")) {
    format(advance_date(start_date, output_frequency, multiplier), 
           if (output_frequency == "year") "%Y" else "%Y-%m")
  } else if (output_frequency == "week") {
    paste(format(advance_date(start_date, "week", multiplier), "%Y-%m-%d"), "w")
  } else if (output_frequency == "day") {
    paste(format(advance_date(start_date, "day", multiplier), "%Y-%m-%d"), "d")
  } else if (output_frequency == "time step" || output_frequency == "every_n_steps") {
    unit <- switch(time_step, "year" = "year", "month" = "month", "week" = "week", "day" = "day")
    advanced_date <- advance_date(start_date, unit, multiplier)
    formatted_date <- switch(unit,
                             "year" = format(advanced_date, "%Y"),
                             "month" = format(advanced_date, "%Y-%m"),
                             "week" = paste(format(advanced_date, "%Y-%m-%d"), "w"),
                             "day" = paste(format(advanced_date, "%Y-%m-%d"), "d"))
  } else {
    stop("Invalid output_frequency or time_step value provided.")
  }
  
  return(formatted_date)
}

# Helper function to update the configuration based on uncertainty settings
update_config <- function(config) {
  config <- draw_parameters(config)  # Draws parameter set for the run
  
  # Update infected and exposed based on uncertainty settings
  if (config$use_initial_condition_uncertainty) {
    config$infected <- matrix_norm_distribution(config$infected_mean, config$infected_sd)
    exposed2 <- matrix_norm_distribution(config$exposed_mean, config$exposed_sd)
    config$exposed[[config$latency_period + 1]] <- exposed2
  } else {
    config$infected <- config$infected_mean
    exposed2 <- config$exposed_mean
    config$exposed[[config$latency_period + 1]] <- exposed2
  }

  # Update host based on uncertainty settings
  if (config$use_host_uncertainty) {
    config$host <- matrix_norm_distribution(config$host_mean, config$host_sd)
  } else {
    config$host <- config$host_mean
  }
  
  # Update susceptible and total populations
  config$susceptible <- config$host - config$infected - exposed2
  config$susceptible[config$susceptible < 0] <- 0
  
  config$total_hosts <- config$host
  config$total_exposed <- exposed2
  
  # Update mortality tracker if mortality is on
  if (config$mortality_on) {
    mortality_tracker2 <- config$mortality_tracker
    mortality_tracker2[[length(mortality_tracker2)]] <- config$infected
    config$mortality_tracker <- mortality_tracker2
  }
  
  config$host <- NULL
  
  return(config)
}
