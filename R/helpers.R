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
  new_matrix <-
    round(matrix(mapply(function(x, y) {rnorm(x, y, n = 1)}, x = mean_matrix, y = sd_matrix),
                 nrow(mean_matrix), ncol(mean_matrix)), digits = 0)
  new_matrix[is.na(new_matrix)] <- 0
  new_matrix[new_matrix < 0] <- 0
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

competency_table_list_creator <- function(competency_table) {
  competency_table2 <- competency_table[, 1:(ncol(competency_table) - 1)]
  competencies <-
    rnorm(n = nrow(competency_table), mean = competency_table$competency_mean,
          sd = competency_table$compentency_sd)
  names(competency_table2)[ncol(competency_table2)] <- "competency"
  while (any(competencies > 1) || any(competencies < 0)) {
    competencies <-
      rnorm(n = nrow(competency_table),
            mean = competency_table$competency_mean, sd = competency_table$compentency_sd)
  }
  competency_table2$competency <- competencies
  competency_table2 <- competency_table2
  competency_table_list <- split(competency_table2, seq_len(nrow(competency_table2)))
  for (i in seq_along(competency_table_list)) {
    competency_table_list[[i]] <- unname(competency_table_list[[i]])
    competency_table_list[[i]] <- as.vector(t(competency_table_list[[i]]))
  }
  return(competency_table_list)
}

# Update host pools when uncertainties are used
host_pool_setup <- function(config) {
  total_infecteds <- config$zero_matrix
  total_exposeds <- config$zero_matrix
  total_hosts <- config$zero_matrix
  for (i in seq_along(config$host_file_list)) {
    host_pool <- config$host_pools[[i]]
    if (config$use_initial_condition_uncertainty) {
      infected <-
        matrix_norm_distribution(config$host_pool_infected_means[[i]],
                                 config$host_pool_infected_sds[[i]])
      while (any(infected < 0)) {
        infected <-
          matrix_norm_distribution(config$host_pool_infected_means[[i]],
                                   config$host_pool_infected_sds[[i]])
      }
      exposed2 <- matrix_norm_distribution(config$host_pool_exposed_means[[i]],
                                           config$host_pool_exposed_sds[[i]])
      while (any(exposed2 < 0)) {
        exposed2 <- matrix_norm_distribution(config$host_pool_exposed_means[[i]],
                                             config$host_pool_exposed_sds[[i]])
      }
      exposed <- host_pool$exposed
      exposed[[config$latency_period + 1]] <- exposed2
      host_pool$infected <- infected
      host_pool$exposed <- exposed
      host_pool$total_exposed <- exposed2

      total_infecteds <- total_infecteds + infected
      total_exposeds <- total_exposeds + exposed2
    }

    if (config$use_host_uncertainty) {
      host <- matrix_norm_distribution(config$host_pool_host_means[[i]],
                                       config$host_pool_host_sds[[i]])
      while (any(host > config$total_populations)) {
        host <- matrix_norm_distribution(config$host_pool_host_means[[i]],
                                         config$host_pool_host_sds[[i]])
      }
      host_pool$total_host <- host
      total_hosts <- total_hosts + host
    }

    susceptible <- host_pool$total_host - host_pool$infected - host_pool$total_exposed
    susceptible[susceptible < 0] <- 0
    host_pool$susceptible <- susceptible

    if (config$mortality_on) {
      mortality_tracker <- host_pool$mortality_tracker
      mortality_tracker[[length(mortality_tracker)]] <- host_pool$infected
      host_pool$mortality_tracker <- mortality_tracker
    }
    config$host_pools[[i]] <- host_pool
  }
  config$total_hosts <- total_hosts
  config$total_exposed <- total_exposeds
  config$total_infecteds <- total_infecteds

  return(config)
}

infected_rast_from_county <- function(county_infections, host, config) {
  infected_rast <- host
  names(host) <- c("host")
  terra::values(infected_rast) <- 0
  county_infections_rast <- terra::rasterize(county_infections, host, field = "FIPS")
  county_infections_rast
  fips_host_cells <- terra::extract(host, county_infections, cells = TRUE)
  fips_host_cells$FIPS <- terra::extract(county_infections_rast, county_infections)$FIPS

  for (id in unique(fips_host_cells$FIPS)) {
    sampler <- fips_host_cells[fips_host_cells$FIPS == id & fips_host_cells$host > 0, ]
    if (config$use_initial_condition_uncertainty) {
      inf_num <-
        round(rnorm(1, county_infections[county_infections$FIPS == id]$infected_mean,
                    county_infections[county_infections$FIPS == id]$infected_sd))
      while (inf_num < 0) {
        inf_num <-
          round(rnorm(1, county_infections[county_infections$FIPS == id]$infected_mean,
                      county_infections[county_infections$FIPS == id]$infected_sd))
      }
    } else {
      inf_num <- county_infections[county_infections$FIPS == id]$infected_mean
    }
    cells <- sample(sampler$cell, inf_num)
    infected_rast[cells] <- 1
  }
  return(infected_rast)
}
