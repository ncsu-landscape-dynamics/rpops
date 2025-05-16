# These functions are designed to reduce code complexity and the need to copy
# and paste code across main functions
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
    format(config$acceptance_rate, digits = 5),
    sep = " "
  )

  if (config$use_quantity) {
    config$acceptance_rate_info <- paste(
      config$acceptance_rate_info,
      "\nquantity:                ",
      config$quantity,
      "\nquantity threshold:      ",
      config$allocation_threshold,
      sep = " "
    )
  }

  if (config$use_allocation) {
    config$acceptance_rate_info <- paste(
      config$acceptance_rate_info,
      "\nallocation:              ",
      config$allocation,
      "\nallocation threshold:    ",
      config$allocation_threshold,
      sep = " "
    )
  }

  if (config$use_configuration) {
    config$acceptance_rate_info <- paste(
      config$acceptance_rate_info,
      "\nconfiguration:           ",
      config$configuration_dis,
      "\nconfiguration threshold: ",
      config$configuration_threshold,
      sep = " "
    )
  }

  if (config$use_accuracy) {
    config$acceptance_rate_info <- paste(
      config$acceptance_rate_info,
      "\naccuracy:                ",
      config$accuracy,
      "\naccuracy threshold:      ",
      config$accuracy_threshold,
      sep = " "
    )
  }

  if (config$use_precision) {
    config$acceptance_rate_info <- paste(
      config$acceptance_rate_info,
      "\nprecision:               ",
      config$precision,
      "\nprecision threshold:     ",
      config$precision_threshold,
      sep = " "
    )

  }

  if (config$use_recall) {
    config$acceptance_rate_info <- paste(
      config$acceptance_rate_info,
      "\nrecall:                  ",
      config$recall,
      "\nrecall threshold:        ",
      config$recall_threshold,
      sep = " "
    )
  }

  if (config$use_specificity) {
    config$acceptance_rate_info <- paste(
      config$acceptance_rate_info,
      "\nspecificity:             ",
      config$specificity,
      "\nspecificity threshold:   ",
      config$specificity_threshold,
      sep = " "
    )
  }

  if (config$use_rmse) {
    config$acceptance_rate_info <- paste(
      config$acceptance_rate_info,
      "\nrmse:                    ",
      config$rmse,
      "\nrmse threshold:          ",
      config$rmse_threshold,
      sep = " "
    )
  }

  if (config$use_distance) {
    config$acceptance_rate_info <- paste(
      config$acceptance_rate_info,
      "\ndistance difference:     ",
      config$distance_difference,
      "\ndistance threshold:      ",
      config$distance_threshold,
      sep = " "
    )
  }

  if (config$use_mcc) {
    config$acceptance_rate_info <- paste(
      config$acceptance_rate_info,
      "\nMCC:                     ",
      config$mcc,
      "\nMCC threshold:           ",
      config$mcc_threshold,
      sep = " "
    )
  }

  config$acceptance_rate_info <- paste(config$acceptance_rate_info, "\n\n", sep = " ")
  return(config)
}

draw_parameters <- function(config) {
  parameters <-
    MASS::mvrnorm(1, config$parameter_means, config$parameter_cov_matrix)
  while (any(
    parameters[1] < 0 |
    parameters[2] <= 0 |
    parameters[3] > 1 |
    parameters[3] <= 0 |
    parameters[4] <= 0 |
    parameters[5] < 0 |
    parameters[6] < 0 |
    parameters[7] < config$res$ew_res / 2 |
    parameters[7] > parameters[8] |
    parameters[8] >
    min(config$rows_cols$num_cols, config$rows_cols$num_rows) * config$res$ew_res
  )) {
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
                        config$res$ew_res])

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
      MASS::mvrnorm(config$number_of_draws,
                    config$parameter_means,
                    config$parameter_cov_matrix)
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
    data.frame(
      disperser_generation = sample(1:999999999, n, replace = FALSE),
      natural_dispersal = sample(1:999999999, n, replace = FALSE),
      anthropogenic_dispersal = sample(1:999999999999, n, replace = FALSE),
      establishment = sample(1:999999999, n, replace = FALSE),
      weather = sample(1:999999999, n, replace = FALSE),
      movement = sample(1:999999999, n, replace = FALSE),
      overpopulation = sample(1:999999999, n, replace = FALSE),
      survival_rate = sample(1:999999999, n, replace = FALSE),
      soil = sample(1:999999999, n, replace = FALSE),
      iteration = seq(1, n, by = 1)
    )

  return(random_seeds)
}

# creates a matrix from a matrix of mean values and a matrix of standard deviations. The two
# matrices must be the same size.
matrix_norm_distribution <- function(mean_matrix, sd_matrix) {
  new_matrix <- matrix(round(rnorm(
    length(mean_matrix), mean = mean_matrix[], sd = sd_matrix[]
  )),
  nrow = nrow(mean_matrix),
  ncol = ncol(mean_matrix))
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
    (((n1 * n2) * (m1 - m2) ^ 2) / ((n1 + n2) * (n1 + n2 - 1)))
}

# Create a unique identifier for the pops_run_lite output filenames
generate_uid <- function() {
  letters_part <- sample(letters[1:26], 6, replace = TRUE)
  numbers_part <- sample(0:9, 6, replace = TRUE)

  uid <- paste0(c(rbind(letters_part, numbers_part)), collapse = "")
  return(uid)
}

# Add time step for the pops_run_lite output filenames
format_output_date <- function(start_date,
                               output_frequency,
                               time_step,
                               multiplier = 0) {
  # Ensure the start date is in Date format and check for validity
  if (is.na(as.Date(start_date, format = "%Y-%m-%d"))) {
    stop("Invalid start_date provided. Please use a valid date in 'YYYY-MM-DD' format.")
  }
  start_date <- as.Date(start_date)

  # Function to increment the date based on a unit and multiplier
  advance_date <- function(date, unit, multiplier) {
    switch(
      unit,
      "year" = date + lubridate::years(multiplier),
      "month" = date + months(multiplier),
      "week" = date + lubridate::weeks(multiplier),
      "day" = date + lubridate::days(multiplier)
    )
  }

  # Determine the correct format and advance the date
  formatted_date <- if (output_frequency %in% c("year", "month")) {
    format(advance_date(start_date, output_frequency, multiplier),
           if (output_frequency == "year")
             "%Y"
           else
             "%Y-%m")
  } else if (output_frequency == "week") {
    paste(format(advance_date(start_date, "week", multiplier), "%Y-%m-%d"), "w")
  } else if (output_frequency == "day") {
    paste(format(advance_date(start_date, "day", multiplier), "%Y-%m-%d"), "d")
  } else if (output_frequency == "time step" ||
             output_frequency == "every_n_steps") {
    unit <- switch(
      time_step,
      "year" = "year",
      "month" = "month",
      "week" = "week",
      "day" = "day"
    )
    advanced_date <- advance_date(start_date, unit, multiplier)
    formatted_date <- switch(
      unit,
      "year" = format(advanced_date, "%Y"),
      "month" = format(advanced_date, "%Y-%m"),
      "week" = paste(format(advanced_date, "%Y-%m-%d"), "w"),
      "day" = paste(format(advanced_date, "%Y-%m-%d"), "d")
    )
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
  gc()

  # Update host based on uncertainty settings
  if (config$use_host_uncertainty) {
    config$host <- matrix_norm_distribution(config$host_mean, config$host_sd)
  } else {
    config$host <- config$host_mean
  }
  gc()

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

# Function to convert  and export infected matrices from pops_lite.R into rasters
write_infected_rasters <- function(config, uid) {
  for (i in seq_len(config$number_of_iterations)) {
    # Read the template raster
    r_template <- rast(config$infected_file)
    # Assign values to the raster
    terra::values(r_template) <- readRDS(file.path(config$output_folder_path,
                                                   paste0(uid, "_", i, ".rds")))$infected[[1]]

    # Write the raster to disk with compression
    terra::writeRaster(
      r_template,
      filename = file.path(config$output_folder_path, paste0(uid, "_", i, ".tif")),
      overwrite = TRUE,
      datatype = "INT2U",
      gdal = c("COMPRESS=ZSTD")
    )
    gc()

    # Print a status message
    cat("Processed iteration",
        i,
        "- output saved to:",
        config$output_folder_path,
        "\n")
  }
}


#' Function to dynamically update the working directory for all file paths in
#' the configuration rds output

update_config_paths <- function(config, root) {
  # Normalize the root directory
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  target_dir <- basename(root)  # Extract the target directory name

  # List of keys to update
  keys_to_update <- c(
    "infected_file", "host_file", "total_populations_file",
    "temperature_coefficient_file", "precipitation_coefficient_file",
    "temperature_file", "survival_rates_file", "treatments_file",
    "movements_file", "quarantine_areas_file", "exposed_file",
    "network_filename", "temperature_coefficient_sd_file",
    "precipitation_coefficient_sd_file", "soil_starting_pest_file",
    "output_folder_path", "reference", "comparison", "ref", "comp"
  )

  # Precompute the replacement pattern
  replace_pattern <- paste0(".*", target_dir)

  # Update paths in place
  for (key in keys_to_update) {
    if (!is.null(config[[key]]) && is.character(config[[key]])) {
      old_path <- normalizePath(config[[key]], winslash = "/", mustWork = FALSE)
      if (grepl(target_dir, old_path, fixed = TRUE)) {
        config[[key]] <- sub(replace_pattern, root, old_path)
      }
    }
  }

  return(config)
}

# Reformat competency_table into list (per host composition) with competency values
# randomly sampled from a normal distribution using mean and sd in competency table

competency_table_list_creator <- function(competency_table) {
  competency_table2 <- competency_table[, 1:(ncol(competency_table) - 1)]
  competencies <-
    rnorm(n = nrow(competency_table), mean = competency_table$competency_mean,
          sd = competency_table$competency_sd)
  names(competency_table2)[ncol(competency_table2)] <- "competency"
  while (any(competencies > 1) || any(competencies < 0)) {
    competencies <-
      rnorm(n = nrow(competency_table),
            mean = competency_table$competency_mean, sd = competency_table$competency_sd)
  }
  competency_table2$competency <- competencies
  competency_table_list <- split(competency_table2, seq_len(nrow(competency_table2)))
  for (i in seq_along(competency_table_list)) {
    competency_table_list[[i]] <- unname(competency_table_list[[i]])
    competency_table_list[[i]] <- as.vector(t(competency_table_list[[i]]))
  }
  return(competency_table_list)
}

# Reformat pest_host_table into list (per host species) with susceptibility and
# mortality rates randomly sampled from a normal distribution using the means and sds i
# in the pest_host_table

pest_host_table_list_creator <- function(pest_host_table) {
  pest_host_table2 <- pest_host_table[, !grepl("_sd", colnames(pest_host_table))]
  susceptibilities <-
    rnorm(n = nrow(pest_host_table), mean = pest_host_table$susceptibility_mean,
          sd = pest_host_table$susceptibility_sd)
  names(pest_host_table2)[1] <- "susceptibility"
  while (any(susceptibilities > 1) || any(susceptibilities < 0)) {
    susceptibilities <- rnorm(n = nrow(pest_host_table),
                              mean = pest_host_table$susceptibility_mean,
                              sd = pest_host_table$susceptibility_sd)
  }
  pest_host_table2$susceptibility <- susceptibilities
  mortality_rates <-
    rnorm(n = nrow(pest_host_table), mean = pest_host_table$mortality_rate_mean,
          sd = pest_host_table$mortality_rate_sd)
  names(pest_host_table2)[2] <- "mortality_rate"
  while (any(mortality_rates > 1) || any(mortality_rates < 0)) {
    mortality_rates <-
      rnorm(n = nrow(pest_host_table), mean = pest_host_table$mortality_rate_mean,
            sd = pest_host_table$mortality_rate_sd)
  }
  pest_host_table2$mortality_rate <- mortality_rates
  pest_host_table_list <- split(pest_host_table2, seq_len(nrow(pest_host_table2)))
  for (i in seq_along(pest_host_table_list)) {
    pest_host_table_list[[i]] <- unname(pest_host_table_list[[i]])
    pest_host_table_list[[i]] <- as.vector(t(pest_host_table_list[[i]]))
  }
  return(pest_host_table_list)
}

# Update host pools when uncertainties are used
host_pool_setup <- function(config) {
  total_infecteds <- config$zero_matrix
  total_exposeds <- config$zero_matrix
  total_hosts <- config$zero_matrix
  for (i in seq_along(config$host_file_list)) {
    host_pool <- config$host_pools[[i]]

    if (config$use_host_uncertainty) {
      host <- matrix_norm_distribution(config$host_pool_host_means[[i]],
                                       config$host_pool_host_sds[[i]])
      while (any(host > config$total_populations, na.rm = TRUE)) {
        host <- matrix_norm_distribution(config$host_pool_host_means[[i]],
                                         config$host_pool_host_sds[[i]])
      }
      host_pool$total_host <- host
      total_hosts <- total_hosts + host
    } else {
      host <- host_pool$total_hosts
      total_hosts <- total_hosts + host_pool$total_hosts
    }

    if (config$use_initial_condition_uncertainty) {
      if (config$county_level_infection_data) {
        host_rast <- terra::rast(config$host_file_list[i])[[1]]
        terra::values(host_rast) <- host
        county_infections <- terra::vect(config$infected_file_list[i])
        infected <- infected_rast_from_county(county_infections, host_rast[[1]], config)
        infected <- terra::as.matrix(infected[[1]], wide = TRUE)
        host_pool$infected <- infected
        if (config$model_type == "SEI" && config$start_exposed) {
          county_exposeds <- terra::vect(config$exposed_file_list[i])
          exposed2 <- infected_rast_from_county(county_exposeds, host_rast[[1]], config)
          exposed2 <- terra::as.matrix(exposed2[[1]], wide = TRUE)
          exposed <- host_pool$exposed
          exposed[[config$latency_period + 1]] <- exposed2
          host_pool$exposed <- exposed
          host_pool$total_exposed <- exposed2
        } else {
          exposed2 <- host_pool$total_exposed
        }

      } else {
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
      }

      total_infecteds <- total_infecteds + infected
      total_exposeds <- total_exposeds + exposed2
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
    if (nrow(sampler) > 0) {
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
  }
  return(infected_rast)
}


calculated_stats_county_level <- function(compare_vect) {
  compare_vect$true_positives <-
    ifelse(compare_vect$reference > 0 & compare_vect$comparison > 0, 1, 0)
  compare_vect$true_negatives <-
    ifelse(compare_vect$reference == 0 & compare_vect$comparison == 0, 1, 0)
  compare_vect$false_negatives <-
    ifelse(compare_vect$reference > 0 & compare_vect$comparison == 0, 1, 0)
  compare_vect$false_positives <-
    ifelse(compare_vect$reference == 0 & compare_vect$comparison > 0, 1, 0)
  compare_vect$unknown_negatives <-
    ifelse(is.na(compare_vect$reference) & compare_vect$comparison == 0, 1, 0)
  compare_vect$unknown_positives <-
    ifelse(is.na(compare_vect$reference) & compare_vect$comparison > 0, 1, 0)
  compare_vect$infected_difference <- compare_vect$reference - compare_vect$comparison

  compare_df <- as.data.frame(compare_vect)
  output <- as.data.frame(compare_vect[1, ])
  output[1, ] <- colSums(compare_df)

  output$total_obs <-
    output$true_negative + output$true_positive + output$false_negative + output$false_positive
  output$accuracy <- (output$true_negative + output$true_positive) / output$total_obs
  output$precision <- output$true_positive / (output$true_positive + output$false_positive)
  output$recall <- output$true_positive / (output$true_positive + output$false_negative)
  output$specificity <- output$true_negative / (output$true_negative + output$false_positive)

  output$tp_fp <- as.double((output$true_positive + output$false_positive))
  output$tp_fn <- as.double((output$true_positive + output$false_negative))
  output$tn_fp <- as.double((output$true_negative + output$false_positive))
  output$tn_fn <- as.double((output$true_negative + output$false_negative))

  if (is.nan(output$tp_fp) || output$tp_fp == 0) {output$tp_fp <- 1}
  if (is.nan(output$tp_fn) || output$tp_fn == 0) {output$tp_fn <- 1}
  if (is.nan(output$tn_fp) || output$tn_fp == 0) {output$tn_fp <- 1}
  if (is.nan(output$tn_fn) || output$tn_fn == 0) {output$tn_fn <- 1}

  output$mcc <-
    ((output$true_positive * output$true_negative) -
       (output$false_positive * output$false_negative)) /
    sqrt(output$tp_fp * output$tp_fn * output$tn_fp * output$tn_fn)

  if (is.nan(output$accuracy)) {output$accuracy <- 0}
  if (is.nan(output$precision)) {output$precision <- 0}
  if (is.nan(output$recall)) {output$recall <- 0}
  if (is.nan(output$specificity)) {output$specificity <- 0}

  output$rmse <- Metrics::rmse(compare_df$reference, compare_df$comparison)

  if (output$false_negative == 0 && output$false_positive == 0) {
    output$odds_ratio <- (output$true_positive * output$true_negative) / 1
  } else if (output$false_negative == 0) {
    output$odds_ratio <- (output$true_positive * output$true_negative) / output$false_positive
  } else if (output$false_positive == 0) {
    output$odds_ratio <- (output$true_positive * output$true_negative) / output$false_negative
  } else {
    output$odds_ratio <- (output$true_positive * output$true_negative) /
      (output$false_negative * output$false_positive)
  }


  output$residual_error <- sum(abs(compare_df$reference -  compare_df$comparison))

  return(output)
}


calculate_all_stats <- function(config, data) {
  all_disagreement <-
    foreach::foreach(
      q = seq_len(length(data$host_pools[[1]]$infected)), .combine = rbind,
      .packages = c("terra", "PoPS")
    ) %do% {
      # need to assign reference, comparison, and mask in inner loop since
      # terra objects are pointers

      comparison <- terra::rast(config$host_file_list[[1]])[[1]]
      terra::values(comparison) <- 0
      reference <- comparison
      mask <- comparison
      infections <- comparison
      for (p in seq_len(length(data$host_pools))) {
        terra::values(infections) <- data$host_pools[[p]]$infected[[q]]
        comparison <- comparison + infections
      }
      terra::values(mask) <- config$mask_matrix
      if (config$county_level_infection_data) {
        reference <- terra::vect(config$infected_years_file[[1]])
        compare_vect <- reference[, c(1, (q + 1))]
        names(compare_vect) <- c("FIPS", "reference")
        compare_vect$comparison <- terra::extract(comparison, reference, fun = "sum")[, 2]
        ad <- calculated_stats_county_level(compare_vect)
        ad <- calculated_stats_county_level(compare_vect)
        ad$quantity_disagreement <- 0
        ad$allocation_disagreement <- 0
        ad$allocation_disagreement <- 0
        ad$configuration_disagreement <- 0
        ad$distance_difference <- 0

      } else {
        terra::values(reference) <- config$infection_years2[[q]]
        ad <-
          quantity_allocation_disagreement(reference,
                                           comparison,
                                           use_configuration = config$use_configuration,
                                           mask = mask,
                                           use_distance = config$use_distance)
        if (file.exists(config$point_file)) {
          obs_data <- terra::vect(config$point_file)
          obs_data <- terra::project(obs_data, comparison)
          s <- extract(comparison, obs_data)
          names(s) <- c("ID", paste("sim_value_output_", q, sep = ""))
          s <- s[2]
          obs_data <- cbind(obs_data, s)
          ## calculate true positive, true negatives, false positives, false
          ## negatives, and other statistics and add them to the data frame
          ## for export
          ad$points_true_positive <-
            nrow(obs_data[obs_data$positive > 0 & obs_data$sim_value_output_1 > 0, ])
          ad$points_false_negative <-
            nrow(obs_data[obs_data$positive > 0 & obs_data$sim_value_output_1 == 0, ])
          ad$points_false_positive <-
            nrow(obs_data[obs_data$positive == 0 & obs_data$sim_value_output_1 > 0, ])
          ad$points_true_negative <-
            nrow(obs_data[obs_data$positive == 0 & obs_data$sim_value_output_1 == 0, ])
          ad$points_total_obs <-
            ad$points_true_negative + ad$points_true_positive +
            ad$points_false_negative + ad$points_false_positive
          ad$points_accuracy <-
            (ad$points_true_negative + ad$points_true_positive) / ad$points_total_obs
          ad$points_precision <-
            ad$points_true_positive / (ad$points_true_positive + ad$points_false_positive)
          ad$points_recall <-
            ad$points_true_positive / (ad$points_true_positive + ad$points_false_negative)
          ad$points_specificiity <-
            ad$points_true_negative / (ad$points_true_negative + ad$points_false_positive)
        }
      }
      ad$output <- q
      ad
    }
  all_disagreement <- data.frame(all_disagreement)
  return(all_disagreement)
}
