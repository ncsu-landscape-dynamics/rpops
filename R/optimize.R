
prior_weight <- function(cost_column, potential_column) {
  minmax_scale <- function(column) {
    if (max(column) - min(column) > 0) {
      return(scales::rescale(column, to = c(0, 1)))
    } else {
      return(0)
    }
  }
  cost_norm <- minmax_scale(cost_column)
  potential_norm <- minmax_scale(potential_column)
  return((potential_norm + 1 - cost_norm) / 2)
}

pixel_treatments <- function(points, treatments_raster) {
  raster <- terra::rasterize(points, treatments_raster, background = 0)
  return(list(raster = raster, cost = sum(points$cost)))
}

buffer_treatments <- function(points, treatments_raster, cost_raster) {
  rasterization_factor <- 10
  buffers <- terra::buffer(points, width = points$buffer_size)
  high_res_raster <- terra::disagg(treatments_raster, fact = rasterization_factor)
  buffers_raster <- terra::rasterize(buffers, high_res_raster)
  count_raster <- terra::aggregate(buffers_raster,
    fact = rasterization_factor,
    fun = function(x) {
      length(x[!is.na(x)])
    }
  )
  treatments_raster <- count_raster / terra::global(count_raster, "max")[[1]]
  treatment_cost_raster <- treatments_raster * cost_raster
  actual_cost <- terra::global(treatment_cost_raster, "sum", na.rm = T)[[1]]
  return(list(raster = treatments_raster, cost = actual_cost))
}

treatments <- function(points, treatments_raster, cost_raster) {
  if ("buffer_size" %in% names(points)) {
    return(buffer_treatments(
      points,
      treatments_raster,
      cost_raster
    ))
  }
  return(pixel_treatments(points, treatments_raster))
}


run_pops <- function(config, treatment_raster = NULL) {
  if (!is.null(treatment_raster)) {
    config$treatment_maps <- list(terra::as.matrix(
      treatment_raster * config$pesticide_efficacy,
      wide = TRUE
    ))
  } else {
    config$treatment_maps <- list(matrix(0, nrow = 1, ncol = 1))
  }

  cl <- parallel::makeCluster(config$core_count)
  doParallel::registerDoParallel(cl)
  infected_stack <-
    foreach::foreach(
      i = seq_len(config$number_of_iterations),
      .combine = c,
      .packages = c("PoPS", "terra")
    ) %dopar% {
      config$random_seed <- round(stats::runif(1, 1, 1000000))
      data <- PoPS::pops_model(
        random_seed = config$random_seed,
        use_lethal_temperature = config$use_lethal_temperature,
        lethal_temperature = config$lethal_temperature,
        lethal_temperature_month = config$lethal_temperature_month,
        use_survival_rates = config$use_survival_rates,
        survival_rate_month = config$survival_rate_month,
        survival_rate_day = config$survival_rate_day,
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
        survival_rates = config$survival_rates,
        weather_coefficient = config$weather_coefficient,
        res = config$res,
        rows_cols = config$rows_cols,
        time_step = config$time_step,
        reproductive_rate = config$reproductive_rate,
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
        percent_natural_dispersal = config$percent_natural_dispersal,
        natural_distance_scale = config$natural_distance_scale,
        anthropogenic_distance_scale = config$anthropogenic_distance_scale,
        natural_dir = config$natural_dir,
        natural_kappa = config$natural_kappa,
        anthropogenic_dir = config$anthropogenic_dir,
        anthropogenic_kappa = config$anthropogenic_kappa,
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
        dispersal_stochasticity = config$dispersal_stochasticity,
        establishment_probability = config$establishment_probability,
        dispersal_percentage = config$dispersal_percentage,
        use_overpopulation_movements = config$use_overpopulation_movements,
        overpopulation_percentage = config$overpopulation_percentage,
        leaving_percentage = config$leaving_percentage,
        leaving_scale_coefficient = config$leaving_scale_coefficient,
        bbox = config$bounding_box,
        network_min_distance = config$network_min_distance,
        network_max_distance = config$network_max_distance,
        network_filename = config$network_filename,
        network_movement = config$network_movement
      )
      run <- c()
      run$infected_area <- data$area_infected
      run$quarantine_escape_distance <- data$quarantine_escape_distance
      run
    }
  stopCluster(cl)

  area_infected_runs <- infected_stack[seq(1, length(infected_stack), 2)]
  infected_area <- data.frame(t(rep(0, length(area_infected_runs[[1]]))))
  if (config$use_quarantine) {
    quarantine_escape_distance_runs <- infected_stack[seq(2, length(infected_stack), 2)]
    quarantine_escape_distances <- data.frame(t(rep(0, length(area_infected_runs[[1]]))))
  }
  for (p in seq_len(length(area_infected_runs))) {
    infected_area[p, ] <- area_infected_runs[[p]]
    if (config$use_quarantine) {
      quarantine_escape_distances[p, ] <- quarantine_escape_distance_runs[[p]]
    }
  }
  results <- list()
  infected_areas <- round(sapply(infected_area, mean, na.rm = TRUE), digits = 0)
  results$infected_area <- infected_areas[[length(infected_areas)]]
  if (config$use_quarantine) {
    quarantine_escape_distances[is.na(quarantine_escape_distances)] <- 0
    quarantine_distance <- round(sapply(quarantine_escape_distances, mean), digits = 0)
    results$quarantine_distance <- quarantine_distance[[length(quarantine_distance)]]
  }
  return(results)
}

pops_init <- function(config) {
  config$function_name <- "multirun"
  config$management <- FALSE
  config$use_initial_condition_uncertainty <- FALSE
  config$use_host_uncertainty <- FALSE
  config <- PoPS::configuration(config)

  if (!is.null(config$failure)) {
    stop(config$failure)
  }
  config$crs <- terra::crs(config$host)
  config$infected <- config$infected_mean
  exposed2 <- config$exposed_mean
  exposed <- config$exposed
  exposed[[config$latency_period + 1]] <- exposed2
  config$exposed <- exposed
  config$host <- config$host_mean
  susceptible <- config$host - config$infected - exposed2
  susceptible[susceptible < 0] <- 0
  config$susceptible <- susceptible
  config$total_hosts <- config$host
  config$total_exposed <- exposed2
  if (config$mortality_on) {
    mortality_tracker2 <- config$mortality_tracker
    mortality_tracker2[[length(mortality_tracker2)]] <- config$infected
    config$mortality_tracker <- mortality_tracker2
  }
  return(config)
}

best_guess <- function(points,
                       weight_column,
                       treatments_raster,
                       cost_raster,
                       budget,
                       config) {
  sorted <- points[order(points[[weight_column]][[1]], decreasing = TRUE), ]
  candidate <- sorted[cumsum(sorted$cost) <= budget, ]

  treatment <- treatments(
    candidate,
    treatments_raster,
    cost_raster
  )
  result <- run_pops(config, treatment$raster)
  return(list(candidate = candidate, result = result))
}

estimate_baseline <- function(config) {
  # run with 0 spread to get initial quarantine distance
  config_no_spread <- config
  config_no_spread$parameter_means[1] <- 0
  config_no_spread$reproductive_rate <- 0
  baseline <- run_pops(config_no_spread)
  quarantine_distance <- baseline$quarantine_distance
  baseline <- run_pops(config)
  baseline$quarantine_distance <- quarantine_distance
  return (baseline)
}

estimate_initial_threshold <- function(points,
                                       weight_column,
                                       treatments_raster,
                                       cost_raster,
                                       budget,
                                       baseline,
                                       score_weights,
                                       config) {
  runs <- 10
  results_list <- list()
  for (run in seq(runs)) {
    candidate <- sample_candidate(points, weight_column, budget)
    treatment <- treatments(
      candidate,
      treatments_raster,
      cost_raster
    )
    results_list[[run]] <- run_pops(config, treatment$raster)
  }
  scores <- sapply(results_list, scoring, baseline, score_weights)
  threshold <- quantile(scores, probs = 0.1)
  threshold_step <- abs(threshold - quantile(scores, probs = 0.2))
  return(list(threshold = threshold, threshold_step = threshold_step))
}

scoring <- function(simulated, baseline, weights=c(1, 1)) {
  scores <- c(NA, NA)

  if (!is.null(simulated$infected_area)) {
    scores[1] <- simulated$infected_area / baseline$infected_area
  }
  if (!is.null(simulated$quarantine_distance)) {
    scores[2] <- (baseline$quarantine_distance - simulated$quarantine_distance) /
      baseline$quarantine_distance
  }
  return (weighted.mean(scores, w=weights, na.rm = TRUE))
}

sample_candidate <- function(points, weight_column, budget) {
  # to allow sample vector size 1
  sample.vec <- function(x, ...) x[sample(length(x), ...)]
  tmp_budget <- 0
  tries <- 0
  count <- 1
  candidate_cats <- c()
  cats <- points$cat
  weights <- points[[weight_column]][[1]]
  costs <- points$cost
  while (length(cats)) {
    sample_cat <- sample.vec(cats, 1, prob = weights)
    selected_index <- which(cats == sample_cat)
    cost <- costs[[selected_index]]
    if (tmp_budget + cost <= budget) {
      tmp_budget <- tmp_budget + cost
      candidate_cats[count] <- sample_cat
      count <- count + 1
      cats <- cats[-selected_index]
      weights <- weights[-selected_index]
      costs <- costs[-selected_index]
    } else {
      tries <- tries + 1
      if (tries > 50) {
        break
      }
    }
  }
  return(points[points$cat %in% candidate_cats, ])
}


generation <- function(points,
                       weight_column,
                       treatments_raster,
                       cost_raster,
                       budget,
                       min_particles,
                       threshold_percentile,
                       threshold,
                       threshold_step,
                       baseline,
                       score_weights,
                       config) {
  score_list <- numeric(min_particles)
  particle_count <- 0
  tested <- 0
  best <- list()
  best$simulation <- baseline
  best$score <- 1
  new_weights <- setNames(as.list(rep(0, length(points$cat))), points$cat)
  while (particle_count < min_particles) {
    tested <- tested + 1
    candidate <- sample_candidate(points, weight_column, budget)
    treatment <- treatments(
      candidate,
      treatments_raster,
      cost_raster
    )
    simulation <- run_pops(config, treatment$raster)
    score <- scoring(simulation, baseline, score_weights)
    # success, add particle
    if (score < threshold) {
      particle_count <- particle_count + 1
      for (cat in as.character(candidate$cat)) {
        new_weights[[cat]] <- new_weights[[cat]] + 1
      }
      score_list[particle_count] <- score
      if (score < best$score) {
        best$simulation <- simulation
        best$candidate <- candidate
        best$score <- score
        best$cost <- treatment$cost
        best$treatment <- treatment$raster
      }
    }
    acceptance_rate <- particle_count / tested
    particles_left <- min_particles - particle_count
    # print intermediate
    if (tested %% 10 == 0) {
      message(
        "Rate: ", acceptance_rate,
        ", particles left: ", particles_left,
        ", best score: ", best$score
      )
    }
    if ((tested == 50) && (acceptance_rate < 0.05 || acceptance_rate > 0.15)) {
      if (acceptance_rate < 0.05) {
        threshold <- threshold + threshold_step
      } else {
        threshold <- threshold - threshold_step
      }
      particle_count <- 0
      tested <- 0
      best$score <- 1
      new_weights <- setNames(as.list(rep(0, length(points$cat))), points$cat)
    }
  }
  new_threshold <- quantile(score_list, probs = threshold_percentile / 100)
  threshold_step <- abs(new_threshold - quantile(score_list, probs = (threshold_percentile + 10) / 100))
  output <- list(
    weights = new_weights,
    acceptance_rate = acceptance_rate,
    threshold = new_threshold,
    threshold_step = threshold_step,
    best = best
  )
  return(output)
}

filter_particles <- function(points, weight_column, iteration, percentile) {
  # filter unsuccessful ones
  filtered <- points[points[[weight_column]] > 0, ]
  weights_percentile <- quantile(filtered[[weight_column]][[1]],
    percentile / 100,
    names = FALSE
  )
  # filter unlikely ones
  points[points[[weight_column]] <= weights_percentile, weight_column] <- 0
  return(points)
}

#' @title Optimize treatments using ABC
#'
#' @description Optimizes treatment locations using ABC based on infestation potential and cost.
#' Treatments are pixel based with the option to specify constant or variable buffer size.
#' The criteria can be infected area, distance to quarantine boundary or both.
#'
#'
#' @inheritParams pops_multirun
#' @param infestation_potential_file Raster file with infestation potential.
#' @param cost_file Raster file with cost of treatment
#'  (constant value or spatially variable). When buffers are used,
#'  the cost must represent the cost of the buffer around a pixel.
#' @param buffer_size_file Raster file with size of buffers
#' (constant values or spatially variable).
#' @param min_particles Number of successful combinations per generation.
#' @param budget Budget to spend. Limits the number of pixels that can be treated.
#' @param score_weights Weights to compute weighted average of infected_area and 
#' distance to quarantine boundary success metrics. For example, c(1, 0) means only
#' infected area is used, c(2, 1) means both metrics are averaged with provided weights.
#' @param filter_percentile Lower value removes fewer pixels from the pool,
#' resulting in more iterations and lower acceptance rate, but potentially better results.
#' @param threshold_percentile Determines threshold for next iteration.

#' @useDynLib PoPS, .registration = TRUE
#' @return results
#' @export
#'
optimize <- function(infestation_potential_file,
                     cost_file,
                     buffer_size_file = NULL,
                     min_particles,
                     budget,
                     score_weights = c(1, 1),
                     filter_percentile = 15,
                     threshold_percentile = 10,
                     infected_file,
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
                     use_survival_rates = FALSE,
                     survival_rate_month = 3,
                     survival_rate_day = 15,
                     survival_rates_file = "",
                     use_lethal_temperature = FALSE,
                     temperature_file = "",
                     lethal_temperature = -12.87,
                     lethal_temperature_month = 1,
                     mortality_on = FALSE,
                     mortality_rate = 0,
                     mortality_time_lag = 0,
                     mortality_frequency = "year",
                     mortality_frequency_n = 1,
                     treatment_dates = c(""),
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
                     output_frequency = "final_step",
                     movements_file = "",
                     use_movements = FALSE,
                     start_exposed = FALSE,
                     generate_stochasticity = TRUE,
                     establishment_stochasticity = TRUE,
                     movement_stochasticity = TRUE,
                     dispersal_stochasticity = TRUE,
                     establishment_probability = 0.5,
                     dispersal_percentage = 0.99,
                     quarantine_areas_file = "",
                     quarantine_directions = "",
                     use_quarantine = FALSE,
                     use_spreadrates = FALSE,
                     use_overpopulation_movements = FALSE,
                     overpopulation_percentage = 0,
                     leaving_percentage = 0,
                     leaving_scale_coefficient = 1,
                     exposed_file = "",
                     mask = NULL,
                     write_outputs = "None",
                     output_folder_path = "",
                     network_filename = "",
                     network_movement = "walk") {
  # parameters for pops
  config <- as.list(environment())

  # extract infected points and associated cost and potential
  infected_raster <- terra::rast(infected_file)
  infected_points <- terra::as.points(infected_raster)
  names(infected_points) <- "infected"
  infected_points <- infected_points[infected_points$infected > 0]
  # cat
  infected_points$cat <- 1:nrow(infected_points)
  # cost
  cost_raster <- terra::rast(cost_file)
  names(cost_raster) <- "cost"
  infected_points <- terra::extract(cost_raster, infected_points, bind = TRUE)
  # infestation potential
  potential_raster <- terra::rast(infestation_potential_file)
  names(potential_raster) <- "potential"
  infected_points <- terra::extract(potential_raster, infected_points, bind = TRUE)

  if (!is.null(buffer_size_file)) {
    buffer_size_raster <- terra::rast(buffer_size_file)
    names(buffer_size_raster) <- "buffer_size"
    infected_points <- terra::extract(buffer_size_raster, infected_points, bind = TRUE)
  }
  buffer_size_raster <- NULL

  # prior weights
  iteration <- 1
  weight_column <- paste0("weight_", iteration)
  infected_points[[weight_column]] <- prior_weight(
    infected_points$cost,
    infected_points$potential
  )

  # prepare treatment raster
  temporary_directory <- tempdir()
  treatments_raster <- terra::rast(terra::ext(infected_raster),
    resolution = terra::res(infected_raster)
  )
  terra::crs(treatments_raster) <- terra::crs(infected_raster)
  config <- pops_init(config)

  # baseline
  baseline <- estimate_baseline(config)
  message("Baseline area:", baseline$infected_area)
  if (score_weights[2] > 0) {
  message(
    "Initial distance to quarantine boundary:",
    baseline$quarantine_distance
  )
  }

  # best guess
  best_guess <- best_guess(
    infected_points,
    weight_column,
    treatments_raster,
    cost_raster,
    budget,
    config
  )
  message("Best guess infected area:", best_guess$result$infected_area)
  if (score_weights[2] > 0) {
  message(
    "Best guess distance to quarantine boundary:",
    best_guess$result$quarantine_distance
  )
  }

  # initial threshold
  thresholds <- c()

  initial_threshold <- estimate_initial_threshold(
    infected_points,
    weight_column,
    treatments_raster,
    cost_raster,
    budget,
    baseline,
    score_weights,
    config
  )
  thresholds[1] <- initial_threshold$threshold
  threshold_step <- initial_threshold$threshold_step

  filtered_points <- infected_points
  tmp_points <- infected_points
  acceptance_rates <- c()
  while (TRUE) {
    message(
      "Iteration ", iteration,
      ", threshold ", thresholds[length(thresholds)]
    )
    results <- generation(
      filtered_points,
      weight_column,
      treatments_raster,
      cost_raster,
      budget,
      min_particles,
      threshold_percentile,
      thresholds[length(thresholds)],
      threshold_step,
      baseline,
      score_weights,
      config
    )
    acceptance_rates <- append(acceptance_rates, results$acceptance_rate)
    thresholds <- append(thresholds, results$threshold)
    threshold_step <- results$threshold_step
    new_weight_column <- paste0("weight_", iteration + 1)
    new_weights <- results$weights
    filtered_points[[new_weight_column]] <- results$weights[match(
      filtered_points$cat,
      as.integer(names(results$weights))
    )]
    # filtering
    tmp <- filter_particles(
      filtered_points,
      new_weight_column,
      iteration,
      filter_percentile
    )
    before <- filtered_points[filtered_points[[new_weight_column]] > 0, ]
    after <- tmp[tmp[[new_weight_column]] > 0, ]
    if (sum(after$cost) >= budget) {
      message(
        "Filtered ", length(before) - length(after),
        ", remaining: ", length(after)
      )
      filtered_points <- tmp
      weight_column <- new_weight_column
    } else {
      break
    }
    iteration <- iteration + 1
  }

  output <- list(
    best_guess = best_guess,
    best = results$best,
    weights = filtered_points,
    acceptance_rates = acceptance_rates,
    thresholds = thresholds
  )

  cat("Baseline area:", baseline$infected_area, "\n")
  cat("Best guess infected area:", best_guess$result$infected_area, "\n")
  cat("Optimized infected area:", results$best$simulation$infected_area, "\n")

  if (score_weights[2] > 0) {
    cat(
      "Initial distance to quarantine boundary:",
      baseline$quarantine_distance, "\n"
    )
    cat(
      "Best guess distance to quarantine boundary:",
      best_guess$result$quarantine_distance, "\n"
    )
    cat(
      "Optimized distance to quarantine boundary:",
      results$best$simulation$quarantine_distance, "\n"
    )
  }
  cat("Acceptance rates: ", acceptance_rates, "\n")
  cat("Thresholds: ", thresholds, "\n")
  if (file.exists(output_folder_path)) {
    file_name <- file.path(output_folder_path, "best_candidate.gpkg")
    terra::writeVector(results$best$candidate, file_name, overwrite = TRUE)
    file_name <- file.path(output_folder_path, "best_treatment.tif")
    terra::writeRaster(results$best$treatment, file_name, overwrite = TRUE)
    file_name <- file.path(output_folder_path, "best_guess_candidate.gpkg")
    terra::writeVector(best_guess$result$candidate, file_name, overwrite = TRUE)
    file_name <- file.path(output_folder_path, "output.rdata")
    save(output, file = file_name)
  }
  return(output)
}
