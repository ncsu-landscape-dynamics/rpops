#' PoPS (Pest or Pathogen Spread) model Multiple Runs
#'
#' A dynamic species distribution model for pest or pathogen spread in forest
#' or agricultural ecosystems. The model is process based meaning that it uses
#' understanding of the effect of weather and other environmental factors on
#' reproduction and survival of the pest/pathogen in order to forecast spread
#' of the pest/pathogen into the future. Run multiple stochastic simulations,
#' propagating uncertainty in parameters, initial conditions, and drivers.
#' The model is process based meaning that it uses understanding of the effect
#' of weather on reproduction and survival of the pest/pathogen in order to
#' forecast spread of the pest/pathogen into the future.
#'
#' @inheritParams pops
#' @param number_of_iterations how many iterations do you want to run to allow the calibration to
#' converge at least 10
#' @param number_of_cores enter how many cores you want to use (default = NA). If not set uses the
#' # of CPU cores - 1. must be an integer >= 1
#' @param write_outputs Either c("summary_outputs", "all_simulations", or "None"). If not
#' "None" output folder path must be provided.
#' @param output_folder_path this is the full path with either / or \\ (e.g.,
#' "C:/user_name/desktop/pops_sod_2020_2023/outputs/")
#'
#' @importFrom terra app rast xres yres classify extract ext as.points ncol nrow project
#' nlyr rowFromCell colFromCell values as.matrix rowFromCell colFromCell crs vect
#' @importFrom stats runif rnorm median sd
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  registerDoSEQ %dopar% %do%
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom lubridate interval time_length mdy %within%
#' @importFrom utils write.csv read.table read.csv
#' @importFrom methods is
#'
#' @return list of infected and susceptible per year
#' @export
#'
pops_multirun <- function(infected_file_list,
                          host_file_list,
                          total_populations_file,
                          parameter_means,
                          parameter_cov_matrix,
                          pest_host_table,
                          competency_table,
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
                          mortality_frequency = "year",
                          mortality_frequency_n = 1,
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
                          dispersal_stochasticity = TRUE,
                          establishment_probability = 0.5,
                          dispersal_percentage = 0.99,
                          quarantine_areas_file = "",
                          use_quarantine = FALSE,
                          use_spreadrates = FALSE,
                          use_overpopulation_movements = FALSE,
                          overpopulation_percentage = 0,
                          leaving_percentage = 0,
                          leaving_scale_coefficient = 1,
                          exposed_file_list = "",
                          mask = NULL,
                          write_outputs = "None",
                          output_folder_path = "",
                          network_filenames = c(""),
                          network_movement_types = c("walk"),
                          network_min_distances = c(0),
                          network_max_distances = c(0),
                          use_initial_condition_uncertainty = FALSE,
                          use_host_uncertainty = FALSE,
                          weather_type = "deterministic",
                          temperature_coefficient_sd_file = "",
                          precipitation_coefficient_sd_file = "",
                          dispersers_to_soils_percentage = 0,
                          quarantine_directions = "",
                          multiple_random_seeds = FALSE,
                          file_random_seeds = NULL,
                          use_soils = FALSE,
                          soil_starting_pest_file = "",
                          start_with_soil_populations = FALSE,
                          county_level_infection_data = FALSE) {
  config <- c()
  config$random_seed <- random_seed
  config$infected_file_list <- infected_file_list
  config$host_file_list <- host_file_list
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
  config$use_survival_rates <- use_survival_rates
  config$survival_rate_month <- survival_rate_month
  config$survival_rate_day <- survival_rate_day
  config$survival_rates_file <- survival_rates_file
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
  config$dispersal_stochasticity <- dispersal_stochasticity
  config$establishment_probability <- establishment_probability
  config$dispersal_percentage <- dispersal_percentage
  config$quarantine_areas_file <- quarantine_areas_file
  config$quarantine_directions <- quarantine_directions
  config$use_quarantine <- use_quarantine
  config$use_spreadrates <- use_spreadrates
  config$use_overpopulation_movements <- use_overpopulation_movements
  config$overpopulation_percentage <- overpopulation_percentage
  config$leaving_percentage <- leaving_percentage
  config$leaving_scale_coefficient <- leaving_scale_coefficient
  config$number_of_iterations <- number_of_iterations
  config$number_of_cores <- number_of_cores
  # add function name for use in configuration function to skip
  # function specific specific configurations namely for validation and
  # calibration.
  config$function_name <- "multirun"
  config$failure <- NULL
  config$exposed_file_list <- exposed_file_list
  config$mask <- mask
  config$write_outputs <- write_outputs
  config$output_folder_path <- output_folder_path
  config$mortality_frequency <- mortality_frequency
  config$mortality_frequency_n <- mortality_frequency_n
  config$network_filenames <- network_filenames
  config$network_movement_types <- network_movement_types
  config$network_min_distances <- network_min_distances
  config$network_max_distances <- network_max_distances
  config$use_initial_condition_uncertainty <- use_initial_condition_uncertainty
  config$use_host_uncertainty <- use_host_uncertainty
  config$weather_type <- weather_type
  config$temperature_coefficient_sd_file <- temperature_coefficient_sd_file
  config$precipitation_coefficient_sd_file <- precipitation_coefficient_sd_file
  config$dispersers_to_soils_percentage <- dispersers_to_soils_percentage
  config$multiple_random_seeds <- multiple_random_seeds
  config$file_random_seeds <- file_random_seeds
  config$use_soils <- use_soils
  config$soil_starting_pest_file <- soil_starting_pest_file
  config$start_with_soil_populations <- start_with_soil_populations
  config$county_level_infection_data <- county_level_infection_data
  config$pest_host_table <- pest_host_table
  config$competency_table <- competency_table

  config <- configuration(config)

  if (!is.null(config$failure)) {
    stop(config$failure)
  }

  if (config$multiple_random_seeds && is.null(config$file_random_seeds) &&
      dir.exists(config$output_folder_path)) {
    write.csv(config$random_seeds, paste0(config$output_folder_path, "forecast_random_seeds.csv"),
              row.names = FALSE)
  }

  config$crs <- terra::crs(config$host)

  i <- NULL
  cl <- parallel::makeCluster(config$core_count)
  doParallel::registerDoParallel(cl)

  infected_stack <-
    foreach::foreach(
      i = seq_len(config$number_of_iterations),
      .combine = c,
      .packages = c("PoPS", "terra")
    ) %dopar% {

      set.seed(config$random_seed[[i]])
      config <- draw_parameters(config) # draws parameter set for the run
      config <- host_pool_setup(config)
      while (any(config$total_hosts > config$total_populations, na.rm = TRUE) ||
            any(config$total_exposed > config$total_populations, na.rm = TRUE) ||
            any(config$total_infecteds > config$total_populations, na.rm = TRUE)) {
        config <- host_pool_setup(config)
      }
      config$competency_table_list <- competency_table_list_creator(config$competency_table)
      config$pest_host_table_list <- pest_host_table_list_creator(config$pest_host_table)

      data <- PoPS::pops_model(
        random_seed = config$random_seed[i],
        multiple_random_seeds = config$multiple_random_seeds,
        random_seeds = as.matrix(config$random_seeds[i, ])[1, ],
        use_lethal_temperature = config$use_lethal_temperature,
        lethal_temperature = config$lethal_temperature,
        lethal_temperature_month = config$lethal_temperature_month,
        use_survival_rates = config$use_survival_rates,
        survival_rate_month = config$survival_rate_month,
        survival_rate_day = config$survival_rate_day,
        host_pools = config$host_pools,
        total_populations = config$total_populations,
        competency_table = config$competency_table_list,
        pest_host_table = config$pest_host_table_list,
        mortality_on = config$mortality_on,
        quarantine_areas = config$quarantine_areas,
        quarantine_directions = config$quarantine_directions,
        treatment_maps = config$treatment_maps,
        treatment_dates = config$treatment_dates,
        pesticide_duration = config$pesticide_duration,
        use_movements = config$use_movements,
        movements = config$movements,
        movements_dates = config$movements_dates,
        weather = config$weather,
        temperature = config$temperature,
        survival_rates = config$survival_rates,
        weather_coefficient = config$weather_coefficient,
        weather_coefficient_sd = config$weather_coefficient_sd,
        res = config$res,
        rows_cols = config$rows_cols,
        time_step = config$time_step,
        reproductive_rate = config$reproductive_rate,
        spatial_indices = config$spatial_indices,
        season_month_start_end = config$season_month_start_end,
        soil_reservoirs = config$soil_reservoirs,
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
        network_min_distances = network_min_distances,
        network_max_distances = network_max_distances,
        network_filenames = config$network_filenames,
        network_movement_types = config$network_movement_types,
        weather_size = config$weather_size,
        weather_type = config$weather_type,
        dispersers_to_soils_percentage = config$dispersers_to_soils_percentage,
        use_soils = config$use_soils)

      outputs <- c()
      outputs$number_infected <- data$number_infected
      outputs$infected_area <- data$area_infected
      outputs$spread_rate <- data$rates
      outputs$quarantine_escape <- data$quarantine_escape
      outputs$quarantine_escape_distance <- data$quarantine_escape_distance
      outputs$quarantine_escape_direction <- data$quarantine_escape_directions
      output_host_pools <- c()

      zero_rast <- terra::rast(config$total_populations_file)[[1]]
      terra::values(zero_rast) <- 0

      total_infecteds_list <- c()
      total_infecteds <- as.matrix(zero_rast, wide = TRUE)

      for (p in seq_len(length(data$host_pools))) {
        output_host_pool <- data$host_pools[[p]]
        output_host_pool$name <- config$host_names[p]
        output_host_pools[[p]] <- output_host_pool

        config$pops_runs_folder_path <- paste(config$output_folder_path, "pops_runs/", sep = "")
        suppressWarnings(dir.create(config$pops_runs_folder_path))
        config$host_pool_folder_path <-
          paste(config$pops_runs_folder_path, config$host_names[p], sep = "")
        suppressWarnings(dir.create(config$host_pool_folder_path))

        infected_out <- zero_rast
        susectible_out <- zero_rast
        exposed_outs <- c()

        for (q in seq_len(length(data$host_pools[[p]]$infected))) {
          total_infecteds_list[[q]] <- total_infecteds
          exposed_out <- zero_rast
          if (q > 1) {
            terra::add(infected_out) <- zero_rast
            terra::add(susectible_out) <- zero_rast
          }
          total_infecteds_list[[q]] <-
            total_infecteds_list[[q]] + data$host_pools[[p]]$infected[[q]]
          terra::values(infected_out[[q]]) <- data$host_pools[[p]]$infected[[q]]
          terra::values(susectible_out[[q]]) <- data$host_pools[[p]]$susceptible[[q]]
          for (k in seq_len(length(data$host_pools[[p]]$exposed[[q]]))) {
            if (k > 1) {
              terra::add(exposed_out) <- zero_rast
            }
            terra::values(exposed_out[[k]]) <- data$host_pools[[p]]$exposed[[q]][[k]]
          }
          if (config$write_outputs == "all_simulations") {
            file_name <-
              paste(config$host_pool_folder_path, "/exposed_", i, "time_step_", q, ".tif", sep = "")
            terra::writeRaster(exposed_out, file_name, overwrite = TRUE)
          }
          exposed_outs[[q]] <- exposed_out
        }
        if (config$write_outputs == "all_simulations") {
          file_name <- paste(config$host_pool_folder_path, "/infected_", i, ".tif", sep = "")
          terra::writeRaster(infected_out, file_name, overwrite = TRUE)
          file_name <-
            paste(config$host_pool_folder_path, "/susectible_", i, ".tif", sep = "")
          terra::writeRaster(susectible_out, file_name, overwrite = TRUE)
        }
        # file_name <- paste(config$host_pool_folder_path, "/exposed_", i, ".tif", sep = "")
        # terra::writeRaster(exposed_out, file_name, overwrite = TRUE)
      }

      outputs$total_infecteds <- total_infecteds_list
      outputs$output_host_pools <- output_host_pools

      outputs
    }

  stopCluster(cl)
  number_infected_runs <- infected_stack[seq(1, length(infected_stack), 8)]
  area_infected_runs <- infected_stack[seq(2, length(infected_stack), 8)]
  spread_rate_runs <- infected_stack[seq(3, length(infected_stack), 8)]
  quarantine_escape_runs <- infected_stack[seq(4, length(infected_stack), 8)]
  quarantine_escape_distance_runs <- infected_stack[seq(5, length(infected_stack), 8)]
  quarantine_escape_directions_runs <- infected_stack[seq(6, length(infected_stack), 8)]
  total_infecteds_runs <- infected_stack[seq(7, length(infected_stack), 8)]
  output_host_pools_runs <- infected_stack[seq(8, length(infected_stack), 8)]

  prediction <- total_infecteds_runs[[1]]
  for (w in seq_len(length(prediction))) {
    prediction[[w]] <- 0
  }
  escape_probability <-
    data.frame(t(rep(0, length(total_infecteds_runs[[1]]))))
  infected_area <- data.frame(t(rep(0, length(total_infecteds_runs[[1]]))))
  infected_number <- data.frame(t(rep(0, length(total_infecteds_runs[[1]]))))
  west_rates <- data.frame(t(rep(0, length(total_infecteds_runs[[1]]))))
  east_rates <- data.frame(t(rep(0, length(total_infecteds_runs[[1]]))))
  south_rates <- data.frame(t(rep(0, length(total_infecteds_runs[[1]]))))
  north_rates <- data.frame(t(rep(0, length(total_infecteds_runs[[1]]))))
  max_values <- data.frame(t(rep(0, length(total_infecteds_runs[[1]]))))
  quarantine_escapes <- data.frame(t(rep(0, length(total_infecteds_runs[[1]]))))
  quarantine_escape_distances <- data.frame(t(rep(0, length(total_infecteds_runs[[1]]))))
  quarantine_escape_directions <- data.frame(t(rep(0, length(total_infecteds_runs[[1]]))))

  for (p in seq_len(length(total_infecteds_runs))) {
    for (w in seq_len(length(prediction))) {
      prob <- total_infecteds_runs[[p]][[w]]
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
    } else {
      west_rates[p, ] <- 0
      east_rates[p, ] <- 0
      south_rates[p, ] <- 0
      north_rates[p, ] <- 0
    }

    if (config$use_quarantine && length(quarantine_escape_runs[[p]]) ==
        length(total_infecteds_runs[[p]])) {
      escape_probability <- escape_probability + quarantine_escape_runs[[p]]
      quarantine_escapes[p, ] <- quarantine_escape_runs[[p]]
      quarantine_escape_distances[p, ] <- quarantine_escape_distance_runs[[p]]
      quarantine_escape_directions[p, ] <- quarantine_escape_directions_runs[[p]]
    }
  }

  probability <- prediction
  for (w in seq_len(length(prediction))) {
    probability[[w]] <- (prediction[[w]] / (length(total_infecteds_runs))) * 100
  }

  infected_areas <-
    round(sapply(infected_area, function(x) {
      c(
        "Mean" = mean(x, na.rm = TRUE),
        "Stand dev" = sd(x)
      )
    }), digits = 0)
  number_infecteds <-
    round(sapply(infected_number, function(x) {
      c(
        "Mean" = mean(x, na.rm = TRUE),
        "Stand dev" = sd(x)
      )
    }),
    digits = 0
    )
  west_rate <-
    round(sapply(west_rates, function(x) {
      c(
        "Mean" = mean(x, na.rm = TRUE),
        "Stand dev" = sd(x)
      )
    }), digits = 0)
  east_rate <-
    round(sapply(east_rates, function(x) {
      c(
        "Mean" = mean(x, na.rm = TRUE),
        "Stand dev" = sd(x)
      )
    }), digits = 0)
  south_rate <-
    round(sapply(south_rates, function(x) {
      c(
        "Mean" = mean(x, na.rm = TRUE),
        "Stand dev" = sd(x)
      )
    }), digits = 0)
  north_rate <-
    round(sapply(north_rates, function(x) {
      c(
        "Mean" = mean(x, na.rm = TRUE),
        "Stand dev" = sd(x)
      )
    }), digits = 0)

  west_rate[is.na(west_rate)] <- 0
  east_rate[is.na(east_rate)] <- 0
  south_rate[is.na(south_rate)] <- 0
  north_rate[is.na(north_rate)] <- 0

  if (use_quarantine) {
    escape_probability <- escape_probability / length(total_infecteds_runs) * 100
    if (
      length(quarantine_escape_distances[quarantine_escape_directions == "N"]) >
        0) {
      north_distance_to_quarantine <-
        round(sapply(
          quarantine_escape_distances[quarantine_escape_directions == "N"],
          function(x) {
            c(
              "Mean" = mean(x, na.rm = TRUE),
              "Stand dev" = sd(x)
            )
          }
        ), digits = 0)
    } else {
      north_distance_to_quarantine <- data.frame(t(rep(NA, length(total_infecteds_runs[[1]]))))
    }

    if (
      length(quarantine_escape_distances[quarantine_escape_directions == "S"]) >
        0) {
      south_distance_to_quarantine <-
        round(sapply(
          quarantine_escape_distances[quarantine_escape_directions == "S"],
          function(x) {
            c(
              "Mean" = mean(x, na.rm = TRUE),
              "Stand dev" = sd(x, na.rm = TRUE)
            )
          }
        ), digits = 0)
    } else {
      south_distance_to_quarantine <- data.frame(t(rep(NA, length(total_infecteds_runs[[1]]))))
    }

    if (
      length(quarantine_escape_distances[quarantine_escape_directions == "E"]) >
        0) {
      east_distance_to_quarantine <-
        round(sapply(
          quarantine_escape_distances[quarantine_escape_directions == "E"],
          function(x) {
            c(
              "Mean" = mean(x, na.rm = TRUE),
              "Stand dev" = sd(x)
            )
          }
        ), digits = 0)
    } else {
      east_distance_to_quarantine <- data.frame(t(rep(NA, length(total_infecteds_runs[[1]]))))
    }

    if (
      length(quarantine_escape_distances[quarantine_escape_directions == "W"]) >
        0) {
      west_distance_to_quarantine <-
        round(sapply(
          quarantine_escape_distances[quarantine_escape_directions == "W"],
          function(x) {
            c(
              "Mean" = mean(x, na.rm = TRUE),
              "Stand dev" = sd(x)
            )
          }
        ), digits = 0)
    } else {
      west_distance_to_quarantine <- data.frame(t(rep(NA, length(total_infecteds_runs[[1]]))))
    }
  } else {
    escape_probability <- data.frame(t(rep(NA, length(total_infecteds_runs[[1]]))))
    north_distance_to_quarantine <- data.frame(t(rep(NA, length(total_infecteds_runs[[1]]))))
    south_distance_to_quarantine <- data.frame(t(rep(NA, length(total_infecteds_runs[[1]]))))
    east_distance_to_quarantine <- data.frame(t(rep(NA, length(total_infecteds_runs[[1]]))))
    west_distance_to_quarantine <- data.frame(t(rep(NA, length(total_infecteds_runs[[1]]))))
  }

  which_median <- function(x) which.min(abs(x - median(x)))

  median_run_index <- which_median(infected_number[[ncol(infected_number)]])
  min_run_index <- which.min(infected_number[[ncol(infected_number)]])
  max_run_index <- which.max(infected_number[[ncol(infected_number)]])

  median_run <- total_infecteds_runs[[median_run_index]]

  min_run <- total_infecteds_runs[[min_run_index]]
  max_run <- total_infecteds_runs[[max_run_index]]

  for (q in seq_len(length(total_infecteds_runs[[1]]))) {
    for (j in seq_len(length(total_infecteds_runs))) {
      if (j == 1) {
        raster_stacks <- list(total_infecteds_runs[[j]][[q]])
      } else {
        raster_stacks[[j]] <- total_infecteds_runs[[j]][[q]]
      }
    }

    raster_stacks2 <- do.call(cbind, raster_stacks)
    raster_stacks2 <- array(raster_stacks2, dim = c(dim(raster_stacks[[1]]), length(raster_stacks)))
    sim_mean <- apply(raster_stacks2, c(1, 2), mean, na.rm = TRUE)
    sim_sd <- apply(raster_stacks2, c(1, 2), sd, na.rm = TRUE)

    simulation_mean <-
      terra::rast(
        nrow = config$rows_cols$num_rows, ncol = config$rows_cols$num_cols,
        xmin = config$xmin, xmax = config$xmax,
        ymin = config$ymin, ymax = config$ymax, crs = config$crs
      )
    simulation_sd <- simulation_mean
    simulation_max <- simulation_mean
    simulation_min <- simulation_mean
    simulation_probability <- simulation_mean
    simulation_median <- simulation_mean
    terra::values(simulation_mean) <- sim_mean
    names(simulation_mean) <- "mean"

    terra::values(simulation_sd) <- sim_sd
    names(simulation_sd) <- "sd"

    terra::values(simulation_max) <- max_run[[q]]
    names(simulation_max) <- "max"

    terra::values(simulation_min) <- min_run[[q]]
    names(simulation_min) <- "min"

    terra::values(simulation_median) <- median_run[[q]]
    names(simulation_median) <- "median"

    terra::values(simulation_probability) <- probability[[q]]
    names(simulation_probability) <- "probability"

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
    simulation_probability_stack <-
      terra::mask(simulation_probability_stack, config$mask, maskvalues = NA, updatevalue = NA)
    simulation_mean_stack <-
      terra::mask(simulation_mean_stack, config$mask, maskvalues = NA, updatevalue = NA)
    simulation_sd_stack <-
      terra::mask(simulation_sd_stack, config$mask, maskvalues = NA, updatevalue = NA)
    simulation_min_stack <-
      terra::mask(simulation_min_stack, config$mask, maskvalues = NA, updatevalue = NA)
    simulation_max_stack <-
      terra::mask(simulation_max_stack, config$mask, maskvalues = NA, updatevalue = NA)
  }

  outputs <-
    list(
      simulation_probability_stack,
      simulation_mean_stack,
      simulation_sd_stack,
      simulation_min_stack,
      simulation_max_stack,
      median_run,
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
      west_distance_to_quarantine,
      output_host_pools_runs
    )

  names(outputs) <-
    c(
      "probability",
      "simulation_mean",
      "simulation_sd",
      "simulation_min",
      "simulation_max",
      "median_run",
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
      "west_distance_to_quarantine",
      "output_host_pools_runs"
    )

  if (config$write_outputs %in% config$output_write_list) {
    file_name <- paste(config$output_folder_path, "simulation_probability.tif", sep = "")
    terra::writeRaster(simulation_probability_stack, file_name, overwrite = TRUE)
    file_name <- paste(config$output_folder_path, "simulation_mean.tif", sep = "")
    terra::writeRaster(simulation_mean_stack, file_name, overwrite = TRUE)
    file_name <- paste(config$output_folder_path, "simulation_sd.tif", sep = "")
    terra::writeRaster(simulation_sd_stack, file_name, overwrite = TRUE)
    file_name <- paste(config$output_folder_path, "simulation_min.tif", sep = "")
    terra::writeRaster(simulation_min_stack, file_name, overwrite = TRUE)
    file_name <- paste(config$output_folder_path, "simulation_max.tif", sep = "")
    terra::writeRaster(simulation_max_stack, file_name, overwrite = TRUE)
    file_name <- paste(config$output_folder_path, "multirun_outputs.rdata", sep = "")
    save(outputs, file = file_name)
  }

  return(outputs)
}
