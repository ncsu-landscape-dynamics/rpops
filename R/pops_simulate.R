#' PoPS Simulate
#'
#' A process-based model for forecasting the spread of forest/agricultural pests
#' or pathogens. This model accounts for weather and environmental effects on
#' reproduction and survival through multiple stochastic simulations. It
#' propagates uncertainty in parameters, initial conditions, and drivers,
#' similar to pops_multirun, but directly exports raw outputs from pops_model.
#'
#' When to Use PoPS Lite:
#' - To export raw simulation data from the `pops_model` function
#' - Extent or resolution constraints make `pops_multirun` impractical due to
#' time or computational constraints.
#'
#' Multiple Random Seed handling:
#' - If use_multiple_random_seeds = TRUE and multiple_random_seeds_file = NULL in the
#' `config_rds_file`, a new set of random seeds will be sampled and exported as
#' [unique_id]_forecast_random_seeds.csv in the output folder.
#' - The same [unique_id] is used for raw output files to link seeds with
#' their respective runs.

#' @param config_rds_file Path to config file produced when calling `configuration`.
#' The config file includes all data necessary used to set up c++ PoPS model
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

pops_simulate <- function(config_rds_file = "") {

  config <- readRDS(config_rds_file)

  if (!is.null(config$failure)) {
    stop(config$failure)
  }

  # update use_multiple_random_seeds if TRUE and multiple_random_seeds_file NULL in config file
  if (config$use_multiple_random_seeds && is.null(config$multiple_random_seeds_file)) {
    config$random_seeds <- create_random_seeds(config$number_of_iterations)
  } else {
    config$random_seeds <- create_random_seeds(1)
  }

  #create output folder if it doesn't already exist
  if (!dir.exists(config$output_path)) {
    suppressWarnings(dir.create(config$output_path, recursive = TRUE))
  }

  # if use_multiple_random_seeds = TRUE and multiple_random_seeds_file = NULL
  # create unique ID to keep track of random seeds used per pops_lite simulation
  if (config$use_multiple_random_seeds &&
      is.null(config$multiple_random_seeds_file) &&
      dir.exists(config$output_path)) {
    uid <- generate_uid()
    write.csv(
      config$random_seeds,
      paste0(
        config$output_path,
        uid,
        "_forecast_random_seeds.csv"
      ),
      row.names = FALSE
    )
  }

  i <- NULL
  cl <- parallel::makeCluster(config$number_of_cores)
  doParallel::registerDoParallel(cl)

  foreach::foreach(
    i = seq_len(config$number_of_iterations),
    .packages = c("PoPS", "terra")
  ) %dopar% {

    set.seed(config$random_seed_list[[i]])
    config <- draw_parameters(config) # draws parameter set for the run
    config <- host_pool_setup(config)
    while (any(config$total_hosts > config$total_populations, na.rm = TRUE) ||
           any(config$total_exposed > config$total_populations, na.rm = TRUE) ||
           any(config$total_infecteds > config$total_populations, na.rm = TRUE)) {
      config <- host_pool_setup(config)
    }
    config$competency_table_list <- competency_table_list_creator(config$competency_table)
    config$pest_host_table_list <- pest_host_table_list_creator(config$pest_host_table)
    config$random_seed <- config$random_seed_list[[i]]
    config$random_seeds <- as.matrix(config$random_seeds_list[i, ])[1, ]

    data <- PoPS::pops_model(config)

    data[c("spatial_indices",
           "soil_reservoirs",
           "total_populations",
           "total_exposed")] <- NULL

    if (config$model_type == "SI") {
      data[c("exposed", "total_exposed", "resistant")] <- NULL
    }

    if (!config$mortality_on) {
      data$mortality <- NULL
    }

    if (!config$use_quarantine) {
      data[c(
        "quarantine_escape",
        "quarantine_escape_directions",
        "quarantine_escape_distance"
      )] <- NULL
    }
    # Remove any null in data
    data <- data[!sapply(data, is.null)]

    if (exists("uid")) {
      fn <- file.path(config$output_path, paste0(uid, "_pops_output_", i, ".rds"))
    } else {
      fn <- file.path(config$output_path, paste0("pops_output_", i, ".rds"))
    }

    saveRDS(
      data,
      file = fn,
      compress = TRUE
    )
    rm(data)
    gc()
  }
  stopCluster(cl)
  return(cat("Raw PoPS runs outputs saved to output_path: "))
}
