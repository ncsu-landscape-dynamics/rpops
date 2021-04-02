context("test-configuration")

config <- c()
config$infected_years_file <-
  system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
config$number_of_observations <- 68
config$prior_number_of_observations <- 0
config$prior_means <- c(0, 0, 0, 0, 0, 0)
config$prior_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
config$params_to_estimate <- c(T, T, T, T, F, F)
config$number_of_generations <- 4
config$generation_size <- 10
config$checks <- c(1200, 100000, 900, 1000)
config$infected_file <-
  system.file("extdata", "simple20x20", "initial_infection.tif",
              package = "PoPS")
config$host_file <-
  system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
config$total_populations_file <-
  system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
config$parameter_means <- c(0, 0, 0, 0, 0, 0)
config$parameter_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
config$temp <- FALSE
config$temperature_coefficient_file <- ""
config$precip <- FALSE
config$precipitation_coefficient_file <- ""
config$model_type <- "SEI"
config$latency_period <- 14
config$time_step <- "day"
config$season_month_start <- 1
config$season_month_end <- 12
config$start_date <- "2003-01-01"
config$end_date <- "2003-02-11"
config$use_lethal_temperature <- FALSE
config$temperature_file <- ""
config$lethal_temperature <- -30
config$lethal_temperature_month <- 1
config$mortality_on <- FALSE
config$mortality_rate <- 0
config$mortality_time_lag <- 0
config$management <- FALSE
config$treatment_dates <- c("2003-01-24")
config$treatments_file <- ""
config$treatment_method <- "ratio"
config$natural_kernel_type <- "exponential"
config$anthropogenic_kernel_type <- "cauchy"
config$natural_dir <- "NONE"
config$natural_kappa <- 0
config$anthropogenic_dir <- "NONE"
config$anthropogenic_kappa <- 0
config$pesticide_duration <- c(0)
config$pesticide_efficacy <- 1.0
config$mask <- NULL
config$success_metric <- "number of locations and total distance"
config$output_frequency <- "year"
config$output_frequency_n <- 1
config$movements_file <- ""
config$use_movements <- FALSE
config$start_exposed <- FALSE
config$generate_stochasticity <- TRUE
config$establishment_stochasticity <- TRUE
config$movement_stochasticity <- TRUE
config$deterministic <- FALSE
config$establishment_probability <- 0.5
config$dispersal_percentage <- 0.99
config$quarantine_areas_file <- ""
config$use_quarantine <- FALSE
config$use_spreadrates <- FALSE
config$percent_natural_dispersal <- 1.0
config$anthropogenic_distance_scale <- 0.0
config$output_frequency <- "week"
config$output_frequency_n <- 1
config$number_of_iterations <- 2
config$number_of_cores <- 1
config$calibration_method <- "ABC"
config$failure <- NULL
config$function_name <- "pops"
config$exposed_file <-
  system.file("extdata", "simple20x20", "initial_infection.tif",
              package = "PoPS")


test_that("Configuration returns proper values when no errors present", {
  config <- configuration(config)
  expect_equal(config$failure, NULL)

  config$start_exposed <- TRUE
  config <- configuration(config)
  expect_equal(config$failure, NULL)

  config$start_exposed <- FALSE
  config$parameter_means <- c(0.2, 20, 0.99, 6000, 0, 0)
  config$parameter_cov_matrix <- matrix(ncol = 6, nrow = 6, 0.1)
  config$function_name <- "multirun"
  config <- configuration(config)
  expect_equal(config$failure, NULL)

  config$infected_files <-
    c(system.file("extdata", "simple20x20", "initial_infection.tif",
                package = "PoPS"))
  config$parameter_means <- list(c(0.2, 20, 0.99, 6000, 0, 0))
  config$parameter_cov_matrix <- list(matrix(ncol = 6, nrow = 6, 0.1))
  config$function_name <- "auto-manage"
  config$species <- c("pest")
  config <- configuration(config)
  expect_equal(config$failure, NULL)
})

test_that("configuration returns proper errors", {
  config$model_type <- "SEID"
  config <- configuration(config)
  expect_equal(config$failure,
               "Model type is not a valid type options are 'SI' or 'SEI'")

  config$model_type <- "SEI"
  config$latency_period <- 0
  config <- configuration(config)
  expect_equal(config$failure,
               "Model type is set to SEI but the latency period is less than 1")

  config$latency_period <- 1
  config$treatment_method <- "everything"
  config <- configuration(config)
  expect_equal(config$failure,
               "treatment method is not one of the valid treatment options")
})
