context("test-configuration")

config <- c()
config$infected_years_file <-
  system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
config$number_of_observations <- 68
config$prior_number_of_observations <- 0
config$prior_means <- c(1, 1, 0.99, 1000, 0, 0, 0, 0)
config$prior_cov_matrix <- matrix(ncol = 8, nrow = 8, 0.1)
config$params_to_estimate <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
config$number_of_generations <- 4
config$generation_size <- 10
config$checks <- c(1200, 100000, 900, 1000)
config$infected_file_list <-
  system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
config$host_file_list <-
  system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
config$total_populations_file <-
  system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
config$parameter_means <- c(0, 1, 0.99, 1000, 0, 0, 0, 0)
config$parameter_cov_matrix <- matrix(ncol = 8, nrow = 8, 0)
config$pest_host_table <-
  system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
config$competency_table <-
  system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")
config$temp <- FALSE
config$temperature_coefficient_file <- ""
config$temperature_coefficient_sd_file <- ""
config$precip <- FALSE
config$precipitation_coefficient_file <- ""
config$precipitation_coefficient_sd_file <- ""
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
config$mortality_frequency <- "Year"
config$mortality_frequency_n <- 1
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
config$exposed_file_list <-
  system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
config$write_outputs <- "None"
config$output_folder_path <- ""
config$network_min_distance <- 0
config$network_max_distance <- 0
config$network_filename <- ""
config$network_movement <- "walk"
config$use_survival_rates <- FALSE
config$survival_rate_month <- 0
config$survival_rate_day <- 0
config$survival_rates_file <- ""
config$use_initial_condition_uncertainty <- FALSE
config$use_host_uncertainty <- FALSE
config$multiple_random_seeds <- FALSE
config$random_seeds <- NULL
config$weather_type <- "deterministic"
config$use_soils <- FALSE
config$soil_starting_pest_file <- ""

test_that("Configuration returns proper values when no errors present", {
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$start_exposed <- TRUE
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$start_exposed <- FALSE
  config$parameter_means <- c(0.2, 20, 0.99, 6000, 0, 0, 0, 0)
  config$parameter_cov_matrix <- matrix(ncol = 8, nrow = 8, 0)
  config$function_name <- "multirun"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "cauchy"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "Cauchy"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "exponential"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "Exponential"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "uniform"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "Uniform"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "deterministic neighbor"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "deterministic-neighbor"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "power law"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "power-law"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "Power-law"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "Power-Law"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "Power Law"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "Power law"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "hyperbolic secant"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "hyperbolic-secant"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "Hyperbolic-secant"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "Hyperbolic-Secant"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "Hyperbolic secant"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "Hyperbolic Secant"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "gamma"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "Gamma"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  # config$natural_kernel_type <- "exponential power"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$natural_kernel_type <- "exponential-power"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$natural_kernel_type <- "Exponential-power"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$natural_kernel_type <- "Exponential-Power"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$natural_kernel_type <- "Exponential power"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "weibull"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "Weibull"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "normal"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  # config$natural_kernel_type <- "log normal"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$natural_kernel_type <- "log-normal"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$natural_kernel_type <- "Log-normal"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$natural_kernel_type <- "Log-Normal"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$natural_kernel_type <- "Log normal"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$natural_kernel_type <- "Log Normal"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "logistic"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$natural_kernel_type <- "Logistic"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "cauchy"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "Cauchy"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "exponential"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "Exponential"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "uniform"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "Uniform"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "deterministic neighbor"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "deterministic-neighbor"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "power law"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "power-law"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "Power-law"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "Power-Law"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "Power Law"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "Power law"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "hyperbolic secant"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "hyperbolic-secant"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "Hyperbolic-secant"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "Hyperbolic-Secant"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "Hyperbolic secant"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "Hyperbolic Secant"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "gamma"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "Gamma"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  # config$anthropogenic_kernel_type <- "exponential power"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$anthropogenic_kernel_type <- "exponential-power"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$anthropogenic_kernel_type <- "Exponential-power"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$anthropogenic_kernel_type <- "Exponential-Power"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$anthropogenic_kernel_type <- "Exponential power"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "weibull"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "Weibull"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "normal"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  # config$anthropogenic_kernel_type <- "log normal"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$anthropogenic_kernel_type <- "log-normal"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$anthropogenic_kernel_type <- "Log-normal"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$anthropogenic_kernel_type <- "Log-Normal"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$anthropogenic_kernel_type <- "Log normal"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)
  #
  # config$anthropogenic_kernel_type <- "Log Normal"
  # config2 <- configuration(config)
  # expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "logistic"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$anthropogenic_kernel_type <- "Logistic"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$network_movement <- "jump"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$network_movement <- "teleport"
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$use_multiple_random_seeds <- TRUE
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$infected_files <-
    c(system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS"))
  config$parameter_means <- list(c(0.2, 20, 0.99, 6000, 0, 0, 0, 0))
  config$parameter_cov_matrix <- list(matrix(ncol = 8, nrow = 8, 0))
  config$function_name <- "auto-manage"
  config$species <- c("pest")
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)
})

test_that("configuration returns proper errors", {
  config$season_month_end <- 15
  config2 <- configuration(config)
  expect_equal(config2$failure, season_month_error)

  config$season_month_end <- 12
  config$model_type <- "SEID"
  config2 <- configuration(config)
  expect_equal(config2$failure, model_type_error)

  config$model_type <- "SEI"
  config$latency_period <- 0
  config2 <- configuration(config)
  expect_equal(config2$failure, latency_period_error)

  config$latency_period <- 1
  config$treatment_method <- "everything"
  config2 <- configuration(config)
  expect_equal(config2$failure, treatment_option_error)

  config$treatment_method <- "ratio"
  config$natural_kernel_type <- "hello"
  config2 <- configuration(config)
  expect_equal(config2$failure, natural_kernel_error)

  config$natural_kernel_type <- "cauchy"
  config$anthropogenic_kernel_type <- "good-bye"
  config2 <- configuration(config)
  expect_equal(config2$failure, anthropogenic_kernel_error)

  config$anthropogenic_kernel_type <- "cauchy"
  config$write_outputs <- "hi"
  config2 <- configuration(config)
  expect_equal(config2$failure, write_outputs_error)

  config$write_outputs <- "summary_outputs"
  config$output_folder_path <- "hi"
  config2 <- configuration(config)
  expect_equal(config2$failure, output_path_error)

  config$network_movement <- "hello"
  config$output_folder_path <- ""
  config2 <- configuration(config)
  expect_equal(config2$failure, network_movement_error)

  config$parameter_means <- c(0, 1, 0.99, 1000, 0, 0, 25, 150)
  config$function_name <- "pops"
  config$write_outputs <- "None"
  config$network_movement <- "walk"
  config$anthropogenic_kernel_type <- "network"
  config2 <- configuration(config)
  expect_equal(config2$failure, network_min_distance_small_error)

  config$parameter_means <- c(0, 1, 0.99, 1000, 0, 0, 175, 150)
  config2 <- configuration(config)
  expect_equal(config2$failure, network_min_distance_large_error)

  config$parameter_means <- c(0, 1, 0.99, 1000, 0, 0, 175, 2200)
  config2 <- configuration(config)
  expect_equal(config2$failure, network_max_distance_large_error)

  config$parameter_means <- c(0, 1, 0.99, 1000, 0, 0, 0, 0)
  config$use_initial_condition_uncertainty <- TRUE
  config2 <- configuration(config)
  expect_equal(config2$failure, initial_cond_uncert_error)

  config$use_initial_condition_uncertainty <- FALSE
  config$use_host_uncertainty <- TRUE
  config2 <- configuration(config)
  expect_equal(config2$failure, host_uncert_error)

  config$use_host_uncertainty <- FALSE
  config$weather_type <- "s"
  config2 <- configuration(config)
  expect_equal(config2$failure, weather_type_error)
})
