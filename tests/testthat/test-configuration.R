context("test-configuration")
config <- yaml::yaml.load_file(system.file("extdata", "test_config.yaml", package = "PoPS"))

config$parameter_means <- as.vector(t(read.csv(system.file("extdata", "parameter_means_default.csv", package = "PoPS"))))
config$parameter_cov_matrix <- as.matrix(read.csv(system.file("extdata", "parameter_cov_matrix_default.csv", package = "PoPS")))
config$infected_years_file <- system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
config$infected_file_list <- system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
config$host_file_list <- system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
config$total_populations_file <- system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
config$pest_host_table <- system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
config$competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")
config$exposed_file_list <- system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")

test_that("Configuration returns proper values when no errors present", {
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$start_exposed <- TRUE
  config2 <- configuration(config)
  expect_equal(config2$failure, NULL)

  config$start_exposed <- FALSE
  config$parameter_means <- c(0.2, 20, 0.99, 6000, 0, 0)
  config$parameter_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
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

  config$function_name <- "pops"
  config$anthropogenic_kernel_type <- "network"
  config$network_min_distances <- c(50)
  config$network_max_distances <- c(150)
  config$network_movement_types <- "hello"
  config$write_outputs <- "None"
  config$output_folder_path <- ""
  config2 <- configuration(config)
  expect_equal(config2$failure, network_movement_error)

  config$network_movement_types <- "walk"
  config$network_min_distances <- c(25)
  config$network_max_distances <- c(150)
  config2 <- configuration(config)
  expect_equal(config2$failure, network_min_distance_small_error)

  config$network_min_distances <- c(175)
  config$network_max_distances <- c(150)
  config2 <- configuration(config)
  expect_equal(config2$failure, network_min_distance_large_error)
  config$network_min_distances <- c(175)
  config$network_max_distances <- c(2200)
  config2 <- configuration(config)
  expect_equal(config2$failure, network_max_distance_large_error)

  config$anthropogenic_kernel_type <- "cauchy"
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
