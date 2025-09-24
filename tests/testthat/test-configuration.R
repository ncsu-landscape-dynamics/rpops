context("test-configuration")

test_that("Configuration returns proper values when no errors present", {
  config_file <- system.file("extdata", "test_config.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(config$failure, NULL)
})

test_that("configuration returns proper errors", {
  config_file <- system.file("extdata", "test_config_season_month_error.yml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(config$failure, season_month_error)

  config_file <- system.file("extdata", "test_config_model_type_error.yml", package = "PoPS")
  config$model_type <- "SEID"
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(config$failure, model_type_error)

  config_file <- system.file("extdata", "test_config_latency_period_error.yml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(config$failure, latency_period_error)

  config_file <- system.file("extdata", "test_config_treatment_option_error.yml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(config$failure, treatment_option_error)

  config_file <- system.file("extdata", "test_config_natural_kernel_error.yml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(config$failure, natural_kernel_error)

  config_file <- system.file("extdata", "test_config_anthropogenic_kernel_error.yml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(config$failure, anthropogenic_kernel_error)

  config_file <- system.file("extdata", "test_config_output_path_error.yml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(config$failure, output_path_error)

  config_file <- system.file("extdata", "test_config_network_movement_error.yml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(config$failure, network_movement_error)

  config_file <- system.file("extdata", "test_config_network_min_distance_small_error.yml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(config$failure, network_min_distance_small_error)

  config_file <- system.file("extdata", "test_config_network_min_distance_large_error.yml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(config$failure, network_min_distance_large_error)

  config_file <- system.file("extdata", "test_config_network_max_distance_large_error.yml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(config$failure, network_max_distance_large_error)

  config_file <- system.file("extdata", "test_config_initial_cond_uncert_error.yml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(config$failure, initial_cond_uncert_error)

  config_file <- system.file("extdata", "test_config_host_uncert_error.yml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(config$failure, host_uncert_error)

  config_file <- system.file("extdata", "test_config_weather_type_error.yml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(config$failure, weather_type_error)
})
