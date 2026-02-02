context("test-calibrate")

test_that(
  "ABC calibration has correctly formatted returns with multiple output comparisons with mask", {
    config_file <-
      system.file("extdata", "configs_calibration", "ts1_config_t1.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    saveRDS(config, file.path(system.file("extdata/", package = "PoPS"),
                              "configs_calibration/ts1_config_t1.rds"))
    config_rds_file <-
      system.file("extdata/configs_calibration/ts1_config_t1.rds", package = "PoPS")
    data <- calibrate(config_rds_file)
    expect_length(data$posterior_means, 6)
    expect_vector(data$posterior_means, ptype = double(), size = 6)
    expect_is(data$posterior_cov_matrix, class = "matrix")
    expect_type(data$posterior_cov_matrix, "double")
    expect_equal(nrow(data$posterior_cov_matrix), 6)
    expect_equal(ncol(data$posterior_cov_matrix), 6)
    expect_gt(data$posterior_means[1], 0)
    expect_gt(data$posterior_means[2], 0)
    expect_lte(data$posterior_means[2], 1000)
    expect_gte(data$posterior_means[3], 0)
    expect_lte(data$posterior_means[3], 1)
    expect_gt(data$posterior_means[4], 0)
    expect_gte(data$posterior_means[5], 0)
    expect_gte(data$posterior_means[6], 0)
    expect_type(data$total_number_of_observations, "double")
    config <- readRDS(config_rds_file)
    expect_equal(data$total_number_of_observations, config$number_of_observations)
    expect_equal(nrow(data$raw_calibration_data),
                 config$number_of_generations * config$generation_size)
  })

test_that(
  "ABC calibration has correctly formatted returns and runs with a single output comparison with
network", {
  config_file <-
    system.file("extdata", "configs_calibration", "ts2_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  saveRDS(config, file.path(system.file("extdata/", package = "PoPS"),
                              "configs_calibration/ts2_config_t1.rds"))
  config_rds_file <-
    system.file("extdata/configs_calibration/ts2_config_t1.rds", package = "PoPS")
  data <- calibrate(config_rds_file)
  expect_length(data$posterior_means, 6)
  expect_vector(data$posterior_means, ptype = double(), size = 6)
  expect_gt(data$posterior_means[1], 0)
  expect_gt(data$posterior_means[2], 0)
  expect_lte(data$posterior_means[2], 1000)
  expect_gte(data$posterior_means[3], 0)
  expect_lte(data$posterior_means[3], 1)
  expect_gt(data$posterior_means[4], 0)
  expect_gte(data$posterior_means[5], 0)
  expect_gte(data$posterior_means[6], 0)
  expect_is(data$posterior_cov_matrix, class = "matrix")
  expect_type(data$posterior_cov_matrix, "double")
  expect_equal(nrow(data$posterior_cov_matrix), 6)
  expect_equal(ncol(data$posterior_cov_matrix), 6)
  expect_type(data$total_number_of_observations, "double")
  config <- readRDS(config_rds_file)
  expect_equal(data$total_number_of_observations, config$number_of_observations)
  expect_equal(nrow(data$raw_calibration_data),
               config$number_of_generations * config$generation_size)
})

test_that(
  "ABC calibration has correctly formatted returns and runs with a single output comparison with
network", {
  config_file <-
    system.file("extdata", "configs_calibration", "ts3_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  saveRDS(config, file.path(system.file("extdata/", package = "PoPS"),
                            "configs_calibration/ts3_config_t1.rds"))
  config_rds_file <-
    system.file("extdata/configs_calibration/ts3_config_t1.rds", package = "PoPS")
  data <- calibrate(config_rds_file)
  expect_length(data$posterior_means, 6)
  expect_vector(data$posterior_means, ptype = double(), size = 6)
  expect_gt(data$posterior_means[1], 0)
  expect_gt(data$posterior_means[2], 0)
  expect_lte(data$posterior_means[2], 1000)
  expect_gte(data$posterior_means[3], 0)
  expect_lte(data$posterior_means[3], 1)
  expect_gt(data$posterior_means[4], 0)
  expect_gte(data$posterior_means[5], 0)
  expect_gte(data$posterior_means[6], 0)
  expect_is(data$posterior_cov_matrix, class = "matrix")
  expect_type(data$posterior_cov_matrix, "double")
  expect_equal(nrow(data$posterior_cov_matrix), 6)
  expect_equal(ncol(data$posterior_cov_matrix), 6)
  expect_type(data$total_number_of_observations, "double")
  config <- readRDS(config_rds_file)
  expect_equal(data$total_number_of_observations, config$number_of_observations)
  expect_equal(nrow(data$raw_calibration_data),
               config$number_of_generations * config$generation_size)
})

test_that(
  "ABC calibration has correctly formatted returns and runs with a single output comparison with
network", {
  config_file <-
    system.file("extdata", "configs_calibration", "ts4_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  saveRDS(config, file.path(system.file("extdata/", package = "PoPS"),
                            "configs_calibration/ts4_config_t1.rds"))
  config_rds_file <-
    system.file("extdata/configs_calibration/ts4_config_t1.rds", package = "PoPS")
  data <- calibrate(config_rds_file)
  expect_length(data$posterior_means, 6)
  expect_vector(data$posterior_means, ptype = double(), size = 6)
  expect_gt(data$posterior_means[1], 0)
  expect_gt(data$posterior_means[2], 0)
  expect_lte(data$posterior_means[2], 1000)
  expect_gte(data$posterior_means[3], 0)
  expect_lte(data$posterior_means[3], 1)
  expect_gt(data$posterior_means[4], 0)
  expect_gte(data$posterior_means[5], 0)
  expect_gte(data$posterior_means[6], 0)
  expect_is(data$posterior_cov_matrix, class = "matrix")
  expect_type(data$posterior_cov_matrix, "double")
  expect_equal(nrow(data$posterior_cov_matrix), 6)
  expect_equal(ncol(data$posterior_cov_matrix), 6)
  expect_type(data$total_number_of_observations, "double")
  config <- readRDS(config_rds_file)
  expect_equal(data$total_number_of_observations, config$number_of_observations)
  expect_equal(nrow(data$raw_calibration_data),
               config$number_of_generations * config$generation_size)
})

test_that(
  "ABC calibration has correctly formatted returns/runs with host and initial condition
uncertainty", {
  config_file <-
    system.file("extdata", "configs_calibration", "ts5_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  saveRDS(config, file.path(system.file("extdata/", package = "PoPS"),
                            "configs_calibration/ts5_config_t1.rds"))
  config_rds_file <-
    system.file("extdata/configs_calibration/ts5_config_t1.rds", package = "PoPS")
  data <- calibrate(config_rds_file)
  expect_length(data$posterior_means, 6)
  expect_vector(data$posterior_means, ptype = double(), size = 6)
  expect_gt(data$posterior_means[1], 0)
  expect_gt(data$posterior_means[2], 0)
  expect_lte(data$posterior_means[2], 1000)
  expect_gte(data$posterior_means[3], 0)
  expect_lte(data$posterior_means[3], 1)
  expect_gt(data$posterior_means[4], 0)
  expect_gte(data$posterior_means[5], 0)
  expect_gte(data$posterior_means[6], 0)
  expect_is(data$posterior_cov_matrix, class = "matrix")
  expect_type(data$posterior_cov_matrix, "double")
  expect_equal(nrow(data$posterior_cov_matrix), 6)
  expect_equal(ncol(data$posterior_cov_matrix), 6)
  expect_type(data$total_number_of_observations, "double")
  config <- readRDS(config_rds_file)
  expect_equal(data$total_number_of_observations, config$number_of_observations)
  expect_equal(nrow(data$raw_calibration_data),
               config$number_of_generations * config$generation_size)
})

test_that("ABC calibration has correctly formatted returns/runs with county level data", {
  config_file <-
    system.file("extdata", "configs_calibration", "ts6_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  saveRDS(config, file.path(system.file("extdata/", package = "PoPS"),
                            "configs_calibration/ts6_config_t1.rds"))
  config_rds_file <-
    system.file("extdata/configs_calibration/ts6_config_t1.rds", package = "PoPS")
  data <- calibrate(config_rds_file)
  expect_length(data$posterior_means, 6)
  expect_vector(data$posterior_means, ptype = double(), size = 6)
  expect_gt(data$posterior_means[1], 0)
  expect_gt(data$posterior_means[2], 0)
  expect_lte(data$posterior_means[2], 1000)
  expect_gte(data$posterior_means[3], 0)
  expect_lte(data$posterior_means[3], 1)
  expect_gt(data$posterior_means[4], 0)
  expect_gte(data$posterior_means[5], 0)
  expect_gte(data$posterior_means[6], 0)
  expect_is(data$posterior_cov_matrix, class = "matrix")
  expect_type(data$posterior_cov_matrix, "double")
  expect_equal(nrow(data$posterior_cov_matrix), 6)
  expect_equal(ncol(data$posterior_cov_matrix), 6)
  expect_type(data$total_number_of_observations, "double")
  config <- readRDS(config_rds_file)
  expect_equal(data$total_number_of_observations, config$number_of_observations)
  expect_equal(nrow(data$raw_calibration_data),
               config$number_of_generations * config$generation_size)
})

test_that("MCMC calibration has correctly formatted returns with multiple output comparisons", {
  config_file <-
    system.file("extdata", "configs_calibration", "ts7_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  saveRDS(config, file.path(system.file("extdata/", package = "PoPS"),
                            "configs_calibration/ts7_config_t1.rds"))
  config_rds_file <-
    system.file("extdata/configs_calibration/ts7_config_t1.rds", package = "PoPS")
  data <- calibrate(config_rds_file)
  expect_length(data$posterior_means, 6)
  expect_vector(data$posterior_means, ptype = double(), size = 6)
  expect_gt(data$posterior_means[1], 0)
  expect_gt(data$posterior_means[2], 0)
  expect_lte(data$posterior_means[2], 1000)
  expect_gte(data$posterior_means[3], 0)
  expect_lte(data$posterior_means[3], 1)
  expect_gt(data$posterior_means[4], 0)
  expect_gte(data$posterior_means[5], 0)
  expect_gte(data$posterior_means[6], 0)
  expect_is(data$posterior_cov_matrix, class = "matrix")
  expect_type(data$posterior_cov_matrix, "double")
  expect_equal(nrow(data$posterior_cov_matrix), 6)
  expect_equal(ncol(data$posterior_cov_matrix), 6)
  expect_type(data$total_number_of_observations, "double")
  config <- readRDS(config_rds_file)
  expect_equal(data$total_number_of_observations, config$number_of_observations)
  expect_equal(nrow(data$raw_calibration_data),
               config$number_of_generations * config$generation_size)
})

test_that(
  "MCMC calibration has correctly formatted returns with multiple output comparisons with mask", {
    config_file <-
      system.file("extdata", "configs_calibration", "ts8_config_t1.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    saveRDS(config, file.path(system.file("extdata/", package = "PoPS"),
                              "configs_calibration/ts8_config_t1.rds"))
    config_rds_file <-
      system.file("extdata/configs_calibration/ts8_config_t1.rds", package = "PoPS")
    data <- calibrate(config_rds_file)
    expect_length(data$posterior_means, 6)
    expect_vector(data$posterior_means, ptype = double(), size = 6)
    expect_gt(data$posterior_means[1], 0)
    expect_gt(data$posterior_means[2], 0)
    expect_lte(data$posterior_means[2], 1000)
    expect_gte(data$posterior_means[3], 0)
    expect_lte(data$posterior_means[3], 1)
    expect_gt(data$posterior_means[4], 0)
    expect_gte(data$posterior_means[5], 0)
    expect_gte(data$posterior_means[6], 0)
    expect_is(data$posterior_cov_matrix, class = "matrix")
    expect_type(data$posterior_cov_matrix, "double")
    expect_equal(nrow(data$posterior_cov_matrix), 6)
    expect_equal(ncol(data$posterior_cov_matrix), 6)
    expect_type(data$total_number_of_observations, "double")
    config <- readRDS(config_rds_file)
    expect_equal(data$total_number_of_observations, config$number_of_observations)
    expect_equal(nrow(data$raw_calibration_data),
                 config$number_of_generations * config$generation_size)
  })

test_that(
  "MCMC calibration has correctly formatted returns with multiple output comparisons with mask", {
    config_file <-
      system.file("extdata", "configs_calibration", "ts9_config_t1.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    saveRDS(config, file.path(system.file("extdata/", package = "PoPS"),
                              "configs_calibration/ts9_config_t1.rds"))
    config_rds_file <-
      system.file("extdata/configs_calibration/ts9_config_t1.rds", package = "PoPS")
    data <- calibrate(config_rds_file)
    expect_length(data$posterior_means, 6)
    expect_vector(data$posterior_means, ptype = double(), size = 6)
    expect_gt(data$posterior_means[1], 0)
    expect_gt(data$posterior_means[2], 0)
    expect_lte(data$posterior_means[2], 1000)
    expect_gte(data$posterior_means[3], 0)
    expect_lte(data$posterior_means[3], 1)
    expect_gt(data$posterior_means[4], 0)
    expect_gte(data$posterior_means[5], 0)
    expect_gte(data$posterior_means[6], 0)
    expect_is(data$posterior_cov_matrix, class = "matrix")
    expect_type(data$posterior_cov_matrix, "double")
    expect_equal(nrow(data$posterior_cov_matrix), 6)
    expect_equal(ncol(data$posterior_cov_matrix), 6)
    expect_type(data$total_number_of_observations, "double")
    config <- readRDS(config_rds_file)
    expect_equal(data$total_number_of_observations, config$number_of_observations)
    expect_equal(nrow(data$raw_calibration_data),
                 config$number_of_generations * config$generation_size)
  })

test_that(
  "MCMC calibration has correctly formatted returns with multiple output comparisons with mask", {
    config_file <-
      system.file("extdata", "configs_calibration", "ts10_config_t1.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    saveRDS(config, file.path(system.file("extdata/", package = "PoPS"),
                              "configs_calibration/ts10_config_t1.rds"))
    config_rds_file <-
      system.file("extdata/configs_calibration/ts10_config_t1.rds", package = "PoPS")
    data <- calibrate(config_rds_file)
    expect_length(data$posterior_means, 6)
    expect_vector(data$posterior_means, ptype = double(), size = 6)
    expect_gt(data$posterior_means[1], 0)
    expect_gt(data$posterior_means[2], 0)
    expect_lte(data$posterior_means[2], 1000)
    expect_gte(data$posterior_means[3], 0)
    expect_lte(data$posterior_means[3], 1)
    expect_gt(data$posterior_means[4], 0)
    expect_gte(data$posterior_means[5], 0)
    expect_gte(data$posterior_means[6], 0)
    expect_is(data$posterior_cov_matrix, class = "matrix")
    expect_type(data$posterior_cov_matrix, "double")
    expect_equal(nrow(data$posterior_cov_matrix), 6)
    expect_equal(ncol(data$posterior_cov_matrix), 6)
    expect_type(data$total_number_of_observations, "double")
    config <- readRDS(config_rds_file)
    expect_equal(data$total_number_of_observations, config$number_of_observations)
    expect_equal(nrow(data$raw_calibration_data),
                 config$number_of_generations * config$generation_size)
  })

test_that(
  "MCMC calibration has correctly formatted returns with host and initial condition uncertainty", {
    config_file <-
      system.file("extdata", "configs_calibration", "ts11_config_t1.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    saveRDS(config, file.path(system.file("extdata/", package = "PoPS"),
                              "configs_calibration/ts11_config_t1.rds"))
    config_rds_file <-
      system.file("extdata/configs_calibration/ts11_config_t1.rds", package = "PoPS")
    data <- calibrate(config_rds_file)
    expect_length(data$posterior_means, 6)
    expect_vector(data$posterior_means, ptype = double(), size = 6)
    expect_gt(data$posterior_means[1], 0)
    expect_gt(data$posterior_means[2], 0)
    expect_lte(data$posterior_means[2], 1000)
    expect_gte(data$posterior_means[3], 0)
    expect_lte(data$posterior_means[3], 1)
    expect_gt(data$posterior_means[4], 0)
    expect_gte(data$posterior_means[5], 0)
    expect_gte(data$posterior_means[6], 0)
    expect_is(data$posterior_cov_matrix, class = "matrix")
    expect_type(data$posterior_cov_matrix, "double")
    expect_equal(nrow(data$posterior_cov_matrix), 6)
    expect_equal(ncol(data$posterior_cov_matrix), 6)
    expect_type(data$total_number_of_observations, "double")
    config <- readRDS(config_rds_file)
    expect_equal(data$total_number_of_observations, config$number_of_observations)
    expect_equal(nrow(data$raw_calibration_data),
                 config$number_of_generations * config$generation_size)
  })

test_that(
  "MCMC calibration has correctly formatted returns with host and initial condition uncertainty", {
    config_file <-
      system.file("extdata", "configs_calibration", "ts12_config_t1.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    saveRDS(config, file.path(system.file("extdata/", package = "PoPS"),
                              "configs_calibration/ts12_config_t1.rds"))
    config_rds_file <-
      system.file("extdata/configs_calibration/ts12_config_t1.rds", package = "PoPS")
    data <- calibrate(config_rds_file)
    expect_length(data$posterior_means, 6)
    expect_vector(data$posterior_means, ptype = double(), size = 6)
    expect_gt(data$posterior_means[1], 0)
    expect_gt(data$posterior_means[2], 0)
    expect_lte(data$posterior_means[2], 1000)
    expect_gte(data$posterior_means[3], 0)
    expect_lte(data$posterior_means[3], 1)
    expect_gt(data$posterior_means[4], 0)
    expect_gte(data$posterior_means[5], 0)
    expect_gte(data$posterior_means[6], 0)
    expect_is(data$posterior_cov_matrix, class = "matrix")
    expect_type(data$posterior_cov_matrix, "double")
    expect_equal(nrow(data$posterior_cov_matrix), 6)
    expect_equal(ncol(data$posterior_cov_matrix), 6)
    expect_type(data$total_number_of_observations, "double")
    config <- readRDS(config_rds_file)
    expect_equal(data$total_number_of_observations, config$number_of_observations)
    expect_equal(nrow(data$raw_calibration_data),
                 config$number_of_generations * config$generation_size)
  })


test_that(
  "ABC calibration has correctly formatted returns with multiple output comparisons with mask", {
  config_file <-
    system.file("extdata", "configs_calibration", "ts13_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  saveRDS(config, file.path(system.file("extdata/", package = "PoPS"),
                            "configs_calibration/ts13_config_t1.rds"))
  config_rds_file <-
    system.file("extdata/configs_calibration/ts13_config_t1.rds", package = "PoPS")
  data <- calibrate(config_rds_file)
  expect_length(data$posterior_means, 6)
  expect_vector(data$posterior_means, ptype = double(), size = 6)
  expect_gt(data$posterior_means[1], 0)
  expect_gt(data$posterior_means[2], 0)
  expect_lte(data$posterior_means[2], 1000)
  expect_gte(data$posterior_means[3], 0)
  expect_lte(data$posterior_means[3], 1)
  expect_gt(data$posterior_means[4], 0)
  expect_gte(data$posterior_means[5], 0)
  expect_gte(data$posterior_means[6], 0)
  expect_is(data$posterior_cov_matrix, class = "matrix")
  expect_type(data$posterior_cov_matrix, "double")
  expect_equal(nrow(data$posterior_cov_matrix), 6)
  expect_equal(ncol(data$posterior_cov_matrix), 6)
  expect_type(data$total_number_of_observations, "double")
  config <- readRDS(config_rds_file)
  expect_equal(data$total_number_of_observations, config$number_of_observations)
  expect_equal(nrow(data$raw_calibration_data),
               config$number_of_generations * config$generation_size)
})
