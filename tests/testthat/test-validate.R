context("test-validate-simulate")

test_that("TS1: Validation has correctly formatted returns with multiple output comparisons", {
    config_file <-
      system.file("extdata", "configs_simulate_validate", "ts1_config_t1.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    saveRDS(config,
            file.path(system.file("extdata/", package = "PoPS"),
                      "configs_simulate_validate/ts1_config_t1.rds"))
    config_rds_file <-
      system.file("extdata/configs_simulate_validate/ts1_config_t1.rds", package = "PoPS")
    pops_simulate(config_rds_file = config_rds_file)
    create_summary_stats_and_stacks(config_rds_file = config_rds_file)
    outputs <- validate(config_rds_file = config_rds_file)

    number_of_timesteps <- 5
    expect_type(outputs, "list")
    expect_length(outputs, 26)
    expect_vector(outputs$quantity_disagreement, size = number_of_timesteps)
    expect_vector(outputs$allocation_disagreement, size = number_of_timesteps)
    expect_vector(outputs$total_disagreement, size = number_of_timesteps)
    expect_vector(outputs$configuration_disagreement, size = number_of_timesteps)
    expect_vector(outputs$false_negatives, size = number_of_timesteps)
    expect_vector(outputs$false_positives, size = number_of_timesteps)
    expect_vector(outputs$true_positives, size = number_of_timesteps)
    expect_vector(outputs$true_negatives, size = number_of_timesteps)
    expect_vector(outputs$unknown_positives, size = number_of_timesteps)
    expect_vector(outputs$unknown_negatives, size = number_of_timesteps)
    expect_vector(outputs$odds_ratio, size = number_of_timesteps)
    expect_vector(outputs$residual_error, size = number_of_timesteps)
    expect_vector(outputs$true_infecteds, size = number_of_timesteps)
    expect_vector(outputs$simulated_infecteds, size = number_of_timesteps)
    expect_vector(outputs$infecteds_difference, size = number_of_timesteps)
  }
)

test_that("Validation has correctly formatted returns and runs with a single output comparison", {
    config_file <-
      system.file("extdata", "configs_simulate_validate", "ts2_config_t1.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    saveRDS(config,
            file.path(system.file("extdata/", package = "PoPS"),
                      "configs_simulate_validate/ts2_config_t1.rds"))
    config_rds_file <-
      system.file("extdata/configs_simulate_validate/ts2_config_t1.rds", package = "PoPS")
    pops_simulate(config_rds_file = config_rds_file)
    create_summary_stats_and_stacks(config_rds_file = config_rds_file)
    outputs <- validate(config_rds_file = config_rds_file)
    number_of_timesteps <- 1
    expect_type(outputs, "list")
    expect_length(outputs, 26)
    expect_vector(outputs$quantity_disagreement, size = number_of_timesteps)
    expect_vector(outputs$allocation_disagreement, size = number_of_timesteps)
    expect_vector(outputs$total_disagreement, size = number_of_timesteps)
    expect_vector(outputs$configuration_disagreement, size = number_of_timesteps)
    expect_vector(outputs$false_negatives, size = number_of_timesteps)
    expect_vector(outputs$false_positives, size = number_of_timesteps)
    expect_vector(outputs$true_positives, size = number_of_timesteps)
    expect_vector(outputs$true_negatives, size = number_of_timesteps)
    expect_vector(outputs$unknown_positives, size = number_of_timesteps)
    expect_vector(outputs$unknown_negatives, size = number_of_timesteps)
    expect_vector(outputs$odds_ratio, size = number_of_timesteps)
    expect_vector(outputs$residual_error, size = number_of_timesteps)
    expect_vector(outputs$true_infecteds, size = number_of_timesteps)
    expect_vector(outputs$simulated_infecteds, size = number_of_timesteps)
    expect_vector(outputs$infecteds_difference, size = number_of_timesteps)
  }
)

test_that(
  "Validation has correctly formatted returns and runs with a single output comparison with mask", {
    config_file <-
      system.file("extdata", "configs_simulate_validate", "ts3_config_t1.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    saveRDS(config,
            file.path(system.file("extdata/", package = "PoPS"),
                      "configs_simulate_validate/ts3_config_t1.rds"))
    config_rds_file <-
      system.file("extdata/configs_simulate_validate/ts3_config_t1.rds", package = "PoPS")
    pops_simulate(config_rds_file = config_rds_file)
    create_summary_stats_and_stacks(config_rds_file = config_rds_file)
    outputs <- validate(config_rds_file = config_rds_file)
    number_of_timesteps <- 1
    expect_type(outputs, "list")
    expect_length(outputs, 26)
    expect_vector(outputs$quantity_disagreement, size = number_of_timesteps)
    expect_vector(outputs$allocation_disagreement, size = number_of_timesteps)
    expect_vector(outputs$total_disagreement, size = number_of_timesteps)
    expect_vector(outputs$configuration_disagreement, size = number_of_timesteps)
    expect_vector(outputs$false_negatives, size = number_of_timesteps)
    expect_vector(outputs$false_positives, size = number_of_timesteps)
    expect_vector(outputs$true_positives, size = number_of_timesteps)
    expect_vector(outputs$true_negatives, size = number_of_timesteps)
    expect_vector(outputs$unknown_positives, size = number_of_timesteps)
    expect_vector(outputs$unknown_negatives, size = number_of_timesteps)
    expect_vector(outputs$odds_ratio, size = number_of_timesteps)
    expect_vector(outputs$residual_error, size = number_of_timesteps)
    expect_vector(outputs$true_infecteds, size = number_of_timesteps)
    expect_vector(outputs$simulated_infecteds, size = number_of_timesteps)
    expect_vector(outputs$infecteds_difference, size = number_of_timesteps)
  }
)

test_that(
  "Validation has correctly formatted returns/runs with host and initial condition uncertainty", {
    config_file <-
      system.file("extdata", "configs_simulate_validate", "ts4_config_t1.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    saveRDS(config,
            file.path(system.file("extdata/", package = "PoPS"),
                      "configs_simulate_validate/ts4_config_t1.rds"))
    config_rds_file <-
      system.file("extdata/configs_simulate_validate/ts4_config_t1.rds", package = "PoPS")
    pops_simulate(config_rds_file = config_rds_file)
    create_summary_stats_and_stacks(config_rds_file = config_rds_file)
    outputs <- validate(config_rds_file = config_rds_file)
    number_of_timesteps <- 1
    expect_type(outputs, "list")
    expect_length(outputs, 26)
    expect_vector(outputs$quantity_disagreement, size = number_of_timesteps)
    expect_vector(outputs$allocation_disagreement, size = number_of_timesteps)
    expect_vector(outputs$total_disagreement, size = number_of_timesteps)
    expect_vector(outputs$configuration_disagreement, size = number_of_timesteps)
    expect_vector(outputs$false_negatives, size = number_of_timesteps)
    expect_vector(outputs$false_positives, size = number_of_timesteps)
    expect_vector(outputs$true_positives, size = number_of_timesteps)
    expect_vector(outputs$true_negatives, size = number_of_timesteps)
    expect_vector(outputs$unknown_positives, size = number_of_timesteps)
    expect_vector(outputs$unknown_negatives, size = number_of_timesteps)
    expect_vector(outputs$odds_ratio, size = number_of_timesteps)
    expect_vector(outputs$residual_error, size = number_of_timesteps)
    expect_vector(outputs$true_infecteds, size = number_of_timesteps)
    expect_vector(outputs$simulated_infecteds, size = number_of_timesteps)
    expect_vector(outputs$infecteds_difference, size = number_of_timesteps)
  }
)

test_that("Validation has correctly formatted returns/runs with polygon level data", {
  config_file <-
    system.file("extdata", "configs_simulate_validate", "ts5_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  saveRDS(config,
          file.path(system.file("extdata/", package = "PoPS"),
                    "configs_simulate_validate/ts5_config_t1.rds"))
  config_rds_file <-
    system.file("extdata/configs_simulate_validate/ts5_config_t1.rds", package = "PoPS")
  pops_simulate(config_rds_file = config_rds_file)
  create_summary_stats_and_stacks(config_rds_file = config_rds_file)
  outputs <- validate(config_rds_file = config_rds_file)
  number_of_timesteps <- 1
  expect_type(outputs, "list")
  expect_length(outputs, 28)
  expect_vector(outputs$false_negatives, size = number_of_timesteps)
  expect_vector(outputs$false_positives, size = number_of_timesteps)
  expect_vector(outputs$true_positives, size = number_of_timesteps)
  expect_vector(outputs$true_negatives, size = number_of_timesteps)
  expect_vector(outputs$unknown_positives, size = number_of_timesteps)
  expect_vector(outputs$unknown_negatives, size = number_of_timesteps)
  expect_vector(outputs$odds_ratio, size = number_of_timesteps)
  expect_vector(outputs$residual_error, size = number_of_timesteps)
  }
)
