context("test-abc_calibration")

test_that("ABC calibration has correctly formatted returns with multiple output comparisons", {
  skip_on_appveyor()
  infected_years_file <- system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  number_of_observations <- 68    ### This is the number of infected cells - just make sure it's consistent across years 
  prior_number_of_observations <- 0
  prior_means <- c(0, 0, 0, 0, 0, 0)    ### leave as 0 for now, means that you are giving them no weight
  prior_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
  params_to_estimate <- c(T, T, T, T, F, F)  ### 1st: reproductive rates, 2nd: natural distance, 3rd: percent natural, 4tH: anthropogenic distance, 5th Natural Kappa, 6th anthropogenic kappa
  number_of_generations <- 6
  generation_size <- 10
  checks = c(1200, 100000, 900, 1000)
  infected_file <- system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
  host_file <- system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
  total_plants_file <- system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
  temp <- FALSE
  temperature_coefficient_file <- ""
  precip <- FALSE
  precipitation_coefficient_file <- ""
  model_type = "SEI"
  latency_period = 14
  time_step <- "day"
  season_month_start <- 1
  season_month_end <- 12
  start_date <- '2003-01-01'
  end_date <- '2003-02-11'
  use_lethal_temperature <- FALSE
  temperature_file <- ""
  lethal_temperature <- -30
  lethal_temperature_month <- 1
  mortality_on <- FALSE
  mortality_rate <- 0
  mortality_time_lag <- 0
  management <- FALSE
  treatment_dates <- c('2003-01-24')
  treatments_file <- ""
  treatment_method <- "ratio"
  natural_kernel_type <- "exponential"
  anthropogenic_kernel_type <- "cauchy"
  natural_dir <- "NONE"
  natural_kappa <- 0
  anthropogenic_dir <- "NONE"
  anthropogenic_kappa <- 0
  pesticide_duration <- c(0)
  pesticide_efficacy <- 1.0
  mask <- NULL
  success_metric <- "number of locations and total distance"
  output_frequency <- "week"
  movements_file = ""
  use_movements = FALSE
  percent_natural_dispersal <- 1.0
  anthropogenic_distance_scale <- 0.0
  
  data <- abc_calibration(infected_years_file, 
                          number_of_observations, 
                          prior_number_of_observations,
                          prior_means, prior_cov_matrix, 
                          params_to_estimate,
                          number_of_generations,
                          generation_size,
                          checks,
                          infected_file, 
                          host_file, 
                          total_plants_file, 
                          temp, 
                          temperature_coefficient_file, 
                          precip, 
                          precipitation_coefficient_file,
                          model_type,
                          latency_period,
                          time_step, 
                          season_month_start, 
                          season_month_end, 
                          start_date, 
                          end_date, 
                          use_lethal_temperature, 
                          temperature_file,
                          lethal_temperature, 
                          lethal_temperature_month,
                          mortality_on, 
                          mortality_rate, 
                          mortality_time_lag, 
                          management, 
                          treatment_dates, 
                          treatments_file,
                          treatment_method,
                          natural_kernel_type, 
                          anthropogenic_kernel_type,
                          natural_dir, 
                          natural_kappa, 
                          anthropogenic_dir, 
                          anthropogenic_kappa,
                          pesticide_duration, 
                          pesticide_efficacy,
                          mask, 
                          success_metric, 
                          output_frequency,
                          movements_file, 
                          use_movements)
  
  expect_length(data$posterior_means, 6)
  expect_vector(data$posterior_means, ptype = double(), size = 6)
  expect_is(data$posterior_cov_matrix, class = "matrix")
  expect_type(data$posterior_cov_matrix, "double")
  expect_equal(nrow(data$posterior_cov_matrix), 6)
  expect_equal(ncol(data$posterior_cov_matrix), 6)
  expect_type(data$total_number_of_observations, "double")
  expect_equal(data$total_number_of_observations, number_of_observations)
  expect_equal(nrow(data$raw_calibration_data), number_of_generations*generation_size)
  
})


# test_that("ABC calibration has correctly formatted returns and runs with a single output comparison", {
#   infected_years_file <- system.file("extdata", "simple20x20", "infected_single.tif", package = "PoPS")
#   number_of_observations <- 68    ### This is the number of infected cells - just make sure it's consistent across years 
#   prior_number_of_observations <- 0
#   prior_means <- c(0, 0, 0, 0, 0, 0)    ### leave as 0 for now, means that you are giving them no weight
#   prior_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
#   params_to_estimate <- c(T, T, T, T, F, F)  ### 1st: reproductive rates, 2nd: natural distance, 3rd: percent natural, 4tH: anthropogenic distance, 5th Natural Kappa, 6th anthropogenic kappa
#   number_of_generations <- 6
#   generation_size <- 10
#   checks = c(1200, 90000, 900, 1000)
#   infected_file <- system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
#   host_file <- system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
#   total_plants_file <- system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
#   temp <- FALSE
#   temperature_coefficient_file <- ""
#   precip <- FALSE
#   precipitation_coefficient_file <- ""
#   model_type = "SI"
#   latency_period = 0
#   time_step <- "month"
#   season_month_start <- 1
#   season_month_end <- 12
#   start_date <- '2003-01-01'
#   end_date <- '2003-12-31'
#   use_lethal_temperature <- FALSE
#   temperature_file <- ""
#   lethal_temperature <- -30
#   lethal_temperature_month <- 1
#   mortality_on <- FALSE
#   mortality_rate <- 0
#   mortality_time_lag <- 0
#   management <- FALSE
#   treatment_dates <- c('2003-01-24')
#   treatments_file <- ""
#   treatment_method <- "ratio"
#   natural_kernel_type <- "exponential"
#   anthropogenic_kernel_type <- "cauchy"
#   natural_dir <- "NONE"
#   natural_kappa <- 0
#   anthropogenic_dir <- "NONE"
#   anthropogenic_kappa <- 0
#   pesticide_duration <- c(0)
#   pesticide_efficacy <- 1.0
#   mask <- NULL
#   success_metric <- "number of locations and total distance"
#   output_frequency <- "year"
#   movements_file = ""
#   use_movements = FALSE
#   percent_natural_dispersal <- 1.0
#   anthropogenic_distance_scale <- 0.0
#   
#   data <- abc_calibration(infected_years_file, 
#                           number_of_observations, 
#                           prior_number_of_observations,
#                           prior_means, prior_cov_matrix, 
#                           params_to_estimate,
#                           number_of_generations,
#                           generation_size,
#                           checks,
#                           infected_file, 
#                           host_file, 
#                           total_plants_file, 
#                           temp, 
#                           temperature_coefficient_file, 
#                           precip, 
#                           precipitation_coefficient_file,
#                           model_type,
#                           latency_period,
#                           time_step, 
#                           season_month_start, 
#                           season_month_end, 
#                           start_date, 
#                           end_date, 
#                           use_lethal_temperature, 
#                           temperature_file,
#                           lethal_temperature, 
#                           lethal_temperature_month,
#                           mortality_on, 
#                           mortality_rate, 
#                           mortality_time_lag, 
#                           management, 
#                           treatment_dates, 
#                           treatments_file,
#                           treatment_method,
#                           natural_kernel_type, 
#                           anthropogenic_kernel_type,
#                           natural_dir, 
#                           natural_kappa, 
#                           anthropogenic_dir, 
#                           anthropogenic_kappa,
#                           pesticide_duration, 
#                           pesticide_efficacy,
#                           mask, 
#                           success_metric, 
#                           output_frequency,
#                           movements_file, 
#                           use_movements)
#   
#   expect_length(data$posterior_means, 6)
#   expect_vector(data$posterior_means, ptype = double(), size = 6)
#   expect_is(data$posterior_cov_matrix, class = "matrix")
#   expect_type(data$posterior_cov_matrix, "double")
#   expect_equal(nrow(data$posterior_cov_matrix), 6)
#   expect_equal(ncol(data$posterior_cov_matrix), 6)
#   expect_type(data$total_number_of_observations, "double")
#   expect_equal(data$total_number_of_observations, number_of_observations)
#   expect_equal(nrow(data$raw_calibration_data), number_of_generations*generation_size)
#   
# })