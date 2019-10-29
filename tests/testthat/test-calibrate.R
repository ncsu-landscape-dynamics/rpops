context("test-calibrate")

test_that("Both reproductive rate and short distance scale are decreasing when starting values are much higher than they should converge to!", {
  infected_years_file <- system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  # infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/test_data/parameter_estimation_test_data/20x20/sim_rr2.2_short_distance_1.5/infected_years.tif"
  num_iterations <- 50
  start_reproductive_rate <- 8
  start_natural_distance_scale <- 50
  sd_reproductive_rate <- 0.2
  sd_natural_distance_scale <- 4
  number_of_cores <- 11
  
  infected_file = system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
  # infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/test_data/parameter_estimation_test_data/20x20/initial_infection.tif"
  host_file = system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
  # host_file <- "H:/Shared drives/APHIS  Projects/PoPS/test_data/parameter_estimation_test_data/20x20/host.tif"
  total_plants_file = system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
  # total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/test_data/parameter_estimation_test_data/20x20/all_plants.tif"
  treatments_file <- ""
  treatment_years <- c(0)
  management <- FALSE
  
  params3 <- calibrate(infected_years_file, num_iterations, start_reproductive_rate, number_of_cores,
                       start_natural_distance_scale, sd_reproductive_rate, sd_natural_distance_scale,
                       infected_file, host_file, total_plants_file, 
                       temp = FALSE, temperature_coefficient_file = "", 
                       precip = FALSE, precipitation_coefficient_file = "", 
                       time_step = 'month', reproductive_rate = 3.0,
                       season_month_start = 1, season_month_end = 12, 
                       start_time =2012, end_time =2016, 
                       use_lethal_temperature = FALSE, temperature_file = "",
                       lethal_temperature = -12.87, lethal_temperature_month = 1,
                       mortality_on = FALSE, mortality_rate = 0, mortality_time_lag = 0, 
                       management, treatment_years = treatment_years, treatments_file = treatments_file,
                       treatment_method = "ratio", treatment_month = 12,
                       percent_natural_dispersal = 1.0,
                       natural_kernel_type = "cauchy", anthropogenic_kernel_type = "cauchy",
                       natural_distance_scale = 59, anthropogenic_distance_scale = 0.0,
                       natural_dir = "NONE", natural_kappa = 0, 
                       anthropogenic_dir = "NONE", anthropogenic_kappa = 0,
                       mask, success_metric = "quantity")


expect_lte(params3$reproductive_rate[[num_iterations]], start_reproductive_rate)
expect_lte(params3$natural_distance_scale[[num_iterations]], start_natural_distance_scale)
})


test_that("Both reproductive rate and short distance scale are decreasing when starting values are much higher than they should converge to!", {
  infected_years_file <- system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  # infected_years_file <- "H:/Shared drives/APHIS  Projects/PoPS/test_data/parameter_estimation_test_data/20x20/sim_rr1.3_short_distance_43/infected_years.tif"
  num_iterations <- 50
  start_reproductive_rate <- 0.1
  start_natural_distance_scale <- 1
  sd_reproductive_rate <- 0.2
  sd_natural_distance_scale <- 1
  
  infected_file = system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
  # infected_file <- "H:/Shared drives/APHIS  Projects/PoPS/test_data/parameter_estimation_test_data/20x20/initial_infection.tif"
  host_file = system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
  # host_file <- "H:/Shared drives/APHIS  Projects/PoPS/test_data/parameter_estimation_test_data/20x20/host.tif"
  total_plants_file = system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
  # total_plants_file <- "H:/Shared drives/APHIS  Projects/PoPS/test_data/parameter_estimation_test_data/20x20/all_plants.tif"
  treatments_file <- ""
  treatment_years <- c(0)
  management <- FALSE
  
  params2 <- calibrate(infected_years_file, num_iterations, start_reproductive_rate, number_of_cores,
                       start_natural_distance_scale, sd_reproductive_rate, sd_natural_distance_scale,
                       infected_file, host_file, total_plants_file, 
                       temp = FALSE, temperature_coefficient_file = "", 
                       precip = FALSE, precipitation_coefficient_file = "", 
                       time_step = 'month', reproductive_rate = 3.0,
                       season_month_start = 1, season_month_end = 12, 
                       start_time =2012, end_time =2016, 
                       use_lethal_temperature = FALSE, temperature_file = "",
                       lethal_temperature = -12.87, lethal_temperature_month = 1,
                       mortality_on = FALSE, mortality_rate = 0, mortality_time_lag = 0, 
                       management, treatment_years = treatment_years, treatments_file = treatments_file,
                       treatment_method = "ratio", treatment_month = 12,
                       percent_natural_dispersal = 1.0,
                       natural_kernel_type = "cauchy", anthropogenic_kernel_type = "cauchy",
                       natural_distance_scale = 59, anthropogenic_distance_scale = 0.0,
                       natural_dir = "NONE", natural_kappa = 0, 
                       anthropogenic_dir = "NONE", anthropogenic_kappa = 0,
                       mask, success_metric = "quantity")
  expect_gte(params2$reproductive_rate[[num_iterations]], start_reproductive_rate)
  expect_gte(params2$short_distance_scale[[num_iterations]], start_short_distance_scale)
})
