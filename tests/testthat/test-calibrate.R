context("test-calibrate")

test_that("Both reproductive rate and short distance scale are decreasing when starting values are much higher than they should converge to!", {
  infected_years_file <- "H:/Team Drives/APHIS  Projects/PoPS/test_data/parameter_estimation_test_data/20x20/sim_rr2.2_short_distance_1.5/infected_years.tif"
  num_iterations <<- 20
  start_reproductive_rate <- 8
  start_short_distance_scale <- 20
  sd_reproductive_rate <- 0.2
  sd_short_distance_scale <- 1
  
  infected_file <- "H:/Team Drives/APHIS  Projects/PoPS/test_data/parameter_estimation_test_data/20x20/initial_infection.tif"
  host_file <- "H:/Team Drives/APHIS  Projects/PoPS/test_data/parameter_estimation_test_data/20x20/host.tif"
  total_plants_file <- "H:/Team Drives/APHIS  Projects/PoPS/test_data/parameter_estimation_test_data/20x20/all_plants.tif"
  treatments_file <- ""
  treatment_years <- c(0)
  management <- FALSE
  
  params3 <- calibrate(infected_years_file, num_interations, start_reproductive_rate, 
                       start_short_distance_scale, sd_reproductive_rate, sd_short_distance_scale,
                       infected_file, host_file, total_plants_file, reproductive_rate =3.0,
                       use_lethal_temperature = FALSE, temp = FALSE, precip = FALSE, management = management, mortality_on = FALSE,
                       temperature_file = "", temperature_coefficient_file = "", 
                       precipitation_coefficient_file ="", treatments_file = treatments_file,
                       season_month_start = 1, season_month_end = 12, time_step = "month",
                       start_time = 2012, end_time = 2016, treatment_years = treatment_years,
                       dispersal_kern = "cauchy", percent_short_distance_dispersal = 1.0,
                       short_distance_scale = 59, long_distance_scale = 0.0,
                       lethal_temperature = -12.87, lethal_temperature_month = 1,
                       mortality_rate = 0, mortality_time_lag = 0,
                       wind_dir = "NONE", kappa = 0)
expect_lte(params3$reproductive_rate[[num_iterations]], start_reproductive_rate)
expect_lte(params3$short_distance_scale[[num_iterations]], start_short_distance_scale)
})


test_that("Both reproductive rate and short distance scale are decreasing when starting values are much higher than they should converge to!", {
  infected_years_file <- "H:/Team Drives/APHIS  Projects/PoPS/test_data/parameter_estimation_test_data/20x20/sim_rr1.3_short_distance_43/infected_years.tif"
  num_iterations <<- 20
  start_reproductive_rate <- 0.2
  start_short_distance_scale <- 5
  sd_reproductive_rate <- 0.2
  sd_short_distance_scale <- 1
  
  infected_file <- "H:/Team Drives/APHIS  Projects/PoPS/test_data/parameter_estimation_test_data/20x20/initial_infection.tif"
  host_file <- "H:/Team Drives/APHIS  Projects/PoPS/test_data/parameter_estimation_test_data/20x20/host.tif"
  total_plants_file <- "H:/Team Drives/APHIS  Projects/PoPS/test_data/parameter_estimation_test_data/20x20/all_plants.tif"
  treatments_file <- ""
  treatment_years <- c(0)
  management <- FALSE
  
  params2 <- calibrate(infected_years_file, num_interations, start_reproductive_rate, 
                       start_short_distance_scale, sd_reproductive_rate, sd_short_distance_scale,
                       infected_file, host_file, total_plants_file, reproductive_rate =3.0,
                       use_lethal_temperature = FALSE, temp = FALSE, precip = FALSE, management = management, mortality_on = FALSE,
                       temperature_file = "", temperature_coefficient_file = "", 
                       precipitation_coefficient_file ="", treatments_file = treatments_file,
                       season_month_start = 1, season_month_end = 12, time_step = "month",
                       start_time = 2012, end_time = 2016, treatment_years = treatment_years,
                       dispersal_kern = "cauchy", percent_short_distance_dispersal = 1.0,
                       short_distance_scale = 59, long_distance_scale = 0.0,
                       lethal_temperature = -12.87, lethal_temperature_month = 1,
                       mortality_rate = 0, mortality_time_lag = 0,
                       wind_dir = "NONE", kappa = 0)
  expect_gte(params2$reproductive_rate[[num_iterations]], start_reproductive_rate)
  expect_gte(params2$short_distance_scale[[num_iterations]], start_short_distance_scale)
})
