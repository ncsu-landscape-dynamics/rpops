# context("test-calibrate")
# 
# test_that("Both reproductive rate and short distance scale are decreasing when starting values are much higher than they should converge to!", {
#   infected_years_file <- system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
#   num_iterations <- 50
#   number_of_observations <- 20
#   prior_number_of_observations <- 5
#   prior_reproductive_rate <- c(8, 0.5)
#   prior_natural_distance_scale <- c(80, 5)
#   prior_percent_natural_dispersal <- c(1.0,0) 
#   prior_anthropogenic_distance_scale <- c(1000,0)
#   number_of_cores <- 11
#   
#   infected_file = system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
#   host_file = system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
#   total_plants_file = system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
#   treatments_file <- ""
#   treatment_years <- c(0)
#   management <- FALSE
#   
#   params3 <- calibrate(infected_years_file, num_iterations,  number_of_cores = NA,
#                         number_of_observations, prior_number_of_observations,
#                         prior_reproductive_rate,
#                         prior_natural_distance_scale,
#                         prior_percent_natural_dispersal, 
#                         prior_anthropogenic_distance_scale,
#                         infected_file, host_file, total_plants_file, 
#                         temp = FALSE, temperature_coefficient_file = "", 
#                         precip = FALSE, precipitation_coefficient_file = "", 
#                         time_step = "month", 
#                         season_month_start = 1, season_month_end = 12, 
#                         start_date = '2012-01-01', end_date = '2016-12-31', 
#                         use_lethal_temperature = FALSE, temperature_file = "",
#                         lethal_temperature = -12.87, lethal_temperature_month = 1,
#                         mortality_on = FALSE, mortality_rate = 0, mortality_time_lag = 0, 
#                         management = FALSE, treatment_dates = c("2012-01-01"), treatments_file = "",
#                         treatment_method = "ratio",
#                         natural_kernel_type = "cauchy", anthropogenic_kernel_type = "cauchy",
#                         natural_dir = "NONE", natural_kappa = 0, 
#                         anthropogenic_dir = "NONE", anthropogenic_kappa = 0,
#                         pesticide_duration = c(0), pesticide_efficacy = 1.0,
#                         mask = NULL, success_metric = "quantity", output_frequency = "year",
#                         movements_file = "", use_movements = FALSE)
#   
# mode_reproductive_rate <- params3$posterior_reproductive_rates  
# mode_reproductive_rate <- mode_reproductive_rate$rate[mode_reproductive_rate$posterior_probability == max(mode_reproductive_rate$posterior_probability)]
# 
# mode_natural_distance <- params3$posterior_natural_distance_scales 
# mode_natural_distance <- mode_natural_distance$rate[mode_natural_distance$posterior_probability == max(mode_natural_distance$posterior_probability)]
# 
# expect_lte(mode_reproductive_rate, prior_reproductive_rate[1])
# expect_lte(mode_natural_distance, prior_natural_distance_scale[1])
# })
# 
# 
# test_that("Both reproductive rate and short distance scale are decreasing when starting values are much higher than they should converge to!", {
#   infected_years_file <- system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
#   num_iterations <- 50
#   number_of_observations <- 20
#   prior_number_of_observations <- 5
#   prior_reproductive_rate <- c(0.1, 0.2)
#   prior_natural_distance_scale <- c(1, 3)
#   prior_percent_natural_dispersal <- c(1.0,0) 
#   prior_anthropogenic_distance_scale <- c(1000,0)
#   number_of_cores <- 11
#   
#   infected_file = system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
#   host_file = system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
#   total_plants_file = system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
#   treatments_file <- ""
#   treatment_years <- c(0)
#   management <- FALSE
#   
#   params2 <- calibrate(infected_years_file, num_iterations,  number_of_cores = NA,
#                        number_of_observations, prior_number_of_observations,
#                        prior_reproductive_rate,
#                        prior_natural_distance_scale,
#                        prior_percent_natural_dispersal, 
#                        prior_anthropogenic_distance_scale,
#                        infected_file, host_file, total_plants_file, 
#                        temp = FALSE, temperature_coefficient_file = "", 
#                        precip = FALSE, precipitation_coefficient_file = "", 
#                        time_step = "month", 
#                        season_month_start = 1, season_month_end = 12, 
#                        start_date = '2012-01-01', end_date = '2016-12-31', 
#                        use_lethal_temperature = FALSE, temperature_file = "",
#                        lethal_temperature = -12.87, lethal_temperature_month = 1,
#                        mortality_on = FALSE, mortality_rate = 0, mortality_time_lag = 0, 
#                        management = FALSE, treatment_dates = c("2012-01-01"), treatments_file = "",
#                        treatment_method = "ratio",
#                        natural_kernel_type = "cauchy", anthropogenic_kernel_type = "cauchy",
#                        natural_dir = "NONE", natural_kappa = 0, 
#                        anthropogenic_dir = "NONE", anthropogenic_kappa = 0,
#                        pesticide_duration = c(0), pesticide_efficacy = 1.0,
#                        mask = NULL, success_metric = "quantity", output_frequency = "year",
#                        movements_file = "", use_movements = FALSE)
#   
#   mode_reproductive_rate <- params2$posterior_reproductive_rates  
#   mode_reproductive_rate <- mode_reproductive_rate$rate[mode_reproductive_rate$posterior_probability == max(mode_reproductive_rate$posterior_probability)]
#   
#   mode_natural_distance <- params2$posterior_natural_distance_scales 
#   mode_natural_distance <- mode_natural_distance$rate[mode_natural_distance$posterior_probability == max(mode_natural_distance$posterior_probability)]
#   
#   expect_gte(mode_reproductive_rate, prior_reproductive_rate[1])
#   expect_gte(mode_natural_distance, prior_natural_distance_scale[1])
# })
