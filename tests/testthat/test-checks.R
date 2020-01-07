context("test-checks")

test_that("Initial raster checks report correct errors and return a raster", {
  infected_file <- ""
  infected <- initial_raster_checks(infected_file)
  expect_equal(infected$checks_passed, FALSE)
  expect_equal(infected$failed_check, "file does not exist")
  
  infected_file <-  system.file("extdata", "simple2x2", "infected.csv", package = "PoPS")
  infected <- initial_raster_checks(infected_file)
  expect_equal(infected$checks_passed, FALSE)
  expect_equal(infected$failed_check, "file is not one of '.grd', '.tif', '.img'")
  
  infected_file <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  infected <- initial_raster_checks(infected_file)
  expect_equal(infected$checks_passed, TRUE)
  expect_equal(infected$raster[], raster(infected_file)[])
})

test_that("Initial raster checks report correct errors and return a raster", {
  infected_file <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  infected <- raster(infected_file)
  host_file <- ""
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check, "file does not exist")
  
  host_file <-  system.file("extdata", "simple2x2", "infected.csv", package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check, "file is not one of '.grd', '.tif', '.img'")
  
  host_file <- system.file("extdata", "simple5x5", "total_plants.tif", package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check, "Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
  
  host_file <- system.file("extdata", "simple2x2", "total_plants_diff_res.tif", package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check, "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  
  host_file <- system.file("extdata", "simple2x2", "total_plants_diff_xres.tif", package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check, "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  
  host_file <- system.file("extdata", "simple2x2", "total_plants_diff_yres.tif", package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check, "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  
  host_file <- system.file("extdata", "simple2x2", "total_plants_with_crs.tif", package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check, "Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
  
  host_file <-  system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, TRUE)
  expect_equal(host$raster[], raster(host_file)[])
  
})

test_that("Treatment checks report correct errors and return a matrices of treatments if all checks pass", {
  treatments_file <- system.file("extdata", "simple2x2", "treatments.tif", package = "PoPS")
  treatment_stack <- stack(treatments_file)
  pesticide_duration <- c(60)
  treatment_dates <- c('2016-01-01', '2016-03-15')
  pesticide_efficacy <- 1.0
  treatment_check <- treatment_checks(treatment_stack, treatments_file, pesticide_duration, treatment_dates, pesticide_efficacy)
  expect_equal(treatment_check$checks_passed, FALSE)
  expect_equal(treatment_check$failed_check, "Length of list for treatment dates and treatments_file must be equal")
  
  treatments_file <- system.file("extdata", "simple2x2", "treatments.tif", package = "PoPS")
  treatment_stack <- stack(treatments_file)
  pesticide_duration <- c(60,45,60)
  treatment_dates <- c('2016-01-01')
  pesticide_efficacy <- 1.0
  treatment_check <- treatment_checks(treatment_stack, treatments_file, pesticide_duration, treatment_dates, pesticide_efficacy)
  expect_equal(treatment_check$checks_passed, FALSE)
  expect_equal(treatment_check$failed_check, "Length of list for treatment dates and pesticide_duration must be equal")
  
  treatments_file <- system.file("extdata", "simple2x2", "treatments.tif", package = "PoPS")
  treatment_stack <- stack(treatments_file)
  pesticide_duration <- c(60)
  treatment_dates <- c('2016-01-01')
  pesticide_efficacy <- 1.0
  treatment_check <- treatment_checks(treatment_stack, treatments_file, pesticide_duration, treatment_dates, pesticide_efficacy)
  expect_equal(treatment_check$checks_passed, TRUE)
  expect_equal(treatment_check$treatment_maps[[1]], as.matrix(treatment_stack[[1]]))
  
  treatments_file <- c(system.file("extdata", "simple2x2", "treatments.tif", package = "PoPS"), system.file("extdata", "simple2x2", "treatments.tif", package = "PoPS"))
  treatment_stack <- stack(treatments_file)
  pesticide_duration <- c(60, 45)
  treatment_dates <- c('2016-01-01', '2016-10-02')
  pesticide_efficacy <- 1.0
  treatment_check <- treatment_checks(treatment_stack, treatments_file, pesticide_duration, treatment_dates, pesticide_efficacy)
  expect_equal(treatment_check$checks_passed, TRUE)
  expect_equal(treatment_check$treatment_maps[[1]], as.matrix(treatment_stack[[1]]))
})


test_that("Treatment method checks report correct errors and TRUE otherwise", {
  treatment_method <- "Test"
  methods_check <- treatment_metric_checks(treatment_method)
  expect_equal(methods_check$checks_passed, FALSE)
  expect_equal(methods_check$failed_check, "treatment method is not one of the valid treatment options")
  treatment_method <- "ratio"
  methods_check <- treatment_metric_checks(treatment_method)
  expect_equal(methods_check$checks_passed, TRUE)
  treatment_method <- "all infected"
  methods_check <- treatment_metric_checks(treatment_method)
  expect_equal(methods_check$checks_passed, TRUE)
  
})

test_that("Success metric checks report correct errors and return configuration disagreement boolean otherwise", {
  success_metric <- "Test"
  metric_check <- metric_checks(success_metric)
  expect_equal(metric_check$checks_passed, FALSE)
  expect_equal(metric_check$failed_check, "Success metric must be one of 'quantity', 'quantity and configuration', 'residual error', or 'odds ratio'")
  success_metric <- "quantity"
  metric_check <- metric_checks(success_metric)
  expect_equal(metric_check$checks_passed, TRUE)
  expect_equal(metric_check$configuration, FALSE)
  success_metric <- "quantity and configuration"
  metric_check <- metric_checks(success_metric)
  expect_equal(metric_check$checks_passed, TRUE)
  expect_equal(metric_check$configuration, TRUE)
  success_metric <- "odds ratio"
  metric_check <- metric_checks(success_metric)
  expect_equal(metric_check$checks_passed, TRUE)
  expect_equal(metric_check$configuration, FALSE)
  success_metric <- "residual error"
  metric_check <- metric_checks(success_metric)
  expect_equal(metric_check$checks_passed, TRUE)
  expect_equal(metric_check$configuration, FALSE)
})

test_that("Percent checks report correct errors and return whether or not to use the anthropogenic distance kernel", {
  percent_natural_dispersal <- 1.01
  percent_check <- percent_checks(percent_natural_dispersal)
  expect_equal(percent_check$checks_passed, FALSE)
  expect_equal(percent_check$failed_check, "Percent natural dispersal must be between 0.0 and 1.0")
  percent_natural_dispersal <- -0.2
  percent_check <- percent_checks(percent_natural_dispersal)
  expect_equal(percent_check$checks_passed, FALSE)
  expect_equal(percent_check$failed_check, "Percent natural dispersal must be between 0.0 and 1.0")
  percent_natural_dispersal <- 1.0
  percent_check <- percent_checks(percent_natural_dispersal)
  expect_equal(percent_check$checks_passed, TRUE)
  expect_equal(percent_check$use_anthropogenic_kernel, FALSE)
  percent_natural_dispersal <- 0.80
  percent_check <- percent_checks(percent_natural_dispersal)
  expect_equal(percent_check$checks_passed, TRUE)
  expect_equal(percent_check$use_anthropogenic_kernel, TRUE)

})

test_that("Time checks report correct errors and return number of time steps, number of years, and number of outputs otherwise", {
  end_date <- 2017-01-01
  start_date <- 2016-01-01
  time_step <- 'hour'
  output_frequency <- 'week'
  time_check <- time_checks(end_date, start_date, time_step, output_frequency)
  expect_equal(time_check$checks_passed, FALSE)
  expect_equal(time_check$failed_check, "Time step must be one of 'week', 'month' or 'day'")
  
  end_date <- 2017-01-01
  start_date <- 2016-01-01
  time_step <- 'week'
  output_frequency <- 'week'
  time_check <- time_checks(end_date, start_date, time_step, output_frequency)
  expect_equal(time_check$checks_passed, FALSE)
  expect_equal(time_check$failed_check, "End time and/or start time not of type numeric and/or in format YYYY-MM-DD")
  
  end_date <- '2017-01-01'
  start_date <- '2016-01-01'
  time_step <- 'month'
  output_frequency <- 'hour'
  time_check <- time_checks(end_date, start_date, time_step, output_frequency)
  expect_equal(time_check$checks_passed, FALSE)
  expect_equal(time_check$failed_check, "Output frequency must be one of 'week', 'month' or 'day'")
  
  end_date <- '2017-01-01'
  start_date <- '2016-01-01'
  time_step <- 'month'
  output_frequency <- 'week'
  time_check <- time_checks(end_date, start_date, time_step, output_frequency)
  expect_equal(time_check$checks_passed, FALSE)
  expect_equal(time_check$failed_check, "Output frequency is more frequent than time_step. The minimum output_frequency you can use is the time_step of your simulation. You can set the output_frequency to 'time_step' to default to most frequent output possible")

  end_date <- '2017-01-01'
  start_date <- '2016-01-01'
  time_step <- 'week'
  output_frequency <- 'week'
  time_check <- time_checks(end_date, start_date, time_step, output_frequency)
  expect_equal(time_check$checks_passed, TRUE)
  expect_equal(time_check$number_of_time_steps, 53)
  expect_equal(time_check$number_of_years, 1)
  expect_equal(time_check$number_of_outputs, 53)
  
  end_date <- '2016-02-01'
  start_date <- '2016-01-01'
  time_step <- 'week'
  output_frequency <- 'week'
  time_check <- time_checks(end_date, start_date, time_step, output_frequency)
  expect_equal(time_check$checks_passed, TRUE)
  expect_equal(time_check$number_of_time_steps, 5)
  expect_equal(time_check$number_of_years, 1)
  expect_equal(time_check$number_of_outputs, 5)
  
  end_date <- '2017-01-01'
  start_date <- '2016-01-01'
  time_step <- 'week'
  output_frequency <- 'month'
  time_check <- time_checks(end_date, start_date, time_step, output_frequency)
  expect_equal(time_check$checks_passed, TRUE)
  expect_equal(time_check$number_of_time_steps, 53)
  expect_equal(time_check$number_of_years, 1)
  expect_equal(time_check$number_of_outputs, 12)
  
  end_date <- '2017-01-01'
  start_date <- '2016-01-01'
  time_step <- 'week'
  output_frequency <- 'year'
  time_check <- time_checks(end_date, start_date, time_step, output_frequency)
  expect_equal(time_check$checks_passed, TRUE)
  expect_equal(time_check$number_of_time_steps, 53)
  expect_equal(time_check$number_of_years, 1)
  expect_equal(time_check$number_of_outputs, 1)
  
})

test_that("Prior checks report correct errors and return reformatted priors, and mean and sd of priors for starting location for MCMC", {
  priors <- 2.1
  prior_check <- prior_checks(priors)
  expect_equal(prior_check$checks_passed, FALSE)
  expect_equal(prior_check$failed_check, "Incorrect format for priors")
  
  priors <- data.frame(data = c(1.5,1.6,1.7,1.8,1.9,2.0,2.1), per = c(0.1,0.12,0.13,0.30,0.13,0.12,0.1))
  prior_check <- prior_checks(priors)
  names(priors) <- c('var', 'prob')
  expect_equal(prior_check$checks_passed, TRUE)
  expect_equal(prior_check$priors, priors)
  expect_equal(prior_check$start_priors, 1.8)
  expect_equal(prior_check$sd_priors, sd(priors$var))
  
  priors <- c(2.1, 0.2)
  prior_check <- prior_checks(priors)
  expect_equal(prior_check$checks_passed, TRUE)
  expect_equal(prior_check$priors, matrix(priors, ncol = 2))
  expect_equal(prior_check$start_priors, 2.1)
  expect_equal(prior_check$sd_priors, 0.2)
  
})

test_that("Bayesian checks report correct errors and return all rates as a data frame and posterior rates as a data frame", {
  priors <- c(2.1, 0.2)
  prior_check <- prior_checks(priors)
  prior <- prior_check$priors
  start_priors <- prior_check$start_priors
  sd_priors <- prior_check$sd_priors
  params <- data.frame(data = c(rep(1.5,10), rep(1.5,10), rep(1.6,15), rep(1.7,15), rep(1.8,30), rep(1.9,15), rep(2.0,15), rep(2.1,15), rep(2.2, 15)))
  count <- 10000000
  prior_weight <- 0.05
  weight <- 1 - prior_weight
  step_size <- 0.1
  bounds <- c(0, Inf)
  round_to <- 1
  round_to_digits <- 1
  bayesian_check <- bayesian_checks(prior, start_priors, sd_priors, params, count, prior_weight, weight, step_size, bounds = c(0, Inf), round_to = 1, round_to_digits = 1)
  expect_equal(bayesian_check$checks_passed, TRUE)
  expect_length(bayesian_check$rates, 4)
  expect_length(bayesian_check$posterior_rates, 2)
  expect_lte(mean(bayesian_check$posterior_rates[,1]), 2.1)
  
  
  
})
