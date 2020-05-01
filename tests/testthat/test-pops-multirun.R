context("test-pops-multirun")

test_that("Multirun model outputs work", {
  infected_file <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file <- system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  total_plants_file <- system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  temperature_file <- ""
  temperature_coefficient_file <- ""
  precipitation_coefficient_file <-""
  use_lethal_temperature <- FALSE
  temp <- FALSE
  precip <- FALSE
  season_month_start <- 5
  season_month_end <- 11
  time_step <- "month"
  start_date <- '2019-01-01'
  end_date <- '2019-12-31'
  lethal_temperature <- -35
  lethal_temperature_month <- 1
  random_seed <- 42
  reproductive_rate <- 0.0
  treatments_file <- ""
  treatment_dates <- c('2019-11-01')
  treatment_method <- "ratio"
  management <- FALSE
  mortality_on <- FALSE
  mortality_rate <- 0
  mortality_time_lag <- 0
  percent_natural_dispersal <- 1.0
  natural_kernel_type <- "cauchy"
  anthropogenic_kernel_type <- "cauchy"
  natural_distance_scale <- 23
  anthropogenic_distance_scale <- 0
  natural_dir <- "NONE"
  natural_kappa <- 0
  anthropogenic_dir <- "NONE"
  anthropogenic_kappa <- 0
  pesticide_duration <- c(0)
  pesticide_efficacy <- 1.0
  random_seed = NULL
  output_frequency = "year"
  movements_file <- ""
  use_movements <- FALSE
  num_iterations <- 5
  number_of_cores <- 2
  model_type <- "SI"
  latency_period <- 0
  
  data <- PoPS::pops_multirun(infected_file, host_file, total_plants_file, 
                              temp, temperature_coefficient_file, 
                              precip, precipitation_coefficient_file, 
                              model_type = model_type,
                              latency_period = latency_period,
                              time_step, reproductive_rate,
                              season_month_start, season_month_end, 
                              start_date, end_date, 
                              use_lethal_temperature, temperature_file,
                              lethal_temperature, lethal_temperature_month,
                              mortality_on, mortality_rate, mortality_time_lag, 
                              management, treatment_dates, treatments_file,
                              treatment_method,
                              percent_natural_dispersal,
                              natural_kernel_type, anthropogenic_kernel_type,
                              natural_distance_scale, anthropogenic_distance_scale,
                              natural_dir, natural_kappa, 
                              anthropogenic_dir, anthropogenic_kappa,
                              num_iterations, number_of_cores,
                              pesticide_duration, pesticide_efficacy,
                              random_seed = NULL, output_frequency,
                              movements_file, use_movements)
  
  expect_equal(length(data), 8)
  expect_equal(as.matrix(data$single_run_out[[1]]), as.matrix(raster(infected_file)))
  expect_equal(as.matrix(data$probability[[1]]), matrix(c(100,0,0,0), nrow = 2, ncol = 2))
  expect_equal(data$number_infecteds[[1]], 5)
  expect_equal(data$number_infecteds[[2]], 0)
  expect_equal(data$infected_areas[[1]], 900)
  expect_equal(data$infected_areas[[2]], 0)
  expect_equal(is.nan(data$west_rate[[1]]), T)
  expect_equal(is.na(data$west_rate[[2]]), T)
  expect_equal(data$east_rate[[1]], 0)
  expect_equal(data$east_rate[[2]], 0)
  expect_equal(data$south_rate[[1]], 0)
  expect_equal(data$south_rate[[2]], 0)
  expect_equal(is.nan(data$north_rate[[1]]), T)
  expect_equal(is.na(data$north_rate[[2]]), T)
  
})