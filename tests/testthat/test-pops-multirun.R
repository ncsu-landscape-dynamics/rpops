context("test-pops-multirun")


test_that("Model stops if files don't exist or aren't the correct extension", {
  infected_file <-
    system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  parameter_means <- c(0, 21, 1, 500, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 6, ncol = 6)

  expect_error(pops_multirun(infected_file = "",
                    host_file =  host_file,
                    total_populations_file =  host_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix),
               "file does not exist")
})

test_that("Multirun model outputs work", {
  skip_on_os("windows")
  infected_file <-
    system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  temperature_file <- ""
  temperature_coefficient_file <- ""
  precipitation_coefficient_file <- ""
  use_lethal_temperature <- FALSE
  temp <- FALSE
  precip <- FALSE
  season_month_start <- 5
  season_month_end <- 11
  time_step <- "month"
  start_date <- "2019-01-01"
  end_date <- "2019-12-31"
  lethal_temperature <- -35
  lethal_temperature_month <- 1
  random_seed <- 42
  treatments_file <- ""
  treatment_dates <- c("2019-11-01")
  treatment_method <- "ratio"
  management <- FALSE
  mortality_on <- FALSE
  mortality_rate <- 0
  mortality_time_lag <- 0
  mortality_frequency <- "Year"
  mortality_frequency_n <- 1
  natural_kernel_type <- "cauchy"
  anthropogenic_kernel_type <- "cauchy"
  natural_dir <- "NONE"
  anthropogenic_dir <- "NONE"
  pesticide_duration <- c(0)
  pesticide_efficacy <- 1.0
  random_seed <- NULL
  output_frequency <- "year"
  movements_file <- ""
  use_movements <- FALSE
  number_of_iterations <- 2
  number_of_cores <- 2
  model_type <- "SI"
  latency_period <- 0
  parameter_means <- c(0, 21, 1, 500, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 6, ncol = 6)
  start_exposed <- FALSE
  generate_stochasticity <- TRUE
  establishment_stochasticity <- TRUE
  movement_stochasticity <- TRUE
  deterministic <- FALSE
  establishment_probability <- 0.5
  dispersal_percentage <- 0.99
  quarantine_areas_file <- ""
  use_quarantine <- FALSE
  output_frequency_n <- 1
  use_spreadrates <- TRUE
  use_overpopulation_movements <- FALSE
  overpopulation_percentage <- 0
  leaving_percentage <- 0
  leaving_scale_coefficient <- 1
  exposed_file <- ""
  mask <- NULL
  write_outputs <- "all_simulations"
  output_folder_path <- tempdir()

  data <- pops_multirun(infected_file,
                        host_file,
                        total_populations_file,
                        parameter_means,
                        parameter_cov_matrix,
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
                        mortality_frequency,
                        mortality_frequency_n,
                        management,
                        treatment_dates,
                        treatments_file,
                        treatment_method,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        number_of_iterations,
                        number_of_cores,
                        pesticide_duration,
                        pesticide_efficacy,
                        random_seed,
                        output_frequency,
                        output_frequency_n,
                        movements_file,
                        use_movements,
                        start_exposed,
                        generate_stochasticity,
                        establishment_stochasticity,
                        movement_stochasticity,
                        deterministic,
                        establishment_probability,
                        dispersal_percentage,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)



  expect_equal(length(data), 19)
  expect_equal(terra::as.matrix(data$single_run[[1]], wide = TRUE),
               terra::as.matrix(terra::rast(infected_file), wide = TRUE))
  expect_equal(terra::as.matrix(data$susceptible_run[[1]], wide = TRUE),
               matrix(c(10, 6, 14, 15), nrow = 2, ncol = 2))
  expect_equal(terra::as.matrix(data$probability[[1]], wide = TRUE),
               matrix(c(100, 0, 0, 0), nrow = 2, ncol = 2))
  expect_equal(data$number_infecteds[[1]], 5)
  expect_equal(data$number_infecteds[[2]], 0)
  expect_equal(data$infected_areas[[1]], 900)
  expect_equal(data$infected_areas[[2]], 0)
  expect_equal(data$west_rate[[1]], 0)
  expect_equal(data$west_rate[[2]], 0)
  expect_equal(data$east_rate[[1]], 0)
  expect_equal(data$east_rate[[2]], 0)
  expect_equal(data$south_rate[[1]], 0)
  expect_equal(data$south_rate[[2]], 0)
  expect_equal(data$north_rate[[1]], 0)
  expect_equal(data$north_rate[[2]], 0)

  output_frequency <- "month"

  data <- pops_multirun(infected_file,
                        host_file,
                        total_populations_file,
                        parameter_means,
                        parameter_cov_matrix,
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
                        mortality_frequency,
                        mortality_frequency_n,
                        management,
                        treatment_dates,
                        treatments_file,
                        treatment_method,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        number_of_iterations,
                        number_of_cores,
                        pesticide_duration,
                        pesticide_efficacy,
                        random_seed,
                        output_frequency,
                        output_frequency_n,
                        movements_file,
                        use_movements,
                        start_exposed,
                        generate_stochasticity,
                        establishment_stochasticity,
                        movement_stochasticity,
                        deterministic,
                        establishment_probability,
                        dispersal_percentage,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)


  expect_equal(length(data), 19)
  expect_equal(terra::as.matrix(data$single_run[[1]], wide = TRUE),
               terra::as.matrix(terra::rast(infected_file), wide = TRUE))
  expect_equal(terra::as.matrix(data$susceptible_run[[1]], wide = TRUE),
               matrix(c(10, 6, 14, 15), nrow = 2, ncol = 2))
  expect_equal(terra::as.matrix(data$probability[[1]], wide = TRUE),
               matrix(c(100, 0, 0, 0), nrow = 2, ncol = 2))
  expect_equal(data$number_infecteds[[1]], 5)
  expect_equal(data$number_infecteds[[2]], 0)
  expect_equal(data$infected_areas[[1]], 900)
  expect_equal(data$infected_areas[[2]], 0)
  expect_equal(data$west_rate[[1]], 0)
  expect_equal(data$west_rate[[2]], 0)
  expect_equal(data$east_rate[[1]], 0)
  expect_equal(data$east_rate[[2]], 0)
  expect_equal(data$south_rate[[1]], 0)
  expect_equal(data$south_rate[[2]], 0)
  expect_equal(data$north_rate[[1]], 0)
  expect_equal(data$north_rate[[2]], 0)

  output_frequency <- "year"
  quarantine_areas_file <-
    system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  use_quarantine <- TRUE

  data <- pops_multirun(infected_file,
                        host_file,
                        total_populations_file,
                        parameter_means,
                        parameter_cov_matrix,
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
                        mortality_frequency,
                        mortality_frequency_n,
                        management,
                        treatment_dates,
                        treatments_file,
                        treatment_method,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        number_of_iterations,
                        number_of_cores,
                        pesticide_duration,
                        pesticide_efficacy,
                        random_seed,
                        output_frequency,
                        output_frequency_n,
                        movements_file,
                        use_movements,
                        start_exposed,
                        generate_stochasticity,
                        establishment_stochasticity,
                        movement_stochasticity,
                        deterministic,
                        establishment_probability,
                        dispersal_percentage,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)


  expect_equal(length(data), 19)
  expect_equal(terra::as.matrix(data$single_run[[1]], wide = TRUE),
               terra::as.matrix(terra::rast(infected_file), wide = TRUE))
  expect_equal(terra::as.matrix(data$susceptible_run[[1]], wide = TRUE),
               matrix(c(10, 6, 14, 15), nrow = 2, ncol = 2))
  expect_equal(terra::as.matrix(data$probability[[1]], wide = TRUE),
               matrix(c(100, 0, 0, 0), nrow = 2, ncol = 2))
  expect_equal(data$number_infecteds[[1]], 5)
  expect_equal(data$number_infecteds[[2]], 0)
  expect_equal(data$infected_areas[[1]], 900)
  expect_equal(data$infected_areas[[2]], 0)
  expect_equal(data$west_rate[[1]], 0)
  expect_equal(data$west_rate[[2]], 0)
  expect_equal(data$east_rate[[1]], 0)
  expect_equal(data$east_rate[[2]], 0)
  expect_equal(data$south_rate[[1]], 0)
  expect_equal(data$south_rate[[2]], 0)
  expect_equal(data$north_rate[[1]], 0)
  expect_equal(data$north_rate[[2]], 0)

  output_frequency <- "month"

  data <- pops_multirun(infected_file,
                        host_file,
                        total_populations_file,
                        parameter_means,
                        parameter_cov_matrix,
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
                        mortality_frequency,
                        mortality_frequency_n,
                        management,
                        treatment_dates,
                        treatments_file,
                        treatment_method,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        number_of_iterations,
                        number_of_cores,
                        pesticide_duration,
                        pesticide_efficacy,
                        random_seed,
                        output_frequency,
                        output_frequency_n,
                        movements_file,
                        use_movements,
                        start_exposed,
                        generate_stochasticity,
                        establishment_stochasticity,
                        movement_stochasticity,
                        deterministic,
                        establishment_probability,
                        dispersal_percentage,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)

  expect_equal(length(data), 19)
  expect_equal(terra::as.matrix(data$single_run[[1]], wide = TRUE),
               terra::as.matrix(terra::rast(infected_file), wide = TRUE))
  expect_equal(terra::as.matrix(data$susceptible_run[[1]], wide = TRUE),
               matrix(c(10, 6, 14, 15), nrow = 2, ncol = 2))
  expect_equal(terra::as.matrix(data$probability[[1]], wide = TRUE),
               matrix(c(100, 0, 0, 0), nrow = 2, ncol = 2))
  expect_equal(data$number_infecteds[[1]], 5)
  expect_equal(data$number_infecteds[[2]], 0)
  expect_equal(data$infected_areas[[1]], 900)
  expect_equal(data$infected_areas[[2]], 0)
  expect_equal(data$west_rate[[1]], 0)
  expect_equal(data$west_rate[[2]], 0)
  expect_equal(data$east_rate[[1]], 0)
  expect_equal(data$east_rate[[2]], 0)
  expect_equal(data$south_rate[[1]], 0)
  expect_equal(data$south_rate[[2]], 0)
  expect_equal(data$north_rate[[1]], 0)
  expect_equal(data$north_rate[[2]], 0)

})

test_that("Multirun model outputs work with mask", {
  skip_on_os("windows")
  infected_file <-
    system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
  host_file <-
    system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
  temperature_file <- ""
  temperature_coefficient_file <- ""
  precipitation_coefficient_file <- ""
  use_lethal_temperature <- FALSE
  temp <- FALSE
  precip <- FALSE
  season_month_start <- 5
  season_month_end <- 11
  time_step <- "month"
  start_date <- "2019-01-01"
  end_date <- "2019-12-31"
  lethal_temperature <- -35
  lethal_temperature_month <- 1
  random_seed <- 42
  treatments_file <- ""
  treatment_dates <- c("2019-11-01")
  treatment_method <- "ratio"
  management <- FALSE
  mortality_on <- FALSE
  mortality_rate <- 0
  mortality_time_lag <- 0
  mortality_frequency <- "Year"
  mortality_frequency_n <- 1
  natural_kernel_type <- "cauchy"
  anthropogenic_kernel_type <- "cauchy"
  natural_dir <- "NONE"
  anthropogenic_dir <- "NONE"
  pesticide_duration <- c(0)
  pesticide_efficacy <- 1.0
  random_seed <- NULL
  output_frequency <- "year"
  movements_file <- ""
  use_movements <- FALSE
  number_of_iterations <- 2
  number_of_cores <- 2
  model_type <- "SI"
  latency_period <- 0
  parameter_means <- c(0, 21, 1, 500, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 6, ncol = 6)
  start_exposed <- FALSE
  generate_stochasticity <- TRUE
  establishment_stochasticity <- TRUE
  movement_stochasticity <- TRUE
  deterministic <- FALSE
  establishment_probability <- 0.5
  dispersal_percentage <- 0.99
  quarantine_areas_file <- ""
  use_quarantine <- FALSE
  output_frequency_n <- 1
  use_spreadrates <- TRUE
  mask <- system.file("extdata", "simple20x20", "mask.tif", package = "PoPS")

  data <- pops_multirun(infected_file,
                        host_file,
                        total_populations_file,
                        parameter_means,
                        parameter_cov_matrix,
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
                        mortality_frequency,
                        mortality_frequency_n,
                        management,
                        treatment_dates,
                        treatments_file,
                        treatment_method,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        number_of_iterations,
                        number_of_cores,
                        pesticide_duration,
                        pesticide_efficacy,
                        random_seed,
                        output_frequency,
                        output_frequency_n,
                        movements_file,
                        use_movements,
                        start_exposed,
                        generate_stochasticity,
                        establishment_stochasticity,
                        movement_stochasticity,
                        deterministic,
                        establishment_probability,
                        dispersal_percentage,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates,
                        use_overpopulation_movements = FALSE,
                        overpopulation_percentage = 0,
                        leaving_percentage = 0,
                        leaving_scale_coefficient = 1,
                        exposed_file = "",
                        mask = mask)

  expect_equal(length(data), 19)
  expect_equal(terra::as.matrix(data$single_run[[1]], wide = TRUE),
               terra::as.matrix(terra::rast(infected_file), wide = TRUE))
  expect_equal(data$number_infecteds[[1]], 1)
  expect_equal(data$number_infecteds[[2]], 0)
  expect_equal(data$infected_areas[[1]], 10000)
  expect_equal(data$infected_areas[[2]], 0)
  expect_equal(data$west_rate[[1]], 0)
  expect_equal(data$west_rate[[2]], 0)
  expect_equal(data$east_rate[[1]], 0)
  expect_equal(data$east_rate[[2]], 0)
  expect_equal(data$south_rate[[1]], 0)
  expect_equal(data$south_rate[[2]], 0)
  expect_equal(data$north_rate[[1]], 0)
  expect_equal(data$north_rate[[2]], 0)
})
