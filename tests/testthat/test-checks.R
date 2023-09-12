context("test-checks")

test_that("Initial raster checks report correct errors and return a raster", {
  infected_file <- ""
  infected <- initial_raster_checks(infected_file)
  expect_equal(infected$checks_passed, FALSE)
  expect_equal(infected$failed_check, PoPS:::file_exists_error)

  infected_file <-
    system.file("extdata", "simple2x2", "infected.csv", package = "PoPS")
  infected <- initial_raster_checks(infected_file)
  expect_equal(infected$checks_passed, FALSE)
  expect_equal(infected$failed_check, PoPS:::raster_type_error)

  infected_file <-
    system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  infected <- initial_raster_checks(infected_file)
  expect_equal(infected$checks_passed, TRUE)
  expect_equal(infected$raster[], terra::rast(infected_file)[])
})

test_that("Initial raster checks report correct errors and return a raster", {
  infected_file <-
    system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  infected <- terra::rast(infected_file)
  host_file <- ""
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check, PoPS:::file_exists_error)

  host_file <-
    system.file("extdata", "simple2x2", "infected.csv", package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check, PoPS:::raster_type_error)

  host_file <-
    system.file("extdata", "simple5x5", "total_plants.tif", package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check, PoPS:::extent_error)

  host_file <-
    system.file("extdata", "simple2x2", "total_plants_diff_res.tif",
                package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check, PoPS:::resolution_error)

  host_file <-
    system.file("extdata", "simple2x2", "total_plants_diff_xres.tif",
                package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check, PoPS:::resolution_error)

  host_file <- system.file("extdata", "simple2x2", "total_plants_diff_yres.tif", package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check, PoPS:::resolution_error)

  host_file <-
    system.file("extdata", "simple2x2", "critical_temp_diff_crs.tif", package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check, PoPS:::crs_error)

  host_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, TRUE)
  expect_equal(host$raster[], terra::rast(host_file)[])

})

test_that(
  "Treatment checks report correct errors and return a matrices of treatments if all checks pass", {
    treatments_file <-
      system.file("extdata", "simple2x2", "treatments.tif", package = "PoPS")
    treatment_stack <- terra::rast(treatments_file)
    pesticide_duration <- c(60)
    treatment_dates <- c("2016-01-01", "2016-03-15")
    pesticide_efficacy <- 1.0
    treatment_check <-
      treatment_checks(treatment_stack,
                       treatments_file,
                       pesticide_duration,
                       treatment_dates,
                       pesticide_efficacy)
    expect_equal(treatment_check$checks_passed, FALSE)
    expect_equal(treatment_check$failed_check, treatment_length_error)

    treatments_file <- system.file("extdata", "simple2x2", "treatments.tif", package = "PoPS")
    treatment_stack <- terra::rast(treatments_file)
    pesticide_duration <- c(60, 45, 60)
    treatment_dates <- c("2016-01-01")
    pesticide_efficacy <- 1.0
    treatment_check <-
      treatment_checks(treatment_stack,
                       treatments_file,
                       pesticide_duration,
                       treatment_dates,
                       pesticide_efficacy)
    expect_equal(treatment_check$checks_passed, FALSE)
    expect_equal(treatment_check$failed_check, pesticide_length_error)

    treatments_file <- system.file("extdata", "simple2x2", "treatments.tif", package = "PoPS")
    treatment_stack <- terra::rast(treatments_file)
    pesticide_duration <- c(60)
    treatment_dates <- c("2016-01-01")
    pesticide_efficacy <- 1.0
    treatment_check <-
      treatment_checks(treatment_stack,
                       treatments_file,
                       pesticide_duration,
                       treatment_dates,
                       pesticide_efficacy)
    expect_equal(treatment_check$checks_passed, TRUE)
    expect_equal(treatment_check$treatment_maps[[1]],
                 terra::as.matrix(treatment_stack[[1]], wide = TRUE))

    treatments_file <-
      c(system.file("extdata", "simple2x2", "treatments.tif", package = "PoPS"),
        system.file("extdata", "simple2x2", "treatments.tif", package = "PoPS"))
    treatment_stack <- terra::rast(treatments_file)
    pesticide_duration <- c(60, 45)
    treatment_dates <- c("2016-01-01", "2016-10-02")
    pesticide_efficacy <- 1.0
    treatment_check <-
      treatment_checks(treatment_stack,
                       treatments_file,
                       pesticide_duration,
                       treatment_dates,
                       pesticide_efficacy)
    expect_equal(treatment_check$checks_passed, TRUE)
    expect_equal(treatment_check$treatment_maps[[1]],
                 terra::as.matrix(treatment_stack[[1]], wide = TRUE))
  })


test_that("Treatment method checks report correct errors and TRUE otherwise", {
  treatment_method <- "Test"
  methods_check <- treatment_metric_checks(treatment_method)
  expect_equal(methods_check$checks_passed, FALSE)
  expect_equal(methods_check$failed_check, treatment_option_error)
  treatment_method <- "ratio"
  methods_check <- treatment_metric_checks(treatment_method)
  expect_equal(methods_check$checks_passed, TRUE)
  treatment_method <- "all infected"
  methods_check <- treatment_metric_checks(treatment_method)
  expect_equal(methods_check$checks_passed, TRUE)
})

test_that(
  "Time checks report correct errors and return number of time steps, number of
  years, and number of outputs otherwise", {
    end_date <- 2017
    start_date <- 2016
    time_step <- "hour"
    output_frequency <- "week"
    output_frequency_n <- 1
    time_check <- time_checks(end_date, start_date, time_step, output_frequency, output_frequency_n)
    expect_equal(time_check$checks_passed, FALSE)
    expect_equal(time_check$failed_check, time_step_error)

    end_date <- 2017
    start_date <- 2016
    time_step <- "week"
    output_frequency <- "week"
    time_check <- time_checks(end_date, start_date, time_step, output_frequency, output_frequency_n)
    expect_equal(time_check$checks_passed, FALSE)
    expect_equal(time_check$failed_check, date_format_error)

    end_date <- "2017-01-01"
    start_date <- "2016-01-01"
    time_step <- "month"
    output_frequency <- "hour"
    time_check <- time_checks(end_date, start_date, time_step, output_frequency, output_frequency_n)
    expect_equal(time_check$checks_passed, FALSE)
    expect_equal(time_check$failed_check, frequency_error)

    time_step <- "month"
    output_frequency <- "week"
    time_check <- time_checks(end_date, start_date, time_step, output_frequency, output_frequency_n)
    expect_equal(time_check$checks_passed, FALSE)
    expect_equal(time_check$failed_check, output_frequency_error)

    time_step <- "week"
    output_frequency <- "week"
    time_check <- time_checks(end_date, start_date, time_step, output_frequency, output_frequency_n)
    expect_equal(time_check$checks_passed, TRUE)
    expect_equal(time_check$number_of_time_steps, 53)
    expect_equal(time_check$number_of_years, 1)
    expect_equal(time_check$number_of_outputs, 53)

    time_step <- "week"
    output_frequency <- "month"
    time_check <- time_checks(end_date, start_date, time_step, output_frequency, output_frequency_n)
    expect_equal(time_check$checks_passed, TRUE)
    expect_equal(time_check$number_of_time_steps, 53)
    expect_equal(time_check$number_of_years, 1)
    expect_equal(time_check$number_of_outputs, 12)

    time_step <- "week"
    output_frequency <- "year"
    time_check <- time_checks(end_date, start_date, time_step, output_frequency, output_frequency_n)
    expect_equal(time_check$checks_passed, TRUE)
    expect_equal(time_check$number_of_time_steps, 53)
    expect_equal(time_check$number_of_years, 1)
    expect_equal(time_check$number_of_outputs, 1)

    end_date <- "2016-02-01"
    start_date <- "2016-01-01"
    time_step <- "week"
    output_frequency <- "week"
    time_check <- time_checks(end_date, start_date, time_step, output_frequency, output_frequency_n)
    expect_equal(time_check$checks_passed, TRUE)
    expect_equal(time_check$number_of_time_steps, 5)
    expect_equal(time_check$number_of_years, 1)
    expect_equal(time_check$number_of_outputs, 5)

    end_date <- "2016-02-01"
    start_date <- "2016-01-01"
    time_step <- "week"
    output_frequency <- "time_step"
    time_check <- time_checks(end_date, start_date, time_step, output_frequency, output_frequency_n)
    expect_equal(time_check$checks_passed, TRUE)
    expect_equal(time_check$number_of_time_steps, 5)
    expect_equal(time_check$number_of_years, 1)
    expect_equal(time_check$number_of_outputs, 5)
  })

test_that("Bayesian MNN checks work properly", {
  prior_means <- c(0.2, 23)
  prior_cov_matrix <- matrix(c(1, 1, 1, 1), ncol = 2, nrow = 2)
  calibrated_means <- c(0.3, 21)
  calibrated_cov_matrix <- matrix(c(1, 2, 1, 1), ncol = 2, nrow = 2)
  prior_weight <- 0.05
  weight <- 1 - prior_weight

  mnn_check <-
    bayesian_mnn_checks(prior_means,
                        prior_cov_matrix,
                        calibrated_means,
                        calibrated_cov_matrix,
                        prior_weight,
                        weight)
  expect_equal(mnn_check$checks_passed, TRUE)
  expect_equal(mnn_check$posterior_means, c(0.295, 21.1))
  expect_equal(mnn_check$posterior_cov_matrix, matrix(c(1, 1.95, 1, 1), ncol = 2, nrow = 2))

  prior_means <- c(0.2, 23)
  prior_cov_matrix <- matrix(c(1, 1, 1, 1), ncol = 2, nrow = 2)
  calibrated_means <- c(0.3, 21)
  calibrated_cov_matrix <- matrix(c(1, 2), ncol = 1, nrow = 2)
  prior_weight <- 0.05
  weight <- 1 - prior_weight

  mnn_check <-
    bayesian_mnn_checks(prior_means,
                        prior_cov_matrix,
                        calibrated_means,
                        calibrated_cov_matrix,
                        prior_weight,
                        weight)
  expect_equal(mnn_check$checks_passed, FALSE)
  expect_equal(mnn_check$failed_check, prior_cov_matrix_error)

  prior_means <- c(0.2, 23)
  prior_cov_matrix <- matrix(c(1, 1, 1, 1), ncol = 2, nrow = 2)
  calibrated_means <- c(0.3)
  calibrated_cov_matrix <- matrix(c(1, 2, 1, 1), ncol = 2, nrow = 2)
  prior_weight <- 0.05
  weight <- 1 - prior_weight

  mnn_check <-
    bayesian_mnn_checks(prior_means,
                        prior_cov_matrix,
                        calibrated_means,
                        calibrated_cov_matrix,
                        prior_weight,
                        weight)
  expect_equal(mnn_check$checks_passed, FALSE)
  expect_equal(mnn_check$failed_check, prior_means_error)
})


test_that("Multispecies checks return erros if any of the inputs aren't the
          same length", {

  infected_files <- c("", "")
  host_file <- c("", "")
  total_populations_file <- c("", "")
  parameter_means <- c("", "")
  parameter_cov_matrix <- c("", "")
  temp <- c(FALSE, FALSE)
  temperature_coefficient_file <- c("", "")
  precip <- c(FALSE, FALSE)
  precipitation_coefficient_file <- c("", "")
  model_type <- c("SI", "SI")
  latency_period <- c(0, 0)
  time_step <- c("month", "month")
  season_month_start <- c(1, 1)
  season_month_end <- c(12, 12)
  use_lethal_temperature <- c(FALSE, FALSE)
  temperature_file <- c("", "")
  lethal_temperature <- c(1, 1)
  lethal_temperature_month <- c(1, 1)
  mortality_on <- c(FALSE, FALSE)
  mortality_rate <- c(0.05, 0.05)
  mortality_time_lag <- c(1, 1)
  natural_kernel_type <- c("cauchy", "cauchy")
  anthropogenic_kernel_type <- c("cauchy", "cauchy")
  natural_dir <- c("NONE", "NONE")
  anthropogenic_dir <- c("NONE", "NONE")
  movements_file <- c("", "")
  use_movements <- c(FALSE, FALSE)
  start_exposed <- c(TRUE, TRUE)
  quarantine_areas_file <- c("", "")
  use_quarantine <- c(TRUE, TRUE)
  use_spreadrates <- c(TRUE, TRUE)
  species <- c("species1", "species2")

  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)

  expect_equal(multispecies_check$checks_passed, TRUE)

  infected_files <- c("")
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)

  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(multispecies_check$failed_check, PoPS:::species_length_error)

  infected_files <- c("", "")
  host_file <- c("")
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(multispecies_check$failed_check, PoPS:::host_file_length_error)

  host_file <- c("", "")
  total_populations_file <- c("")
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(multispecies_check$failed_check, PoPS:::total_population_length_error)

  total_populations_file <- c("", "")
  parameter_means <- c("")
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::parameter_length_error)

  parameter_means <- c("", "")
  parameter_cov_matrix <- c("")
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::cov_matrix_length_error)

  parameter_cov_matrix <- c("", "")
  temp <- c(FALSE)
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(multispecies_check$failed_check, PoPS:::temp_length_error)

  temp <- c(FALSE, FALSE)
  temperature_coefficient_file <- c("")
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(multispecies_check$failed_check, PoPS:::temperature_coefficient_length_error)

  temperature_coefficient_file <- c("", "")
  precip <- c(FALSE)
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(multispecies_check$failed_check, PoPS:::precip_length_error)

  precip <- c(FALSE, FALSE)
  precipitation_coefficient_file <- c("")
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::precipitation_coefficient_length_error)

  precipitation_coefficient_file <- c("", "")
  model_type <- c("SI")
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(multispecies_check$failed_check, PoPS:::model_type_length_error)

  model_type <- c("SI", "SI")
  latency_period <- c(0)
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::latency_period_length_error)

  latency_period <- c(0, 0)
  time_step <- c("month")
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(multispecies_check$failed_check, PoPS:::time_step_length_error)

  time_step <- c("month", "month")
  season_month_start <- c(1)
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::season_month_start_length_error)

  season_month_start <- c(1, 1)
  season_month_end <- c(12)
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::season_month_end_length_error)

  season_month_end <- c(12, 12)
  use_lethal_temperature <- c(FALSE)
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::use_lethal_length_error)

  use_lethal_temperature <- c(FALSE, FALSE)
  temperature_file <- c("")
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::temperature_file_length_error)

  temperature_file <- c("", "")
  lethal_temperature <- c(1)
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::lethal_temperature_length_error)

  lethal_temperature <- c(1, 1)
  lethal_temperature_month <- c(1)
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::lethal_temperature_month_length_error)

  lethal_temperature_month <- c(1, 1)
  mortality_on <- c(FALSE)
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::mortality_on_length_error)

  mortality_on <- c(FALSE, FALSE)
  mortality_rate <- c(0.05)
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::mortality_rate_length_error)

  mortality_rate <- c(0.05, 0.05)
  mortality_time_lag <- c(1)
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::mortality_time_lag_length_error)

  mortality_time_lag <- c(1, 1)
  natural_kernel_type <- c("cauchy")
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::natural_kernel_length_error)

  natural_kernel_type <- c("cauchy", "cauchy")
  anthropogenic_kernel_type <- c("cauchy")
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::anthropogenic_kernel_length_error)

  anthropogenic_kernel_type <- c("cauchy", "cauchy")
  natural_dir <- c("NONE")
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::natural_dir_length_error)

  natural_dir <- c("NONE", "NONE")
  anthropogenic_dir <- c("NONE")
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::anthropogenic_dir_length_error)

  anthropogenic_dir <- c("NONE", "NONE")
  movements_file <- c("")
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::movements_file_length_error)

  movements_file <- c("", "")
  use_movements <- c(FALSE)
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::use_movements_length_error)

  use_movements <- c(FALSE, FALSE)
  start_exposed <- c(TRUE)
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::start_exposed_length_error)

  start_exposed <- c(TRUE, TRUE)
  quarantine_areas_file <- c("")
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::quarantine_areas_length_error)

  quarantine_areas_file <- c("", "")
  use_quarantine <- c(TRUE)
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::use_quarantine_length_error)

  use_quarantine <- c(TRUE, TRUE)
  use_spreadrates <- c(TRUE)
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(
    multispecies_check$failed_check, PoPS:::use_spreadrates_length_error)

  use_spreadrates <- c(TRUE, TRUE)
  species <- c("species1")
  multispecies_check <-
    multispecies_checks(species,
                        infected_files,
                        parameter_means,
                        parameter_cov_matrix,
                        natural_kernel_type,
                        anthropogenic_kernel_type,
                        natural_dir,
                        anthropogenic_dir,
                        model_type,
                        host_file,
                        total_populations_file,
                        temp,
                        temperature_coefficient_file,
                        precip,
                        precipitation_coefficient_file,
                        latency_period,
                        time_step,
                        season_month_start,
                        season_month_end,
                        use_lethal_temperature,
                        temperature_file,
                        lethal_temperature,
                        lethal_temperature_month,
                        mortality_on,
                        mortality_rate,
                        mortality_time_lag,
                        movements_file,
                        use_movements,
                        start_exposed,
                        quarantine_areas_file,
                        use_quarantine,
                        use_spreadrates)
  expect_equal(multispecies_check$checks_passed, FALSE)
  expect_equal(multispecies_check$failed_check, PoPS:::species_length_error)
})

test_that("Multiple random seed draws and checks work", {
  randoms_file <- "NULL.fs"
  randoms <- random_seeds_file_checks(randoms_file)
  expect_equal(randoms$checks_passed, FALSE)
  expect_equal(randoms$failed_check, PoPS:::file_exists_error)
  
  randoms_file <-
    system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  randoms <- random_seeds_file_checks(randoms_file)
  expect_equal(randoms$checks_passed, FALSE)
  expect_equal(randoms$failed_check, PoPS:::file_type_error)
  
  randoms_file <-
    system.file("extdata", "simple2x2", "randoms.csv", package = "PoPS")
  randoms <- random_seeds_file_checks(randoms_file)
  expect_equal(randoms$checks_passed, TRUE)
})