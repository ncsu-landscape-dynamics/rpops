context("test-checks")

test_that("Initial raster checks report correct errors and return a raster", {
  infected_file <- ""
  infected <- initial_raster_checks(infected_file)
  expect_equal(infected$checks_passed, FALSE)
  expect_equal(infected$failed_check, "file does not exist")

  infected_file <-
    system.file("extdata", "simple2x2", "infected.csv", package = "PoPS")
  infected <- initial_raster_checks(infected_file)
  expect_equal(infected$checks_passed, FALSE)
  expect_equal(infected$failed_check,
               "file is not one of '.grd', '.tif', '.img', or '.vrt'")

  infected_file <-
    system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  infected <- initial_raster_checks(infected_file)
  expect_equal(infected$checks_passed, TRUE)
  expect_equal(infected$raster[], terra::rast(infected_file)[])
})

test_that("Initial raster checks report correct errors and return a raster", {
  extent_error <-
    "Extents of input rasters do not match. Ensure that all of your input
    rasters have the same extent"
  resolution_error <-
    "Resolution of input rasters do not match. Ensure that all of your input
    rasters have the same resolution"
  crs_error <-
    "Coordinate reference system (crs) of input rasters do not match. Ensure
    that all of your input rasters have the same crs"

  infected_file <-
    system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  infected <- terra::rast(infected_file)
  host_file <- ""
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check, "file does not exist")

  host_file <-
    system.file("extdata", "simple2x2", "infected.csv", package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check,
               "file is not one of '.grd', '.tif', '.img', or '.vrt'")

  host_file <-
    system.file("extdata", "simple5x5", "total_plants.tif", package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check,
               extent_error)

  host_file <-
    system.file("extdata", "simple2x2", "total_plants_diff_res.tif",
                package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check,
               resolution_error)

  host_file <-
    system.file("extdata", "simple2x2", "total_plants_diff_xres.tif",
                package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check,
               resolution_error)

  host_file <- system.file("extdata", "simple2x2",
                           "total_plants_diff_yres.tif", package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check,
               resolution_error)

  host_file <-
    system.file("extdata", "simple2x2", "critical_temp_diff_crs.tif",
                package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, FALSE)
  expect_equal(host$failed_check,
               crs_error)

  host_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  host <- secondary_raster_checks(host_file, infected)
  expect_equal(host$checks_passed, TRUE)
  expect_equal(host$raster[], terra::rast(host_file)[])

})

test_that(
  "Treatment checks report correct errors and return a matrices of treatments
  if all checks pass", {
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
    expect_equal(
      treatment_check$failed_check,
      "Length of list for treatment dates and treatments_file must be equal")

    treatments_file <-
      system.file("extdata", "simple2x2", "treatments.tif", package = "PoPS")
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
    expect_equal(
      treatment_check$failed_check,
      "Length of list for treatment dates and pesticide_duration must be equal")

    treatments_file <-
      system.file("extdata", "simple2x2", "treatments.tif", package = "PoPS")
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
                 terra::as.matrix(treatment_stack[[1]],
                                  wide = TRUE))

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
                 terra::as.matrix(treatment_stack[[1]],
                                  wide = TRUE))
  })


test_that("Treatment method checks report correct errors and TRUE otherwise", {
  treatment_method <- "Test"
  methods_check <- treatment_metric_checks(treatment_method)
  expect_equal(methods_check$checks_passed, FALSE)
  expect_equal(methods_check$failed_check,
               "treatment method is not one of the valid treatment options")
  treatment_method <- "ratio"
  methods_check <- treatment_metric_checks(treatment_method)
  expect_equal(methods_check$checks_passed, TRUE)
  treatment_method <- "all infected"
  methods_check <- treatment_metric_checks(treatment_method)
  expect_equal(methods_check$checks_passed, TRUE)

})

test_that(
  "Success metric checks report correct errors and return configuration
  disagreement boolean otherwise", {
    success_metric <- "Test"
    metric_check <- metric_checks(success_metric)
    expect_equal(metric_check$checks_passed, FALSE)
    expect_equal(
      metric_check$failed_check,
      "Success metric must be one of 'quantity','quantity and configuration',
    'residual error', or 'odds ratio'")
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

test_that(
  "Percent checks report correct errors and return whether or not to use the
  anthropogenic distance kernel", {
    percent_natural_dispersal <- 1.01
    percent_check <- percent_checks(percent_natural_dispersal)
    expect_equal(percent_check$checks_passed, FALSE)
    expect_equal(percent_check$failed_check,
                 "Percent natural dispersal must be between 0.0 and 1.0")
    percent_natural_dispersal <- -0.2
    percent_check <- percent_checks(percent_natural_dispersal)
    expect_equal(percent_check$checks_passed, FALSE)
    expect_equal(percent_check$failed_check,
                 "Percent natural dispersal must be between 0.0 and 1.0")
    percent_natural_dispersal <- 1.0
    percent_check <- percent_checks(percent_natural_dispersal)
    expect_equal(percent_check$checks_passed, TRUE)
    expect_equal(percent_check$use_anthropogenic_kernel, FALSE)
    percent_natural_dispersal <- 0.80
    percent_check <- percent_checks(percent_natural_dispersal)
    expect_equal(percent_check$checks_passed, TRUE)
    expect_equal(percent_check$use_anthropogenic_kernel, TRUE)

  })

test_that(
  "Time checks report correct errors and return number of time steps, number of
  years, and number of outputs otherwise", {
    output_frequency_error <-
      "Output frequency is more frequent than time_step. The minimum
      output_frequency you can use is the time_step of your simulation. You can
      set the output_frequency to 'time_step' to default to most frequent
      output possible"
    time_format_error <-
    "End time and/or start time not of type numeric and/or in format YYYY-MM-DD"
    end_date <- 2017
    start_date <- 2016
    time_step <- "hour"
    output_frequency <- "week"
    time_check <- time_checks(end_date, start_date, time_step, output_frequency)
    expect_equal(time_check$checks_passed, FALSE)
    expect_equal(time_check$failed_check,
                 "Time step must be one of 'week', 'month' or 'day'")

    end_date <- 2017
    start_date <- 2016
    time_step <- "week"
    output_frequency <- "week"
    time_check <- time_checks(end_date, start_date, time_step, output_frequency)
    expect_equal(time_check$checks_passed, FALSE)
    expect_equal(
      time_check$failed_check,
      time_format_error)

    end_date <- "2017-01-01"
    start_date <- "2016-01-01"
    time_step <- "month"
    output_frequency <- "hour"
    time_check <- time_checks(end_date, start_date, time_step, output_frequency)
    expect_equal(time_check$checks_passed, FALSE)
    expect_equal(
      time_check$failed_check,
      "Output frequency must be either 'week', 'month', 'day', 'year', 'time_step', or 'every_n_steps'")

    end_date <- "2017-01-01"
    start_date <- "2016-01-01"
    time_step <- "month"
    output_frequency <- "week"
    time_check <- time_checks(end_date, start_date, time_step, output_frequency)
    expect_equal(time_check$checks_passed, FALSE)
    expect_equal(
      time_check$failed_check,
      output_frequency_error)

    end_date <- "2017-01-01"
    start_date <- "2016-01-01"
    time_step <- "week"
    output_frequency <- "week"
    time_check <- time_checks(end_date, start_date, time_step, output_frequency)
    expect_equal(time_check$checks_passed, TRUE)
    expect_equal(time_check$number_of_time_steps, 53)
    expect_equal(time_check$number_of_years, 1)
    expect_equal(time_check$number_of_outputs, 53)

    end_date <- "2016-02-01"
    start_date <- "2016-01-01"
    time_step <- "week"
    output_frequency <- "week"
    time_check <- time_checks(end_date, start_date, time_step, output_frequency)
    expect_equal(time_check$checks_passed, TRUE)
    expect_equal(time_check$number_of_time_steps, 5)
    expect_equal(time_check$number_of_years, 1)
    expect_equal(time_check$number_of_outputs, 5)

    end_date <- "2017-01-01"
    start_date <- "2016-01-01"
    time_step <- "week"
    output_frequency <- "month"
    time_check <- time_checks(end_date, start_date, time_step, output_frequency)
    expect_equal(time_check$checks_passed, TRUE)
    expect_equal(time_check$number_of_time_steps, 53)
    expect_equal(time_check$number_of_years, 1)
    expect_equal(time_check$number_of_outputs, 12)

    end_date <- "2017-01-01"
    start_date <- "2016-01-01"
    time_step <- "week"
    output_frequency <- "year"
    time_check <- time_checks(end_date, start_date, time_step, output_frequency)
    expect_equal(time_check$checks_passed, TRUE)
    expect_equal(time_check$number_of_time_steps, 53)
    expect_equal(time_check$number_of_years, 1)
    expect_equal(time_check$number_of_outputs, 1)

  })

test_that("Prior checks report correct errors and return reformatted priors,
          and mean and sd of priors for starting location for MCMC", {
            skip_on_appveyor()
            priors <- 2.1
            prior_check <- prior_checks(priors)
            expect_equal(prior_check$checks_passed, FALSE)
            expect_equal(prior_check$failed_check,
                         "Incorrect format for priors")

            priors <-
              data.frame(data = c(1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1),
                         per = c(0.1, 0.12, 0.13, 0.30, 0.13, 0.12, 0.1))
            prior_check <- prior_checks(priors)
            names(priors) <- c("var", "prob")
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

test_that(
  "Bayesian checks report correct errors and return all rates as a data frame
  and posterior rates as a data frame", {
    skip_on_appveyor()
    priors <- 2.1
    prior_check <- prior_checks(priors)
    prior <- 2.1
    start_priors <- 2.1
    sd_priors <- 0.2
    params <-
      data.frame(data = c(rep(1.5, 10),
                          rep(1.5, 10),
                          rep(1.6, 15),
                          rep(1.7, 15),
                          rep(1.8, 30),
                          rep(1.9, 15),
                          rep(2.0, 15),
                          rep(2.1, 15),
                          rep(2.2, 15)))
    count <- 10000000
    prior_weight <- 0.05
    weight <- 1 - prior_weight
    step_size <- 0.1
    bounds <- c(0, Inf)
    round_to <- 1
    round_to_digits <- 1
    bayesian_check <-
      bayesian_checks(prior,
                      start_priors,
                      sd_priors,
                      params,
                      count,
                      prior_weight,
                      weight,
                      step_size,
                      bounds = c(0, Inf),
                      round_to = 1,
                      round_to_digits = 1)
    expect_equal(bayesian_check$checks_passed, TRUE)
    expect_length(bayesian_check$rates, 4)
    expect_length(bayesian_check$posterior_rates, 2)
    expect_lte(weighted.mean(bayesian_check$posterior_rates[, 1],
                             bayesian_check$posterior_rates[, 2]),
               2.1)
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
  expect_equal(mnn_check$posterior_cov_matrix,
               matrix(c(1, 1.95, 1, 1), ncol = 2, nrow = 2))
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
  expect_equal(multispecies_check$failed_check,
               "Length of list for species and infected_files must be equal")

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
  expect_equal(multispecies_check$failed_check,
               "Length of list for infected_files and host_file must be equal")

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
  expect_equal(
    multispecies_check$failed_check,
    "Length of list for infected_files and total_populations_file must be equal")

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
    multispecies_check$failed_check,
    "Length of list for parameter_means and infected_files must be equal")

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
    multispecies_check$failed_check,
    "Length of list for parameter_cov_matrix and infected_files must be equal")

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
  expect_equal(multispecies_check$failed_check,
               "Length of list for infected_files and temp must be equal")

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
  expect_equal(
    multispecies_check$failed_check,
    "Length of list for infected_files and temperature_coefficient_file must be equal")

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
  expect_equal(multispecies_check$failed_check,
               "Length of list for infected_files and precip must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and precipitation_coefficient_file must be equal")

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
  expect_equal(multispecies_check$failed_check,
               "Length of list for model_type and infected_files must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and latency_period must be equal")

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
  expect_equal(multispecies_check$failed_check,
               "Length of list for infected_files and time_step must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and season_month_start must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and season_month_end must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and use_lethal_temperature must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and temperature_file must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and lethal_temperature must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and lethal_temperature_month must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and mortality_on must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and mortality_rate must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and mortality_time_lag must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and natural_kernel_type must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and anthropogenic_kernel_type must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and natural_dir must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and anthropogenic_dir must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and movements_file must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and use_movements must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and start_exposed must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and quarantine_areas_file must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and use_quarantine must be equal")

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
    multispecies_check$failed_check,
    "Length of list for infected_files and use_spreadrates must be equal")

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
  expect_equal(multispecies_check$failed_check,
               "Length of list for species and infected_files must be equal")

})
