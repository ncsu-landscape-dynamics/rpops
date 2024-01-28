context("test-pops")

test_that("Model stops if files don't exist or aren't the correct extension", {
  infected_file_list <- c("")
  host_file_list <- c(system.file("extdata", "simple2x2", "host.tif", package = "PoPS"))
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file <-
    system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file <-
    system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2010-12-31"
  parameter_means <- c(0, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  pest_host_table <- system.file("extdata", "pest_host_table_singlehost.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table))

  infected_file_list <- system.file("extdata", "simple2x2", "infected.csv", package = "PoPS")
  expect_error(pops(
    infected_file_list =  infected_file_list,
    host_file_list =  host_file_list,
    total_populations_file =  host_file_list,
    parameter_means = parameter_means,
    parameter_cov_matrix = parameter_cov_matrix,
    pest_host_table = pest_host_table,
    competency_table = competency_table),
    raster_type_error, fixed = TRUE)

  host_file_list <- ""
  infected_file_list <- c(system.file("extdata", "simple2x2", "infected.tif", package = "PoPS"))
  expect_error(pops(infected_file_list =  infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file =  total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table))

  host_file_list <- c(system.file("extdata", "simple2x2", "infected.csv", package = "PoPS"))
  expect_error(pops(infected_file_list =  infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file =  host_file_list,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               raster_type_error, fixed = TRUE)

  host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
  total_populations_file <- ""
  expect_error(pops(infected_file_list =  infected_file_list,
                    host_file_list =  host_file_list,
                    total_populations_file =  total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               detailed_file_exists_error(total_populations_file))

  expect_error(pops(infected_file_list =  infected_file_list,
                    host_file_list =  host_file_list,
                    total_populations_file =
                      system.file("extdata", "simple2x2", "infected.csv", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               raster_type_error, fixed = TRUE)

  temperature_file <- ""
  expect_error(pops(infected_file_list =  infected_file_list,
                    host_file_list =  host_file_list,
                    total_populations_file =  host_file_list,
                    use_lethal_temperature = TRUE,
                    temperature_file = temperature_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               detailed_file_exists_error(temperature_file))

  expect_error(pops(infected_file_list =  infected_file_list,
                    host_file_list =  host_file_list,
                    total_populations_file =  host_file_list,
                    use_lethal_temperature = TRUE,
                    temperature_file =
                      system.file("extdata", "simple2x2", "infected.csv", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               raster_type_error, fixed = TRUE)

  temperature_file <-
    system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  temperature_coefficient_file <- ""
  expect_error(pops(infected_file_list =  infected_file_list,
                    host_file_list =  host_file_list,
                    total_populations_file =  host_file_list,
                    temp = TRUE,
                    temperature_coefficient_file = temperature_coefficient_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               detailed_file_exists_error(temperature_coefficient_file))

  expect_error(pops(infected_file_list =  infected_file_list,
                    host_file_list =  host_file_list,
                    total_populations_file =  host_file_list,
                    temp = TRUE,
                    temperature_coefficient_file =
                      system.file("extdata", "simple2x2", "infected.csv", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               raster_type_error, fixed = TRUE)

  precipitation_coefficient_file <- ""
  expect_error(pops(infected_file_list =  infected_file_list,
                    host_file_list =  host_file_list,
                    total_populations_file =  host_file_list,
                    precip = TRUE,
                    precipitation_coefficient_file = "",
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               detailed_file_exists_error(precipitation_coefficient_file))

  expect_error(pops(infected_file_list =  infected_file_list,
                    host_file_list =  host_file_list,
                    total_populations_file =  host_file_list,
                    precip = TRUE,
                    precipitation_coefficient_file =
                      system.file("extdata", "simple2x2", "infected.csv", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               raster_type_error, fixed = TRUE)

  treatments_file <- ""
  expect_error(pops(infected_file_list =  infected_file_list,
                    host_file_list =  host_file_list,
                    total_populations_file =  host_file_list,
                    management = TRUE,
                    treatments_file = treatments_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               detailed_file_exists_error(treatments_file))

  expect_error(pops(infected_file_list =  infected_file_list,
                    host_file_list =  host_file_list,
                    total_populations_file =  host_file_list,
                    management = TRUE,
                    treatments_file =
                      system.file("extdata", "simple2x2", "infected.csv", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               raster_type_error, fixed = TRUE)

  exposed_file_list <- ""
  expect_error(pops(infected_file_list =  infected_file_list,
                    host_file_list =  host_file_list,
                    total_populations_file =  host_file_list,
                    model_type = "SEI",
                    exposed_file_list = exposed_file_list,
                    latency_period = 2,
                    start_exposed = TRUE,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table))

  quarantine_areas_file <- ""
  expect_error(pops(infected_file_list =  infected_file_list,
                    host_file_list =  host_file_list,
                    total_populations_file =  host_file_list,
                    quarantine_areas_file = quarantine_areas_file,
                    use_quarantine = TRUE,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               detailed_file_exists_error(quarantine_areas_file))

  expect_error(pops(infected_file_list =  infected_file_list,
                    host_file_list =  host_file_list,
                    total_populations_file =  host_file_list,
                    parameter_means = c(0, 0, 0, 0),
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               paramter_means_error, fixed = TRUE)

  expect_error(pops(infected_file_list =  infected_file_list,
                    host_file_list =  host_file_list,
                    total_populations_file =  host_file_list,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = matrix(0, nrow = 5, ncol = 6),
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               covariance_mat_error, fixed = TRUE)
})

test_that("Model stops if treatments don't have correct dimenisions", {
  infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2010-12-31"
  parameter_means <- c(0, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  treatments_file <- system.file("extdata", "simple2x2", "treatments.tif", package = "PoPS")
  treatment_dates <- c("2008-01-01", "2008-05-01")
  pest_host_table <- system.file("extdata", "pest_host_table_singlehost.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list =  host_file_list,
                    total_populations_file =  host_file_list,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table,
                    management = TRUE,
                    treatments_file = treatments_file,
                    treatment_dates = treatment_dates),
               treatment_length_error, fixed = TRUE)
})

test_that("Model stops if time and date parameters are of the wrong type and/or dimension", {
  infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file <-
    system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file <-
    system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2010-12-31"
  parameter_means <- c(0, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  pest_host_table <- system.file("extdata", "pest_host_table_singlehost.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    time_step = "two",
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               time_step_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    end_date = "two",
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               date_format_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    end_date = 156,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               date_format_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    start_date = "five",
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               date_format_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    start_date = 19,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               date_format_error, fixed = TRUE)
})

test_that("Model stops if kernel is of the wrong type and/or dimension", {
  infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file <-
    system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file <- system.file("extdata", "simple2x2",
                                  "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2010-12-31"
  parameter_means <- c(0, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  pest_host_table <- system.file("extdata", "pest_host_table_singlehost.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  expect_error(
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         natural_kernel_type = "none",
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table),
    natural_kernel_error, fixed = TRUE)
  expect_error(
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         anthropogenic_kernel_type = "none",
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table),
    anthropogenic_kernel_error, fixed = TRUE)
})

test_that("Input raster resolutions, extents, and crs all match", {
  infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file <-
    system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file <-
    system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2010-12-31"
  parameter_means <- c(0, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  pest_host_table <- system.file("extdata", "pest_host_table_singlehost.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list =
                      system.file("extdata", "simple5x5", "total_plants.tif", package = "PoPS"),
                    total_populations_file = total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               extent_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file =
                      system.file("extdata", "simple5x5", "total_plants.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               extent_error, fixed = TRUE)
  expect_error(pops(infected_file_list =
                      system.file("extdata", "simple5x5", "total_plants.tif", package = "PoPS"),
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               extent_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    use_lethal_temperature = TRUE,
                    temperature_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_extent.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               extent_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    temp = TRUE,
                    temperature_coefficient_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_extent.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               extent_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    precip = TRUE,
                    precipitation_coefficient_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_extent.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               extent_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    temp = TRUE,
                    temperature_coefficient_file =
                      system.file("extdata", "simple2x2", "critical_temp.tif", package = "PoPS"),
                    precip = TRUE,
                    precipitation_coefficient_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_extent.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               extent_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    management = TRUE,
                    treatments_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_extent.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               extent_error, fixed = TRUE)

  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list =
                      system.file("extdata", "simple2x2",
                                  "total_plants_diff_res.tif", package = "PoPS"),
                    total_populations_file = total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list =
                      system.file("extdata", "simple2x2",
                                  "total_plants_diff_xres.tif", package = "PoPS"),
                    total_populations_file = total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list =
                      system.file("extdata", "simple2x2",
                                  "total_plants_diff_yres.tif", package = "PoPS"),
                    total_populations_file = total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file =
                      system.file("extdata", "simple2x2",
                                  "total_plants_diff_res.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file =
                      system.file("extdata", "simple2x2",
                                  "total_plants_diff_xres.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file =
                      system.file("extdata", "simple2x2",
                                  "total_plants_diff_yres.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list =
                      system.file("extdata", "simple2x2",
                                  "total_plants_diff_res.tif", package = "PoPS"),
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list =
                      system.file("extdata", "simple2x2",
                                  "total_plants_diff_xres.tif", package = "PoPS"),
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list =
                      system.file("extdata", "simple2x2",
                                  "total_plants_diff_yres.tif", package = "PoPS"),
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    use_lethal_temperature = TRUE,
                    temperature_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_res.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    use_lethal_temperature = TRUE,
                    temperature_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_xres.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    use_lethal_temperature = TRUE,
                    temperature_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_yres.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    temp = TRUE,
                    temperature_coefficient_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_res.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    temp  = TRUE,
                    temperature_coefficient_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_xres.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    temp = TRUE,
                    temperature_coefficient_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_yres.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    precip = TRUE,
                    precipitation_coefficient_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_res.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    precip  = TRUE,
                    precipitation_coefficient_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_xres.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    precip = TRUE,
                    precipitation_coefficient_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_yres.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    temp = TRUE,
                    temperature_coefficient_file =
                      system.file("extdata", "simple2x2", "critical_temp.tif", package = "PoPS"),
                    precip = TRUE,
                    precipitation_coefficient_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_res.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    temp = TRUE,
                    temperature_coefficient_file =
                      system.file("extdata", "simple2x2", "critical_temp.tif", package = "PoPS"),
                    precip = TRUE,
                    precipitation_coefficient_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_xres.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    temp = TRUE,
                    temperature_coefficient_file =
                      system.file("extdata", "simple2x2", "critical_temp.tif", package = "PoPS"),
                    precip = TRUE,
                    precipitation_coefficient_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_yres.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    management = TRUE,
                    treatments_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_res.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    management = TRUE,
                    treatments_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_xres.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    management = TRUE,
                    treatments_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_yres.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               resolution_error, fixed = TRUE)

  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_crs.tif", package = "PoPS"),
                    total_populations_file = total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               crs_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_crs.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               crs_error, fixed = TRUE)
  expect_error(pops(infected_file_list =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_crs.tif", package = "PoPS"),
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               crs_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    use_lethal_temperature = TRUE,
                    temperature_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_crs.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               crs_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    temp = TRUE,
                    temperature_coefficient_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_crs.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               crs_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    precip = TRUE,
                    precipitation_coefficient_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_crs.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               crs_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    temp = TRUE,
                    temperature_coefficient_file =
                      system.file("extdata", "simple2x2", "critical_temp.tif",
                                  package = "PoPS"),
                    precip = TRUE,
                    precipitation_coefficient_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_crs.tif",
                                  package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               crs_error, fixed = TRUE)
  expect_error(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    management = TRUE,
                    treatments_file =
                      system.file("extdata", "simple2x2",
                                  "critical_temp_diff_crs.tif",
                                  package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table),
               crs_error, fixed = TRUE)

})

test_that("Infected results return initial infected if reproductive rate is set to 0", {
  infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file <-
    system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file <-
    system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2010-12-31"
  parameter_means <- c(0, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  expect_equal(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table)$host_pools[[1]]$infected[[1]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  expect_equal(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    use_lethal_temperature = TRUE,
                    lethal_temperature = -16,
                    lethal_temperature_month = 1,
                    temperature_file =
                      system.file("extdata", "simple2x2", "critical_temp.tif", package = "PoPS"),
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table)$host_pools[[1]]$infected[[1]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  expect_equal(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table,
                    temp = TRUE,
                    temperature_coefficient_file = coefficient_file)$host_pools[[1]]$infected[[1]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  expect_equal(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table,
                    precip = TRUE,
                    precipitation_coefficient_file =
                      coefficient_file)$host_pools[[1]]$infected[[1]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  expect_equal(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table,
                    temp = TRUE,
                    temperature_coefficient_file = coefficient_file,
                    precip = TRUE,
                    precipitation_coefficient_file =
                      coefficient_file)$host_pools[[1]]$infected[[1]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  expect_equal(pops(infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table,
                    use_lethal_temperature = TRUE,
                    lethal_temperature = -16,
                    lethal_temperature_month = 1,
                    temperature_file =
                      system.file("extdata", "simple2x2", "critical_temp.tif", package = "PoPS"),
                    temp = TRUE,
                    temperature_coefficient_file = coefficient_file,
                    precip = TRUE,
                    precipitation_coefficient_file =
                      coefficient_file)$host_pools[[1]]$infected[[1]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))

  expect_equal(
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         temp = TRUE,
         temperature_coefficient_file =
           system.file("extdata", "simple2x2",
                       "temperature_coefficient_weeks.tif", package = "PoPS"),
         time_step = "week")$host_pools[[1]]$infected[[1]],
    terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  expect_equal(
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         precip = TRUE,
         precipitation_coefficient_file =
           system.file("extdata", "simple2x2",
                       "temperature_coefficient_weeks.tif", package = "PoPS"),
         time_step = "week")$host_pools[[1]]$infected[[1]],
    terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  expect_equal(
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         temp = TRUE,
         temperature_coefficient_file =
           system.file("extdata", "simple2x2",
                       "temperature_coefficient_weeks.tif", package = "PoPS"),
         precip = TRUE,
         precipitation_coefficient_file =
           system.file("extdata", "simple2x2",
                       "temperature_coefficient_weeks.tif", package = "PoPS"),
         time_step = "week")$host_pools[[1]]$infected[[1]],
    terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  expect_equal(
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         use_lethal_temperature = TRUE,
         lethal_temperature = -16,
         lethal_temperature_month = 1,
         temperature_file =
           system.file("extdata", "simple2x2", "critical_temp.tif",
                       package = "PoPS"),
         temp = TRUE,
         temperature_coefficient_file =
           system.file("extdata", "simple2x2",
                       "temperature_coefficient_weeks.tif", package = "PoPS"),
         precip = TRUE,
         precipitation_coefficient_file =
           system.file("extdata", "simple2x2",
                       "temperature_coefficient_weeks.tif", package = "PoPS"),
         time_step = "week")$host_pools[[1]]$infected[[1]],
    terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))

  expect_equal(
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         temp = TRUE,
         time_step = "day",
         temperature_coefficient_file =
           system.file("extdata", "simple2x2", "temperature_coefficient_days.tif",
                       package = "PoPS"))$host_pools[[1]]$infected[[1]],
    terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  expect_equal(
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         precip = TRUE,
         precipitation_coefficient_file =
           system.file("extdata", "simple2x2",
                       "temperature_coefficient_days.tif", package = "PoPS"),
         time_step = "day")$host_pools[[1]]$infected[[1]],
    terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  expect_equal(
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         temp = TRUE,
         temperature_coefficient_file =
           system.file("extdata", "simple2x2",
                       "temperature_coefficient_days.tif", package = "PoPS"),
         precip = TRUE,
         precipitation_coefficient_file =
           system.file("extdata", "simple2x2",
                       "temperature_coefficient_days.tif", package = "PoPS"),
         time_step = "day")$host_pools[[1]]$infected[[1]],
    terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  expect_equal(
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         use_lethal_temperature = TRUE,
         lethal_temperature = -16,
         lethal_temperature_month = 1,
         temperature_file =
           system.file("extdata", "simple2x2", "critical_temp.tif", package = "PoPS"),
         temp = TRUE,
         temperature_coefficient_file =
           system.file("extdata", "simple2x2",
                       "temperature_coefficient_days.tif", package = "PoPS"),
         precip = TRUE,
         precipitation_coefficient_file =
           system.file("extdata", "simple2x2",
                       "temperature_coefficient_days.tif", package = "PoPS"),
         time_step = "day")$host_pools[[1]]$infected[[1]],
    terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
})

test_that(
  "Infected results returns all 0's if minimum temp drops below lethal temperature", {
    infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
    host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
    total_populations_file <-
      system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
    coefficient_file <-
      system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
    temperature_file <-
      system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
    start_date <- "2008-01-01"
    end_date <- "2010-12-31"
    parameter_means <- c(1, 21, 1, 500, 0, 0, 0, 0)
    parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
    pest_host_table <-
      system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
    competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

    expect_equal(pops(infected_file_list = infected_file_list,
                      host_file_list = host_file_list,
                      total_populations_file = total_populations_file,
                      use_lethal_temperature = TRUE,
                      temperature_file = temperature_file,
                      parameter_means = parameter_means,
                      parameter_cov_matrix = parameter_cov_matrix,
                      pest_host_table = pest_host_table,
                      competency_table = competency_table)$host_pools[[1]]$infected[[1]],
                 matrix(0, ncol = 2, nrow = 2))
    expect_equal(pops(infected_file_list = infected_file_list,
                      host_file_list = host_file_list,
                      total_populations_file = total_populations_file,
                      use_lethal_temperature = TRUE,
                      temperature_file = temperature_file,
                      precip = TRUE,
                      precipitation_coefficient_file = coefficient_file,
                      parameter_means = parameter_means,
                      parameter_cov_matrix = parameter_cov_matrix,
                      pest_host_table = pest_host_table,
                      competency_table = competency_table)$host_pools[[1]]$infected[[1]],
                 matrix(0, ncol = 2, nrow = 2))
    expect_equal(pops(infected_file_list = infected_file_list,
                      host_file_list = host_file_list,
                      total_populations_file = total_populations_file,
                      use_lethal_temperature = TRUE,
                      temperature_file = temperature_file,
                      temp = TRUE,
                      temperature_coefficient_file = coefficient_file,
                      parameter_means = parameter_means,
                      parameter_cov_matrix = parameter_cov_matrix,
                      pest_host_table = pest_host_table,
                      competency_table = competency_table)$host_pools[[1]]$infected[[1]],
                 matrix(0, ncol = 2, nrow = 2))
    expect_equal(pops(infected_file_list = infected_file_list,
                      host_file_list = host_file_list,
                      total_populations_file = total_populations_file,
                      use_lethal_temperature = TRUE,
                      temperature_file = temperature_file,
                      temp = TRUE,
                      temperature_coefficient_file = coefficient_file,
                      precip = TRUE,
                      precipitation_coefficient_file = coefficient_file,
                      parameter_means = parameter_means,
                      parameter_cov_matrix = parameter_cov_matrix,
                      pest_host_table = pest_host_table,
                      competency_table = competency_table)$host_pools[[1]]$infected[[1]],
                 matrix(0, ncol = 2, nrow = 2))

  })

test_that(
  "Infected results returns less infection after survival rates than before", {
    infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
    host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
    total_populations_file <-
      system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
    coefficient_file <-
      system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
    survival_rates_file <-
      system.file("extdata", "simple2x2", "survival_rates.tif", package = "PoPS")
    start_date <- "2008-01-01"
    end_date <- "2010-12-31"
    parameter_means <- c(0, 21, 1, 500, 0, 0, 0, 0)
    parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
    pest_host_table <-
      system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
    competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

    reduced_inf <- matrix(0, ncol = 2, nrow = 2)
    reduced_inf[1, 1] <- 3

    expect_equal(pops(infected_file_list = infected_file_list,
                      host_file_list = host_file_list,
                      total_populations_file = total_populations_file,
                      use_survival_rates = TRUE,
                      survival_rates_file = survival_rates_file,
                      parameter_means = parameter_means,
                      parameter_cov_matrix = parameter_cov_matrix,
                      pest_host_table = pest_host_table,
                      competency_table = competency_table)$host_pools[[1]]$infected[[1]],
                 reduced_inf)
    expect_equal(pops(infected_file_list = infected_file_list,
                      host_file_list = host_file_list,
                      total_populations_file = total_populations_file,
                      use_survival_rates = TRUE,
                      survival_rates_file = survival_rates_file,
                      precip = TRUE,
                      precipitation_coefficient_file = coefficient_file,
                      parameter_means = parameter_means,
                      parameter_cov_matrix = parameter_cov_matrix,
                      pest_host_table = pest_host_table,
                      competency_table = competency_table)$host_pools[[1]]$infected[[1]],
                 reduced_inf)
    expect_equal(pops(infected_file_list = infected_file_list,
                      host_file_list = host_file_list,
                      total_populations_file = total_populations_file,
                      use_survival_rates = TRUE,
                      survival_rates_file = survival_rates_file,
                      temp = TRUE,
                      temperature_coefficient_file = coefficient_file,
                      parameter_means = parameter_means,
                      parameter_cov_matrix = parameter_cov_matrix,
                      pest_host_table = pest_host_table,
                      competency_table = competency_table)$host_pools[[1]]$infected[[1]],
                 reduced_inf)
    expect_equal(pops(infected_file_list = infected_file_list,
                      host_file_list = host_file_list,
                      total_populations_file = total_populations_file,
                      use_survival_rates = TRUE,
                      survival_rates_file = survival_rates_file,
                      temp = TRUE,
                      temperature_coefficient_file = coefficient_file,
                      precip = TRUE,
                      precipitation_coefficient_file = coefficient_file,
                      parameter_means = parameter_means,
                      parameter_cov_matrix = parameter_cov_matrix,
                      pest_host_table = pest_host_table,
                      competency_table = competency_table)$host_pools[[1]]$infected[[1]],
                 reduced_inf)

  })

test_that("Infected and Susceptible results return all 0's if treatments file is all 1's but
          leaves a proportion of susceptibles if treatment method is ratio", {
            infected_file_list <-
              system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
            host_file_list <-
              system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
            total_populations_file <-
              system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
            coefficient_file <-
              system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
            temperature_file <-
              system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif",
                          package = "PoPS")
            start_date <- "2008-01-01"
            end_date <- "2010-12-31"
            treatments_file <- system.file("extdata", "simple2x2", "treatments.tif",
                                           package = "PoPS")
            parameter_means <- c(0, 21, 1, 500, 0, 0, 0, 0)
            parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
            pest_host_table <-
              system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
            competency_table <-
              system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

            data <-
              pops(infected_file_list = infected_file_list,
                   host_file_list = host_file_list,
                   total_populations_file = total_populations_file,
                   management  = TRUE,
                   treatment_dates = c("2008-12-01"),
                   treatments_file = treatments_file,
                   parameter_means = parameter_means,
                   parameter_cov_matrix = parameter_cov_matrix,
                   pest_host_table = pest_host_table,
                   competency_table = competency_table,
                   start_date = start_date,
                   end_date = end_date)

            expect_equal(data$host_pools[[1]]$infected[[1]], matrix(0, ncol = 2, nrow = 2))
            expect_equal(data$host_pools[[1]]$susceptible[[1]], matrix(0, ncol = 2, nrow = 2))

            data <-
              pops(infected_file_list = infected_file_list,
                   host_file_list = host_file_list,
                   treatment_method = "all infected",
                   total_populations_file = total_populations_file,
                   management  = TRUE,
                   treatment_dates = c("2008-12-01"),
                   treatments_file = treatments_file,
                   parameter_means = parameter_means,
                   parameter_cov_matrix = parameter_cov_matrix,
                   pest_host_table = pest_host_table,
                   competency_table = competency_table,
                   start_date = start_date,
                   end_date = end_date)

            expect_equal(data$host_pools[[1]]$infected[[1]], matrix(0, ncol = 2, nrow = 2))
            expect_equal(data$host_pools[[1]]$susceptible[[1]], matrix(0, ncol = 2, nrow = 2))

            treatments_file <-
              system.file("extdata", "simple2x2", "treatmentshalf.tif", package = "PoPS")

            data <-
              pops(infected_file_list = infected_file_list,
                   host_file_list = host_file_list,
                   treatment_method = "ratio",
                   total_populations_file = total_populations_file,
                   management  = TRUE,
                   treatment_dates = c("2008-12-01"),
                   treatments_file = treatments_file,
                   parameter_means = parameter_means,
                   parameter_cov_matrix = parameter_cov_matrix,
                   pest_host_table = pest_host_table,
                   competency_table = competency_table,
                   start_date = start_date,
                   end_date = end_date)

            expect_equal(data$host_pools[[1]]$infected[[1]],
                         matrix(c(2, 0, 0, 0), ncol = 2, nrow = 2))
            expect_equal(data$host_pools[[1]]$susceptible[[1]],
                         matrix(c(6, 7, 3, 7), ncol = 2, nrow = 2))

            data <-
              pops(infected_file_list = infected_file_list,
                   host_file_list = host_file_list,
                   treatment_method = "all infected",
                   total_populations_file = total_populations_file,
                   management  = TRUE,
                   treatment_dates = c("2008-12-01"),
                   treatments_file = treatments_file,
                   parameter_means = parameter_means,
                   parameter_cov_matrix = parameter_cov_matrix,
                   pest_host_table = pest_host_table,
                   competency_table = competency_table,
                   start_date = start_date,
                   end_date = end_date)

            expect_equal(data$host_pools[[1]]$infected[[1]],
                         matrix(c(0, 0, 0, 0), ncol = 2, nrow = 2))
            expect_equal(data$host_pools[[1]]$susceptible[[1]],
                         matrix(c(6, 7, 3, 7), ncol = 2, nrow = 2))
          })

test_that("Infected results are greater than initial infected", {
  infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file <-
    system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file <-
    system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2010-12-31"
  parameter_means <- c(1, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  expect_equal(all(pops(infected_file_list = infected_file_list,
                        host_file_list = host_file_list,
                        total_populations_file = total_populations_file,
                        parameter_means = parameter_means,
                        parameter_cov_matrix = parameter_cov_matrix,
                        pest_host_table = pest_host_table,
                        competency_table = competency_table
  )$host_pools[[1]]$infected[[1]] >=
    terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)
  expect_equal(all(
    pops(infected_file_list = infected_file_list,
         host_file_list =
           system.file("extdata", "simple2x2",
                       "total_plants_host_greater_than_infected.tif", package = "PoPS"),
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table)$host_pools[[1]]$infected[[1]] >=
      terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)

})

test_that("All kernel types lead to spread", {
  infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file <-
    system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file <-
    system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2008-12-31"
  time_step <- "month"
  parameter_means <- c(3.0, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  data <- pops(infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               time_step = time_step,
               natural_kernel_type = "exponential")

  infecteds <- data$host_pools[[1]]$infected[[1]]
  expect_equal(all(infecteds >=
                     terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)
  expect_gt(infecteds[1, 2] + infecteds[2, 1] + infecteds[2, 2], 0)

  data <- pops(infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               natural_kernel_type = "cauchy")
  infecteds <- data$host_pools[[1]]$infected[[1]]
  expect_equal(all(infecteds >=
                     terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)
  expect_gt(infecteds[1, 2] + infecteds[2, 1] + infecteds[2, 2], 0)

  data <- pops(infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               natural_kernel_type = "uniform")
  infecteds <- data$host_pools[[1]]$infected[[1]]
  expect_equal(all(infecteds >=
                     terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)
  expect_gt(infecteds[1, 2] + infecteds[2, 1] + infecteds[2, 2], 0)

  data <- pops(infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               natural_kernel_type = "hyperbolic secant")
  infecteds <- data$host_pools[[1]]$infected[[1]]
  expect_equal(all(infecteds >=
                     terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)
  expect_gt(infecteds[1, 2] + infecteds[2, 1] + infecteds[2, 2], 0)

  data <- pops(infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               natural_kernel_type = "weibull")
  infecteds <- data$host_pools[[1]]$infected[[1]]
  expect_equal(all(infecteds >=
                     terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)
  expect_gt(infecteds[1, 2] + infecteds[2, 1] + infecteds[2, 2], 0)

  data <- pops(infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               natural_kernel_type = "logistic")
  infecteds <- data$host_pools[[1]]$infected[[1]]
  expect_equal(all(infecteds >=
                     terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)
  expect_gt(infecteds[1, 2] + infecteds[2, 1] + infecteds[2, 2], 0)

  data <- pops(infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               natural_kernel_type = "gamma")
  infecteds <- data$host_pools[[1]]$infected[[1]]
  expect_equal(all(infecteds >=
                     terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)
  expect_gt(infecteds[1, 2] + infecteds[2, 1] + infecteds[2, 2], 0)

  parameter_means <- c(0.4, 2, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)

  data <- pops(infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               natural_kernel_type = "power law")
  infecteds <- data$host_pools[[1]]$infected[[1]]
  expect_equal(all(infecteds >=
                     terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)
  expect_gte(infecteds[1, 2] + infecteds[2, 1] + infecteds[2, 2], 0)

  ## currently not working
  # doesn't disperse outside of originally infected cell
  # parameter_means <- c(0.4, 1000, 1, 500, 0, 0, 0, 0)
  # parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)

  # data <- pops(infected_file_list = infected_file_list,
  #              host_file_list = host_file_list,
  #              total_populations_file = total_populations_file,
  #              parameter_means = parameter_means,
  #              parameter_cov_matrix = parameter_cov_matrix,
  #              pest_host_table = pest_host_table,
  #              competency_table = competency_table,
  #              natural_kernel_type = "exponential-power")
  # infecteds <- data$host_pools[[1]]$infected[[1]]
  # expect_equal(all(infecteds >=
  #                    terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)),
  #              TRUE)
  # expect_gt(infecteds[1,2] + infecteds[2,1] + infecteds[2,2], 0)
  #
  # # bad array when icdf is available
  # data <- pops(infected_file_list = infected_file_list,
  #              host_file_list = host_file_list,
  #              total_populations_file = total_populations_file,
  #              parameter_means = parameter_means,
  #              parameter_cov_matrix = parameter_cov_matrix,
  #              pest_host_table = pest_host_table,
  #              competency_table = competency_table,
  #              natural_kernel_type = "log normal")
  # infecteds <- data$host_pools[[1]]$infected[[1]]
  # expect_equal(all(infecteds >=
  #                    terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)),
  #              TRUE)
  # expect_gt(infecteds[1,2] + infecteds[2,1] + infecteds[2,2], 0)

  # checks for anthropogenic kernel type
  data <- pops(infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               anthropogenic_kernel_type = "exponential")

  expect_equal(all(data$host_pools[[1]]$infected[[1]] >=
                     terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)

  data <- pops(infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               anthropogenic_kernel_type = "cauchy")
  expect_equal(all(data$host_pools[[1]]$infected[[1]] >=
                     terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)

  data <- pops(infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               anthropogenic_kernel_type = "uniform")
  expect_equal(all(data$host_pools[[1]]$infected[[1]] >=
                     terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)

  data <- pops(infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               anthropogenic_kernel_type = "hyperbolic secant")
  expect_equal(all(data$host_pools[[1]]$infected[[1]] >=
                     terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)

  data <- pops(infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               anthropogenic_kernel_type = "logistic")
  expect_equal(all(data$host_pools[[1]]$infected[[1]] >=
                     terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)

  data <- pops(infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               anthropogenic_kernel_type = "weibull")
  expect_equal(all(data$host_pools[[1]]$infected[[1]] >=
                     terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)

  data <- pops(infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               anthropogenic_kernel_type = "power law")
  expect_equal(all(data$host_pools[[1]]$infected[[1]] >=
                     terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)

  data <- pops(infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               anthropogenic_kernel_type = "gamma")
  expect_equal(all(data$host_pools[[1]]$infected[[1]] >=
                     terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)), TRUE)

  #
  #   data <- pops(infected_file_list = infected_file_list,
  #                host_file_list = host_file_list,
  #                total_populations_file = total_populations_file,
  #                parameter_means = parameter_means,
  #                parameter_cov_matrix = parameter_cov_matrix,
  #                pest_host_table = pest_host_table,
  #                competency_table = competency_table,
  #                anthropogenic_kernel_type = "exponential-power")
  #   expect_equal(all(data$host_pools[[1]]$infected[[1]] >=
  #                      terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)),
  #                TRUE)

  ## currently not working
  #
  #
  # data <- pops(infected_file_list = infected_file_list,
  #              host_file_list = host_file_list,
  #              total_populations_file = total_populations_file,
  #              parameter_means = parameter_means,
  #              parameter_cov_matrix = parameter_cov_matrix,
  #              pest_host_table = pest_host_table,
  #              competency_table = competency_table,
  #              anthropogenic_kernel_type = "log normal")
  # expect_equal(all(data$host_pools[[1]]$infected[[1]] >=
  #                    terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)),
  #              TRUE)

})

test_that("Susceptibles are never negative", {
  infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file <-
    system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file <-
    system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2010-12-31"
  parameter_means <- c(0.4, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  data <- pops(infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               random_seed = 42,
               start_date = start_date,
               end_date = end_date)

  expect_equal(all(data$host_pools[[1]]$susceptible[[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data$host_pools[[1]]$susceptible[[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data$host_pools[[1]]$susceptible[[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)

  parameter_means <- c(0.5, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  data <-
    pops(infected_file_list = infected_file_list,
         host_file_list =
           system.file("extdata", "simple2x2",
                       "total_plants_host_greater_than_infected.tif", package = "PoPS"),
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         random_seed = 42,
         start_date = start_date,
         end_date = end_date)

  expect_equal(all(data$host_pools[[1]]$susceptible[[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data$host_pools[[1]]$susceptible[[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data$host_pools[[1]]$susceptible[[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)

})

test_that("SEI model works as intended", {
  infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file <-
    system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file <-
    system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2008-12-31"
  model_type <- "SI"
  latency_period <- 2
  time_step <- "month"
  output_frequency <- "month"
  treatment_dates <- "2008-02-25"
  parameter_means <- c(0.4, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  data <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         random_seed = 42,
         start_date = start_date,
         end_date = end_date,
         model_type = model_type,
         latency_period = latency_period,
         output_frequency = output_frequency,
         time_step = time_step,
         treatment_dates = treatment_dates)
  model_type <- "SEI"
  data2 <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         random_seed = 42,
         start_date = start_date,
         end_date = end_date,
         model_type = model_type,
         latency_period = latency_period,
         output_frequency = output_frequency,
         time_step = time_step,
         treatment_dates = treatment_dates)

  expect_equal(all(data2$exposed[[1]][[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[1]][[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[1]][[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[2]][[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[2]][[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[2]][[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[3]][[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[3]][[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[3]][[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[4]][[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[4]][[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[4]][[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[5]][[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[5]][[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[5]][[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[6]][[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[6]][[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[6]][[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[7]][[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[7]][[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[7]][[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[8]][[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[8]][[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[8]][[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[9]][[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[9]][[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[9]][[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[10]][[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[10]][[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[10]][[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[11]][[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[11]][[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[11]][[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[12]][[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[12]][[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data2$exposed[[12]][[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)

  expect_equal(all(data$host_pools[[1]]$susceptible[[1]] <=
                     data2$host_pools[[1]]$susceptible[[1]]), TRUE)
  expect_equal(all(data$host_pools[[1]]$susceptible[[2]] <=
                     data2$host_pools[[1]]$susceptible[[1]]), TRUE)
  expect_equal(all(data$host_pools[[1]]$susceptible[[3]] <=
                     data2$host_pools[[1]]$susceptible[[1]]), TRUE)

  expect_equal(all(data$host_pools[[1]]$infected[[1]] >= data2$host_pools[[1]]$infected[[1]]), TRUE)
  expect_equal(all(data$infected[[2]] >= data2$host_pools[[1]]$infected[[1]]), TRUE)
  expect_equal(all(data$infected[[3]] >= data2$host_pools[[1]]$infected[[1]]), TRUE)

  start_exposed <- TRUE
  exposed_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  model_type <- "SEI"
  data3 <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         random_seed = 42,
         start_date = start_date,
         end_date = end_date,
         model_type = model_type,
         latency_period = latency_period,
         output_frequency = output_frequency,
         time_step = time_step,
         treatment_dates = treatment_dates,
         start_exposed = start_exposed,
         exposed_file_list = exposed_file_list)

  expect_equal(all(data3$host_pools[[1]]$susceptible[[1]] <=
                     data2$host_pools[[1]]$susceptible[[1]]), TRUE)
  expect_equal(all(data3$host_pools[[1]]$susceptible[[2]] <=
                     data2$host_pools[[1]]$susceptible[[1]]), TRUE)
  expect_equal(all(data3$host_pools[[1]]$susceptible[[3]] <=
                     data2$host_pools[[1]]$susceptible[[1]]), TRUE)

  expect_equal(all(data3$host_pools[[1]]$infected[[1]] >=
                     data2$host_pools[[1]]$infected[[1]]), TRUE)
  expect_equal(all(data3$infected[[2]] >= data2$host_pools[[1]]$infected[[1]]), TRUE)
  expect_equal(all(data3$infected[[3]] >= data2$host_pools[[1]]$infected[[1]]), TRUE)

})

test_that("Infected results with weather are less than those without weather", {
  infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file <-
    system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file <-
    system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2010-12-31"
  parameter_means <- c(2.0, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  coefficient_sd_file <- system.file("extdata", "simple2x2", "coefficient_sd.tif", package = "PoPS")
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  data <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         random_seed = 42,
         start_date = start_date,
         end_date = end_date)

  data_temp <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         temp = TRUE,
         temperature_coefficient_file = coefficient_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         random_seed = 42,
         start_date = start_date,
         end_date = end_date)

  data_precip <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         precip = TRUE,
         precipitation_coefficient_file = coefficient_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         random_seed = 42,
         start_date = start_date,
         end_date = end_date)

  data_weather <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         temp = TRUE,
         temperature_coefficient_file = coefficient_file,
         precip = TRUE,
         precipitation_coefficient_file = coefficient_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         random_seed = 42,
         start_date = start_date,
         end_date = end_date)

  data_temp_wsd <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         temp = TRUE,
         temperature_coefficient_file = coefficient_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         random_seed = 42,
         start_date = start_date,
         end_date = end_date,
         weather_type = "probabilistic",
         temperature_coefficient_sd_file = coefficient_sd_file)

  data_precip_wsd <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         precip = TRUE,
         precipitation_coefficient_file = coefficient_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         random_seed = 42,
         start_date = start_date,
         end_date = end_date,
         weather_type = "probabilistic",
         precipitation_coefficient_sd_file = coefficient_sd_file)

  data_weather_wsd <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         temp = TRUE,
         temperature_coefficient_file = coefficient_file,
         precip = TRUE,
         precipitation_coefficient_file = coefficient_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         random_seed = 42,
         start_date = start_date,
         end_date = end_date,
         weather_type = "probabilistic",
         temperature_coefficient_sd_file = coefficient_sd_file,
         precipitation_coefficient_sd_file = coefficient_sd_file)

  expect_gte(sum(data$host_pools[[1]]$infected[[1]]), sum(data_temp$host_pools[[1]]$infected[[1]]))
  expect_gte(sum(data$host_pools[[1]]$infected[[2]]), sum(data_temp$host_pools[[1]]$infected[[2]]))
  expect_gte(sum(data$host_pools[[1]]$infected[[3]]), sum(data_temp$host_pools[[1]]$infected[[3]]))

  expect_gte(sum(data$host_pools[[1]]$infected[[1]]),
             sum(data_precip$host_pools[[1]]$infected[[1]]))
  expect_gte(sum(data$host_pools[[1]]$infected[[2]]),
             sum(data_precip$host_pools[[1]]$infected[[2]]))
  expect_gte(sum(data$host_pools[[1]]$infected[[3]]),
             sum(data_precip$host_pools[[1]]$infected[[3]]))

  expect_gte(sum(data$host_pools[[1]]$infected[[1]]),
             sum(data_weather$host_pools[[1]]$infected[[1]]))
  expect_gte(sum(data$host_pools[[1]]$infected[[2]]),
             sum(data_weather$host_pools[[1]]$infected[[2]]))
  expect_gte(sum(data$host_pools[[1]]$infected[[3]]),
             sum(data_weather$host_pools[[1]]$infected[[3]]))

  expect_gte(sum(data$host_pools[[1]]$infected[[2]]),
             sum(data_temp_wsd$host_pools[[1]]$infected[[2]]))
  expect_gte(sum(data$host_pools[[1]]$infected[[3]]),
             sum(data_temp_wsd$host_pools[[1]]$infected[[3]]))
  expect_gte(sum(data$host_pools[[1]]$infected[[1]]),
             sum(data_temp_wsd$host_pools[[1]]$infected[[1]]))

  expect_gte(sum(data$host_pools[[1]]$infected[[1]]),
             sum(data_precip_wsd$host_pools[[1]]$infected[[1]]))
  expect_gte(sum(data$host_pools[[1]]$infected[[2]]),
             sum(data_precip_wsd$host_pools[[1]]$infected[[2]]))
  expect_gte(sum(data$host_pools[[1]]$infected[[3]]),
             sum(data_precip_wsd$host_pools[[1]]$infected[[3]]))

  expect_gte(sum(data$host_pools[[1]]$infected[[1]]),
             sum(data_weather_wsd$host_pools[[1]]$infected[[1]]))
  expect_gte(sum(data$host_pools[[1]]$infected[[2]]),
             sum(data_weather_wsd$host_pools[[1]]$infected[[2]]))
  expect_gte(sum(data$host_pools[[1]]$infected[[3]]),
             sum(data_weather_wsd$host_pools[[1]]$infected[[3]]))
})

test_that(
  "Infected results are greater with same parameters for weekly spread vs. monthly", {
    infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
    host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
    total_populations_file <-
      system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
    coefficient_file <-
      system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
    temperature_file <-
      system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
    start_date <- "2008-01-01"
    end_date <- "2010-12-31"
    parameter_means <- c(0.2, 21, 1, 500, 0, 0, 0, 0)
    parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
    pest_host_table <-
      system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
    competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

    data_week <-
      pops(infected_file_list = infected_file_list,
           host_file_list = host_file_list,
           total_populations_file = total_populations_file,
           parameter_means = parameter_means,
           parameter_cov_matrix = parameter_cov_matrix,
           pest_host_table = pest_host_table,
           competency_table = competency_table,
           time_step = "week",
           random_seed = 42,
           start_date = start_date,
           end_date = end_date)
    data_month <-
      pops(infected_file_list = infected_file_list,
           host_file_list = host_file_list,
           total_populations_file = total_populations_file,
           parameter_means = parameter_means,
           parameter_cov_matrix = parameter_cov_matrix,
           pest_host_table = pest_host_table,
           competency_table = competency_table,
           time_step = "month",
           random_seed = 42,
           start_date = start_date,
           end_date = end_date)

    expect_equal(all(data_week$host_pools[[1]]$infected[[1]] >=
                       data_month$host_pools[[1]]$infected[[1]]), TRUE)
    expect_equal(all(data_week$host_pools[[1]]$infected[[2]] >=
                       data_month$host_pools[[1]]$infected[[2]]), TRUE)

  })

test_that("Infected results are greater with same parameters for daily spread vs. monthly and
          weekly", {
            infected_file_list <-
              system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
            host_file_list <-
              system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
            total_populations_file <-
              system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
            coefficient_file <-
              system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
            temperature_file <-
              system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif",
                          package = "PoPS")
            start_date <- "2008-01-01"
            end_date <- "2010-12-31"
            parameter_means <- c(0.1, 21, 1, 500, 0, 0, 0, 0)
            parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
            pest_host_table <-
              system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
            competency_table <-
              system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

            data_day <-
              pops(infected_file_list = infected_file_list,
                   host_file_list = host_file_list,
                   total_populations_file = total_populations_file,
                   parameter_means = parameter_means,
                   parameter_cov_matrix = parameter_cov_matrix,
                   pest_host_table = pest_host_table,
                   competency_table = competency_table,
                   time_step = "day",
                   random_seed = 42)
            data_week <-
              pops(infected_file_list = infected_file_list,
                   host_file_list = host_file_list,
                   total_populations_file = total_populations_file,
                   parameter_means = parameter_means,
                   parameter_cov_matrix = parameter_cov_matrix,
                   pest_host_table = pest_host_table,
                   competency_table = competency_table,
                   time_step = "week",
                   random_seed = 42,
                   start_date = start_date,
                   end_date = end_date)
            data_month <-
              pops(infected_file_list = infected_file_list,
                   host_file_list = host_file_list,
                   total_populations_file = total_populations_file,
                   parameter_means = parameter_means,
                   parameter_cov_matrix = parameter_cov_matrix,
                   pest_host_table = pest_host_table,
                   competency_table = competency_table,
                   time_step = "month",
                   random_seed = 42,
                   start_date = start_date,
                   end_date = end_date)

            expect_equal(all(data_day$host_pools[[1]]$infected[[1]] >=
                               data_month$host_pools[[1]]$infected[[1]]), TRUE)
            expect_equal(all(data_day$host_pools[[1]]$infected[[1]] >=
                               data_week$host_pools[[1]]$infected[[1]]), TRUE)
            expect_equal(all(data_week$host_pools[[1]]$infected[[1]] >=
                               data_month$host_pools[[1]]$infected[[1]]), TRUE)
          })

test_that(
  "Infected results are greater without treatment than with treatment", {
    infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
    host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
    total_populations_file <-
      system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
    treatments_file <- system.file("extdata", "simple2x2", "treatments_1_1.tif", package = "PoPS")
    treatment_dates <- c("2008-03-05")
    start_date <- "2008-01-01"
    end_date <- "2009-12-31"
    parameter_means <- c(0.8, 21, 1, 500, 0, 0, 0, 0)
    parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
    pest_host_table <-
      system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
    competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

    data <-
      pops(infected_file_list = infected_file_list,
           host_file_list = host_file_list,
           total_populations_file = total_populations_file,
           parameter_means = parameter_means,
           parameter_cov_matrix = parameter_cov_matrix,
           pest_host_table = pest_host_table,
           competency_table = competency_table,
           random_seed = 42,
           start_date = start_date,
           end_date = end_date)
    data_treat <-
      pops(infected_file_list = infected_file_list,
           host_file_list = host_file_list,
           total_populations_file = total_populations_file,
           management = TRUE,
           treatment_dates = treatment_dates,
           treatments_file = treatments_file,
           parameter_means = parameter_means,
           parameter_cov_matrix = parameter_cov_matrix,
           pest_host_table = pest_host_table,
           competency_table = competency_table,
           random_seed = 42,
           start_date = start_date,
           end_date = end_date)

    expect_equal(all(data$host_pools[[1]]$infected[[1]] >=
                       data_treat$host_pools[[1]]$infected[[1]]), TRUE)
    expect_equal(all(data$host_pools[[1]]$infected[[2]] >=
                       data_treat$host_pools[[1]]$infected[[2]]), TRUE)
  })

test_that("Infected results are greater with higher reproductive rate", {
  infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2010-12-31"
  parameter_means <- c(1.0, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  data_1 <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         time_step = "month",
         random_seed = 42,
         start_date = start_date,
         end_date = end_date)
  parameter_means <- c(0.75, 21, 1, 500, 0, 0, 0, 0)
  data_075 <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         time_step = "month",
         random_seed = 42,
         start_date = start_date,
         end_date = end_date)
  parameter_means <- c(0.5, 21, 1, 500, 0, 0, 0, 0)
  data_050 <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         time_step = "month",
         random_seed = 42,
         start_date = start_date,
         end_date = end_date)
  parameter_means <- c(0.25, 21, 1, 500, 0, 0, 0, 0)
  data_025 <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         time_step = "month",
         random_seed = 42,
         start_date = start_date,
         end_date = end_date)
  parameter_means <- c(0.1, 21, 1, 500, 0, 0, 0, 0)
  data_010 <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         time_step = "month",
         random_seed = 42,
         start_date = start_date,
         end_date = end_date)

  expect_gte(sum(data_1$host_pools[[1]]$infected[[1]]), sum(data_075$host_pools[[1]]$infected[[1]]))
  expect_gte(sum(data_1$host_pools[[1]]$infected[[1]]), sum(data_050$host_pools[[1]]$infected[[1]]))
  expect_gte(sum(data_1$host_pools[[1]]$infected[[1]]), sum(data_025$host_pools[[1]]$infected[[1]]))
  expect_gte(sum(data_1$host_pools[[1]]$infected[[1]]), sum(data_010$host_pools[[1]]$infected[[1]]))

  expect_gte(sum(data_075$host_pools[[1]]$infected[[1]]),
             sum(data_050$host_pools[[1]]$infected[[1]]))
  expect_gte(sum(data_075$host_pools[[1]]$infected[[1]]),
             sum(data_025$host_pools[[1]]$infected[[1]]))
  expect_gte(sum(data_075$host_pools[[1]]$infected[[1]]),
             sum(data_010$host_pools[[1]]$infected[[1]]))

  expect_gte(sum(data_050$host_pools[[1]]$infected[[1]]),
             sum(data_025$host_pools[[1]]$infected[[1]]))
  expect_gte(sum(data_050$host_pools[[1]]$infected[[2]]),
             sum(data_025$host_pools[[1]]$infected[[2]]))
  expect_gte(sum(data_050$host_pools[[1]]$infected[[1]]),
             sum(data_010$host_pools[[1]]$infected[[1]]))
  expect_gte(sum(data_050$host_pools[[1]]$infected[[2]]),
             sum(data_010$host_pools[[1]]$infected[[2]]))

  expect_gte(sum(data_025$host_pools[[1]]$infected[[1]]),
             sum(data_010$host_pools[[1]]$infected[[1]]))
  expect_gte(sum(data_025$host_pools[[1]]$infected[[2]]),
             sum(data_010$host_pools[[1]]$infected[[2]]))
})

test_that("Treatments apply no matter what time step", {
  infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file <-
    system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file <-
    system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date <- "2009-01-01"
  end_date <- "2009-12-31"
  treatments_file <- system.file("extdata", "simple2x2", "treatments.tif", package = "PoPS")
  parameter_means <- c(0, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  dates <- seq.Date(as.Date(start_date), as.Date(end_date), by = "days")
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")
  for (i in seq_len(length(dates))) {
    data <-
      pops(infected_file_list = infected_file_list,
           host_file_list = host_file_list,
           total_populations_file = total_populations_file,
           management  = TRUE,
           treatment_dates = c(as.character(dates[i])),
           treatments_file = treatments_file,
           parameter_means = parameter_means,
           parameter_cov_matrix = parameter_cov_matrix,
           pest_host_table = pest_host_table,
           competency_table = competency_table,
           start_date = start_date,
           end_date = end_date)
    expect_equal(data$host_pools[[1]]$infected[[1]], matrix(0, ncol = 2, nrow = 2))
    expect_equal(data$host_pools[[1]]$susceptible[[1]], matrix(0, ncol = 2, nrow = 2))
  }
})

test_that("Pesticide treatments apply no matter what time step", {
  infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file <-
    system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file <-
    system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date <- "2009-01-01"
  end_date <- "2009-12-31"
  treatments_file <- system.file("extdata", "simple2x2", "treatments.tif", package = "PoPS")
  pesticide_duration <- c(120)
  pesticide_efficacy <- 1.0
  parameter_means <- c(0, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  dates <- seq.Date(as.Date(start_date), as.Date("2009-06-30"), by = "days")
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  for (i in seq_len(length(dates))) {
    data <-
      pops(infected_file_list = infected_file_list,
           host_file_list = host_file_list,
           total_populations_file = total_populations_file,
           management  = TRUE,
           treatment_dates = c(as.character(dates[i])),
           treatments_file = treatments_file,
           parameter_means = parameter_means,
           parameter_cov_matrix = parameter_cov_matrix,
           pest_host_table = pest_host_table,
           competency_table = competency_table,
           start_date = start_date,
           end_date = end_date,
           pesticide_duration = pesticide_duration,
           pesticide_efficacy = pesticide_efficacy)
    expect_equal(data$host_pools[[1]]$infected[[1]], matrix(0, ncol = 2, nrow = 2))
    expect_equal(data$host_pools[[1]]$susceptible[[1]],
                 terra::as.matrix(terra::rast(host_file_list), wide = TRUE))
  }

  pesticide_duration <- c(120)
  pesticide_efficacy <- 0.5

  for (i in seq_len(length(dates))) {
    data <-
      pops(infected_file_list = infected_file_list,
           host_file_list = host_file_list,
           total_populations_file = total_populations_file,
           management  = TRUE,
           treatment_dates = c(as.character(dates[i])),
           treatments_file = treatments_file,
           parameter_means = parameter_means,
           parameter_cov_matrix = parameter_cov_matrix,
           pest_host_table = pest_host_table,
           competency_table = competency_table,
           start_date = start_date,
           end_date = end_date,
           pesticide_duration = pesticide_duration,
           pesticide_efficacy = pesticide_efficacy)
    expect_equal(data$host_pools[[1]]$infected[[1]], matrix(c(3, 0, 0, 0), ncol = 2, nrow = 2))
    expect_equal(data$host_pools[[1]]$susceptible[[1]],
                 matrix(c(14, 14, 6, 15), ncol = 2, nrow = 2))
  }
})

test_that("Changing the output frequency returns the correct number of outputs and output
          statistics", {
            infected_file_list <-
              system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
            host_file_list <-
              system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
            total_populations_file <-
              system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
            start_date <- "2009-01-01"
            end_date <- "2009-12-31"
            treatment_dates <- c(start_date)
            parameter_means <- c(0, 21, 1, 500, 0, 0, 0, 0)
            parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
            pest_host_table <-
              system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
            competency_table <-
              system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

            data <-
              pops(output_frequency = "year",
                   time_step = "month",
                   treatment_dates = treatment_dates,
                   infected_file_list = infected_file_list,
                   host_file_list = host_file_list,
                   total_populations_file = total_populations_file,
                   parameter_means = parameter_means,
                   parameter_cov_matrix = parameter_cov_matrix,
                   pest_host_table = pest_host_table,
                   competency_table = competency_table,
                   start_date = start_date,
                   end_date = end_date)
            expect_equal(length(data$host_pools[[1]]$infected), 1)

            data <-
              pops(output_frequency = "year",
                   time_step = "week",
                   treatment_dates = start_date,
                   infected_file_list = infected_file_list,
                   host_file_list = host_file_list,
                   total_populations_file = total_populations_file,
                   parameter_means = parameter_means,
                   parameter_cov_matrix = parameter_cov_matrix,
                   pest_host_table = pest_host_table,
                   competency_table = competency_table,
                   start_date = start_date,
                   end_date = end_date)
            expect_equal(length(data$host_pools[[1]]$infected), 1)

            data <-
              pops(output_frequency = "year",
                   time_step = "day",
                   treatment_dates = start_date,
                   infected_file_list = infected_file_list,
                   host_file_list = host_file_list,
                   total_populations_file = total_populations_file,
                   parameter_means = parameter_means,
                   parameter_cov_matrix = parameter_cov_matrix,
                   pest_host_table = pest_host_table,
                   competency_table = competency_table,
                   start_date = start_date,
                   end_date = end_date)
            expect_equal(length(data$host_pools[[1]]$infected), 1)

            data <-
              pops(output_frequency = "month",
                   time_step = "week",
                   treatment_dates = start_date,
                   infected_file_list = infected_file_list,
                   host_file_list = host_file_list,
                   total_populations_file = total_populations_file,
                   parameter_means = parameter_means,
                   parameter_cov_matrix = parameter_cov_matrix,
                   pest_host_table = pest_host_table,
                   competency_table = competency_table,
                   start_date = start_date,
                   end_date = end_date)
            expect_equal(length(data$host_pools[[1]]$infected), 12)

            data <-
              pops(output_frequency = "month",
                   time_step = "day",
                   treatment_dates = start_date,
                   infected_file_list = infected_file_list,
                   host_file_list = host_file_list,
                   total_populations_file = total_populations_file,
                   parameter_means = parameter_means,
                   parameter_cov_matrix = parameter_cov_matrix,
                   pest_host_table = pest_host_table,
                   competency_table = competency_table,
                   start_date = start_date,
                   end_date = end_date)
            expect_equal(length(data$host_pools[[1]]$infected), 12)

            data <-
              pops(output_frequency = "week",
                   time_step = "week",
                   treatment_dates = start_date,
                   infected_file_list = infected_file_list,
                   host_file_list = host_file_list,
                   total_populations_file = total_populations_file,
                   parameter_means = parameter_means,
                   parameter_cov_matrix = parameter_cov_matrix,
                   pest_host_table = pest_host_table,
                   competency_table = competency_table,
                   start_date = start_date,
                   end_date = end_date)
            expect_equal(length(data$host_pools[[1]]$infected), 52)

            data <-
              pops(output_frequency = "week",
                   time_step = "day",
                   treatment_dates = start_date,
                   infected_file_list = infected_file_list,
                   host_file_list = host_file_list,
                   total_populations_file = total_populations_file,
                   parameter_means = parameter_means,
                   parameter_cov_matrix = parameter_cov_matrix,
                   pest_host_table = pest_host_table,
                   competency_table = competency_table,
                   start_date = start_date,
                   end_date = end_date)
            expect_equal(length(data$host_pools[[1]]$infected), 52)

            expect_error(pops(output_frequency = "day",
                              time_step = "week",
                              treatment_dates = start_date,
                              infected_file_list = infected_file_list,
                              host_file_list = host_file_list,
                              total_populations_file = total_populations_file,
                              parameter_means = parameter_means,
                              parameter_cov_matrix = parameter_cov_matrix,
                              pest_host_table = pest_host_table,
                              competency_table = competency_table,
                              start_date = start_date,
                              end_date = end_date), output_frequency_error)

            expect_error(pops(output_frequency = "day",
                              time_step = "month",
                              treatment_dates = start_date,
                              infected_file_list = infected_file_list,
                              host_file_list = host_file_list,
                              total_populations_file = total_populations_file,
                              parameter_means = parameter_means,
                              parameter_cov_matrix = parameter_cov_matrix,
                              pest_host_table = pest_host_table,
                              competency_table = competency_table,
                              start_date = start_date,
                              end_date = end_date), output_frequency_error)

            expect_error(pops(output_frequency = "week",
                              time_step = "month",
                              treatment_dates = start_date,
                              infected_file_list = infected_file_list,
                              host_file_list = host_file_list,
                              total_populations_file = total_populations_file,
                              parameter_means = parameter_means,
                              parameter_cov_matrix = parameter_cov_matrix,
                              pest_host_table = pest_host_table,
                              competency_table = competency_table,
                              start_date = start_date,
                              end_date = end_date), output_frequency_error)

            data <- pops(output_frequency = "day",
                         time_step = "day",
                         treatment_dates = start_date,
                         infected_file_list = infected_file_list,
                         host_file_list = host_file_list,
                         total_populations_file = total_populations_file,
                         parameter_means = parameter_means,
                         parameter_cov_matrix = parameter_cov_matrix,
                         pest_host_table = pest_host_table,
                         competency_table = competency_table,
                         start_date = start_date,
                         end_date = end_date)
            expect_equal(length(data$host_pools[[1]]$infected), 364)

            data <- pops(output_frequency = "every_n_steps",
                         output_frequency_n = 5,
                         time_step = "day",
                         treatment_dates = start_date,
                         infected_file_list = infected_file_list,
                         host_file_list = host_file_list,
                         total_populations_file = total_populations_file,
                         parameter_means = parameter_means,
                         parameter_cov_matrix = parameter_cov_matrix,
                         pest_host_table = pest_host_table,
                         competency_table = competency_table,
                         start_date = start_date,
                         end_date = end_date)
            expect_equal(length(data$host_pools[[1]]$infected), 72)
          })

test_that(
  "Outputs occur with non-full year date range for all time step output frequency combinations", {
    infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
    host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
    total_populations_file <-
      system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
    start_date <- "2009-05-01"
    end_date <- "2009-10-29"
    treatment_dates <- start_date
    parameter_means <- c(0, 21, 1, 500, 0, 0, 0, 0)
    parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
    pest_host_table <-
      system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
    competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

    data <- pops(output_frequency = "year",
                 time_step = "month",
                 treatment_dates = treatment_dates,
                 infected_file_list = infected_file_list,
                 host_file_list = host_file_list,
                 total_populations_file = total_populations_file,
                 parameter_means = parameter_means,
                 parameter_cov_matrix = parameter_cov_matrix,
                 pest_host_table = pest_host_table,
                 competency_table = competency_table,
                 start_date = start_date,
                 end_date = end_date)
    expect_equal(length(data$host_pools[[1]]$infected), 1)

    data <- pops(output_frequency = "year",
                 time_step = "week",
                 treatment_dates = start_date,
                 infected_file_list = infected_file_list,
                 host_file_list = host_file_list,
                 total_populations_file = total_populations_file,
                 parameter_means = parameter_means,
                 parameter_cov_matrix = parameter_cov_matrix,
                 pest_host_table = pest_host_table,
                 competency_table = competency_table,
                 start_date = start_date,
                 end_date = end_date)
    expect_equal(length(data$host_pools[[1]]$infected), 1)

    data <- pops(output_frequency = "year",
                 time_step = "day",
                 treatment_dates = start_date,
                 infected_file_list = infected_file_list,
                 host_file_list = host_file_list,
                 total_populations_file = total_populations_file,
                 parameter_means = parameter_means,
                 parameter_cov_matrix = parameter_cov_matrix,
                 pest_host_table = pest_host_table,
                 competency_table = competency_table,
                 start_date = start_date,
                 end_date = end_date)
    expect_equal(length(data$host_pools[[1]]$infected), 1)

    data <- pops(output_frequency = "month",
                 time_step = "week",
                 treatment_dates = start_date,
                 infected_file_list = infected_file_list,
                 host_file_list = host_file_list,
                 total_populations_file = total_populations_file,
                 parameter_means = parameter_means,
                 parameter_cov_matrix = parameter_cov_matrix,
                 pest_host_table = pest_host_table,
                 competency_table = competency_table,
                 start_date = start_date,
                 end_date = end_date)
    expect_equal(length(data$host_pools[[1]]$infected), 5)

    data <- pops(output_frequency = "month",
                 time_step = "day",
                 treatment_dates = start_date,
                 infected_file_list = infected_file_list,
                 host_file_list = host_file_list,
                 total_populations_file = total_populations_file,
                 parameter_means = parameter_means,
                 parameter_cov_matrix = parameter_cov_matrix,
                 pest_host_table = pest_host_table,
                 competency_table = competency_table,
                 start_date = start_date,
                 end_date = end_date)
    expect_equal(length(data$host_pools[[1]]$infected), 5)

    data <- pops(output_frequency = "week",
                 time_step = "week",
                 treatment_dates = start_date,
                 infected_file_list = infected_file_list,
                 host_file_list = host_file_list,
                 total_populations_file = total_populations_file,
                 parameter_means = parameter_means,
                 parameter_cov_matrix = parameter_cov_matrix,
                 pest_host_table = pest_host_table,
                 competency_table = competency_table,
                 start_date = start_date,
                 end_date = end_date)
    expect_equal(length(data$host_pools[[1]]$infected), 26)

    data <- pops(output_frequency = "week",
                 time_step = "day",
                 treatment_dates = start_date,
                 infected_file_list = infected_file_list,
                 host_file_list = host_file_list,
                 total_populations_file = total_populations_file,
                 parameter_means = parameter_means,
                 parameter_cov_matrix = parameter_cov_matrix,
                 pest_host_table = pest_host_table,
                 competency_table = competency_table,
                 start_date = start_date,
                 end_date = end_date)
    expect_equal(length(data$host_pools[[1]]$infected), 26)

    data <-
      pops(output_frequency = "day"
           , time_step = "day",
           treatment_dates = start_date,
           infected_file_list = infected_file_list,
           host_file_list = host_file_list,
           total_populations_file = total_populations_file,
           parameter_means = parameter_means,
           parameter_cov_matrix = parameter_cov_matrix,
           pest_host_table = pest_host_table,
           competency_table = competency_table,
           start_date = start_date,
           end_date = end_date)
    expect_equal(length(data$host_pools[[1]]$infected), 182)

    data <-
      pops(output_frequency = "time_step",
           time_step = "day",
           treatment_dates = start_date,
           infected_file_list = infected_file_list,
           host_file_list = host_file_list,
           total_populations_file = total_populations_file,
           parameter_means = parameter_means,
           parameter_cov_matrix = parameter_cov_matrix,
           pest_host_table = pest_host_table,
           competency_table = competency_table,
           start_date = start_date,
           end_date = end_date)
    expect_equal(length(data$host_pools[[1]]$infected), 182)
  })

test_that("Quarantine and spread rates work at all timings", {
  infected_file_list <-
    system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
  start_date <- "2009-01-01"
  end_date <- "2009-12-31"
  treatment_dates <- start_date
  parameter_means <- c(0, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  quarantine_areas_file <-
    system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  data <- pops(output_frequency = "year",
               time_step = "month",
               treatment_dates = start_date,
               infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               start_date = start_date,
               end_date = end_date,
               use_quarantine = TRUE,
               use_spreadrates = TRUE,
               quarantine_areas_file = quarantine_areas_file)
  expect_equal(length(data$host_pools[[1]]$infected), 1)
  expect_equal(length(data$quarantine_escape), 1)
  expect_equal(length(data$quarantine_escape_distance), 1)
  expect_equal(length(data$quarantine_escape_directions), 1)
  expect_equal(length(data$rates), 1)

  data <- pops(output_frequency = "year",
               time_step = "week",
               treatment_dates = start_date,
               infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               start_date = start_date,
               end_date = end_date,
               use_quarantine = TRUE,
               use_spreadrates = TRUE,
               quarantine_areas_file = quarantine_areas_file)
  expect_equal(length(data$host_pools[[1]]$infected), 1)
  expect_equal(length(data$quarantine_escape), 1)
  expect_equal(length(data$quarantine_escape_distance), 1)
  expect_equal(length(data$quarantine_escape_directions), 1)
  expect_equal(length(data$rates), 1)

  data <- pops(output_frequency = "year",
               time_step = "day",
               treatment_dates = start_date,
               infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               start_date = start_date,
               end_date = end_date,
               use_quarantine = TRUE,
               use_spreadrates = TRUE,
               quarantine_areas_file = quarantine_areas_file)
  expect_equal(length(data$host_pools[[1]]$infected), 1)
  expect_equal(length(data$quarantine_escape), 1)
  expect_equal(length(data$quarantine_escape_distance), 1)
  expect_equal(length(data$quarantine_escape_directions), 1)
  expect_equal(length(data$rates), 1)

  data <- pops(output_frequency = "month",
               time_step = "week",
               treatment_dates = start_date,
               infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               start_date = start_date,
               end_date = end_date,
               use_quarantine = TRUE,
               use_spreadrates = TRUE,
               quarantine_areas_file = quarantine_areas_file)
  expect_equal(length(data$host_pools[[1]]$infected), 12)
  expect_equal(length(data$quarantine_escape), 12)
  expect_equal(length(data$quarantine_escape_distance), 12)
  expect_equal(length(data$quarantine_escape_directions), 12)
  expect_equal(length(data$rates), 12)

  data <- pops(output_frequency = "month",
               time_step = "day",
               treatment_dates = start_date,
               infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               start_date = start_date,
               end_date = end_date,
               use_quarantine = TRUE,
               use_spreadrates = TRUE,
               quarantine_areas_file = quarantine_areas_file)
  expect_equal(length(data$host_pools[[1]]$infected), 12)
  expect_equal(length(data$quarantine_escape), 12)
  expect_equal(length(data$quarantine_escape_distance), 12)
  expect_equal(length(data$quarantine_escape_directions), 12)
  expect_equal(length(data$rates), 12)

  data <- pops(output_frequency = "week",
               time_step = "week",
               treatment_dates = start_date,
               infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               start_date = start_date,
               end_date = end_date,
               use_quarantine = TRUE,
               use_spreadrates = TRUE,
               quarantine_areas_file = quarantine_areas_file)
  expect_equal(length(data$host_pools[[1]]$infected), 52)
  expect_equal(length(data$quarantine_escape), 52)
  expect_equal(length(data$quarantine_escape_distance), 52)
  expect_equal(length(data$quarantine_escape_directions), 52)
  expect_equal(length(data$rates), 52)

  data <- pops(output_frequency = "week",
               time_step = "day",
               treatment_dates = start_date,
               infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               start_date = start_date,
               end_date = end_date,
               use_quarantine = TRUE,
               use_spreadrates = TRUE,
               quarantine_areas_file = quarantine_areas_file)
  expect_equal(length(data$host_pools[[1]]$infected), 52)
  expect_equal(length(data$quarantine_escape), 52)
  expect_equal(length(data$quarantine_escape_distance), 52)
  expect_equal(length(data$quarantine_escape_directions), 52)
  expect_equal(length(data$rates), 52)

  data <- pops(output_frequency = "day",
               time_step = "day",
               treatment_dates = start_date,
               infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               start_date = start_date,
               end_date = end_date, use_quarantine = TRUE,
               use_spreadrates = TRUE,
               quarantine_areas_file = quarantine_areas_file,
               quarantine_directions = "N,E,S,W")
  expect_equal(length(data$host_pools[[1]]$infected), 364)
  expect_equal(length(data$quarantine_escape), 364)
  expect_equal(length(data$quarantine_escape_distance), 364)
  expect_equal(length(data$quarantine_escape_directions), 364)
  expect_equal(length(data$rates), 364)

  data <- pops(output_frequency = "time_step",
               time_step = "day",
               treatment_dates = start_date,
               infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               start_date = start_date,
               end_date = end_date,
               use_quarantine = TRUE,
               use_spreadrates = TRUE,
               quarantine_areas_file = quarantine_areas_file,
               quarantine_directions = "N")
  expect_equal(length(data$host_pools[[1]]$infected), 364)
  expect_equal(length(data$quarantine_escape), 364)
  expect_equal(length(data$quarantine_escape_distance), 364)
  expect_equal(length(data$quarantine_escape_directions), 364)
  expect_equal(length(data$rates), 364)
})

test_that("Mortality works as expected with multiple ", {
  infected_file_list <-
    system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
  start_date <- "2009-01-01"
  end_date <- "2009-12-31"
  treatment_dates <- start_date
  parameter_means <- c(0, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  pest_host_table <- system.file("extdata", "pest_host_table_singlehost.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  data <- pops(output_frequency = "month",
               time_step = "month",
               treatment_dates = start_date,
               infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               start_date = start_date,
               end_date = end_date,
               mortality_frequency = "month",
               mortality_frequency_n = 1)

  expect_equal(length(data$host_pools[[1]]$mortality), 12)
  expect_equal(data$host_pools[[1]]$mortality[[1]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[2]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[3]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))

  data <- pops(output_frequency = "week",
               time_step = "week",
               treatment_dates = start_date,
               infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               start_date = start_date,
               end_date = end_date,
               mortality_frequency = "month",
               mortality_frequency_n = 1)

  expect_equal(length(data$host_pools[[1]]$mortality), 12)
  expect_equal(data$host_pools[[1]]$mortality[[1]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[2]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[3]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))

  pest_host_table <- system.file("extdata", "pest_host_table_singlehost025.csv", package = "PoPS")
  data <- pops(output_frequency = "week",
               time_step = "week",
               treatment_dates = start_date,
               infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               start_date = start_date,
               end_date = end_date,
               mortality_frequency = "month",
               mortality_frequency_n = 1)

  expect_equal(length(data$host_pools[[1]]$mortality), 12)
  expect_equal(data$host_pools[[1]]$mortality[[1]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[2]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[3]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[4]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[5]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))

  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost025tl3.csv", package = "PoPS")
  data <- pops(output_frequency = "week",
               time_step = "week",
               treatment_dates = start_date,
               infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               start_date = start_date,
               end_date = end_date,
               mortality_frequency = "month",
               mortality_frequency_n = 1)

  expect_equal(length(data$host_pools[[1]]$mortality), 12)
  expect_equal(data$host_pools[[1]]$mortality[[1]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[2]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[3]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[4]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[5]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[6]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[7]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))


  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost010tl1.csv", package = "PoPS")
  data <- pops(output_frequency = "week",
               time_step = "week",
               treatment_dates = start_date,
               infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               start_date = start_date,
               end_date = end_date,
               mortality_frequency = "month",
               mortality_frequency_n = 1)

  expect_equal(length(data$host_pools[[1]]$mortality), 12)
  expect_equal(data$host_pools[[1]]$mortality[[1]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[2]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[3]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[4]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[5]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[6]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[7]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[8]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[9]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[10]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[11]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
})

test_that("Movements works as expected", {
  infected_file_list <-
    system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
  start_date <- "2009-01-01"
  end_date <- "2009-12-31"
  treatment_dates <- start_date
  parameter_means <- c(0, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  use_movements <- TRUE
  movements_file <- system.file("extdata", "simple20x20", "movements.tif", package = "PoPS")
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  expect_error(pops(output_frequency = "month",
                    time_step = "month",
                    treatment_dates = start_date,
                    infected_file_list = infected_file_list,
                    host_file_list = host_file_list,
                    total_populations_file = total_populations_file,
                    parameter_means = parameter_means,
                    parameter_cov_matrix = parameter_cov_matrix,
                    pest_host_table = pest_host_table,
                    competency_table = competency_table,
                    start_date = start_date,
                    end_date = end_date,
                    use_movements = use_movements,
                    movements_file = movements_file,
                    random_seed = 42),
               detailed_file_exists_error(movements_file))

  movements_file <- system.file("extdata", "simple20x20", "movements.csv", package = "PoPS")
  data <- pops(output_frequency = "month",
               time_step = "month",
               treatment_dates = start_date,
               infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               start_date = start_date,
               end_date = end_date,
               use_movements = use_movements,
               movements_file = movements_file,
               random_seed = 42)

  expect_equal(length(data$host_pools[[1]]$infected), 12)
  expect_equal(data$host_pools[[1]]$infected[[1]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  expect_equal(data$host_pools[[1]]$infected[[2]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  expect_equal(data$host_pools[[1]]$infected[[3]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  expect_equal(data$host_pools[[1]]$infected[[4]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  infected_move <- matrix(0, ncol = 20, nrow = 20)
  infected_move[2, 1] <- 1
  expect_equal(data$host_pools[[1]]$infected[[5]], infected_move)
  sus <- terra::rast(host_file_list) - terra::rast(infected_file_list)
  sus <- terra::as.matrix(sus, wide = TRUE)
  sus5 <- sus
  sus5[1, 1] <- sus5[1, 1] - 199
  sus5[2, 1] <- sus5[2, 1] + 199
  sus6 <- sus5
  sus6[1, 2] <- sus6[1, 2] - 50
  sus6[2, 2] <- sus6[2, 2] + 50
  expect_equal(data$host_pools[[1]]$susceptible[[1]], sus)
  expect_equal(data$host_pools[[1]]$susceptible[[2]], sus)
  expect_equal(data$host_pools[[1]]$susceptible[[3]], sus)
  expect_equal(data$host_pools[[1]]$susceptible[[4]], sus)
  expect_equal(data$host_pools[[1]]$susceptible[[5]], sus5)
  expect_equal(data$host_pools[[1]]$susceptible[[6]], sus6)

  data <- pops(output_frequency = "month",
               time_step = "month",
               treatment_dates = start_date,
               infected_file_list = infected_file_list,
               host_file_list = host_file_list,
               total_populations_file = total_populations_file,
               parameter_means = parameter_means,
               parameter_cov_matrix = parameter_cov_matrix,
               pest_host_table = pest_host_table,
               competency_table = competency_table,
               start_date = start_date,
               end_date = end_date,
               use_movements = use_movements,
               movements_file = movements_file,
               random_seed = 45)

  expect_equal(length(data$host_pools[[1]]$infected), 12)
  expect_equal(data$host_pools[[1]]$infected[[1]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  expect_equal(data$host_pools[[1]]$infected[[2]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  expect_equal(data$host_pools[[1]]$infected[[3]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  expect_equal(data$host_pools[[1]]$infected[[4]],
               terra::as.matrix(terra::rast(infected_file_list), wide = TRUE))
  infected_move <- matrix(0, ncol = 20, nrow = 20)
  infected_move[2, 1] <- 1
  expect_equal(data$host_pools[[1]]$infected[[5]], infected_move)
  sus <- terra::rast(host_file_list) - terra::rast(infected_file_list)
  sus <- terra::as.matrix(sus, wide = TRUE)
  sus5 <- sus
  sus5[1, 1] <- sus5[1, 1] - 199
  sus5[2, 1] <- sus5[2, 1] + 199
  sus6 <- sus5
  sus6[1, 2] <- sus6[1, 2] - 50
  sus6[2, 2] <- sus6[2, 2] + 50
  expect_equal(data$host_pools[[1]]$susceptible[[1]], sus)
  expect_equal(data$host_pools[[1]]$susceptible[[2]], sus)
  expect_equal(data$host_pools[[1]]$susceptible[[3]], sus)
  expect_equal(data$host_pools[[1]]$susceptible[[4]], sus)
  expect_equal(data$host_pools[[1]]$susceptible[[5]], sus5)
  expect_equal(data$host_pools[[1]]$susceptible[[6]], sus6)
})

test_that(
  "Overpopulation dispersal works as expected with directionality to prevent dispersers from
  leaving the simulated area", {
    infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
    host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
    total_populations_file <-
      system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
    start_date <- "2008-01-01"
    end_date <- "2008-12-31"
    parameter_means <- c(0, 21, 1, 500, 0, 0, 0, 0)
    parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
    pest_host_table <-
      system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
    competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

    data <-
      pops(infected_file_list = infected_file_list,
           host_file_list = host_file_list,
           total_populations_file = total_populations_file,
           parameter_means = parameter_means,
           parameter_cov_matrix = parameter_cov_matrix,
           pest_host_table = pest_host_table,
           competency_table = competency_table,
           start_date = start_date,
           end_date = end_date,
           use_overpopulation_movements = TRUE,
           overpopulation_percentage = 0.2,
           leaving_percentage = 0.5,
           leaving_scale_coefficient = 0.5,
           natural_dir = "SE")
    test_mat <- terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)
    expect_lte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
    expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
    expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
    expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
  })

test_that("Deterministic dispersal works as expected", {
  infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2008-12-31"
  parameter_means <- c(2, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  data <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         start_date = start_date,
         end_date = end_date,
         generate_stochasticity = FALSE,
         establishment_stochasticity = FALSE,
         movement_stochasticity = FALSE,
         dispersal_stochasticity  = TRUE)
  test_mat <- terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
})

test_that("Network dispersal works as expected", {
  infected_file_list <-
    system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2008-03-31"
  parameter_means <- c(2, 21, 1, 500, 0, 0, 100, 1000)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  network_filename <-  system.file("extdata", "simple20x20", "segments.csv", package = "PoPS")
  anthropogenic_kernel_type <- "network"
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  data <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         start_date = start_date,
         end_date = end_date,
         anthropogenic_kernel_type = anthropogenic_kernel_type,
         network_filename = network_filename)

  test_mat <- terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
})

test_that("uncertainty propogation works as expected", {
  infected_file_list <- system.file("extdata", "simple20x20", "infected_wsd.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple20x20", "host_w_sd2.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2008-03-31"
  parameter_means <- c(2, 21, 1, 500, 0, 0, 100, 1000)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  anthropogenic_kernel_type <- "cauchy"
  use_initial_condition_uncertainty <- TRUE
  use_host_uncertainty <- TRUE
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  data <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         start_date = start_date,
         end_date = end_date,
         anthropogenic_kernel_type = anthropogenic_kernel_type,
         use_initial_condition_uncertainty = use_initial_condition_uncertainty,
         use_host_uncertainty = use_host_uncertainty)

  test_mat <- terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])

  infected_file_list <- system.file("extdata", "simple20x20", "infected_wsd.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
  total_populations_file <- system.file("extdata", "simple20x20", "all_plants.tif",
                                        package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2008-03-31"
  parameter_means <- c(2, 21, 1, 500, 0, 0, 100, 1000)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  anthropogenic_kernel_type <- "cauchy"
  use_initial_condition_uncertainty <- TRUE
  use_host_uncertainty <- FALSE

  data <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         start_date = start_date,
         end_date = end_date,
         anthropogenic_kernel_type = anthropogenic_kernel_type,
         use_initial_condition_uncertainty = use_initial_condition_uncertainty,
         use_host_uncertainty = use_host_uncertainty)

  test_mat <- terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])

  infected_file_list <-
    system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple20x20", "host_w_sd2.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2008-03-31"
  parameter_means <- c(2, 21, 1, 500, 0, 0, 100, 1000)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  anthropogenic_kernel_type <- "cauchy"
  use_initial_condition_uncertainty <- FALSE
  use_host_uncertainty <- TRUE
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  data <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         start_date = start_date,
         end_date = end_date,
         anthropogenic_kernel_type = anthropogenic_kernel_type,
         use_initial_condition_uncertainty = use_initial_condition_uncertainty,
         use_host_uncertainty = use_host_uncertainty)

  test_mat <- terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
})

test_that("multiple_random seeds works and returns expected results", {
  infected_file_list <-
    system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2008-03-31"
  parameter_means <- c(5, 21, 1, 500, 0, 0, 100, 1000)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  anthropogenic_kernel_type <- "cauchy"
  multiple_random_seeds <- TRUE
  file_random_seeds <- NULL
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  data <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         start_date = start_date,
         end_date = end_date,
         anthropogenic_kernel_type = anthropogenic_kernel_type,
         multiple_random_seeds = multiple_random_seeds,
         file_random_seeds = file_random_seeds)

  test_mat <- terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])

  file_random_seeds <- system.file("extdata", "simple2x2", "randoms.csv", package = "PoPS")

  data <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         start_date = start_date,
         end_date = end_date,
         anthropogenic_kernel_type = anthropogenic_kernel_type,
         multiple_random_seeds = multiple_random_seeds,
         file_random_seeds = file_random_seeds)

  test_mat <- terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
})


test_that("Using soils returns expected results", {
  infected_file_list <-
    system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple2x2", "host.tif", package = "PoPS")
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2009-12-31"
  parameter_means <- c(5, 21, 1, 500, 0, 0, 100, 1000)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  use_soils <- TRUE
  dispersers_to_soils_percentage <- 0.05
  coefficient_file <-
    system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  data <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         start_date = start_date,
         end_date = end_date,
         temp = TRUE,
         temperature_coefficient_file = coefficient_file,
         use_soils = use_soils,
         dispersers_to_soils_percentage = dispersers_to_soils_percentage)

  test_mat <- terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
  expect_equal(length(data$soil_reservoirs[[1]]), 20)
  expect_equal(length(data$soil_reservoirs[[2]]), 20)

  data <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         start_date = start_date,
         end_date = end_date,
         temp = TRUE,
         temperature_coefficient_file = coefficient_file,
         use_soils = use_soils,
         dispersers_to_soils_percentage = dispersers_to_soils_percentage,
         soil_starting_pest_file = infected_file_list,
         start_with_soil_populations = TRUE)

  test_mat <- terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
  expect_equal(length(data$soil_reservoirs[[1]]), 20)
  expect_equal(length(data$soil_reservoirs[[2]]), 20)
})

test_that("Using multiple hosts works as expected", {
  infected_file_list <-
    c(system.file("extdata", "simple2x2", "infected_oak.tif", package = "PoPS"),
      system.file("extdata", "simple2x2", "infected_tanoak.tif", package = "PoPS"),
      system.file("extdata", "simple2x2", "infected_baylaurel.tif", package = "PoPS"))
  host_file_list <-
    c(system.file("extdata", "simple2x2", "host_oak.tif", package = "PoPS"),
      system.file("extdata", "simple2x2", "host_tanoak.tif", package = "PoPS"),
      system.file("extdata", "simple2x2", "host_baylaurel.tif", package = "PoPS"))
  total_populations_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  start_date <- "2008-01-01"
  end_date <- "2009-12-31"
  parameter_means <- c(5, 21, 1, 500, 0, 0, 100, 1000)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  coefficient_file <-
    system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  pest_host_table <-
    system.file("extdata", "pest_host_table.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_multihost.csv", package = "PoPS")

  data <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         start_date = start_date,
         end_date = end_date,
         temp = TRUE,
         temperature_coefficient_file = coefficient_file)

  test_mat <- terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
  expect_equal(length(data$soil_reservoirs[[1]]), 20)
  expect_equal(length(data$soil_reservoirs[[2]]), 20)

  data <-
    pops(infected_file_list = infected_file_list,
         host_file_list = host_file_list,
         total_populations_file = total_populations_file,
         parameter_means = parameter_means,
         parameter_cov_matrix = parameter_cov_matrix,
         pest_host_table = pest_host_table,
         competency_table = competency_table,
         start_date = start_date,
         end_date = end_date,
         temp = TRUE,
         temperature_coefficient_file = coefficient_file,
         use_soils = use_soils,
         dispersers_to_soils_percentage = dispersers_to_soils_percentage,
         soil_starting_pest_file = infected_file_list,
         start_with_soil_populations = TRUE)

  test_mat <- terra::as.matrix(terra::rast(infected_file_list), wide = TRUE)
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
  expect_equal(length(data$soil_reservoirs[[1]]), 20)
  expect_equal(length(data$soil_reservoirs[[2]]), 20)
})
