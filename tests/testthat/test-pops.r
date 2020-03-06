context("test-pops")

test_that("Model stops if files don't exist or aren't the correct extension", {
  infected_file = system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file = system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file = system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date = "2008-01-01"
  end_date = "2010-12-31"
  
  expect_equal(pops(infected_file = "", host_file =  host_file, total_plants_file =  host_file), "file does not exist")
  expect_equal(pops(infected_file =  system.file("extdata", "simple2x2", "infected.csv", package = "PoPS"), host_file =  host_file, total_plants_file =  host_file), "file is not one of '.grd', '.tif', '.img'")
  expect_equal(pops(infected_file =  infected_file, host_file = "", total_plants_file =  host_file), "file does not exist")
  expect_equal(pops(infected_file =  infected_file, host_file =  system.file("extdata", "simple2x2", "infected.csv", package = "PoPS"), total_plants_file =  host_file), "file is not one of '.grd', '.tif', '.img'")
  expect_equal(pops(infected_file =  infected_file, host_file =  host_file, total_plants_file = ""), "file does not exist")
  expect_equal(pops(infected_file =  infected_file, host_file =  host_file, total_plants_file =  system.file("extdata", "simple2x2", "infected.csv", package = "PoPS")), "file is not one of '.grd', '.tif', '.img'")
  expect_equal(pops(infected_file =  infected_file, host_file =  host_file, total_plants_file =  host_file, use_lethal_temperature = TRUE, temperature_file = ""), "file does not exist")
  expect_equal(pops(infected_file =  infected_file, host_file =  host_file, total_plants_file =  host_file, use_lethal_temperature = TRUE, temperature_file =  system.file("extdata", "simple2x2", "infected.csv", package = "PoPS")), "file is not one of '.grd', '.tif', '.img'")
  expect_equal(pops(infected_file =  infected_file, host_file =  host_file, total_plants_file =  host_file, temp = TRUE, temperature_coefficient_file = ""), "file does not exist")
  expect_equal(pops(infected_file =  infected_file, host_file =  host_file, total_plants_file =  host_file, temp = TRUE, temperature_coefficient_file =  system.file("extdata", "simple2x2", "infected.csv", package = "PoPS")), "file is not one of '.grd', '.tif', '.img'")
  expect_equal(pops(infected_file =  infected_file, host_file =  host_file, total_plants_file =  host_file, precip = TRUE, precipitation_coefficient_file = ""), "file does not exist")
  expect_equal(pops(infected_file =  infected_file, host_file =  host_file, total_plants_file =  host_file, precip = TRUE, precipitation_coefficient_file =  system.file("extdata", "simple2x2", "infected.csv", package = "PoPS")), "file is not one of '.grd', '.tif', '.img'")
  expect_equal(pops(infected_file =  infected_file, host_file =  host_file, total_plants_file =  host_file, management = TRUE, treatments_file = ""), "file does not exist")
  expect_equal(pops(infected_file =  infected_file, host_file =  host_file, total_plants_file =  host_file, management = TRUE, treatments_file =  system.file("extdata", "simple2x2", "infected.csv", package = "PoPS")), "file is not one of '.grd', '.tif', '.img'")
  
})

test_that("Model stops if parameters are of the wrong type and/or dimension", {
  infected_file = system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file = system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file = system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date = "2008-01-01"
  end_date = "2010-12-31"
  
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, time_step = "two"), "Time step must be one of 'week', 'month' or 'day'")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, end_date = "two"), "End time and/or start time not of type numeric and/or in format YYYY-MM-DD")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, end_date = 156), "End time and/or start time not of type numeric and/or in format YYYY-MM-DD")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, start_date = "five"), "End time and/or start time not of type numeric and/or in format YYYY-MM-DD")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, start_date = 19), "End time and/or start time not of type numeric and/or in format YYYY-MM-DD")

})

test_that("Input raster resolutions, extents, and crs all match", {
  infected_file = system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file = system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file = system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date = "2008-01-01"
  end_date = "2010-12-31"
  
  expect_equal(pops(infected_file = infected_file, host_file = system.file("extdata", "simple5x5", "total_plants.tif", package = "PoPS"), total_plants_file = host_file), "Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = system.file("extdata", "simple5x5", "total_plants.tif", package = "PoPS")), "Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
  expect_equal(pops(infected_file = system.file("extdata", "simple5x5", "total_plants.tif", package = "PoPS"), host_file = host_file, total_plants_file = host_file), "Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, use_lethal_temperature = TRUE, temperature_file = system.file("extdata", "simple2x2", "critical_temp_diff_extent.tif", package = "PoPS")), "Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, temp = TRUE, temperature_coefficient_file = system.file("extdata", "simple2x2", "critical_temp_diff_extent.tif", package = "PoPS")), "Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, precip = TRUE, precipitation_coefficient_file = system.file("extdata", "simple2x2", "critical_temp_diff_extent.tif", package = "PoPS")), "Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, temp = TRUE, temperature_coefficient_file = system.file("extdata", "simple2x2", "critical_temp.tif", package = "PoPS"), precip = TRUE, precipitation_coefficient_file = system.file("extdata", "simple2x2", "critical_temp_diff_extent.tif", package = "PoPS")), "Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, management = TRUE, treatments_file = system.file("extdata", "simple2x2", "critical_temp_diff_extent.tif", package = "PoPS")), "Extents of input rasters do not match. Ensure that all of your input rasters have the same extent")
  
  expect_equal(pops(infected_file = infected_file, host_file = system.file("extdata", "simple2x2", "total_plants_diff_res.tif", package = "PoPS"), total_plants_file = host_file), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = system.file("extdata", "simple2x2", "total_plants_diff_xres.tif", package = "PoPS"), total_plants_file = host_file), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = system.file("extdata", "simple2x2", "total_plants_diff_yres.tif", package = "PoPS"), total_plants_file = host_file), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = system.file("extdata", "simple2x2", "total_plants_diff_res.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = system.file("extdata", "simple2x2", "total_plants_diff_xres.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = system.file("extdata", "simple2x2", "total_plants_diff_yres.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = system.file("extdata", "simple2x2", "total_plants_diff_res.tif", package = "PoPS"), host_file = host_file, total_plants_file = host_file), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = system.file("extdata", "simple2x2", "total_plants_diff_xres.tif", package = "PoPS"), host_file = host_file, total_plants_file = host_file), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = system.file("extdata", "simple2x2", "total_plants_diff_yres.tif", package = "PoPS"), host_file = host_file, total_plants_file = host_file), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, use_lethal_temperature = TRUE, temperature_file = system.file("extdata", "simple2x2", "critical_temp_diff_res.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, use_lethal_temperature = TRUE, temperature_file = system.file("extdata", "simple2x2", "critical_temp_diff_xres.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, use_lethal_temperature = TRUE, temperature_file = system.file("extdata", "simple2x2", "critical_temp_diff_yres.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, temp = TRUE, temperature_coefficient_file = system.file("extdata", "simple2x2", "critical_temp_diff_res.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, temp  = TRUE, temperature_coefficient_file = system.file("extdata", "simple2x2", "critical_temp_diff_xres.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, temp = TRUE, temperature_coefficient_file = system.file("extdata", "simple2x2", "critical_temp_diff_yres.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, precip = TRUE, precipitation_coefficient_file = system.file("extdata", "simple2x2", "critical_temp_diff_res.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, precip  = TRUE, precipitation_coefficient_file = system.file("extdata", "simple2x2", "critical_temp_diff_xres.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, precip = TRUE, precipitation_coefficient_file = system.file("extdata", "simple2x2", "critical_temp_diff_yres.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, temp = TRUE, temperature_coefficient_file = system.file("extdata", "simple2x2", "critical_temp.tif", package = "PoPS"), precip = TRUE, precipitation_coefficient_file = system.file("extdata", "simple2x2", "critical_temp_diff_res.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, temp = TRUE, temperature_coefficient_file = system.file("extdata", "simple2x2", "critical_temp.tif", package = "PoPS"), precip = TRUE, precipitation_coefficient_file = system.file("extdata", "simple2x2", "critical_temp_diff_xres.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, temp = TRUE, temperature_coefficient_file = system.file("extdata", "simple2x2", "critical_temp.tif", package = "PoPS"), precip = TRUE, precipitation_coefficient_file = system.file("extdata", "simple2x2", "critical_temp_diff_yres.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, management = TRUE, treatments_file = system.file("extdata", "simple2x2", "critical_temp_diff_res.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, management = TRUE, treatments_file = system.file("extdata", "simple2x2", "critical_temp_diff_xres.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, management = TRUE, treatments_file = system.file("extdata", "simple2x2", "critical_temp_diff_yres.tif", package = "PoPS")), "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
  
  expect_equal(pops(infected_file = infected_file, host_file = system.file("extdata", "simple2x2", "total_plants_with_crs.tif", package = "PoPS"), total_plants_file = host_file), "Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = system.file("extdata", "simple2x2", "total_plants_with_crs.tif", package = "PoPS")), "Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
  expect_equal(pops(infected_file = system.file("extdata", "simple2x2", "total_plants_with_crs.tif", package = "PoPS"), host_file = host_file, total_plants_file = host_file), "Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, use_lethal_temperature = TRUE, temperature_file = system.file("extdata", "simple2x2", "critical_temp_diff_crs.tif", package = "PoPS")), "Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, temp = TRUE, temperature_coefficient_file = system.file("extdata", "simple2x2", "critical_temp_diff_crs.tif", package = "PoPS")), "Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, precip = TRUE, precipitation_coefficient_file = system.file("extdata", "simple2x2", "critical_temp_diff_crs.tif", package = "PoPS")), "Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, temp = TRUE, temperature_coefficient_file = system.file("extdata", "simple2x2", "critical_temp.tif", package = "PoPS"), precip = TRUE, precipitation_coefficient_file = system.file("extdata", "simple2x2", "critical_temp_diff_crs.tif", package = "PoPS")), "Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, management = TRUE, treatments_file = system.file("extdata", "simple2x2", "critical_temp_diff_crs.tif", package = "PoPS")), "Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
  
})

test_that("Infected results return initial infected if reproductive rate is set to 0", {
  infected_file = system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file = system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file = system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date = "2008-01-01"
  end_date = "2010-12-31"
  
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.0)$infected[[1]], raster::as.matrix(raster::raster(infected_file)))
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, use_lethal_temperature = TRUE, temperature_file = system.file("extdata", "simple2x2", "critical_temp.tif", package = "PoPS"), reproductive_rate = 0.0)$infected[[1]], raster::as.matrix(raster::raster(infected_file)))
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.0, temp = TRUE, temperature_coefficient_file = coefficient_file)[[1]][[1]], raster::as.matrix(raster::raster(infected_file)))
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.0, precip = TRUE, precipitation_coefficient_file = coefficient_file)[[1]][[1]], raster::as.matrix(raster::raster(infected_file)))
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.0, temp = TRUE, temperature_coefficient_file = coefficient_file, precip = TRUE, precipitation_coefficient_file = coefficient_file)[[1]][[1]], raster::as.matrix(raster::raster(infected_file)))
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.0, use_lethal_temperature = TRUE, temperature_file = system.file("extdata", "simple2x2", "critical_temp.tif", package = "PoPS"), temp = TRUE, temperature_coefficient_file = coefficient_file, precip = TRUE, precipitation_coefficient_file = coefficient_file)[[1]][[1]], raster::as.matrix(raster::raster(infected_file)))

  # expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.0, temp = TRUE, temperature_coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient_weeks.tif", package = "PoPS"), time_step = "week")[[1]][[1]], raster::as.matrix(raster::raster(infected_file)))
  # expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.0, precip = TRUE, precipitation_coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient_weeks.tif", package = "PoPS"), time_step = "week")[[1]][[1]], raster::as.matrix(raster::raster(infected_file)))
  # expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.0, temp = TRUE, temperature_coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient_weeks.tif", package = "PoPS"), precip = TRUE, precipitation_coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient_weeks.tif", package = "PoPS"), time_step = "week")[[1]][[1]], raster::as.matrix(raster::raster(infected_file)))
  # expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.0, use_lethal_temperature = TRUE, temperature_file = system.file("extdata", "simple2x2", "critical_temp.tif", package = "PoPS"), temp = TRUE, temperature_coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient_weeks.tif", package = "PoPS"), precip = TRUE, precipitation_coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient_weeks.tif", package = "PoPS"), time_step = "week")[[1]][[1]], raster::as.matrix(raster::raster(infected_file)))
  # 
  # expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.0, temp = TRUE, temperature_coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient_days.tif", package = "PoPS"))[[1]][[1]], raster::as.matrix(raster::raster(infected_file, time_step = "day")))
  # expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.0, precip = TRUE, precipitation_coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient_days.tif", package = "PoPS"))[[1]][[1]], raster::as.matrix(raster::raster(infected_file, time_step = "day")))
  # expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.0, temp = TRUE, temperature_coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient_days.tif", package = "PoPS"), precip = TRUE, precipitation_coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient_days.tif", package = "PoPS"), time_step = "day")[[1]][[1]], raster::as.matrix(raster::raster(infected_file)))
  # expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.0, use_lethal_temperature = TRUE, temperature_file = system.file("extdata", "simple2x2", "critical_temp.tif", package = "PoPS"), temp = TRUE, temperature_coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient_days.tif", package = "PoPS"), precip = TRUE, precipitation_coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient_days.tif", package = "PoPS"), time_step = "day")[[1]][[1]], raster::as.matrix(raster::raster(infected_file)))

})

test_that("Infected results returns all 0's if minimum temp drops below lethal temperature", {
  infected_file = system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file = system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file = system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date = "2008-01-01"
  end_date = "2010-12-31"
  
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, use_lethal_temperature = TRUE, temperature_file = temperature_file)$infected[[1]], matrix(0,ncol = 2, nrow = 2))
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, use_lethal_temperature = TRUE, temperature_file = temperature_file, precip = TRUE, precipitation_coefficient_file = coefficient_file)$infected[[1]], matrix(0,ncol = 2, nrow = 2))
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, use_lethal_temperature = TRUE, temperature_file = temperature_file, temp = TRUE, temperature_coefficient_file = coefficient_file)$infected[[1]], matrix(0,ncol = 2, nrow = 2))
  expect_equal(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, use_lethal_temperature = TRUE, temperature_file = temperature_file, temp = TRUE, temperature_coefficient_file = coefficient_file, precip = TRUE, precipitation_coefficient_file = coefficient_file)$infected[[1]], matrix(0,ncol = 2, nrow = 2))

  })

test_that("Infected and Susceptible results return all 0's if treatments file is all 1's", {
  infected_file = system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file = system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file = system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date = "2008-01-01"
  end_date = "2010-12-31"
  treatments_file = system.file("extdata", "simple2x2", "treatments.tif", package = "PoPS")
  
  data <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, management  = TRUE, treatment_dates = c("2008-12-01"), treatments_file = treatments_file, reproductive_rate = 1.0, start_date = start_date, end_date = end_date)
  
  expect_equal(data$infected[[1]], matrix(0,ncol = 2, nrow = 2))
  expect_equal(data$susceptible[[1]], matrix(0,ncol = 2, nrow = 2))
  
})

test_that("Infected results are greater than initial infected", {
  infected_file = system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file = system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file = system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date = "2008-01-01"
  end_date = "2010-12-31"
  
  expect_equal(all(pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file)[[1]][[1]] >= raster::as.matrix(raster::raster(infected_file))), TRUE)
  expect_equal(all(pops(infected_file = infected_file, host_file = system.file("extdata", "simple2x2", "total_plants_host_greater_than_infected.tif", package = "PoPS"), total_plants_file = host_file)[[1]][[1]] >= raster::as.matrix(raster::raster(infected_file))), TRUE)

})

test_that("Susceptibles are never negative", {
  infected_file = system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file = system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file = system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date = "2008-01-01"
  end_date = "2010-12-31"
  
  data <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.4, random_seed = 42, start_date = start_date, end_date = end_date)
  
  expect_equal(all(data[[2]][[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data[[2]][[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data[[2]][[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)

  data <- pops(infected_file = infected_file, host_file = system.file("extdata", "simple2x2", "total_plants_host_greater_than_infected.tif", package = "PoPS"), total_plants_file = host_file, reproductive_rate = 0.5, random_seed = 42, start_date = start_date, end_date = end_date)
  
  expect_equal(all(data[[2]][[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data[[2]][[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data[[2]][[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)

})

test_that("Infected results with weather are less than those without weather", {
  infected_file = system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file = system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file = system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date = "2008-01-01"
  end_date = "2010-12-31"
  
  data <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.4, random_seed = 42, start_date = start_date, end_date = end_date)
  data_temp <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, temp = TRUE, temperature_coefficient_file = coefficient_file, reproductive_rate = 0.4, random_seed = 42, start_date = start_date, end_date = end_date)
  data_precip <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, precip = TRUE, precipitation_coefficient_file = coefficient_file, reproductive_rate = 0.4, random_seed = 42, start_date = start_date, end_date = end_date)
  data_weather <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, temp = TRUE, temperature_coefficient_file = coefficient_file, precip = TRUE, precipitation_coefficient_file = coefficient_file, reproductive_rate = 0.4, random_seed = 42, start_date = start_date, end_date = end_date)
  
  expect_gte(sum(data[[1]][[1]]), sum(data_temp[[1]][[1]]))
  expect_gte(sum(data[[1]][[2]]), sum(data_temp[[1]][[2]]))
  expect_gte(sum(data[[1]][[3]]), sum(data_temp[[1]][[3]]))

  expect_gte(sum(data[[1]][[1]]), sum(data_precip[[1]][[1]]))
  expect_gte(sum(data[[1]][[2]]), sum(data_precip[[1]][[2]]))
  expect_gte(sum(data[[1]][[3]]), sum(data_precip[[1]][[3]]))

  expect_gte(sum(data[[1]][[1]]), sum(data_weather[[1]][[1]]))
  expect_gte(sum(data[[1]][[2]]), sum(data_weather[[1]][[2]]))
  expect_gte(sum(data[[1]][[3]]), sum(data_weather[[1]][[3]]))

})

test_that("Infected results are greater with same parameters for weekly spread vs. monthly", {
  infected_file = system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file = system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file = system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date = "2008-01-01"
  end_date = "2010-12-31"
  
  data_week <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.2, time_step = "week", random_seed = 42, start_date = start_date, end_date = end_date)
  data_month <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.2, time_step = "month", random_seed = 42, start_date = start_date, end_date = end_date)
  
  expect_equal(all(data_week[[1]][[1]] >= data_month[[1]][[1]]), TRUE)
  expect_equal(all(data_week[[1]][[2]] >= data_month[[1]][[2]]), TRUE)

})

test_that("Infected results are greater with same parameters for daily spread vs. monthly and weekly", {
  infected_file = system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file = system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file = system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date = "2008-01-01"
  end_date = "2010-12-31"
  
  
  data_day <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.1, time_step = "day", random_seed = 42)
  data_week <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.1, time_step = "week", random_seed = 42, start_date = start_date, end_date = end_date)
  data_month <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.1, time_step = "month", random_seed = 42, start_date = start_date, end_date = end_date)
  
  expect_equal(all(data_day[[1]][[1]] >= data_month[[1]][[1]]), TRUE)

  expect_equal(all(data_day[[1]][[1]] >= data_week[[1]][[1]]), TRUE)
})

test_that("Infected results are greater without treatment than with treatment",{
  infected_file = system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file = system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  treatments_file = system.file("extdata", "simple2x2", "treatments_1_1.tif", package = "PoPS")
  treatment_dates = c("2008-03-05")
  start_date = "2008-01-01"
  end_date = "2010-12-31"
  
  data <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.8, random_seed = 44, start_date = start_date, end_date = end_date)
  data_treat <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, management = TRUE, treatment_dates = treatment_dates, treatments_file = treatments_file, reproductive_rate = 0.8, random_seed = 44, start_date = start_date, end_date = end_date)
  
  expect_equal(all(data[[1]][[1]] >= data_treat[[1]][[1]]), TRUE)
  expect_equal(all(data[[1]][[2]] >= data_treat[[1]][[2]]), TRUE)
})

test_that("Infected results are greater with higher reproductive rate", {
  infected_file = system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file = system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  start_date = "2008-01-01"
  end_date = "2010-12-31"
  
  data_1 <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 1.00, time_step = "month", random_seed = 42, start_date = start_date, end_date = end_date)
  data_075 <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.75, time_step = "month", random_seed = 42, start_date = start_date, end_date = end_date)
  data_050 <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.50, time_step = "month", random_seed = 42, start_date = start_date, end_date = end_date)
  data_025 <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.25, time_step = "month", random_seed = 42, start_date = start_date, end_date = end_date)
  data_010 <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, reproductive_rate = 0.10, time_step = "month", random_seed = 42, start_date = start_date, end_date = end_date)
  
  expect_gte(sum(data_1[[1]][[1]]), sum(data_075[[1]][[1]]))

  expect_gte(sum(data_050[[1]][[1]]), sum(data_025[[1]][[1]]))
  expect_gte(sum(data_050[[1]][[2]]), sum(data_025[[1]][[2]]))

  expect_gte(sum(data_025[[1]][[1]]), sum(data_010[[1]][[1]]))
  expect_gte(sum(data_025[[1]][[2]]), sum(data_010[[1]][[2]]))

})

test_that("Treatments apply no matter what time step", {
  infected_file = system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  host_file = system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  coefficient_file = system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
  temperature_file = system.file("extdata", "simple2x2", "critical_temp_all_below_threshold.tif", package = "PoPS")
  start_date = "2009-01-01"
  end_date = "2009-12-31"
  treatments_file = system.file("extdata", "simple2x2", "treatments.tif", package = "PoPS")
  dates <- seq.Date(as.Date(start_date), as.Date(end_date), by = "days")
  for (i in 1:length(dates)) {
    data <- pops(infected_file = infected_file, host_file = host_file, total_plants_file = host_file, management  = TRUE, treatment_dates = c( as.character(dates[i])), treatments_file = treatments_file, reproductive_rate = 0, start_date = start_date, end_date = end_date)
    expect_equal(data$infected[[1]], matrix(0,ncol = 2, nrow = 2))
    expect_equal(data$susceptible[[1]], matrix(0,ncol = 2, nrow = 2))
  }
})
