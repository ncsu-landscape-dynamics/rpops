context("test-helpers")

test_that("Raster mean and sd returns a raster from the mean and sd", {

  host_file <-
    system.file("extdata", "simple20x20", "host_w_sd.tif", package = "PoPS")
  host <- terra::rast(host_file)
  hosts <- output_from_raster_mean_and_sd(host)
  hosts2 <- output_from_raster_mean_and_sd(host)
  expect_true(any(hosts[hosts > host[[1]]] > 0))
  expect_true(any(hosts[hosts == host[[1]]] > 0))
  expect_true(any(hosts[hosts < host[[1]]] > 0))
  expect_true(any(hosts2[hosts2 > host[[1]]] > 0))
  expect_true(any(hosts2[hosts2 == host[[1]]] > 0))
  expect_true(any(hosts2[hosts2 < host[[1]]] > 0))

})
