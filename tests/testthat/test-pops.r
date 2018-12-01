context("test-pops")

test_that("Model stops if files don't exist", {
  expect_equal(pops(infected_file = "", host_file = "C:/Users/Chris/Desktop/rpops/data/simple2x2/total_plants.tif", total_plants_file = "C:/Users/Chris/Desktop/rpops/data/simple2x2/total_plants.tif"), "Infected file does not exist")
  expect_equal(pops(infected_file = "C:/Users/Chris/Desktop/rpops/data/simple2x2/infected.tif", host_file = "", total_plants_file = "C:/Users/Chris/Desktop/rpops/data/simple2x2/total_plants.tif"), "Host file does not exist")
  expect_equal(pops(infected_file = "C:/Users/Chris/Desktop/rpops/data/simple2x2/infected.tif", host_file = "C:/Users/Chris/Desktop/rpops/data/simple2x2/total_plants.tif", total_plants_file = ""), "Total plants file does not exist")
  
  })
