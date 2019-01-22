context("test-quantity_allocation_disagreement")

test_that("Quantity allocation is 0 when comparison and reference are the exact same!", {
  comp <- raster(matrix(1, nrow = 2, ncol = 2))
  ref <- raster(matrix(1, nrow = 2, ncol = 2))
  data <- quantity_allocation_disagreement(ref, comp)
  expect_equal(data$quantity_disagreement, 0)
  expect_equal(data$allocation_disagreement, 0)
  expect_equal(data$total_disagreement,0)
  expect_equal(data$omission, 0)
  expect_equal(data$commission, 0)
  expect_equal(data$directional_disagreement, 0)
  expect_equal(data$number_of_infected_comp, 4)
})

test_that("Check that quantity disagreement, total disagreement, ommision, and directional disagreement are 4 when ref is all 1's and comp is all 0's!", {
  comp <- raster(matrix(0, nrow = 2, ncol = 2))
  ref <- raster(matrix(1, nrow = 2, ncol = 2))
  data <- quantity_allocation_disagreement(ref, comp)
  expect_equal(data$quantity_disagreement, 4)
  expect_equal(data$allocation_disagreement, 0)
  expect_equal(data$total_disagreement,4)
  expect_equal(data$omission, 4)
  expect_equal(data$commission, 0)
  expect_equal(data$directional_disagreement, 4)
  expect_equal(data$number_of_infected_comp, 0)
})

test_that("Check that quantity disagreement, total disagreement, number_of_infected_comp, and commision are 4 and directional disagreement is -4 when ref is all 0's and comp is all 1's!", {
  comp <- raster(matrix(1, nrow = 2, ncol = 2))
  ref <- raster(matrix(0, nrow = 2, ncol = 2))
  data <- quantity_allocation_disagreement(ref, comp)
  expect_equal(data$quantity_disagreement, 4)
  expect_equal(data$allocation_disagreement, 0)
  expect_equal(data$total_disagreement,4)
  expect_equal(data$omission, 0)
  expect_equal(data$commission, 4)
  expect_equal(data$directional_disagreement, -4)
  expect_equal(data$number_of_infected_comp, 4)
})

test_that("Check that allocation disgreement and total disagreement are 4 and number_of_infected_comp, ommision and commision are 4 and directional disagreement is 0 when ref has 1's at [2,1] and [2,2] and comp  1's at [1,1] and [1,2]!", {
  comp <- raster(matrix(0, nrow = 2, ncol = 2))
  comp[1,1] <- 1
  comp[1,2] <- 1
  ref <- raster(matrix(0, nrow = 2, ncol = 2))
  ref[2,1] <- 1
  ref[2,2] <- 1
  data <- quantity_allocation_disagreement(ref, comp)
  expect_equal(data$quantity_disagreement, 0)
  expect_equal(data$allocation_disagreement, 4)
  expect_equal(data$total_disagreement,4)
  expect_equal(data$omission, 2)
  expect_equal(data$commission, 2)
  expect_equal(data$directional_disagreement, 0)
  expect_equal(data$number_of_infected_comp, 2)
})