context("test-quantity_allocation_disagreement")

test_that(
  "Quantity allocation is 0 when comparison and reference are the exact
  same!", {
  comparison <- terra::rast(matrix(1, nrow = 2, ncol = 2), crs = "epsg:4326")
  reference <- terra::rast(matrix(1, nrow = 2, ncol = 2), crs = "epsg:4326")
  use_configuration <- TRUE
  mask <- NULL
  use_distance <- TRUE
  data <-
    quantity_allocation_disagreement(reference, comparison, use_configuration, mask, use_distance)

  expect_equal(data$quantity_disagreement, 0)
  expect_equal(data$allocation_disagreement, 0)
  expect_equal(data$total_disagreement, 0)
  expect_equal(data$false_negative, 0)
  expect_equal(data$false_positive, 0)
  expect_equal(data$odds_ratio, 0)
  expect_equal(data$accuracy, 1)
  expect_equal(data$precision, 1)
  expect_equal(data$recall, 1)
  expect_equal(data$specificity, 0)
  expect_equal(data$configuration_disagreement, 0)
  expect_equal(data$distance_difference, 0)
})

test_that("Check that quantity disagreement, total disagreement, ommision, and
          directional disagreement are 4 when ref is all 1's and comp is all
          0's!", {
  comp <- terra::rast(matrix(0, nrow = 2, ncol = 2), crs = "epsg:4326")
  ref <- terra::rast(matrix(1, nrow = 2, ncol = 2), crs = "epsg:4326")
  comparison <- comp
  reference <- ref
  use_configuration <- TRUE
  mask <- NULL
  use_distance <- FALSE
  data <- quantity_allocation_disagreement(ref, comp, use_configuration, mask, use_distance)
  expect_equal(data$quantity_disagreement, 4)
  expect_equal(data$allocation_disagreement, 0)
  expect_equal(data$total_disagreement, 4)
  expect_equal(data$false_negative, 4)
  expect_equal(data$false_positive, 0)
  expect_equal(data$odds_ratio, 0)
  expect_equal(data$accuracy, 0)
  expect_equal(data$precision, 0)
  expect_equal(data$recall, 0)
  expect_equal(data$specificity, 0)
  expect_equal(data$configuration_disagreement, 0.75)
  expect_equal(data$distance_difference, 0)
})

test_that("Check that quantity disagreement, total disagreement,
          number_of_infected_comp, and commision are 4 and directional
          disagreement is -4 when ref is all 0's and comp is all 1's!", {
  comp <- terra::rast(matrix(1, nrow = 2, ncol = 2))
  ref <- terra::rast(matrix(0, nrow = 2, ncol = 2))
  comparison <- comp
  reference <- ref
  use_configuration <- FALSE
  mask <- NULL
  data <- quantity_allocation_disagreement(ref, comp, use_configuration, mask)
  expect_equal(data$quantity_disagreement, 4)
  expect_equal(data$allocation_disagreement, 0)
  expect_equal(data$total_disagreement, 4)
  expect_equal(data$false_negative, 0)
  expect_equal(data$false_positive, 4)
  expect_equal(data$odds_ratio, 0)
  expect_equal(data$accuracy, 0)
  expect_equal(data$precision, 0)
  expect_equal(data$recall, 0)
  expect_equal(data$specificity, 0)
})

test_that("Check that allocation disgreement and total disagreement are 4 and
          number_of_infected_comp, ommision and commision are 4 and directional
          disagreement is 0 when ref has 1's at [2,1] and [2,2] and comp  1's
          at [1,1] and [1,2]!", {
  comp <- terra::rast(matrix(0, nrow = 2, ncol = 2), crs = "epsg:4326")
  comp[1, 1] <- 1
  comp[1, 2] <- 1
  ref <- terra::rast(matrix(0, nrow = 2, ncol = 2), crs = "epsg:4326")
  ref[2, 1] <- 1
  ref[2, 2] <- 1
  comparison <- comp
  reference <- ref
  use_configuration <- TRUE
  mask <- NULL
  use_distance <- TRUE
  data <- quantity_allocation_disagreement(ref, comp, use_configuration, mask, use_distance)
  expect_equal(data$quantity_disagreement, 0)
  expect_equal(data$allocation_disagreement, 4)
  expect_equal(data$total_disagreement, 4)
  expect_equal(data$false_negative, 2)
  expect_equal(data$false_positive, 2)
  expect_equal(data$odds_ratio, 0)
  expect_equal(data$accuracy, 0)
  expect_equal(data$precision, 0)
  expect_equal(data$recall, 0)
  expect_equal(data$specificity, 0)
  expect_gt(data$distance_difference, 0)
  expect_equal(data$configuration_disagreement, 0)
})

test_that(
  "Check that allocation disgreement and total disagreement are 4 and
  number_of_infected_comp, ommision and commision are 4 and directional
  disagreement is 0 when ref has 1's at [2,1] and [2,2] and comp  1's at [1,1]
  and [1,2]!", {
  comp <- terra::rast(matrix(0, nrow = 2, ncol = 2))
  comp[2, 1] <- 1
  comp[1, 2] <- 1
  ref <- terra::rast(matrix(0, nrow = 2, ncol = 2))
  ref[2, 1] <- 1
  ref[2, 2] <- 1
  comparison <- comp
  reference <- ref
  use_configuration <- FALSE
  mask <- NULL
  data <- quantity_allocation_disagreement(ref, comp, use_configuration, mask)
  expect_equal(data$quantity_disagreement, 0)
  expect_equal(data$allocation_disagreement, 2)
  expect_equal(data$total_disagreement, 2)
  expect_equal(data$false_negative, 1)
  expect_equal(data$false_positive, 1)
  expect_equal(data$odds_ratio, 1)
  expect_equal(data$accuracy, 0.5)
  expect_equal(data$precision, 0.5)
  expect_equal(data$recall, 0.5)
  expect_equal(data$specificity, 0.5)
})

test_that("Check that configuration disagreement works", {
  comp <- terra::rast(matrix(0, nrow = 2, ncol = 2))
  comp[2, 1] <- 1
  comp[1, 2] <- 1
  ref <- terra::rast(matrix(0, nrow = 2, ncol = 2))
  ref[2, 1] <- 1
  ref[1, 2] <- 1
  comparison <- comp
  reference <- ref
  use_configuration <- TRUE
  mask <- NULL
  data <- quantity_allocation_disagreement(ref, comp, use_configuration, mask)
  expect_equal(data$quantity_disagreement, 0)
  expect_equal(data$allocation_disagreement, 0)
  expect_equal(data$total_disagreement, 0)
  expect_equal(data$false_negative, 0)
  expect_equal(data$false_positive, 0)
  expect_equal(data$odds_ratio, 4)
  expect_equal(data$configuration_disagreement, 0)

  comp <- terra::rast(matrix(0, nrow = 2, ncol = 2))
  comp[2, 1] <- 1
  comp[1, 2] <- 1
  ref <- terra::rast(matrix(0, nrow = 2, ncol = 2))
  ref[2, 1] <- 1
  ref[2, 2] <- 1
  mask <- terra::rast(matrix(0, nrow = 2, ncol = 2))
  mask[1, 2] <- NA
  mask[2, 2] <- NA
  comparison <- comp
  reference <- ref
  use_configuration <- TRUE

  data <- quantity_allocation_disagreement(ref, comp, use_configuration, mask)
  expect_equal(data$quantity_disagreement, 0)
  expect_equal(data$allocation_disagreement, 0)
  expect_equal(data$total_disagreement, 0)
  expect_equal(data$false_negative, 0)
  expect_equal(data$false_positive, 0)
  expect_equal(data$odds_ratio, 1)
  expect_equal(data$configuration_disagreement, 0)


  data <- quantity_allocation_disagreement(ref, comp, use_configuration, mask = NULL)
  expect_equal(data$quantity_disagreement, 0)
  expect_equal(data$allocation_disagreement, 2)
  expect_equal(data$total_disagreement, 2)
  expect_equal(data$false_negative, 1)
  expect_equal(data$false_positive, 1)
  expect_equal(data$odds_ratio, 1)
  expect_gte(data$configuration_disagreement, 0)
})
