context("test-helpers")

test_that("Get all infected returns all infected locations", {
  infected_file <-
    system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  infected <- raster(infected_file)
  test <- get_all_infected(infected, direction = 4)
  expect_equal(nrow(test), 1)
  expect_equal(test$cell_number, 1)
  expect_equal(test$detections, 5)
  expect_equal(test$i, 1)
  expect_equal(test$j, 1)
  expect_equal(test$group, 1)
  expect_equal(test$group_size, 1)

  infected_file <-
    system.file("extdata", "simple20x20", "infected_single.tif", package = "PoPS")
  infected <- raster(infected_file)
  test <- get_all_infected(infected, direction = 4)
  expect_equal(nrow(test), 393)
  expect_equal(unique(test$group), 1)
  expect_equal(unique(test$group_size), 393)

  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[1]]
  test <- get_all_infected(infected, direction = 4)
  expect_equal(nrow(test), 1)
  expect_equal(test$cell_number, 1)
  expect_equal(test$detections, 4)
  expect_equal(test$i, 1)
  expect_equal(test$j, 1)
  expect_equal(test$group, 1)
  expect_equal(test$group_size, 1)

  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[2]]
  test <- get_all_infected(infected, direction = 4)
  expect_equal(nrow(test), 4)
  expect_equal(unique(test$group), 4)
  expect_equal(max(test$group_size), 1)

  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[3]]
  test <- get_all_infected(infected, direction = 4)
  expect_equal(nrow(test), sum(infected[infected > 0] > 0))
  expect_equal(sum(test$detections), sum(infected[infected >0]))
  expect_equal(unique(test$group), 3)
  expect_equal(max(test$group_size), 16)

  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[3]]
  test <- get_all_infected(infected, direction = 8)
  expect_equal(nrow(test), sum(infected[infected > 0] > 0))
  expect_equal(sum(test$detections), sum(infected[infected >0]))
  expect_equal(unique(test$group), 2)
  expect_equal(max(test$group_size), 16)

})


test_that("Get all infected returns all infected locations", {
  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[1]]
  foci <- get_foci(infected)
  expect_equal(foci$i, 1)
  expect_equal(foci$j, 1)

  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[2]]
  foci <- get_foci(infected)
  expect_equal(foci$i, 5)
  expect_equal(foci$j, 3)

  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[3]]
  foci <- get_foci(infected)
  expect_equal(foci$i, 9)
  expect_equal(foci$j, 5)

})

# treatment_auto

# get_infection_distances

# get_infection_border
test_that("Get all infected returns all infected locations", {
  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[1]]
  border <- get_infection_border(infected)
  expect_equal(border$i, 1)
  expect_equal(border$j, 1)
  expect_equal(nrow(border), 1)

  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[2]]
  border <- get_infection_border(infected)
  expect_equal(nrow(border), 4)

  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[3]]
  border <- get_infection_border(infected)
  expect_equal(nrow(border), 22)

})





# output_from_raster_mean_and_sd