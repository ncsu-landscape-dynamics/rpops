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
  expect_equal(length(unique(test$group)), 4)
  expect_equal(max(test$group_size), 1)

  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[3]]
  test <- get_all_infected(infected, direction = 4)
  expect_equal(nrow(test), sum(infected[infected > 0] > 0))
  expect_equal(sum(test$detections), sum(infected[infected >0]))
  expect_equal(length(unique(test$group)), 3)
  expect_equal(max(test$group_size), 16)

  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[3]]
  test <- get_all_infected(infected, direction = 8)
  expect_equal(nrow(test), sum(infected[infected > 0] > 0))
  expect_equal(sum(test$detections), sum(infected[infected >0]))
  expect_equal(length(unique(test$group)), 2)
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

test_that("Get all infected returns all infected locations", {
  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[1]]
  distances <- get_infection_distances(infected, method = "Foci", points = c())
  expect_equal(distances$distance, 0)
  expect_equal(nrow(distances), 1)

  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[1]]
  distances <- get_infection_distances(infected, method = "Border", points = c())
  expect_equal(distances$distance, 0)
  expect_equal(nrow(distances), 1)

  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[2]]
  distances <- get_infection_distances(infected, method = "Foci", points = c())
  expect_equal(nrow(distances), 4)
  expect_equal(all(distances$distance > 1), TRUE)

  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[2]]
  distances <- get_infection_distances(infected, method = "Border", points = c())
  expect_equal(nrow(distances), 4)
  expect_equal(all(distances$distance == 0), TRUE)

  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[3]]
  distances <- get_infection_distances(infected, method = "Foci", points = c())
  expect_equal(nrow(distances), sum(infected[infected > 0] > 0))
  expect_equal(all(distances$distance > 0), TRUE)

  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[3]]
  distances <- get_infection_distances(infected, method = "Border", points = c())
  expect_equal(nrow(distances), sum(infected[infected > 0] > 0))
  expect_equal(nrow(distances[distances$distance == 0,]), 22)
  expect_equal(nrow(distances[distances$distance > 0,]), 3)

  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[3]]
  distances <-
    get_infection_distances(infected,
                            method = "Points",
                            points = data.frame(i = 1, j = 1))
  expect_equal(nrow(distances), sum(infected[infected > 0] > 0))
  expect_equal(nrow(distances[distances$distance > 0,]), 25)

})

test_that("Raster mean and sd returns a raster from the mean and sd", {

  host_file <-
    system.file("extdata", "simple20x20", "host_w_sd.tif", package = "PoPS")
  host <- stack(host_file)
  hosts <- output_from_raster_mean_and_sd(host)
  hosts2 <- output_from_raster_mean_and_sd(host)
  expect_true(any(hosts[hosts > host[[1]]] > 0))
  expect_true(any(hosts[hosts == host[[1]]] > 0))
  expect_true(any(hosts[hosts < host[[1]]] > 0))
  expect_true(any(hosts2[hosts2 > host[[1]]] > 0))
  expect_true(any(hosts2[hosts2 == host[[1]]] > 0))
  expect_true(any(hosts2[hosts2 < host[[1]]] > 0))

})

# treatment_auto

test_that("Automated treatment location selection", {

  host_file <-
    system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
  host <- stack(host_file)
  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected <- infected[[3]]
  number_of_locations <- 2

  treatments <- treatment_auto(infected,
                               host,
                               method = "Foci",
                               priority = "group size",
                               number_of_locations = number_of_locations,
                               points = c(),
                               treatment_efficacy = 1,
                               buffer_cells = 1.5,
                               direction_first = TRUE,
                               treatment_priority = "equal",
                               treatment_rank = c(0))

  number_of_locations <- 5
  treatments <- treatment_auto(infected,
                               host,
                               method = "Points",
                               priority = "host",
                               number_of_locations = number_of_locations,
                               points = data.frame(i = 1, j = 1),
                               treatment_efficacy = 1,
                               buffer_cells = 1.5,
                               direction_first = TRUE,
                               treatment_priority = "equal",
                               treatment_rank = c(0))


  host_file <-
    system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
  host <- stack(host_file)
  infected_file <-
    system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
  infected <- stack(infected_file)
  infected1 <- infected[[3]]
  infected2 <- infected[[2]]
  infecteds <- stack(infected1, infected2)
  hosts <- stack(host, host)
  number_of_locations <- 7

  treatments <- treatment_auto(infecteds,
                               hosts,
                               method = "Border",
                               priority = "infected",
                               number_of_locations = number_of_locations,
                               points = c(),
                               treatment_efficacy = 1,
                               buffer_cells = 1,
                               direction_first = TRUE,
                               treatment_priority = "ranked",
                               treatment_rank = c(2,1))
  expect_equal()

})

