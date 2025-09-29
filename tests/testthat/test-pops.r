context("test-pops")

test_that("Test Set 1: Infected results return initial infected if reproductive rate is set to 0", {
  config_file <- system.file("extdata", "configs_pops", "ts1_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(pops(config)$host_pools[[1]]$infected[[1]], config$host_pools[[1]]$infected)

  config_file <- system.file("extdata", "configs_pops", "ts1_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(pops(config)$host_pools[[1]]$infected[[1]], config$host_pools[[1]]$infected)

  config_file <- system.file("extdata", "configs_pops", "ts1_config_t3.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(pops(config)$host_pools[[1]]$infected[[1]], config$host_pools[[1]]$infected)

  config_file <- system.file("extdata", "configs_pops", "ts1_config_t4.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(pops(config)$host_pools[[1]]$infected[[1]], config$host_pools[[1]]$infected)

  config_file <- system.file("extdata", "configs_pops", "ts1_config_t5.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(pops(config)$host_pools[[1]]$infected[[1]], config$host_pools[[1]]$infected)

  config_file <- system.file("extdata", "configs_pops", "ts1_config_t6.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(pops(config)$host_pools[[1]]$infected[[1]], config$host_pools[[1]]$infected)

  config_file <- system.file("extdata", "configs_pops", "ts1_config_t7.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(pops(config)$host_pools[[1]]$infected[[1]], config$host_pools[[1]]$infected)

  config_file <- system.file("extdata", "configs_pops", "ts1_config_t8.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(pops(config)$host_pools[[1]]$infected[[1]], config$host_pools[[1]]$infected)

  config_file <- system.file("extdata", "configs_pops", "ts1_config_t9.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(pops(config)$host_pools[[1]]$infected[[1]], config$host_pools[[1]]$infected)

  config_file <- system.file("extdata", "configs_pops", "ts1_config_t10.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(pops(config)$host_pools[[1]]$infected[[1]], config$host_pools[[1]]$infected)

  config_file <- system.file("extdata", "configs_pops", "ts1_config_t11.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(pops(config)$host_pools[[1]]$infected[[1]], config$host_pools[[1]]$infected)

  config_file <- system.file("extdata", "configs_pops", "ts1_config_t12.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(pops(config)$host_pools[[1]]$infected[[1]], config$host_pools[[1]]$infected)

  config_file <- system.file("extdata", "configs_pops", "ts1_config_t13.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(pops(config)$host_pools[[1]]$infected[[1]], config$host_pools[[1]]$infected)
})

test_that(
  "Test set 2: Infected results returns all 0's if minimum temp drops below lethal temperature", {
    config_file <- system.file("extdata", "configs_pops", "ts2_config_t1.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    expect_equal(pops(config)$host_pools[[1]]$infected[[1]], matrix(0, ncol = 2, nrow = 2))

    config_file <- system.file("extdata", "configs_pops", "ts2_config_t2.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    expect_equal(pops(config)$host_pools[[1]]$infected[[1]], matrix(0, ncol = 2, nrow = 2))

    config_file <- system.file("extdata", "configs_pops", "ts2_config_t3.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    expect_equal(pops(config)$host_pools[[1]]$infected[[1]], matrix(0, ncol = 2, nrow = 2))

    config_file <- system.file("extdata", "configs_pops", "ts2_config_t4.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    expect_equal(pops(config)$host_pools[[1]]$infected[[1]], matrix(0, ncol = 2, nrow = 2))
  })

test_that(
  "Test set 3: Infected results returns less infection after survival rates than before", {
    reduced_inf <- matrix(0, ncol = 2, nrow = 2)
    reduced_inf[1, 1] <- 3
    config_file <- system.file("extdata", "configs_pops", "ts3_config_t1.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    expect_equal(pops(config)$host_pools[[1]]$infected[[1]], reduced_inf)

    config_file <- system.file("extdata", "configs_pops", "ts3_config_t2.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    expect_equal(pops(config)$host_pools[[1]]$infected[[1]], reduced_inf)

    config_file <- system.file("extdata", "configs_pops", "ts3_config_t3.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    expect_equal(pops(config)$host_pools[[1]]$infected[[1]], reduced_inf)

    config_file <- system.file("extdata", "configs_pops", "ts3_config_t4.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    expect_equal(pops(config)$host_pools[[1]]$infected[[1]], reduced_inf)
  })

test_that(
  "Test set 4: Infected and Susceptible results return all 0's if treatments file is all 1's
but leaves a proportion of susceptibles if treatment method is ratio", {
  config_file <- system.file("extdata", "configs_pops", "ts4_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(data$host_pools[[1]]$infected[[1]], matrix(0, ncol = 2, nrow = 2))
  expect_equal(data$host_pools[[1]]$susceptible[[1]], matrix(0, ncol = 2, nrow = 2))

  config_file <- system.file("extdata", "configs_pops", "ts4_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(data$host_pools[[1]]$infected[[1]], matrix(0, ncol = 2, nrow = 2))
  expect_equal(data$host_pools[[1]]$susceptible[[1]], matrix(0, ncol = 2, nrow = 2))

  config_file <- system.file("extdata", "configs_pops", "ts4_config_t3.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(data$host_pools[[1]]$infected[[1]], matrix(c(2, 0, 0, 0), ncol = 2, nrow = 2))
  expect_equal(data$host_pools[[1]]$susceptible[[1]], matrix(c(6, 7, 3, 7), ncol = 2, nrow = 2))

  config_file <- system.file("extdata", "configs_pops", "ts4_config_t4.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(data$host_pools[[1]]$infected[[1]], matrix(c(0, 0, 0, 0), ncol = 2, nrow = 2))
  expect_equal(data$host_pools[[1]]$susceptible[[1]], matrix(c(6, 7, 3, 7), ncol = 2, nrow = 2))
})

test_that("Test set 5: Infected results are greater than initial infected", {
  config_file <- system.file("extdata", "configs_pops", "ts5_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(all(pops(config)$host_pools[[1]]$infected[[1]] >=
                     config$host_pools[[1]]$infected), TRUE)

  config_file <- system.file("extdata", "configs_pops", "ts5_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  expect_equal(all(pops(config)$host_pools[[1]]$infected[[1]] >=
                     config$host_pools[[1]]$infected), TRUE)
})

## this one is random
test_that("Test set 6: All kernel types lead to spread", {
  config_file <- system.file("extdata", "configs_pops", "ts6_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  infecteds <- data$host_pools[[1]]$infected[[1]]
  expect_equal(all(infecteds >= config$host_pools[[1]]$infected), TRUE)
  expect_gt(infecteds[1, 2] + infecteds[2, 1] + infecteds[2, 2], 0)

  config_file <- system.file("extdata", "configs_pops", "ts6_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  infecteds <- data$host_pools[[1]]$infected[[1]]
  expect_equal(all(infecteds >= config$host_pools[[1]]$infected), TRUE)
  expect_gt(infecteds[1, 2] + infecteds[2, 1] + infecteds[2, 2], 0)

  config_file <- system.file("extdata", "configs_pops", "ts6_config_t3.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  infecteds <- data$host_pools[[1]]$infected[[1]]
  expect_equal(all(infecteds >= config$host_pools[[1]]$infected), TRUE)
  expect_gt(infecteds[1, 2] + infecteds[2, 1] + infecteds[2, 2], 0)

  config_file <- system.file("extdata", "configs_pops", "ts6_config_t4.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  infecteds <- data$host_pools[[1]]$infected[[1]]
  expect_equal(all(infecteds >= config$host_pools[[1]]$infected), TRUE)
  expect_gt(infecteds[1, 2] + infecteds[2, 1] + infecteds[2, 2], 0)

  config_file <- system.file("extdata", "configs_pops", "ts6_config_t5.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  infecteds <- data$host_pools[[1]]$infected[[1]]
  expect_equal(all(infecteds >= config$host_pools[[1]]$infected), TRUE)
  expect_gt(infecteds[1, 2] + infecteds[2, 1] + infecteds[2, 2], 0)

  config_file <- system.file("extdata", "configs_pops", "ts6_config_t6.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  infecteds <- data$host_pools[[1]]$infected[[1]]
  expect_equal(all(infecteds >= config$host_pools[[1]]$infected), TRUE)
  expect_gt(infecteds[1, 2] + infecteds[2, 1] + infecteds[2, 2], 0)

  config_file <- system.file("extdata", "configs_pops", "ts6_config_t7.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  infecteds <- data$host_pools[[1]]$infected[[1]]
  expect_equal(all(infecteds >= config$host_pools[[1]]$infected), TRUE)
  expect_gt(infecteds[1, 2] + infecteds[2, 1] + infecteds[2, 2], 0)

  config_file <- system.file("extdata", "configs_pops", "ts6_config_t8.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  infecteds <- data$host_pools[[1]]$infected[[1]]
  expect_equal(all(infecteds >= config$host_pools[[1]]$infected), TRUE)
  expect_gte(infecteds[1, 2] + infecteds[2, 1] + infecteds[2, 2], 0)

  # checks for anthropogenic kernel type
  config_file <- system.file("extdata", "configs_pops", "ts6_config_t9.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(all(data$host_pools[[1]]$infected[[1]] >= config$host_pools[[1]]$infected), TRUE)

  config_file <- system.file("extdata", "configs_pops", "ts6_config_t9.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(all(data$host_pools[[1]]$infected[[1]] >= config$host_pools[[1]]$infected), TRUE)

  config_file <- system.file("extdata", "configs_pops", "ts6_config_t10.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(all(data$host_pools[[1]]$infected[[1]] >= config$host_pools[[1]]$infected), TRUE)

  config_file <- system.file("extdata", "configs_pops", "ts6_config_t11.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(all(data$host_pools[[1]]$infected[[1]] >= config$host_pools[[1]]$infected), TRUE)

  config_file <- system.file("extdata", "configs_pops", "ts6_config_t12.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(all(data$host_pools[[1]]$infected[[1]] >= config$host_pools[[1]]$infected), TRUE)

  config_file <- system.file("extdata", "configs_pops", "ts6_config_t13.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(all(data$host_pools[[1]]$infected[[1]] >= config$host_pools[[1]]$infected), TRUE)

  config_file <- system.file("extdata", "configs_pops", "ts6_config_t14.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(all(data$host_pools[[1]]$infected[[1]] >= config$host_pools[[1]]$infected), TRUE)

  config_file <- system.file("extdata", "configs_pops", "ts6_config_t15.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(all(data$host_pools[[1]]$infected[[1]] >= config$host_pools[[1]]$infected), TRUE)
})

test_that("Test set 7: Susceptibles are never negative", {
  config_file <- system.file("extdata", "configs_pops", "ts7_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(all(data$host_pools[[1]]$susceptible[[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data$host_pools[[1]]$susceptible[[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data$host_pools[[1]]$susceptible[[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)

  config_file <- system.file("extdata", "configs_pops", "ts7_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(all(data$host_pools[[1]]$susceptible[[1]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data$host_pools[[1]]$susceptible[[2]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
  expect_equal(all(data$host_pools[[1]]$susceptible[[3]] >= matrix(0, ncol = 2, nrow = 2)), TRUE)
})

test_that("Test set 8: SEI model works as intended", {
  config_file <- system.file("extdata", "configs_pops", "ts8_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)

  config_file <- system.file("extdata", "configs_pops", "ts8_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data2 <- pops(config)

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

  config_file <- system.file("extdata", "configs_pops", "ts8_config_t3.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data3 <- pops(config)
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

test_that("Test set 9: Infected results with weather are less than those without weather", {
  config_file <- system.file("extdata", "configs_pops", "ts9_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)

  config_file <- system.file("extdata", "configs_pops", "ts9_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data_temp <- pops(config)

  config_file <- system.file("extdata", "configs_pops", "ts9_config_t3.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data_precip <- pops(config)

  config_file <- system.file("extdata", "configs_pops", "ts9_config_t4.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data_weather <- pops(config)

  config_file <- system.file("extdata", "configs_pops", "ts9_config_t5.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data_temp_wsd <- pops(config)

  config_file <- system.file("extdata", "configs_pops", "ts9_config_t6.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data_precip_wsd <- pops(config)

  config_file <- system.file("extdata", "configs_pops", "ts9_config_t7.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data_weather_wsd <- pops(config)

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
  "Test set 10: Infected results are greater with same parameters for weekly spread vs. monthly", {
    config_file <- system.file("extdata", "configs_pops", "ts10_config_t1.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    data_week <- pops(config)
    config_file <- system.file("extdata", "configs_pops", "ts10_config_t2.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    data_month <- pops(config)
    expect_equal(all(data_week$host_pools[[1]]$infected[[1]] >=
                       data_month$host_pools[[1]]$infected[[1]]), TRUE)
    expect_equal(all(data_week$host_pools[[1]]$infected[[2]] >=
                       data_month$host_pools[[1]]$infected[[2]]), TRUE)
  })

test_that(
  "Test set 11: Infected results are greater with same parameters for daily spread vs. monthly and
  weekly", {
    config_file <- system.file("extdata", "configs_pops", "ts11_config_t1.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    data_day <- pops(config)
    config_file <- system.file("extdata", "configs_pops", "ts11_config_t2.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    data_week <- pops(config)
    config_file <- system.file("extdata", "configs_pops", "ts11_config_t3.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    data_month <- pops(config)
    expect_equal(all(data_day$host_pools[[1]]$infected[[1]] >=
                       data_month$host_pools[[1]]$infected[[1]]), TRUE)
    expect_equal(all(data_day$host_pools[[1]]$infected[[1]] >=
                       data_week$host_pools[[1]]$infected[[1]]), TRUE)
    expect_equal(all(data_week$host_pools[[1]]$infected[[1]] >=
                       data_month$host_pools[[1]]$infected[[1]]), TRUE)
  })

test_that("Test set 12: Infected results are greater without treatment than with treatment", {
  config_file <- system.file("extdata", "configs_pops", "ts12_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  config_file <- system.file("extdata", "configs_pops", "ts12_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data_treat <- pops(config)
  expect_equal(all(data$host_pools[[1]]$infected[[1]] >=
                     data_treat$host_pools[[1]]$infected[[1]]), TRUE)
  expect_equal(all(data$host_pools[[1]]$infected[[2]] >=
                     data_treat$host_pools[[1]]$infected[[2]]), TRUE)
})

test_that("Test set 13: Infected results are greater with higher reproductive rate", {
  config_file <- system.file("extdata", "configs_pops", "ts13_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data_1 <- pops(config)
  config_file <- system.file("extdata", "configs_pops", "ts13_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data_075 <- pops(config)
  config_file <- system.file("extdata", "configs_pops", "ts13_config_t3.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data_050 <- pops(config)
  config_file <- system.file("extdata", "configs_pops", "ts13_config_t4.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data_025 <- pops(config)
  config_file <- system.file("extdata", "configs_pops", "ts13_config_t5.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data_010 <- pops(config)

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

test_that("Test set 14: Treatments apply no matter what time step", {
  config_file <- system.file("extdata", "configs_pops", "ts14_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  dates <- seq.Date(as.Date(config$start_date), as.Date(config$end_date), by = "days")
  for (i in seq_len(length(dates))) {
    config$treatment_dates <- c(paste(dates[[i]]))
    data <- pops(config)
    expect_equal(data$host_pools[[1]]$infected[[1]], matrix(0, ncol = 2, nrow = 2))
    expect_equal(data$host_pools[[1]]$susceptible[[1]], matrix(0, ncol = 2, nrow = 2))
  }
})

test_that("Test set 15: Pesticide treatments apply no matter what time step", {
  config_file <- system.file("extdata", "configs_pops", "ts15_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  dates <- seq.Date(as.Date(config$start_date), as.Date(config$end_date), by = "days")
  # only go through the day in which pesticide duration end doesn't go beyond the end_date
  for (i in 1:245) {
    config$treatment_dates <- c(paste(dates[[i]]))
    data <-pops(config)
    expect_equal(data$host_pools[[1]]$infected[[1]], matrix(0, ncol = 2, nrow = 2))
    expect_equal(data$host_pools[[1]]$susceptible[[1]],
                 config$host_pools[[1]]$susceptible + config$host_pools[[1]]$infected)
  }

  config_file <- system.file("extdata", "configs_pops", "ts15_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  for (i in 1:245) {
    config$treatment_dates <- c(paste(dates[[i]]))
    data <-pops(config)
    expect_equal(data$host_pools[[1]]$infected[[1]], matrix(c(3, 0, 0, 0), ncol = 2, nrow = 2))
    expect_equal(data$host_pools[[1]]$susceptible[[1]],
                 matrix(c(14, 14, 6, 15), ncol = 2, nrow = 2))
  }
})

test_that(
  "Test set 16: Changing the output frequency returns the correct number of outputs and output
statistics", {
  config_file <- system.file("extdata", "configs_pops", "ts16_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 1)

  config_file <- system.file("extdata", "configs_pops", "ts16_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 1)

  config_file <- system.file("extdata", "configs_pops", "ts16_config_t3.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 1)

  config_file <- system.file("extdata", "configs_pops", "ts16_config_t4.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 12)

  config_file <- system.file("extdata", "configs_pops", "ts16_config_t5.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 12)

  config_file <- system.file("extdata", "configs_pops", "ts16_config_t6.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 52)

  config_file <- system.file("extdata", "configs_pops", "ts16_config_t7.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 52)

  config_file <- system.file("extdata", "configs_pops", "ts16_config_t8.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 364)

  config_file <- system.file("extdata", "configs_pops", "ts16_config_t9.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 72)
})

test_that(
  "Test set 17: Outputs occur with non-full year date range for all time step output frequency
  combinations", {
    config_file <- system.file("extdata", "configs_pops", "ts17_config_t1.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    data <- pops(config)
    expect_equal(length(data$host_pools[[1]]$infected), 1)

    config_file <- system.file("extdata", "configs_pops", "ts17_config_t2.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    data <- pops(config)
    expect_equal(length(data$host_pools[[1]]$infected), 1)

    config_file <- system.file("extdata", "configs_pops", "ts17_config_t3.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    data <- pops(config)
    expect_equal(length(data$host_pools[[1]]$infected), 1)

    config_file <- system.file("extdata", "configs_pops", "ts17_config_t4.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    data <- pops(config)
    expect_equal(length(data$host_pools[[1]]$infected), 5)

    config_file <- system.file("extdata", "configs_pops", "ts17_config_t5.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    data <- pops(config)
    expect_equal(length(data$host_pools[[1]]$infected), 5)

    config_file <- system.file("extdata", "configs_pops", "ts17_config_t6.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    data <- pops(config)
    expect_equal(length(data$host_pools[[1]]$infected), 26)

    config_file <- system.file("extdata", "configs_pops", "ts17_config_t7.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    data <- pops(config)
    expect_equal(length(data$host_pools[[1]]$infected), 26)

    config_file <- system.file("extdata", "configs_pops", "ts17_config_t8.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    data <- pops(config)
    expect_equal(length(data$host_pools[[1]]$infected), 182)

    config_file <- system.file("extdata", "configs_pops", "ts17_config_t9.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    data <- pops(config)
    expect_equal(length(data$host_pools[[1]]$infected), 182)
  })

test_that("Test set 18: Quarantine and spread rates work at all timings", {
  config_file <- system.file("extdata", "configs_pops", "ts18_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 1)
  expect_equal(length(data$quarantine_escape), 1)
  expect_equal(length(data$quarantine_escape_distance), 1)
  expect_equal(length(data$quarantine_escape_directions), 1)
  expect_equal(length(data$rates), 1)

  config_file <- system.file("extdata", "configs_pops", "ts18_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 1)
  expect_equal(length(data$quarantine_escape), 1)
  expect_equal(length(data$quarantine_escape_distance), 1)
  expect_equal(length(data$quarantine_escape_directions), 1)
  expect_equal(length(data$rates), 1)

  config_file <- system.file("extdata", "configs_pops", "ts18_config_t3.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 12)
  expect_equal(length(data$quarantine_escape), 12)
  expect_equal(length(data$quarantine_escape_distance), 12)
  expect_equal(length(data$quarantine_escape_directions), 12)
  expect_equal(length(data$rates), 12)

  config_file <- system.file("extdata", "configs_pops", "ts18_config_t4.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 12)
  expect_equal(length(data$quarantine_escape), 12)
  expect_equal(length(data$quarantine_escape_distance), 12)
  expect_equal(length(data$quarantine_escape_directions), 12)
  expect_equal(length(data$rates), 12)

  config_file <- system.file("extdata", "configs_pops", "ts18_config_t5.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 52)
  expect_equal(length(data$quarantine_escape), 52)
  expect_equal(length(data$quarantine_escape_distance), 52)
  expect_equal(length(data$quarantine_escape_directions), 52)
  expect_equal(length(data$rates), 52)

  config_file <- system.file("extdata", "configs_pops", "ts18_config_t6.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 52)
  expect_equal(length(data$quarantine_escape), 52)
  expect_equal(length(data$quarantine_escape_distance), 52)
  expect_equal(length(data$quarantine_escape_directions), 52)
  expect_equal(length(data$rates), 52)

  config_file <- system.file("extdata", "configs_pops", "ts18_config_t7.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 364)
  expect_equal(length(data$quarantine_escape), 364)
  expect_equal(length(data$quarantine_escape_distance), 364)
  expect_equal(length(data$quarantine_escape_directions), 364)
  expect_equal(length(data$rates), 364)

  config_file <- system.file("extdata", "configs_pops", "ts18_config_t8.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 364)
  expect_equal(length(data$quarantine_escape), 364)
  expect_equal(length(data$quarantine_escape_distance), 364)
  expect_equal(length(data$quarantine_escape_directions), 364)
  expect_equal(length(data$rates), 364)
})

test_that("Test set 19: Mortality works as expected", {
  config_file <- system.file("extdata", "configs_pops", "ts19_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$mortality), 12)
  expect_equal(data$host_pools[[1]]$mortality[[1]], config$host_pools[[1]]$infected)
  expect_equal(data$host_pools[[1]]$mortality[[2]], config$host_pools[[1]]$infected)
  expect_equal(data$host_pools[[1]]$mortality[[3]], config$host_pools[[1]]$infected)
  expect_equal(data$host_pools[[1]]$mortality[[4]], config$host_pools[[1]]$infected)

  config_file <- system.file("extdata", "configs_pops", "ts19_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$mortality), 12)
  expect_equal(data$host_pools[[1]]$mortality[[1]], config$host_pools[[1]]$infected)
  expect_equal(data$host_pools[[1]]$mortality[[2]], config$host_pools[[1]]$infected)
  expect_equal(data$host_pools[[1]]$mortality[[3]], config$host_pools[[1]]$infected)

  config_file <- system.file("extdata", "configs_pops", "ts19_config_t3.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$mortality), 12)
  expect_equal(data$host_pools[[1]]$mortality[[1]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[2]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[3]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[4]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[5]], config$host_pools[[1]]$infected)

  config_file <- system.file("extdata", "configs_pops", "ts19_config_t4.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$mortality), 12)
  expect_equal(data$host_pools[[1]]$mortality[[1]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[2]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[3]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[4]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[5]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[6]], matrix(0, ncol = 20, nrow = 20))
  expect_equal(data$host_pools[[1]]$mortality[[7]], config$host_pools[[1]]$infected)

  config_file <- system.file("extdata", "configs_pops", "ts19_config_t5.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
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
  expect_equal(data$host_pools[[1]]$mortality[[11]], config$host_pools[[1]]$infected)
})

test_that("Test set 20: Movements works as expected", {
  config_file <- system.file("extdata", "configs_pops", "ts20_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 12)
  expect_equal(data$host_pools[[1]]$infected[[1]], config$host_pools[[1]]$infected)
  expect_equal(data$host_pools[[1]]$infected[[2]], config$host_pools[[1]]$infected)
  expect_equal(data$host_pools[[1]]$infected[[3]], config$host_pools[[1]]$infected)
  expect_equal(data$host_pools[[1]]$infected[[4]], config$host_pools[[1]]$infected)
  infected_move <- matrix(0, ncol = 20, nrow = 20)
  infected_move[2, 1] <- 1
  expect_equal(data$host_pools[[1]]$infected[[5]], infected_move)
  sus <- config$host_pools[[1]]$susceptible
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

  config_file <- system.file("extdata", "configs_pops", "ts20_config_t3.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  expect_equal(length(data$host_pools[[1]]$infected), 12)
  expect_equal(data$host_pools[[1]]$infected[[1]], config$host_pools[[1]]$infected)
  expect_equal(data$host_pools[[1]]$infected[[2]], config$host_pools[[1]]$infected)
  expect_equal(data$host_pools[[1]]$infected[[3]], config$host_pools[[1]]$infected)
  expect_equal(data$host_pools[[1]]$infected[[4]], config$host_pools[[1]]$infected)
  infected_move <- matrix(0, ncol = 20, nrow = 20)
  infected_move[2, 1] <- 1
  expect_equal(data$host_pools[[1]]$infected[[5]], infected_move)
  sus <- config$host_pools[[1]]$susceptible
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
  "Test set 21: Overpopulation dispersal works as expected with directionality to prevent
  dispersers from leaving the simulated area", {
    config_file <- system.file("extdata", "configs_pops", "ts21_config_t1.yaml", package = "PoPS")
    config <- configuration(config_file = config_file, testing = TRUE)
    data <- pops(config)
    test_mat <- config$host_pools[[1]]$infected
    expect_lte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
    expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
    expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
    expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
  })

test_that("Test set 22: Network dispersal works as expected", {
  config_file <- system.file("extdata", "configs_pops", "ts22_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_mat <- config$host_pools[[1]]$infected
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
})

test_that("Test set 23: Network dispersal works with multiple networks", {
  config_file <- system.file("extdata", "configs_pops", "ts23_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_mat <- config$host_pools[[1]]$infected
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])

  config_file <- system.file("extdata", "configs_pops", "ts23_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_mat <- config$host_pools[[1]]$infected
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
})

test_that("Test set 24: uncertainty propogation works as expected", {
  config_file <- system.file("extdata", "configs_pops", "ts24_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_mat <- config$host_pools[[1]]$infected
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])

  config_file <- system.file("extdata", "configs_pops", "ts24_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_mat <- config$host_pools[[1]]$infected
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])

  config_file <- system.file("extdata", "configs_pops", "ts24_config_t3.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_mat <- config$host_pools[[1]]$infected
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
})

test_that("Test that 25: multiple_random seeds works and returns expected results", {
  config_file <- system.file("extdata", "configs_pops", "ts25_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_mat <- config$host_pools[[1]]$infected
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])

  config_file <- system.file("extdata", "configs_pops", "ts25_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_mat <- config$host_pools[[1]]$infected
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
})

test_that("Test set 26: Using soils returns expected results", {
  config_file <- system.file("extdata", "configs_pops", "ts26_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_mat <- config$host_pools[[1]]$infected
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
  expect_equal(length(data$soil_reservoirs[[1]]), 20)
  expect_equal(length(data$soil_reservoirs[[2]]), 20)

  config_file <- system.file("extdata", "configs_pops", "ts26_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_mat <- config$host_pools[[1]]$infected
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
  expect_equal(length(data$soil_reservoirs[[1]]), 20)
  expect_equal(length(data$soil_reservoirs[[2]]), 20)
})

test_that("Test set 27: Using multiple hosts works as expected", {
  config_file <- system.file("extdata", "configs_pops", "ts27_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_mat <- config$host_pools[[1]]$infected
  test_mat[[1]] <- 4
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
  test_mat <- config$host_pools[[2]]$infected
  expect_lte(data$host_pools[[2]]$infected[[1]][[1]], test_mat[[1]])
  expect_lte(data$host_pools[[2]]$infected[[1]][[2]], test_mat[[2]])
  expect_lte(data$host_pools[[2]]$infected[[1]][[3]], test_mat[[3]])
  expect_lte(data$host_pools[[2]]$infected[[1]][[4]], test_mat[[4]])

  config_file <- system.file("extdata", "configs_pops", "ts27_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_mat <- config$host_pools[[1]]$infected
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
  test_mat <- config$host_pools[[2]]$infected
  expect_gte(data$host_pools[[2]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[2]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[2]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[2]]$infected[[1]][[4]], test_mat[[4]])

  config_file <- system.file("extdata", "configs_pops", "ts27_config_t3.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_mat <- config$host_pools[[1]]$infected
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
  test_mat <- config$host_pools[[2]]$infected
  expect_gte(data$host_pools[[2]]$infected[[1]][[1]], test_mat[[1]])
  expect_gte(data$host_pools[[2]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[2]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[2]]$infected[[1]][[4]], test_mat[[4]])
})

test_that("Test set 28: Using multiple hosts with uncertainty works as expected", {
  config_file <- system.file("extdata", "configs_pops", "ts28_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_mat <- config$host_pools[[1]]$infected
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]] + data$host_pools[[1]]$mortality[[1]][[1]],
             test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
  test_mat <- config$host_pools[[2]]$infected
  test_mat[[3]] <- 0
  expect_gte(data$host_pools[[2]]$infected[[1]][[1]] + data$host_pools[[2]]$infected[[1]][[1]],
             test_mat[[1]])
  expect_gte(data$host_pools[[2]]$infected[[1]][[2]] + data$host_pools[[2]]$infected[[1]][[2]],
             test_mat[[2]])
  expect_gte(data$host_pools[[2]]$infected[[1]][[3]] + data$host_pools[[2]]$infected[[1]][[3]],
             test_mat[[3]])
  expect_gte(data$host_pools[[2]]$infected[[1]][[4]] + data$host_pools[[2]]$infected[[1]][[4]],
             test_mat[[4]])

  config_file <- system.file("extdata", "configs_pops", "ts28_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_mat <- config$host_pools[[1]]$infected
  expect_gte(data$host_pools[[1]]$infected[[1]][[1]] + data$host_pools[[1]]$mortality[[1]][[1]],
             test_mat[[1]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[1]]$infected[[1]][[4]], test_mat[[4]])
  test_mat <- config$host_pools[[2]]$infected
  expect_gte(data$host_pools[[2]]$infected[[1]][[1]] + data$host_pools[[2]]$mortality[[1]][[1]]
             , test_mat[[1]])
  expect_gte(data$host_pools[[2]]$infected[[1]][[2]], test_mat[[2]])
  expect_gte(data$host_pools[[2]]$infected[[1]][[3]], test_mat[[3]])
  expect_gte(data$host_pools[[2]]$infected[[1]][[4]], test_mat[[4]])
})

test_that("Test set 29: county level infection works as expected", {
  config_file <- system.file("extdata", "configs_pops", "ts29_config_t1.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_infected <- terra::vect(system.file(config$starting_infected_files[[1]], package = "PoPS"))
  expect_equal(data$number_infected[[1]], sum(test_infected$infected_mean))

  config_file <- system.file("extdata", "configs_pops", "ts29_config_t2.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_infected <- terra::vect(system.file(config$starting_infected_files[[1]], package = "PoPS"))
  expect_gte(data$number_infected[[1]], sum(test_infected$infected_mean))

  config_file <- system.file("extdata", "configs_pops", "ts29_config_t3.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_infected <- terra::vect(system.file(config$starting_infected_files[[1]], package = "PoPS"))
  expect_gte(data$number_infected[[1]], 0)

  config_file <- system.file("extdata", "configs_pops", "ts29_config_t4.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_infected <- terra::vect(system.file(config$starting_infected_files[[1]], package = "PoPS"))
  expect_gte(data$number_infected[[1]], 0)

  config_file <- system.file("extdata", "configs_pops", "ts29_config_t5.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_infected <- terra::vect(system.file(config$starting_infected_files[[1]], package = "PoPS"))
  expect_gte(data$number_infected[[1]], 0)
  expect_gte(sum(data$host_pools[[1]]$total_exposed[[1]]), 0)

  config_file <- system.file("extdata", "configs_pops", "ts29_config_t6.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_infected <- terra::vect(system.file(config$starting_infected_files[[1]], package = "PoPS"))
  expect_equal(data$number_infected[[1]], sum(test_infected$infected_mean))
  expect_equal(sum(data$host_pools[[1]]$total_exposed[[1]]), sum(test_infected$infected_mean))

  config_file <- system.file("extdata", "configs_pops", "ts29_config_t7.yaml", package = "PoPS")
  config <- configuration(config_file = config_file, testing = TRUE)
  data <- pops(config)
  test_infected <- terra::vect(system.file(config$starting_infected_files[[1]], package = "PoPS"))
  expect_equal(data$number_infected[[1]], sum(test_infected$infected_mean))
  expect_equal(sum(data$host_pools[[1]]$total_exposed[[1]]), sum(test_infected$infected_mean))
})
