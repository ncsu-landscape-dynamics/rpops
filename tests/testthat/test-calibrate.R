context("test-calibrate")

test_that("Model stops if files don't exist or aren't the correct extension", {
  infected_file <-
    system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  infected_years_file <-
    system.file("extdata", "simple20x20", "infected_years.tif",
                package = "PoPS"
    )
  host_file <-
    system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  prior_means <- c(0, 21, 1, 500, 0, 0)
  prior_cov_matrix <- matrix(0, nrow = 6, ncol = 6)
  number_of_observations = 1
  prior_number_of_observations = 0

  expect_error(calibrate(infected_years_file = infected_years_file,
                         number_of_observations = number_of_observations,
                         prior_number_of_observations = prior_number_of_observations,
                         prior_means = prior_means,
                         prior_cov_matrix = prior_cov_matrix,
                         infected_file = "",
                         host_file =  host_file,
                         total_populations_file =  host_file),
               "file does not exist")
})

test_that("ABC calibration has correctly formatted returns with multiple output
          comparisons with mask", {
            skip_on_ci()
            skip_on_covr()
            infected_years_file <-
              system.file("extdata", "simple20x20", "infected_years.tif",
                          package = "PoPS")
            number_of_observations <- 68
            prior_number_of_observations <- 0
            prior_means <- c(0, 0, 0, 0, 0, 0)
            prior_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
            params_to_estimate <- c(T, T, T, T, F, F)
            number_of_generations <- 2
            generation_size <- 2
            infected_file <-
              system.file("extdata", "simple20x20", "initial_infection.tif",
                          package = "PoPS")
            host_file <-
              system.file("extdata", "simple20x20", "host.tif",
                          package = "PoPS")
            total_populations_file <-
              system.file("extdata", "simple20x20", "all_plants.tif",
                          package = "PoPS")
            temp <- FALSE
            temperature_coefficient_file <- ""
            precip <- FALSE
            precipitation_coefficient_file <- ""
            model_type <- "SI"
            latency_period <- 14
            time_step <- "week"
            season_month_start <- 1
            season_month_end <- 12
            start_date <- "2003-01-01"
            end_date <- "2003-02-11"
            use_lethal_temperature <- FALSE
            temperature_file <- ""
            lethal_temperature <- -30
            lethal_temperature_month <- 1
            mortality_on <- FALSE
            mortality_rate <- 0
            mortality_time_lag <- 0
            mortality_frequency <- "Year"
            mortality_frequency_n <- 1
            management <- FALSE
            treatment_dates <- c("2003-01-24")
            treatments_file <- ""
            treatment_method <- "ratio"
            natural_kernel_type <- "exponential"
            anthropogenic_kernel_type <- "cauchy"
            natural_dir <- "NONE"
            natural_kappa <- 0
            anthropogenic_dir <- "NONE"
            anthropogenic_kappa <- 0
            pesticide_duration <- c(0)
            pesticide_efficacy <- 1.0
            mask <- system.file("extdata", "simple20x20", "mask.tif",
                                package = "PoPS")
            output_frequency <- "week"
            output_frequency_n <- 1
            movements_file <- ""
            use_movements <- FALSE
            percent_natural_dispersal <- 1.0
            anthropogenic_distance_scale <- 0.1
            number_of_iterations <- 10
            number_of_cores <- 2
            start_exposed <- FALSE
            generate_stochasticity <- TRUE
            establishment_stochasticity <- TRUE
            movement_stochasticity <- TRUE
            deterministic <- FALSE
            establishment_probability <- 0.5
            dispersal_percentage <- 0.99
            quarantine_areas_file <- ""
            use_quarantine <- FALSE
            use_spreadrates <- FALSE
            calibration_method <- "ABC"
            use_overpopulation_movements <- FALSE
            overpopulation_percentage <- 0
            leaving_percentage <- 0
            leaving_scale_coefficient <- 1
            exposed_file <- ""
            verbose <- TRUE
            write_outputs <- "None"
            output_folder_path <- ""
            use_distance <- FALSE
            use_rmse <- FALSE

            data <- calibrate(infected_years_file,
                              number_of_observations,
                              prior_number_of_observations,
                              prior_means,
                              prior_cov_matrix,
                              params_to_estimate,
                              number_of_generations,
                              generation_size,
                              infected_file,
                              host_file,
                              total_populations_file,
                              temp,
                              temperature_coefficient_file,
                              precip,
                              precipitation_coefficient_file,
                              model_type,
                              latency_period,
                              time_step,
                              season_month_start,
                              season_month_end,
                              start_date,
                              end_date,
                              use_lethal_temperature,
                              temperature_file,
                              lethal_temperature,
                              lethal_temperature_month,
                              mortality_on,
                              mortality_rate,
                              mortality_time_lag,
                              mortality_frequency,
                              mortality_frequency_n,
                              management,
                              treatment_dates,
                              treatments_file,
                              treatment_method,
                              natural_kernel_type,
                              anthropogenic_kernel_type,
                              natural_dir,
                              natural_kappa,
                              anthropogenic_dir,
                              anthropogenic_kappa,
                              pesticide_duration,
                              pesticide_efficacy,
                              mask,
                              output_frequency,
                              output_frequency_n,
                              movements_file,
                              use_movements,
                              start_exposed,
                              generate_stochasticity,
                              establishment_stochasticity,
                              movement_stochasticity,
                              deterministic,
                              establishment_probability,
                              dispersal_percentage,
                              quarantine_areas_file,
                              use_quarantine,
                              use_spreadrates,
                              use_overpopulation_movements,
                              overpopulation_percentage,
                              leaving_percentage,
                              leaving_scale_coefficient,
                              calibration_method,
                              number_of_iterations,
                              exposed_file,
                              verbose,
                              write_outputs,
                              output_folder_path,
                              use_distance,
                              use_rmse)

            expect_length(data$posterior_means, 6)
            expect_vector(data$posterior_means, ptype = double(), size = 6)
            expect_is(data$posterior_cov_matrix, class = "matrix")
            expect_type(data$posterior_cov_matrix, "double")
            expect_equal(nrow(data$posterior_cov_matrix), 6)
            expect_equal(ncol(data$posterior_cov_matrix), 6)
            expect_type(data$total_number_of_observations, "double")
            expect_equal(data$total_number_of_observations,
                         number_of_observations)
            expect_equal(nrow(data$raw_calibration_data),
                         number_of_generations * generation_size)
          })


test_that("ABC calibration has correctly formatted returns and runs with a
          single output comparison", {
            infected_years_file <-
              system.file("extdata", "simple20x20", "infected_single.tif",
                          package = "PoPS")
            number_of_observations <- 68
            prior_number_of_observations <- 0
            prior_means <- c(0, 0, 0, 0, 0, 0)
            prior_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
            params_to_estimate <- c(T, T, T, T, F, F)
            number_of_generations <- 2
            generation_size <- 2
            infected_file <-
              system.file("extdata", "simple20x20", "initial_infection.tif",
                          package = "PoPS")
            host_file <-
              system.file("extdata", "simple20x20", "host.tif",
                          package = "PoPS")
            total_populations_file <-
              system.file("extdata", "simple20x20", "all_plants.tif",
                          package = "PoPS")
            temp <- FALSE
            temperature_coefficient_file <- ""
            precip <- FALSE
            precipitation_coefficient_file <- ""
            model_type <- "SI"
            latency_period <- 0
            time_step <- "month"
            season_month_start <- 1
            season_month_end <- 12
            start_date <- "2003-01-01"
            end_date <- "2003-02-11"
            use_lethal_temperature <- FALSE
            temperature_file <- ""
            lethal_temperature <- -30
            lethal_temperature_month <- 1
            mortality_on <- FALSE
            mortality_rate <- 0
            mortality_time_lag <- 0
            mortality_frequency <- "Year"
            mortality_frequency_n <- 1
            management <- FALSE
            treatment_dates <- c("2003-01-24")
            treatments_file <- ""
            treatment_method <- "ratio"
            natural_kernel_type <- "exponential"
            anthropogenic_kernel_type <- "cauchy"
            natural_dir <- "NONE"
            natural_kappa <- 0
            anthropogenic_dir <- "NONE"
            anthropogenic_kappa <- 0
            pesticide_duration <- c(0)
            pesticide_efficacy <- 1.0
            mask <- NULL
            output_frequency <- "year"
            output_frequency_n <- 1
            movements_file <- ""
            use_movements <- FALSE
            percent_natural_dispersal <- 1.0
            anthropogenic_distance_scale <- 0.1
            number_of_iterations <- 10
            number_of_cores <- 2
            start_exposed <- FALSE
            generate_stochasticity <- TRUE
            establishment_stochasticity <- TRUE
            movement_stochasticity <- TRUE
            deterministic <- FALSE
            establishment_probability <- 0.5
            dispersal_percentage <- 0.99
            quarantine_areas_file <- ""
            use_quarantine <- FALSE
            output_frequency_n <- 1
            use_spreadrates <- FALSE
            calibration_method <- "ABC"
            number_of_iterations <- 10
            use_overpopulation_movements <- FALSE
            overpopulation_percentage <- 0
            leaving_percentage <- 0
            leaving_scale_coefficient <- 1
            exposed_file <- ""
            verbose <- TRUE
            write_outputs <- "None"
            output_folder_path <- ""
            use_distance <- FALSE
            use_rmse <- FALSE

            data <- calibrate(infected_years_file,
                              number_of_observations,
                              prior_number_of_observations,
                              prior_means,
                              prior_cov_matrix,
                              params_to_estimate,
                              number_of_generations,
                              generation_size,
                              infected_file,
                              host_file,
                              total_populations_file,
                              temp,
                              temperature_coefficient_file,
                              precip,
                              precipitation_coefficient_file,
                              model_type,
                              latency_period,
                              time_step,
                              season_month_start,
                              season_month_end,
                              start_date,
                              end_date,
                              use_lethal_temperature,
                              temperature_file,
                              lethal_temperature,
                              lethal_temperature_month,
                              mortality_on,
                              mortality_rate,
                              mortality_time_lag,
                              mortality_frequency,
                              mortality_frequency_n,
                              management,
                              treatment_dates,
                              treatments_file,
                              treatment_method,
                              natural_kernel_type,
                              anthropogenic_kernel_type,
                              natural_dir,
                              natural_kappa,
                              anthropogenic_dir,
                              anthropogenic_kappa,
                              pesticide_duration,
                              pesticide_efficacy,
                              mask,
                              output_frequency,
                              output_frequency_n,
                              movements_file,
                              use_movements,
                              start_exposed,
                              generate_stochasticity,
                              establishment_stochasticity,
                              movement_stochasticity,
                              deterministic,
                              establishment_probability,
                              dispersal_percentage,
                              quarantine_areas_file,
                              use_quarantine,
                              use_spreadrates,
                              use_overpopulation_movements,
                              overpopulation_percentage,
                              leaving_percentage,
                              leaving_scale_coefficient,
                              calibration_method,
                              number_of_iterations,
                              exposed_file,
                              verbose,
                              write_outputs,
                              output_folder_path,
                              use_distance,
                              use_rmse)

            expect_length(data$posterior_means, 6)
            expect_vector(data$posterior_means, ptype = double(), size = 6)
            expect_is(data$posterior_cov_matrix, class = "matrix")
            expect_type(data$posterior_cov_matrix, "double")
            expect_equal(nrow(data$posterior_cov_matrix), 6)
            expect_equal(ncol(data$posterior_cov_matrix), 6)
            expect_type(data$total_number_of_observations, "double")
            expect_equal(data$total_number_of_observations,
                         number_of_observations)
            expect_equal(nrow(data$raw_calibration_data),
                         number_of_generations * generation_size)

          })

test_that("MCMC calibration has correctly formatted returns with multiple output
          comparisons", {
            infected_years_file <-
              system.file("extdata", "simple20x20", "infected_years.tif",
                          package = "PoPS")
            number_of_observations <- 68
            prior_number_of_observations <- 0
            prior_means <- c(0, 0, 0, 0, 0, 0)
            prior_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
            params_to_estimate <- c(T, T, T, T, F, F)
            number_of_generations <- 2
            generation_size <- 2
            infected_file <-
              system.file("extdata", "simple20x20", "initial_infection.tif",
                          package = "PoPS")
            host_file <-
              system.file("extdata", "simple20x20", "host.tif",
                          package = "PoPS")
            total_populations_file <-
              system.file("extdata", "simple20x20", "all_plants.tif",
                          package = "PoPS")
            temp <- FALSE
            temperature_coefficient_file <- ""
            precip <- FALSE
            precipitation_coefficient_file <- ""
            model_type <- "SI"
            latency_period <- 14
            time_step <- "week"
            season_month_start <- 1
            season_month_end <- 12
            start_date <- "2003-01-01"
            end_date <- "2003-02-11"
            use_lethal_temperature <- FALSE
            temperature_file <- ""
            lethal_temperature <- -30
            lethal_temperature_month <- 1
            mortality_on <- FALSE
            mortality_rate <- 0
            mortality_time_lag <- 0
            mortality_frequency <- "Year"
            mortality_frequency_n <- 1
            management <- FALSE
            treatment_dates <- c("2003-01-24")
            treatments_file <- ""
            treatment_method <- "ratio"
            natural_kernel_type <- "exponential"
            anthropogenic_kernel_type <- "cauchy"
            natural_dir <- "NONE"
            natural_kappa <- 0
            anthropogenic_dir <- "NONE"
            anthropogenic_kappa <- 0
            pesticide_duration <- c(0)
            pesticide_efficacy <- 1.0
            mask <- NULL
            output_frequency <- "week"
            output_frequency_n <- 1
            movements_file <- ""
            use_movements <- FALSE
            percent_natural_dispersal <- 1.0
            anthropogenic_distance_scale <- 0.1
            number_of_iterations <- 10
            number_of_cores <- 2
            start_exposed <- FALSE
            generate_stochasticity <- TRUE
            establishment_stochasticity <- TRUE
            movement_stochasticity <- TRUE
            deterministic <- FALSE
            establishment_probability <- 0.5
            dispersal_percentage <- 0.99
            quarantine_areas_file <- ""
            use_quarantine <- FALSE
            use_spreadrates <- FALSE
            calibration_method <- "MCMC"
            use_overpopulation_movements <- FALSE
            overpopulation_percentage <- 0
            leaving_percentage <- 0
            leaving_scale_coefficient <- 1
            exposed_file <- ""
            verbose <- TRUE
            write_outputs <- "None"
            output_folder_path <- ""
            use_distance <- FALSE
            use_rmse <- FALSE

            data <- calibrate(infected_years_file,
                              number_of_observations,
                              prior_number_of_observations,
                              prior_means,
                              prior_cov_matrix,
                              params_to_estimate,
                              number_of_generations,
                              generation_size,
                              infected_file,
                              host_file,
                              total_populations_file,
                              temp,
                              temperature_coefficient_file,
                              precip,
                              precipitation_coefficient_file,
                              model_type,
                              latency_period,
                              time_step,
                              season_month_start,
                              season_month_end,
                              start_date,
                              end_date,
                              use_lethal_temperature,
                              temperature_file,
                              lethal_temperature,
                              lethal_temperature_month,
                              mortality_on,
                              mortality_rate,
                              mortality_time_lag,
                              mortality_frequency,
                              mortality_frequency_n,
                              management,
                              treatment_dates,
                              treatments_file,
                              treatment_method,
                              natural_kernel_type,
                              anthropogenic_kernel_type,
                              natural_dir,
                              natural_kappa,
                              anthropogenic_dir,
                              anthropogenic_kappa,
                              pesticide_duration,
                              pesticide_efficacy,
                              mask,
                              output_frequency,
                              output_frequency_n,
                              movements_file,
                              use_movements,
                              start_exposed,
                              generate_stochasticity,
                              establishment_stochasticity,
                              movement_stochasticity,
                              deterministic,
                              establishment_probability,
                              dispersal_percentage,
                              quarantine_areas_file,
                              use_quarantine,
                              use_spreadrates,
                              use_overpopulation_movements,
                              overpopulation_percentage,
                              leaving_percentage,
                              leaving_scale_coefficient,
                              calibration_method,
                              number_of_iterations,
                              exposed_file,
                              verbose,
                              write_outputs,
                              output_folder_path,
                              use_distance,
                              use_rmse)

            expect_length(data$posterior_means, 6)
            expect_vector(data$posterior_means, ptype = double(), size = 6)
            expect_is(data$posterior_cov_matrix, class = "matrix")
            expect_type(data$posterior_cov_matrix, "double")
            expect_equal(nrow(data$posterior_cov_matrix), 6)
            expect_equal(ncol(data$posterior_cov_matrix), 6)
            expect_type(data$total_number_of_observations, "double")
            expect_equal(data$total_number_of_observations,
                         number_of_observations)
            expect_equal(nrow(data$raw_calibration_data),
                         number_of_iterations)

          })

test_that("MCMC calibration has correctly formatted returns with multiple output
          comparisons with mask", {
            skip_on_ci()
            skip_on_covr()
            infected_years_file <-
              system.file("extdata", "simple20x20", "infected_years.tif",
                          package = "PoPS")
            number_of_observations <- 68
            prior_number_of_observations <- 0
            prior_means <- c(0, 0, 0, 0, 0, 0)
            prior_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
            params_to_estimate <- c(T, T, T, T, F, F)
            number_of_generations <- 2
            generation_size <- 2
            infected_file <-
              system.file("extdata", "simple20x20", "initial_infection.tif",
                          package = "PoPS")
            host_file <-
              system.file("extdata", "simple20x20", "host.tif",
                          package = "PoPS")
            total_populations_file <-
              system.file("extdata", "simple20x20", "all_plants.tif",
                          package = "PoPS")
            temp <- FALSE
            temperature_coefficient_file <- ""
            precip <- FALSE
            precipitation_coefficient_file <- ""
            model_type <- "SI"
            latency_period <- 14
            time_step <- "week"
            season_month_start <- 1
            season_month_end <- 12
            start_date <- "2003-01-01"
            end_date <- "2003-02-11"
            use_lethal_temperature <- FALSE
            temperature_file <- ""
            lethal_temperature <- -30
            lethal_temperature_month <- 1
            mortality_on <- FALSE
            mortality_rate <- 0
            mortality_time_lag <- 0
            mortality_frequency <- "Year"
            mortality_frequency_n <- 1
            management <- FALSE
            treatment_dates <- c("2003-01-24")
            treatments_file <- ""
            treatment_method <- "ratio"
            natural_kernel_type <- "exponential"
            anthropogenic_kernel_type <- "cauchy"
            natural_dir <- "NONE"
            natural_kappa <- 0
            anthropogenic_dir <- "NONE"
            anthropogenic_kappa <- 0
            pesticide_duration <- c(0)
            pesticide_efficacy <- 1.0
            mask <- system.file("extdata", "simple20x20", "mask.tif",
                                package = "PoPS")
            output_frequency <- "week"
            output_frequency_n <- 1
            movements_file <- ""
            use_movements <- FALSE
            percent_natural_dispersal <- 1.0
            anthropogenic_distance_scale <- 0.1
            number_of_iterations <- 10
            number_of_cores <- 2
            start_exposed <- FALSE
            generate_stochasticity <- TRUE
            establishment_stochasticity <- TRUE
            movement_stochasticity <- TRUE
            deterministic <- FALSE
            establishment_probability <- 0.5
            dispersal_percentage <- 0.99
            quarantine_areas_file <- ""
            use_quarantine <- FALSE
            use_spreadrates <- FALSE
            calibration_method <- "MCMC"
            use_overpopulation_movements <- FALSE
            overpopulation_percentage <- 0
            leaving_percentage <- 0
            leaving_scale_coefficient <- 1
            exposed_file <- ""
            verbose <- TRUE
            write_outputs <- "None"
            output_folder_path <- ""
            use_distance <- FALSE
            use_rmse <- FALSE

            data <- calibrate(infected_years_file,
                              number_of_observations,
                              prior_number_of_observations,
                              prior_means,
                              prior_cov_matrix,
                              params_to_estimate,
                              number_of_generations,
                              generation_size,
                              infected_file,
                              host_file,
                              total_populations_file,
                              temp,
                              temperature_coefficient_file,
                              precip,
                              precipitation_coefficient_file,
                              model_type,
                              latency_period,
                              time_step,
                              season_month_start,
                              season_month_end,
                              start_date,
                              end_date,
                              use_lethal_temperature,
                              temperature_file,
                              lethal_temperature,
                              lethal_temperature_month,
                              mortality_on,
                              mortality_rate,
                              mortality_time_lag,
                              mortality_frequency,
                              mortality_frequency_n,
                              management,
                              treatment_dates,
                              treatments_file,
                              treatment_method,
                              natural_kernel_type,
                              anthropogenic_kernel_type,
                              natural_dir,
                              natural_kappa,
                              anthropogenic_dir,
                              anthropogenic_kappa,
                              pesticide_duration,
                              pesticide_efficacy,
                              mask,
                              output_frequency,
                              output_frequency_n,
                              movements_file,
                              use_movements,
                              start_exposed,
                              generate_stochasticity,
                              establishment_stochasticity,
                              movement_stochasticity,
                              deterministic,
                              establishment_probability,
                              dispersal_percentage,
                              quarantine_areas_file,
                              use_quarantine,
                              use_spreadrates,
                              use_overpopulation_movements,
                              overpopulation_percentage,
                              leaving_percentage,
                              leaving_scale_coefficient,
                              calibration_method,
                              number_of_iterations,
                              exposed_file,
                              verbose,
                              write_outputs,
                              output_folder_path,
                              use_distance,
                              use_rmse)

            expect_length(data$posterior_means, 6)
            expect_vector(data$posterior_means, ptype = double(), size = 6)
            expect_is(data$posterior_cov_matrix, class = "matrix")
            expect_type(data$posterior_cov_matrix, "double")
            expect_equal(nrow(data$posterior_cov_matrix), 6)
            expect_equal(ncol(data$posterior_cov_matrix), 6)
            expect_type(data$total_number_of_observations, "double")
            expect_equal(data$total_number_of_observations,
                         number_of_observations)
            expect_equal(nrow(data$raw_calibration_data),
                         number_of_iterations)

          })

