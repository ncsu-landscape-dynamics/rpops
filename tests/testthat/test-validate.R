context("test-validate")

test_that("Model stops if files don't exist or aren't the correct extension", {
  infected_file_list <- system.file("extdata", "simple2x2", "infected.tif", package = "PoPS")
  infected_years_file <-
    system.file("extdata", "simple20x20", "infected_single.tif", package = "PoPS")
  host_file_list <- system.file("extdata", "simple2x2", "total_plants.tif", package = "PoPS")
  parameter_means <- c(0, 21, 1, 500, 0, 0, 0, 0)
  parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
  pest_host_table <-
    system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
  competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")

  expect_error(validate(infected_years_file = infected_years_file,
                        infected_file_list = "",
                        host_file_list =  host_file_list,
                        total_populations_file =  host_file_list,
                        parameter_means = parameter_means,
                        parameter_cov_matrix = parameter_cov_matrix,
                        pest_host_table = pest_host_table,
                        competency_table = competency_table),
               PoPS:::file_exists_error, fixed = TRUE)

  expect_error(validate(infected_years_file = infected_years_file,
                        infected_file_list = infected_file_list,
                        host_file_list =  host_file_list,
                        total_populations_file =  host_file_list,
                        start_date = "2008-01-01",
                        end_date = "2009-12-31",
                        output_frequency = "year",
                        parameter_means = parameter_means,
                        parameter_cov_matrix = parameter_cov_matrix,
                        pest_host_table = pest_host_table,
                        competency_table = competency_table),
               PoPS:::infection_years_length_error(1, 2), fixed = TRUE)
})

test_that(
  "Validation has correctly formatted returns with multiple output comparisons", {
    skip_on_os("windows")
    infected_years_file <-
      system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
    parameter_means <- c(1.8, 16.4, 0.973, 7803, 0, 0, 0, 0)
    parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
    infected_file_list <-
      system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
    host_file_list <- system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
    total_populations_file <-
      system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
    temp <- FALSE
    temperature_coefficient_file <- ""
    precip <- FALSE
    precipitation_coefficient_file <- ""
    model_type <- "SEI"
    latency_period <- 14
    time_step <- "day"
    season_month_start <- 1
    season_month_end <- 12
    start_date <- "2003-01-01"
    end_date <- "2003-02-11"
    use_lethal_temperature <- FALSE
    temperature_file <- ""
    lethal_temperature <- -30
    lethal_temperature_month <- 1
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
    anthropogenic_distance_scale <- 0.0
    number_of_iterations <- 4
    number_of_cores <- 2
    start_exposed <- FALSE
    generate_stochasticity <- TRUE
    establishment_stochasticity <- TRUE
    movement_stochasticity <- TRUE
    dispersal_stochasticity <- FALSE
    establishment_probability <- 0.5
    dispersal_percentage <- 0.99
    quarantine_areas_file <- ""
    use_quarantine <- FALSE
    use_spreadrates <- FALSE
    use_overpopulation_movements <- FALSE
    overpopulation_percentage <- 0
    leaving_percentage <- 0
    leaving_scale_coefficient <- 1
    exposed_file_list <- ""
    write_outputs <- "all_simulations"
    output_folder_path <- tempdir()
    point_file <- ""
    use_survival_rates <- FALSE
    survival_rate_month <- 3
    survival_rate_day <- 15
    survival_rates_file <- ""
    network_filename <- ""
    network_movement <- "walk"
    use_distance <- FALSE
    use_configuration <- FALSE
    use_initial_condition_uncertainty <- FALSE
    use_host_uncertainty <- FALSE
    weather_type <- "deterministic"
    temperature_coefficient_sd_file <- ""
    precipitation_coefficient_sd_file <- ""
    dispersers_to_soils_percentage <- 0
    quarantine_directions <- ""
    multiple_random_seeds <- TRUE
    file_random_seeds <- NULL
    use_soils <- FALSE
    soil_starting_pest_file <- ""
    start_with_soil_populations <- FALSE

    outputs <- validate(
      infected_years_file,
      number_of_iterations,
      number_of_cores,
      parameter_means,
      parameter_cov_matrix,
      pest_host_table = pest_host_table,
      competency_table = competency_table,
      infected_file_list,
      host_file_list,
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
      use_survival_rates,
      survival_rate_month,
      survival_rate_day,
      survival_rates_file,
      use_lethal_temperature,
      temperature_file,
      lethal_temperature,
      lethal_temperature_month,
      mortality_frequency,
      mortality_frequency_n,
      management,
      treatment_dates,
      treatments_file,
      treatment_method,
      natural_kernel_type,
      anthropogenic_kernel_type,
      natural_dir,
      anthropogenic_dir,
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
      dispersal_stochasticity,
      establishment_probability,
      dispersal_percentage,
      quarantine_areas_file,
      use_quarantine,
      use_spreadrates,
      use_overpopulation_movements,
      overpopulation_percentage,
      leaving_percentage,
      leaving_scale_coefficient,
      exposed_file_list,
      write_outputs,
      output_folder_path,
      point_file,
      network_filename,
      network_movement,
      use_distance,
      use_configuration,
      use_initial_condition_uncertainty,
      use_host_uncertainty,
      weather_type,
      temperature_coefficient_sd_file,
      precipitation_coefficient_sd_file,
      dispersers_to_soils_percentage,
      quarantine_directions,
      multiple_random_seeds,
      file_random_seeds,
      use_soils,
      soil_starting_pest_file,
      start_with_soil_populations)

    expect_type(outputs, "list")
    expect_length(outputs, 12)
    data <- outputs[[1]]
    expect_length(data, 26)
    expect_vector(data$quantity_disagreement, size = number_of_iterations)
    expect_vector(data$allocation_disagreement, size = number_of_iterations)
    expect_vector(data$total_disagreement, size = number_of_iterations)
    expect_vector(data$configuration_disagreement, size = number_of_iterations)
    expect_vector(data$false_negatives, size = number_of_iterations)
    expect_vector(data$false_positives, size = number_of_iterations)
    expect_vector(data$true_positives, size = number_of_iterations)
    expect_vector(data$true_negatives, size = number_of_iterations)
    expect_vector(data$unknown_positives, size = number_of_iterations)
    expect_vector(data$unknown_negatives, size = number_of_iterations)
    expect_vector(data$odds_ratio, size = number_of_iterations)
    expect_vector(data$residual_error, size = number_of_iterations)
    expect_vector(data$true_infecteds, size = number_of_iterations)
    expect_vector(data$simulated_infecteds, size = number_of_iterations)
    expect_vector(data$infecteds_difference, size = number_of_iterations)
  }
)

test_that(
  "Validation has correctly formatted returns and runs with a single output comparison", {
    skip_on_os("windows")
    infected_years_file <-
      system.file("extdata", "simple20x20", "infected_single.tif", package = "PoPS")
    number_of_observations <- 68
    parameter_means <- c(1.8, 16.4, 0.973, 7803, 0, 0, 0, 0)
    parameter_cov_matrix <- matrix(ncol = 8, nrow = 8, 0)
    pest_host_table <-
      system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
    competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")
    checks <- c(500, 60000, 900, 1000)
    infected_file_list <-
      system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
    host_file_list <-
      system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
    total_populations_file <-
      system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
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
    end_date <- "2003-12-31"
    use_lethal_temperature <- FALSE
    temperature_file <- ""
    lethal_temperature <- -30
    lethal_temperature_month <- 1
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
    movements_file <- ""
    use_movements <- FALSE
    percent_natural_dispersal <- 1.0
    anthropogenic_distance_scale <- 0.0
    number_of_iterations <- 4
    number_of_cores <- 2
    start_exposed <- FALSE
    generate_stochasticity <- TRUE
    establishment_stochasticity <- TRUE
    movement_stochasticity <- TRUE
    dispersal_stochasticity <- FALSE
    establishment_probability <- 0.5
    dispersal_percentage <- 0.99
    quarantine_areas_file <- ""
    use_quarantine <- FALSE
    output_frequency_n <- 1
    use_spreadrates <- FALSE
    use_overpopulation_movements <- FALSE
    overpopulation_percentage <- 0
    leaving_percentage <- 0
    leaving_scale_coefficient <- 1
    exposed_file_list <- ""
    write_outputs <- "None"
    output_folder_path <- ""
    point_file <- system.file("extdata", "simple20x20", "points.gpkg", package = "PoPS")
    use_survival_rates <- FALSE
    survival_rate_month <- 3
    survival_rate_day <- 15
    survival_rates_file <- ""
    network_filename <- ""
    network_movement <- "walk"
    use_distance <- FALSE
    use_configuration <- FALSE
    use_initial_condition_uncertainty <- FALSE
    use_host_uncertainty <- FALSE
    weather_type <- "deterministic"
    temperature_coefficient_sd_file <- ""
    precipitation_coefficient_sd_file <- ""
    dispersers_to_soils_percentage <- 0
    quarantine_directions <- ""
    multiple_random_seeds <- FALSE
    file_random_seeds <- NULL
    use_soils <- FALSE
    soil_starting_pest_file <- ""
    start_with_soil_populations <- FALSE

    outputs <- validate(
      infected_years_file,
      number_of_iterations,
      number_of_cores,
      parameter_means,
      parameter_cov_matrix,
      pest_host_table = pest_host_table,
      competency_table = competency_table,
      infected_file_list,
      host_file_list,
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
      use_survival_rates,
      survival_rate_month,
      survival_rate_day,
      survival_rates_file,
      use_lethal_temperature,
      temperature_file,
      lethal_temperature,
      lethal_temperature_month,
      mortality_frequency,
      mortality_frequency_n,
      management,
      treatment_dates,
      treatments_file,
      treatment_method,
      natural_kernel_type,
      anthropogenic_kernel_type,
      natural_dir,
      anthropogenic_dir,
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
      dispersal_stochasticity,
      establishment_probability,
      dispersal_percentage,
      quarantine_areas_file,
      use_quarantine,
      use_spreadrates,
      use_overpopulation_movements,
      overpopulation_percentage,
      leaving_percentage,
      leaving_scale_coefficient,
      exposed_file_list,
      write_outputs,
      output_folder_path,
      point_file,
      network_filename,
      network_movement,
      use_distance,
      use_configuration,
      use_initial_condition_uncertainty,
      use_host_uncertainty,
      weather_type,
      temperature_coefficient_sd_file,
      precipitation_coefficient_sd_file,
      dispersers_to_soils_percentage,
      quarantine_directions,
      multiple_random_seeds,
      file_random_seeds,
      use_soils,
      soil_starting_pest_file,
      start_with_soil_populations)

    expect_type(outputs, "list")
    expect_length(outputs, 2)
    data <- outputs[[1]]
    expect_length(data, 35)
    expect_vector(data$quantity_disagreement, size = number_of_iterations)
    expect_vector(data$allocation_disagreement, size = number_of_iterations)
    expect_vector(data$total_disagreement, size = number_of_iterations)
    expect_vector(data$configuration_disagreement, size = number_of_iterations)
    expect_vector(data$false_negatives, size = number_of_iterations)
    expect_vector(data$false_positives, size = number_of_iterations)
    expect_vector(data$true_positives, size = number_of_iterations)
    expect_vector(data$true_negatives, size = number_of_iterations)
    expect_vector(data$unknown_positives, size = number_of_iterations)
    expect_vector(data$unknown_negatives, size = number_of_iterations)
    expect_vector(data$odds_ratio, size = number_of_iterations)
    expect_vector(data$residual_error, size = number_of_iterations)
    expect_vector(data$true_infecteds, size = number_of_iterations)
    expect_vector(data$simulated_infecteds, size = number_of_iterations)
    expect_vector(data$infecteds_difference, size = number_of_iterations)
  }
)

test_that(
  "Validation has correctly formatted returns and runs with a single output comparison with mask", {
    skip_on_os("windows")
    infected_years_file <-
      system.file("extdata", "simple20x20", "infected_single.tif", package = "PoPS")
    number_of_observations <- 68
    parameter_means <- c(1.8, 16.4, 0.973, 7803, 0, 0, 0, 0)
    parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
    pest_host_table <-
      system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
    competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")
    checks <- c(500, 60000, 900, 1000)
    infected_file_list <-
      system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
    host_file_list <-
      system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
    total_populations_file <-
      system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
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
    end_date <- "2003-12-31"
    use_lethal_temperature <- FALSE
    temperature_file <- ""
    lethal_temperature <- -30
    lethal_temperature_month <- 1
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
    mask <- system.file("extdata", "simple20x20", "mask.tif", package = "PoPS")
    output_frequency <- "year"
    movements_file <- ""
    use_movements <- FALSE
    percent_natural_dispersal <- 1.0
    anthropogenic_distance_scale <- 0.0
    number_of_iterations <- 4
    number_of_cores <- 2
    start_exposed <- FALSE
    generate_stochasticity <- TRUE
    establishment_stochasticity <- TRUE
    movement_stochasticity <- TRUE
    dispersal_stochasticity <- FALSE
    establishment_probability <- 0.5
    dispersal_percentage <- 0.99
    quarantine_areas_file <- ""
    use_quarantine <- FALSE
    output_frequency_n <- 1
    use_spreadrates <- FALSE
    use_overpopulation_movements <- FALSE
    overpopulation_percentage <- 0
    leaving_percentage <- 0
    leaving_scale_coefficient <- 1
    exposed_file_list <- ""
    write_outputs <- "None"
    output_folder_path <- ""
    point_file <- ""
    use_survival_rates <- FALSE
    survival_rate_month <- 3
    survival_rate_day <- 15
    survival_rates_file <- ""
    network_filename <- ""
    network_movement <- "walk"
    use_distance <- FALSE
    use_configuration <- FALSE
    use_initial_condition_uncertainty <- FALSE
    use_host_uncertainty <- FALSE
    weather_type <- "deterministic"
    temperature_coefficient_sd_file <- ""
    precipitation_coefficient_sd_file <- ""
    dispersers_to_soils_percentage <- 0
    quarantine_directions <- ""
    multiple_random_seeds <- FALSE
    file_random_seeds <- NULL
    use_soils <- FALSE
    soil_starting_pest_file <- ""
    start_with_soil_populations <- FALSE

    outputs <- validate(
      infected_years_file,
      number_of_iterations,
      number_of_cores,
      parameter_means,
      parameter_cov_matrix,
      pest_host_table = pest_host_table,
      competency_table = competency_table,
      infected_file_list,
      host_file_list,
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
      use_survival_rates,
      survival_rate_month,
      survival_rate_day,
      survival_rates_file,
      use_lethal_temperature,
      temperature_file,
      lethal_temperature,
      lethal_temperature_month,
      mortality_frequency,
      mortality_frequency_n,
      management,
      treatment_dates,
      treatments_file,
      treatment_method,
      natural_kernel_type,
      anthropogenic_kernel_type,
      natural_dir,
      anthropogenic_dir,
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
      dispersal_stochasticity,
      establishment_probability,
      dispersal_percentage,
      quarantine_areas_file,
      use_quarantine,
      use_spreadrates,
      use_overpopulation_movements,
      overpopulation_percentage,
      leaving_percentage,
      leaving_scale_coefficient,
      exposed_file_list,
      write_outputs,
      output_folder_path,
      point_file,
      network_filename,
      network_movement,
      use_distance,
      use_configuration,
      use_initial_condition_uncertainty,
      use_host_uncertainty,
      weather_type,
      temperature_coefficient_sd_file,
      precipitation_coefficient_sd_file,
      dispersers_to_soils_percentage,
      quarantine_directions,
      multiple_random_seeds,
      file_random_seeds,
      use_soils,
      soil_starting_pest_file,
      start_with_soil_populations)

    expect_type(outputs, "list")
    expect_length(outputs, 2)
    data <- outputs[[1]]
    expect_length(data, 26)
    expect_vector(data$quantity_disagreement, size = number_of_iterations)
    expect_vector(data$allocation_disagreement, size = number_of_iterations)
    expect_vector(data$total_disagreement, size = number_of_iterations)
    expect_vector(data$configuration_disagreement, size = number_of_iterations)
    expect_vector(data$false_negatives, size = number_of_iterations)
    expect_vector(data$false_positives, size = number_of_iterations)
    expect_vector(data$true_positives, size = number_of_iterations)
    expect_vector(data$true_negatives, size = number_of_iterations)
    expect_vector(data$unknown_positives, size = number_of_iterations)
    expect_vector(data$unknown_negatives, size = number_of_iterations)
    expect_vector(data$odds_ratio, size = number_of_iterations)
    expect_vector(data$residual_error, size = number_of_iterations)
    expect_vector(data$true_infecteds, size = number_of_iterations)
    expect_vector(data$simulated_infecteds, size = number_of_iterations)
    expect_vector(data$infecteds_difference, size = number_of_iterations)
  }
)


test_that(
  "Validation has correctly formatted returns/runs with host and initial condition uncertainty", {
    skip_on_os("windows")
    infected_years_file <-
      system.file("extdata", "simple20x20", "infected_single.tif", package = "PoPS")
    number_of_observations <- 68
    parameter_means <- c(1.8, 16.4, 0.973, 7803, 0, 0, 0, 0)
    parameter_cov_matrix <- matrix(0, nrow = 8, ncol = 8)
    pest_host_table <-
      system.file("extdata", "pest_host_table_singlehost_nomort.csv", package = "PoPS")
    competency_table <- system.file("extdata", "competency_table_singlehost.csv", package = "PoPS")
    checks <- c(500, 60000, 900, 1000)
    infected_file_list <-
      system.file("extdata", "simple20x20", "infected_wsd.tif", package = "PoPS")
    host_file_list <-
      system.file("extdata", "simple20x20", "host_w_sd2.tif", package = "PoPS")
    total_populations_file <-
      system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
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
    end_date <- "2003-12-31"
    use_lethal_temperature <- FALSE
    temperature_file <- ""
    lethal_temperature <- -30
    lethal_temperature_month <- 1
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
    mask <- system.file("extdata", "simple20x20", "mask.tif", package = "PoPS")
    output_frequency <- "year"
    movements_file <- ""
    use_movements <- FALSE
    percent_natural_dispersal <- 1.0
    anthropogenic_distance_scale <- 0.0
    number_of_iterations <- 4
    number_of_cores <- 2
    start_exposed <- FALSE
    generate_stochasticity <- TRUE
    establishment_stochasticity <- TRUE
    movement_stochasticity <- TRUE
    dispersal_stochasticity <- FALSE
    establishment_probability <- 0.5
    dispersal_percentage <- 0.99
    quarantine_areas_file <- ""
    use_quarantine <- FALSE
    output_frequency_n <- 1
    use_spreadrates <- FALSE
    use_overpopulation_movements <- FALSE
    overpopulation_percentage <- 0
    leaving_percentage <- 0
    leaving_scale_coefficient <- 1
    exposed_file_list <- ""
    write_outputs <- "None"
    output_folder_path <- ""
    point_file <- ""
    use_survival_rates <- FALSE
    survival_rate_month <- 3
    survival_rate_day <- 15
    survival_rates_file <- ""
    network_filename <- ""
    network_movement <- "walk"
    use_distance <- FALSE
    use_configuration <- FALSE
    use_initial_condition_uncertainty <- TRUE
    use_host_uncertainty <- TRUE
    weather_type <- "deterministic"
    temperature_coefficient_sd_file <- ""
    precipitation_coefficient_sd_file <- ""
    dispersers_to_soils_percentage <- 0
    quarantine_directions <- ""
    multiple_random_seeds <- FALSE
    file_random_seeds <- NULL
    use_soils <- FALSE
    soil_starting_pest_file <- ""
    start_with_soil_populations <- FALSE

    outputs <- validate(
      infected_years_file,
      number_of_iterations,
      number_of_cores,
      parameter_means,
      parameter_cov_matrix,
      pest_host_table = pest_host_table,
      competency_table = competency_table,
      infected_file_list,
      host_file_list,
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
      use_survival_rates,
      survival_rate_month,
      survival_rate_day,
      survival_rates_file,
      use_lethal_temperature,
      temperature_file,
      lethal_temperature,
      lethal_temperature_month,
      mortality_frequency,
      mortality_frequency_n,
      management,
      treatment_dates,
      treatments_file,
      treatment_method,
      natural_kernel_type,
      anthropogenic_kernel_type,
      natural_dir,
      anthropogenic_dir,
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
      dispersal_stochasticity,
      establishment_probability,
      dispersal_percentage,
      quarantine_areas_file,
      use_quarantine,
      use_spreadrates,
      use_overpopulation_movements,
      overpopulation_percentage,
      leaving_percentage,
      leaving_scale_coefficient,
      exposed_file_list,
      write_outputs,
      output_folder_path,
      point_file,
      network_filename,
      network_movement,
      use_distance,
      use_configuration,
      use_initial_condition_uncertainty,
      use_host_uncertainty,
      weather_type,
      temperature_coefficient_sd_file,
      precipitation_coefficient_sd_file,
      dispersers_to_soils_percentage,
      quarantine_directions,
      multiple_random_seeds,
      file_random_seeds,
      use_soils,
      soil_starting_pest_file,
      start_with_soil_populations)

    expect_type(outputs, "list")
    expect_length(outputs, 2)
    data <- outputs[[1]]
    expect_length(data, 26)
    expect_vector(data$quantity_disagreement, size = number_of_iterations)
    expect_vector(data$allocation_disagreement, size = number_of_iterations)
    expect_vector(data$total_disagreement, size = number_of_iterations)
    expect_vector(data$configuration_disagreement, size = number_of_iterations)
    expect_vector(data$false_negatives, size = number_of_iterations)
    expect_vector(data$false_positives, size = number_of_iterations)
    expect_vector(data$true_positives, size = number_of_iterations)
    expect_vector(data$true_negatives, size = number_of_iterations)
    expect_vector(data$unknown_positives, size = number_of_iterations)
    expect_vector(data$unknown_negatives, size = number_of_iterations)
    expect_vector(data$odds_ratio, size = number_of_iterations)
    expect_vector(data$residual_error, size = number_of_iterations)
    expect_vector(data$true_infecteds, size = number_of_iterations)
    expect_vector(data$simulated_infecteds, size = number_of_iterations)
    expect_vector(data$infecteds_difference, size = number_of_iterations)
  }
)
