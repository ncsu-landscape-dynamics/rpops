context("test-validate")

test_that(
  "Validation has correctly formatted returns with multiple output
  comparisons", {
    skip_on_os("windows")
    infected_years_file <-
      system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
    parameter_means <- c(1.8, 16.4, 0.973, 7803, 0, 0)
    parameter_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
    infected_file <-
      system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
    host_file <-
      system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
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
    anthropogenic_distance_scale <- 0.0
    number_of_iterations <- 4
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
    use_overpopulation_movements <- FALSE
    overpopulation_percentage <- 0
    leaving_percentage <- 0
    leaving_scale_coefficient <- 1
    exposed_file <- ""
    write_outputs <- "None"
    output_folder_path <- ""
    point_file <- ""

    outputs <- validate(
      infected_years_file,
      number_of_iterations,
      number_of_cores,
      parameter_means,
      parameter_cov_matrix,
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
      exposed_file,
      write_outputs,
      output_folder_path,
      point_file)

    expect_type(outputs, "list")
    expect_length(outputs, 12)
    data <- outputs[[1]]
    expect_length(data, 22)
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
    expect_vector(data$true_infected, size = number_of_iterations)
    expect_vector(data$simulated_infected, size = number_of_iterations)
    expect_vector(data$infected_difference, size = number_of_iterations)
  }
)

test_that(
  "Validation has correctly formatted returns and runs with a
  single output comparison", {
    skip_on_os("windows")
    infected_years_file <-
      system.file("extdata", "simple20x20", "infected_single.tif", package = "PoPS")
    number_of_observations <- 68
    parameter_means <- c(1.8, 16.4, 0.973, 7803, 0, 0)
    parameter_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
    checks <- c(500, 60000, 900, 1000)
    infected_file <-
      system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
    host_file <-
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
    deterministic <- FALSE
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
    exposed_file <- ""
    write_outputs <- "None"
    output_folder_path <- ""

    outputs <- validate(
      infected_years_file,
      number_of_iterations,
      number_of_cores,
      parameter_means,
      parameter_cov_matrix,
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
      exposed_file,
      write_outputs,
      output_folder_path
    )

    expect_type(outputs, "list")
    expect_length(outputs, 2)
    data <- outputs[[1]]
    expect_length(data, 22)
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
    expect_vector(data$true_infected, size = number_of_iterations)
    expect_vector(data$simulated_infected, size = number_of_iterations)
    expect_vector(data$infected_difference, size = number_of_iterations)
  }
)

test_that(
  "Validation has correctly formatted returns and runs with a
  single output comparison with mask", {
    skip_on_os("windows")
    infected_years_file <-
      system.file("extdata", "simple20x20", "infected_single.tif", package = "PoPS")
    number_of_observations <- 68
    parameter_means <- c(1.8, 16.4, 0.973, 7803, 0, 0)
    parameter_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
    checks <- c(500, 60000, 900, 1000)
    infected_file <-
      system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
    host_file <-
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
    deterministic <- FALSE
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
    exposed_file <- ""
    write_outputs <- "None"
    output_folder_path <- ""

    outputs <- validate(
      infected_years_file,
      number_of_iterations,
      number_of_cores,
      parameter_means,
      parameter_cov_matrix,
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
      exposed_file,
      write_outputs,
      output_folder_path
    )

    expect_type(outputs, "list")
    expect_length(outputs, 2)
    data <- outputs[[1]]
    expect_length(data, 22)
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
    expect_vector(data$true_infected, size = number_of_iterations)
    expect_vector(data$simulated_infected, size = number_of_iterations)
    expect_vector(data$infected_difference, size = number_of_iterations)
  }
)
