# These functions are designed to improve data format checks and reduce copy
# and pasting of code across functions
Sys.setenv("R_TESTS" = "")

initial_raster_checks <- function(x, use_s3 = FALSE, bucket = "") {
  checks_passed <- TRUE

  if (use_s3) {
    if (!aws.s3::head_object(x, bucket)) {
      checks_passed <- FALSE
      failed_check <- file_exists_error
    }
  } else {
    if (!all(file.exists(x))) {
      checks_passed <- FALSE
      failed_check <- file_exists_error
    }
  }

  if (checks_passed && !all((tools::file_ext(x) %in% raster_list))) {
    checks_passed <- FALSE
    failed_check <- raster_type_error
  }

  if (checks_passed) {
    if (use_s3) {
      aws.s3::save_object(object = x, bucket = bucket, file = x, check_region = FALSE)
      r <- terra::rast(x)
    } else {
      r <- terra::rast(x)
    }
  }

  if (checks_passed) {
    outs <- list(checks_passed, r)
    names(outs) <- c("checks_passed", "raster")
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- failed_check_list
    return(outs)
  }
}

# adds checks to test for raster extent, resolution, and crs x2 being the
# raster already through initial checks for comparison
secondary_raster_checks <- function(x, x2, use_s3 = FALSE, bucket = "") {
  checks_passed <- TRUE

  if (use_s3) {
    if (!aws.s3::head_object(x, bucket)) {
      checks_passed <- FALSE
      failed_check <- file_exists_error
    }
  } else {
    if (!all(file.exists(x))) {
      checks_passed <- FALSE
      failed_check <- file_exists_error
    }
  }

  if (checks_passed && !all((tools::file_ext(x) %in% raster_list))) {
    checks_passed <- FALSE
    failed_check <- raster_type_error
  }

  if (checks_passed) {
    if (use_s3) {
      aws.s3::save_object(object = x, bucket = bucket, file = x, check_region = FALSE)
      r <- terra::rast(x)
    } else {
      r <- terra::rast(x)
    }
  }

  if (checks_passed && !(terra::ext(x2) == terra::ext(r))) {
    checks_passed <- FALSE
    failed_check <- extent_error
  }

  if (checks_passed && !(terra::xres(x2) == terra::xres(r) &&
    terra::yres(x2) == terra::yres(r))) {
    checks_passed <- FALSE
    failed_check <- resolution_error
  }

  if (checks_passed) {
    crs1 <- terra::crs(r, describe = TRUE)
    crs2 <- terra::crs(x2, describe = TRUE)
    if (is.na(crs1$code)) {
      crs1$code <- "1"
      }
    if (is.na(crs2$code)) {
      crs2$code <- "1"
      }
    if (!(crs1$code == crs2$code)) {
      checks_passed <- FALSE
      failed_check <- crs_error
    }
  }

  if (checks_passed) {
    outs <- list(checks_passed, r)
    names(outs) <- c("checks_passed", "raster")
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- failed_check_list
    return(outs)
  }
}

treatment_checks <- function(treatment_stack,
                             treatments_file,
                             pesticide_duration,
                             treatment_dates,
                             pesticide_efficacy) {
  checks_passed <- TRUE

  if (checks_passed && length(treatments_file) != length(treatment_dates)) {
    checks_passed <- FALSE
    failed_check <- treatment_length_error
  }

  if (checks_passed && length(pesticide_duration) != length(treatment_dates)) {
    checks_passed <- FALSE
    failed_check <- pesticide_length_error
  }

  if (checks_passed) {
    if (pesticide_duration[1] > 0) {
      treatment_maps <-
        list(terra::as.matrix(treatment_stack[[1]] * pesticide_efficacy, wide = TRUE))
    } else {
      treatment_maps <- list(terra::as.matrix(treatment_stack[[1]], wide = TRUE))
    }

    if (terra::nlyr(treatment_stack) >= 2) {
      for (i in 2:terra::nlyr(treatment_stack)) {
        if (pesticide_duration[i] > 0) {
          treatment_maps[[i]] <-
            list(terra::as.matrix(treatment_stack[[i]] * pesticide_efficacy, wide = TRUE))
        } else {
          treatment_maps[[i]] <- list(terra::as.matrix(treatment_stack[[i]], wide = TRUE))
        }
      }
    }
  }

  if (checks_passed) {
    outs <- list(checks_passed, treatment_maps)
    names(outs) <- c("checks_passed", "treatment_maps")
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- failed_check_list
    return(outs)
  }
}

treatment_metric_checks <- function(treatment_method) {
  checks_passed <- TRUE

  if (!treatment_method %in% c("ratio", "all infected")) {
    checks_passed <- FALSE
    failed_check <- treatment_method_error
  }

  if (checks_passed) {
    outs <- list(checks_passed)
    names(outs) <- c("checks_passed")
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- failed_check_list
    return(outs)
  }
}

time_checks <- function(end_date, start_date, time_step, output_frequency, output_frequency_n) {
  checks_passed <- TRUE

  if (checks_passed && !(time_step %in% list("week", "month", "day"))) {
    checks_passed <- FALSE
    failed_check <- time_step_error
  }

  if (checks_passed && (!is(end_date, "character") ||
    !is(start_date, "character") ||
    !is(as.Date(end_date, format = "%Y-%m-%d"), "Date") ||
    !is(as.Date(start_date, format = "%Y-%m-%d"), "Date") ||
    is.na(as.Date(end_date, format = "%Y-%m-%d")) ||
    is.na(as.Date(start_date, format = "%Y-%m-%d")))) {
    checks_passed <- FALSE
    failed_check <- date_format_error
  }

  if (checks_passed && !(output_frequency %in% output_frequency_list)) {
    checks_passed <- FALSE
    failed_check <- output_type_error
  }

  if (checks_passed && output_frequency == "day") {
    if (time_step == "week" || time_step == "month") {
      checks_passed <- FALSE
      failed_check <- output_frequency_error
    }
  }

  if (checks_passed && output_frequency == "week") {
    if (time_step == "month") {
      checks_passed <- FALSE
      failed_check <- output_frequency_error
    }
  }

  if (checks_passed) {
    duration <- lubridate::interval(start_date, end_date)

    if (time_step == "week") {
      number_of_time_steps <- ceiling(lubridate::time_length(duration, "week"))
    } else if (time_step == "month") {
      number_of_time_steps <- ceiling(lubridate::time_length(duration, "month"))
    } else if (time_step == "day") {
      number_of_time_steps <- ceiling(lubridate::time_length(duration, "day"))
    }

    number_of_years <- ceiling(lubridate::time_length(duration, "year"))

    if (output_frequency == "time_step") {
      output_frequency <- time_step
    }

    if (output_frequency == "week") {
      number_of_outputs <- ceiling(lubridate::time_length(duration, "week"))
    } else if (output_frequency == "month") {
      number_of_outputs <- ceiling(lubridate::time_length(duration, "month"))
    } else if (output_frequency == "day") {
      number_of_outputs <- ceiling(lubridate::time_length(duration, "day"))
    } else if (output_frequency == "year") {
      number_of_outputs <- ceiling(lubridate::time_length(duration, "year"))
    } else if (output_frequency == "time_step") {
      number_of_outputs <- number_of_time_steps
    } else if (output_frequency == "every_n_steps") {
      number_of_outputs <- number_of_time_steps / output_frequency_n
    }

    if (output_frequency == "year" && time_step == "day" &&
      number_of_time_steps < 365) {
      output_frequency <- "final_step"
    } else if (output_frequency == "year" && time_step == "week" &&
      number_of_time_steps < 52) {
      output_frequency <- "final_step"
    } else if (output_frequency == "year" && time_step == "month" &&
      number_of_time_steps < 12) {
      output_frequency <- "final_step"
    } else {
      output_frequency <- output_frequency
    }
  }

  if (checks_passed) {
    outs <- list(
      checks_passed, number_of_time_steps, number_of_years,
      number_of_outputs, output_frequency
    )
    names(outs) <- c("checks_passed", "number_of_time_steps", "number_of_years",
                     "number_of_outputs", "output_frequency")
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- failed_check_list
    return(outs)
  }
}

bayesian_mnn_checks <- function(prior_means,
                                prior_cov_matrix,
                                calibrated_means,
                                calibrated_cov_matrix,
                                prior_weight, weight) {
  checks_passed <- TRUE
  if (length(prior_means) == length(calibrated_means) && prior_weight > 0) {
    posterior_means <- prior_means * prior_weight + calibrated_means * weight
  } else if (prior_weight == 0) {
    posterior_means <- calibrated_means
  } else {
    checks_passed <- FALSE
    failed_check <- prior_means_error
  }

  if (nrow(prior_cov_matrix) == nrow(calibrated_cov_matrix) &&
    ncol(prior_cov_matrix) == ncol(calibrated_cov_matrix) &&
    prior_weight > 0) {
    posterior_cov_matrix <- prior_cov_matrix * prior_weight +
      calibrated_cov_matrix * weight
  } else if (prior_weight == 0) {
    posterior_cov_matrix <- calibrated_cov_matrix
  } else {
    checks_passed <- FALSE
    failed_check <- prior_cov_matrix_error
  }

  if (checks_passed) {
    outs <- list(checks_passed, posterior_means, posterior_cov_matrix)
    names(outs) <- c("checks_passed", "posterior_means", "posterior_cov_matrix")
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- failed_check_list
    return(outs)
  }
}

multispecies_checks <- function(species,
                                infected_files,
                                parameter_means,
                                parameter_cov_matrix,
                                natural_kernel_type,
                                anthropogenic_kernel_type,
                                natural_dir,
                                anthropogenic_dir,
                                model_type,
                                host_file,
                                total_populations_file,
                                temp,
                                temperature_coefficient_file,
                                precip,
                                precipitation_coefficient_file,
                                latency_period,
                                time_step,
                                season_month_start,
                                season_month_end,
                                use_lethal_temperature,
                                temperature_file,
                                lethal_temperature,
                                lethal_temperature_month,
                                mortality_on,
                                mortality_rate,
                                mortality_time_lag,
                                movements_file,
                                use_movements,
                                start_exposed,
                                quarantine_areas_file,
                                use_quarantine,
                                use_spreadrates) {
  checks_passed <- TRUE

  if (checks_passed && length(species) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- species_length_error
  }

  if (checks_passed && length(parameter_means) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- parameter_length_error
  }

  if (checks_passed && length(parameter_cov_matrix) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- cov_matrix_length_error
  }

  if (checks_passed && length(model_type) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- model_type_length_error
  }

  if (checks_passed && length(natural_kernel_type) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- natural_kernel_length_error
  }

  if (checks_passed && length(anthropogenic_kernel_type) !=
    length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- anthropogenic_kernel_length_error
  }

  if (checks_passed && length(natural_dir) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- natural_dir_length_error
  }

  if (checks_passed && length(anthropogenic_dir) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- anthropogenic_dir_length_error
  }

  if (checks_passed && length(host_file) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- host_file_length_error
  }

  if (checks_passed && length(total_populations_file) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- total_population_length_error
  }

  if (checks_passed && length(temp) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- temp_length_error
  }

  if (checks_passed && length(temperature_coefficient_file) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- temperature_coefficient_length_error
  }

  if (checks_passed && length(precip) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- precip_length_error
  }

  if (checks_passed && length(precipitation_coefficient_file) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- precipitation_coefficient_length_error
  }

  if (checks_passed && length(latency_period) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- latency_period_length_error
  }

  if (checks_passed && length(time_step) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- time_step_length_error
  }

  if (checks_passed && length(season_month_start) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- season_month_start_length_error
  }

  if (checks_passed && length(season_month_end) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- season_month_end_length_error
  }

  if (checks_passed && length(use_lethal_temperature) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- use_lethal_length_error
  }

  if (checks_passed && length(temperature_file) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- temperature_file_length_error
  }

  if (checks_passed && length(lethal_temperature) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- lethal_temperature_length_error
  }

  if (checks_passed && length(lethal_temperature_month) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- lethal_temperature_month_length_error
  }

  if (checks_passed && length(mortality_on) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- mortality_on_length_error
  }

  if (checks_passed && length(mortality_rate) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- mortality_rate_length_error
  }

  if (checks_passed && length(mortality_time_lag) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- mortality_time_lag_length_error
  }

  if (checks_passed && length(movements_file) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- movements_file_length_error
  }

  if (checks_passed && length(use_movements) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- use_movements_length_error
  }

  if (checks_passed && length(start_exposed) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- start_exposed_length_error
  }

  if (checks_passed && length(quarantine_areas_file) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- quarantine_areas_length_error
  }

  if (checks_passed && length(use_quarantine) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- use_quarantine_length_error
  }

  if (checks_passed && length(use_spreadrates) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- use_spreadrates_length_error
  }

  if (checks_passed) {
    outs <- list(checks_passed)
    names(outs) <- c("checks_passed")
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- failed_check_list
    return(outs)
  }
}

movement_checks <- function(x, rast, start_date, end_date) {
  checks_passed <- TRUE

  if (!all(file.exists(x))) {
    checks_passed <- FALSE
    failed_check <- file_exists_error
  }

  if (checks_passed && !all((tools::file_ext(x) %in% csv_list))) {
    checks_passed <- FALSE
    failed_check <- file_type_error
  }

  if (checks_passed) {
    moves <- read.csv(x, header = TRUE)
    movement_from <- sp::SpatialPointsDataFrame(moves[, 1:2],
      data = moves[, c(1:2, 5:6)],
      proj4string = sp::CRS("+init=epsg:4326")
    )
    movement_to <- sp::SpatialPointsDataFrame(moves[, 3:4],
      data = moves[, 3:6],
      proj4string = sp::CRS("+init=epsg:4326")
    )
    movement_from <-
      suppressWarnings(sp::spTransform(movement_from, CRSobj = terra::crs(rast)))
    movement_to <-
      suppressWarnings(sp::spTransform(movement_to, CRSobj = terra::crs(rast)))
    move_from <- terra::vect(movement_from)
    move_to <- terra::vect(movement_to)
    cell_from <- terra::extract(rast, move_from, cells = TRUE)
    cell_to <- terra::extract(rast, move_to, cells = TRUE)
    rowcol_from <- terra::rowColFromCell(rast, cell_from[, 3])
    rowcol_to <- terra::rowColFromCell(rast, cell_to[, 3])
    movements <- data.frame(
      row_from = rowcol_from[, 1],
      col_from = rowcol_from[, 2],
      row_to = rowcol_to[, 1],
      col_to = rowcol_to[, 2],
      num_animals = moves$animals,
      date = moves$date
    )
    movements <- movements[!is.na(movements$row_from) &
      !is.na(movements$col_from) &
      !is.na(movements$row_to) &
      !is.na(movements$col_to) &
      !is.na(movements$num_animals) &
      !is.na(movements$date), ]
    movements <- movements[movements$num_animals > 0, ]
    movements$date <- paste("'", movements$date, "'", sep = "")
    movements$date <- lubridate::mdy(movements$date)
    duration <- lubridate::interval(start_date, end_date)
    movements <- movements[movements$date %within% duration, ]
    movements <- movements[order(movements$date, decreasing = FALSE), ]
    movements_dates <- as.character(movements$date)
    movements_r <- movements
    # subtract 1 from the movement index to account for r indexing starts at 1
    # and C++ starts at 0
    movements[, 1:4] <- movements[, 1:4] - 1
    movements <- unname(movements)
    movement2 <- as.matrix(movements[, 1:5])
    movement2 <- unname(movement2, force = TRUE)
    movement <- list()
    # movements_date
    for (i in seq_len(nrow(movement2))) {
      movement[[i]] <- movement2[i, 1:5]
    }

    movement <- unname(movement)
  }

  if (checks_passed) {
    outs <- list(checks_passed, movement, movements_dates, movements_r)
    names(outs) <-
      c("checks_passed", "movements", "movements_dates", "movements_r")
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- failed_check_list
    return(outs)
  }
}

random_seeds_file_checks <- function(x, number_of_iterations = 1) {
  checks_passed <- TRUE
  if (!all(file.exists(x))) {
    checks_passed <- FALSE
    failed_check <- file_exists_error
  }

  if (checks_passed && !all((tools::file_ext(x) %in% csv_list))) {
    checks_passed <- FALSE
    failed_check <- file_type_error
  }

  if (checks_passed) {
    random_seeds <- read.csv(x)
    if (NCOL(random_seeds) != 9 && NROW(random_seeds != number_of_iterations)) {
      checks_passed <- FALSE
      failed_check <- random_seeds_dimensions_error
    }
  }

  if (checks_passed) {
    outs <- list(checks_passed, random_seeds)
    names(outs) <- c("checks_passed", "random_seeds")
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- failed_check_list
    return(outs)
  }
}
