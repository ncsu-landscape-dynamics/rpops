## These functions are designed to improve data format checks and reduce copy and pasting of code across functions
Sys.setenv("R_TESTS"="")

initial_raster_checks <- function(x) {
  checks_passed <- TRUE
  
  if (!all(file.exists(x))) {
    checks_passed <- FALSE
    failed_check <- "file does not exist" 
  }
  
  if (checks_passed && !all((raster::extension(x) %in% c(".grd", ".tif", ".img")))) {
    checks_passed <- FALSE
    failed_check <- "file is not one of '.grd', '.tif', '.img'"
  }
  
  if (checks_passed) {
    r <- raster::stack(x)
    r <- raster::reclassify(r, matrix(c(NA,0), ncol = 2, byrow = TRUE), right = NA)
  }
  
  
  if (checks_passed) {
    outs <- list(checks_passed, r)
    names(outs) <- c('checks_passed', 'raster')
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- c('checks_passed', 'failed_check')
    return(outs)
  }
}

## adds checks to test for raster extent, resolution, and crs x2 being the raster already through initial checks for comparison
secondary_raster_checks <- function(x, x2) {
  checks_passed <- TRUE
  
  if (!all(file.exists(x))) {
    checks_passed <- FALSE
    failed_check <- "file does not exist" 
  }
  
  if (checks_passed && !all((raster::extension(x) %in% c(".grd", ".tif", ".img")))) {
    checks_passed <- FALSE
    failed_check <- "file is not one of '.grd', '.tif', '.img'"
  }
  
  if (checks_passed) {
    r <- raster::stack(x)
    r2 <- raster::reclassify(r, matrix(c(NA,0), ncol = 2, byrow = TRUE), right = NA)
    if (!(raster::extent(r2) == raster::extent(r))) {
      raster::extent(r2) <- raster::extent(r)
    }
    r <- r2
  }
  
  if (checks_passed && !(raster::extent(x2) == raster::extent(r))) {
    checks_passed <- FALSE
    failed_check <- "Extents of input rasters do not match. Ensure that all of your input rasters have the same extent"
  }
  
  if (checks_passed && !(raster::xres(x2) == raster::xres(r) && raster::yres(x2) == raster::yres(r))) {
    checks_passed <- FALSE
    failed_check <- "Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution"
  }
  
  if (checks_passed && !raster::compareCRS(r,x2)) {
    checks_passed <- FALSE
    failed_check <- "Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs"
  }
  
  if (checks_passed) {
    outs <- list(checks_passed, r)
    names(outs) <- c('checks_passed', 'raster')
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- c('checks_passed', 'failed_check')
    return(outs)
  }
}

treatment_checks <- function(treatment_stack, treatments_file, pesticide_duration, treatment_dates, pesticide_efficacy) {
  checks_passed <- TRUE
  
  if (checks_passed && length(treatments_file) != length(treatment_dates)) {
    checks_passed <- FALSE
    failed_check <- "Length of list for treatment dates and treatments_file must be equal"
  }
  
  if (checks_passed && length(pesticide_duration) != length(treatment_dates)) {
    checks_passed <- FALSE
    failed_check <- "Length of list for treatment dates and pesticide_duration must be equal"
  }
  
  if (checks_passed) {
    if (pesticide_duration[1] > 0) {
      treatment_maps <- list(raster::as.matrix(treatment_stack[[1]] * pesticide_efficacy))
    } else {
      treatment_maps <- list(raster::as.matrix(treatment_stack[[1]]))
    }
    
    if (raster::nlayers(treatment_stack) >= 2) {
      for(i in 2:raster::nlayers(treatment_stack)) {
        if (pesticide_duration[i] > 0) {
          treatment_maps[[i]] <- raster::as.matrix(treatment_stack[[i]] * pesticide_efficacy)
        } else {
          treatment_maps[[i]] <- raster::as.matrix(treatment_stack[[i]])
          
        }
      }
    }
  }
  
  if (checks_passed) {
    outs <- list(checks_passed, treatment_maps)
    names(outs) <- c('checks_passed', 'treatment_maps')
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- c('checks_passed', 'failed_check')
    return(outs)
  }
}

treatment_metric_checks <- function(treatment_method) {
  checks_passed <- TRUE
  
  if (!treatment_method %in% c("ratio", "all infected")) {
    checks_passed <- FALSE
    failed_check <- "treatment method is not one of the valid treatment options"
  }
  
  if (checks_passed) {
    outs <- list(checks_passed)
    names(outs) <- c('checks_passed')
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- c('checks_passed', 'failed_check')
    return(outs)
  }
}

metric_checks <- function(success_metric) {
  checks_passed <- TRUE
  
  if (success_metric == "quantity") {
    configuration <- FALSE
  } else if (success_metric == "quantity and configuration") {
    configuration <- TRUE
  } else if (success_metric == "odds ratio") {
    configuration <- FALSE
  } else if (success_metric == "residual error") {
    configuration <- FALSE
  } else {
    checks_passed <- FALSE
    failed_check <- "Success metric must be one of 'quantity', 'quantity and configuration', 'residual error', or 'odds ratio'"
  }
  
  if (checks_passed) {
    outs <- list(checks_passed, configuration)
    names(outs) <- c('checks_passed', 'configuration')
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- c('checks_passed', 'failed_check')
    return(outs)
  }
  
}

percent_checks <- function(percent_natural_dispersal) {
  checks_passed <- TRUE
  
  if(percent_natural_dispersal == 1.0) {
    use_anthropogenic_kernel = FALSE
  } else if (percent_natural_dispersal < 1.0  && percent_natural_dispersal >= 0.0) {
    use_anthropogenic_kernel = TRUE
  } else {
    checks_passed <- FALSE
    failed_check <- "Percent natural dispersal must be between 0.0 and 1.0"
  }
  
  if (checks_passed) {
    outs <- list(checks_passed, use_anthropogenic_kernel)
    names(outs) <- c('checks_passed', 'use_anthropogenic_kernel')
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- c('checks_passed', 'failed_check')
    return(outs)
  }
}

time_checks <- function(end_date, start_date, time_step, output_frequency) {
  checks_passed <- TRUE
  
  if (checks_passed && !(time_step %in% list("week", "month", "day"))) {
    checks_passed <- FALSE
    failed_check <- "Time step must be one of 'week', 'month' or 'day'"
  }
  
  if (checks_passed && (!is(end_date, "character") || !is(start_date, "character") || !is(as.Date(end_date, format="%Y-%m-%d"), "Date") || !is(as.Date(start_date, format="%Y-%m-%d"), "Date") || is.na(as.Date(end_date, format="%Y-%m-%d")) || is.na(as.Date(start_date, format="%Y-%m-%d")))){
    checks_passed <- FALSE
    failed_check <- "End time and/or start time not of type numeric and/or in format YYYY-MM-DD"
  }
  
  if (checks_passed && !(output_frequency %in% list("week", "month", "day", "year", "time_step"))) {
    checks_passed <- FALSE
    failed_check <- "Output frequency must be one of 'week', 'month' or 'day'"
  }
  
  if (checks_passed && output_frequency == "day") {
    if (time_step == "week" || time_step == "month") {
      checks_passed <- FALSE
      failed_check <- "Output frequency is more frequent than time_step. The minimum output_frequency you can use is the time_step of your simulation. You can set the output_frequency to 'time_step' to default to most frequent output possible"
    }
  }
  
  if (checks_passed && output_frequency == "week") {
    if (time_step == "month") {
      checks_passed <- FALSE
      failed_check <- "Output frequency is more frequent than time_step. The minimum output_frequency you can use is the time_step of your simulation. You can set the output_frequency to 'time_step' to default to most frequent output possible"
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
    }
  }
  
  if (checks_passed) {
    outs <- list(checks_passed, number_of_time_steps, number_of_years, number_of_outputs)
    names(outs) <- c('checks_passed', 'number_of_time_steps', 'number_of_years', 'number_of_outputs')
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- c('checks_passed', 'failed_check')
    return(outs)
  }

}

## check for making sure priors are in the proper format and output the mean where the mean serves as the starting point for the calibration
prior_checks <- function(priors) {
  checks_passed <- TRUE
  
  if (is(priors, "numeric") && length(priors) == 2) {
    priors <- matrix(priors, ncol = 2)
  } 
  
  if ((is(priors, "data.frame") | is(priors,"matrix")) && ncol(priors) == 2) {
    if (is(priors,"matrix") && base::nrow(priors) == 1) {
      start_priors <- priors[1]
      sd_priors <- priors[2]
    } else if (is(priors, "data.frame") && nrow(priors) == 1) {
      start_priors <- priors[[1]]
      sd_priors <- 0
    } else if ((is(priors, "data.frame") | is(priors,"matrix")) && base::nrow(priors) > 1) {
      names(priors) <- c('var', 'prob')
      start_priors <- priors$var[priors$prob == max(priors$prob)]
      if (length(start_priors) > 1) {
        start_priors <- mean(start_priors)
      } 
      sd_priors <- sd(priors$var)
    }
  } else {
    checks_passed <- FALSE
    failed_check <- "Incorrect format for priors"
  }
  
  if (checks_passed) {
    outs <- list(checks_passed, priors, start_priors, sd_priors)
    names(outs) <- c('checks_passed', 'priors', 'start_priors', 'sd_priors')
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- c('checks_passed', 'failed_check')
    return(outs)
  }
}

## helper for taking priors and calibartion to posteriors
bayesian_checks <- function(prior, start_priors, sd_priors, params, count, prior_weight, weight, step_size, bounds = c(0, Inf), round_to = 1, round_to_digits = 1) {
  checks_passed <- TRUE
  
  if ((is(prior, "matrix") && nrow(prior) == 1) || (is(prior, "numeric") && length(prior) == 1)) {
    priors <- round(rnorm(count, start_priors, sd_priors)/round_to, digits = round_to_digits)*round_to
    priors <- as.data.frame(table(priors))
    priors$priors <- as.numeric(as.character(priors$priors))
    priors$prob <- round(priors$Freq/count, digits = 3)
    priors$prob[priors$priors == bounds[2]] <- sum(priors$prob[priors$priors >= bounds[2]])
    priors$prob[priors$priors == bounds[1]] <- sum(priors$prob[priors$priors <= bounds[1]])
    priors <- priors[priors$prob > 0.000,]
    priors <- priors[priors$priors <= bounds[2] & priors$priors >= bounds[1],]
  } else if ((is(priors, "data.frame") | is(priors,"matrix")) && nrow(prior) > 1) {
    priors <- prior[ , 1:2]
    names(priors) <- c('priors', 'prob')
  }
  
  if (is(params, 'matrix')) {
    params <- data.frame(params)
    names(params) <- c('params', 'prob')
  }
  calibration_count <- length(params)
  calibrated_rates <- base::as.data.frame(table(params))
  calibrated_rates$params <- as.numeric(as.character(calibrated_rates$params))
  calibrated_rates$prob <- round(calibrated_rates$Freq/calibration_count, digits = 3)
  
  min_rate <- min(min(priors$priors), min(calibrated_rates$params))
  max_rate <- max(max(priors$priors), max(calibrated_rates$params))
  
  rates <- data.frame(rate = round(seq(min_rate, max_rate, step_size), digits = round_to_digits), 
                                   prior_probability = rep(0, length(seq(min_rate, max_rate, step_size))),
                                   calibrated_probability = rep(0, length(seq(min_rate, max_rate, step_size))),
                                   posterior_probability = rep(0, length(seq(min_rate, max_rate, step_size))))
  
  for (i in 1:base::nrow(rates)) {
    if (length(priors$prob[priors$priors == rates$rate[i]]) > 0) {
      rates$prior_probability[i] <- priors$prob[priors$priors == rates$rate[i]]
    }
    if (length(calibrated_rates$prob[calibrated_rates$params == rates$rate[i]]) > 0) {
      rates$calibrated_probability[i] <- calibrated_rates$prob[calibrated_rates$params == rates$rate[i]]
    }
  }
  rates$posterior_probability <- round(rates$prior_probability*prior_weight + rates$calibrated_probability * weight, digits = 3)
  posterior_rates <- rates[,c(1,4)]
  
  if (checks_passed) {
    outs <- list(checks_passed, rates, posterior_rates)
    names(outs) <- c('checks_passed', 'rates', 'posterior_rates')
    return(outs)
  } 
}


bayesian_MNN_checks <- function(prior_means, prior_cov_matrix, calibrated_means, calibrated_cov_matrix, prior_weight, weight) {
  checks_passed <- TRUE
  if (length(prior_means) == length(calibrated_means) && prior_weight > 0) {
    posterior_means <- prior_means * prior_weight + calibrated_means * weight
  } else if (prior_weight == 0){
    posterior_means <- calibrated_means
  } else {
    checks_passed <- FALSE 
    failed_check <- "There are not enough prior_means to communte the posterior means for all paramaters"
  }
  
  if (nrow(prior_cov_matrix) == nrow(calibrated_cov_matrix) && ncol(prior_cov_matrix) == ncol(calibrated_cov_matrix) && prior_weight > 0) {
    posterior_cov_matrix <- prior_cov_matrix * prior_weight + calibrated_cov_matrix * weight
  } else if (prior_weight == 0) {
    posterior_cov_matrix <- calibrated_cov_matrix
  } else { 
    checks_passed <- FALSE 
    failed_check <- "The prior covariance matrix is not the correct dimension to match the calibrated covariance matrix in order to computer the posterior covariance matrix"
  }
  
  if (checks_passed) {
    outs <- list(checks_passed, posterior_means, posterior_cov_matrix)
    names(outs) <- c('checks_passed', 'posterior_means', 'posterior_cov_matrix')
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- c('checks_passed', 'failed_check')
    return(outs)
  }
}

multispecies_checks <- function(species, infected_files, reproductive_rate, percent_natural_dispersal, natural_kernel_type, anthropogenic_kernel_type, 
                                natural_distance_scale, anthropogenic_distance_scale, natural_dir, natural_kappa, anthropogenic_dir, anthropogenic_kappa) {
  checks_passed <- TRUE
  
  if (checks_passed && length(species) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- "Length of list for species and infected_files must be equal"
  }
  
  if (checks_passed && length(reproductive_rate) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- "Length of list for infected_files and reproductive_rate must be equal"
  }
  
  if (checks_passed && length(percent_natural_dispersal) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- "Length of list for infected_files and percent_natural_dispersal must be equal"
  }
  
  if (checks_passed && length(natural_kernel_type) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- "Length of list for infected_files and natural_kernel_type must be equal"
  }
  
  if (checks_passed && length(anthropogenic_kernel_type) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- "Length of list for infected_files and anthropogenic_kernel_type must be equal"
  }
  
  if (checks_passed && length(natural_distance_scale) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- "Length of list for infected_files and natural_distance_scale must be equal"
  }
  
  if (checks_passed && length(anthropogenic_distance_scale) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- "Length of list for infected_files and anthropogenic_distance_scale must be equal"
  }
  
  if (checks_passed && length(natural_dir) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- "Length of list for infected_files and natural_dir must be equal"
  }
  
  if (checks_passed && length(natural_kappa) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- "Length of list for infected_files and natural_kappa must be equal"
  }
  
  if (checks_passed && length(anthropogenic_dir) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- "Length of list for infected_files and anthropogenic_dir must be equal"
  }
  
  if (checks_passed && length(anthropogenic_kappa) != length(infected_files)) {
    checks_passed <- FALSE
    failed_check <- "Length of list for infected_files and anthropogenic_kappa must be equal"
  }
  
  if (checks_passed) {
    outs <- list(checks_passed)
    names(outs) <- c('checks_passed')
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- c('checks_passed', 'failed_check')
    return(outs)
  }
}

movement_checks <- function(x, rast, start_date, end_date) {
  checks_passed <- TRUE
  
  if (!all(file.exists(x))) {
    checks_passed <- FALSE
    failed_check <- "file does not exist" 
  }
  
  if (checks_passed && !all((raster::extension(x) %in% c(".csv", ".txt")))) {
    checks_passed <- FALSE
    failed_check <- "file is not one of '.csv' or '.txt'"
  }
  
  if (checks_passed) {
    moves <- read.csv(x, header = TRUE)
    movement_from <- SpatialPointsDataFrame(moves[,1:2], data = moves, proj4string = CRS("+init=epsg:4326"))
    movement_to <- SpatialPointsDataFrame(moves[,3:4], data = moves, proj4string = CRS("+init=epsg:4326"))
    movement_from <- spTransform(movement_from, CRSobj = crs(rast))
    movement_to <- spTransform(movement_to, CRSobj = crs(rast))
    cell_from <-  raster::extract(rast, movement_from, cellnumbers = TRUE)
    cell_to <-  raster::extract(rast, movement_to, cellnumbers = TRUE)
    rowcol_from <- rowColFromCell(rast, cell_from[,1])
    rowcol_to <- rowColFromCell(rast, cell_to[,1])
    movements <- data.frame(row_from = rowcol_from[,1], col_from = rowcol_from[,2], row_to = rowcol_to[,1], col_to = rowcol_to[,2], num_animals = moves$animals, date = moves$date)
    movements <- movements[!is.na(movements$row_from) & !is.na(movements$col_from) & !is.na(movements$row_to) & !is.na(movements$col_to) & !is.na(movements$num_animals) & !is.na(movements$date),]
    movements <- movements[movements$num_animals > 0, ]
    movements$date <- paste("'", movements$date, "'", sep = "")
    movements$date <- lubridate::mdy(movements$date)
    duration <- lubridate::interval(start_date, end_date)
    movements <- movements[movements$date %within% duration,]
    movements <- movements[order(movements$date, decreasing = FALSE), ]
    movements_dates <- as.character(movements$date)
    movements_r <- movements
    movements[,1:4] <- movements[, 1:4] - 1 # subtract 1 from the movement index to account for r indexing starts at 1 and C++ starts at 0
    movements <- unname(movements)
    movement2 <- as.matrix(movements[, 1:5])
    movement2 <- unname(movement2, force = TRUE)
    movement <- list()
    # movements_date
    for (i in 1:nrow(movement2)){
      movement[[i]] <- movement2[i,1:5]
    }

    movement <- unname(movement)
  }
  
  if (checks_passed) {
    outs <- list(checks_passed, movement, movements_dates, movements_r)
    names(outs) <- c('checks_passed', 'movements', 'movements_dates', 'movements_r')
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- c('checks_passed', 'failed_check')
    return(outs)
  }
}

parameter_checks <- function(n, parameter_means, parameter_cov_matrix) {
  checks_passed <- TRUE
  
  if (nrow(parameter_cov_matrix) != 6 | ncol(parameter_cov_matrix) != 6) {
    checks_passed <- FALSE
    failed_check <- "parameter covariance matrix is not 6 x 6" 
  }
  
  if (length(parameter_means) != 6) {
    checks_passed <- FALSE
    failed_check <- "parameter means is not a vector of length 6" 
  }
  
  if(checks_passed) {
    parameters <- data.frame(MASS::mvrnorm(n, parameter_means, parameter_cov_matrix))
  }
  
  if (checks_passed) {
    outs <- list(checks_passed, parameters)
    names(outs) <- c('checks_passed', 'parameters')
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- c('checks_passed', 'failed_check')
    return(outs)
  }
}

quarantine_checks <- function() {
  checks_passed <- TRUE
  
  if (nrow(parameter_cov_matrix) != 6 | ncol(parameter_cov_matrix) != 6) {
    checks_passed <- FALSE
    failed_check <- "parameter covariance matrix is not 6 x 6" 
  }
  
  if(checks_passed) {
    parameters <- data.frame(MASS::mvrnorm(n, parameter_means, parameter_cov_matrix))
  }
  
  if (checks_passed) {
    outs <- list(checks_passed, parameters)
    names(outs) <- c('checks_passed', 'parameters')
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    names(outs) <- c('checks_passed', 'failed_check')
    return(outs)
  }
}