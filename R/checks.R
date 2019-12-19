## These functions are designed to improve data format checks and reduce copy and pasting of code across functions

initial_raster_check <- function(x) {
  checks_passed <- TRUE
  
  if (!file.exists(x)) {
    checks_passed <- FALSE
    failed_check <- "file does not exist" 
  }
  
  if (!(raster::extension(x) %in% c(".grd", ".tif", ".img"))) {
    checks_passed <- FALSE
    failed_check <- "file is not one of '.grd', '.tif', '.img'"
  }
  
  if (checks_passed) {
    r<- raster::raster(x)
    r <- raster::reclassify(r, matrix(c(NA,0), ncol = 2, byrow = TRUE), right = NA)
  }
  
  if (checks_passed) {
    outs <- list(checks_passed, r)
    return(outs)
  } else {
    outs <- list(checks_passed, failed_check)
    return(outs)
  }
}





time_checks <- function(end_date, start_date, time_step, output_frequency) {
  checks_passed <- TRUE
  
  if (!(time_step %in% list("week", "month", "day"))) {
    checks_passed <- FALSE
    failed_check <- "Time step must be one of 'week', 'month' or 'day'"
  }
  
  if (class(end_date) != "character" || class(start_date) != "character" || class(as.Date(end_date, format="%Y-%m-%d")) != "Date" || class(as.Date(start_date, format="%Y-%m-%d")) != "Date" || is.na(as.Date(end_date, format="%Y-%m-%d")) || is.na(as.Date(start_date, format="%Y-%m-%d"))){
    checks_passed <- FALSE
    failed_check <- "End time and/or start time not of type numeric and/or in format YYYY-MM-DD"
  }
  
  if (!(output_frequency %in% list("week", "month", "day", "year", "time_step"))) {
    checks_passed <- FALSE
    failed_check <- "Time step must be one of 'week', 'month' or 'day'"
  }
  
  if (output_frequency == "day") {
    if (time_step == "week" || time_step == "month") {
      checks_passed <- FALSE
      failed_check <- "Output frequency is more frequent than time_step. The minimum output_frequency you can use is the time_step of your simulation. You can set the output_frequency to 'time_step' to default to most frequent output possible"
    }
  }
  
  if (output_frequency == "week") {
    if (time_step == "month") {
      checks_passed <- FALSE
      failed_check <- "Output frequency is more frequent than time_step. The minimum output_frequency you can use is the time_step of your simulation. You can set the output_frequency to 'time_step' to default to most frequent output possible"
    }
  }
  
  
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

