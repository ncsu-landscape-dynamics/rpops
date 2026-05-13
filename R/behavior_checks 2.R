# Behavior module input validation helpers
# Mirrors the pattern of checks.R

#' Check behavior type parameters
#'
#' @param type_names character vector of grower type names
#' @param type_probs numeric vector of probabilities (must sum to 1)
#' @param spatial_structure one of 'random', 'clustered', 'empirical'
#' @param cluster_model one of 'gaussian', 'exponential' (only checked if clustered)
#' @return list with checks_passed (logical) and failed_check (character) if FALSE
behavior_type_checks <- function(type_names, type_probs, spatial_structure,
                                   cluster_model = "gaussian") {
  checks_passed <- TRUE
  failed_check <- NULL

  if (length(type_names) == 0) {
    checks_passed <- FALSE
    failed_check <- type_names_empty_error
    return(list(checks_passed = checks_passed, failed_check = failed_check))
  }

  if (length(type_names) != length(type_probs)) {
    checks_passed <- FALSE
    failed_check <- type_names_length_error
    return(list(checks_passed = checks_passed, failed_check = failed_check))
  }

  if (abs(sum(type_probs) - 1) > 1e-6) {
    checks_passed <- FALSE
    failed_check <- type_probs_sum_error
    return(list(checks_passed = checks_passed, failed_check = failed_check))
  }

  if (!spatial_structure %in% behavior_structure_list) {
    checks_passed <- FALSE
    failed_check <- behavior_structure_error
    return(list(checks_passed = checks_passed, failed_check = failed_check))
  }

  if (spatial_structure == "clustered" && !cluster_model %in% cluster_model_list) {
    checks_passed <- FALSE
    failed_check <- cluster_model_error
    return(list(checks_passed = checks_passed, failed_check = failed_check))
  }

  return(list(checks_passed = checks_passed, failed_check = failed_check))
}

#' Check decision dates
#'
#' @param decision_dates character vector of YYYY-MM-DD dates
#' @param start_date simulation start date string
#' @param end_date simulation end date string
#' @return list with checks_passed (logical) and failed_check (character) if FALSE
decision_date_checks <- function(decision_dates, start_date, end_date) {
  checks_passed <- TRUE
  failed_check <- NULL

  parsed <- tryCatch(
    as.Date(decision_dates, format = "%Y-%m-%d"),
    error = function(e) NULL
  )

  if (is.null(parsed) || any(is.na(parsed))) {
    checks_passed <- FALSE
    failed_check <- decision_dates_format_error
    return(list(checks_passed = checks_passed, failed_check = failed_check))
  }

  sim_start <- as.Date(start_date, format = "%Y-%m-%d")
  sim_end   <- as.Date(end_date,   format = "%Y-%m-%d")

  if (any(parsed < sim_start) || any(parsed > sim_end)) {
    checks_passed <- FALSE
    failed_check <- decision_dates_bounds_error
    return(list(checks_passed = checks_passed, failed_check = failed_check))
  }

  return(list(checks_passed = checks_passed, failed_check = failed_check))
}

#' Check management unit matrix integrity
#'
#' @param management_unit_matrix integer matrix of unit IDs
#' @param host_matrix integer matrix of the same dimensions
#' @return list with checks_passed (logical) and failed_check (character) if FALSE
management_unit_checks <- function(management_unit_matrix, host_matrix) {
  checks_passed <- TRUE
  failed_check <- NULL

  if (!identical(dim(management_unit_matrix), dim(host_matrix))) {
    checks_passed <- FALSE
    failed_check <- management_unit_dim_error
    return(list(checks_passed = checks_passed, failed_check = failed_check))
  }

  vals <- as.integer(management_unit_matrix[!is.na(management_unit_matrix)])
  if (any(vals < 0) || !all(vals == floor(vals))) {
    checks_passed <- FALSE
    failed_check <- management_unit_id_error
    return(list(checks_passed = checks_passed, failed_check = failed_check))
  }

  return(list(checks_passed = checks_passed, failed_check = failed_check))
}

#' Check trust/learning dynamics configuration
#'
#' @param trust_config list with learning_rate, success_threshold, memory_length
#' @return list with checks_passed (logical) and failed_check (character) if FALSE
trust_config_checks <- function(trust_config) {
  checks_passed <- TRUE
  failed_check <- NULL

  lr <- trust_config$learning_rate
  if (is.null(lr) || !is.numeric(lr) || lr <= 0 || lr > 1) {
    checks_passed <- FALSE
    failed_check <- learning_rate_error
    return(list(checks_passed = checks_passed, failed_check = failed_check))
  }

  st <- trust_config$success_threshold
  if (is.null(st) || !is.numeric(st) || st < 0 || st > 1) {
    checks_passed <- FALSE
    failed_check <- success_threshold_behavior_error
    return(list(checks_passed = checks_passed, failed_check = failed_check))
  }

  ml <- trust_config$memory_length
  if (is.null(ml) || !is.numeric(ml) || ml < 1 || ml != floor(ml)) {
    checks_passed <- FALSE
    failed_check <- memory_length_error
    return(list(checks_passed = checks_passed, failed_check = failed_check))
  }

  return(list(checks_passed = checks_passed, failed_check = failed_check))
}
