#' Validates the accuracy of the calibrated model on out of sample data.
#'
#' This function uses the quantity, allocation, and configuration disagreement
#' to validate the model across the landscape using the parameters from the
#' calibrate function. Ideally the model is calibrated with 2 or more years of
#' data and validated for the last year or if you have 6 or more years of data
#' then the model can be validated for the final 2 years.
#'
#' @inheritParams pops_simulate
#'
#' @importFrom terra app classify extract as.points project values crs vect
#' @importFrom Metrics rmse
#' @importFrom utils write.csv read.table read.csv
#'
#' @return a data frame of statistical measures of model performance.
#' @export
#'
validate <- function(config_rds_file) {
  config <- readRDS(config_rds_file)
  if (!is.null(config$failure)) {
    stop(config$failure)
  }
  filelist <- list.files(file.path(config$output_path), pattern = "pops_ensemble_ts*")
  base_rast <- terra::rast(file.path(config$input_path, config$host_files[[1]]))
  if (config$use_mask) {
    mask <- base_rast[[1]]
    terra::values(mask) <- config$mask_matrix
  } else {
    mask <- NULL
  }

  if (config$county_level_infection_data) {
    reference <- terra::vect(file.path(config$input_path, config$infected_val_cal_file))
  } else {
    reference <- terra::rast(file.path(config$input_path, config$infected_val_cal_file))
  }

  vals <- lapply(seq_len(length(filelist)), function(j) {
    if (config$county_level_infection_data) {
      compare_vect <- reference[, c(1, (j + 1))]
      names(compare_vect) <- c("FIPS", "reference")
    }
    comparison <- terra::rast(file.path(config$output_path, filelist[j]))
    y <- lapply(1:config$number_of_iterations, function(i) {
      if (config$county_level_infection_data) {
        compare_vect$comparison <- terra::extract(comparison, reference, fun = "sum")[, 2]
        ad <- calculated_stats_county_level(compare_vect)
        ad <- calculated_stats_county_level(compare_vect)
        ad$quantity_disagreement <- 0
        ad$allocation_disagreement <- 0
        ad$allocation_disagreement <- 0
        ad$configuration_disagreement <- 0
        ad$distance_difference <- 0
        ad
      } else {
        quantity_allocation_disagreement(reference[[j]], comparison[[i]],
                                         use_configuration = config$use_configuration,
                                         mask = mask,
                                         use_distance = config$use_distance)
      }
    })
    s <- as.data.frame(do.call(rbind, y))
    s$output_ts <- j
    write.csv(s, file = paste(config$output_path, "val_metrics_output_step_", j, ".csv", sep = ""),
              row.names = FALSE)
    means <- colMeans(s)
    return(means)
  })

  means <- as.data.frame(do.call(rbind, vals))
  write.csv(means,
            file = paste(config$output_path, "validation_means_per_timestep.csv", sep = ""),
            row.names = FALSE)
  return(means)
}
