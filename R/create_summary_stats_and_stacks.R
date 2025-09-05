#' Create summary stats and raster summaries
#'
#' This function takes the outputs from PoPS lite and calculates the mean, max, min, and standard
#' deviation of each cell and calculates the summary statistics for the study area.
#'
#' @param config the config used for running PoPS lite for the simulations.
#' @return creates and writes raster stacks from all pops_lite runs in the outputs folder. Also
#' reates and writes mean, standard deviation, median, min, and max runs, and summary statistics.
#' @export
#'

create_summary_stats_and_stacks <- function(config) {
  raster_template <- rast(config$host_file_list[[1]])[[1]]
  filelist <- list.files(file.path(config$output_folder_path), pattern = "pops_output*")

  inf_indices <- lapply(seq_len(config$number_of_outputs), function(i) {
    seq(i, 2 * config$number_of_outputs *config$number_of_iterations, 2 * config$number_of_outputs)
  })
  area_indices <- lapply((3 + seq_len(config$number_of_outputs)), function(i) {
    seq(i, 2 * config$number_of_outputs *config$number_of_iterations, 2 * config$number_of_outputs)
  })
  all_indices <-  c(inf_indices, area_indices)

  inf_names <-  unlist(lapply(seq_len(config$number_of_outputs), function(i) {
    paste0("number_infecteds_y", i)
  }))
  area_names <-  unlist(lapply(seq_len(config$number_of_outputs), function(i) {
    paste0("area_infected_y", i)
  }))
  all_names <-  c(inf_names, area_names)

  rasts <- lapply(1:config$number_of_outputs, function(j) {
    y <- lapply(1:config$number_of_iterations, function(i){
      file <- readRDS(file.path(config$output_folder_path, filelist[i]))
      values(raster_template) <- file$host_pools[[1]]$infected[[j]]
      raster_template
    })
    return(rast(y))
  })

  x <- unlist(lapply(1:config$number_of_iterations, function(r) {
    file <- readRDS(file.path(config$output_folder_path, filelist[r]))
    c(file$number_infected, file$area_infected)
  }))

  all_stats <- data.frame(matrix(ncol = config$number_of_outputs * 2, nrow = config$number_of_iterations))
  for (i in seq_len(length(all_indices))) {
    all_stats[, i] <- x[all_indices[[i]]]
  }
  names(all_stats) <- all_names

  all_means <- rast(lapply(rasts, terra::mean))
  sd_s <- rast(lapply(rasts, terra::stdev))
  which_median <- function(x) which.min(abs(x - median(x)))
  median_run_index <- which_median(all_stats$number_infecteds_y1)
  min_run_index <- which.min(all_stats$number_infecteds_y1)
  max_run_index <- which.max(all_stats$number_infecteds_y1)
  for (i in seq_len(config$number_of_outputs)) {
    if (i ==1) {
      median_run <- rasts[[i]][[median_run_index]]
      min_run <- rasts[[i]][[min_run_index]]
      max_run <- rasts[[i]][[max_run_index]]
    } else {
      median_run <- c(median_run, rasts[[i]][[median_run_index]])
      min_run <- c(min_run, rasts[[i]][[min_run_index]])
      max_run <- c(max_run, rasts[[i]][[max_run_index]])
    }
  }

  writeRaster(median_run, file.path(config$output_folder_path, paste0("pops_median.tif")), overwrite = TRUE, gdal = c("COMPRESS=NONE"))
  writeRaster(min_run, file.path(config$output_folder_path, paste0("pops_min.tif")), overwrite = TRUE, gdal = c("COMPRESS=NONE"))
  writeRaster(max_run, file.path(config$output_folder_path, paste0("pops_max.tif")), overwrite = TRUE, gdal = c("COMPRESS=NONE"))
  writeRaster(all_means, file.path(config$output_folder_path, paste0("pops_mean.tif")), overwrite = TRUE, gdal = c("COMPRESS=NONE"))
  writeRaster(sd_s, file.path(config$output_folder_path, paste0("pops_sd.tif")), overwrite = TRUE, gdal = c("COMPRESS=NONE"))
  write.csv(all_stats, (file.path(config$output_folder_path, paste0("all_stats.csv"))))
  summary_stats <-  data.frame(t(c(colMeans(all_stats), sapply(all_stats, sd))))
  write.csv(summary_stats, file.path(config$output_folder_path, paste0("all_stats.csv")))
  return(all_means)
}