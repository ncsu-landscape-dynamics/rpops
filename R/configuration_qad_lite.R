#' Create PoPS configuration file for quantity allocation and disagreement lite function
#'
#' Pre-process input data for the quantity, allocation, and disagreement lite function.
#' This function separates pre-processing steps from validation to optimize memory usage.
#'
#' @param reference Path to the raster with ground truth data. Values > 1 are
#' reclassified to 1, values < 1 to 0, and NA remains unchanged (except for RMSE).
#' @param comparison Path to the raster with simulated data. Values > 1 are
#' reclassified to 1, and values < 1 to 0 (except for RMSE).
#' @param mask Path to the raster mask.
#' @param use_configuration Boolean to use configuration disagreement for model
#' comparison. Default is FALSE.
#' @param use_distance Boolean to compare distances between simulations and
#' observations. Default is FALSE.
#' @param output_folder_path this is the full path with either / or \
#' (e.g., "C:/user_name/desktop/pops_sod_2020_2023/outputs/"). If not provided,
#' the config_qad_rds, masked and reclassified rasters will be exported to the
#' same directories as the input reference and comparison folders.
#'
#' @export config_qad_lite list of all data necessary used to set up the PoPS
#' `quantity_allocation_disagreement_lite` function
#'
#' @importFrom terra app rast xres yres classify extract ext as.points ncol nrow
#' nlyr rowFromCell colFromCell values as.matrix rowFromCell colFromCell crs
#' rowColFromCell global vect


reference <- "C:/Users/blaginh/Desktop/pops_runs/input_data/res_1000m/slf/slf_2014_2020.tif"
comparison <- "C:/Users/blaginh/Desktop/pops_runs/results/6o1v3i7c_1000m_2020.tif"
mask <- "C:/Users/blaginh/Desktop/pops_runs/input_data/res_1000m/base/mask_1000m.tif"
use_configuration <- TRUE
use_distance <- TRUE
use_rmse <- TRUE
output_folder_path <- "C:/Users/blaginh/Desktop/pops_runs/results/qad_outputs/"

configuration_qad_lite <- function(reference,
                                   comparison,
                                   mask,
                                   use_configuration = FALSE,
                                   use_distance = FALSE,
                                   use_rmse = FALSE,
                                   output_folder_path = "") {
  
  # Helper function to adjust file paths based on output_folder_path
  adjust_path <- function(file, folder) {
    if (folder != "")
      file.path(folder, basename(file))
    else
      file
  }
  
  # Create an empty list for configuration
  config_qad_lite <- list()
  
  # Generate filenames for intermediate outputs
  config_qad_rds_fn <- adjust_path(gsub(".tif", "_qad.rds", comparison), output_folder_path)
  config_qad_lite$ref <- adjust_path(gsub(".tif", "_qad_masked.tif", reference),
                                                     output_folder_path)
  positives_in_ref_fn <- adjust_path(gsub(".tif", "_qad_masked_positives.rds", reference),
                                                     output_folder_path)
  config_qad_lite$reference <- adjust_path(gsub(".tif", "_qad_reclassed.tif", reference),
                                                         output_folder_path)
  ref_values_fn <- adjust_path(gsub(".tif", "_qad_reclassed_vals.rds", reference),
                                                        output_folder_path)
  config_qad_lite$ref_points <- adjust_path(gsub(".tif", "_qad_masked_positive_pts.gpkg", reference),
                                               output_folder_path)
  positives_in_reference_fn <- adjust_path(gsub(".tif", "_qad_reclassed_positives.rds", reference),
                                                           output_folder_path)
  config_qad_lite$comp <- adjust_path(gsub(".tif", "_qad_masked.tif", comparison),
                                                      output_folder_path)
  config_qad_lite$comparison <- adjust_path(gsub(".tif", "_qad_reclassed.tif", comparison),
                                                          output_folder_path)
  config_qad_lite$use_configuration <- use_configuration
  config_qad_lite$use_distance <- use_distance
  config_qad_lite$use_rmse <- use_rmse
  
  # Check for existing config file
  if (file.exists(config_qad_rds_fn)) {
    message(
      sprintf(
        "Configuration file '%s' already exists. Skipping processing.",
        config_qad_rds_fn
      )
    )
    return(NULL)
  }
  
  # Process reference (masked only)
  if (file.exists(config_qad_lite$ref)) {
    ref <- terra::rast(config_qad_lite$ref)
    config_qad_lite$positives_in_ref <- readRDS(positives_in_ref_fn)
  } else {
    ref <- terra::mask(terra::rast(reference), terra::rast(mask))
    terra::writeRaster(
      ref,
      config_qad_lite$ref,
      overwrite = TRUE,
      gdal = "COMPRESS=ZSTD",
      datatype = "INT1U"
    )
    config_qad_lite$positives_in_ref <- terra::global(ref, fun = "sum", na.rm = TRUE)$sum
    saveRDS(
      config_qad_lite$positives_in_ref,
      positives_in_ref_fn,
      compress = TRUE
    )
  }
  
  comp <- terra::mask(terra::rast(comparison), terra::rast(mask))
  terra::writeRaster(
    comp,
    config_qad_lite$comp,
    overwrite = TRUE,
    datatype = "INT2U",
    gdal = "COMPRESS=ZSTD"
  )
  config_qad_lite$positives_in_comp <- terra::global(comp, fun = "sum", na.rm = TRUE)$sum
  
  if (use_distance) {
    # Directly classify `comparison` raster for point extraction
    sim_points <- terra::as.points(terra::classify(comp, matrix(
      c(-Inf, 0, NA), ncol = 3, byrow = TRUE
    )), na.rm = TRUE)
    names(sim_points) <- "data"
    
    # Directly classify `reference` raster for point extraction
    if (!file.exists(config_qad_lite$ref_points)) {
      ref_points <- terra::as.points(terra::classify(ref, matrix(
        c(-Inf, 0, NA), ncol = 3, byrow = TRUE
      )), na.rm = TRUE)
      names(ref_points) <- "data"
      terra::writeVector(ref_points, config_qad_lite$ref_points, overwrite = TRUE)
    }
    ref_points <- terra::vect(config_qad_lite$ref_points)
    dist <- terra::distance(ref_points, sim_points)
    rm(ref_points, sim_points); gc()
    config_qad_lite$distance_differences <- apply(dist, 2, min)
    rm(dist); gc()
    all_distances <- function(distance_differences) {
      distance_differences <- round(sqrt(sum(distance_differences^2)), digits = 0)
      return(distance_differences)
      }
    config_qad_lite$distance_differences <- lapply(config_qad_lite$distance_differences, all_distances)
    config_qad_lite$distance_differences <- unlist(config_qad_lite$distance_differences, recursive = TRUE, use.names = TRUE)
  }
  
  if (use_rmse) {
    totals <- ref + comp
    obs_points <- terra::as.points(terra::classify(totals, matrix(c(-Inf, 0, NA),
                                                                  ncol = 3,
                                                                  byrow = TRUE)),
                                   na.rm = TRUE)
    names(obs_points) <- "data"
    config_qad_lite$actual_predicted <- terra::extract(c(ref, comp), obs_points)[, -1]
    names(config_qad_lite$actual_predicted) <- c("actual", "predicted")
    rm(obs_points, totals); gc()
  }
  
  # Process reference (masked + reclassified)
  if (file.exists(config_qad_lite$reference)) {
    reference <- terra::rast(config_qad_lite$reference)
    config_qad_lite$ref_values <- readRDS(config_qad_lite$ref_values)
    config_qad_lite$positives_in_reference <- readRDS(positives_in_reference_fn)
  } else {
    reference <- terra::classify(ref, matrix(
      c(NA, 0, 2, 1, Inf, 1, 0, 0.99, 0),
      ncol = 3,
      byrow = TRUE
    ), right = FALSE)
    terra::writeRaster(
      reference,
      config_qad_lite$reference,
      overwrite = TRUE,
      gdal = "COMPRESS=ZSTD",
      datatype = "FLT4S"
    )
    config_qad_lite$ref_values <- terra::values(reference, na.rm = TRUE)
    saveRDS(
      config_qad_lite$ref_values,
      ref_values_fn,
      compress = TRUE
    )
    config_qad_lite$positives_in_reference <- sum(config_qad_lite$ref_values == 1, na.rm = TRUE)
    saveRDS(
      config_qad_lite$positives_in_reference,
      positives_in_reference_fn,
      compress = TRUE
    )
  }
  rm(ref); gc()
  
  # Process comparison (masked + reclassified)
  if (file.exists(config_qad_lite$comparison)) {
    comparison <- terra::rast(config_qad_lite$comparison)
  } else {
    comparison <- terra::classify(comp, matrix(
      c(1, Inf, 1, 0, 1, 0, NA, 0, 0),
      ncol = 3,
      byrow = TRUE
    ), right = FALSE)
    terra::writeRaster(
      comparison,
      config_qad_lite$comparison,
      overwrite = TRUE,
      gdal = "COMPRESS=ZSTD",
      datatype = "INT1U"
    )
  }
  
  rm(comp); gc()
  
  # Extract values and calculate totals
  config_qad_lite$comp_values <- terra::values(comparison, na.rm = TRUE)
  config_qad_lite$positives_in_comparison <- sum(config_qad_lite$comp_values == 1, na.rm = TRUE)
  
  rm(reference, comparison); gc()
  config_qad_lite$output_folder_path <- output_folder_path
  
  # Save the configuration object
  saveRDS(config_qad_lite,
          config_qad_rds_fn,
          compress = TRUE)
  message(
    sprintf(
      "Configuration file '%s' has been created.",
      config_qad_rds_fn
    )
  )
}
