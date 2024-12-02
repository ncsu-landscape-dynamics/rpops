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
#' @export config_qad_list list of all data necessary used to set up the PoPS
#' `quantity_allocation_disagreement_lite` function
#'
#' @importFrom terra app rast xres yres classify extract ext as.points ncol nrow
#' nlyr rowFromCell colFromCell values as.matrix rowFromCell colFromCell crs
#' rowColFromCell global vect
#' 
config_qad_lite <- function(reference,
                         comparison,
                         mask,
                         use_configuration = FALSE,
                         use_distance = FALSE,
                         output_folder_path = "") {
  # Helper function to adjust file paths based on output_folder_path
  adjust_path <- function(file, folder) {
    if (folder != "")
      file.path(folder, basename(file))
    else
      file
  }
  
  # Generate filenames for intermediate outputs
  config_qad_rds_fn <- adjust_path(gsub(".tif", "_qad.rds", comparison), output_folder_path)
  reference_masked_fn <- adjust_path(gsub(".tif", "_qad_masked.tif", reference), output_folder_path)
  comparison_masked_fn <- adjust_path(gsub(".tif", "_qad_masked.tif", comparison), output_folder_path)
  reference_classified_fn <- adjust_path(gsub(".tif", "_qad_reclassed.tif", reference),
                                         output_folder_path)
  comparison_classified_fn <- adjust_path(gsub(".tif", "_qad_reclassed.tif", comparison),
                                          output_folder_path)
  
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
  if (file.exists(reference_masked_fn)) {
    ref <- terra::rast(reference_masked_fn)
  } else {
    ref <- terra::mask(terra::rast(reference), terra::rast(mask))
  }
  
  # Process comparison (masked only)
  if (file.exists(comparison_masked_fn)) {
    comp <- terra::rast(comparison_masked_fn)
  } else {
    comp <- terra::mask(terra::rast(comparison), terra::rast(mask))
  }
  
  # Process reference (masked + reclassified)
  if (file.exists(reference_classified_fn)) {
    reference <- terra::rast(reference_classified_fn)
  } else {
    reference <- terra::classify(ref, matrix(
      c(NA, 0, 2, 1, Inf, 1, 0, 0.99, 0),
      ncol = 3,
      byrow = TRUE
    ), right = FALSE)
  }
  
  # Process comparison (masked + reclassified)
  if (file.exists(comparison_classified_fn)) {
    comparison <- terra::rast(comparison_classified_fn)
  } else {
    comparison <- terra::classify(comp, matrix(
      c(1, Inf, 1, 0, 1, 0, NA, 0, 0),
      ncol = 3,
      byrow = TRUE
    ), right = FALSE)
  }
  
  # Extract values and calculate totals
  ref_values <- terra::values(reference)
  comp_values <- terra::values(comparison)
  
  positives_in_reference <- sum(ref_values == 1, na.rm = TRUE)
  positives_in_comparison <- sum(comp_values == 1, na.rm = TRUE)
  
  positives_in_ref <- sum(terra::values(ref), na.rm = TRUE)
  positives_in_comp <- sum(terra::values(comp), na.rm = TRUE)
  
  # Save the configuration object
  config_qad_lite <- list(
    reference_fn = reference_classified_fn,
    comparison_fn = comparison_classified_fn,
    use_distance = use_distance,
    use_configuration = use_configuration,
    ref_fn = reference_masked_fn,
    comp_fn = comparison_masked_fn,
    ref_values = ref_values,
    comp_values = comp_values,
    positives_in_reference = positives_in_reference,
    positives_in_comparison = positives_in_comparison,
    positives_in_ref = positives_in_ref,
    positives_in_comp = positives_in_comp
  )
  
  saveRDS(config_qad_lite, config_qad_rds_fn, compress = TRUE)
  message(sprintf("Configuration file '%s' has been created.", config_qad_rds_fn))
}

