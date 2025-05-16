#' Calculates quantity and allocation disagreement between two raster datasets.
#' `quantity_allocation_disagreement_lite` functions similarly to
#' `quantity_allocation_disagreement` but uses pre-processed inputs from
#' `config_qad_lite`.
#'
#' `config_qad_lite` handles pre-processing, preparing the reference and
#' comparison raster files, and saves them as an RDS file. This RDS is then
#' used by `quantity_allocation_disagreement_lite` to compute metrics efficiently.
#'
#' Separating pre-processing from metric calculation minimizes memory usage
#' in `quantity_allocation_disagreement_lite`, making it suitable for
#' large extent, fine-resolution rasters.

#'
#' Uses quantity and allocation disagreement metrics by Pontius and Millones
#' (2014) and omission and commission errors of comparing a modeled raster
#' data set to a reference raster data set.
#'
#' @param config_qad_lite_file the raster with ground truth data. For all metrics expect
#' RMSE these values are reclassified with values > 1 becoming 1, values < 1
#' @param use_configuration Boolean if you want to use configuration
#' disagreement for comparing model runs. NOTE: overrides boolean provided in
#' the quantity, allocation, and disagreement config file.
#' @param use_distance Boolean if you want to compare distance between
#' simulations and observations. NOTE: overrides boolean provided in
#' the quantity, allocation, and disagreement config file.
#' @param new_dirs_path (Default = `NULL`) Specify a new root directory to update
#' input and output file paths in the `config_qad_lite_file`. This helps adapt the file paths
#' when running the configuration on a different workstation with a different root path,
#' without needing to recreate the entire `config_qad_lite_file`.
#' If `new_dirs_path` is `NULL` (default), the original paths in `config_qad_lite_file` are used.
#' When specified, `new_dirs_path` should point to the top-level folder containing
#' the input files and output folder. The folder structure under this top-level
#' directory must match the structure in the original `config_qad_lite_file` If no match is
#' found, the original input file paths and output folder remain unchanged.

#' @importFrom landscapemetrics lsm_c_np lsm_c_enn_mn lsm_c_para_mn lsm_c_lpi
#' @importFrom terra cells xres ncol nrow yres ext compareGeom vect
#' @importFrom Metrics rmse
#'
#' @return A data frame with spatial configuration metrics. Particularly
#' quantity, allocation, and total disagreement,  omission and commission, and
#' directional disagreement where directional disagreement.
#'
#' @export

quantity_allocation_disagreement_lite <-
  function(config_qad_lite_file, new_dirs_path = NULL, use_configuration = NULL,
           use_distance = NULL) {

    config_qad_lite <- readRDS(config_qad_lite_file)

    if (!is.null(use_configuration)) {
    config_qad_lite$use_configuration <- use_configuration
    }

    if (!is.null(use_distance)) {
      config_qad_lite$use_distance <- use_distance
    }

    if (!is.null(new_dirs_path)) {
      config_qad_lite <- update_config_paths(config_qad_lite, new_dirs_path)
    }

    reference <- terra::rast(config_qad_lite$reference)
    comparison <- terra::rast(config_qad_lite$comparison)
    ref <- terra::rast(config_qad_lite$ref)
    comp <- terra::rast(config_qad_lite$comp)

    if (config_qad_lite$use_configuration) {
      # Helper function to extract metric value for class 1 or return 0
      get_metric_value <- function(metric, class = 1) {
        if (any(unique(metric$class) %in% class)) {
          metric$value[metric$class == class]
        } else {
          0
        }
      }

      # Calculate metrics for the comparison data
      np_comp <- if (config_qad_lite$positives_in_comparison == 0)
        0
      else {
        np <- landscapemetrics::lsm_c_np(comparison, directions = 8)
        get_metric_value(np)
      }

      enn_mn_comp <- if (np_comp > 1) {
        enn_mn <-
          suppressWarnings(landscapemetrics::lsm_c_enn_mn(comparison, directions = 8,
                                                          verbose = TRUE))
        get_metric_value(enn_mn)
      } else {
        0
      }

      para_mn_comp <- if (np_comp > 0) {
        para_mn <- landscapemetrics::lsm_c_para_mn(comparison, directions = 8)
        get_metric_value(para_mn)
      } else {
        0
      }

      lpi_comp <- if (np_comp > 0) {
        lpi <- landscapemetrics::lsm_c_lpi(comparison, directions = 8)
        get_metric_value(lpi)
      } else {
        0
      }

      # Calculate metrics for the reference data
      np_ref <- {
        np <- landscapemetrics::lsm_c_np(reference, directions = 8)
        get_metric_value(np)
      }

      enn_mn_ref <- if (np_ref > 1) {
        enn_mn <- landscapemetrics::lsm_c_enn_mn(reference, directions = 8, verbose = TRUE)
        get_metric_value(enn_mn)
      } else {
        0
      }

      para_mn_ref <- {
        para_mn <- landscapemetrics::lsm_c_para_mn(reference, directions = 8)
        get_metric_value(para_mn)
      }

      lpi_ref <- {
        lpi <- landscapemetrics::lsm_c_lpi(reference, directions = 8)
        get_metric_value(lpi)
      }

      # Calculate differences
      calculate_change <- function(comp, ref) {
        if (ref == 0) {
          if (comp == 0)
            0
          else
            1
        } else {
          min(abs((comp - ref) / ref), 1)
        }
      }

      change_np <- calculate_change(np_comp, np_ref)
      change_enn_mn <- calculate_change(enn_mn_comp, enn_mn_ref)
      change_para_mn <- calculate_change(para_mn_comp, para_mn_ref)
      change_lpi <- calculate_change(lpi_comp, lpi_ref)

      # Overall configuration disagreement
      configuration_disagreement <- (change_np + change_enn_mn + change_para_mn + change_lpi) / 4
    } else {
      configuration_disagreement <- NA
    }

    # Calculate confusion matrix
    true_negative <- sum(config_qad_lite$ref_values == 0 &
                           config_qad_lite$comp_values == 0,
                         na.rm = TRUE)
    false_positive <- sum(config_qad_lite$ref_values == 0 &
                            config_qad_lite$comp_values == 1,
                          na.rm = TRUE)
    true_positive <- sum(config_qad_lite$ref_values == 1 &
                           config_qad_lite$comp_values == 1,
                         na.rm = TRUE)
    false_negative <- sum(config_qad_lite$ref_values == 1 &
                            config_qad_lite$comp_values == 0,
                          na.rm = TRUE)
    unknown_positive <- sum(config_qad_lite$ref_values == 2 &
                              config_qad_lite$comp_values == 1,
                            na.rm = TRUE)
    unknown_negative <- sum(config_qad_lite$ref_values == 2 &
                              config_qad_lite$comp_values == 0,
                            na.rm = TRUE)

    config_qad_lite$ref_values <- NULL
    config_qad_lite$comp_values <- NULL
    gc()

    # Calculate derived metrics
    total_obs <- true_negative + true_positive + false_negative + false_positive
    accuracy <- (true_negative + true_positive) / total_obs
    precision <- true_positive / (true_positive + false_positive)
    recall <- true_positive / (true_positive + false_negative)
    specificity <- true_negative / (true_negative + false_positive)

    # Calculate MCC value
    tp_fp <- ifelse(
      is.nan(true_positive + false_positive) ||
        (true_positive + false_positive) == 0,
      1,
      true_positive + false_positive
    )
    tp_fn <- ifelse(
      is.nan(true_positive + false_negative) ||
        (true_positive + false_negative) == 0,
      1,
      true_positive + false_negative
    )
    tn_fp <- ifelse(
      is.nan(true_negative + false_positive) ||
        (true_negative + false_positive) == 0,
      1,
      true_negative + false_positive
    )
    tn_fn <- ifelse(
      is.nan(true_negative + false_negative) ||
        (true_negative + false_negative) == 0,
      1,
      true_negative + false_negative
    )

    # Convert values to numeric to avoid integer overflow
    tp_fp <- as.numeric(tp_fp)
    tp_fn <- as.numeric(tp_fn)
    tn_fp <- as.numeric(tn_fp)
    tn_fn <- as.numeric(tn_fn)

    mcc <-
      ((true_positive * true_negative) - (false_positive * false_negative)) /
      sqrt(tp_fp * tp_fn * tn_fp * tn_fn)
    norm_mcc <- (mcc + 1) / 2

    # Replace NaN values for metrics with 0
    accuracy <- ifelse(is.nan(accuracy), 0, accuracy)
    precision <- ifelse(is.nan(precision), 0, precision)
    recall <- ifelse(is.nan(recall), 0, recall)
    specificity <- ifelse(is.nan(specificity), 0, specificity)

    # Calculate disagreements
    quantity_disagreement <- abs(
      config_qad_lite$positives_in_comparison - config_qad_lite$positives_in_reference
    )
    allocation_disagreement <- 2 * min(false_negative, false_positive)
    total_disagreement <- quantity_disagreement + allocation_disagreement

    # Calculate RMSE if applicable
    RMSE <- if (config_qad_lite$use_rmse) {
      Metrics::rmse(config_qad_lite$actual_predicted[, 1],
                    config_qad_lite$actual_predicted[, 2])
    } else {
      NA
    }

    config_qad_lite$actual_predicted <- NULL
    gc()

    # Calculate distance differences if applicable

    if (config_qad_lite$use_distance) {
      distance_difference <- sum(config_qad_lite$distance_differences)
      config_qad_lite$distance_differences <- NULL
      gc()
    } else {
      distance_difference <- NA
    }


    # Calculate odds ratio with adjustments to avoid NA or Inf
    odds_ratio <- if (false_negative == 0 && false_positive == 0) {
      (true_positive * true_negative) / 1
    } else if (false_negative == 0) {
      (true_positive * true_negative) / false_positive
    } else if (false_positive == 0) {
      (true_positive * true_negative) / false_negative
    } else {
      (true_positive * true_negative) / (false_negative * false_positive)
    }

    # Create the output data frame directly
    output <- data.frame(
      false_negatives = false_negative,
      false_positives = false_positive,
      true_positives = true_positive,
      true_negatives = true_negative,
      unknown_positives = unknown_positive,
      unknown_negatives = unknown_negative,
      accuracy = accuracy,
      precision = precision,
      recall = recall,
      specificity = specificity,
      quantity_disagreement = quantity_disagreement,
      allocation_disagreement = allocation_disagreement,
      total_disagreement = total_disagreement,
      configuration_disagreement = configuration_disagreement,
      odds_ratio = odds_ratio,
      residual_error = terra::global(abs(ref - comp), "sum", na.rm = TRUE)[[1]],
      true_infected_locations = config_qad_lite$positives_in_reference,
      simulated_infected_locations = config_qad_lite$positives_in_comparison,
      infected_locations_difference = config_qad_lite$positives_in_comparison -
        config_qad_lite$positives_in_reference,
      true_infecteds = config_qad_lite$positives_in_ref,
      simulated_infecteds = config_qad_lite$positives_in_comp,
      infecteds_difference = config_qad_lite$positives_in_comp -
        config_qad_lite$positives_in_ref,
      rmse = RMSE,
      distance_difference = distance_difference,
      mcc = mcc,
      norm_mcc = norm_mcc
    )
    dir.create(file.path(config_qad_lite$output_folder_path, "results"),
               showWarnings = FALSE, recursive = TRUE)

    results_fn <- file.path(config_qad_lite$output_folder_path, "results",
                            gsub("reclassed.tif", "results.csv",
                                 basename(config_qad_lite$comparison)))
    write.csv(output, results_fn, row.names = FALSE)
    return(cat("Completed & exported results:", basename(results_fn), "\n"))
  }
