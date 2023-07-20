#' Compares quantity and allocation disagreement of two raster data sets
#'
#' Uses quantity and allocation disagreement metrics by Pontius and Millones
#' (2014) and omission and commission errors of comparing a modeled raster
#' data set to a reference raster data set.
#'
#' @param reference the raster with ground truth data. For all metrics expect
#' RMSE these values are reclassified with values > 1 becoming 1, values < 1
#' going to 0, and NA values staying NA.
#' @param comparison the raster with simulated data. For all metrics expect
#' RMSE these values are reclassified with values > 1 becoming 1, values < 1
#' going to 0.
#' @param use_configuration Boolean if you want to use configuration
#' disagreement for comparing model runs. Default is FALSE.
#' @param mask Used to provide a mask to remove 0's that are not true
#' negatives from comparisons.
#' @param use_distance Boolean if you want to compare distance between
#' simulations and observations. Default is FALSE.
#'
#' @importFrom landscapemetrics lsm_c_np lsm_c_enn_mn lsm_c_para_mn lsm_c_lpi
#' @importFrom terra cells xres ncol nrow yres ext compareGeom vect
#' @importFrom Metrics rmse
#'
#' @return A data frame with spatial configuration metrics. Particularly
#' quantity, allocation, and total disagreement,  omission and commission, and
#' directional disagreement where directional disagreement.
#'
#' @export
#'
quantity_allocation_disagreement <-
  function(reference, comparison, use_configuration = FALSE,
           mask = NULL, use_distance = FALSE) {
    if (!is.null(mask)) {
      reference <- terra::mask(reference, mask)
      comparison <- terra::mask(comparison, mask)
    }
    # test that the comparison raster is the same extent, resolution, and crs as
    # the reference to ensure that they can be compared accurately
    # save initial reference and comparison to use for residual error
    # calculation and then reclassify the original reference and comparison to
    # binary values for other classifications as we are only concerned with
    # locational accuracy of infection not exact population predictions
    # (residual error is a comparison of exact population numbers)
    comp <- comparison
    ref <- reference
    rcl_comp <- c(1, Inf, 1, 0, 1, 0, NA, 0, 0)
    rclmat_comp <- matrix(rcl_comp, ncol = 3, byrow = TRUE)
    ## use 2 to indicate areas that aren't sampled in the reference data. This
    ## allows for the calculation of pure non-inflated accuracy statistics and
    ## to examine where the model is predicting
    rcl_ref <- c(NA, 0, 2, 1, Inf, 1, 0, 0.99, 0)
    rclmat_ref <- matrix(rcl_ref, ncol = 3, byrow = TRUE)
    reference <- terra::classify(reference, rclmat_ref, right = FALSE)
    comparison <- terra::classify(comparison, rclmat_comp, right = FALSE)

    if (use_configuration) {
      if (sum(terra::values(comparison) > 0, na.rm = TRUE) == 0) {
        np_comp <- 0
        enn_mn_comp <- 0
        lpi_comp <- 0
        para_mn_comp <- 0
      } else {
        # calculate number of infected patches
        np_comps <- landscapemetrics::lsm_c_np(comparison, directions = 8)
        if (any(unique(np_comps$class) %in% 1)) {
          np_comp <- np_comps$value[np_comps$class == 1]
        } else {
          np_comp <- 0
        }

        # calculate the mean euclidean distance between patches
        if (np_comp > 1) {
          enn_mn_comps <- landscapemetrics::lsm_c_enn_mn(comparison, directions = 8, verbose = TRUE)
          if (any(unique(enn_mn_comps$class) %in% 1)) {
            enn_mn_comp <- enn_mn_comps$value[enn_mn_comps$class == 1]
          } else {
            enn_mn_comp <- 0
          }
        } else if (np_comp <= 1) {
          enn_mn_comp <- 0
        }

        # calculate the mean perimeter-area ratio of patches and the difference
        para_mn_comps <- landscapemetrics::lsm_c_para_mn(comparison, directions = 8)
        if (any(unique(para_mn_comps$class) %in% 1)) {
          para_mn_comp <- para_mn_comps$value[para_mn_comps$class == 1]
        } else {
          para_mn_comp <- 0
        }

        # calculate the largest patch index and difference
        lpi_comps <- landscapemetrics::lsm_c_lpi(comparison, directions = 8)
        if (any(unique(lpi_comps$class) %in% 1)) {
          lpi_comp <- lpi_comps$value[lpi_comps$class == 1]
        } else {
          lpi_comp <- 0
        }
      }
      # calculate number of infected patches in reference (observed) and comparison (simulated)
      # data
      np_refs <- landscapemetrics::lsm_c_np(reference, directions = 8)
      if (any(unique(np_refs$class) %in% 1)) {
        np_ref <- np_refs$value[np_refs$class == 1]
      } else {
        np_ref <- 0
      }
      change_np <- abs((np_comp - np_ref) / (np_ref))
      if (change_np >= 1) {
        change_np <- 1
      }

      # calculate the mean euclidean distance between patches
      if (np_ref > 1) {
        enn_mn_refs <- landscapemetrics::lsm_c_enn_mn(reference, directions = 8, verbose = TRUE)
        if (any(unique(enn_mn_refs$class) %in% 1)) {
          enn_mn_ref <- enn_mn_refs$value[enn_mn_refs$class == 1]
        } else {
          enn_mn_ref <- 0
        }
      } else if (np_ref == 1) {
        enn_mn_ref <- 0
      }

      if (enn_mn_ref != 0) {
        change_enn_mn <- abs((enn_mn_comp - enn_mn_ref) / (enn_mn_ref))
        if (change_enn_mn >= 1) {
          change_enn_mn <- 1
        }
      } else if (enn_mn_comp == 0 && enn_mn_ref == 0) {
        change_enn_mn <- 0
      } else {
        change_enn_mn <- 1
      }

      # calculate the mean perimeter-area ratio of patches and the difference
      para_mn_refs <- landscapemetrics::lsm_c_para_mn(reference, directions = 8)
      if (any(unique(para_mn_refs$class) %in% 1)) {
        para_mn_ref <- para_mn_refs$value[para_mn_refs$class == 1]
      } else {
        para_mn_ref <- 0
      }

      change_para_mn <- abs((para_mn_comp - para_mn_ref) / (para_mn_ref))
      if (change_para_mn > 1) {
        change_para_mn <- 1
      }

      # calculate the largest patch index and difference
      lpi_refs <- landscapemetrics::lsm_c_lpi(reference, directions = 8)
      if (any(unique(lpi_refs$class) %in% 1)) {
        lpi_ref <- lpi_refs$value[lpi_refs$class == 1]
      } else {
        lpi_ref <- 0
      }

      change_lpi <- abs((lpi_comp - lpi_ref) / (lpi_ref))
      if (change_lpi >= 1) {
        change_lpi <- 1
      }

      configuration_disagreement <- ((change_np + change_enn_mn + change_para_mn + change_lpi) / 4)
    } else {
      configuration_disagreement <- 0
    }

    # calculate reference and comparison totals to use for creation of
    # probabilities
    positives_in_reference <- sum(terra::values(reference) == 1, na.rm = TRUE)
    positives_in_comparison <- sum(terra::values(comparison) == 1, na.rm = TRUE)

    positives_in_ref <- sum(terra::values(ref), na.rm = TRUE)
    positives_in_comp <- sum(terra::values(comp), na.rm = TRUE)

    ## calculate confusion matrix for accuracy assessment
    true_negative <- sum(terra::values(reference) == 0 &
                           terra::values(comparison) == 0, na.rm = TRUE)
    false_positive <- sum(terra::values(reference) == 0 &
                            terra::values(comparison) == 1, na.rm = TRUE)
    true_positive <- sum(terra::values(reference) == 1 &
                           terra::values(comparison) == 1, na.rm = TRUE)
    false_negative <- sum(terra::values(reference) == 1 &
                            terra::values(comparison) == 0, na.rm = TRUE)
    unknown_positive <- sum(terra::values(reference) == 2 &
                              terra::values(comparison) == 1, na.rm = TRUE)
    unknown_negative <- sum(terra::values(reference) == 2 &
                              terra::values(comparison) == 0, na.rm = TRUE)
    total_obs <- true_negative + true_positive + false_negative + false_positive
    accuracy <- (true_negative + true_positive) / total_obs
    precision <- true_positive / (true_positive + false_positive)
    recall <- true_positive / (true_positive + false_negative)
    specificity <- true_negative / (true_negative + false_positive)
    
    ## calculate MCC value
    tp_fp <- as.double((true_positive + false_positive))
    tp_fn <- as.double((true_positive + false_negative))
    tn_fp <- as.double((true_negative + false_positive))
    tn_fn <- as.double((true_negative + false_negative))
    
    if (is.nan(tp_fp) || tp_fp == 0) {tp_fp <- 1}
    if (is.nan(tp_fn) || tp_fn == 0) {tp_fn <- 1}
    if (is.nan(tn_fp) || tn_fp == 0) {tn_fp <- 1}
    if (is.nan(tn_fn) || tn_fn == 0) {tn_fn <- 1}
    
    mcc <- ((true_positive * true_negative) - (false_positive * false_negative)) /
      sqrt(tp_fp * tp_fn * tn_fp * tn_fn)

    if (is.nan(accuracy)) {accuracy <- 0}
    if (is.nan(precision)) {precision <- 0}
    if (is.nan(recall)) {recall <- 0}
    if (is.nan(specificity)) {specificity <- 0}

    # calculate quantity and allocation disagreements for infected/infested from
    # probabilities based on Death to Kappa (Pontius et al. 2011)
    quantity_disagreement <- abs(positives_in_comparison - positives_in_reference)
    allocation_disagreement <- 2 * min(false_negative, false_positive)
    total_disagreement <- quantity_disagreement + allocation_disagreement

    # calculate RMSE for comparison (only accounts for areas that are infected or predicted to
    # be infected)
    totals <- ref + comp
    obs_points <- terra::as.points(totals)
    names(obs_points) <- "data"
    obs_points <- obs_points[obs_points$data > 0]
    actual <- terra::extract(ref, obs_points)
    predicted <- terra::extract(comp, obs_points)
    RMSE <- Metrics::rmse(actual[, 2], predicted[, 2])
    distance_difference <- 0

    if (use_distance) {
      sim_points <- terra::as.points(comp)
      names(sim_points) <- "data"
      sim_points <- sim_points[sim_points$data > 0]
      ref_points <- terra::as.points(ref)
      names(ref_points) <- "data"
      ref_points <- ref_points[ref_points$data > 0]
      dist <- terra::distance(ref_points, sim_points)
      if (is(dist, "matrix")) {
        distance_differences <- apply(dist, 2, min)
      }

      all_distances <- function(distance_differences) {
        distance_differences <- round(sqrt(sum(distance_differences^2)), digits = 0)
        return(distance_differences)
      }
      distance_differences <- lapply(distance_differences, all_distances)
      distance_differences <- unlist(distance_differences, recursive = TRUE, use.names = TRUE)

      distance_difference <- sum(distance_differences)
    }

    # calculate odds ratio with adjustments so can never be NA or INF
    if (false_negative == 0 && false_positive == 0) {
      odds_ratio <- (true_positive * true_negative) / 1
    } else if (false_negative == 0) {
      odds_ratio <- (true_positive * true_negative) / false_positive
    } else if (false_positive == 0) {
      odds_ratio <- (true_positive * true_negative) / false_negative
    } else {
      odds_ratio <- (true_positive * true_negative) / (false_negative * false_positive)
    }
    # create data frame for outputs and add calculated values to it
    output <-
      data.frame(
        quantity_disagreement = 0,
        allocation_disagreement = 0,
        total_disagreement = 0,
        configuration_disagreement = 0,
        false_negatives = 0,
        false_positives = 0,
        true_positives = 0,
        true_negatives = 0,
        odds_ratio = 0
      )
    output$false_negatives <- false_negative
    output$false_positives <- false_positive
    output$true_positives <- true_positive
    output$true_negatives <- true_negative
    output$unknown_positives <- unknown_positive
    output$unknown_negatives <- unknown_negative
    output$accuracy <- accuracy
    output$precision <- precision
    output$recall <- recall
    output$specificity <- specificity
    output$quantity_disagreement <- quantity_disagreement
    output$allocation_disagreement <- allocation_disagreement
    output$total_disagreement <- total_disagreement
    output$configuration_disagreement <- configuration_disagreement
    output$odds_ratio <- odds_ratio
    output$residual_error <- terra::global(abs(ref - comp), "sum", na.rm = TRUE)[[1]]
    output$true_infected_locations <- positives_in_reference
    output$simulated_infected_locations <- positives_in_comparison
    output$infected_locations_difference <- positives_in_comparison - positives_in_reference
    output$true_infecteds <- positives_in_ref
    output$simulated_infecteds <- positives_in_comp
    output$infecteds_difference <- positives_in_comp - positives_in_ref
    output$rmse <- RMSE
    output$distance_difference <- distance_difference
    output$mcc <- mcc

    return(output)
  }
