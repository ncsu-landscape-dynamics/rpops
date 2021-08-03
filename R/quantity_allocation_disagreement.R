#' Compares quantity and allocation disagreement of two raster datasets
#'
#' Uses quantity and allocation disagreement metrics by Pontius and Millones
#' (2014) and omission and comission errors of comparing a modeled raster
#' dataset to a reference raster datatset.
#'
#' @param reference the binary (0 or 1) raster with ground truth data
#' @param comparison the binary (0 or 1) raster with simulated data to be
#' compared to ground truth
#' @param configuration Set to true if you want to use configuration
#' disagreement for comparing model runs. Default is FALSE.
#' @param mask Used to provide a mask to remove 0's that are not true
#' negatives from comparisons.
#'
#' @importFrom landscapemetrics lsm_c_np lsm_c_enn_mn lsm_c_para_mn lsm_c_lpi
#' @importFrom terra cells xres ncol nrow yres ext compareGeom
#'
#' @return A data frame with spatial configuration metrics. Particularly
#' quantity, allocation, and total disagreement,  omission and comission, and
#' directional disagreement where directional disagreement.
#'
#' @export
#'
quantity_allocation_disagreement <-
  function(reference, comparison, configuration = FALSE, mask = NULL) {
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
    rcl_comp <- c(1, Inf, 1, 0, 0.99, 0)
    rclmat_comp <- matrix(rcl_comp, ncol = 3, byrow = TRUE)
    ## use 2 to indicate areas that aren't sampled in the reference data. This
    ## allows for the calculation of pure non-inflated accuracy statistics and
    ## to examine where the model is predicting
    rcl_ref <- c(NA, 0, 2, 1, Inf, 1, 0, 0.99, 0)
    rclmat_ref <- matrix(rcl_ref, ncol = 3, byrow = TRUE)
    reference <- terra::classify(reference, rclmat_ref)
    comparison <- terra::classify(comparison, rclmat_comp)

    if (configuration == TRUE) {
      # calculate number of infected patches
      np_ref <- landscapemetrics::lsm_c_np(reference, directions = 8)$value[2]
      if (sum(terra::values(comparison) > 0, na.rm = TRUE) == 0) {
        np_comp <- 0
        enn_mn_comp <- 0
        lpi_comp <- 0
        para_mn_comp <- 0
      } else {
        np_comp <-
          landscapemetrics::lsm_c_np(comparison, directions = 8)$value[2]
      }

      change_np <- abs((np_comp - np_ref) / (np_ref))
      if (change_np > 1) {
        change_np <- 1
      }
      if (change_np >= 1) {
        change_np <- 1
      }

      # calculate the mean euclidean distance between patches
      if (np_ref > 1) {
        enn_mn_ref <-
          landscapemetrics::lsm_c_enn_mn(reference,
            directions = 8,
            verbose = TRUE
          )$value[2]
      } else if (np_ref == 1) {
        enn_mn_ref <- 0
      }

      if (sum(terra::values(comparison) > 0, na.rm = TRUE) != 0 &&
          np_comp > 1) {
        enn_mn_comp <-
          landscapemetrics::lsm_c_enn_mn(comparison,
            directions = 8,
            verbose = TRUE
          )$value[2]
      } else if (sum(terra::values(comparison) > 0, na.rm = TRUE) != 0 &&
                 np_comp <= 1) {
        enn_mn_comp <- 0
      }

      if (enn_mn_ref != 0) {
        change_enn_mn <- abs((enn_mn_comp - enn_mn_ref) / (enn_mn_ref))
        if (change_enn_mn > 1) {
          change_enn_mn <- 1
        }
      } else if (enn_mn_comp == 0 && enn_mn_ref == 0) {
        change_enn_mn <- 0
      } else {
        change_enn_mn <- 1
      }

      # calculate the mean perimeter-area ratio of patches and the difference
      para_mn_ref <-
        landscapemetrics::lsm_c_para_mn(reference, directions = 8)$value[2]
      if (sum(terra::values(comparison) > 0, na.rm = TRUE) == 0) {
        para_mn_comp <- 0
      } else if (sum(terra::values(comparison) > 0, na.rm = TRUE) != 0) {
        para_mn_comp <-
          landscapemetrics::lsm_c_para_mn(comparison, directions = 8)$value[2]
      }

      change_para_mn <- abs((para_mn_comp - para_mn_ref) / (para_mn_ref))
      if (change_para_mn > 1) {
        change_para_mn <- 1
      }

      # calculate the largest patch index and difference
      lpi_ref <- landscapemetrics::lsm_c_lpi(reference, directions = 8)$value[2]
      if (sum(terra::values(comparison) > 0, na.rm = TRUE) == 0) {
        lpi_comp <- 0
      } else if (sum(terra::values(comparison) > 0, na.rm = TRUE) != 0) {
        lpi_comp <-
          landscapemetrics::lsm_c_lpi(comparison, directions = 8)$value[2]
      }

      change_lpi <- abs((lpi_comp - lpi_ref) / (lpi_ref))
      if (change_lpi > 1) {
        change_lpi <- 1
      }
      configuration_disagreement <-
        ((change_np + change_enn_mn + change_para_mn + change_lpi) / 4)
    } else {
      configuration_disagreement <- 0
    }

    # calculate reference and comparison totals to use for creation of
    # probabilities
    positives_in_reference <- sum(terra::values(reference) == 1, na.rm = TRUE)
    positives_in_comparison <- sum(terra::values(comparison) == 1, na.rm = TRUE)

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
    specificiity <- true_negative / (true_negative + false_positive)

    # calculate quantity and allocation disagreements for infected/infested from
    # probabilities based on Death to Kappa (Pontius et al. 2011)
    quantity_disagreement <-
      abs(positives_in_comparison - positives_in_reference)
    allocation_disagreement <- 2 * min(false_negative, false_positive)
    total_disagreement <- quantity_disagreement + allocation_disagreement
    # calculate odds ratio with adjustments so can never be NA or INF
    if (false_negative == 0 && false_positive == 0) {
      odds_ratio <- (true_positive * true_negative) / 1
    } else if (false_negative == 0) {
      odds_ratio <- (true_positive * true_negative) / false_positive
    } else if (false_positive == 0) {
      odds_ratio <- (true_positive * true_negative) / false_negative
    } else {
      odds_ratio <-
        (true_positive * true_negative) / (false_negative * false_positive)
    }
    # create data frame for outputs and add calculated values to it
    output <-
      data.frame(
        quantity_disagreement = 0,
        allocation_disagreement = 0,
        total_disagreement = 0,
        configuration_disagreement = 0,
        false_negative = 0,
        false_positive = 0,
        true_positives = 0,
        true_negatives = 0,
        odds_ratio = 0
      )
    output$false_negative <- false_negative
    output$false_positive <- false_positive
    output$true_positives <- true_positive
    output$true_negatives <- true_negative
    output$unknown_positive <- unknown_positive
    output$unknown_negative <- unknown_negative
    output$accuracy <- accuracy
    output$precision <- precision
    output$recall <- recall
    output$specificiity <- specificiity
    output$quantity_disagreement <- quantity_disagreement
    output$allocation_disagreement <- allocation_disagreement
    output$total_disagreement <- total_disagreement
    output$configuration_disagreement <- configuration_disagreement
    output$odds_ratio <- odds_ratio
    output$residual_error <-
      terra::global(abs(ref - comp), "sum", na.rm = TRUE)[[1]]
    output$true_infected <- positives_in_reference
    output$simulated_infected <- positives_in_comparison
    output$infected_difference <-
      positives_in_comparison - positives_in_reference

    return(output)
  }
