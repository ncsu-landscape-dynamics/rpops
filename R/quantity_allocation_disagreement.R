#' Compares quantity and allocation disagreement of two raster datasets
#' 
#' Uses quantity and allocation disagreement metrics by Pontius and Millones (2014) and 
#' omission and comission errors of comparing a modeled raster dataset to a reference raster 
#' datatset.
#'
#' @param reference the binary (0 or 1) raster with ground truth data
#' @param comparison the binary (0 or 1) raster with simulated data to be compared to ground truth
#' @param configuration Set to true if you want to use configuration disagreement for comparing model runs. Default is FALSE.
#' @param mask Used to provide a mask to remove 0's that are not true negatives from comparisons. 
#'
#' @importFrom landscapemetrics lsm_c_np lsm_c_enn_mn lsm_c_para_mn lsm_c_lpi 
#' @importFrom raster cellsFromExtent xres ncol nrow yres extent
#' 
#' @return A data frame with spatial configuration metrics. Particularly quantity, allocation, and 
#' total disagreement,  omission and comission, and directional disagreement where directional disagreement.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' library(raster)
#' reference <- raster(matrix(1, ncol = 2, nrow = 2))
#' reference[1,1] <- 0
#' comparison <- raster(matrix(1, ncol = 2, nrow = 2))
#' quantity_allocation_disagreement(reference, comparison)
#' }
quantity_allocation_disagreement <- function(reference, comparison, configuration = FALSE, mask = NULL){
  if (is.null(mask)) {
    
  } else {
    reference[is.na(mask)] <- NA
    comparison[is.na(mask)] <- NA
  }
  ## test that the comparison raster is the same extent, resolution, and crs as the reference to ensure that they can be compared accurately
  raster::compareRaster(reference, comparison)
  # save initial reference and comparison to use for residual error calculation and then reclassify the original reference and comparison to binary values for other
  # classifications as we are only concerned with locational accuracy of infection not exact population predictions (residual error is a comparison of exact population numbers)
  comp <- comparison
  ref <- reference
  rcl <- c(1, Inf, 1, 0, 0.99, 0)
  rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
  reference <- reclassify(reference, rclmat)
  comparison <- reclassify(comparison, rclmat)
  
  if (configuration == TRUE) {
    # calculate number of infected patches
    NP_ref <- landscapemetrics::lsm_c_np(reference, directions = 8)$value[2]
    if (sum(comparison[comparison > 0]) == 0) {
      NP_comp <- 0
      ENN_MN_comp <- 0
      LPI_comp <- 0
      PARA_MN_comp <- 0
    } else {
      NP_comp <- landscapemetrics::lsm_c_np(comparison, directions = 8)$value[2]
    }

    change_NP <- abs((NP_comp - NP_ref)/(NP_ref))
    if (change_NP > 1) {change_NP <- 1}
    if (change_NP >= 1) {change_NP <- 1}

    # calculate the mean euclidean distance between patches
    if (NP_ref > 1) {
      ENN_MN_ref <- landscapemetrics::lsm_c_enn_mn(reference, directions = 8, verbose = TRUE)$value[2]
    } else  if (NP_ref == 1) {
      ENN_MN_ref <- 0
    }
    
    if (sum(comparison[comparison > 0]) != 0 && NP_comp > 1) {
      ENN_MN_comp <- landscapemetrics::lsm_c_enn_mn(comparison, directions = 8, verbose = TRUE)$value[2]
    } else  if (sum(comparison[comparison > 0]) != 0 && NP_comp <= 1) {
      ENN_MN_comp <- 0
    }
    
    if (ENN_MN_ref != 0) {
      change_ENN_MN <- abs((ENN_MN_comp - ENN_MN_ref)/(ENN_MN_ref))
      if (change_ENN_MN > 1) {change_ENN_MN <- 1}
    } else if (ENN_MN_comp == 0 && ENN_MN_ref == 0) {
      change_ENN_MN <- 0
    } else {
      change_ENN_MN <- 1
    }
    
    # calculate the mean perimeter-area ratio of patches and the difference
    PARA_MN_ref <- landscapemetrics::lsm_c_para_mn(reference, directions = 8)$value[2]
    if (sum(comparison[comparison >0]) == 0) {
      PARA_MN_comp <- 0
    } else if (sum(comparison[comparison >0]) != 0) {
      PARA_MN_comp <- landscapemetrics::lsm_c_para_mn(comparison, directions = 8)$value[2]
    }
    
    change_PARA_MN <- abs((PARA_MN_comp - PARA_MN_ref)/(PARA_MN_ref))
    if (change_PARA_MN > 1) {change_PARA_MN <- 1}
    
    # calculate the largest patch index and difference
    LPI_ref <- landscapemetrics::lsm_c_lpi(reference, directions = 8)$value[2]
    if (sum(comparison[comparison >0]) == 0) {
      LPI_comp <- 0
    } else if (sum(comparison[comparison >0]) != 0) {
      LPI_comp <- landscapemetrics::lsm_c_lpi(comparison, directions = 8)$value[2]
    }
    
    change_LPI <- abs((LPI_comp - LPI_ref)/(LPI_ref))
    if (change_LPI > 1) {change_LPI <- 1}
    configuration_disagreement <- ((change_NP + change_ENN_MN + change_PARA_MN + change_LPI) / 4)

  } else {
    configuration_disagreement <- 0
  }

  ## calculate reference and comparison totals to use for creation of probabilities 
  positives_in_reference <- sum(reference[!is.na(reference)] == 1)
  negatives_in_reference <- sum(reference[!is.na(reference)] == 0)
  total_in_reference <- sum(reference[!is.na(reference)] >= 0)
  positives_in_comparison <- sum(comparison[!is.na(reference)] == 1)
  negatives_in_comparison <- sum(comparison[!is.na(reference)] == 0)
  
  ## calculate confusion matrix for accurracy assessment
  true_positive <- sum(comparison[reference == 1] == 1)
  false_positive <- sum(comparison[reference == 0] == 1)
  false_negative <- sum(comparison[reference == 1] == 0)
  true_negative <- sum(comparison[reference == 0] == 0)
  
  ## calculate quantity and allocation disagreements for infected/infested from probabilities based on Death to Kappa (Pontius et al. 2011)
  # quantity_disagreement <- abs((probability_11 + probability_01) - (probability_11 + probability_10))
  # allocation_disagreement <- min((probability_11 + probability_10) - probability_11, (probability_11 + probability_01) - probability_11)
  quantity_disagreement <- abs(positives_in_comparison - positives_in_reference)
  allocation_disagreement <- 2 * min(false_negative, false_positive)
  total_disagreement <- quantity_disagreement + allocation_disagreement
  ## calculate odds ratio with adjustments so can never be NA or INF
  if(false_negative == 0 && false_positive == 0) {
    odds_ratio <- (true_positive * true_negative) / 1
  } else if (false_negative == 0) {
    odds_ratio <- (true_positive * true_negative) / false_positive
  } else if (false_positive == 0) {
    odds_ratio <- (true_positive * true_negative) / false_negative
  } else {
    odds_ratio <- (true_positive * true_negative) / (false_negative * false_positive)
  }
  ## create data frame for outputs and add calculated values to it
  output <- data.frame(quantity_disagreement = 0, allocation_disagreement = 0, total_disagreement = 0, configuration_disagreement = 0, omission = 0, commission = 0,  true_positives = 0, true_negatives = 0, odds_ratio = 0)
  output$omission <- false_negative
  output$commission <- false_positive
  output$true_positives <- true_positive
  output$true_negatives <- true_negative
  output$quantity_disagreement <- quantity_disagreement
  output$allocation_disagreement <- allocation_disagreement
  output$total_disagreement <- total_disagreement
  output$configuration_disagreement <- configuration_disagreement
  output$odds_ratio <- odds_ratio
  output$residual_error <- cellStats(abs(ref - comp), 'sum')
  output$true_infected <- positives_in_reference
  output$simulated_infected <- positives_in_comparison
  output$infected_difference <- positives_in_comparison - positives_in_reference
  
  return(output)
}