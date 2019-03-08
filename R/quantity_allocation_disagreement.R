#' Compares quantity and allocation disagreement of two raster datasets
#' 
#' Uses quantity and allocation disagreement metrics by Pontius and Millones (2014) and 
#' omission and comission errors of comparing a modeled raster dataset to a reference raster 
#' datatset.
#'
#' @param reference the binary (0 or 1) raster with ground truth data
#' @param comparison the binary (0 or 1) raster with simulated data to be compared to ground truth
#'
#' @importFrom landscapemetrics lsm_c_np lsm_c_enn_mn lsm_c_para_mn lsm_c_lpi 
#' 
#' @return A data frame with spatial configuration metrics. Particularly quantity, allocation, and 
#' total disagreement,  omission and comission, and directional disagreement where directional disagreement.
#' 
#' @export
#'
#' @examples
#' library(raster)
#' reference <- raster(matrix(1, ncol = 2, nrow = 2))
#' reference[1,1] <- 0
#' comparison <- raster(matrix(1, ncol = 2, nrow = 2))
#' quantity_allocation_disagreement(reference, comparison)
#' 
quantity_allocation_disagreement <- function(reference, comparison){
  # test that the comparison raster is the same extent, resolution, and crs as the reference (if not end)
  raster::compareRaster(reference, comparison)
  compare <- reference - comparison
  # compare3 <- reference + comparison
  # extent = extent(reference)
  # num_of_cells = max(cellsFromExtent(reference, extent(reference)))
  
  # calculate number of infected patches
  reference_no0 <- reference
  comparison_no0 <- comparison
  reference_no0[reference_no0 == 0] <- NA
  comparison_no0[comparison_no0 == 0] <- NA
  NP_ref <- landscapemetrics::lsm_c_np(reference_no0, directions = 8)$value
  if (sum(comparison_no0[comparison_no0 >0]) == 0) {
    NP_comp <- 0
    ENN_MN_comp <- 0
    PARA_MN_comp <- 0
    LPI_comp <- 0
  } else {
    NP_comp <- landscapemetrics::lsm_c_np(comparison_no0, directions = 8)$value
  }
  
  change_NP <- abs((NP_comp - NP_ref)/NP_ref)
  
  # calculate the mean euclidean distance between patches
  if (NP_ref > 1) {
    ENN_MN_ref <- landscapemetrics::lsm_c_enn_mn(reference_no0, directions = 8, verbose = TRUE)$value
  } else  if (NP_ref == 1) {
    ENN_MN_ref <- 0
  }
  
  if (sum(comparison_no0[comparison_no0 > 0]) != 0 && NP_comp > 1) {
    ENN_MN_comp <- landscapemetrics::lsm_c_enn_mn(comparison_no0, directions = 8, verbose = TRUE)$value
  } else  if (sum(comparison_no0[comparison_no0 > 0]) != 0 && NP_comp <= 1) {
    ENN_MN_comp <- 0
  }
  
  if (ENN_MN_ref != 0) {
    change_ENN_MN <- abs((ENN_MN_comp - ENN_MN_ref)/ENN_MN_ref)
  } else if (ENN_MN_comp == 0 && ENN_MN_ref == 0) {
    change_ENN_MN <- 0
  } else {
    change_ENN_MN <- 1
  }

  # calculate the mean perimeter-area ratio of patches and the difference
  PARA_MN_ref <- landscapemetrics::lsm_c_para_mn(reference_no0, directions = 8)$value
  if (sum(comparison_no0[comparison_no0 >0]) == 0) {
    PARA_MN_comp <- 0
  } else if (sum(comparison_no0[comparison_no0 >0]) != 0) {
    PARA_MN_comp <- landscapemetrics::lsm_c_para_mn(comparison_no0, directions = 8)$value
  }
    
  change_PARA_MN <- abs((PARA_MN_comp - PARA_MN_ref)/PARA_MN_ref) 
  
  # calculate the largest patch index and difference
  LPI_ref <- landscapemetrics::lsm_c_lpi(reference_no0, directions = 8)$value
  if (sum(comparison_no0[comparison_no0 >0]) == 0) {
    LPI_comp <- 0
  } else if (sum(comparison_no0[comparison_no0 >0]) != 0) {
    LPI_comp <- landscapemetrics::lsm_c_lpi(comparison_no0, directions = 8)$value
  }
  
  change_LPI <- abs(LPI_comp - LPI_ref) / 100
  
  # calculate landscape similarity index between reference and comparison
  LSI <- 1 - ((change_NP + change_ENN_MN + change_PARA_MN + change_LPI) / 4)
  if (LSI < 0) { LSI <- 0 }
  
  ## create data frame for comparison
  output <- data.frame(quantity_disagreement = 0, allocation_disagreement = 0, total_disagreement = 0, omission = 0, commission = 0 , number_of_infected_comp = 0, directional_disagreement = 0, landscape_similarity = 0)
  output$total_disagreement <- sum(compare[compare == 1]) + abs(sum(compare[compare == -1]))
  output$quantity_disagreement <- abs(sum(compare[compare == 1]) + sum(compare[compare == -1]))
  output$allocation_disagreement <- output$total_disagreement - output$quantity_disagreement
  output$omission <- abs(sum(compare[compare == 1]))
  output$commission <- abs(sum(compare[compare == -1]))
  output$number_of_infected_comp <- sum(comparison[comparison == 1])
  output$directional_disagreement <- sum(compare[compare == 1]) + sum(compare[compare == -1])
  output$landscape_similarity <- LSI
  # output$true_positives <- abs(sum(compare3[compare3 ==2]))
  # output$true_negatives <- num_of_cells - output$omission - output$commission - output$true_positives
  # output$odds_ratio = (output$true_positives*output$true_negatives)/(output$omission*output$commission)
  
  return(output)
}
