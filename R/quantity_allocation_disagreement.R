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
  ## test that the comparison raster is the same extent, resolution, and crs as the reference (if not end)
  raster::compareRaster(reference, comparison)
  compare <- reference - comparison
  compare2 <- raster::as.matrix(compare)
  reference2 <- raster::as.matrix(reference)
  comparison2 <- raster::as.matrix(comparison)
  
  #calculate number of infected patches
  reference_no0 <- reference
  comparison_no0 <- comparison
  reference_no0[reference_no0==0]<-NA
  comparison_no0[comparison_no0==0]<-NA
  NP_ref_full <- landscapemetrics::lsm_c_np(reference_no0, directions=8)
  NP_comp_full <- landscapemetrics::lsm_c_np(comparison_no0, directions=8)
  NP_ref <- NP_ref_full$value
  NP_comp <- NP_comp_full$value
  change_NP <- abs((NP_comp - NP_ref)/NP_ref)
  
  #calculate the mean euclidean distance between patches
  if (NP_comp > 1) {
    ENN_MN_ref_full <- landscapemetrics::lsm_c_enn_mn(reference_no0, directions=8, verbose=TRUE)
    ENN_MN_comp_full <- landscapemetrics::lsm_c_enn_mn(comparison_no0, directions=8, verbose=TRUE)
    ENN_MN_ref <- ENN_MN_ref_full$value
    ENN_MN_comp <- ENN_MN_comp_full$value
    change_ENN_MN <- abs((ENN_MN_comp - ENN_MN_ref)/ENN_MN_ref)
  } else  if (NP_comp == 1) {
    ENN_MN_comp <- NA
    change_ENN_MN <- NA
  }
  
  # calculate the mean perimeter-area ratio of patches and the difference
  PARA_MN_ref_full <- landscapemetrics::lsm_c_para_mn(reference_no0, directions=8)
  PARA_MN_comp_full <- landscapemetrics::lsm_c_para_mn(comparison_no0, directions=8)
  PARA_MN_ref <- PARA_MN_ref_full$value
  PARA_MN_comp <- PARA_MN_comp_full$value
  change_PARA_MN <- abs((PARA_MN_comp - PARA_MN_ref)/PARA_MN_ref) 
  
  # calculate the largest patch index and difference
  LPI_ref_full <- landscapemetrics::lsm_c_lpi(reference_no0, directions=8)
  LPI_comp_full <- landscapemetrics::lsm_c_lpi(comparison_no0, directions=8)
  LPI_ref <- LPI_ref_full$value
  LPI_comp <- LPI_comp_full$value
  change_LPI <- abs(LPI_comp - LPI_ref) / 100
  
  # calculate landscape similarity index between reference and comparison
  LSI <- 1 - ((change_NP + change_ENN_MN + change_PARA_MN + change_LPI)/8)
  
  ## create data frame for comparison
  output <- data.frame(quantity_disagreement = 0, allocation_disagreement = 0, total_disagreement = 0 , omission = 0, commission = 0 ,number_of_infected_comp =0,directional_disagreement = 0, landscape_similarity=0)
  output$total_disagreement <- sum(abs(compare2))
  output$quantity_disagreement <- abs(sum(reference2)-sum(comparison2))
  output$allocation_disagreement <- output$total_disagreement - output$quantity_disagreement
  output$omission <- abs(sum(compare2[compare2== 1]))
  output$commission <- abs(sum(compare2[compare2== -1]))
  output$number_of_infected_comp <- sum(comparison[comparison==1])
  output$directional_disagreement <- sum(compare2)
  output$landscape_similarity <- LSI
  
  return(output)
}
