#' Quantity Allocation Disagreement
#'
#' @param reference the binary (0 or 1) raster with ground truth data
#' @param comparison the binary (0 or 1) raster with simulated data to be compared to ground truth
#'
#' @return data frame of with disagreement and omission/comimission
#' @export
#'
#' @examples
quantity_allocation_disagreement <- function(reference, comparison){
  ## test that the comparison raster is the same extent, resolution, and crs as the reference (if not end)
  compareRaster(reference, comparison)
  compare <- reference - comparison
  compare2 <- as.matrix(compare)
  reference2 <- as.matrix(reference)
  comparison2 <- as.matrix(comparison)
  
  ## create data frame for comparison
  output <- data.frame(quantity_disagreement = 0, allocation_disagreement = 0, total_disagreement = 0 , omission = 0, commission = 0 ,number_of_infected_comp =0,total_allocation_distance = 0)
  output$total_disagreement <- sum(abs(compare2))
  output$quantity_disagreement <- abs(sum(reference2)-sum(comparison2))
  output$allocation_disagreement <- output$total_disagreement - output$quantity_disagreement
  output$omission <- abs(sum(compare2[compare2== -1]))
  output$commission <- abs(sum(compare2[compare2== 1]))
  output$number_of_infected_comp <- sum(comparison[comparison==1])
  
  return(output)
}
