#' @title PoPS (Pest or Pathogen Spread) model
#'
#' @description A dynamic species distribution model for pest or pathogen spread
#' in forest or agricultural ecosystems. The model is process based meaning that
#' it uses understanding of the effect of weather and other environmental
#' factors on reproduction and survival of the pest/pathogen in order to
#' forecast spread of the pest/pathogen into the future. This function performs
#' a single stochastic realization of the model and is predominantly used for
#' automated tests of model features.
#' @param config Path to config file produced when calling `configuration`.
#' The config file includes all data necessary used to set up c++ PoPS model
#'
#' @useDynLib PoPS, .registration = TRUE
#' @importFrom terra app rast xres yres classify extract ext as.points ncol nrow project
#' nlyr rowFromCell colFromCell values as.matrix rowFromCell colFromCell crs
#' @importFrom Rcpp sourceCpp evalCpp
#' @importFrom stats runif
#' @importFrom lubridate interval time_length mdy %within%
#' @importFrom utils read.csv read.table
#' @importFrom methods is
#' @return list of infected and susceptible per year
#' @export
#'

pops <- function(config) {

  if (!is.null(config$failure)) {
    stop(config$failure)
  }
  set.seed(config$random_seed_list[[1]])
  config <- draw_parameters(config) # draws parameter set for the run
  config <- host_pool_setup(config)
  while (any(config$total_hosts > config$total_populations, na.rm = TRUE) ||
         any(config$total_exposed > config$total_populations, na.rm = TRUE) ||
         any(config$total_infecteds > config$total_populations, na.rm = TRUE)) {
    config <- host_pool_setup(config)
  }
  config$competency_table_list <- competency_table_list_creator(config$competency_table)
  config$pest_host_table_list <- pest_host_table_list_creator(config$pest_host_table)
  config$random_seed <- config$random_seed_list[[1]]
  config$random_seeds <- as.matrix(config$random_seeds_list[1, ])[1, ]

  data <- PoPS::pops_model(config)

  return(data)
}
