
# add error variables that can be returned from functions and used in tests
file_exists_error <- "does not exist."

detailed_file_exists_error <- function(file_name) {
  error_message <-
    paste(deparse(substitute(file_name)), "does not exist. Current path is", file_name, sep = " ")
  return(error_message)
}
frequency_error <-
  "Output frequency must be either 'week', 'month', 'day', 'year', 'time_step', or 'every_n_steps'"
raster_type_error <-
  "file is not one of '.grd', '.tif', '.img', or '.vrt'"
extent_error <-
  "Extents of input rasters do not match. Ensure that all of your input rasters have the same
extent"
resolution_error <-
  "Resolution of input rasters do not match. Ensure that all of your input rasters have the same
resolution"
crs_error <-
  "Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input
rasters have the same crs"
treatment_length_error <-
  "Length of list for treatment dates and treatments_file must be equal"
pesticide_length_error <-
  "Length of list for treatment dates and pesticide_duration must be equal"
treatment_method_error <-
  "treatment method is not one of the valid treatment options"
time_step_error <-
  "Time step must be one of 'week', 'month' or 'day'"
output_frequency_error <-
  "Output frequency is more frequent than time_step. The minimum output_frequency you can use is
  the time_step of your simulation. You can set the output_frequency to 'time_step' to default to
  the most frequent output possible"
date_format_error <-
  "End time and/or start time not of type numeric and/or in format YYYY-MM-DD"
output_type_error <-
  "Output frequency must be either 'week', 'month', 'day', 'year', 'time_step', or 'every_n_steps'"
prior_means_error <-
  "There are not enough prior_means to communte the posterior means for all paramaters"
prior_cov_matrix_error <-
  "The prior covariance matrix is not the correct dimension to match the calibrated covariance
  matrix in order to compute the posterior covariance matrix"
# multi-species errors based on matching lengths of parameters
species_length_error <-
  "Length of list for species and infected_files must be equal"
parameter_length_error <-
  "Length of list for parameter_means and infected_files must be equal"
cov_matrix_length_error <-
  "Length of list for parameter_cov_matrix and infected_files must be equal"
model_type_length_error <-
  "Length of list for model_type and infected_files must be equal"
natural_kernel_length_error <-
  "Length of list for infected_files and natural_kernel_type must be equal"
anthropogenic_kernel_length_error <-
  "Length of list for infected_files and anthropogenic_kernel_type must be equal"
natural_dir_length_error <-
  "Length of list for infected_files and natural_dir must be equal"
anthropogenic_dir_length_error <-
  "Length of list for infected_files and anthropogenic_dir must be equal"
host_file_length_error <-
  "Length of list for infected_files and host_file must be equal"
total_population_length_error <-
  "Length of list for infected_files and total_populations_file must be equal"
temp_length_error <-
  "Length of list for infected_files and temp must be equal"
temperature_coefficient_length_error <-
  "Length of list for infected_files and temperature_coefficient_file must be equal"
precip_length_error <-
  "Length of list for infected_files and precip must be equal"
precipitation_coefficient_length_error <-
  "Length of list for infected_files and precipitation_coefficient_file must be equal"
latency_period_length_error <-
  "Length of list for infected_files and latency_period must be equal"
time_step_length_error <-
  "Length of list for infected_files and time_step must be equal"
season_month_start_length_error <-
  "Length of list for infected_files and season_month_start must be equal"
season_month_end_length_error <-
  "Length of list for infected_files and season_month_end must be equal"
use_lethal_length_error <-
  "Length of list for infected_files and use_lethal_temperature must be equal"
temperature_file_length_error <-
  "Length of list for infected_files and temperature_file must be equal"
lethal_temperature_length_error <-
  "Length of list for infected_files and lethal_temperature must be equal"
lethal_temperature_month_length_error <-
  "Length of list for infected_files and lethal_temperature_month must be equal"
mortality_on_length_error <-
  "Length of list for infected_files and mortality_on must be equal"
mortality_rate_length_error <-
  "Length of list for infected_files and mortality_rate must be equal"
mortality_time_lag_length_error <-
  "Length of list for infected_files and mortality_time_lag must be equal"
movements_file_length_error <-
  "Length of list for infected_files and movements_file must be equal"
use_movements_length_error <-
  "Length of list for infected_files and use_movements must be equal"
start_exposed_length_error <-
  "Length of list for infected_files and start_exposed must be equal"
quarantine_areas_length_error <-
  "Length of list for infected_files and quarantine_areas_file must be equal"
use_quarantine_length_error <-
  "Length of list for infected_files and use_quarantine must be equal"
use_spreadrates_length_error <-
  "Length of list for infected_files and use_spreadrates must be equal"
file_type_error <-
  "file is not one of '.csv' or '.txt'"
natural_kernel_error <-
  "Natural kernel type not one of 'cauchy', 'exponential', 'uniform', 'deterministic neighbor',
  'power law', 'hyperbolic secant', 'gamma', 'weibull', 'logistic'"
anthropogenic_kernel_error <-
  "Anthropogenic kernel type not one of 'cauchy', 'exponential', 'uniform', 'deterministic
  neighbor', 'power law', 'hyperbolic secant', 'gamma', 'weibull', 'logistic', 'network'"
covariance_mat_error <-
  "parameter covariance matrix is not 6 x 6"
paramter_means_error <-
  "parameter means is not a vector of length 6"
write_outputs_error <-
  "write_outputs is not one of c('all simulations', 'summary_outputs', 'None')"
output_path_error <-
  "output path doesn't exist"
model_type_error <-
  "Model type is not a valid type options are 'SI' or 'SEI'"
season_month_error <-
  "Season month start or end not between 1 and 12"
latency_period_error <-
  "Model type is set to SEI but the latency period is less than 1"
treatment_option_error <-
  "treatment method is not one of the valid treatment options"
network_min_distance_small_error <-
  "network min distance is less than half the cell resolution"
network_min_distance_large_error <-
  "network min distance is greater than the network max distance"
network_max_distance_large_error <-
  "network max distance is greater than the resoultion times the minimum NS or EW extent"
network_movement_error <- "network movement is not of type 'walk', 'jump', or 'teleport'"

infection_years_length_error <- function(num_layers_infected_years, number_of_time_steps) {
  error_message <-
    paste("The infection years file must have enough layers to match the number of outputs from the
    model. The number of layers of your infected year file is", num_layers_infected_years,
          "and the number of outputs is", number_of_time_steps, sep = " ")
  return(error_message)
}

success_metric_error <- "success_metric is not one of the listed options."

initial_cond_uncert_error <-
  "use_initial_condition_uncertainty is TRUE but the number of layers in the infected file is not 2.
  This should be a raster file with 2 layers the first being the mean value and the second the
  stadard deviation"

host_uncert_error <-
  "use_host_uncertainty is TRUE but the number of layers in the host file is not 2. This should be
  a raster file with 2 layers the first being the mean value and the second the stadard deviation."

random_seeds_dimensions_error <-
  "Either number of rows and columns in random_seeds_file does not equal the number_of_iterations
  set in the model or the number of columns does not equal the number of unique random seeds"

weather_type_error <- "Weather type is not one of 'probabilistic', 'deterministic', or 'none'"

weather_size_deterministic_error <-
  "Weather coeeficient number of layers with deterministic is not equal to the total number of time
  steps."

weather_size_probabilitic_error <-
  "Weather coefficient number of layers with probablisitc is not equal to the total number of time
  steps annual."

weather_sd_layer_error <-
  "weather coefficient sd file number of layers not equal to number of layers in weather coefficient
  file"

multihost_file_length_error <-
  "infected_file_list and host_file_list are not of the same length, ensure both of these files are
  of the length of the number of host species you want to simulate"

competency_table_column_length_error <-
  "competency_table doesn't have the same number of columns as number of files in host_file_list"

competency_table_row_length_error <-
  "competency_table needs to have at least 1 more row than the number of hosts being modeled which
  is represented by the number of file in the host_file_list"

competency_table_wrong_columns <-
  "Check column order and headings. The competency table requires a column for each
  host species, followed by a competency_mean column and competency_sd column"

competency_value_error <-
  "competency_table competency_mean and competency_sd values must be between 0 and 1"

pest_host_table_wrong_columns <-
  "pest_host_table must the 6 columns named and order: host, susceptibility_mean,
  susceptibility_sd, mortality_rate_mean, mortality_rate_sd, mortality_time_lag"

pest_host_susceptbility_value_error <-
  "pest_host_table susceptiblity_mean and susceptibility_sd values must be between 0 and 1"

pest_host_mortality_rate_value_error <-
  "pest_host_table mortality_rate_mean and mortality_rate_sd values must be between 0 and 1"

pest_host_table_row_length_error <-
  "pest_host_table doesn't have the same number of rows as number of files in host_file_list"

multihosts_gt_totpop_error <-
  "All hosts sum to more than the total populations in some cells. Check rasters to ensure that
  combined summed host layers are not greater than total populations raster."

multiinfected_gt_totpop_error <-
  "All infecteds sum to more than the total populations in some cells. Check rasters to ensure that
  combined summed infected layers are not greater than total populations raster."

multiexposed_gt_totpop_error <-
  "All exposeds sum to more than the total populations in some cells. Check rasters to ensure that
  combined summed exposed layers are not greater than total populations raster."

crs_infected_county_error <-
  "Coordinate reference system (crs) of input infected vector does not match. Ensure that all of
  your input rasters and vectors have the same crs"
