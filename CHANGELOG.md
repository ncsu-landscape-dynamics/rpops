# Change Log

All notable changes to this project should be documented in this file.

The _Unreleased_ section should become the release once the release is ready
and the text can be used as part of the release description.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
with an additional rule that interfaces to the core library follow the same numbering
although they are in separate repositories, so increasing version number there
increases version here although there is no specific tag or release for that in
this repository.

## [Unreleased]

### Added

- Final step as an option for output frequency in all functions (@petrasovaa, #163).

- Added ability to use mean and sd for temp and precip coefficient. Uses 3 new parameters: weather_type, temperature_coefficient_sd_file, and precipitation_coefficient_sd_file. This is added to help understand uncertainty in model predictions due to weather drivers (@ChrisJones687, #168).

- Added the ability to use different seeds for all processes within pops-core. Adds 2 new parameters: multiple_random_seeds (boolean for using this functionality) and random_seeds (set to NULL to allow internal model selection of kernels or pass a CSV with the number of rows being the number of model runs and columns being the kernels in order (@ChrisJones687, #168 and @petrasovaa).

- Added ability for model to simulate pathogen/pest survival in soil and then emergence/rain splash in an area (@ChrisJones687, #178).

- Model can now write out forecast ensemble members (@ChrisJones687, #179)

- Added ability for model to use multi-host pops-core API. The model now takes 2 new parameter tables for hosts. It also takes a list of rasters for each host allowing testing of hypothesis related to different host combinations being more conducive to spread on the landscape (@petrasovaa, @ChrisJones687, @wenzeslaus, #186 and #188).

- Added the ability to calibrate, validate, and simulate from county level infection data (@ChrisJones687, #196 and #198).

- Added in PoPS lite features which make the model faster and less memory intensive (@cyborginhas and @ChrisJones687, #207).

- Added in ability to use multiple networks as part of the anthropogenic kernel (@ChrisJones687, #226).
    
### Changed

- Quarantine directions can now be calculated only for directions of interest (@petrasovaa, #168).

- In PoPS Core we are more explicit with types and round to the nearest int for model aspects that require it (@wenzeslaus, #203, #204, #205, #206, and #210 and @ChrisJones687, #208 and #209).

### Fixed

- Fixed MCC calculations that returned NaN with very large datasets (@ChrisJones687, #161).

- Fixed quarantine directions outputs in pops-multirun (@petrasovaa, #162).

- Fixed probabilistic weather not having the right number of time steps to check against (@ChrisJones687, #200)

- Fixed error with checking for any when NAs are present in areas that are water (@ChrisJones687, #201)

### Removed

- Removed dependencies: packages rgdal, raster, sp (@petrasovaa, #187)

## [2.0.1] 2022-11-16

### Added

- `calibrate` can now use the Mathews Correlation Coefficient as the summary statistic of 
  to keep or reject parameter sets (@ChrisJones687, #145).
  
- `validate`, `calibrate`, `pops_multirun`,  and `pops` now take use_survival_rates, survival_rate_month, 
  survival_rate_day, and survival_rates_file (@ChrisJones687, #147).
  
 - `validate`, `calibrate`, `pops_multirun`,  and `pops` now take network_movement as a parameter. This 
    parameter controls how dispersal occurs along the network (@ChrisJones687, #147).
    
### Changed
- `validate`, `calibrate`, `pops_multirun`, `auto_manage` and `pops` no longer take
  network_min_distance and network_max_distance as these parameters are now passed
  in through the parameter_means and parameter_cov_matrix parameters and are calibrated
  as part of the calibration if network kernel is selected (@ChrisJones687, #140).
  
- `calibrate` now calibrates the network_min_distance and network_max_distance parameters
  during calibration and they are now part of the parameter_means and parameter_cov_matrix
  that are exported from the calibration (@ChrisJones687, #140).
  
- `calibrate` has more flexible success metric options removes use_distance, use_rmse, and use_mcc
  parameters and replaces it with the more flexible success_metrics parameter. Users can now
  select multiple combinations of different success metrics for the calibration 
  (@ChrisJones687, #150).
  
- `pops_multirun` removed the ability to write all simulations. (@ChrisJones687, #144).

- the deterministic parameter has been renamed to dispersal_stochasticity. This 
  was done to be more consistent with generate_stochasticity, movement_stochasticity, and 
  establishment_stochasticity parameters (@ChrisJones687, #147).
  
- `validate`, `calibrate`, `pops_multirun`, and `pops` now propogate uncertainty from host and 
  initial conditions. This adds the parameters use_initial_condition_uncertainty and
  use_host_uncertainty. If use_initial_condition_uncertainty is TRUE the infected_file and/or
  exposed_file need to have 2 layers a mean and standard deviation. If use_host_uncertainty is
  is TRUE the host_file needs to have 2 layers a mean and standard deviation. For each model run
  a host and/or initial conditions are drawn from the mean and sd layers so that each run has a
  unique host and/or initial conditions this will allow for both propogation of uncertainty from 
  these sources but also partitioning (@ChrisJones687, #151).

### Fixed

- fixed error on parameter draws in select situations (@ChrisJones687, #146).
  
## [2.0.0] - 2021-12-14

### Added

- `quantity_allocation_disagreement` now takes in the use_distance parameter which is FALSE
  by default. This allows the model to commute the minimum total distance between observed
  and simulated infestations (@ChrisJones687, #130).
  
- `validate`, `calibrate`, `pops_multirun`, `auto_manage` and `pops` now take in 
  network_min_distance, network_max_distance, and network_filename parameters these are 
  used when the anthropogenic_kernel_type = "network". This allows directed spread along
  a network such as a railroad (@ChrisJones687 and @wenzeslaus, #131)

### Changed

- `validate` now takes the variable point_file and uses it to calculate statistics based on the
  point_file in addition to the raster file also calculates new measures of model performance
  accuracy, precision, recall, and specificity (@ChrisJones687, #124).
  
- `calibrate` no longer uses success_metrics and checks parameters as these are both handled
  internally by auto updating if the values are too far of in the first generation. This 
  makes for a simpler and faster user experience (@ChrisJones687, #130).
  
- `validate` no longer uses success_metrics parameter but now has the added parameters 
  use_rmse and use_distance that are used in `quantity_allocation_disagreement`.
  This allows for a more intuitive user interface (@ChrisJones687, #130).
  
- `quantity_allocation_disagreement` variable configuration changed to use_configuration
  to be more consistent with variable names. (@ChrisJones687, #130).
  
- `pops_multirun` can now use variable write_output = "all_simulations" to write out
  susceptible, exposed, and infected rasters from all simulations. This will allow for a 
  future update where simulations can easily be started from previous simulation outputs.
  (@ChrisJones687, #133)
  
## [1.1.0] - 2021-06-22

### Added

- New dispersal kernels added: Uniform, Power-law, Deterministic neighbor, 
  Hyperbolic-Secant, Gamma, Weibull, Normal, and Logistic (@ChrisJones687, #73).
  * These kernels are usable for both (`natural_dispersal_kernel`, `anthropogenic_dispersal_kernel`).
  * Used in the following functions (`pops`, `pops_multirun`, `auto_manage`, `calibrate`, and `validate`).

- Overpopulation module added to pops-core (@wenzeslaus, #83).

- Spatial Index to increase computational speed (@ChrisJones687, #67)

### Changed

- Exposed and Infected populations can both be present at the start of a simulation (@ChrisJones687, #92).

- Internal functions for data handling have switched from the raster package to the terra package (@ChrisJones687, #79).
  * Adds terra dependency.
  * Greatly speeds up data handling and preparation.

- Mask parameter now used in pops-multirun for post processing data for visualization (@ChrisJones687, #104).

- Calibration function now takes verbose parameter (@nfkruska, #99).

- All functions now support vrt data types (@ChrisJones687, #97).

- Raster files can now be read from S3 buckets (@ChrisJones687, #75).
  * needed for model-api for dashboard
  
- Outputs can now be saved with `write_outputs` and `output_folder_path` parameters (@ChrisJones687, #111).

- Host map now used as initial mask for post processing and validation calculations (@ChrisJones687, #112).

- Movements now moves exposed, resistant, and mortality tracked populations (@ChrisJones687, #118).

- Mortality can now occur at various timesteps not just yearly (@ChrisJones687, #118).
  * Can be either "day", "week", "month", "year", or "every_n_steps".
  * adds parameters `mortality_frequency` and `mortality_frequency_n`.

- Validation now exports the statistics for each output and the cumulative statistics for each year (@ChrisJones687, #121).

- Treatments now update total_hosts (@ChrisJones687, #122)

### Fixed

- Mask parameter works as intended in validate function after terra update (@ChrisJones687, #104).

- Improved pops_model documentation updating (@ChrisJones687, #94).

- `Output_frequency` can now be `every_n_steps` (@ChrisJones, #118).

## [1.0.2] - 2020-10-09

### Changed

- Exposed populations now exported each time-step (@ChrisJones687, #64).

### Fixed

- Multirun and validate draw a parameter set from the distribution for each run as intended

- Fixed error when calibrating for less than 4 parameters

## [1.0.1] - 2020-09-16

- [Patch release of r.pops.spread](https://github.com/ncsu-landscape-dynamics/r.pops.spread/releases/tag/v1.0.1) (no changes in rpops)


## [1.0.0] - 2020-09-16

Version 1.0.0 of the _PoPS Core_ C++ library and its interfaces: _rpops_ R package and _r.pops.spread_ GRASS GIS module. The release of rpops includes:

- Susceptible-infected (`SI`) and susceptible-exposed-infected (`SEI`) host phases (`model_type`, `latency_period`).

- Host mortality tracking (`mortality_rate`, `mortality`).

- Host removal and pesticide application treatments (`treatments`, `treatment_date`, `pesticide_duration`).

- Host resistance based on pesticide application treatments (`pesticide_duration` > 0).

- Treatments applied only to a ratio of hosts (`treatment_application`).

- Yearly pest removal based on lethal temperature (`lethal_temperature`, `lethal_month`).

- Two different dispersal kernels (`natural_dispersal_kernel`, `anthropogenic_dispersal_kernel`).

- Cauchy and exponential radial dispersal kernels use Von Mises distribution.

- Seasonal spread (`seasonality` in months).

- Multiple stochastic runs (`pops_multirun`).

- Parallel execution of multiple runs (`number_of_cores` in `pops_multirun`).

- Output of average and standard deviation of infected hosts across multiple runs and for a single stochastic run (`simulation_mean`,  `single_run`,  `simulation_sd`).

- Average and standard deviations for output averages (`number_infecteds`, `infected_areas`).

- Infection probability output in percent (`probability`).

- Spread rate measurement in 4 cardinal directions (`west_rate`, `east_rate`, `south_rate`, `north_rate` ).

- Distance to quarantine in 4 cardinal directions (`north_distance_to_quarantine`, `south_distance_to_quarantine`, `east_distance_to_quarantine`, `west_distance_to_quarantine` ).

- Probability of quarantine escape (`escape_probability`).


[unreleased]: https://github.com/ncsu-landscape-dynamics/rpops/compare/main...v2.0.1
[2.0.1]: https://github.com/ncsu-landscape-dynamics/rpops/compare/v2.0.0...v2.0.1
[2.0.0]: https://github.com/ncsu-landscape-dynamics/rpops/compare/1.1.0...v2.0.0
[1.1.0]: https://github.com/ncsu-landscape-dynamics/rpops/compare/v1.0.2...1.1.0
[1.0.2]: https://github.com/ncsu-landscape-dynamics/rpops/compare/v1.0.0...v1.0.2
[1.0.1]: https://github.com/ncsu-landscape-dynamics/rpops/releases/tag/v1.0.0
[1.0.0]: https://github.com/ncsu-landscape-dynamics/rpops/releases/tag/v1.0.0
