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

- New dispersal kernels added: Uniform, Power-law, Deterministic neighbor, 
  Hyperbolic-Secant, Gamma, Weibull, Normal, and Logistic (@ChrisJones687, #73).
  * These kernels are usable for both (`natural_dispersal_kernel`, `anthropogenic_dispersal_kernel`).
  * Used in the following functions (`pops`, `pops_multirun`, `auto_manage`, `calibrate`, and `validate`).

- Overpopulation module added to pops-core (@wenzeslaus, #83).

- Spatial Index to increase computational speed (@ChrisJones, #67)

### Changed
- Exposed and Infected populations can both be present at the start of a simulation (@ChrisJones687, #92).

- Internal functions for data handling have switched from the raster package to the terra package (@ChrisJones687, #79).
  * Adds terra dependency.
  * Greatly speeds up data handling and preparation.

- Mask parameter now used in pops-multirun for post processing data for visualization (@ChrisJones687, #104).

- Calibration function now takes verbose parameter (@nfkruska, #99).

- All functions now support vrt data types (@ChrisJones687, #97).

- Raster files can now be read from S3 buckets (@chrisJones687, #75).
  * needed for model-api for dashboard
  
- Outputs can now be saved with `write_outputs` and `output_folder_path` parameters (@ChrisJones687, #111)

- Host map now used as initial mask for postprocessing and validation calculations (@ChrisJones687, #112)

### Fixed

- Mask parameter works as intended in validate function after terra update (@ChrisJones687, #104).

- Improved pops_model documentation updating (@ChrisJones687, #94).

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


[unreleased]: https://github.com/ncsu-landscape-dynamics/rpops/compare/v1.0.2...master
[1.0.2]: https://github.com/ncsu-landscape-dynamics/rpops/compare/v1.0.0...v1.0.2
[1.0.1]: https://github.com/ncsu-landscape-dynamics/rpops/releases/tag/v1.0.0
[1.0.0]: https://github.com/ncsu-landscape-dynamics/rpops/releases/tag/v1.0.0
