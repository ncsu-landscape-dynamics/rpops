# PoPS (Pest or Pathogen Spread) R Package <img src="man/PoPS_GIF_transparent.gif" align="right" width="15%" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/ncsu-landscape-dynamics/rpops/workflows/R-CMD-check/badge.svg)](https://github.com/ncsu-landscape-dynamics/rpops/actions?query=workflow%3AR-CMD-check)
[![lint](https://github.com/ncsu-landscape-dynamics/rpops/workflows/lint/badge.svg)](https://github.com/ncsu-landscape-dynamics/rpops/actions?query=workflow%3Alint)
[![Build Status](https://travis-ci.org/ncsu-landscape-dynamics/rpops.svg?branch=master)](https://travis-ci.org/ncsu-landscape-dynamics/rpops)
[![codecov](https://codecov.io/gh/ncsu-landscape-dynamics/rpops/branch/master/graph/badge.svg)](https://codecov.io/gh/ncsu-landscape-dynamics/rpops)
  <!-- badges: end -->

## Overview

This is the R package for simulating spread of pests and pathogens. The package is an R package with multiple functions built around the PoPS (Pest or Pathogen Spread) model implemented in the C++ library maintained in the [PoPS Core Repository](https://github.com/ncsu-landscape-dynamics/pops-core). 

PoPs is a stochastic spread model of pests and pathogens in forest and agricultural landscapes to learn more visit [popsmodel.org](https://popsmodel.org/). The R package provides an easy way for researchers to calibrate, validate, and test what if scenarios of treatment interventions. The model is also available in GRASS GIS you can install and use [r.pops.spread](https://github.com/ncsu-landscape-dynamics/r.pops.spread) to run the model in GRASS GIS.

## Installation
Requires terra version between 1.1.17 and 1.2.5. terra 1.2.10 has a breaking bug in the extract function.

If you are on Windows, you need to first install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) to get a C++ compiler that R can use.

### Version Specific

If you want to install a specific version just change the version number.

```R
install.packages("remotes")
remotes::install_github("ncsu-landscape-dynamics/rpops", ref = "v1.0.2")
library(PoPS)

```
### Development Version

```R
install.packages("remotes")
remotes::install_github("ncsu-landscape-dynamics/rpops")
library(PoPS)

```
## Features
the PoPS package in R is built on top of PoPS Core C++ library includes:

* Susceptible-infected (`SI`) and susceptible-exposed-infected (`SEI`) model types (`model_type`, `latency_period`).
* Host mortality tracking (`mortality_rate`, `mortality_time_lag`, `mortality`).
* Host removal and pesticide application treatments (`treatments`, `treatment_date`, `pesticide_duration`).
* Host resistance based on pesticide application treatments (`pesticide_duration` > 0).
* Treatments applied only to a ratio of hosts (`treatment_application`).
* Yearly pest removal based on lethal temperature (`lethal_temperature`, `lethal_month`).
* Two different dispersal kernels (`natural_dispersal_kernel`, `anthropogenic_dispersal_kernel`).
* Cauchy, Exponential, Uniform, Power-law, Deterministic neighbor, Hyperbolic-Secant, Gamma, Weibull, Normal, and Logistic radial dispersal kernels use Von Mises distribution.
* Seasonal spread (`seasonality` in months).
* Host movement for animals moving from farm to farm or plants via nursery trade (`use_movements`, `movements_file`).
* Reduced stochasticity options (`generate_stochasticity`, `establishment_stochasticity`, `movement_stochasticity`) and deterministic versions of kernels (`deterministic`) other required parameters if reducing stochasticity are (`dispersal_percentage`, `establishment_probability`).
* Spread rate measurement in 4 cardinal directions (`west_rate`, `east_rate`, `south_rate`, `north_rate`) when (`use_spreadrates`) is true.
* Distance to quarantine in 4 cardinal directions (`north_distance_to_quarantine`, `south_distance_to_quarantine`, `east_distance_to_quarantine`, `west_distance_to_quarantine`) when (`use_quarantine`) is true and (`quarantine_areas_file`) is provided.
* Probability of quarantine escape (`escape_probability`).
* Overpopulation function (individuals in areas of high population leave the area and disperse longer distances on average) (`use_overpopulation`, `overpopulation_percentage`, `leaving_percentage`, `leaving_scale_coefficient`).
* Flexible output frequency with n number of days, weeks, months, or years as options (`output_frequency`, `output_frequency_n`).

### Functions in rpops
* `calibrate:` Calibration of the model parameters using either MCMC (markov chain monte carlo) or ABC (approximate bayesian computation). 
* `validate:` Validation of the model using quantiy, allocation, and configuration disagreement.
* `pops_multirun:` Parallel execution of multiple stochastic runs (`number_of_cores` used to set cores used if left NULL defaults to using n - 1 cores on the machine). Outputs statistics of infected/infested hosts across multiple runs (`simulation_mean`,  `single_run`,  `simulation_sd`, `simulation_min`, `simulation_max`), current state using the median (`infected`,  `exposed`, and `susceptible`), average and standard deviations whole area statistics (`number_infecteds`, `infected_areas`), and probability of infection (`probability`) which is the percent of model runs that a cell has at least one infestation/infection.
* `pops:` Runs a single stochastic run of the model. This function is primarily used for automated testing of model functionality.

## How to cite

If you use this software or code, please cite the following papers:

* Jones, Chris, Jones, Shannon, Petrasova, Anna, Petras, Vaclav, Gaydos, Devon, 
  Skrip, Megan, Yu, Takeuchi, Bigsby, Kevin, and Meentemeyer, Ross, 2021.
  Iteratively forecasting biological invasions with PoPS and a little help from 
  our friends.
  *Frontiers in Ecology and the Environment* 
  [DOI: 10.1002/fee.2357](https://doi.org/10.1002/fee.2357)

In case you are using the automatic management feature in rpops or the
steering version of r.pops.spread (from the branch steering), please
cite also:

* Petrasova, A., Gaydos, D.A., Petras, V., Jones, C.M., Mitasova, H. and
  Meentemeyer, R.K., 2020.
  Geospatial simulation steering for adaptive management.
  *Environmental Modelling & Software* 133: 104801.
  [DOI: 10.1016/j.envsoft.2020.104801](https://doi.org/10.1016/j.envsoft.2020.104801)

In addition to citing the above paper, we also encourage you to
reference, link, and/or acknowledge specific version of the software
you are using for example:

* *We have used rpops R package version 1.0.0 from
  <https://github.com/ncsu-landscape-dynamics/rpops>*.

## Contributing

Please see the [pops-core](https://github.com/ncsu-landscape-dynamics/pops-core#readme) repository for contributing best practices and release policies. Other than that, just open pull requests against this repo. We suggest following the [Style Guide from Hadley](http://adv-r.had.co.nz/Style.html) for R code.

## Authors and contributors

### Authors

_(alphabetical order)_

* Chris Jones
* Margaret Lawrimore
* Vaclav Petras
* Anna Petrasova

### Contributors

_(alphabetical order)_

* Zexi Chen
* Devon Gaydos
* Nick Kruskamp
* Francesco Tonini

See Git commit history, GitHub insights, or CHANGELOG.md file for details about
contributions.

## License

Permission to use, copy, modify, and distribute this software and its documentation
under the terms of the GNU General Public License version 2 or higher is hereby
granted. No representations are made about the suitability of this software for any
purpose. It is provided "as is" without express or implied warranty.
See the GNU General Public License for more details.
