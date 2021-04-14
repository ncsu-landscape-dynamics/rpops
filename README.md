# PoPS (Pest or Pathogen Spread) R <img src="man/PoPS_GIF_transparent.gif" align="right" width="15%" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/ncsu-landscape-dynamics/rpops/workflows/R-CMD-check/badge.svg)](https://github.com/ncsu-landscape-dynamics/rpops/actions?query=workflow%3AR-CMD-check)
[![lint](https://github.com/ncsu-landscape-dynamics/rpops/workflows/lint/badge.svg)](https://github.com/ncsu-landscape-dynamics/rpops/actions?query=workflow%3Alint)
[![Build Status](https://travis-ci.org/ncsu-landscape-dynamics/rpops.svg?branch=master)](https://travis-ci.org/ncsu-landscape-dynamics/rpops)
[![codecov](https://codecov.io/gh/ncsu-landscape-dynamics/rpops/branch/master/graph/badge.svg)](https://codecov.io/gh/ncsu-landscape-dynamics/rpops)
  <!-- badges: end -->

## Overview

This is the R package for simulating spread of pests and pathogens. The package is an R package with multiple functions built around the PoPS (Pest or Pathogen Spread) model implemented in the C++ library maintained in the [PoPS Core Repository](https://github.com/ncsu-landscape-dynamics/pops-core). 

PoPs is a stochastic spread model of pests and pathogens in forest and agricultural landscapes to learn more visit [popsmodel.org](https://popsmodel.org/). The R package provides an easy way for researchers to calibrate, validate, and test what if scenarios of treatment interventions. The model is also available in GRASS GIS you can install and use [r.pops.spread](https://github.com/ncsu-landscape-dynamics/r.pops.spread) to run the model in GRASS GIS.

## Features
The PoPS Core C++ library and its interfaces: rpops R package and r.pops.spread GRASS GIS module. The release of rpops includes:

* Susceptible-infected (`SI`) and susceptible-exposed-infected (`SEI`) host phases (`model_type`, `latency_period`),
* Host mortality tracking (`mortality_rate`, `mortality`),
* Host removal and pesticide application treatments (`treatments`, `treatment_date`, `pesticide_duration`),
* Host resistance based on pesticide application treatments (`pesticide_duration` > 0),
* Treatments applied only to a ratio of hosts (`treatment_application`),
* Yearly pest removal based on lethal temperature (`lethal_temperature`, `lethal_month`),
* Two different dispersal kernels (`natural_dispersal_kernel`, `anthropogenic_dispersal_kernel`),
* Cauchy, Exponential, Uniform, Power-law, Deterministic neighbor, Hyperbolic-Secant, Gamma, Weibull, Normal, and Logistic radial dispersal kernels use Von Mises distribution,
* Seasonal spread (`seasonality` in months),
* Multiple stochastic runs (`pops_multirun`),
* Parallel execution of multiple runs (`number_of_cores` in `pops_multirun`),
* Output of average and standard deviation of infected hosts across multiple runs and for a single stochastic run (`simulation_mean`,  `single_run`,  `simulation_sd`, `simulation_min`, `simulation_max`),
* Average and standard deviations for output averages (`number_infecteds`, `infected_areas`),
* Infection probability output in percent (this is the percent of model runs that a cell has at least one infestation/infection (`probability`),
* Spread rate measurement in 4 cardinal directions (`west_rate`, `east_rate`, `south_rate`, `north_rate` ),
* Distance to quarantine in 4 cardinal directions (`north_distance_to_quarantine`, `south_distance_to_quarantine`, `east_distance_to_quarantine`, `west_distance_to_quarantine` ),
* Probability of quarantine escape (`escape_probability`).
* Overpopulation function (individuals in areas of high population leave the area and disperse longer distances on average)
* Calibration of the model parameters using either MCMC (markov chain monte carlo) or ABC (approximate bayesian computation) (`calibrate`)
* Validation of the model using quantiy, allocation, and configuration disagreement (`validate`)
## How to cite

If you use this software or code, please cite the following papers:

* Ross K. Meentemeyer, Nik J. Cunniffe, Alex R. Cook, Joao A. N. Filipe,
  Richard D. Hunter, David M. Rizzo, and Christopher A. Gilligan, 2011.
  Epidemiological modeling of invasion in heterogeneous landscapes:
  spread of sudden oak death in California (1990–2030).
  *Ecosphere* 2:art17.
  [DOI: 10.1890/ES10-00192.1](https://doi.org/10.1890/ES10-00192.1)

* Tonini, Francesco, Douglas Shoemaker, Anna Petrasova, Brendan Harmon,
  Vaclav Petras, Richard C. Cobb, Helena Mitasova,
  and Ross K. Meentemeyer, 2017.
  Tangible geospatial modeling for collaborative solutions
  to invasive species management.
  *Environmental Modelling & Software* 92: 176-188.
  [DOI: 10.1016/j.envsoft.2017.02.020](https://doi.org/10.1016/j.envsoft.2017.02.020)

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

## How to install

Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) if not already installed. Once Rtools is installed you can the latest version of PoPS by 
running the code below and changing the ref to match the latest semantic version
or the specific version that you are looking to install. PoPS requires [terra](https://github.com/rspatial/terra) version 1.1-17 or higher.

```R
install.packages("devtools")
library(devtools)
devtools::install_github("ncsu-landscape-dynamics/rpops", ref = "v1.0.2")
library(PoPS)

```

## Contributing


Please see the [pops-core](https://github.com/ncsu-landscape-dynamics/pops-core#readme) repository for contributing best practices and release policies. Other than that, just open pull requests against this repo. We suggest following the [Style Guide from Hadley](http://adv-r.had.co.nz/Style.html) for R code.

## Authors and contributors

### Authors

_(alphabetical order)_

* Chris Jones
* Margaret Lawrimore
* Vaclav Petras
* Anna Petrasova

### Previous contributors

_(alphabetical order)_

* Zexi Chen
* Devon Gaydos
* Francesco Tonini

See Git commit history, GitHub insights, or CHANGELOG.md file (if present)
for details about contributions.

## License

Permission to use, copy, modify, and distribute this software and its documentation
under the terms of the GNU General Public License version 2 or higher is hereby
granted. No representations are made about the suitability of this software for any
purpose. It is provided "as is" without express or implied warranty.
See the GNU General Public License for more details.
