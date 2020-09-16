# PoPS (Pest or Pathogen Spread) Model R Wrapper <img src="man/PoPS_Logo.png" align="right" width="14%" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/ncsu-landscape-dynamics/rpops/workflows/R-CMD-check/badge.svg)](https://github.com/ncsu-landscape-dynamics/rpops/actions?query=workflow%3AR-CMD-check)
[![lint](https://github.com/ncsu-landscape-dynamics/rpops/workflows/lint/badge.svg)](https://github.com/ncsu-landscape-dynamics/rpops/actions?query=workflow%3Alint)
[![Build Status](https://travis-ci.org/ncsu-landscape-dynamics/rpops.svg?branch=master)](https://travis-ci.org/ncsu-landscape-dynamics/rpops)
[![codecov](https://codecov.io/gh/ncsu-landscape-dynamics/rpops/branch/master/graph/badge.svg)](https://codecov.io/gh/ncsu-landscape-dynamics/rpops)
  <!-- badges: end -->

## Overview

This is the R package for simulating spread of pests and pathogens. The package is an R package with multiple functions built around the PoPS (Pest or Pathogen Spread) model implemented in the C++ library maintained in the [PoPS Core Repository](https://github.com/ncsu-landscape-dynamics/pops-core). 

PoPs is a stochastic spread model of pests and pathogens in forest and agricultural landscapes. The R package provides an easy way for researchers to calibrate, validate, and test what if scenarios of treatment interventions. The model is also available in GRASS GIS you can install and use [r.pops.spread](https://github.com/ncsu-landscape-dynamics/r.pops.spread) to run the model in GRASS GIS.

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

Install Rtools(https://cran.r-project.org/bin/windows/Rtools/) if not already installed. 
```R
install.packages("devtools")
library(devtools)
devtools::install_github("ncsu-landscape-dynamics/rpops")
library(PoPS)

## if you get an error that says it failed check use:
Sys.getenv("BINPREF")
## should return "C:/Rtools/mingw_$(WIN)/bin/"
## if not run the command below
cat('Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")',
    file = file.path(Sys.getenv("HOME"), ".Rprofile"), 
    sep = "\n", append = FALSE)
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
