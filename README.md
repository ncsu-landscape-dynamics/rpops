# PoPS (Pest or Pathogen Spread) Model R Wrapper <img src="man/PoPS_Logo.png" align="right" width="14%" />

[![Build Status](https://travis-ci.org/ncsu-landscape-dynamics/rpops.svg?branch=master)](https://travis-ci.org/ncsu-landscape-dynamics/rpops)
[![Build status](https://ci.appveyor.com/api/projects/status/pixncc1jo7gtu0wx/branch/master?svg=true)](https://ci.appveyor.com/project/ChrisJones687/rpops/branch/master)
[![codecov](https://codecov.io/gh/ncsu-landscape-dynamics/rpops/branch/master/graph/badge.svg)](https://codecov.io/gh/ncsu-landscape-dynamics/rpops)

## Overview

An R wrapper using RCPP for the [PoPS C++ library](https://github.com/ncsu-landscape-dynamics/PoPShttps://github.com/ncsu-landscape-dynamics/PoPS). PoPs is a stochastic spread model of pests and pathogens in forest and agricultural landscapes. It has been generalized and new features added but was originally developed for *Phytophthora ramorum* and the original version of the model was based on this reference paper: Ross K. Meentemeyer, Nik J. Cunniffe, Alex R. Cook, Joao A. N. Filipe, Richard D. Hunter, David M. Rizzo, and Christopher A. Gilligan 2011. Epidemiological modeling of invasion in heterogeneous landscapes: spread of sudden oak death in California (1990–2030). *Ecosphere* 2:art17. [http://dx.doi.org/10.1890/ES10-00192.1] (http://www.esajournals.org/doi/abs/10.1890/ES10-00192.1) 

## How to install

Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) if not already installed. 

```R
install.packages("devtools")
library(devtools)
devtools::install_github("ncsu-landscape-dynamics/rpops")
# if you want a specific branch use where branch_name is the name of the branch you want to use.
devtools::install_github("ncsu-landscape-dynamics/rpops", ref = "branch_name")
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

This section is designed to clarify the branch structure of this repository and where new features and bug fixes should go.

### Branch Structure

1. **master** is the stable version of the model that is used for official releases. 
2. **bugfix/thingnotworking** are branched off of master then merged back once bug is fixed.
3. **feature/new_feature** is where new features are developed before they are merged into master. For example, infect and vector are currently being developed and will be merged together prior to being merged to master for an official major version release.

### Bug Fixes

Most bugs/issues will be found in the **master** branch as it is the branch being used in the R package and Grass Module. Thus bug fixes should be merged into **master** once tested. Bug fixes should be released as minor versions (e.g. if major release is 1.0 then the first bug fix would be released as version 1.1).

### New Features

When creating new features create a branch from **master** using the following syntax **feature/new_feature**. For example, we want to add a transportation network model for human assisted dispersal, the branch created would be named feature/transportation_network_model (or similar). New features will be merged into **master** once tested based on the priorities of our stakeholders first. Once new features are tested with the latest bug fixes and any other new features being included in the next major release we will merge them into **master** and create an official major release version (e.g. update from version 1.1 to version 2.0). 

If you are interested in contributing to PoPS development and are not a core developer on the model, please take a look at following
documents to make the process as seamless as possible.

1. [Contributor Code of Conduct](contributing_docs/CODE_OF_CONDUCT.md)
1. [Style Guide from hadley](http://adv-r.had.co.nz/Style.html)
1. [Contributor Guide](contributing_docs/CONTRIBUTING.md)


## Authors

* Chris Jones
* Anna Petrasova
* Vaclav Petras
* Devon Gaydos

## License

Permission to use, copy, modify, and distribute this software and its documentation under the terms of the GNU General Public License is hereby granted. No representations are made about the suitability of this software for any purpose. It is provided "as is" without express or implied warranty. See the GNU General Public License for more details.
