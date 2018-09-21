## PoPS R wrapper
library(Rcpp)
library(raster)
library(rgdal)
Sys.setenv("PKF_CXXFLAGS"="-std=c++11")
sourceCpp("pops.cpp")

## Read in data

infected = raster(matrix(c(5,0,0,0), ncol=2, nrow=2))
mortality_tracker = raster(matrix(c(0,0,0,0), ncol=2, nrow=2))
temperature = raster(matrix(c(5,0,0,0), ncol=2, nrow=2))

lethal_temperature = -4.5
random_seed = 42
reproductive_rate = 3.0
dispersal_kernel = "CAUCHY"
weather = TRUE
short_distance_scale = 0.0

pops_model(random_seed = random_seed, lethal_temperature = lethal_temperature, 
           reproductive_rate = reproductive_rate, 
           weather = weather, short_distance_scale = short_distance_scale)
