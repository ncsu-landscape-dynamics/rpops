## PoPS R wrapper
library(Rcpp)
library(raster)
library(rgdal)
Sys.setenv("PKF_CXXFLAGS"="-std=c++11")
sourceCpp("pops.cpp")
pops_model(42)
