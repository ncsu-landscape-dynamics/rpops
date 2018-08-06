## PoPS R wrapper
library(Rcpp)
Sys.setenv("PKF_CXXFLAGS"="-std=c++11")
sourceCpp("pops.cpp")
timesTwo(42)
