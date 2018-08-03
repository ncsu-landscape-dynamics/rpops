## PoPS R wrapper
library(Rcpp)
Sys.setenv("PKF_CXXFLAGS"="-std=c++11")
sourceCpp("popss.cpp")
timesTwo(42)
