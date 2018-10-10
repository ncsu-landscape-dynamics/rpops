## PoPS R wrapper
library(Rcpp)
library(raster)
library(rgdal)
Sys.setenv("PKF_CXXFLAGS"="-std=c++11")
sourceCpp("pops.cpp")

## Read in data
infected = matrix(c(5,0,0,0), ncol=2, nrow=2)
infected = raster(infected, xmn = 0, ymn = 0, xmx = 60, ymx = 60)
susceptible = matrix(c(10,6,14,15), ncol=2, nrow=2)
susceptible = raster(susceptible, xmn = 0, ymn = 0, xmx = 60, ymx = 60)
total_plants = matrix(c(15,6,14,15), ncol=2, nrow=2)
total_plants = raster(total_plants, xmn = 0, ymn = 0, xmx = 60, ymx = 60)
mortality_tracker = matrix(c(0,0,0,0), ncol=2, nrow=2)
mortality_tracker = raster(mortality_tracker, xmn = 0, ymn = 0, xmx = 60, ymx = 60)
temperature = matrix(c(5,0,0,0), ncol=2, nrow=2)
temperature = raster(temperature, xmn = 0, ymn = 0, xmx = 60, ymx = 60)
weather_coefficient = matrix(c(0.8,0.5,0.9,0.2), ncol=2, nrow=2)
weather_coefficient = raster(weather_coefficient, xmn = 0, ymn = 0, xmx = 60, ymx = 60)

cols = as.numeric(ncol(susceptible))
rows = as.numeric(nrow(susceptible))
ew_res = xres(susceptible)
ns_res = yres(susceptible)
lethal_temperature = -4.5
random_seed = 42
reproductive_rate = 3.0
dispersal_kernel = "CAUCHY"
weather = TRUE
short_distance_scale = 1.5
infected = as.matrix(infected)
susceptible = as.matrix(susceptible)
total_plants = as.matrix(total_plants)
mortality_tracker = as.matrix(mortality_tracker)
temperature = as.matrix(temperature)
weather_coefficient = as.matrix(weather_coefficient)

pops_model(random_seed = random_seed, lethal_temperature = lethal_temperature, 
           reproductive_rate = reproductive_rate, 
           weather = weather, short_distance_scale = short_distance_scale, infected = infected,
           susceptible = susceptible, mortality_tracker =mortality_tracker,
           total_plants = total_plants, temperature = temperature,
           weather_coefficient = weather_coefficient, ew_res = ew_res, ns_res = ns_res)
