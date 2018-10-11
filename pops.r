## PoPS R wrapper
library(Rcpp)
library(raster)
library(rgdal)
Sys.setenv("PKF_CXXFLAGS"="-std=c++11")
sourceCpp("pops.cpp")

## Read in data
infected_file = 

infected = raster(infected_file)
susceptible = raster(susceptible_file)
total_plants = raster(total_plants_file)
temperature = raster(temperature_file)
weather_coefficient = raster(weather_coefficient_file)

mortality_tracker = raster(mortality_tracker_file)

cols = as.numeric(ncol(susceptible))
rows = as.numeric(nrow(susceptible))
ew_res = xres(susceptible)
ns_res = yres(susceptible)
lethal_temperature = -4.5
random_seed = 42
reproductive_rate = 8.0
dispersal_kernel = "CAUCHY"
weather = TRUE
short_distance_scale = 1.5
use_lethal_temperature = TRUE
lethal_temperature_month = 1
infected = as.matrix(infected)
susceptible = as.matrix(susceptible)
total_plants = as.matrix(total_plants)
mortality_tracker = as.matrix(mortality_tracker)
temperature = as.matrix(temperature)
weather_coefficient = as.matrix(weather_coefficient)

pops_model(random_seed = random_seed, 
           lethal_temperature = lethal_temperature, use_lethal_temperature, lethal_temperature_month,
           reproductive_rate = reproductive_rate, 
           weather = weather, short_distance_scale = short_distance_scale, infected = infected,
           susceptible = susceptible, mortality_tracker =mortality_tracker,
           total_plants = total_plants, temperature = temperature,
           weather_coefficient = weather_coefficient, ew_res = ew_res, ns_res = ns_res)
