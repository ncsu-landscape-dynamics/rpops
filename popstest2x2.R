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

season_month_start = 6
season_month_end = 11
time_step = "month"
start_time = 2018
end_time = 2020

number_of_years = end_time-start_time+1
temperature_stack <- stack(lapply(1:number_of_years, function(i) temperature))
temperature = list(as.matrix(temperature_stack[[1]]))
for(i in 2:number_of_years) {
   temperature[[i]] <- as.matrix(temperature_stack[[i]])
}

if (time_step == "week") {
  number_of_time_steps = (end_time-start_time+1)*52 +1
} else if (time_step == "month") {
  number_of_time_steps = (end_time-start_time+1)*12
} else if (time_step == "day") {
  number_of_time_steps = (end_time-start_time+1)*365
}
weather_coefficient_stack <- stack(lapply(1:number_of_time_steps, function(i) weather_coefficient))
weather_coefficient <- list(as.matrix(weather_coefficient_stack[[1]]))
for(i in 2:number_of_time_steps) {
  weather_coefficient[[i]] <- as.matrix(weather_coefficient_stack[[i]])
}

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
use_lethal_temperature = TRUE
lethal_temperature_month = 1
infected = as.matrix(infected)
susceptible = as.matrix(susceptible)
total_plants = as.matrix(total_plants)
mortality_tracker = as.matrix(mortality_tracker)
weather_coefficient = as.matrix(weather_coefficient)
dispersal_kern = "cauchy"
percent_short_distance_dispersal = 1.0
long_distance_scale = 0.0
wind_dir = "NONE"
kappa = 0


data <- pops_model(random_seed = random_seed, 
           lethal_temperature = lethal_temperature, use_lethal_temperature, lethal_temperature_month,
           reproductive_rate = reproductive_rate, 
           weather = weather, short_distance_scale = short_distance_scale, infected = infected,
           susceptible = susceptible, mortality_tracker =mortality_tracker,
           total_plants = total_plants, temperature = temperature,
           weather_coefficient = weather_coefficient, 
           ew_res = ew_res, ns_res = ns_res,
           time_step = time_step,
           season_month_start = season_month_start, season_month_end = season_month_end,
           start_time = start_time, end_time = end_time,
           dispersal_kern = dispersal_kern, percent_short_distance_dispersal = percent_short_distance_dispersal,
           long_distance_scale = long_distance_scale,
           wind_dir = wind_dir, kappa = kappa)
