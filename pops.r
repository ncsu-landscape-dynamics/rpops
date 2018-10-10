## PoPS R wrapper
library(Rcpp)
library(raster)
library(rgdal)
# library(velox)
# library(rgrass7)
# library(link2GI)
##link = findGRASS()
Sys.setenv("PKF_CXXFLAGS"="-std=c++11")
sourceCpp("pops.cpp")

# initGRASS(gisBase = "C:/Program Files/GRASS GIS 7.4.1",
#           home = getwd(), gisDbase = "GRASS_TEMP", override = TRUE)

## Read in data
# sus = "C:/Users/Chris/Desktop/sus.tif"
infected = raster(matrix(c(5,0,0,0), ncol=2, nrow=2))
sus = matrix(c(10,6,14,15), ncol=2, nrow=2)
#writeRaster(x= sus, filename = "C:/Users/Chris/Desktop/sus.tif", overwrite = TRUE, format = 'GTiff')
total_plants = raster(matrix(c(15,6,14,15), ncol=2, nrow=2))
mortality_tracker = raster(matrix(c(0,0,0,0), ncol=2, nrow=2))
weather_coefficient = raster(matrix(c(5,0,0,0), ncol=2, nrow=2))

cols = as.numeric(ncol(sus))
rows = as.numeric(nrow(sus))
ew_res = xres(sus)
ns_res = yres(sus)
lethal_temperature = -4.5
random_seed = 42
reproductive_rate = 7.0
dispersal_kernel = "CAUCHY"
weather = TRUE
short_distance_scale = 0.5

pops_model(random_seed = random_seed, lethal_temperature = lethal_temperature, 
           reproductive_rate = reproductive_rate, 
           weather = weather, short_distance_scale = short_distance_scale,
           sus = sus)
