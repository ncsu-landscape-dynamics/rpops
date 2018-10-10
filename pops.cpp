// #define POPS_RASTER_WITH_GRASS_GIS

#include <Rcpp.h>
#include "pops/simulation.hpp"
#include "pops/raster.hpp"
#include "pops/date.hpp"
#include <iostream>
#include <vector>
#include <map>
#include <memory>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>


// extern "C"{
// #include "grass/gis.h"
// #include "grass/glocale.h"
// #include "grass/raster.h"
// }

using std::string;
using std::cout;
using std::cerr;
using std::endl;

using namespace Rcpp;
using namespace pops;

// typedef Raster<int> IntegerRaster;
// typedef Raster<double> DoubleRaster;
// #define POPS_RASTER_WITH_GRASS_GIS

class Season
{
public:
    Season(int start, int end)
        : m_start_month(start), m_end_month(end)
    {}
    inline bool month_in_season(int month)
    {
        return month >= m_start_month && month <= m_end_month;
    }
private:
    int m_start_month;
    int m_end_month;
};

// RCPP_MODULE(rast) {
//   class_<Raster<int>>("int_rast");
//   
// }

// RcppExport SEXP Raster__new(SEXP rast_) {
//   Raster<int> rast = as<Raster<int>>(rast_);
//   Rcpp::XPtr<Raster<int>>
//   ptr( new Raster<int>(rast), true);
//   return ptr;
// }

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List pops_model(int random_seed, double lethal_temperature,
                double reproductive_rate,
                bool weather,
                double short_distance_scale,
                IntegerMatrix sus,
                int cols,
                int rows,
                int ew_res,
                int ns_res,
                double start_time = 2018, double end_time = 2018,
                double percent_short_distance_dispersal = 0.0,
                double long_distance_scale = 0.0
                )
  {
  // Environment raster("package:raster");
  // Function rast = raster["raster"];
  
  //Raster<int> suscept = Raster<int>::from_grass_raster(sus);
  // infected = rast(infected);
  Raster<int> susceptible = Raster<int>(cols, rows, ew_res, ns_res);
  // mortality_tracker = rast(mortality_tracker);
  // total_plants = rast(total_plants);
  // temperature = rast(temperature);
  // weather_coefficient = rast(weather_coefficient);
  
  Raster<int> infected = {{5, 0}, {0, 0}};
  Raster<int> mortality_tracker = {{0, 0}, {0, 0}};
  // Raster<int> susceptible = {{10, 6}, {14, 15}};
  Raster<int> total_plants = {{15, 6}, {14, 15}};
  Raster<double> temperature = {{5, 0}, {0, 0}};
  Raster<double> weather_coefficient = {{0.8, 0.8}, {0.2, 0.8}};
  std::vector<std::tuple<int, int>> outside_dispersers;
  DispersalKernel dispersal_kernel = CAUCHY;
  
  pops::Date dd_start(start_time, 01, 01);
  pops::Date dd_end(end_time, 12, 31);
  Season season(6,11);
  string step = "month";
  pops::Date dd_current(dd_start);

  int counter = 0;
  
for (int current_time_step = 0; ; current_time_step++, step == "month" ? dd_current.increased_by_month() : dd_current.increased_by_week()) {
    if (season.month_in_season(dd_current.month()))
    counter += 1;
    if (dd_current >= dd_end)
            break;
}

  // Simulation<IntegerMatrix, NumericMatrix> simulation(random_seed, infected);
  Simulation<Raster<int>, Raster<double>> simulation(random_seed, infected);
  simulation.remove(infected, susceptible,
                    temperature, lethal_temperature);
  simulation.generate(infected, weather, weather_coefficient, reproductive_rate);
  simulation.disperse(susceptible, infected,
                      mortality_tracker, total_plants,
                      outside_dispersers, weather, weather_coefficient, 
                      dispersal_kernel, short_distance_scale);
  // this is a test of vector of rasters
  // std::vector<IntegerMatrix> v;
  std::vector<Raster<int>> v;
  v = {infected,susceptible};
   // this is a test of vector of rasters
   
  cout << infected;
  cout << outside_dispersers.size() << endl;
  for(Raster<int> n : v) {
    cout << n;  // this is a test of vector of rasters
  }
  cout << counter;
  cout << dd_start;
  cout << dd_end;
  cout << dd_current;
  //cout << sus;
  
  return 0;
}

