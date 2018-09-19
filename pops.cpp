#include <Rcpp.h>
using namespace Rcpp;
#include "pops/simulation.hpp"
#include "pops/raster.hpp"
#include <iostream>

using std::cout;
using namespace pops;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
NumericVector pops_model(NumericVector x) {
  Raster<int> infected = {{5, 0}, {0, 0}};
  Raster<int> mortality_tracker = {{0, 0}, {0, 0}};
  Raster<int> susceptible = {{10, 6}, {14, 15}};
  Raster<int> total_plants = {{15, 6}, {14, 15}};
  Raster<double> temperature = {{5, 0}, {0, 0}};
  Raster<double> weather_coefficient = {{0.8, 0.8}, {0.2, 0.8}};
  std::vector<std::tuple<int, int>> outside_dispersers;
  DispersalKernel dispersal_kernel = CAUCHY;
  bool weather = true;
  double lethal_temperature = -4.5;
  double reproductive_rate = 4.5;
  double short_distance_scale = 0.0;
  Simulation<Raster<int>, Raster<double>> simulation(42, infected);
  simulation.remove(infected, susceptible,
                    temperature, lethal_temperature);
  simulation.generate(infected, weather, weather_coefficient, reproductive_rate);
  simulation.disperse(susceptible, infected,
                      mortality_tracker, total_plants,
                      outside_dispersers, weather, weather_coefficient, 
                      dispersal_kernel, short_distance_scale);
  cout << infected;
  cout << outside_dispersers.size() << endl;
  return 0;
}


