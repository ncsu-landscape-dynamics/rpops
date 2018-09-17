#include <Rcpp.h>
using namespace Rcpp;
#include "pops/simulation.hpp"
#include "pops/raster.hpp"
#include <iostream>

using std::cout;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
NumericVector pops(NumericVector x) {
  Raster<int> infected = {{2, 0}, {0, 0}};
  Raster<int> exposed = {{0, 0}, {0, 0}};
  Raster<int> diseased = {{0, 0}, {0, 0}};
  Raster<int> infected_cohort = {{0, 0}, {0, 0}};
  Raster<int> susceptible = {{5, 6}, {14, 15}};
  Raster<int> total_plants = {{10, 6}, {14, 15}};
  Raster<double> temperature = {{5, 0}, {0, 0}};
  // Raster<double> weather_coef_value = {{0.5,0.7},{0.2,0.8}};
  std::vector<std::tuple<int, int>> outside_dispersers;
  Dispersal_kernel dispersal_kernel = CAUCHY;
  double lethal_temperature = -4.5;
  double reproductive_rate = 4.5;
  double short_distance_scale = 18.7;
  Simulation<Raster<int>, Raster<double>> simulation(42, infected);
  simulation.remove(infected, susceptible, exposed, diseased,
                    temperature, lethal_temperature);
  simulation.generate(infected, 0, reproductive_rate);
  simulation.disperse(susceptible, infected,
                      infected_cohort, total_plants,
                      outside_dispersers, dispersal_kernel, 0,
                      0, short_distance_scale);
  cout << infected;
  cout << outside_dispersers.size() << endl;
  return 0;
}


