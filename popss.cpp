#include <Rcpp.h>
using namespace Rcpp;
#include "popss/simulation.hpp"
#include "popss/raster.hpp"
#include <iostream>

using std::cout;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  Raster<int> infected = {{5, 0}, {0, 0}};
  Raster<int> infected_cohort = {{0, 0}, {0, 0}};
  Raster<int> susceptible = {{5, 6}, {14, 15}};
  Raster<int> all = {{5, 0}, {0, 0}};
  Raster<double> temperature = {{5, 0}, {0, 0}};
  std::vector<std::tuple<int, int>> escaping;
  Rtype rtype = CAUCHY;
  double lethal_temperature = -4.5;
  double weather_coef_value = 5.8;
  double spore_rate = 4.5;
  double scale1 = 18.7;
  Simulation<Raster<int>, Raster<double>> simulation(42, infected);
  simulation.SporeRemove(infected, susceptible,
                         temperature, lethal_temperature);
  simulation.SporeGen(infected, 0, weather_coef_value, spore_rate);
  simulation.SporeSpreadDisp_singleSpecies(susceptible, infected,
                                           infected_cohort, all,
                                           escaping, rtype, 0,
                                           weather_coef_value, scale1);
  cout << infected;
  cout << escaping.size() << endl;
  return 0;
}


