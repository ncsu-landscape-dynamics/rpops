#include <Rcpp.h>
#include "simulation.hpp"
#include "raster.hpp"
#include "date.hpp"
#include "treatments.hpp"
#include <iostream>
#include <vector>
#include <map>
#include <memory>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>

using std::string;
using std::cout;
using std::cerr;
using std::endl;

using namespace Rcpp;
using namespace pops;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

bool all_infected(IntegerMatrix susceptible)
{
  bool allInfected = true;
  for (int j = 0; j < susceptible.rows(); j++) {
    for (int k = 0; k < susceptible.cols(); k++) {
      if (susceptible(j, k) > 0)
        allInfected = false;
    }
  }
  return allInfected;
}

Direction direction_enum_from_string(const std::string& text)
{
  std::map<std::string, Direction> mapping{
    {"N", N}, {"NE", NE}, {"E", E}, {"SE", SE}, {"S", S},
    {"SW", SW}, {"W", W}, {"NW", NW}, {"NONE", NONE}
  };
  try {
    return mapping.at(text);
  }
  catch (const std::out_of_range&) {
    throw std::invalid_argument("direction_enum_from_string: Invalid"
                                  " value '" + text +"' provided");
  }
}

DispersalKernel radial_type_from_string(const std::string& text)
{
  if (text == "cauchy")
    return CAUCHY;
  else if (text == "cauchy_double_scale")
    return CAUCHY_DOUBLE_SCALE;
  else
    throw std::invalid_argument("radial_type_from_string: Invalid"
                                  " value '" + text +"' provided");
}

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List pops_model(int random_seed, 
                double lethal_temperature, bool use_lethal_temperature, int lethal_temperature_month,
                double reproductive_rate,
                bool weather, bool mortality_on,
                double short_distance_scale,
                IntegerMatrix infected,
                IntegerMatrix susceptible,
                IntegerMatrix mortality_tracker,
                IntegerMatrix mortality,
                IntegerMatrix total_plants,
                std::vector<NumericMatrix> treatment_maps,
                std::vector<int> treatment_years,
                std::vector<NumericMatrix> temperature,
                std::vector<NumericMatrix> weather_coefficient,
                int ew_res, int ns_res,
                std::string time_step,
                double mortality_rate = 0.0, int mortality_time_lag = 2,
                int season_month_start = 1, int season_month_end = 12,
                double start_time = 2018, double end_time = 2018,
                std::string dispersal_kern = "cauchy", double percent_short_distance_dispersal = 0.0,
                double long_distance_scale = 0.0,
                std::string wind_dir = "NONE", double kappa = 0
)
{
  
  std::vector<std::tuple<int, int>> outside_dispersers;
  DispersalKernel dispersal_kernel = radial_type_from_string(dispersal_kern);
  pops::Date dd_start(start_time, 01, 01);
  pops::Date dd_end(end_time, 12, 31);
  Direction wind_direction = direction_enum_from_string(wind_dir);
  Season season(season_month_start,season_month_end);
  pops::Date dd_current(dd_start);
  Simulation<IntegerMatrix, NumericMatrix> simulation(random_seed, infected, ew_res, ns_res);
  int counter = 0;
  
  std::vector<IntegerMatrix> infected_vector;
  std::vector<IntegerMatrix> susceptible_vector;
  std::vector<int> simulated_weeks;
  
  for (unsigned current_time_step = 0; ; current_time_step++, time_step == "month" ? dd_current.increased_by_month() : dd_current.increased_by_week()) {
      
      if (dd_current.year() > dd_end.year()) {
        break;
      }
      
    if (all_infected(susceptible)) {
      Rcerr << "At timestep " << dd_current << " all suspectible hosts are infected!" << std::endl;
      infected_vector.push_back(Rcpp::clone(infected));
      susceptible_vector.push_back(Rcpp::clone(susceptible));
      break;
    }
    
    if (use_lethal_temperature && dd_current.month() == lethal_temperature_month && dd_current.year() <= dd_end.year()) {
      unsigned simulation_year = dd_current.year() - dd_start.year();
      if (simulation_year >= temperature.size()){
        Rcerr << "Not enough years of temperature data" << std::endl;
      }
      simulation.remove(infected, susceptible, temperature[simulation_year], lethal_temperature);
    }
    
    if (season.month_in_season(dd_current.month())) {
      counter += 1;
      simulated_weeks.push_back(current_time_step);
      
      if (current_time_step >= weather_coefficient.size()  && weather == TRUE) {
        Rcerr << "Not enough time steps of weather coefficient data" << std::endl;
      }

      simulation.generate(infected, weather, weather_coefficient[current_time_step], reproductive_rate);
      simulation.disperse(susceptible, infected,
                          mortality_tracker, total_plants,
                          outside_dispersers, 
                          weather, weather_coefficient[current_time_step], 
                          dispersal_kernel, short_distance_scale,
                          percent_short_distance_dispersal, long_distance_scale,
                          wind_direction, kappa);
    
    }
    
    if ((time_step == "month" ? dd_current.is_last_month_of_year() : dd_current.is_last_week_of_year())) {
      infected_vector.push_back(Rcpp::clone(infected));
      susceptible_vector.push_back(Rcpp::clone(susceptible));
    }
    
    if (dd_current >= dd_end) {
      break;
    }
  
  }

  return List::create(
    _["infected_vector"] = infected_vector,
    _["susceptible_vector"] = susceptible_vector
  );
  
}
