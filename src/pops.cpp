#include <Rcpp.h>
// change to use pops::model
#include "model.hpp"
//#include "simulation.hpp"
#include "raster.hpp"
#include "date.hpp"
//#include "treatments.hpp"
//#include "kernel.hpp"
#include "kernel_types.hpp"
#include "radial_kernel.hpp"
#include "short_long_kernel.hpp"
//#include "spread_rate.hpp"
#include "statistics.hpp"
#include "switch_kernel.hpp"
#include "uniform_kernel.hpp"
//#include "scheduling.hpp"
#include <iostream>
#include <vector>
#include <tuple>
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

template<int... Indices>
struct indices {
  using next = indices<Indices..., sizeof...(Indices)>;
};

template<int Size>
struct build_indices {
  using type = typename build_indices<Size - 1>::type::next;
};

template<>
struct build_indices<0> {
  using type = indices<>;
};

template<typename T>
using Bare = typename std::remove_cv<typename std::remove_reference<T>::type>::type;

template<typename Tuple>
constexpr
  typename build_indices<std::tuple_size<Bare<Tuple>>::value>::type
  make_indices()
  { return {}; }

template<typename Tuple, int... Indices>
std::array<
  typename std::tuple_element<0, Bare<Tuple>>::type,
  std::tuple_size<Bare<Tuple>>::value
>
to_array(Tuple&& tuple, indices<Indices...>)
{
  using std::get;
  return {{ get<Indices>(std::forward<Tuple>(tuple))... }};
}

template<typename Tuple>
auto to_array(Tuple&& tuple)
  -> decltype( to_array(std::declval<Tuple>(), make_indices<Tuple>()) )
  {
    return to_array(std::forward<Tuple>(tuple), make_indices<Tuple>());
  }

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
// TODO add in reduced stochasticity options?
List pops_model(int random_seed, 
                bool use_lethal_temperature, double lethal_temperature, int lethal_temperature_month,
                IntegerMatrix infected,
                std::vector<IntegerMatrix> exposed,
                IntegerMatrix susceptible,
                IntegerMatrix total_plants,
                bool mortality_on,
                IntegerMatrix mortality_tracker,
                IntegerMatrix mortality,
                std::vector<NumericMatrix> treatment_maps,
                std::vector<std::string> treatment_dates,
                std::vector<int> pesticide_duration,
                IntegerMatrix resistant,
                bool use_movements, std::vector<std::vector<int>> movements,
                std::vector<std::string> movements_dates,
                bool weather,
                std::vector<NumericMatrix> temperature,
                std::vector<NumericMatrix> weather_coefficient,
                double ew_res, double ns_res, int num_rows, int num_cols,
                std::string time_step, double reproductive_rate,
                double mortality_rate = 0.0, int mortality_time_lag = 2,
                int season_month_start = 1, int season_month_end = 12,
                std::string start_date = "2018-01-01", std::string end_date = "2018-12-31",
                std::string treatment_method = "ratio",
                std::string natural_kernel_type = "cauchy", 
                std::string anthropogenic_kernel_type = "cauchy", 
                bool use_anthropogenic_kernel = false, double percent_natural_dispersal = 0.0,
                double natural_distance_scale = 21, double anthropogenic_distance_scale = 0.0, 
                std::string natural_dir = "NONE", double natural_kappa = 0,
                std::string anthropogenic_dir = "NONE", double anthropogenic_kappa = 0,
                std::string output_frequency = "year", std::string model_type_ = "SI",
                int latency_period = 0
)
{
  // TODO no establishment_probability, establishment_stochasticity, or establishment_probability
  // steps and use_treatments set later (set to 0 and false for now)
  Config config = {random_seed, num_rows, num_cols, ew_res, ns_res, 0, true, true, 0.5, 
  use_lethal_temperature, lethal_temperature, weather, reproductive_rate, model_type_, latency_period,
  natural_kernel_type, natural_distance_scale, natural_dir, natural_kappa, use_anthropogenic_kernel,
  percent_natural_dispersal, anthropogenic_kernel_type, anthropogenic_distance_scale, anthropogenic_dir,
  anthropogenic_kappa, false, mortality_on, mortality_rate, mortality_time_lag};
  
  std::vector<std::tuple<int, int>> outside_dispersers;
  TreatmentApplication treatment_application = treatment_app_enum_from_string(treatment_method);
  pops::Date dd_start(start_date);
  pops::Date dd_end(end_date);
  Season season(season_month_start, season_month_end);
  pops::Date dd_current(dd_start);

  std::vector<std::array<double,4>> spread_rates_vector;
  std::tuple<double,double,double,double> spread_rates;
  IntegerMatrix total_dispersers(config.rows, config.cols);

  int num_infected;
  std::vector<int> number_infected;
  double area_infect;
  std::vector<double> area_infected;

  if (output_frequency == "time_step") {
    output_frequency = time_step;
  }

  std::vector<IntegerMatrix> infected_vector;
  std::vector<IntegerMatrix> susceptible_vector;
  std::vector<IntegerMatrix> mortality_tracker_vector;
  std::vector<IntegerMatrix> mortality_vector;
  std::vector<IntegerMatrix> resistant_vector;
  std::vector<IntegerMatrix> total_host_vector;
  std::vector<IntegerMatrix> dispersers_vector;

  StepUnit step_unit = step_unit_enum_from_string(time_step);

  // Define simulation time step
  Scheduler scheduler(dd_start, dd_end, step_unit, 1);
  // set config.steps
  config.steps = scheduler.get_num_steps();
  // Define spread schedule
  std::vector<bool> spread_schedule = scheduler.schedule_spread(season);
  // Define spread rate schedule
  std::vector<bool> spread_rate_schedule = scheduler.schedule_action_end_of_year();
  // Define mortality schedule
  std::vector<bool> mortality_schedule = scheduler.schedule_action_end_of_year();
  // Define lethality schedule
  std::vector<bool> lethality_schedule = scheduler.schedule_action_yearly(lethal_temperature_month, 1);
  // Define output schedule
  std::vector<bool> output_schedule;
  if (output_frequency == "year") {
    output_schedule = scheduler.schedule_action_end_of_year();
  } else if (output_frequency == "month") {
    output_schedule = scheduler.schedule_action_monthly();
  } else if (output_frequency == "week") {
    if (time_step == "day") {
      output_schedule = scheduler.schedule_action_nsteps(7);
    } else if (time_step == "week") {
      output_schedule = scheduler.schedule_action_nsteps(1);
    }
  } else if (output_frequency == "day") {
    output_schedule = scheduler.schedule_action_nsteps(1);
  }

  Treatments<IntegerMatrix, NumericMatrix> treatments(scheduler);
  //config.use_treatments = false; set to false in constructor
  for (unsigned t = 0; t < treatment_maps.size(); t++) {
    treatments.add_treatment(treatment_maps[t], pops::Date(treatment_dates[t]), pesticide_duration[t], treatment_application);
    config.use_treatments = true;
  }

  unsigned count_lethal = get_number_of_scheduled_actions(lethality_schedule);
  if (config.use_lethal_temperature && count_lethal > temperature.size()) {
    Rcerr << "Not enough years of temperature data" << std::endl;
  }

  unsigned count_weather = get_number_of_scheduled_actions(spread_schedule);
  if (config.weather && count_weather > weather_coefficient.size()) {
    Rcerr << "Not enough indices of weather coefficient data" << std::endl;
  }

  unsigned spread_rate_outputs = get_number_of_scheduled_actions(spread_rate_schedule);
  SpreadRate<IntegerMatrix> spreadrate(infected, config.ew_res, config.ns_res, spread_rate_outputs);
  // Define movement schedule  
  unsigned last_index = 0;
  unsigned move_scheduled;
  std::vector<unsigned> movement_schedule;
  if (use_movements) {
    for (unsigned move = 0; move < movements_dates.size(); ++move) {
        pops::Date movement_date(movements_dates[move]);
        move_scheduled = unsigned(scheduler.schedule_action_date(movement_date));
        movement_schedule.push_back(move_scheduled);
    }
  }

  Model<IntegerMatrix, NumericMatrix, NumericMatrix> model(config);
  IntegerMatrix dispersers;
  for (unsigned current_index = 0; current_index < config.steps; ++current_index) {

    // if (all_infected(susceptible)) {
    //   Rcerr << "All suspectible hosts are infected!" << std::endl;
    //   infected_vector.push_back(Rcpp::clone(infected));
    //   susceptible_vector.push_back(Rcpp::clone(susceptible));
    //   resistant_vector.push_back(Rcpp::clone(resistant));
    //   break;
    // }
    
    if (spread_schedule[current_index]) {
      dispersers(config.rows, config.cols);
    }
    if (config.use_mortality && mortality_schedule[current_index]) {
      mortality_tracker_vector.push_back(Rcpp::clone(mortality_tracker));
      std::fill(mortality_tracker.begin(), mortality_tracker.end(), 0);
    }
    // TODO pass in current_index as step to run_step?
    model.run_step(current_index, spread_schedule, mortality_schedule, lethality_schedule,
    spread_rate_schedule, count_weather, infected, susceptible, total_plants, dispersers, exposed, mortality_tracker_vector,
    mortality, temperature, weather_coefficient, treatments, resistant, outside_dispersers,
    spreadrate);

   if (spread_schedule[current_index]) {
      total_dispersers += dispersers;
      if (use_movements) {
        last_index = simulation.movement(infected, susceptible, mortality_tracker, total_plants, current_index, last_index, movements, movement_schedule);
      }
    }

    if (config.use_mortality && mortality_schedule[current_index]) {
      mortality_vector.push_back(Rcpp::clone(mortality));
    }

    if (output_schedule[current_index]) {
      infected_vector.push_back(Rcpp::clone(infected));
      susceptible_vector.push_back(Rcpp::clone(susceptible));
      resistant_vector.push_back(Rcpp::clone(resistant));
      total_host_vector.push_back(Rcpp::clone(total_plants));
      dispersers_vector.push_back(Rcpp::clone(total_dispersers));
      // exposed_vector = Rcpp::clone(exposed);

      num_infected = sum_of_infected(infected);
      number_infected.push_back(num_infected);
      area_infect = area_of_infected(infected, config.ew_res, config.ns_res);
      area_infected.push_back(area_infect);
      // reinitialize total dispersers so each output isn't an accumulation of the previous output
      total_dispersers(config.rows, config.cols);
    }

    if (spread_rate_schedule[current_index]) {
      // leaving this line here (also in model.hpp) because without it simulation_year is unknown
       unsigned simulation_year = simulation_step_to_action_step(spread_rate_schedule, current_index);
      spread_rates = spreadrate.yearly_rate(simulation_year);
      auto sr = to_array(spread_rates);
      spread_rates_vector.push_back(sr);
    }
  }

  return List::create(
    _["infected"] = infected_vector,
    _["exposed"] = exposed,
    _["susceptible"] = susceptible_vector,
    _["resistant"] = resistant_vector,
    _["mortality"] = mortality_vector,
    _["rates"] = spread_rates_vector,
    _["number_infected"] = number_infected,
    _["area_infected"] = area_infected,
    _["total_hosts"] = total_host_vector,
    _["propogules"] = dispersers_vector
  );
  
}
