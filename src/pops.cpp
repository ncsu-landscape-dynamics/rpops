#include <Rcpp.h>
#include "simulation.hpp"
#include "raster.hpp"
#include "date.hpp"
#include "treatments.hpp"
#include "kernel.hpp"
#include "kernel_types.hpp"
#include "radial_kernel.hpp"
#include "short_long_kernel.hpp"
#include "spread_rate.hpp"
#include "statistics.hpp"
#include "switch_kernel.hpp"
#include "uniform_kernel.hpp"
#include "scheduling.hpp"
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
  
  std::vector<std::tuple<int, int>> outside_dispersers;
  DispersalKernelType natural_dispersal_kernel_type = kernel_type_from_string(natural_kernel_type);
  DispersalKernelType anthropogenic_dispersal_kernel_type = kernel_type_from_string(anthropogenic_kernel_type);
  TreatmentApplication treatment_application = treatment_app_enum_from_string(treatment_method);
  ModelType model_type = model_type_from_string(model_type_);
  pops::Date dd_start(start_date);
  pops::Date dd_end(end_date);
  Direction natural_direction = direction_from_string(natural_dir);
  Direction anthropogenic_direction = direction_from_string(anthropogenic_dir);
  Season season(season_month_start, season_month_end);
  pops::Date dd_current(dd_start);
  Simulation<IntegerMatrix, NumericMatrix> simulation(random_seed, num_rows, num_cols, model_type, latency_period);
  RadialDispersalKernel natural_radial_dispersal_kernel(ew_res, ns_res, natural_dispersal_kernel_type, natural_distance_scale, natural_direction, natural_kappa);
  RadialDispersalKernel anthropogenic_radial_dispersal_kernel(ew_res, ns_res, anthropogenic_dispersal_kernel_type, anthropogenic_distance_scale, anthropogenic_direction, anthropogenic_kappa);
  UniformDispersalKernel uniform_kernel(num_rows, num_cols);
  SwitchDispersalKernel natural_dispersal_kernel(natural_dispersal_kernel_type, natural_radial_dispersal_kernel, uniform_kernel);
  SwitchDispersalKernel anthropogenic_dispersal_kernel(anthropogenic_dispersal_kernel_type, anthropogenic_radial_dispersal_kernel, uniform_kernel);
  DispersalKernel kernel(natural_dispersal_kernel, anthropogenic_dispersal_kernel, use_anthropogenic_kernel, percent_natural_dispersal);
  std::vector<std::array<double,4>> spread_rates_vector;
  std::tuple<double,double,double,double> spread_rates;
  IntegerMatrix total_dispersers(num_rows, num_cols);
  
  int num_infected;
  std::vector<int> number_infected;
  double area_infect;
  std::vector<double> area_infected;
  
  if (output_frequency == "time_step") {
    output_frequency = time_step;
  }
  
  int first_mortality_year = mortality_time_lag;
  
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
  bool use_treatments = false;
  for (unsigned t = 0; t < treatment_maps.size(); t++) {
    treatments.add_treatment(treatment_maps[t], pops::Date(treatment_dates[t]), pesticide_duration[t], treatment_application);
    use_treatments = true;
  }
  
  unsigned count_lethal = get_number_of_scheduled_actions(lethality_schedule);
    if (use_lethal_temperature && count_lethal > temperature.size()){
      Rcerr << "Not enough years of temperature data" << std::endl;
    }
    
  unsigned count_weather = get_number_of_scheduled_actions(spread_schedule);
      if (weather && count_weather > weather_coefficient.size()) {
        Rcerr << "Not enough indices of weather coefficient data" << std::endl;
      }
  
  unsigned mortality_simulation_year;
  
  unsigned spread_rate_outputs = get_number_of_scheduled_actions(spread_rate_schedule);
  SpreadRate<IntegerMatrix> spreadrate(infected, ew_res, ns_res, spread_rate_outputs);
  // Define movement schedule  
  unsigned last_index = 0;
  unsigned move_scheduled;
  std::vector<unsigned> movement_schedule;
  if (use_movements) {
    for (unsigned move = 0; move < movements_dates.size(); ++move){
        pops::Date movement_date(movements_dates[move]);
        move_scheduled = unsigned(scheduler.schedule_action_date(movement_date));
        movement_schedule.push_back(move_scheduled);
    }
  }

  for (unsigned current_index = 0; current_index < scheduler.get_num_steps(); ++current_index) {
      
      // if (all_infected(susceptible)) {
      //   Rcerr << "All suspectible hosts are infected!" << std::endl;
      //   infected_vector.push_back(Rcpp::clone(infected));
      //   susceptible_vector.push_back(Rcpp::clone(susceptible));
      //   resistant_vector.push_back(Rcpp::clone(resistant));
      //   break;
      // }
    
      if (lethality_schedule[current_index]  && use_lethal_temperature) {
        int lethal_step = simulation_step_to_action_step(lethality_schedule, current_index);
        simulation.remove(infected, susceptible, temperature[lethal_step], lethal_temperature);
      }
      
      if (use_treatments) {
        bool managed = treatments.manage(current_index, infected, susceptible, resistant);
        
        if (mortality_on && managed) {
          for (unsigned l = 0; l < mortality_tracker_vector.size(); l++) {
            treatments.manage_mortality(current_index, mortality_tracker_vector[l]);
          }
        }
        
        if (model_type == ModelType::SusceptibleExposedInfected && managed) {
          for (unsigned e = 0; e < exposed.size(); e++){
            treatments.manage_mortality(current_index, exposed[e]);
          }
        }
      }

      if (spread_schedule[current_index]) {
        IntegerMatrix dispersers(num_rows, num_cols);
        simulation.generate(dispersers, infected, weather, weather_coefficient[current_index], reproductive_rate);
        total_dispersers += dispersers;
        // dispersers_vector.push_back(Rcpp::clone(dispersers));
        simulation.disperse_and_infect(current_index, dispersers, susceptible, exposed, infected, mortality_tracker, total_plants,
                            outside_dispersers, weather, weather_coefficient[current_index], kernel);
        if (use_movements) {
          last_index = simulation.movement(infected, susceptible, mortality_tracker, total_plants, current_index, last_index, movements, movement_schedule);
        }
      }
      
      if (mortality_on && mortality_schedule[current_index]) {
        mortality_tracker_vector.push_back(Rcpp::clone(mortality_tracker));
        std::fill(mortality_tracker.begin(), mortality_tracker.end(), 0);
        mortality_simulation_year = simulation_step_to_action_step(mortality_schedule, current_index);
        simulation.mortality(infected, mortality_rate, mortality_simulation_year, first_mortality_year, mortality, mortality_tracker_vector);
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
        area_infect = area_of_infected(infected, ew_res, ns_res);
        area_infected.push_back(area_infect);
        // reinitialize total dispersers so each output isn't an accumulation of the previous output
        total_dispersers(num_rows, num_cols);
      }
      
      if (spread_rate_schedule[current_index]) {
        unsigned simulation_year = simulation_step_to_action_step(spread_rate_schedule, current_index);
        spreadrate.compute_yearly_spread_rate(infected, simulation_year);
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
