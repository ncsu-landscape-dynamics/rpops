#include "date.hpp"
#include "kernel.hpp"
#include "kernel_types.hpp"
#include "model.hpp"
#include "natural_anthropogenic_kernel.hpp"
#include "radial_kernel.hpp"
#include "raster.hpp"
#include "scheduling.hpp"
#include "spread_rate.hpp"
#include "statistics.hpp"
#include "switch_kernel.hpp"
#include "treatments.hpp"
#include "uniform_kernel.hpp"
#include "quarantine.hpp"
#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::to_string;

using namespace Rcpp;
using namespace pops;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List pops_model_cpp(
    int random_seed,
    bool multiple_random_seeds,
    std::vector<unsigned> random_seeds,
    double lethal_temperature,
    int lethal_temperature_month,
    IntegerMatrix infected,
    IntegerMatrix total_exposed,
    std::vector<IntegerMatrix> exposed,
    IntegerMatrix susceptible,
    IntegerMatrix total_populations,
    IntegerMatrix total_hosts,
    std::vector<IntegerMatrix> mortality_tracker,
    IntegerMatrix mortality,
    IntegerMatrix quarantine_areas,
    std::string quarantine_directions,
    std::vector<NumericMatrix> treatment_maps,
    std::vector<std::string> treatment_dates,
    std::vector<int> pesticide_duration,
    IntegerMatrix resistant,
    std::vector<std::vector<int>> movements,
    std::vector<std::string> movements_dates,
    std::vector<NumericMatrix> temperature,
    std::vector<NumericMatrix> survival_rates,
    std::vector<NumericMatrix> weather_coefficient,
    List bbox,
    List res,
    List rows_cols,
    double reproductive_rate,
    std::vector<std::vector<int>> spatial_indices,
    List season_month_start_end,
    List frequency_config,
    List bool_config,
    double mortality_rate = 0.0,
    int mortality_time_lag = 2,
    std::string start_date = "2018-01-01",
    std::string end_date = "2018-12-31",
    std::string treatment_method = "ratio",
    std::string natural_kernel_type = "cauchy",
    std::string anthropogenic_kernel_type = "cauchy",
    double percent_natural_dispersal = 0.0,
    double natural_distance_scale = 21,
    double anthropogenic_distance_scale = 0.0,
    std::string natural_dir = "NONE",
    double natural_kappa = 0,
    std::string anthropogenic_dir = "NONE",
    double anthropogenic_kappa = 0,
    Nullable<List> frequencies_n_config = R_NilValue,
    std::string model_type_ = "SI",
    int latency_period = 0,
    double establishment_probability = 0,
    double dispersal_percentage = 0.99,
    int survival_rate_month = 0,
    int survival_rate_day = 0,
    Nullable<List> overpopulation_config = R_NilValue,
    Nullable<List> network_config = R_NilValue,
    Nullable<List> network_data_config = R_NilValue,
    int weather_size = 0,
    std::string weather_type = "deterministic",
    double dispersers_to_soils_percentage = 0)
{
    Config config;
    config.random_seed = random_seed;
    config.rows = rows_cols["num_rows"];
    config.cols = rows_cols["num_cols"];
    config.ew_res = res["ew_res"];
    config.ns_res = res["ns_res"];

    config.generate_stochasticity = bool_config["generate_stochasticity"];
    config.establishment_stochasticity = bool_config["establishment_stochasticity"];
    config.movement_stochasticity = bool_config["movement_stochasticity"];
    config.dispersal_stochasticity = bool_config["dispersal_stochasticity"];
    config.establishment_probability = establishment_probability;

    config.use_lethal_temperature = bool_config["use_lethal_temperature"];
    config.lethal_temperature = lethal_temperature;
    config.lethal_temperature_month = lethal_temperature_month;
    config.weather = bool_config["weather"];
    config.weather_size = weather_size;
    config.weather_type = weather_type;

    config.reproductive_rate = reproductive_rate;
    config.model_type = model_type_;
    config.latency_period_steps = latency_period;
    config.natural_kernel_type = natural_kernel_type;
    config.natural_scale = natural_distance_scale;
    config.natural_direction = natural_dir;
    config.natural_kappa = natural_kappa;

    config.use_survival_rate = bool_config["use_survival_rate"];
    config.survival_rate_day = survival_rate_day;
    config.survival_rate_month = survival_rate_month;

    config.use_anthropogenic_kernel = bool_config["use_anthropogenic_kernel"];
    config.percent_natural_dispersal = percent_natural_dispersal;
    config.anthro_kernel_type = anthropogenic_kernel_type;
    config.anthro_scale = anthropogenic_distance_scale;
    config.anthro_direction = anthropogenic_dir;
    config.anthro_kappa = anthropogenic_kappa;

    std::string time_step = frequency_config["time_step"];
    std::string output_frequency = frequency_config["output_frequency"];
    std::string quarantine_frequency = frequency_config["quarantine_frequency"];
    std::string spreadrate_frequency = frequency_config["spreadrate_frequency"];
    std::string mortality_frequency = frequency_config["mortality_frequency"];
    // use_treatment set later
    config.use_mortality = bool_config["mortality_on"];
    config.mortality_rate = mortality_rate;
    config.mortality_time_lag = mortality_time_lag;
    if (output_frequency == "time_step") {
        output_frequency = time_step;
    }
    config.use_movements = bool_config["use_movements"];
    // movement_schedule used later

    config.dispersal_percentage = dispersal_percentage;

    if (frequencies_n_config.isNotNull()) {
        List freq_n_config(frequencies_n_config);
        config.output_frequency_n = freq_n_config["output_frequency_n"];
        config.quarantine_frequency_n = freq_n_config["quarantine_frequency_n"];
        config.spreadrate_frequency_n = freq_n_config["spreadrate_frequency_n"];
        config.mortality_frequency_n = freq_n_config["mortality_frequency_n"];
    }

    config.output_frequency = output_frequency;
    config.quarantine_frequency = quarantine_frequency;
    config.use_quarantine = bool_config["use_quarantine"];
    config.quarantine_directions = quarantine_directions;
    config.spreadrate_frequency = spreadrate_frequency;
    config.mortality_frequency = mortality_frequency;
    config.use_spreadrates = bool_config["use_spreadrates"];
    config.use_overpopulation_movements = bool_config["use_overpopulation_movements"];
    if (config.use_overpopulation_movements && overpopulation_config.isNotNull()) {
        List over_config(overpopulation_config);
        config.overpopulation_percentage = over_config["overpopulation_percentage"];
        config.leaving_percentage = over_config["leaving_percentage"];
        config.leaving_scale_coefficient = over_config["leaving_scale_coefficient"];
    }

    std::vector<std::tuple<int, int>> outside_dispersers;
    TreatmentApplication treatment_application =
        treatment_app_enum_from_string(treatment_method);
    config.set_date_start(start_date);
    config.set_date_end(end_date);
    config.set_step_unit(time_step);
    config.set_step_num_units(1);
    int start_month = season_month_start_end["start_month"];
    int end_month = season_month_start_end["end_month"];
    config.set_season_start_end_month(start_month, end_month);
    config.dispersers_to_soils_percentage = dispersers_to_soils_percentage;

    std::vector<std::array<double, 4>> spread_rates_vector;
    std::tuple<double, double, double, double> spread_rates;
    IntegerMatrix total_dispersers(config.rows, config.cols);
    IntegerMatrix established_dispersers(config.rows, config.cols);

    int num_infected;
    std::vector<int> number_infected;
    double area_infect;
    std::vector<double> area_infected;

    std::vector<IntegerMatrix> infected_vector;
    std::vector<IntegerMatrix> susceptible_vector;
    std::vector<IntegerMatrix> mortality_vector;
    std::vector<IntegerMatrix> resistant_vector;
    std::vector<IntegerMatrix> total_populations_vector;
    std::vector<IntegerMatrix> total_exposed_vector;
    std::vector<IntegerMatrix> dispersers_vector;
    std::vector<IntegerMatrix> exposed_v;
    std::vector<std::vector<IntegerMatrix>> exposed_vector;

    config.create_schedules();

    Treatments<IntegerMatrix, NumericMatrix> treatments(config.scheduler());
    for (unsigned t = 0; t < treatment_maps.size(); t++) {
        treatments.add_treatment(
            treatment_maps[t],
            pops::Date(treatment_dates[t]),
            pesticide_duration[t],
            treatment_application);
        config.use_treatments = true;
    }

    if (config.use_lethal_temperature) {
        if (config.num_lethal() > temperature.size()) {
            Rcerr << "Not enough years of temperature data" << std::endl;
        }
    }

    unsigned count_weather = get_number_of_scheduled_actions(config.spread_schedule());
    if (config.weather && count_weather > weather_coefficient.size()) {
        Rcerr << "Not enough indices of weather coefficient data" << std::endl;
    }

    unsigned spread_rate_outputs;
    if (config.use_spreadrates) {
        spread_rate_outputs = config.rate_num_steps();
    }
    else {
        spread_rate_outputs = 0;
    }
    SpreadRate<IntegerMatrix> spreadrate(
        infected, config.ew_res, config.ns_res, spread_rate_outputs, spatial_indices);
    unsigned move_scheduled;
    if (config.use_movements) {
        for (unsigned move = 0; move < movements_dates.size(); ++move) {
            pops::Date movement_date(movements_dates[move]);
            move_scheduled =
                unsigned(config.scheduler().schedule_action_date(movement_date));
            config.movement_schedule.push_back(move_scheduled);
        }
    }

    unsigned quarantine_outputs;
    if (config.use_quarantine) {
        quarantine_outputs = config.quarantine_num_steps();
    }
    else {
        quarantine_outputs = 0;
    }

    QuarantineEscape<IntegerMatrix> quarantine(
        quarantine_areas,
        config.ew_res,
        config.ns_res,
        quarantine_outputs,
        config.quarantine_directions);
    bool quarantine_escape;
    std::vector<bool> quarantine_escapes;
    int escape_dist;
    std::vector<int> escape_dists;
    Direction escape_direction;
    std::vector<std::string> escape_directions;

    std::unique_ptr<Network<int>> network{nullptr};
    if (network_config.isNotNull() && network_data_config.isNotNull()) {
        // The best place for bbox handling would be with rows, cols, and
        // resolution, but since it is required only for network, it is here.
        config.bbox.north = bbox["north"];
        config.bbox.south = bbox["south"];
        config.bbox.east = bbox["east"];
        config.bbox.west = bbox["west"];
        List net_config(network_config);
        config.network_min_distance = net_config["network_min_distance"];
        config.network_max_distance = net_config["network_max_distance"];
        std::string network_movement = net_config["network_movement"];
        config.network_movement = network_movement;
        network.reset(new Network<int>(config.bbox, config.ew_res, config.ns_res));
        List net_data_config(network_data_config);
        std::ifstream network_stream{
            Rcpp::as<std::string>(net_data_config["network_filename"])};
        network->load(network_stream);
    }

    ModelType mt = model_type_from_string(config.model_type);
    Simulation<IntegerMatrix, NumericMatrix> simulation(
        config.rows, config.cols, mt, config.latency_period_steps);

    Model<IntegerMatrix, NumericMatrix, int> model(config);
    for (unsigned current_index = 0; current_index < config.scheduler().get_num_steps();
         ++current_index) {

        IntegerMatrix dispersers(config.rows, config.cols);
        model.run_step(
            current_index,
            infected,
            susceptible,
            total_populations,
            total_hosts,
            dispersers,
            established_dispersers,
            total_exposed,
            exposed,
            mortality_tracker,
            mortality,
            temperature,
            survival_rates,
            treatments,
            resistant,
            outside_dispersers,
            spreadrate,
            quarantine,
            quarantine_areas,
            movements,
            network ? *network : Network<int>::null_network(),
            spatial_indices);

        // keeps track of cumulative dispersers or propagules from a site.
        if (config.spread_schedule()[current_index]) {
            total_dispersers += dispersers;
        }

        if (config.use_mortality && config.mortality_schedule()[current_index]) {
            mortality_vector.push_back(Rcpp::clone(mortality));
        }

        if (config.output_schedule()[current_index]) {
            infected_vector.push_back(Rcpp::clone(infected));
            susceptible_vector.push_back(Rcpp::clone(susceptible));
            resistant_vector.push_back(Rcpp::clone(resistant));
            total_populations_vector.push_back(Rcpp::clone(total_populations));
            total_exposed_vector.push_back(Rcpp::clone(total_exposed));
            dispersers_vector.push_back(Rcpp::clone(total_dispersers));

            if (config.model_type == "SEI") {
                exposed_v.clear();

                for (unsigned e = 0; e < exposed.size(); e++) {
                    exposed_v.push_back(Rcpp::clone(exposed[e]));
                }
            }
            else {
                exposed_v = exposed;
            }

            // exposed_v = exposed;
            exposed_vector.push_back(exposed_v);

            num_infected = sum_of_infected(infected, spatial_indices);
            number_infected.push_back(num_infected);
            area_infect = area_of_infected(
                infected, config.ew_res, config.ns_res, spatial_indices);
            area_infected.push_back(area_infect);
            total_dispersers(config.rows, config.cols);
        }

        // update spread rate outputs if they are used and scheduled for that time step
        if (config.use_spreadrates && config.spread_rate_schedule()[current_index]) {
            unsigned simulation_step = simulation_step_to_action_step(
                config.spread_rate_schedule(), current_index);
            spread_rates = spreadrate.step_rate(simulation_step);
            auto sr = to_array(spread_rates);
            spread_rates_vector.push_back(sr);
        }

        // update quarantine outputs if they are used and scheduled for that time step
        if (config.use_quarantine && config.quarantine_schedule()[current_index]) {
            unsigned quarantine_step = simulation_step_to_action_step(
                config.quarantine_schedule(), current_index);
            quarantine_escape = quarantine.escaped(quarantine_step);
            escape_dist = quarantine.distance(quarantine_step);
            escape_direction = quarantine.direction(quarantine_step);
            quarantine_escapes.push_back(quarantine_escape);
            escape_dists.push_back(escape_dist);
            escape_directions.push_back(quarantine_enum_to_string(escape_direction));
        }
    }

    return List::create(
        _["infected"] = infected_vector,
        _["exposed"] = exposed_vector,
        _["susceptible"] = susceptible_vector,
        _["resistant"] = resistant_vector,
        _["mortality"] = mortality_vector,
        _["rates"] = spread_rates_vector,
        _["number_infected"] = number_infected,
        _["area_infected"] = area_infected,
        _["total_populations"] = total_populations_vector,
        _["total_exposed"] = total_exposed_vector,
        _["propogules"] = dispersers_vector,
        _["quarantine_escape"] = quarantine_escapes,
        _["quarantine_escape_distance"] = escape_dists,
        _["quarantine_escape_directions"] = escape_directions,
        _["spatial_indices"] = spatial_indices);
}
