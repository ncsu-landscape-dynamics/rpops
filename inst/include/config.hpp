/*
 * Tests for the PoPS Config class.
 *
 * Copyright (C) 2020-2022 by the authors.
 *
 * Authors: Vaclav Petras <wenzeslaus gmail com>
 *
 * This file is part of PoPS.

 * PoPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * PoPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with PoPS. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef POPS_CONFIG_HPP
#define POPS_CONFIG_HPP

#include "scheduling.hpp"
#include "utils.hpp"

#include <vector>
#include <iostream>
#include <regex>
#include <string>
#include <sstream>

namespace pops {

/**
 * Read key-value pairs from text into a map
 *
 * Text can be, e.g., comma-separated pairs of key and value where
 * key and value are separated by equal (key=value,key2=value2) or
 * YAML-style lines of `key: value`.
 *
 * Both separators are a single character and need to be different from
 * each other. Common separators such as comma, semicolor, or colon will
 * work. Special characters for regular expressions such as bracket or asterisk
 * will confuse the parser.
 *
 * @param stream Text as stream
 * @param record_separator Character which separates individual records
 * @param key_value_separator Character which separates the key and value
 * @param conversion Function to convert string to value (use lambda)
 *
 * @see Other overloads.
 */
template<typename Value, typename Conversion>
std::map<std::string, Value> read_key_value_pairs(
    std::istream& stream,
    char record_separator,
    char key_value_separator,
    Conversion conversion)
{
    std::map<std::string, Value> config;
    std::string line;
    std::string expression_string(R"(\s*([^= ]+)[\s]*=[\s]*([^ ]+))");
    std::replace(
        expression_string.begin(), expression_string.end(), '=', key_value_separator);
    std::regex expression(expression_string);
    while (std::getline(stream, line, record_separator)) {
        std::smatch match;
        if (regex_search(line, match, expression) && match.size() == 3) {
            std::string value(match[2]);
            config[match[1]] = conversion(value.c_str());
        }
        else {
            throw std::invalid_argument(std::string("Incorrect format of: ") + line);
        }
    }
    return config;
}

/**
 * Read key-value pairs from text into a map
 *
 * @see Other overloads.
 */
template<typename Value, typename Conversion>
std::map<std::string, Value> read_key_value_pairs(
    const std::string& text,
    char record_separator,
    char key_value_separator,
    Conversion conversion)
{
    std::istringstream stream(text);
    return read_key_value_pairs<Value>(
        stream, record_separator, key_value_separator, conversion);
}

/**
 * Read key-value pairs from text into a map
 *
 * @see Other overloads.
 */
template<typename Value, typename Conversion>
std::map<std::string, Value> read_key_value_pairs(
    const char* text,
    char record_separator,
    char key_value_separator,
    Conversion conversion)
{
    std::string std_text(text);
    return read_key_value_pairs<Value>(
        std_text, record_separator, key_value_separator, conversion);
}

/** Configuration for Model */
class Config
{
public:
    /**
     * @brief Row of a table with pest-host-use data
     *
     * One row is stored for each host.
     */
    struct PestHostUseTableDataRow
    {
        double susceptibility;  ///< Susceptibility for a given host
        double mortality_rate;  ///< Mortality rate for a given host
        double mortality_time_lag;  ///< Time lag of mortality in mortality steps
    };
    /**
     * @brief Row of a table with competency data
     */
    struct CompetencyTableDataRow
    {
        std::vector<bool> presence_absence;
        double competency;
    };

    // Seed
    int random_seed{0};
    bool multiple_random_seeds{false};
    std::map<std::string, unsigned> random_seeds;
    // Size
    int rows{0};
    int cols{0};
    double ew_res{0};
    double ns_res{0};
    BBox<double> bbox;
    // Reduced stochasticity
    bool generate_stochasticity{true};
    bool establishment_stochasticity{true};
    bool movement_stochasticity{true};
    bool dispersal_stochasticity{true};
    double establishment_probability{0};
    // Temperature
    bool use_lethal_temperature{false};
    double lethal_temperature{-273.15};  // 0 K
    int lethal_temperature_month{0};
    bool weather{false};
    int weather_size{0};  ///< Number of weather steps (size of weather time series)
    std::string weather_type;  ///< probabilistic, deterministic
    double reproductive_rate{0};
    // survival rate
    bool use_survival_rate{false};
    int survival_rate_month{0};
    int survival_rate_day{0};
    // SI/SEI
    std::string model_type;
    int latency_period_steps;
    // Kernels
    std::string natural_kernel_type;
    double natural_scale{0};
    std::string natural_direction;
    double natural_kappa{0};
    bool use_anthropogenic_kernel{false};
    double percent_natural_dispersal{1};
    std::string anthro_kernel_type;
    double anthro_scale{0};
    std::string anthro_direction;
    std::string network_movement;  ///< walk, jump, teleport
    double network_min_distance{0};
    double network_max_distance{0};
    double anthro_kappa{0};
    double shape{1.0};
    // Treatments
    bool use_treatments{false};
    // Mortality
    bool use_mortality{false};
    std::string mortality_frequency;
    unsigned mortality_frequency_n;
    /** Mortality rate if used without pest-host-use table
     *
     * @see read_pest_host_use_table()
     */
    double mortality_rate{0};
    /** Time lag of mortality in simulation steps if used without pest-host-use table */
    int mortality_time_lag{0};
    // Quarantine
    bool use_quarantine{false};
    std::string quarantine_frequency;
    unsigned quarantine_frequency_n;
    std::string quarantine_directions;
    // Movements
    bool use_movements{false};
    std::vector<unsigned> movement_schedule;
    double dispersal_percentage{0.99};
    std::string output_frequency;
    unsigned output_frequency_n;
    bool use_spreadrates{true};
    std::string spreadrate_frequency;
    unsigned spreadrate_frequency_n;
    bool use_overpopulation_movements{false};
    double overpopulation_percentage{0};
    double leaving_percentage{0};
    double leaving_scale_coefficient{1};
    double dispersers_to_soils_percentage{0};  ///< Ratio of dispersers going into soil

    void create_schedules()
    {
        scheduler_ = Scheduler(date_start_, date_end_, step_unit_, step_num_units_);
        spread_schedule_ =
            scheduler_.schedule_spread(Season(season_start_month_, season_end_month_));
        output_schedule_ =
            schedule_from_string(scheduler_, output_frequency, output_frequency_n);
        if (use_mortality)
            mortality_schedule_ = schedule_from_string(
                scheduler_, mortality_frequency, mortality_frequency_n);
        if (use_lethal_temperature)
            lethal_schedule_ =
                scheduler_.schedule_action_yearly(lethal_temperature_month, 1);
        if (use_survival_rate)
            survival_rate_schedule_ = scheduler_.schedule_action_yearly(
                survival_rate_month, survival_rate_day);
        if (use_spreadrates)
            spread_rate_schedule_ = schedule_from_string(
                scheduler_, spreadrate_frequency, spreadrate_frequency_n);
        if (use_quarantine)
            quarantine_schedule_ = schedule_from_string(
                scheduler_, quarantine_frequency, quarantine_frequency_n);
        if (weather_size)
            weather_table_ = scheduler_.schedule_weather(weather_size);
        schedules_created_ = true;
    }

    const Scheduler& scheduler() const
    {
        if (!schedules_created_)
            throw std::logic_error(
                "Schedules were not created before calling scheduler()");
        return scheduler_;
    }

    const std::vector<bool>& spread_schedule() const
    {
        if (!schedules_created_)
            throw std::logic_error(
                "Schedules were not created before calling spread_schedule()");
        return spread_schedule_;
    }

    const std::vector<bool>& mortality_schedule() const
    {
        if (!schedules_created_)
            throw std::logic_error(
                "Schedules were not created before calling mortality_schedule()");
        return mortality_schedule_;
    }

    const std::vector<bool>& lethal_schedule() const
    {
        if (!use_lethal_temperature)
            throw std::logic_error(
                "lethal_schedule() not available when use_lethal_temperature is false");
        if (!schedules_created_)
            throw std::logic_error(
                "Schedules were not created before calling lethal_schedule()");
        return lethal_schedule_;
    }

    const std::vector<bool>& survival_rate_schedule() const
    {
        if (!use_survival_rate)
            throw std::logic_error(
                "survival_rate_schedule() not available when use_survival_rate is false");
        if (!schedules_created_)
            throw std::logic_error(
                "Schedules were not created before calling survival_rate_schedule()");
        return survival_rate_schedule_;
    }

    const std::vector<bool>& spread_rate_schedule() const
    {
        if (!use_spreadrates)
            throw std::logic_error(
                "spread_rate_schedule() not available when use_spreadrates is false");
        if (!schedules_created_)
            throw std::logic_error(
                "Schedules were not created before calling spread_rate_schedule()");
        return spread_rate_schedule_;
    }

    const std::vector<bool>& quarantine_schedule() const
    {
        if (!use_quarantine)
            throw std::logic_error(
                "quarantine_schedule() not available when use_quarantine is false");
        if (!schedules_created_)
            throw std::logic_error(
                "Schedules were not created before calling quarantine_schedule()");
        return quarantine_schedule_;
    }

    const std::vector<bool>& output_schedule() const
    {
        if (!schedules_created_)
            throw std::logic_error(
                "Schedules were not created before calling output_schedule()");
        return output_schedule_;
    }

    /**
     * @brief Get weather table for converting simulation steps to weather steps
     * @return Weather table as a vector of weather steps by reference
     */
    const std::vector<unsigned>& weather_table() const
    {
        if (!schedules_created_)
            throw std::logic_error(
                "Schedules were not created before calling weather_table()");
        if (!weather_size)
            throw std::logic_error(
                "weather_table() is not available when weather_size is zero");
        return weather_table_;
    }

    /**
     * @brief Convert simulation step to weather step
     * @param step Simulation step
     * @return Weather step (usable as an index of the weather array)
     */
    unsigned simulation_step_to_weather_step(unsigned step)
    {
        if (!schedules_created_)
            throw std::logic_error(
                "Schedules were not created before calling simulation_step_to_weather_step()");
        if (!weather_size)
            throw std::logic_error(
                "simulation_step_to_weather_step() is not available when weather_size is zero");
        return weather_table_.at(step);
    }

    unsigned num_mortality_steps()
    {
        if (!schedules_created_)
            throw std::logic_error(
                "Schedules were not created before calling num_mortality_steps()");
        return get_number_of_scheduled_actions(mortality_schedule_);
    }

    unsigned num_lethal()
    {
        if (!use_lethal_temperature)
            throw std::logic_error(
                "num_lethal() not available when use_lethal_temperature is false");
        if (!schedules_created_)
            throw std::logic_error(
                "Schedules were not created before calling num_lethal()");
        return get_number_of_scheduled_actions(lethal_schedule_);
    }

    unsigned num_survival_rate()
    {
        if (!use_survival_rate)
            throw std::logic_error(
                "num_survival_rate() not available when use_survival_rate is false");
        if (!schedules_created_)
            throw std::logic_error(
                "Schedules were not created before calling num_survival_rate()");
        return get_number_of_scheduled_actions(survival_rate_schedule_);
    }

    unsigned rate_num_steps()
    {
        if (!use_spreadrates)
            throw std::logic_error(
                "rate_num_steps() not available when use_spreadrates is false");
        if (!schedules_created_)
            throw std::logic_error(
                "Schedules were not created before calling rate_num_steps()");
        return get_number_of_scheduled_actions(spread_rate_schedule_);
    }

    unsigned quarantine_num_steps()
    {
        if (!use_quarantine)
            throw std::logic_error(
                "quarantine_num_steps() not available when use_quarantine is false");
        if (!schedules_created_)
            throw std::logic_error(
                "Schedules were not created before calling quarantine_num_steps()");
        return get_number_of_scheduled_actions(quarantine_schedule_);
    }

    const Date& date_start() const
    {
        return date_start_;
    }

    template<typename... Args>
    void set_date_start(Args&&... args)
    {
        date_start_ = Date(std::forward<Args>(args)...);
    }

    const Date& date_end() const
    {
        return date_end_;
    }

    template<typename... Args>
    void set_date_end(Args&&... args)
    {
        date_end_ = Date(std::forward<Args>(args)...);
    }

    StepUnit step_unit() const
    {
        return step_unit_;
    }

    void set_step_unit(StepUnit step_unit)
    {
        step_unit_ = step_unit;
    }

    void set_step_unit(const std::string& text)
    {
        step_unit_ = step_unit_enum_from_string(text);
    }

    unsigned step_num_units() const
    {
        return step_num_units_;
    }

    void set_step_num_units(unsigned step_num_units)
    {
        step_num_units_ = step_num_units;
    }

    // TODO: move to Season?
    void set_season_start_end_month(int start, int end)
    {
        season_start_month_ = start;
        season_end_month_ = end;
    }

    void set_season_start_end_month(const std::string& start, const std::string& end)
    {
        season_start_month_ = std::stoi(start);
        season_end_month_ = std::stoi(end);
    }

    /**
     * @brief Set disperser arrival behavior for model
     * @param value Arrival behavior string "infect" or "land"
     */
    void set_arrival_behavior(std::string value)
    {
        if (value != "infect" && value != "land") {
            throw std::invalid_argument(
                "arrival behavior can be 'infect' or 'land' but not: " + value);
        }
        arrival_behavior_ = value;
    }

    /**
     * @brief Get disperser arrival behavior for model
     * @return Arrival behavior as string "infect" or "land"
     */
    const std::string& arrival_behavior() const
    {
        return arrival_behavior_;
    }

    /**
     * Read seeds from text.
     *
     * @note All seeds are mandatory regardless of the other configuration value.
     *
     * @see read_key_value_pairs() for parameters and behavior.
     */
    void
    read_seeds(const std::string& text, char record_separator, char key_value_separator)
    {
        this->random_seeds = read_key_value_pairs<unsigned>(
            text, record_separator, key_value_separator, [](std::string text) {
                return std::stoul(text);
            });
        this->multiple_random_seeds = true;
    }

    /**
     * Read seeds from vector unsigned ints (list of integers).
     *
     * @note All seeds are mandatory regardless of the other configuration value.
     */
    void read_seeds(const std::vector<unsigned>& seeds)
    {
        static const std::vector<std::string> names{
            "disperser_generation",
            "natural_dispersal",
            "anthropogenic_dispersal",
            "establishment",
            "weather",
            "lethal_temperature",
            "movement",
            "overpopulation",
            "survival_rate",
            "soil"};
        if (names.size() != seeds.size()) {
            throw std::invalid_argument(
                "read_seeds: wrong number of seeds (" + std::to_string(seeds.size())
                + " instead of " + std::to_string(names.size()) + ")");
        }
        size_t i = 0;
        for (const auto& name : names) {
            this->random_seeds[name] = seeds.at(i);
            ++i;
        }
        this->multiple_random_seeds = true;
    }

    /**
     * @brief Get data for the pest-host-use table
     * @return Reference to the internal table
     *
     * @see PestHostUseTableDataRow
     */
    const std::vector<PestHostUseTableDataRow>& pest_host_use_table_data() const
    {
        return pest_host_use_table_data_;
    }

    /**
     * @brief Get data for the competency table
     * @return Reference to the internal table
     *
     * @see CompetencyTableDataRow
     */
    const std::vector<CompetencyTableDataRow>& competency_table_data() const
    {
        return competency_table_data_;
    }

    /**
     * @brief Read pest-host-use table data from vector of vectors of doubles
     *
     * The nested vectors need to be of size 3. The order of values is susceptibility,
     * mortality rate, and mortality time lag.
     *
     * @param values Table data
     */
    void read_pest_host_use_table(const std::vector<std::vector<double>>& values)
    {
        for (const auto& row : values) {
            if (row.size() < 3) {
                throw std::invalid_argument(
                    "3 values are required for each pest-host-use table row");
            }
            PestHostUseTableDataRow resulting_row;
            resulting_row.susceptibility = row[0];
            resulting_row.mortality_rate = row[1];
            resulting_row.mortality_time_lag = row[2];
            pest_host_use_table_data_.push_back(std::move(resulting_row));
        }
    }

    /**
     * @brief Use existing config parameters to create pest-host-use table
     *
     * This will create table with date for the given number of hosts with values for
     * all hosts being the same. Susceptibility is set to 1 and mortality is taken from
     * existing config attributes *mortality_rate* and *mortality_time_lag*.
     *
     * @param num_of_hosts Number of hosts
     *
     * @see #mortality_rate
     * @see #mortality_time_lag
     */
    void create_pest_host_use_table_from_parameters(int num_of_hosts)
    {
        for (int i = 0; i < num_of_hosts; ++i) {
            PestHostUseTableDataRow resulting_row;
            resulting_row.susceptibility = 1;
            resulting_row.mortality_rate = this->mortality_rate;
            resulting_row.mortality_time_lag = this->mortality_time_lag;
            pest_host_use_table_data_.push_back(std::move(resulting_row));
        }
    }

    /**
     * @brief Read competency table from vector of vectors of doubles
     *
     * The nested vectors are rows of the table which need to have size 2 or higher.
     * Each vector contains the combination of hosts and competency score.
     *
     * First n-1 items are the host presence and absence data which will be converted
     * from double to bool, i.e., use 0 for absence, 1 for presence. The number of these
     * items should be the number of hosts. Last (the nth) item in each vector is the
     * competency score for the given combination and will be used as double.
     *
     * For 1 host with competency 1, the table is `{{1, 1}, {0, 0}}`. For 2 hosts, the
     * table may look like this:
     *
     * ```
     * {
     *   {1, 0, 0.1},
     *   {0, 1, 0.4},
     *   {1, 1, 0.8},
     *   {0, 0, 0}
     * }
     * ```
     *
     * @param values Table data
     */
    void read_competency_table(const std::vector<std::vector<double>>& values)
    {
        for (const auto& row : values) {
            if (row.size() < 2) {
                throw std::invalid_argument(
                    "At least 2 values are required for each competency table row");
            }
            CompetencyTableDataRow resulting_row;
            for (auto it = row.begin(); it < std::prev(row.end()); ++it) {
                resulting_row.presence_absence.push_back(bool(*it));
            }
            resulting_row.competency = row.back();
            competency_table_data_.push_back(std::move(resulting_row));
        }
    }

private:
    Date date_start_{"0-01-01"};
    Date date_end_{"0-01-02"};

    int season_start_month_{1};
    int season_end_month_{12};

    StepUnit step_unit_{StepUnit::Day};
    unsigned step_num_units_{1};

    Scheduler scheduler_{date_start_, date_end_, step_unit_, step_num_units_};
    bool schedules_created_{false};

    std::string arrival_behavior_{"infect"};  ///< Disperser arrival behavior

    std::vector<bool> spread_schedule_;
    std::vector<bool> output_schedule_;
    std::vector<bool> mortality_schedule_;
    std::vector<bool> lethal_schedule_;
    std::vector<bool> survival_rate_schedule_;
    std::vector<bool> spread_rate_schedule_;
    std::vector<bool> quarantine_schedule_;
    std::vector<unsigned> weather_table_;

    /** Storage for the pest-host-use table data */
    std::vector<PestHostUseTableDataRow> pest_host_use_table_data_;
    /** Storage for the competency table data */
    std::vector<CompetencyTableDataRow> competency_table_data_;
};

}  // namespace pops

#endif  // POPS_CONFIG_HPP
