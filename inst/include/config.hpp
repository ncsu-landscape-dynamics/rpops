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
    double mortality_rate{0};
    int mortality_time_lag{0};  // Time lag of mortality in mortality steps
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

private:
    Date date_start_{"0-01-01"};
    Date date_end_{"0-01-02"};

    int season_start_month_{1};
    int season_end_month_{12};

    StepUnit step_unit_{StepUnit::Day};
    unsigned step_num_units_{1};

    Scheduler scheduler_{date_start_, date_end_, step_unit_, step_num_units_};
    bool schedules_created_{false};

    std::vector<bool> spread_schedule_;
    std::vector<bool> output_schedule_;
    std::vector<bool> mortality_schedule_;
    std::vector<bool> lethal_schedule_;
    std::vector<bool> survival_rate_schedule_;
    std::vector<bool> spread_rate_schedule_;
    std::vector<bool> quarantine_schedule_;
    std::vector<unsigned> weather_table_;
};

}  // namespace pops

#endif  // POPS_CONFIG_HPP
