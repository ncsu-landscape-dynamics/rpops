/*
 * PoPS model - environment for hosts and pests
 *
 * Copyright (C) 2022 by the authors.
 *
 * Authors: Vaclav Petras (wenzeslaus gmail com)
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#ifndef POPS_ENVIRONMENT_HPP
#define POPS_ENVIRONMENT_HPP

#include <cmath>
#include <memory>
#include <tuple>
#include <vector>
#include <random>
#include <string>
#include <stdexcept>

#include "environment_interface.hpp"
#include "normal_distribution_with_uniform_fallback.hpp"
#include "utils.hpp"
#include "host_pool_interface.hpp"
#include "generator_provider.hpp"

namespace pops {

/**
 * @brief Type of weather
 *
 * This includes all the ways of how weather coefficient for one step can be obtained.
 */
enum class WeatherType
{
    Deterministic,  ///< Weather is taken from a time series
    Probabilistic,  ///< Weather is generated from a distribution
    None,  ///< No weather
};

/*! Get a corresponding enum value for a string which represents a weather type.
 *
 * Throws an std::invalid_argument exception if the values was not
 * found or is not supported (which is the same thing).
 */
inline WeatherType weather_type_from_string(const std::string& text)
{
    if (text == "deterministic" || text == "Deterministic")
        return WeatherType::Deterministic;
    if (text == "probabilistic" || text == "Probabilistic")
        return WeatherType::Probabilistic;
    if (text.empty() || text == "none" || text == "None" || text == "NONE")
        return WeatherType::None;
    throw std::invalid_argument(
        "weather_type_from_string: Invalid value '" + text + "' provided");
}

/**
 * Encapsulates surrounding environment
 *
 * Currently, only handles weather coefficient for soils. Holds only the current state.
 */
template<
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex,
    typename Generator>
class Environment
    : public EnvironmentInterface<IntegerRaster, FloatRaster, RasterIndex, Generator>
{
public:
    Environment() {}

    /**
     * @brief Update the current weather coefficient
     *
     * @param raster Raster with the weather coefficient.
     */
    void update_weather_coefficient(const FloatRaster& raster) override
    {
        current_weather_coefficient = &raster;
        weather_ = true;
    }

    /**
     * @brief Update the current weather coefficient using mean and standard deviation
     *
     * Normal distribution is used to generate new value for each cell using *mean*
     * and *stddev*. Generated values which would fall out of the range for weather
     * coefficient, i.e., outside of 0-1 interval, are replaced by random value from
     * a uniform distribution.
     *
     * The values in *mean* are checked to be in the 0-1 interval because it is assumed
     * that the range for mean should be the same as for the actual coefficient value.
     * The values in *stddev* are not checked.
     *
     * @param mean Raster of mean weather coefficient for each cell
     * @param stddev Raster of standard deviation of weather coefficient for each cell
     * @param generator Random number generator provider
     *
     * @throw std::invalid_argument when dimensions of *mean* and *stddev* differ or
     * when mean is out of range
     */
    void update_weather_from_distribution(
        const FloatRaster& mean,
        const FloatRaster& stddev,
        Generator& generator) override
    {
        if (mean.rows() != stddev.rows()) {
            throw std::invalid_argument(
                "Mean and stddev need to have the same number of rows ("
                + std::to_string(mean.rows()) + " != " + std::to_string(stddev.rows())
                + ")");
        }
        if (mean.cols() != stddev.cols()) {
            throw std::invalid_argument(
                "Mean and stddev need to have the same number of columns ("
                + std::to_string(mean.cols()) + " != " + std::to_string(stddev.cols())
                + ")");
        }
        stored_weather_coefficient = FloatRaster(mean.rows(), mean.cols());
        // possibly use suitable cells here
        for (RasterIndex i = 0; i < mean.rows(); ++i) {
            for (RasterIndex j = 0; j < mean.cols(); ++j) {
                auto mean_value = mean(i, j);
                // In general, to get a specific shape of the distribution, mean can be
                // anything, but we limit that assuming that mean which is out of the
                // desired coefficient range is an erroneous value.
                // Notably, this test is not perfomed for deterministic weather.
                if (mean_value < weather_coefficient_min
                    || mean_value > weather_coefficient_max) {
                    throw std::invalid_argument(
                        std::string("Weather coefficient mean is expected to be ")
                        + "between " + std::to_string(weather_coefficient_min) + " and "
                        + std::to_string(weather_coefficient_max) + ", but is "
                        + std::to_string(mean_value) + " at (" + std::to_string(i)
                        + ", " + std::to_string(j) + ")");
                }
                NormalDistributionWithUniformFallback<double> distribution{
                    mean_value,
                    stddev(i, j),
                    weather_coefficient_min,
                    weather_coefficient_max};
                stored_weather_coefficient(i, j) = distribution(generator.weather());
            }
        }
        current_weather_coefficient = &stored_weather_coefficient;
        weather_ = true;
    }

    /**
     * @brief Get weather coefficient at a given cell
     *
     * @param row Cell row number
     * @param col Cell column number
     * @return Current value at the given cell
     *
     * @throw std::logic_error when coefficient is not set
     */
    double weather_coefficient_at(RasterIndex row, RasterIndex col) const override
    {
        if (!current_weather_coefficient) {
            throw std::logic_error("Weather coefficient used, but not provided");
        }
        return current_weather_coefficient->operator()(row, col);
    }

    double influence_reproductive_rate_at(
        RasterIndex row, RasterIndex col, double value) const override
    {
        if (!weather_)
            return value;
        return value * weather_coefficient_at(row, col);
    }

    double influence_suitability_at(
        RasterIndex row, RasterIndex col, double value) const override
    {
        if (!weather_)
            return value;
        return value * weather_coefficient_at(row, col);
    }

    int total_population_at(RasterIndex row, RasterIndex col) const override
    {
        // If total population is used, use that instead of computing it.
        if (total_population_)
            return total_population_->operator()(row, col);
        int sum = 0;
        if (other_individuals_)
            sum += other_individuals_->operator()(row, col);
        for (const auto& host : hosts_)
            sum += host->total_hosts_at(row, col);
        return sum;
    }

    /**
     * @copydoc EnvironmentInterface::host_presence_at()
     *
     * Hosts are in the order of how host pools were registered to the environment.
     */
    std::vector<bool> host_presence_at(RasterIndex row, RasterIndex col) const override
    {
        std::vector<bool> presence;
        presence.reserve(hosts_.size());
        for (const auto& host : hosts_)
            presence.push_back(host->total_hosts_at(row, col));
        return presence;
    }

    void set_other_individuals(const IntegerRaster* individuals) override
    {
        other_individuals_ = individuals;
    }

    void set_total_population(const IntegerRaster* individuals) override
    {
        total_population_ = individuals;
    }

    // While copydoc generates correct doc here, it triggers a
    // warning with Doxygen 1.8.17, so using see here for now.
    /**
     * @see EnvironmentInterface::add_host()
     *
     * @note The function is no-op if host pool already registered. This may throw an
     * exception in the future.
     */
    void add_host(const HostPoolInterface<RasterIndex>* host) override
    {
        if (container_contains(hosts_, host))
            return;
        hosts_.push_back(host);
    }

    /**
     * @brief Remove all hosts from the environment.
     *
     * This function is useful for reusing an environment object in different contexts,
     * especially in tests, but it does not have an epidemiological meaning.
     */
    void remove_hosts()
    {
        hosts_.clear();
    }

    /**
     * @copydoc EnvironmentInterface::host_index()
     */
    size_t host_index(const HostPoolInterface<RasterIndex>* host) const override
    {
        auto it = std::find(hosts_.begin(), hosts_.end(), host);
        if (it == hosts_.end())
            throw std::invalid_argument(
                "Environment::host_index: Host is not in the environment");
        return std::distance(hosts_.begin(), it);
    }

    /**
     * @brief Get weather coefficient raster
     *
     * @return Reference to the current weather coefficient raster
     *
     * @throw std::logic_error when coefficient is not set
     */
    const FloatRaster& weather_coefficient() const override
    {
        if (!current_weather_coefficient) {
            throw std::logic_error("Weather coefficient used, but not provided");
        }
        return *current_weather_coefficient;
    }

    void update_temperature(const FloatRaster& raster) override
    {
        temperature_ = &raster;
    }

    double temperature_at(RasterIndex row, RasterIndex col) const override
    {
        if (!temperature_) {
            throw std::logic_error("Temperature used, but not provided");
        }
        return temperature_->operator()(row, col);
    }

protected:
    static constexpr double weather_coefficient_min = 0;
    static constexpr double weather_coefficient_max = 1;

    /**
     * Current weather coefficient
     *
     * Value may not be set and these cases should produce exceptions.
     */
    const FloatRaster* current_weather_coefficient{nullptr};
    FloatRaster stored_weather_coefficient;
    bool weather_{false};

    std::vector<const HostPoolInterface<RasterIndex>*> hosts_;  // host, non-owning
    const IntegerRaster* other_individuals_{nullptr};  // non-hosts, non-owning
    const IntegerRaster* total_population_{nullptr};  // non-hosts, non-owning

    const FloatRaster* temperature_{nullptr};
};

}  // namespace pops

#endif  // POPS_ENVIRONMENT_HPP
