/*
 * PoPS model - pest or pathogen spread simulation
 *
 * Copyright (C) 2023 by the authors.
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

#ifndef POPS_ENVIRONMENT_INTERFACE_HPP
#define POPS_ENVIRONMENT_INTERFACE_HPP

#include "host_pool_interface.hpp"

namespace pops {

template<
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex,
    typename Generator>
class EnvironmentInterface
{
public:
    virtual ~EnvironmentInterface() {}
    virtual void update_weather_coefficient(const FloatRaster& raster) = 0;
    virtual void update_weather_from_distribution(
        const FloatRaster& mean, const FloatRaster& stddev, Generator& generator) = 0;
    virtual double weather_coefficient_at(RasterIndex row, RasterIndex col) const = 0;
    virtual double influence_reproductive_rate_at(
        RasterIndex row, RasterIndex col, double value) const = 0;
    virtual double influence_probability_of_establishment_at(
        RasterIndex row, RasterIndex col, double value) const = 0;
    virtual int total_population_at(RasterIndex row, RasterIndex col) const = 0;
    virtual std::vector<bool>
    host_presence_at(RasterIndex row, RasterIndex col) const = 0;
    // in general, we may have other individuals-non const
    virtual void set_other_individuals(const IntegerRaster* individuals) = 0;
    virtual void set_total_population(const IntegerRaster* individuals) = 0;
    virtual void add_host(const HostPoolInterface<RasterIndex>* host) = 0;
    virtual size_t host_index(const HostPoolInterface<RasterIndex>* host) const = 0;
    virtual const FloatRaster& weather_coefficient() const = 0;
    virtual void update_temperature(const FloatRaster& raster) = 0;
    virtual double temperature_at(RasterIndex row, RasterIndex col) const = 0;
};

}  // namespace pops

#endif  // POPS_ENVIRONMENT_INTERFACE_HPP
