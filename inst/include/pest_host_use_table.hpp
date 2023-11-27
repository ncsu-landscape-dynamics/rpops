/*
 * PoPS model - Pest-host-use table for hosts and pest
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

#ifndef POPS_PEST_HOST_USE_TABLE_HPP
#define POPS_PEST_HOST_USE_TABLE_HPP

#include <vector>
#include <stdexcept>
#include <string>

#include "config.hpp"

namespace pops {

template<typename HostPool>
class PestHostUseTable
{
public:
    using Environment = typename HostPool::Environment;

    PestHostUseTable(const Environment& environment) : environment_(environment) {}
    PestHostUseTable(const Config& config, const Environment& environment)
        : environment_(environment)
    {
        for (const auto& row : config.pest_host_use_table_data()) {
            susceptibilities_.push_back(row.susceptibility);
            mortality_rates_.push_back(row.mortality_rate);
            mortality_time_lags_.push_back(row.mortality_time_lag);
        }
    }

    void
    add_host_info(double susceptibility, double mortality_rate, int mortality_time_lag)
    {
        susceptibilities_.push_back(susceptibility);
        mortality_rates_.push_back(mortality_rate);
        mortality_time_lags_.push_back(mortality_time_lag);
    }

    double susceptibility(const HostPool* host) const
    {
        // This is using index because the environment is part of competency table,
        // otherwise a map which would use pointer to host would work here, too.
        auto host_index = environment_.host_index(host);
        return susceptibilities_.at(host_index);
    }

    double mortality_rate(const HostPool* host) const
    {
        auto host_index = environment_.host_index(host);
        return mortality_rates_.at(host_index);
    }

    double mortality_time_lag(const HostPool* host) const
    {
        auto host_index = environment_.host_index(host);
        return mortality_time_lags_.at(host_index);
    }

private:
    std::vector<double> susceptibilities_;
    std::vector<double> mortality_rates_;
    std::vector<int> mortality_time_lags_;
    const Environment& environment_;
};

}  // namespace pops

#endif  // POPS_PEST_HOST_USE_TABLE_HPP
