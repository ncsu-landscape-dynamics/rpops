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

/**
 * Pest-host-use table holding susceptibilities, mortality rates, and mortality time
 * lags for multiple hosts.
 */
template<typename HostPool>
class PestHostUseTable
{
public:
    using Environment = typename HostPool::Environment;

    /**
     * @brief Create an empty pest-host-use table
     *
     * @param environment Reference to the environment
     */
    PestHostUseTable(const Environment& environment) : environment_(environment) {}

    /**
     * @brief Create a pest-host-use table using values in config
     *
     * @param config Configuration with pest-host-use table data
     * @param environment Reference to the environment
     */
    PestHostUseTable(const Config& config, const Environment& environment)
        : environment_(environment)
    {
        for (const auto& row : config.pest_host_use_table_data()) {
            susceptibilities_.push_back(row.susceptibility);
            mortality_rates_.push_back(row.mortality_rate);
            mortality_time_lags_.push_back(row.mortality_time_lag);
        }
    }

    /**
     * @brief Add info for one host to the table
     *
     * Order of addition matters and should be the same as additions to the environment.
     *
     * @param susceptibility Host susceptibility
     * @param mortality_rate Host mortality rate
     * @param mortality_time_lag Host mortality time lag
     */
    void
    add_host_info(double susceptibility, double mortality_rate, int mortality_time_lag)
    {
        susceptibilities_.push_back(susceptibility);
        mortality_rates_.push_back(mortality_rate);
        mortality_time_lags_.push_back(mortality_time_lag);
    }

    /**
     * @brief Get susceptibility for the given host
     * @param host Pointer to the host to get the information for
     * @return Susceptibility score
     */
    double susceptibility(const HostPool* host) const
    {
        // This is using index because the environment is part of competency table,
        // otherwise a map which would use pointer to host would work here, too.
        auto host_index = environment_.host_index(host);
        return susceptibilities_.at(host_index);
    }

    /**
     * @brief Get mortality rate for the given host
     * @param host Pointer to the host to get the information for
     * @return Mortality rate value
     */
    double mortality_rate(const HostPool* host) const
    {
        auto host_index = environment_.host_index(host);
        return mortality_rates_.at(host_index);
    }

    /**
     * @brief Get mortality time lag for the given host
     * @param host Pointer to the host to get the information for
     * @return Mortality time lag value
     */
    double mortality_time_lag(const HostPool* host) const
    {
        auto host_index = environment_.host_index(host);
        return mortality_time_lags_.at(host_index);
    }

private:
    std::vector<double> susceptibilities_;  ///< List of susceptibilities for hosts
    std::vector<double> mortality_rates_;  ///< List of mortality_rates for hosts
    std::vector<int> mortality_time_lags_;  ///< List of mortality time lags for hosts
    const Environment& environment_;  ///< Environment used for host indexing
};

}  // namespace pops

#endif  // POPS_PEST_HOST_USE_TABLE_HPP
