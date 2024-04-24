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

#ifndef POPS_MULTI_HOST_POOL_HPP
#define POPS_MULTI_HOST_POOL_HPP

#include <vector>
#include <algorithm>
#include <random>

#include "competency_table.hpp"
#include "config.hpp"
#include "pest_host_table.hpp"
#include "utils.hpp"

namespace pops {

/**
 * Host pool for multiple hosts
 *
 * Keeps the same interface as (single) HostPool given that HostPool has methods to
 * behave in a multi-host way if needed. This allows most external operations such as
 * actions to work without distinguishing single- and multi-host pools.
 */
template<
    typename HostPool,
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex,
    typename GeneratorProvider>
class MultiHostPool
{
public:
    /**
     * Standard random number generator to be passed directly to the methods.
     */
    using Generator = typename GeneratorProvider::Generator;

    /**
     * @brief Create MultiHostPool object with given host pools and configuration
     *
     * @param host_pools List of host pools to use
     * @param config Configuration to use
     */
    MultiHostPool(const std::vector<HostPool*>& host_pools, const Config& config)
        : host_pools_(host_pools), config_(config)
    {}

    /**
     * @brief Set pest-host table for all hosts
     *
     * The existing object will be used (not copy is performed).
     *
     * @param table Reference to the table object
     */
    void set_pest_host_table(const PestHostTable<HostPool>& table)
    {
        for (auto& host_pool : host_pools_) {
            host_pool->set_pest_host_table(table);
        }
    }

    /**
     * @brief Set competency table for all hosts
     *
     * The existing object will be used (not copy is performed).
     *
     * @param table Reference to the table object
     */
    void set_competency_table(const CompetencyTable<HostPool>& table)
    {
        for (auto& host_pool : host_pools_) {
            host_pool->set_competency_table(table);
        }
    }

    /**
     * @brief Get suitable cells spatial index
     *
     * Always uses the index only from the first single-host pool
     * assuming that the indices are the same.
     * Assumes that at least one single-host pool was added.
     *
     * @return Const reference to the index
     */
    const std::vector<std::vector<int>>& suitable_cells() const
    {
        return host_pools_.at(0)->suitable_cells();
    }

    /**
     * @brief Check whether the cell is outside of the raster extent
     *
     * Always uses the only the first single-host pool for the check
     * assuming that the extents are the same.
     * Assumes that at least one single-host pool was added.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @return true if outside of the raster, false if inside
     */
    bool is_outside(RasterIndex row, RasterIndex col)
    {
        return host_pools_.at(0)->is_outside(row, col);
    }

    /**
     * @brief Step all hosts forward
     * @param step Step in the simulation (>=0)
     *
     * @see HostPool::step_forward()
     */
    void step_forward(unsigned step)
    {
        for (auto& host_pool : host_pools_) {
            host_pool->step_forward(step);
        }
    }

    /**
     * @copydoc HostPool::remove_all_infected_at()
     */
    void remove_all_infected_at(RasterIndex row, RasterIndex col, Generator& generator)
    {
        for (auto& host_pool : host_pools_) {
            host_pool->remove_all_infected_at(row, col, generator);
        }
    }

    /**
     * @copydoc HostPool::remove_infection_by_ratio_at()
     */
    void remove_infection_by_ratio_at(
        RasterIndex row, RasterIndex col, double ratio, Generator& generator)
    {
        for (auto& host_pool : host_pools_) {
            host_pool->remove_infection_by_ratio_at(row, col, ratio, generator);
        }
    }

    /**
     * @copydoc HostPool::apply_mortality_at(RasterIndex, RasterIndex, double, int)
     */
    void apply_mortality_at(
        RasterIndex row, RasterIndex col, double mortality_rate, int mortality_time_lag)
    {
        for (auto& host_pool : host_pools_) {
            host_pool->apply_mortality_at(row, col, mortality_rate, mortality_time_lag);
        }
    }

    /**
     * @copydoc HostPool::apply_mortality_at(RasterIndex, RasterIndex)
     */
    void apply_mortality_at(RasterIndex row, RasterIndex col)
    {
        for (auto& host_pool : host_pools_) {
            host_pool->apply_mortality_at(row, col);
        }
    }

    /**
     * @copydoc HostPool::step_forward_mortality()
     */
    void step_forward_mortality()
    {
        for (auto& host_pool : host_pools_) {
            host_pool->step_forward_mortality();
        }
    }

    /**
     * @brief Total number of all hosts over all host pools
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     *
     * @return Number of hosts
     */
    int total_hosts_at(RasterIndex row, RasterIndex col) const
    {
        int sum = 0;
        for (auto& host_pool : host_pools_) {
            sum += host_pool->total_hosts_at(row, col);
        }
        return sum;
    }

    /**
     * @brief Get number of infected hosts over all host pools at a given cell
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     *
     * @return Number of infected hosts
     */
    int infected_at(RasterIndex row, RasterIndex col) const
    {
        int infected = 0;
        for (auto& host_pool : host_pools_) {
            infected += host_pool->infected_at(row, col);
        }
        return infected;
    }

    /**
     * @copydoc HostPool::dispersers_from()
     */
    int dispersers_from(RasterIndex row, RasterIndex col, Generator& generator) const
    {
        int sum{0};
        for (auto& host_pool : host_pools_) {
            sum += host_pool->dispersers_from(row, col, generator);
        }
        return sum;
    }

    /**
     * @brief Move pests from a cell
     *
     * This is a multi-host version of single host pests_from method.
     * It randomly distributes the number to un-infect between multiple
     * hosts.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param count Number of pests requested to move from the cell
     * @param generator Random number generator
     *
     * @return Number of pests actually moved from the cell
     *
     * @note For consitency with the previous implementation, this does not modify
     * mortality cohorts nor touches the exposed cohorts.
     */
    int pests_from(RasterIndex row, RasterIndex col, int count, Generator& generator)
    {
        std::vector<int> infected;
        int index = 0;
        for (auto& host_pool : host_pools_) {
            infected.insert(infected.end(), host_pool->infected_at(row, col), index);
            index++;
        }

        index = 0;
        int collect_count = 0;
        std::vector<int> draw = draw_n_from_v(infected, count, generator);
        for (auto& host_pool : host_pools_) {
            count = std::count(draw.begin(), draw.end(), index);
            collect_count += host_pool->pests_from(row, col, count, generator);
            index++;
        }

        return collect_count;
    }
    /**
     * @brief Move pests to a cell
     *
     * This is a multi-host version of single host pests_to method.
     * It randomly distributes the number to infect between multiple
     * hosts.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param count Number of pests requested to move to the cell
     * @param generator Random number generator
     *
     * @return Number of accepted pests
     *
     * @note For consistency with the previous implementation, this does not make hosts
     * exposed in the SEI model. Instead, the hosts are infected right away. This may
     * become a feature in the future.
     *
     * @note For consistency with the previous implementation, this does not modify the
     * mortality cohorts. This wil need to be fixed in the future.
     *
     * @note This may be merged with add_disperser_at() in the future.
     */
    int pests_to(RasterIndex row, RasterIndex col, int count, Generator& generator)
    {
        std::vector<int> susceptible;
        int index = 0;
        for (auto& host_pool : host_pools_) {
            susceptible.insert(
                susceptible.end(), host_pool->susceptible_at(row, col), index);
            index++;
        }
        index = 0;
        int collect_count = 0;
        std::vector<int> draw = draw_n_from_v(susceptible, count, generator);
        for (auto& host_pool : host_pools_) {
            count = std::count(draw.begin(), draw.end(), index);
            collect_count += host_pool->pests_to(row, col, count, generator);
            index++;
        }

        return collect_count;
    }

    /**
     * @brief Move a disperser to a cell with establishment test
     *
     * Config::arrival_behavior() is used to determine if the dispersers should be
     * landing in a cell (`"land"`) or infecting a specific host directly (`"infect"`).
     * Landing performs the establishment test on a combined suitability from all hosts
     * in a cell and then uses one host to add the  new disperser to. Infecting picks a
     * host and lets the host accept or reject the disperser based on its own
     * establishment test.
     *
     * Uses static function HostPool::can_disperser_establish() from HostPool class to
     * do the test for multiple hosts for landing, so any host pool class needs to
     * implement that besides its disperser methods. It assumes that the configuration
     * for multi-host is applicable for the combined suitability of individual hosts.
     *
     * For one single host pool only and host pool which implements its `disperser_to()`
     * using `can_disperser_establish()` and `add_disperser_at()`, this gives identical
     * results for both landing and infecting arrival behavior.
     *
     * @see HostPool::can_disperser_establish()
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param generator Random number generator
     *
     * @return 1 if established, 0 otherwise
     */
    int disperser_to(RasterIndex row, RasterIndex col, Generator& generator)
    {
        std::vector<double> suitabilities;
        double total_suitability_score = 0;
        for (auto& host_pool : host_pools_) {
            // The suitability accounts for weather and, importantly, number of
            // susceptible hosts, so host pool without available hosts in a given cell
            // is less likely to be selected over a host pool with available hosts
            // (however, it can happen in case zeros are not exactly zeros in the
            // discrete distribution used later and code checks can't tell, so we need
            // to account for that case later anyway).
            double suitability = host_pool->suitability_at(row, col);
            // The resulting individual suitability can be 0-1. The individual
            // suitabilities are used as weights for picking the host, so their absolute
            // range does not matter. The total is used in a stochastic test and should
            // be <=1 which should be ensured in the input data.
            suitabilities.push_back(suitability);
            total_suitability_score += suitability;
        }
        if (total_suitability_score <= 0) {
            // While the score should always be >= 0, it may be == 0 if no hosts (host
            // individuals) are present. No hosts present cause all suitabilities to be
            // zero which is not permissible for the host picking later and it is enough
            // information for us to know there won't be any establishment.
            return 0;
        }
        if (total_suitability_score > 1) {
            throw std::invalid_argument(
                "Total suitability score is " + std::to_string(total_suitability_score)
                + " but it needs to be <=1");
        }
        auto host = pick_host_by_weight(host_pools_, suitabilities, generator);
        if (config_.arrival_behavior() == "land") {
            // The operations are ordered so that for single host, this gives an
            // identical results to the infect behavior (influenced by usage of random
            // numbers and presence of susceptible hosts).
            if (host->susceptible_at(row, col) <= 0)
                return 0;
            bool establish = HostPool::can_disperser_establish(
                total_suitability_score,
                config_.establishment_stochasticity,
                config_.establishment_probability,
                generator);
            if (establish)
                return host->add_disperser_at(row, col);  // simply increases the counts
            return 0;
        }
        else if (config_.arrival_behavior() == "infect") {
            return host->disperser_to(row, col, generator);  // with establishment test
        }
        throw std::invalid_argument(
            "Unknown value for arrival_behavior: " + config_.arrival_behavior());
    }
    /**
     * @brief Move hosts from a cell to a cell (multi-host)
     *
     * Currently just works for first host.
     *
     * @param row_from Row index of the source cell
     * @param col_from Column index of the source cell
     * @param row_to Row index of the target cell
     * @param col_to Column index of the target cell
     * @param count Number of pests to move
     * @param generator Random number generator for distribution of host among pools
     *
     * @return Number of pests actually moved between cells
     *
     * @note Mortality is not supported.
     */
    int move_hosts_from_to(
        RasterIndex row_from,
        RasterIndex col_from,
        RasterIndex row_to,
        RasterIndex col_to,
        int count,
        Generator& generator)
    {
        return host_pools_[0]->move_hosts_from_to(
            row_from, col_from, row_to, col_to, count, generator);
    }

    /**
     * @brief Get list of host pools
     * @return Reference to host pools
     */
    const std::vector<HostPool*>& host_pools()
    {
        return host_pools_;
    }

private:
    /**
     * @brief Pick a host given a weight for each host
     *
     * The weights don't need to be 0-1 nor they need to add up to 1. However, their sum
     * needs to be >0.
     *
     * If there is only one host, it returns that host without using the random number
     * generator.
     *
     * @param hosts List of pointers to host pools
     * @param weights Weight values for each host pool
     * @param generator Random number generator
     *
     * @return Pointer to selected host pool
     *
     * @throw std::invalid_argument if *hosts* is empty.
     */
    static HostPool* pick_host_by_weight(
        std::vector<HostPool*>& hosts,
        const std::vector<double>& weights,
        Generator& generator)
    {
        if (!hosts.size()) {
            std::invalid_argument("List of hosts is empty in multi host pool");
        }
        if (hosts.size() == 1) {
            return hosts[0];
        }
        std::discrete_distribution<int> distribution{weights.begin(), weights.end()};
        return hosts.at(distribution(generator));
    }

    /**
     * List of non-owning pointers to individual host pools.
     */
    std::vector<HostPool*> host_pools_;
    /**
     * Reference to configuration
     */
    const Config& config_;
};

}  // namespace pops

#endif  // POPS_MULTI_HOST_POOL_HPP
