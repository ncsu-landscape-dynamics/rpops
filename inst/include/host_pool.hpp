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

#ifndef POPS_HOST_POOL_HPP
#define POPS_HOST_POOL_HPP

#include <vector>
#include <random>
#include <stdexcept>
#include <algorithm>

#include "host_pool_interface.hpp"
#include "model_type.hpp"
#include "environment_interface.hpp"
#include "competency_table.hpp"
#include "pest_host_table.hpp"
#include "utils.hpp"

namespace pops {

/**
 * Host pool managing susceptible, exposed, infected, and resistant hosts.
 *
 * @tparam IntegerRaster Integer raster type
 * @tparam FloatRaster Floating point raster type
 * @tparam RasterIndex Type for indexing the rasters
 * @tparam GeneratorProvider Provider of random number generators
 *
 * GeneratorProvider needs to provide Generator member which is the type of the
 * underlying random number generators.
 */
template<
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex,
    typename GeneratorProvider>
class HostPool : public HostPoolInterface<RasterIndex>
{
public:
    /**
     * Type of environment object providing information about weather and other
     * environmental properties.
     */
    using Environment = EnvironmentInterface<
        IntegerRaster,
        FloatRaster,
        RasterIndex,
        GeneratorProvider>;
    /**
     * Standard random number generator to be passed directly to the methods.
     */
    using Generator = typename GeneratorProvider::Generator;

    /**
     * @brief Creates an object with stored references and host properties.
     *
     * The *exposed* vector is a list of hosts exposed in the previous steps.
     * The length of the vector is the number of steps of the latency
     * period plus one. See the step_forward() function for details on how
     * the vector is used.
     *
     * The *total_populations* can be total number of hosts in the basic case
     * or it can be the total size of population of all relevant species
     * both host and non-host if dilution effect should be applied.
     *
     * If establishment stochasticity is disabled,
     * *establishment_probability* is used to decide whether or not
     * a disperser is established in a cell. Value 1 means that all
     * dispersers will establish and value 0 means that no dispersers
     * will establish.
     *
     * The *mortality_tracker_vector* is a vector of matrices for tracking infected
     * host infection over time. Expectation is that mortality tracker is of
     * length (1/mortality_rate + mortality_time_lag).
     *
     * Host is added to the environment by the constructor. Afterwards, the environment
     * is not modified.
     *
     * @param model_type Type of the model (SI or SEI)
     * @param susceptible Raster of susceptible hosts
     * @param exposed Raster of exposed or infected hosts
     * @param latency_period Length of the latency period in steps
     * @param infected Infected hosts
     * @param total_exposed Raster tracking all exposed hosts
     * @param resistant Resistant hosts
     * @param mortality_tracker_vector Raster tracking hosts for mortality
     * @param died Raster tracking all hosts which died
     * @param total_hosts Total number of hosts
     * @param environment Environment which influences the processes
     * @param dispersers_stochasticity Enable stochasticity in generating of dispersers
     * @param reproductive_rate Reproductive rate in ideal conditions
     * @param establishment_stochasticity Enable stochasticity in establishment step
     * @param establishment_probability Fixed probability disperser establishment
     * @param rows Number of rows for the study area
     * @param cols Number of columns for the study area
     * @param suitable_cells Cells where hosts are known to occur
     */
    HostPool(
        ModelType model_type,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& exposed,
        unsigned latency_period,
        IntegerRaster& infected,
        IntegerRaster& total_exposed,
        IntegerRaster& resistant,
        std::vector<IntegerRaster>& mortality_tracker_vector,
        IntegerRaster& died,
        IntegerRaster& total_hosts,
        Environment& environment,
        bool dispersers_stochasticity,
        double reproductive_rate,
        bool establishment_stochasticity,
        double establishment_probability,
        RasterIndex rows,
        RasterIndex cols,
        std::vector<std::vector<int>>& suitable_cells)
        : susceptible_(susceptible),
          infected_(infected),
          exposed_(exposed),
          latency_period_(latency_period),
          total_exposed_(total_exposed),
          resistant_(resistant),
          mortality_tracker_vector_(mortality_tracker_vector),
          died_(died),
          total_hosts_(total_hosts),
          environment_(environment),
          model_type_(model_type),
          dispersers_stochasticity_(dispersers_stochasticity),
          reproductive_rate_(reproductive_rate),
          establishment_stochasticity_(establishment_stochasticity),
          deterministic_establishment_probability_(establishment_probability),
          rows_(rows),
          cols_(cols),
          suitable_cells_(suitable_cells)
    {
        environment.add_host(this);
    }

    /**
     * @brief Set pest-host table
     *
     * If set, apply_mortality_at(RasterIndex, RasterIndex) can be used instead of the
     * version which takes mortality parameters directly. Susceptibility is used
     * automatically if the table is set.
     *
     * Pointer to the existing object is stored and used. So, the table can be modified,
     * but the table object needs to exists during the lifetime of this object.
     *
     * @param table Reference to the table
     */
    void set_pest_host_table(const PestHostTable<HostPool>& table)
    {
        this->pest_host_table_ = &table;
    }

    /**
     * @brief Set competency table
     *
     * Competency is used automatically if the table is set. No competency is considered
     * if the table is not set.
     *
     * Pointer to the existing object is stored and used. So, the table can be modified,
     * but the table object needs to exists during the lifetime of this object.
     *
     * @param table Reference to the table
     */
    void set_competency_table(const CompetencyTable<HostPool>& table)
    {
        this->competency_table_ = &table;
    }

    /**
     * @brief Move disperser to a cell in the host pool
     *
     * Processes event when a disperser lands in a cell potentially establishing on a
     * host. The disperser may or may not establish a based on host availability,
     * weather, establishment probability, and stochasticity.
     *
     * Any new dispersers targeting host in the host pool should be processed using this
     * function.
     *
     * @param row Row number of the target cell
     * @param col Column number of the target cell
     * @param generator Random number generator
     *
     * @return true if disperser has established in the cell, false otherwise
     *
     * @throw std::runtime_error if model type is unsupported (i.e., not SI or SEI)
     */
    int disperser_to(RasterIndex row, RasterIndex col, Generator& generator)
    {
        if (susceptible_(row, col) <= 0)
            return 0;
        double probability_of_establishment = suitability_at(row, col);
        bool establish = can_disperser_establish(
            probability_of_establishment,
            establishment_stochasticity_,
            deterministic_establishment_probability_,
            generator);
        if (establish)
            return add_disperser_at(row, col);
        return 0;
    }

    /**
     * @brief Test whether a disperser can establish
     *
     * This static (object-independent) function to test disperser establishment allows
     * code reuse between a single host and a multi-host case.
     *
     * @param probability_of_establishment Probability of establishment
     * @param establishment_stochasticity true if establishment stochasticity is enabled
     * @param deterministic_establishment_probability Establishment probability for
     * deterministic establishment
     * @param generator Random number generator (used with enabled stochasticity)
     * @return true if disperser can establish, false otherwise
     */
    static bool can_disperser_establish(
        double probability_of_establishment,
        bool establishment_stochasticity,
        double deterministic_establishment_probability,
        Generator& generator)
    {
        std::uniform_real_distribution<double> distribution_uniform(0.0, 1.0);
        double establishment_tester = 1 - deterministic_establishment_probability;
        if (establishment_stochasticity)
            establishment_tester = distribution_uniform(generator);
        if (establishment_tester < probability_of_establishment)
            return true;
        return false;
    }

    /**
     * @brief Add disperser to a cell
     *
     * Turns disperser into infection considering model type (SI, SEI).
     *
     * Unlike disperser_to(), this is transforming a disperser into infection right away
     * without any further evaluation of establishment or stochasticity.
     *
     * @param row Row number of the target cell
     * @param col Column number of the target cell
     *
     * @return 1 if disperser was added, 0 if there was not available host.
     *
     * @throw std::runtime_error if model type is unsupported (i.e., not SI or SEI)
     *
     * @note This may be merged with pests_to() in the future.
     */
    int add_disperser_at(RasterIndex row, RasterIndex col)
    {
        if (susceptible_(row, col) <= 0)
            return 0;
        susceptible_(row, col) -= 1;
        if (model_type_ == ModelType::SusceptibleInfected) {
            infected_(row, col) += 1;
            mortality_tracker_vector_.back()(row, col) += 1;
        }
        else if (model_type_ == ModelType::SusceptibleExposedInfected) {
            exposed_.back()(row, col) += 1;
            total_exposed_(row, col) += 1;
        }
        else {
            throw std::runtime_error(
                "Unknown ModelType value in HostPool::add_disperser_at()");
        }
        return 1;
    }

    /**
     * @brief Get dispersers produced in a cell
     *
     * Each time the function is called it generates number of dispersers based on the
     * current infection, host attributes, and the environment.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param generator Random number generator
     * @return Number of generated dispersers
     */
    int dispersers_from(RasterIndex row, RasterIndex col, Generator& generator) const
    {
        if (infected_at(row, col) <= 0)
            return 0;
        double lambda =
            environment_.influence_reproductive_rate_at(row, col, reproductive_rate_);
        if (competency_table_) {
            lambda *= competency_table_->competency_at(row, col, this);
        }
        int dispersers_from_cell = 0;
        if (dispersers_stochasticity_) {
            std::poisson_distribution<int> distribution(lambda);
            for (int k = 0; k < infected_at(row, col); k++) {
                dispersers_from_cell += distribution(generator);
            }
        }
        else {
            dispersers_from_cell =
                static_cast<int>(std::floor(lambda * infected_at(row, col)));
        }
        return dispersers_from_cell;
    }

    /**
     * @brief Get suitability score for a cell
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     *
     * @return suitability score
     */
    double suitability_at(RasterIndex row, RasterIndex col) const
    {
        double suitability = (double)(susceptible_(row, col))
                             / environment_.total_population_at(row, col);
        if (pest_host_table_) {
            suitability *= pest_host_table_->susceptibility(this);
        }
        suitability = environment_.influence_suitability_at(row, col, suitability);
        if (suitability < 0 || suitability > 1) {
            throw std::invalid_argument(
                "Suitability should be >=0 and <=1, not " + std::to_string(suitability)
                + " (susceptible: " + std::to_string(susceptible_(row, col))
                + ", total population: "
                + std::to_string(environment_.total_population_at(row, col))
                + ", susceptibility: "
                + std::to_string(pest_host_table_->susceptibility(this)) + ")");
        }
        return suitability;
    }

    /**
     * @brief Move pests from a cell
     *
     * Moves pest from a cell reducing number of infected and increasing number of
     * susceptible hosts in the cell.
     *
     * This is a request to move the number of pests given by *count* from a cell.
     * The function may reduce the number based on the availability of pests in the cell
     * or internal mechanics of the pool. However, the current implementation simply
     * moves the given number of pests without any checks.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param count Number of pests requested to move from the cell
     * @param generator Random number generator (for compatibility with multi-host API)
     *
     * @return Number of pests actually moved from the cell
     *
     * @note For consistency with the previous implementation, this does not modify
     * mortality cohorts nor touches the exposed cohorts.
     */
    int pests_from(RasterIndex row, RasterIndex col, int count, Generator& generator)
    {
        UNUSED(generator);
        susceptible_(row, col) += count;
        infected_(row, col) -= count;
        return count;
    }

    /**
     * @brief Move pests to a cell
     *
     * This directly modifies the number of infected hosts. Unlike disperser_to(), this
     * is assuming the pests will establish right away without any further evaluation of
     * establishment or stochasticity.
     *
     * When there is more pests than the target cell can accept based on the number of
     * susceptible hosts, all susceptible hosts are infected. The number of actually
     * infected hosts (accepted pests) is returned.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param count Number of pests requested to move to the cell
     * @param generator Random number generator (for compatibility with multi-host API)
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
        UNUSED(generator);
        // The target cell can accept all.
        if (susceptible_(row, col) >= count) {
            susceptible_(row, col) -= count;
            infected_(row, col) += count;
        }
        // More pests than the target cell can accept.
        else {
            count = susceptible_(row, col);
            susceptible_(row, col) -= count;
            infected_(row, col) += count;
        }
        return count;
    }

    //
    /**
     * @brief Move hosts from a cell to a cell
     *
     * Moves specified number of hosts from a cell to another cell. The distribution
     * of the hosts among susceptible, exposed, and other pools is stochastic.
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
        int total_hosts_moved{count};
        if (count > total_hosts_(row_from, col_from)) {
            total_hosts_moved = total_hosts_(row_from, col_from);
        }
        int total_infecteds = infected_(row_from, col_from);
        int suscepts = susceptible_(row_from, col_from);
        int expose = total_exposed_(row_from, col_from);
        int resist = resistant_(row_from, col_from);
        // Set up vector of numeric categories (infected = 1, susceptible = 2,
        // exposed = 3) for drawing number of moved in each category.
        std::vector<int> categories(total_infecteds, 1);
        categories.insert(categories.end(), suscepts, 2);
        categories.insert(categories.end(), expose, 3);
        categories.insert(categories.end(), resist, 4);

        std::vector<int> draw = draw_n_from_v(categories, total_hosts_moved, generator);
        int infected_moved = std::count(draw.begin(), draw.end(), 1);
        int susceptible_moved = std::count(draw.begin(), draw.end(), 2);
        int exposed_moved = std::count(draw.begin(), draw.end(), 3);
        int resistant_moved = std::count(draw.begin(), draw.end(), 4);

        if (exposed_moved > 0) {
            std::vector<int> exposed_draw = draw_n_from_cohorts(
                exposed_, exposed_moved, row_from, col_from, generator);
            int index = 0;
            for (auto& raster : exposed_) {
                raster(row_from, col_from) -= exposed_draw[index];
                raster(row_to, col_to) += exposed_draw[index];
                index += 1;
            }
        }
        if (infected_moved > 0) {
            std::vector<int> mortality_draw = draw_n_from_cohorts(
                mortality_tracker_vector_,
                infected_moved,
                row_from,
                col_from,
                generator);
            int index = 0;
            for (auto& raster : mortality_tracker_vector_) {
                raster(row_from, col_from) -= mortality_draw[index];
                raster(row_to, col_to) += mortality_draw[index];
                index += 1;
            }
        }
        // Ensure that the target cell of host movement is in suitable cells.
        // Since suitable cells originally comes from the total hosts, check first total
        // hosts and proceed only if there was no host.
        if (total_hosts_(row_to, col_to) == 0) {
            std::vector<int> new_index = {row_to, col_to};
            if (!container_contains(suitable_cells_, new_index)) {
                suitable_cells_.push_back(new_index);
            }
        }

        infected_(row_from, col_from) -= infected_moved;
        susceptible_(row_from, col_from) -= susceptible_moved;
        total_hosts_(row_from, col_from) -= total_hosts_moved;
        total_exposed_(row_from, col_from) -= exposed_moved;
        resistant_(row_from, col_from) -= resistant_moved;
        infected_(row_to, col_to) += infected_moved;
        susceptible_(row_to, col_to) += susceptible_moved;
        total_hosts_(row_to, col_to) += total_hosts_moved;
        total_exposed_(row_to, col_to) += exposed_moved;
        resistant_(row_to, col_to) += resistant_moved;

        // Returned total hosts actually moved is based only on the total host and no
        // other checks are performed. This assumes that the counts are correct in the
        // object (precondition).
        return total_hosts_moved;
    }

    /**
     * @brief Completely remove any hosts
     *
     * Removes hosts completely (as opposed to moving them to another pool).
     * If mortality is not active, the *mortality* parameter is ignored.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param susceptible Number of susceptible hosts to remove.
     * @param exposed Number of exposed hosts to remove by cohort.
     * @param infected Number of infected hosts to remove.
     * @param mortality Number of infected hosts in each mortality cohort.
     *
     * @note This does not remove resistant just like the original implementation in
     * treatments.
     */
    void completely_remove_hosts_at(
        RasterIndex row,
        RasterIndex col,
        int susceptible,
        std::vector<int> exposed,
        int infected,
        const std::vector<int>& mortality)
    {
        if (susceptible > 0)
            susceptible_(row, col) = susceptible_(row, col) - susceptible;

        if (exposed.size() != exposed_.size()) {
            throw std::invalid_argument(
                "counts is not the same size as the internal list of exposed ("
                + std::to_string(exposed.size())
                + " != " + std::to_string(exposed_.size()) + ") for cell ("
                + std::to_string(row) + ", " + std::to_string(col) + ")");
        }

        // no simple zip in C++, falling back to indices
        for (size_t i = 0; i < exposed.size(); ++i) {
            exposed_[i](row, col) -= exposed[i];
        }

        // Possibly reuse in the I->S removal.
        if (infected <= 0)
            return;
        if (!mortality_tracker_vector_.size()) {
            infected_(row, col) -= infected;
            reset_total_host(row, col);
            return;
        }
        if (mortality_tracker_vector_.size() != mortality.size()) {
            throw std::invalid_argument(
                "mortality is not the same size as the internal mortality tracker ("
                + std::to_string(mortality_tracker_vector_.size())
                + " != " + std::to_string(mortality.size()) + ") for cell ("
                + std::to_string(row) + ", " + std::to_string(col) + ")");
        }

        int mortality_total = 0;
        for (size_t i = 0; i < mortality_tracker_vector_.size(); ++i) {
            if (mortality_tracker_vector_[i](row, col) < mortality[i]) {
                throw std::invalid_argument(
                    "Mortality value [" + std::to_string(i) + "] is too high ("
                    + std::to_string(mortality[i]) + " > "
                    + std::to_string(mortality_tracker_vector_[i](row, col))
                    + ") for cell (" + std::to_string(row) + ", " + std::to_string(col)
                    + ")");
            }
            mortality_tracker_vector_[i](row, col) =
                mortality_tracker_vector_[i](row, col) - mortality[i];
            mortality_total += mortality[i];
        }
        // These two values will only match if we actually compute one from another
        // and once we don't need to keep the exact same double to int results for
        // tests. First condition always fails the tests. The second one may potentially
        // fail.
        if (infected != mortality_total) {
            throw std::invalid_argument(
                "Total of removed mortality values differs from removed infected "
                "count ("
                + std::to_string(mortality_total) + " != " + std::to_string(infected)
                + ") for cell (" + std::to_string(row) + ", " + std::to_string(col)
                + ")");
        }
        if (infected_(row, col) < mortality_total) {
            throw std::invalid_argument(
                "Total of removed mortality values is higher than current number "
                "of infected hosts ("
                + std::to_string(mortality_total) + " > " + std::to_string(infected)
                + ") for cell (" + std::to_string(row) + ", " + std::to_string(col)
                + ")");
        }
        infected_(row, col) -= infected;
        reset_total_host(row, col);
    }

    /**
     * @brief Remove infected hosts and make the hosts susceptible
     *
     * Distribution of removed infected hosts among mortality groups is stochastic.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param count Number of host to remove infection from
     * @param generator Random number generator to provide stochasticity for mortality
     *
     * @note Using this method either assumes the SI model or it needs to be used
     * together with remove_exposed_at() to handle SEI model.
     */
    void remove_infected_at(
        RasterIndex row, RasterIndex col, int count, Generator& generator)
    {
        // remove percentage of infestation/infection in the infected class
        infected_(row, col) -= count;
        // remove the removed infected from mortality cohorts
        if (count > 0) {
            std::vector<int> mortality_draw = draw_n_from_cohorts(
                mortality_tracker_vector_, count, row, col, generator);
            int index = 0;
            for (auto& raster : mortality_tracker_vector_) {
                raster(row, col) -= mortality_draw[index];
                index += 1;
            }
        }
        // move infested/infected host back to susceptible pool
        susceptible_(row, col) += count;
    }

    /**
     * @brief Make all infected hosts susceptible at the given cell
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param generator Random number generator to provide stochasticity for mortality
     */
    void remove_all_infected_at(RasterIndex row, RasterIndex col, Generator& generator)
    {
        auto count = this->infected_at(row, col);
        this->remove_infected_at(row, col, count, generator);
    }

    /**
     * @brief Remove percentage of infestation/infection (I->S)
     *
     * Besides removing percentage of infected, it removes the same percentage for total
     * exposed and remove individuals randomly from each cohort.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param ratio Ratio of removed infection
     * @param generator Random number generator to provide stochasticity for mortality
     */
    void remove_infection_by_ratio_at(
        RasterIndex row, RasterIndex col, double ratio, Generator& generator)
    {
        auto infected = this->infected_at(row, col);
        int removed_infected = infected - std::lround(infected * ratio);
        this->remove_infected_at(row, col, removed_infected, generator);
        auto exposed = this->exposed_at(row, col);
        int total_removed_exposed = exposed - std::lround(exposed * ratio);
        this->remove_exposed_at(row, col, total_removed_exposed, generator);
    }

    /**
     * @brief Remove exposed hosts and make the hosts susceptible
     *
     * Distribution of removed hosts among exposed groups is stochastic.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param count Number of host to remove infection from
     * @param generator Random number generator to provide stochasticity for mortality
     */
    void
    remove_exposed_at(RasterIndex row, RasterIndex col, int count, Generator& generator)
    {
        // remove the same percentage for total exposed and remove randomly from
        // each cohort
        total_exposed_(row, col) -= count;
        if (count > 0) {
            std::vector<int> exposed_draw =
                draw_n_from_cohorts(exposed_, count, row, col, generator);
            int index = 0;
            for (auto& raster : exposed_) {
                raster(row, col) -= exposed_draw[index];
                index += 1;
            }
        }
        // move infested/infected host back to susceptible pool
        susceptible_(row, col) += count;
    }

    /**
     * @brief Make hosts resistant in a given cell
     *
     * If mortality is not active, the *mortality* parameter is ignored.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param susceptible Number of susceptible hosts to make resistant
     * @param exposed Number of exposed hosts in each cohort to make resistant
     * @param infected Number of infected hosts to make resistant
     * @param mortality Number of infected hosts in each mortality cohort to make
     * resistant
     *
     * @throws std::invalid_argument if counts exceed the maximum amounts possible given
     * the current state (currently checked only for susceptible)
     * @throws std::invalid_argument if sizes don't equal to number of cohorts
     */
    void make_resistant_at(
        RasterIndex row,
        RasterIndex col,
        int susceptible,
        const std::vector<int>& exposed,
        int infected,
        const std::vector<int>& mortality)
    {
        int total_resistant = 0;

        if (susceptible_(row, col) < susceptible) {
            throw std::invalid_argument(
                "Total of newly resistant is higher than current number ("
                + std::to_string(susceptible) + " > "
                + std::to_string(susceptible_(row, col)) + ") for cell ("
                + std::to_string(row) + ", " + std::to_string(col) + ")");
        }

        susceptible_(row, col) -= susceptible;
        total_resistant += susceptible;

        if (exposed.size() != exposed_.size()) {
            throw std::invalid_argument(
                "exposed is not the same size as the internal list of exposed ("
                + std::to_string(exposed.size())
                + " != " + std::to_string(exposed_.size()) + ") for cell ("
                + std::to_string(row) + ", " + std::to_string(col) + ")");
        }
        // no simple zip in C++, falling back to indices
        for (size_t i = 0; i < exposed.size(); ++i) {
            exposed_[i](row, col) -= exposed[i];
            total_resistant += exposed[i];
        }
        infected_(row, col) -= infected;
        total_resistant += infected;
        resistant_(row, col) += total_resistant;
        if (!mortality_tracker_vector_.size()) {
            reset_total_host(row, col);
            return;
        }
        if (mortality_tracker_vector_.size() != mortality.size()) {
            throw std::invalid_argument(
                "mortality is not the same size as the internal mortality tracker ("
                + std::to_string(mortality_tracker_vector_.size())
                + " != " + std::to_string(mortality.size()) + ") for cell ("
                + std::to_string(row) + ", " + std::to_string(col) + ")");
        }
        int mortality_total = 0;
        // no simple zip in C++, falling back to indices
        for (size_t i = 0; i < mortality_tracker_vector_.size(); ++i) {
            mortality_tracker_vector_[i](row, col) -= mortality[i];
            mortality_total += mortality[i];
        }
        // These two values will only match if we actually compute one from another
        // and once we don't need to keep the exact same double to int results for
        // tests. First condition always fails the tests. The second one may potentially
        // fail.
        if (false && infected != mortality_total) {
            throw std::invalid_argument(
                "Total of mortality values differs from formerly infected, now resistant "
                "count ("
                + std::to_string(mortality_total) + " != " + std::to_string(infected)
                + " for cell (" + std::to_string(row) + ", " + std::to_string(col)
                + "))");
        }
        reset_total_host(row, col);
    }

    /**
     * @brief Remove resistance of host at a cell
     *
     * All resistant hosts are moved to the susceptible group.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     */
    void remove_resistance_at(RasterIndex row, RasterIndex col)
    {
        susceptible_(row, col) += resistant_(row, col);
        resistant_(row, col) = 0;
    }

    /**
     * @brief Apply mortality at a given cell
     *
     * Each cohort in exposed to mortality rate except those inside the mortality time
     * lag. In indexes of the cohort vector that are in the mortality_time_lag, no
     * mortality occurs. In the last year of mortality tracking, the all remaining
     * tracked infected hosts are removed. In all other indexes the number of tracked
     * individuals is multiplied by the mortality rate to calculate the number of hosts
     * that die that time step.
     *
     * If mortality rate is zero (<=0), no mortality is applied and mortality tracker
     * vector stays as is, i.e., no hosts die.
     *
     * To be used together with step_forward_mortality().
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param mortality_rate Percent of infected hosts that die each time period
     * @param mortality_time_lag Time lag prior to mortality beginning
     */
    void apply_mortality_at(
        RasterIndex row, RasterIndex col, double mortality_rate, int mortality_time_lag)
    {
        if (mortality_rate <= 0)
            return;
        int max_index = mortality_tracker_vector_.size() - mortality_time_lag - 1;
        for (int index = 0; index <= max_index; index++) {
            int mortality_in_index = 0;
            if (mortality_tracker_vector_[index](row, col) > 0) {
                // used to ensure that all infected hosts in the last year of
                // tracking mortality
                if (index == 0) {
                    mortality_in_index = mortality_tracker_vector_[index](row, col);
                }
                else {
                    mortality_in_index = static_cast<int>(std::floor(
                        mortality_rate * mortality_tracker_vector_[index](row, col)));
                }
                mortality_tracker_vector_[index](row, col) -= mortality_in_index;
                died_(row, col) += mortality_in_index;
                if (mortality_in_index > infected_(row, col)) {
                    throw std::runtime_error(
                        "Mortality[" + std::to_string(index)
                        + "] is higher than current number of infected hosts ("
                        + std::to_string(mortality_in_index) + " > "
                        + std::to_string(infected_(row, col)) + ") for cell ("
                        + std::to_string(row) + ", " + std::to_string(col) + ")");
                }
                if (mortality_in_index > total_hosts_(row, col)) {
                    throw std::runtime_error(
                        "Mortality[" + std::to_string(index)
                        + "] is higher than current number of total hosts ("
                        + std::to_string(mortality_in_index) + " > "
                        + std::to_string(total_hosts_(row, col)) + ") for cell ("
                        + std::to_string(row) + ", " + std::to_string(col) + ")");
                }
                if (infected_(row, col) > 0) {
                    infected_(row, col) -= mortality_in_index;
                }
                if (total_hosts_(row, col) > 0) {
                    total_hosts_(row, col) -= mortality_in_index;
                }
            }
        }
    }

    /**
     * @brief Apply mortality at a given cell
     *
     * Uses pest-host table for mortality parameters.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     *
     * @see apply_mortality_at(RasterIndex, RasterIndex, double, int)
     */
    void apply_mortality_at(RasterIndex row, RasterIndex col)
    {
        if (!pest_host_table_) {
            throw std::invalid_argument(
                "Set pest-host table before calling apply_mortality_at "
                "or provide parameters in the function call");
        }
        this->apply_mortality_at(
            row,
            col,
            pest_host_table_->mortality_rate(this),
            pest_host_table_->mortality_time_lag(this));
    }

    /**
     * @brief Make a step forward in tracking mortality cohorts
     *
     * This makes the mortality cohorts age.
     *
     * To be used together with apply_mortality_at().
     */
    void step_forward_mortality()
    {
        rotate_left_by_one(mortality_tracker_vector_);
    }

    /**
     * @brief Get number of infected hosts at a given cell
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     *
     * @return Number of infected hosts
     */
    int infected_at(RasterIndex row, RasterIndex col) const
    {
        return infected_(row, col);
    }

    /**
     * @brief Get number of susceptible hosts at a given cell
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     *
     * @return Number of susceptible hosts
     */
    int susceptible_at(RasterIndex row, RasterIndex col) const
    {
        return susceptible_(row, col);
    }

    /**
     * @brief Get number of exposed hosts at a given cell
     *
     * The number is the total number of exposed hosts in all exposed cohorts.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     *
     * @return Total number of exposed hosts
     *
     * @note The number is taken from the internally managed total. This will be merged
     * with computed_exposed_at() in the future.
     */
    int exposed_at(RasterIndex row, RasterIndex col) const
    {
        // Future code could remove total exposed and compute that on the fly.
        return total_exposed_(row, col);
    }

    /**
     * @brief Get number of resistant hosts at a given cell
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     *
     * @return Number of resistant hosts
     *
     * @note The number is computed from all cohorts. This will be merged with
     * exposed_at() in the future.
     */
    int computed_exposed_at(RasterIndex row, RasterIndex col) const
    {
        int sum = 0;
        for (const auto& raster : exposed_)
            sum += raster(row, col);
        return sum;
    }

    /**
     * @brief Get exposed hosts in each cohort at a given cell
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     *
     * @return Vector with number of exposed hosts per cohort
     */
    std::vector<int> exposed_by_group_at(RasterIndex row, RasterIndex col) const
    {
        std::vector<int> all;
        all.reserve(exposed_.size());
        for (const auto& raster : exposed_)
            all.push_back(raster(row, col));
        return all;
    }

    /**
     * @brief Get infected hosts in each mortality cohort at a given cell
     *
     * If mortality is not active, it returns number of all infected individuals
     * in the first and only item of the vector.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     *
     * @return Vector with number of infected hosts per mortality cohort
     */
    std::vector<int> mortality_by_group_at(RasterIndex row, RasterIndex col) const
    {
        std::vector<int> all;

        if (!mortality_tracker_vector_.size()) {
            all.push_back(infected_at(row, col));
            return all;
        }

        all.reserve(mortality_tracker_vector_.size());
        for (const auto& raster : mortality_tracker_vector_)
            all.push_back(raster(row, col));
        return all;
    }

    /**
     * @brief Get number of resistant hosts at a given cell
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     *
     * @return Number of resistant hosts
     */
    int resistant_at(RasterIndex row, RasterIndex col) const
    {
        return resistant_(row, col);
    }

    /**
     * @copydoc pops::HostPoolInterface::total_hosts_at()
     *
     * @note Computes the total host from susceptible and infected and does not consider
     * exposed and resistant.
     *
     * @note Computes the value on the fly and does not use the raster storage for total
     * host.
     */
    int total_hosts_at(RasterIndex row, RasterIndex col) const override
    {
        return susceptible_at(row, col) + infected_at(row, col);
    }

    /**
     * @brief Return true if cell is outside of the rectangle covered by hosts
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     *
     * @return true if cell is outside, false otherwise
     */
    bool is_outside(RasterIndex row, RasterIndex col)
    {
        return row < 0 || row >= rows_ || col < 0 || col >= cols_;
    }

    /** Infect exposed hosts (E to I transition in the SEI model)
     *
     * Applicable to SEI model, no-operation otherwise, i.e., the state
     * is left intact for SI.
     *
     * Before the first latency period is over,
     * the E to I transition won't happen because no item in the exposed
     * vector is old enough to become infected.
     *
     * The position of the items in the exposed vector determines their
     * age, i.e., for how long the hosts are exposed. The oldest item
     * is at the front and youngest at the end.
     * Before the the first latency period is over, items in the front
     * are still empty (unused) because no hosts were exposed for the
     * given time period.
     * After the first latency
     * period, this needs to be true before the function is called and
     * it is true after the function
     * finished with the difference that after the function is called,
     * the last item is empty in the sense that it does not contain any
     * hosts.
     *
     * When the E to I transition happens, hosts from the oldest item
     * in the exposed vector are moved to the infected (and mortality
     * tracker). They are removed from the exposed item and this item
     * is moved to the back of the vector.
     *
     * Like in disperse(), there is no distinction between *infected*
     * and *mortality_tracker*, but different usage is expected outside
     * of this function.
     *
     * The raster class used with the simulation class needs to support
     * `.fill(value)` method for this function to work.
     *
     * Step is used to evaluate the latency period.
     *
     * @param step Step in the simulation (>=0)
     */
    void step_forward(unsigned step)
    {
        if (model_type_ == ModelType::SusceptibleExposedInfected) {
            if (step >= latency_period_) {
                // Oldest item needs to be in the front
                auto& oldest = exposed_.front();
                // Move hosts to infected raster
                infected_ += oldest;
                mortality_tracker_vector_.back() += oldest;
                total_exposed_ += (oldest * (-1));
                // Reset the raster
                // (hosts moved from the raster)
                oldest.fill(0);
            }
            // Age the items and the used one to the back
            // elements go one position to the left
            // new oldest goes to the front
            // old oldest goes to the back
            rotate_left_by_one(exposed_);
        }
        else if (model_type_ == ModelType::SusceptibleInfected) {
            // no-op
        }
        else {
            throw std::runtime_error(
                "Unknown ModelType value in Simulation::infect_exposed()");
        }
    }

    /**
     * @brief Get suitable cells index
     *
     * Suitable cells are all cells which need to be modified when all cells with host
     * need to be modified.
     *
     * @return List of cell indices
     */
    const std::vector<std::vector<int>>& suitable_cells() const
    {
        return suitable_cells_;
    }

    /**
     * @brief Get list which contains this host pool
     *
     * This is for compatibility with multi-host pool. In case of this host pool, it
     * returns one item which is a pointer to this host pool.
     *
     * @return Read-write reference to a vector of size 1.
     */
    std::vector<HostPool*>& host_pools()
    {
        return host_pools_;
    }

private:
    std::vector<HostPool*> host_pools_ = {this};  // non-owning
    /**
     * @brief Reset total host value using the individual pools.
     *
     * This considers susceptible, exposed, infected, and resistant.
     *
     * This function is needed when the user needs total host as a result and everything
     * is provided through individual rasters (otherwise the simplest implementation
     * would just compute it on the fly always.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     */
    void reset_total_host(RasterIndex row, RasterIndex col)
    {
        total_hosts_(row, col) = susceptible_(row, col) + computed_exposed_at(row, col)
                                 + infected_(row, col) + resistant_(row, col);
    }

    IntegerRaster& susceptible_;
    IntegerRaster& infected_;

    std::vector<IntegerRaster>& exposed_;
    unsigned latency_period_{0};
    IntegerRaster& total_exposed_;

    IntegerRaster& resistant_;

    std::vector<IntegerRaster>& mortality_tracker_vector_;
    IntegerRaster& died_;

    IntegerRaster& total_hosts_;
    const Environment& environment_;

    ModelType model_type_;

    bool dispersers_stochasticity_{false};
    double reproductive_rate_{0};
    bool establishment_stochasticity_{true};
    double deterministic_establishment_probability_{0};

    /** pest-host table */
    const PestHostTable<HostPool>* pest_host_table_{nullptr};
    /** Competency table */
    const CompetencyTable<HostPool>* competency_table_{nullptr};

    RasterIndex rows_{0};
    RasterIndex cols_{0};

    std::vector<std::vector<int>>& suitable_cells_;
};

}  // namespace pops

#endif  // POPS_HOST_POOL_HPP
