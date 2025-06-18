/*
 * Tests for the PoPS Model class.
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

#ifndef POPS_MODEL_HPP
#define POPS_MODEL_HPP

#include "environment.hpp"
#include "config.hpp"
#include "treatments.hpp"
#include "spread_rate.hpp"
#include "host_pool.hpp"
#include "pest_pool.hpp"
#include "actions.hpp"
#include "multi_host_pool.hpp"
#include "switch_kernel.hpp"
#include "kernel.hpp"
#include "scheduling.hpp"
#include "quarantine.hpp"
#include "soils.hpp"
#include "generator_provider.hpp"

#include <vector>

namespace pops {

template<
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex,
    typename Generator = std::default_random_engine,
    typename KernelFactory = DispersalKernel<Generator>(
        const Config&, const IntegerRaster&, const Network<RasterIndex>&)>
class Model
{
protected:
    Config config_;
    RandomNumberGeneratorProvider<Generator> generator_provider_;
    DispersalKernelType natural_kernel;
    DispersalKernelType anthro_kernel;
    UniformDispersalKernel uniform_kernel;
    DeterministicNeighborDispersalKernel natural_neighbor_kernel;
    DeterministicNeighborDispersalKernel anthro_neighbor_kernel;
    KernelFactory& kernel_factory_;
    /**
     * Surrounding environment
     */
    Environment<
        IntegerRaster,
        FloatRaster,
        RasterIndex,
        RandomNumberGeneratorProvider<Generator>>
        environment_;
    /**
     * Optionally created soil pool
     */
    std::shared_ptr<SoilPool<
        IntegerRaster,
        FloatRaster,
        RasterIndex,
        RandomNumberGeneratorProvider<Generator>>>
        soil_pool_{nullptr};
    unsigned last_index{0};

    /**
     * @brief Create overpopulation movement kernel
     *
     * Same as the natural kernel. The natural kernel parameters are used,
     * but the scale for radial and deterministic kernel is multiplied
     * by the leaving scale coefficient.
     *
     * @param dispersers The disperser raster (reference, for deterministic kernel)
     * @param network Network (initialized or not)
     * @return Created kernel
     */
    SwitchDispersalKernel<IntegerRaster, RasterIndex>
    create_overpopulation_movement_kernel(
        const IntegerRaster& dispersers, const Network<RasterIndex>& network)
    {
        RadialDispersalKernel<IntegerRaster> radial_kernel(
            config_.ew_res,
            config_.ns_res,
            natural_kernel,
            config_.natural_scale * config_.leaving_scale_coefficient,
            direction_from_string(config_.natural_direction),
            config_.natural_kappa,
            config_.shape);
        DeterministicDispersalKernel<IntegerRaster> deterministic_kernel(
            natural_kernel,
            dispersers,
            config_.dispersal_percentage,
            config_.ew_res,
            config_.ns_res,
            config_.natural_scale * config_.leaving_scale_coefficient,
            config_.shape);
        NetworkDispersalKernel<RasterIndex> network_kernel(
            network, config_.network_min_distance, config_.network_max_distance);
        SwitchDispersalKernel<IntegerRaster, RasterIndex> selectable_kernel(
            natural_kernel,
            radial_kernel,
            deterministic_kernel,
            uniform_kernel,
            network_kernel,
            natural_neighbor_kernel,
            config_.dispersal_stochasticity);
        return selectable_kernel;
    }

public:
    /** Type for single-host pool */
    using StandardSingleHostPool = HostPool<
        IntegerRaster,
        FloatRaster,
        RasterIndex,
        RandomNumberGeneratorProvider<Generator>>;
    /** Type for multi-host pool */
    using StandardMultiHostPool = MultiHostPool<
        StandardSingleHostPool,
        IntegerRaster,
        FloatRaster,
        RasterIndex,
        RandomNumberGeneratorProvider<Generator>>;
    /** Type for pest pool */
    using StandardPestPool = PestPool<IntegerRaster, FloatRaster, RasterIndex>;

    Model(
        const Config& config,
        KernelFactory& kernel_factory =
            create_dynamic_kernel<Generator, IntegerRaster, RasterIndex>)
        : config_(config),
          generator_provider_(config),
          natural_kernel(kernel_type_from_string(config.natural_kernel_type)),
          anthro_kernel(kernel_type_from_string(config.anthro_kernel_type)),
          uniform_kernel(config.rows, config.cols),
          natural_neighbor_kernel(direction_from_string(config.natural_direction)),
          anthro_neighbor_kernel(direction_from_string(config.anthro_direction)),
          kernel_factory_(kernel_factory)
    {}

    /**
     * @brief Run one step of the simulation.
     *
     * The *total_populations* can be total number of hosts in the basic case
     * or it can be the total size of population of all relevant species
     * both host and non-host if dilution effect should be applied.
     * When movement is applied, *total_populations* needs to be only the
     * total number of hosts because movement does not support non-host
     * individuals.
     *
     * *dispersers* is for internal use and for tracking dispersers creation.
     * The input values are ignored and the output is not the current existing
     * dispersers, but only the number of dispersers generated (and subsequently
     * used) in this step. There are no dispersers in between simulation steps.
     *
     * @param step Step number in the simulation
     * @param[in,out] infected Infected hosts
     * @param[in,out] susceptible Susceptible hosts
     * @param[in,out] total_hosts All host individuals in the area. Is equal to
     * infected + exposed + susceptible in the cell
     * @param[in,out] total_populations All host and non-host individuals in the area
     * @param[out] dispersers Dispersing individuals (used internally)
     * @param[out] established_dispersers Dispersers originating from a given cell which
     * established elsewhere
     * @param[in,out] total_exposed Sum of all exposed hosts (if SEI model is active)
     * @param[in,out] exposed Exposed hosts (if SEI model is active)
     * @param[in,out] mortality_tracker Mortality tracker used to generate *died*.
     * Expectation is that mortality tracker is of length (1/mortality_rate +
     * mortality_time_lag)
     * @param[out] died Infected hosts which died this step based on the mortality
     * schedule
     * @param[in] temperatures Vector of temperatures used to evaluate lethal
     * temperature
     * @param[in] survival_rates Pest survival rates
     * @param[in,out] resistant Resistant hosts (host temporarily removed from
     * susceptible hosts)
     * @param[in,out] outside_dispersers Dispersers escaping the rasters (adds to the
     * vector)
     * @param[in,out] quarantine Quarantine escape tracker
     * @param[in] quarantine_areas Quarantine areas
     * @param[in] movements Table of host movements
     * @param network Network (initialized or Network::null_network() if unused)
     * @param[in,out] suitable_cells List of indices of cells with hosts
     *
     * @note The parameters roughly correspond to Simulation::disperse()
     * and Simulation::disperse_and_infect() functions, so these can be used
     * for further reference.
     */
    void run_step(
        int step,
        IntegerRaster& infected,
        IntegerRaster& susceptible,
        IntegerRaster& total_populations,
        IntegerRaster& total_hosts,
        IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        IntegerRaster& total_exposed,
        std::vector<IntegerRaster>& exposed,
        std::vector<IntegerRaster>& mortality_tracker,
        IntegerRaster& died,
        const std::vector<FloatRaster>& temperatures,
        const std::vector<FloatRaster>& survival_rates,
        IntegerRaster& resistant,
        std::vector<std::tuple<int, int>>& outside_dispersers,  // out
        QuarantineEscapeAction<IntegerRaster>& quarantine,  // out
        const IntegerRaster& quarantine_areas,
        const std::vector<std::vector<int>> movements,
        const Network<RasterIndex>& network,
        std::vector<std::vector<int>>& suitable_cells)
    {
        StandardSingleHostPool host_pool(
            config_,
            susceptible,
            exposed,
            infected,
            total_exposed,
            resistant,
            mortality_tracker,
            died,
            total_hosts,
            environment_,
            suitable_cells);
        std::vector<StandardSingleHostPool*> host_pools = {&host_pool};
        StandardMultiHostPool multi_host_pool(host_pools, config_);
        StandardPestPool pest_pool{
            dispersers, established_dispersers, outside_dispersers};
        SpreadRateAction<StandardMultiHostPool, RasterIndex> spread_rate(
            multi_host_pool,
            config_.rows,
            config_.cols,
            config_.ew_res,
            config_.ns_res,
            0);
        Treatments<StandardSingleHostPool, Raster<double>> treatments(
            config_.scheduler());
        run_step(
            step,
            multi_host_pool,
            pest_pool,
            total_populations,
            treatments,
            temperatures,
            survival_rates,
            spread_rate,
            quarantine,
            quarantine_areas,
            movements,
            network);
    }

    /**
     * @brief Run one step of the simulation.
     *
     * @param step Step number in the simulation
     * @param host_pool Host pool
     * @param pest_pool Pest pool
     * @param[in,out] total_populations All host and non-host individuals in the area
     * @param[in,out] treatments Treatments to be applied (also tracks use of
     * treatments)
     * @param[in] temperatures Vector of temperatures used to evaluate lethal
     * temperature
     * @param[in] survival_rates Pest survival rates
     * @param spread_rate Spread rate action
     * @param[in,out] quarantine Quarantine escape tracker
     * @param[in] quarantine_areas Quarantine areas
     * @param[in] movements Table of host movements
     * @param network Network (initialized or Network::null_network() if unused)
     */
    void run_step(
        int step,
        StandardMultiHostPool& host_pool,
        StandardPestPool& pest_pool,
        IntegerRaster& total_populations,
        Treatments<StandardSingleHostPool, FloatRaster>& treatments,
        const std::vector<FloatRaster>& temperatures,
        const std::vector<FloatRaster>& survival_rates,
        SpreadRateAction<StandardMultiHostPool, RasterIndex>& spread_rate,
        QuarantineEscapeAction<IntegerRaster>& quarantine,
        const IntegerRaster& quarantine_areas,
        const std::vector<std::vector<int>> movements,
        const Network<RasterIndex>& network)
    {
        // Soil step is the same as simulation step.
        if (soil_pool_)
            soil_pool_->next_step(step);
        // removal of dispersers due to lethal temperatures
        if (config_.use_lethal_temperature && config_.lethal_schedule()[step]) {
            int lethal_step =
                simulation_step_to_action_step(config_.lethal_schedule(), step);
            this->environment().update_temperature(temperatures[lethal_step]);
            RemoveByTemperature<
                StandardMultiHostPool,
                IntegerRaster,
                FloatRaster,
                RasterIndex,
                RandomNumberGeneratorProvider<Generator>>
                remove(this->environment(), config_.lethal_temperature);
            remove.action(host_pool, generator_provider_);
        }
        // removal of percentage of dispersers
        if (config_.use_survival_rate && config_.survival_rate_schedule()[step]) {
            int survival_step =
                simulation_step_to_action_step(config_.survival_rate_schedule(), step);
            SurvivalRateAction<StandardMultiHostPool, IntegerRaster, FloatRaster>
                survival(survival_rates[survival_step]);
            survival.action(host_pool, generator_provider_);
        }
        // actual spread
        if (config_.spread_schedule()[step]) {
            auto dispersal_kernel =
                kernel_factory_(config_, pest_pool.dispersers(), network);
            auto overpopulation_kernel =
                create_overpopulation_movement_kernel(pest_pool.dispersers(), network);
            SpreadAction<
                StandardMultiHostPool,
                StandardPestPool,
                IntegerRaster,
                FloatRaster,
                RasterIndex,
                decltype(dispersal_kernel),
                RandomNumberGeneratorProvider<Generator>>
                spread_action{dispersal_kernel};

            environment_.set_total_population(&total_populations);
            // Soils are activated by an independent function call for model, but spread
            // action is temporary, so it is activated for every step.
            if (this->soil_pool_) {
                spread_action.activate_soils(
                    soil_pool_, config_.dispersers_to_soils_percentage);
            }
            spread_action.action(host_pool, pest_pool, generator_provider_);
            host_pool.step_forward(step);
            if (config_.use_overpopulation_movements) {
                MoveOverpopulatedPests<
                    StandardMultiHostPool,
                    StandardPestPool,
                    IntegerRaster,
                    FloatRaster,
                    RasterIndex,
                    decltype(overpopulation_kernel)>
                    move_pest{
                        overpopulation_kernel,
                        config_.overpopulation_percentage,
                        config_.leaving_percentage,
                        config_.rows,
                        config_.cols};
                move_pest.action(host_pool, pest_pool, generator_provider_);
            }
            if (config_.use_movements) {
                HostMovement<
                    StandardMultiHostPool,
                    IntegerRaster,
                    FloatRaster,
                    RasterIndex>
                    host_movement{
                        static_cast<unsigned>(
                            step),  // Step is int, but indexing happens with unsigned.
                        last_index,
                        movements,
                        config_.movement_schedule};
                last_index = host_movement.action(host_pool, generator_provider_);
            }
        }
        // treatments
        if (config_.use_treatments) {
            for (auto& host : host_pool.host_pools()) {
                treatments.manage(step, *host);
            }
        }
        if (config_.use_mortality && config_.mortality_schedule()[step]) {
            // expectation is that mortality tracker is of length (1/mortality_rate
            // + mortality_time_lag).
            // TODO: died.zero(); should be done by the caller if needed, document!
            Mortality<StandardMultiHostPool, IntegerRaster, FloatRaster> mortality;
            mortality.action(host_pool);
        }
        // compute spread rate
        if (config_.use_spreadrates && config_.spread_rate_schedule()[step]) {
            unsigned rates_step =
                simulation_step_to_action_step(config_.spread_rate_schedule(), step);
            spread_rate.action(host_pool, rates_step);
        }
        // compute quarantine escape
        if (config_.use_quarantine && config_.quarantine_schedule()[step]) {
            unsigned action_step =
                simulation_step_to_action_step(config_.quarantine_schedule(), step);
            quarantine.action(host_pool, quarantine_areas, action_step);
        }
    }

    /**
     * @brief Get the associated random number generator provider
     * @return Reference to the generator provider
     */
    RandomNumberGeneratorProvider<Generator>& random_number_generator()
    {
        return generator_provider_;
    }

    /**
     * @brief Get surrounding environment
     * @return Environment object by reference
     */
    Environment<
        IntegerRaster,
        FloatRaster,
        RasterIndex,
        RandomNumberGeneratorProvider<Generator>>&
    environment()
    {
        return environment_;
    }

    /**
     * @brief Activate movement to and from soil pool
     *
     * The size of vector of rasters for cohorts is the number of simulation steps
     * dispersers stay in the soil.
     *
     * @param rasters Vector of rasters for cohorts
     */
    void activate_soils(std::vector<IntegerRaster>& rasters)
    {
        // The soil pool is created again for every new activation.
        this->soil_pool_.reset(new SoilPool<
                               IntegerRaster,
                               FloatRaster,
                               RasterIndex,
                               RandomNumberGeneratorProvider<Generator>>(
            rasters,
            this->environment_,
            config_.generate_stochasticity,
            config_.establishment_stochasticity));
    }
};

}  // namespace pops

#endif  // POPS_MODEL_HPP
