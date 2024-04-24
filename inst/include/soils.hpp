/*
 * PoPS model - soils pool
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

#ifndef POPS_SOILS_HPP
#define POPS_SOILS_HPP

#include <cmath>
#include <memory>
#include <tuple>
#include <vector>
#include <random>
#include <string>
#include <stdexcept>

#include "utils.hpp"
#include "environment.hpp"

namespace pops {

/** Handles disperser (pathogen) storage in soils.
 *
 * Takes care of adding dispersers to the pool and of taking them out.
 */
template<
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex,
    typename GeneratorProvider>
class SoilPool
{
public:
    /**
     * @brief Create a fully functioning soil pool object.
     *
     * The size of the vector of provided rasters is number of simulation steps the
     * dispersers persist in the soil.
     *
     * @param rasters List of soil rasters (cannot be empty)
     * @param environment Surrounding environment (weather)
     * @param generate_stochasticity Use stochasticity when releasing from the pool
     * @param establishment_stochasticity Use stochasticity when adding to the pool
     * @param fixed_establishment_probability Non-stochastic establishment probability
     */
    SoilPool(
        std::vector<IntegerRaster>& rasters,
        const Environment<IntegerRaster, FloatRaster, RasterIndex, GeneratorProvider>&
            environment,
        bool generate_stochasticity = true,
        bool establishment_stochasticity = true,
        double fixed_establishment_probability = 0)
        : rasters_(&rasters),
          environment_(&environment),
          generate_stochasticity_(generate_stochasticity),
          establishment_stochasticity_(establishment_stochasticity),
          fixed_establishment_probability_(fixed_establishment_probability)
    {
        if (rasters.empty()) {
            throw std::logic_error(
                "List of rasters of SoilPool needs to have at least one item");
        }
    }

    /**
     * Release (generate) dispersers from the pool
     *
     * Dispersers are released from each cohort randomly.
     */
    template<typename Generator>
    int dispersers_from(RasterIndex row, RasterIndex col, Generator& generator)
    {
        auto count = this->total_at(row, col);
        double lambda = environment_->weather_coefficient_at(row, col);
        int dispersers = 0;
        if (this->generate_stochasticity_) {
            std::poisson_distribution<int> distribution(lambda);
            for (int k = 0; k < count; k++) {
                dispersers += distribution(generator);
            }
        }
        else {
            dispersers = lambda * count;
        }
        auto draw = draw_n_from_cohorts(*rasters_, dispersers, row, col, generator);
        size_t index = 0;
        for (auto count : draw) {
            (*rasters_)[index](row, col) -= count;
            ++index;
        }
        return dispersers;
    }

    /**
     * Move (add, store) one disperser to the pool
     *
     * Disperser may or many not establish in the soil.
     */
    template<typename Generator>
    void disperser_to(RasterIndex row, RasterIndex col, Generator& generator)
    {
        double current_probability = environment_->weather_coefficient_at(row, col);
        double tester;
        if (this->establishment_stochasticity_)
            tester = distribution_uniform_(generator);
        else
            tester = 1 - fixed_establishment_probability_;
        if (tester < current_probability) {
            this->add_at(row, col);
        }
    }

    /**
     * Add (store) dispersers to the pool
     *
     * See disperser_to().
     */
    template<typename Generator>
    void dispersers_to(
        int dispersers, RasterIndex row, RasterIndex col, Generator& generator)
    {
        for (int i = 0; i < dispersers; i++)
            this->disperser_to(row, col, generator);
    }

    /**
     * Get total number of dispersers at a specific location
     *
     * All cohorts for one cell are combined.
     */
    int total_at(RasterIndex row, RasterIndex col) const
    {
        int total = 0;
        for (auto& raster : *rasters_) {
            total += raster(row, col);
        }
        return total;
    }

    /**
     * Advance to the next simulation step
     *
     * Makes the cohorts age one time step. The oldest cohort disappears. The actual
     * time of the step is driven by how often this is called and the size of the soil
     * raster vector.
     *
     * Internally, this rotates the cohorts and clears what becomes the youngest cohort.
     */
    void next_step(int step)
    {
        UNUSED(step);
        rotate_left_by_one(*rasters_);
        rasters_->back().fill(0);
    }

protected:
    std::vector<IntegerRaster>* rasters_{nullptr};  ///< Disperser cohorts
    /**
     * Surrounding environment
     */
    const Environment<IntegerRaster, FloatRaster, RasterIndex, GeneratorProvider>*
        environment_{nullptr};
    bool generate_stochasticity_{false};  ///< Outgoing dispersers stochasticity
    bool establishment_stochasticity_{false};  ///< Incoming dispersers
    /**
     * Probability for establishment when stochastic is disabled.
     */
    double fixed_establishment_probability_{0};
    /**
     * Distribution driving stochastic establishment.
     */
    std::uniform_real_distribution<double> distribution_uniform_{0.0, 1.0};

    /**
     * Add disperser or dispersers at a specific place
     *
     * Dispersers are always added.
     */
    void add_at(RasterIndex row, RasterIndex col, int value = 1)
    {
        rasters_->back()(row, col) += value;
    }
};

}  // namespace pops

#endif  // POPS_SOILS_HPP
