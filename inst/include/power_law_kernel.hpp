/*
 * PoPS model - power law dispersal kernel
 *
 * Copyright (C) 2015-2020 by the authors.
 *
 * Authors: Margaret Lawrimore malawrim ncsu edu
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#ifndef POPS_POWER_LAW_KERNEL_HPP
#define POPS_POWER_LAW_KERNEL_HPP

#include "kernel_types.hpp"
#include "raster.hpp"

#include <random>

namespace pops {

using std::pow;

/*! Dispersal kernel for power law distribution
 *  class utilized by RadialKernel and DeterministicKernel
 */
class PowerLawKernel
{
protected:
    double alpha;
    double xmin;
    std::uniform_real_distribution<double> distribution;

public:
    PowerLawKernel(double a, double xm) : alpha(a), xmin(xm), distribution(0.0, 1.0)
    {
        //         if (alpha < 1.0) {
        //             throw std::invalid_argument("alpha must be greater than or equal
        //             to 1.0");
        //         }
        if (xmin == 0) {
            throw std::invalid_argument("xmin cannot equal 0.0");
        }
    }

    /*!
     *  Returns random value from power law distribution
     *  Used by RadialKernel to determine location of spread
     *  @param generator uniform random number generator
     *  @return value from power law distribution
     */
    template<class Generator>
    double random(Generator& generator)
    {
        double x = distribution(generator);
        return icdf(x);
    }

    /*!
     *  Power law probability density function
     *  Used by DeterministicKernel to determine location of spread
     *  @param x point within same space of distribution
     *  @return relative likelihood that a random variable would equal x
     *  @note only works with alpha < 1 so center square in matrix is always zero
     */
    double pdf(double x)
    {
        if (x < 0) {
            throw std::invalid_argument("x cannot be less than 0.0");
        }
        // x = distance to center and is always 0 for the center cell
        // so have to add xmin to account for minimum x values
        x = x + xmin;
        if (x < xmin) {
            throw std::invalid_argument("x must be greater than or equal to xmin");
        }
        return ((alpha - 1.0) / xmin) * pow(x / xmin, -alpha);
    }

    /*!
     *  Power law inverse cumulative distribution (quantile) function
     *  Used by DeterministicKernel to determine maximum distance of spread
     *  @param x proportion of the distribution
     *  @return value in distribution that is less than or equal to probability (x)
     */
    double icdf(double x)
    {
        if (x <= 0 || x >= 1) {
            throw std::invalid_argument("icdf: x must be between 0 and 1.0");
        }
        return pow(x / xmin, (-alpha + 1.0));
    }
};

}  // namespace pops

#endif  // POPS_POWER_LAW_KERNEL_HPP
