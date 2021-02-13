/*
 * PoPS model - exponential dispersal kernel
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

#ifndef POPS_EXPONENTIAL_KERNEL_HPP
#define POPS_EXPONENTIAL_KERNEL_HPP

#include "kernel_types.hpp"
#include "raster.hpp"
#include "utils.hpp"

#include <random>

namespace pops {

using std::exp;
using std::log;

/*! Dispersal kernel for exponential distribution
 *  class utilized by RadialKernel and DeterministicKernel
 */
class ExponentialKernel
{
protected:
    double beta;
    std::exponential_distribution<double> exponential_distribution;

public:
    ExponentialKernel(double b) : beta(b), exponential_distribution(1.0 / beta)
    {
        if (beta <= 0) {
            throw std::invalid_argument("beta must be greater than 0.0");
        }
    }

    /*!
     *  Returns random value from exponential distribution
     *  Used by RadialKernel to determine location of spread
     *  @param generator uniform random number generator
     *  @return value from exponential distribution
     */
    template<class Generator>
    double random(Generator& generator)
    {
        return std::abs(exponential_distribution(generator));
    }

    /*!
     *  Exponential probability density function
     *  Used by DeterministicKernel to determine location of spread
     *  @param x point within same space of distribution
     *  @return relative likelihood that a random variable would equal x
     */
    double pdf(double x)
    {
        if (x < 0) {
            throw std::invalid_argument("x must be greater than or equal to zero");
        }
        return (1.0 / beta) * (exp(-x / beta));
    }

    /*!
     *  Exponential inverse cumulative distribution (quantile) function
     *  Used by DeterministicKernel to determine maximum distance of spread
     *  @param x proportion of the distribution
     *  @return value in distribution that is less than or equal to probability (x)
     */
    double icdf(double x)
    {
        if (x <= 0 || x >= 1) {
            throw std::invalid_argument("icdf: x must be between 0.0 and 1.0");
        }
        else {
            return -beta * log(1 - x);
        }
    }
};

}  // namespace pops

#endif  // POPS_EXPONENTIAL_KERNEL_HPP
