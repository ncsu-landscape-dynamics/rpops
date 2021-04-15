/*
 * PoPS model - logistic dispersal kernel
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

#ifndef POPS_LOGISTIC_KERNEL_HPP
#define POPS_LOGISTIC_KERNEL_HPP

#include "kernel_types.hpp"
#include "raster.hpp"

#include <random>

namespace pops {

using std::pow;
using std::exp;
using std::log;
/*! Dispersal kernel for logistic distribution
 *  class utilized by RadialKernel and DeterministicKernel
 */
class LogisticKernel
{
protected:
    double s;
    std::uniform_real_distribution<double> distribution;

public:
    LogisticKernel(double scale) : s(scale), distribution(0.0, 1.0)
    {
        if (s <= 0) {
            throw std::invalid_argument(
                "LogisticKernel: scale (s) must be greater than 0.0");
        }
    }

    /*!
     *  Returns random value from logistic distribution
     *  Used by RadialKernel to determine location of spread
     *  @param generator uniform random number generator
     *  @return value from logistic distribution
     */
    template<class Generator>
    double random(Generator& generator)
    {
        double x = distribution(generator);
        return icdf(x);
    }

    /*!
     *  Logistic probability density function
     *  Used by DeterministicKernel to determine location of spread
     *  @param x point within same space of distribution
     *  @return relative likelihood that a random variable would equal x
     */
    double pdf(double x)
    {
        if (s == 1) {
            return exp(-x) / pow(1 + exp(-x), 2);
        }
        return (exp(-x / s)) / (s * pow(1 + exp(-x / s), 2));
    }

    /*!
     *  Logistic inverse cumulative distribution (quantile) function
     *  Used by DeterministicKernel to determine maximum distance of spread
     *  @param x proportion of the distribution
     *  @return value in distribution that is less than or equal to probability (x)
     */
    double icdf(double x)
    {
        if (x <= 0 || x >= 1) {
            throw std::invalid_argument("icdf: x must be between 0.0 and 1.0");
        }
        return s * log(x / (1.0 - x));
    }
};

}  // namespace pops

#endif  // POPS_LOGISTIC_KERNEL_HPP
