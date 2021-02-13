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

#ifndef POPS_EXPONENTIAL_POWER_KERNEL_HPP
#define POPS_EXPONENTIAL_POWER_KERNEL_HPP

#include "kernel_types.hpp"
#include "gamma_kernel.hpp"
#include "raster.hpp"

#include <random>

namespace pops {

using std::pow;
using std::exp;

/*! Dispersal kernel for exponential power distribution
 *  class utilized by RadialKernel and DeterministicKernel
 */
class ExponentialPowerKernel
{
protected:
    double alpha;
    double beta;
    std::uniform_real_distribution<double> distribution;

public:
    ExponentialPowerKernel(double a, double b)
        : alpha(a), beta(b), distribution(0.0, 1.0)
    {
        if (alpha <= 0 || beta <= 0) {
            throw std::invalid_argument(
                "alpha and beta must greater than or equal to 0");
        }
    }

    /*!
     *  Returns random value from exponential power distribution
     *  Used by RadialKernel to determine location of spread
     *  @param generator uniform random number generator
     *  @return value from exponential power distribution
     */
    template<class Generator>
    double random(Generator& generator)
    {
        double x = distribution(generator);
        return icdf(x);
    }

    /*!
     *  Exponential Power probability density function
     *  Used by DeterministicKernel to determine location of spread
     *  @param x point within same space of distribution
     *  @return relative likelihood that a random variable would equal x
     */
    double pdf(double x)
    {
        return (beta / (2 * alpha * std::tgamma(1.0 / beta)))
               * pow(exp(-x / alpha), beta);
    }

    /*!
     *  Exponential Power inverse cumulative distribution (quantile) function
     *  Used by DeterministicKernel to determine maximum distance of spread
     *  @param x proportion of the distribution
     *  @return value in distribution that is less than or equal to probability (x)
     */
    double icdf(double x)
    {
        if (x <= 0 || x >= 1) {
            throw std::invalid_argument("icdf: x must be between 0.0 and 1.0");
        }
        GammaKernel gamma_distribution(1.0 / beta, 1.0 / pow(alpha, beta));
        double gamma = gamma_distribution.icdf(2 * std::abs(x - 0.5));
        return (x - 0.5) * pow(gamma, 1.0 / beta);
    }
};

}  // namespace pops

#endif  // POPS_EXPONENTIAL_POWER_KERNEL_HPP
