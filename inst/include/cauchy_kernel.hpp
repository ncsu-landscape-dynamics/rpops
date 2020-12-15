/*
 * PoPS model - cauchy dispersal kernel
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

#ifndef POPS_CAUCHY_KERNEL_HPP
#define POPS_CAUCHY_KERNEL_HPP

#include "utils.hpp"
#include "kernel_types.hpp"
#include "raster.hpp"

#include <random>

namespace pops {

using std::pow;
using std::tan;

/*! Dispersal kernel for cauchy distribution
 *  class utilized by RadialKernel and DeterministicKernel
 */
class CauchyKernel
{
protected:
    double s;
    std::cauchy_distribution<double> cauchy_distribution;

public:
    CauchyKernel(double scale) : s(scale), cauchy_distribution(0, s)
    {
        if (scale <= 0) {
            throw std::invalid_argument("scale (s) must be greater than 0.0");
        }
    }

    /*!
     *  Returns random value from cauchy distribution
     *  Used by RadialKernel to determine location of spread
     *  @param generator uniform random number generator
     *  @return value from cauchy distribution
     */
    template<class Generator>
    double random(Generator& generator)
    {
        return std::abs(cauchy_distribution(generator));
    }

    /*!
     *  Cauchy probability density function
     *  Used by DeterministicKernel to determine location of spread
     *  @param x point within same space of distribution
     *  @return relative likelihood that a random variable would equal x
     */
    double pdf(double x)
    {
        return 1 / ((s * M_PI) * (1 + (pow(x / s, 2))));
    }

    /*!
     *  Cauchy inverse cumulative distribution (quantile) function
     *  Used by DeterministicKernel to determine maximum distance of spread
     *  @param x proportion of the distribution
     *  @return value in distribution that is less than or equal to probability (x)
     */
    double icdf(double x)
    {
        if (x <= 0 || x >= 1) {
            throw std::invalid_argument("icdf: x must be between 0.0 and 1.0");
        }
        return s * tan(M_PI * (x - 0.5));
    }
};

}  // namespace pops

#endif  // POPS_CAUCHY_KERNEL_HPP
