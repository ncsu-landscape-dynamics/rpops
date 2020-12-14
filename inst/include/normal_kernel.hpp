/*
 * PoPS model - normal dispersal kernel
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

#ifndef POPS_NORMAL_KERNEL_HPP
#define POPS_NORMAL_KERNEL_HPP

#include "kernel_types.hpp"
#include "raster.hpp"

#include <random>

namespace pops {

using std::pow;
using std::exp;
using std::log;
using std::sqrt;

/*! Dispersal kernel for normal distribution
 *  class utilized by RadialKernel and DeterministicKernel
 */
class NormalKernel
{
protected:
    double sigma;
    std::normal_distribution<double> normal_distribution;

public:
    NormalKernel(double s) : sigma(s), normal_distribution(0.0, sigma)
    {
        if (sigma == 0) {
            throw std::invalid_argument("sigma cannot be zero");
        }
    }

    /*!
     *  Returns random value from normal distribution
     *  Used by RadialKernel to determine location of spread
     *  @param generator uniform random number generator
     *  @return value from normal distribution
     */
    template<class Generator>
    double random(Generator& generator)
    {
        return std::abs(normal_distribution(generator));
    }

    /*!
     *  Normal probability density function
     *  Used by DeterministicKernel to determine location of spread
     *  @param x point within same space of distribution
     *  @return relative likelihood that a random variable would equal x
     */
    double pdf(double x)
    {
        if (sigma == 1) {
            return 1.0 / (sqrt(2 * M_PI)) * exp(-0.5 * pow(x, 2));
        }
        return 1.0 / (sigma * sqrt(2 * M_PI)) * exp(-0.5 * pow(x / sigma, 2));
    }

    /*!
     *  Normal inverse cumulative distribution (quantile) function
     *  Used by DeterministicKernel to determine maximum distance of spread
     *  @param x proportion of the distribution
     *  @return value in distribution that is less than or equal to probability (x)
     *  @note uses approximation for inverse error function from "A handy
     *    approximation for the error function and its inverse"  by Sergei Winitzki
     */
    double icdf(double x)
    {
        if (x <= 0 || x >= 1) {
            throw std::invalid_argument("icdf: x must be between 0.0 and 1.0");
        }
        //  approximation for inverse error function
        double y = (2 * x) - 1;
        float sign = (y < 0) ? -1.0f : 1.0f;
        //  0.147 used for a relative error of about 2*10^-3
        float b = 2.0 / (M_PI * 0.147) + 0.5f * log(1 - pow(y, 2));
        double inverf =
            (sign * sqrt(-b + sqrt(pow(b, 2) - (1.0 / (0.147) * log(1 - pow(y, 2))))));

        return sigma * std::sqrt(2) * inverf;
    }
};

}  // namespace pops

#endif  // POPS_NORMAL_KERNEL_HPP
