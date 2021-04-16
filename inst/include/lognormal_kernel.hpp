/*
 * PoPS model - log normal dispersal kernel
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

#ifndef POPS_LOGNORMAL_KERNEL_HPP
#define POPS_LOGNORMAL_KERNEL_HPP

#include "kernel_types.hpp"
#include "raster.hpp"
#include "utils.hpp"

#include <random>

namespace pops {

using std::pow;
using std::exp;
using std::log;
using std::sqrt;

/*! Dispersal kernel for log normal distribution
 *  class utilized by RadialKernel and DeterministicKernel
 */
class LogNormalKernel
{
protected:
    double sigma;
    std::lognormal_distribution<double> lognormal_distribution;

public:
    LogNormalKernel(double s) : sigma(s), lognormal_distribution(0.0, sigma)
    {
        if (sigma <= 0) {
            throw std::invalid_argument("sigma must be greater than 0.0");
        }
    }

    /*!
     *  Returns random value from log normal distribution
     *  Used by RadialKernel to determine location of spread
     *  @param generator uniform random number generator
     *  @return value from log normal distribution
     */
    template<class Generator>
    double random(Generator& generator)
    {
        return std::abs(lognormal_distribution(generator));
    }

    /*!
     *  Log normal probability density function
     *  Used by DeterministicKernel to determine location of spread
     *  @param x point within same space of distribution
     *  @return relative likelihood that a random variable would equal x
     */
    double pdf(double x)
    {
        if (x < 0) {
            throw std::invalid_argument("x must be greater than or equal to 0.0");
        }
        if (x == 0) {
            return 0;
        }
        return (1 / (x * sigma * sqrt(2 * M_PI)))
               * exp(-(pow(log(x), 2)) / (2 * pow(sigma, 2)));
    }

    /*!
     *  Log normal inverse cumulative distribution (quantile) function
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
        // float b = 2 / (M_PI * 0.147) + 0.5f * log(1 - pow(y, 2));
        // double inverf =
        //    (sign * sqrt(-b + sqrt(pow(b, 2) - (1 / (0.147) * log(1 - pow(y, 2))))));
        double a = 0.140012;
        double t = 2.0 / (M_PI * a);
        double l = log(1 - pow(y, 2));
        double inverf =
            sign * sqrt(sqrt(pow(t + (l / 2.0), 2) - (l / a)) - (t + (l / 2.0)));
        return exp(sqrt(2 * pow(sigma, 2)) * inverf);
    }
};

}  // namespace pops

#endif  // POPS_LOGNORMAL_KERNEL_HPP
