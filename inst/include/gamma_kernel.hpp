/*
 * PoPS model - gamma dispersal kernel
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

#ifndef POPS_GAMMA_KERNEL_HPP
#define POPS_GAMMA_KERNEL_HPP

#include "kernel_types.hpp"
#include "raster.hpp"
#include "lognormal_kernel.hpp"

#include <random>

namespace pops {

using std::pow;
using std::exp;

/*! Dispersal kernel for gamma distribution
 *  class utilized by RadialKernel and DeterministicKernel
 */
class GammaKernel
{
protected:
    double alpha;
    double theta;
    std::gamma_distribution<double> gamma_distribution;

public:
    GammaKernel(double a, double t)
        : alpha(a), theta(t), gamma_distribution(alpha, 1.0 / theta)
    {
        if (alpha <= 0 || theta <= 0) {
            throw std::invalid_argument(
                "alpha and theta must greater than or equal to 0");
        }
    }

    /*!
     *  Returns random value from gamma distribution
     *  Used by RadialKernel to determine location of spread
     *  @param generator uniform random number generator
     *  @return value from gamma distribution
     */
    template<class Generator>
    double random(Generator& generator)
    {
        return std::abs(gamma_distribution(generator));
    }

    /*!
     *  Gamma probability density function
     *  Used by DeterministicKernel to determine location of spread
     *  @param x point within same space of distribution
     *  @return relative likelihood that a random variable would equal x
     */
    double pdf(double x)
    {
        if (x < 0) {
            throw std::invalid_argument("x must greater than or equal to 0");
        }
        return 1.0 / (std::tgamma(alpha) * pow(theta, alpha)) * pow(x, (alpha - 1))
               * exp(-x / theta);
    }

    /*!
     *  Gamma cumulative distribution function used by gamma icdf
     *  @param x value in distribution
     *  @return probability of x
     */
    double cdf(double x)
    {
        double sum = 0.0;
        double beta = 1.0 / theta;
        for (int i = 0; i < alpha; i++) {
            // tgamma = (i-1)! used since c++ has no factorial in std lib
            sum += (1.0 / std::tgamma(i + 1)) * exp(-beta * x) * pow(beta * x, i);
        }
        return 1 - sum;
    }

    /*!
     *  Gamma inverse cumulative distribution (quantile) function
     *  Used by DeterministicKernel to determine maximum distance of spread
     *  There is no known closed-form solution for Gamma icdf, used Newton's method
     *  References include: Abramowitz, M. and Stegun, I.A. (1964) Handbook
     *  of Mathematical Functions, Dover, New York, section 26.1.
     *  Evans, M., Hastings, N., and Peacock, B. (1993) Statistical Distributions,
     *  2nd ed., Wiley.
     *  @param x proportion of the distribution
     *  @return value in distribution that is less than or equal to probability (x)
     */
    double icdf(double x)
    {
        if (x <= 0 || x >= 1) {
            throw std::invalid_argument("icdf: x must be between 0.0 and 1.0");
        }
        // pick starting approximation using lognormal icdf
        LogNormalKernel lognormal(1);
        double guess = lognormal.icdf(x);
        double check = cdf(guess);
        double numiterations = 1000;  // will need to adjust this
        double precision = 0.001;  // will need to adjust this
        for (int i = 0; i < numiterations; i++) {
            if (check < (x - precision) || check > (x + precision)) {
                double dif = check - x;
                // if dif is positive guess is greater than needed
                // if dif is negative guess is less than needed
                double past_guess = guess;
                double derivative = dif / pdf(guess);
                // limit size of next guess
                guess = std::max(guess / 10, std::min(guess * 10, guess - derivative));
                check = cdf(guess);
                // Check if we went to far and need to backtrack
                int count = 0;
                bool run = true;
                while ((std::abs(dif) < std::abs(check - x)) && run) {
                    guess = (guess + past_guess) / 2.0;
                    check = cdf(guess);
                    count++;
                    if (count > 20) {
                        run = false;
                    }
                }
            }
            else {
                return guess;
            }
        }
        throw std::invalid_argument("unable to find solution to gamma icdf ");
        return -1;
    }
};

}  // namespace pops

#endif  // POPS_GAMMA_KERNEL_HPP
