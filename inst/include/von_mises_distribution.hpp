/*
 * PoPS model - radial disperal kernel
 *
 * Copyright (C) 2015-2020 by the authors.
 *
 * Authors: Vaclav Petras (wenzeslaus gmail com)
 *          Chris Jones (cjones1688 gmail com)
 *          Anna Petrasova (kratochanna gmail com)
 *          Zexi Chen (zchen22 ncsu edu)
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#ifndef POPS_VON_MISES_DISTRIBUTION_HPP
#define POPS_VON_MISES_DISTRIBUTION_HPP

#include "utils.hpp"

#include <cmath>
#include <random>

namespace pops {

// we need to bring sqrt to the namespace
// otherwise only candidates are raster-related
using std::sqrt;

/*! Von Mises Distribution (Circular data distribution)

    mu is the mean angle, expressed in radians between 0 and 2*pi,
    and kappa is the concentration parameter, which must be greater
    than or equal to zero. If kappa is equal to zero, this distribution
    reduces to a uniform random angle over the range 0 to 2*pi.
*/
class VonMisesDistribution
{
public:
    VonMisesDistribution(double mu, double kappa)
        : mu(mu), kappa(kappa), distribution(0.0, 1.0)
    {}
    template<class Generator>
    double operator()(Generator& generator)
    {
        double a, b, c, f, r, theta, u1, u2, u3, z;

        if (kappa <= 1.e-06)
            return 2 * PI * distribution(generator);

        a = 1.0 + sqrt(1.0 + 4.0 * kappa * kappa);
        b = (a - sqrt(2.0 * a)) / (2.0 * kappa);
        r = (1.0 + b * b) / (2.0 * b);

        while (true) {
            u1 = distribution(generator);
            z = cos(PI * u1);
            f = (1.0 + r * z) / (r + z);
            c = kappa * (r - f);
            u2 = distribution(generator);
            if (u2 <= c * (2.0 - c) || u2 < c * exp(1.0 - c))
                break;
        }

        u3 = distribution(generator);
        if (u3 > 0.5) {
            theta = fmod(mu + acos(f), 2 * PI);
        }
        else {
            theta = fmod(mu - acos(f), 2 * PI);
        }
        return theta;
    }

private:
    double mu;
    double kappa;
    std::uniform_real_distribution<double> distribution;
};

}  // namespace pops

#endif  // POPS_RADIAL_KERNEL_HPP
