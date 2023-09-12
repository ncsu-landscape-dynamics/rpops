/*
 * PoPS model - Normal distribution with uniform distribution fallback
 *
 * Copyright (C) 2023 by the authors.
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
#ifndef NORMAL_DISTRIBUTION_WITH_UNIFORM_FALLBACK_HPP
#define NORMAL_DISTRIBUTION_WITH_UNIFORM_FALLBACK_HPP

#include <random>

namespace pops {

/**
 * Normal distribution with uniform distribution fallback when out of range
 *
 * This is a combination of normal distribution and a uniform distribution.
 * If value from the normal distribution is out of range, a uniform distibution is
 * used to supply a new value which is guaranteed to be in the specified interval.
 * The resulting PDF has the shape of the original normal distribution but is clipped
 * (or trimmed) according to the range and is shifted up in the direction of y axis.
 * Unlike truncated normal distribution, this spreads the trunkated
 * values equally in the interval.
 */
template<typename Real>
class NormalDistributionWithUniformFallback
{
public:
    /**
     * @brief Create so-called frozen distribution object.
     * @param mean mean of the normal distribution
     * @param stddev standard deviation of the normal distribution
     * @param low lower limit of the interval
     * @param high upper limit of the interval
     *
     * The interval is closed on both sides, i.e., [low, high] (both values can be
     * included in the output). The implemenation currently relies on a common bug
     * across compilers, but it can be adjusted in the future.
     */
    NormalDistributionWithUniformFallback(Real mean, Real stddev, Real low, Real high)
        : low_(low),
          high_(high),
          normal_distribution_(mean, stddev),
          uniform_distribution_(low, high)
    {}
    /**
     * Get a random value from the distribution.
     */
    template<class Generator>
    Real operator()(Generator& generator)
    {
        Real value = normal_distribution_(generator);
        if (value < low_ || value > high_) {
            // Value is out of range, get a random value in range.
            return uniform_distribution_(generator);
        }
        return value;
    }

private:
    Real low_;
    Real high_;
    std::normal_distribution<Real> normal_distribution_;
    /**
     * While the standard specifies that the range is [a, b), most
     * implementations return [a, b] as of February 2023 which is
     * exactly what we desire (and doing [a, b] using std::nextafter
     * would cause out-of-range for us).
     */
    std::uniform_real_distribution<Real> uniform_distribution_;
};

}  // namespace pops

#endif  // NORMAL_DISTRIBUTION_WITH_UNIFORM_FALLBACK_HPP
