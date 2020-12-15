/*
 * PoPS model - deterministic dispersal kernel
 *
 * Copyright (C) 2015-2020 by the authors.
 *
 * Authors: Margaret Lawrimore (malawrim ncsu edu)
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#ifndef POPS_DETERMINISTIC_KERNEL_HPP
#define POPS_DETERMINISTIC_KERNEL_HPP

#include <vector>
#include <tuple>

#include "raster.hpp"
#include "kernel_types.hpp"
#include "utils.hpp"
#include "hyperbolic_secant_kernel.hpp"
#include "logistic_kernel.hpp"
#include "exponential_power_kernel.hpp"
#include "exponential_kernel.hpp"
#include "cauchy_kernel.hpp"
#include "gamma_kernel.hpp"
#include "lognormal_kernel.hpp"
#include "normal_kernel.hpp"
#include "weibull_kernel.hpp"
#include "power_law_kernel.hpp"

namespace pops {

using std::pow;
using std::tan;
using std::exp;
using std::log;
using std::ceil;
using std::abs;
using std::sqrt;

/*!
 * Dispersal kernel for deterministic spread to cell with highest probability of
 * spread
 *
 * Dispersal Kernel type determines use of Exponential or Cauchy distribution
 * to find probability.
 *
 * dispersal_percentage is the percent of all possible dispersal to be included
 * in the moving window size (e.g for 99% input 0.99).
 *
 * Useful for testing as it is deterministic and provides fully replicable results
 */
template<typename IntegerRaster>
class DeterministicDispersalKernel
{
protected:
    const IntegerRaster& dispersers_;
    // row/col position of middle cell
    int mid_row = 0;
    int mid_col = 0;
    // position of cell from previous call
    int prev_row = -1;
    int prev_col = -1;
    // number of rows/cols in the probability window
    int number_of_rows = 0;
    int number_of_columns = 0;
    // maximum distance from center cell to outer cells
    double max_distance{0};
    Raster<double> probability;
    Raster<double> probability_copy;
    CauchyKernel cauchy;
    ExponentialKernel exponential;
    WeibullKernel weibull;
    LogNormalKernel log_normal;
    NormalKernel normal;
    HyperbolicSecantKernel hyperbolic_secant;
    PowerLawKernel power_law;
    LogisticKernel logistic;
    GammaKernel gamma;
    ExponentialPowerKernel exponential_power;

    DispersalKernelType kernel_type_;
    double proportion_of_dispersers;
    // the west-east resolution of the pixel
    double east_west_resolution;
    // the north-south resolution of the pixel
    double north_south_resolution;

public:
    DeterministicDispersalKernel(
        DispersalKernelType dispersal_kernel,
        const IntegerRaster& dispersers,
        double dispersal_percentage,
        double ew_res,
        double ns_res,
        double distance_scale,
        double shape = 1.0)
        : dispersers_(dispersers),
          cauchy(distance_scale),
          exponential(distance_scale),
          weibull(distance_scale, shape),
          log_normal(distance_scale),
          normal(distance_scale),
          hyperbolic_secant(distance_scale),
          power_law(distance_scale, shape),
          logistic(distance_scale),
          gamma(distance_scale, shape),
          exponential_power(distance_scale, shape),
          kernel_type_(dispersal_kernel),
          east_west_resolution(ew_res),
          north_south_resolution(ns_res)
    {
        // We initialize max distance only for the supported kernels.
        // For the others, we report the error only when really called
        // to allow use of this class in initialization phase.
        if (kernel_type_ == DispersalKernelType::Cauchy) {
            max_distance = cauchy.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::Exponential) {
            max_distance = exponential.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::Weibull) {
            max_distance = weibull.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::Normal) {
            max_distance = normal.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::LogNormal) {
            max_distance = log_normal.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::HyperbolicSecant) {
            max_distance = hyperbolic_secant.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::PowerLaw) {
            max_distance = power_law.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::Logistic) {
            max_distance = logistic.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::Gamma) {
            max_distance = gamma.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::ExponentialPower) {
            max_distance = exponential_power.icdf(dispersal_percentage);
        }
        number_of_columns = ceil(max_distance / east_west_resolution) * 2 + 1;
        number_of_rows = ceil(max_distance / north_south_resolution) * 2 + 1;
        Raster<double> prob_size(number_of_rows, number_of_columns, 0);
        probability = prob_size;
        probability_copy = prob_size;
        mid_row = number_of_rows / 2;
        mid_col = number_of_columns / 2;
        double sum = 0.0;
        for (int i = 0; i < number_of_rows; i++) {
            for (int j = 0; j < number_of_columns; j++) {
                double distance_to_center = std::sqrt(
                    pow((abs(mid_row - i) * east_west_resolution), 2)
                    + pow((abs(mid_col - j) * north_south_resolution), 2));
                // determine probability based on distance
                if (kernel_type_ == DispersalKernelType::Cauchy) {
                    probability(i, j) = abs(cauchy.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::Exponential) {
                    probability(i, j) = abs(exponential.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::Weibull) {
                    probability(i, j) = abs(weibull.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::Normal) {
                    probability(i, j) = abs(normal.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::LogNormal) {
                    probability(i, j) = abs(log_normal.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::PowerLaw) {
                    probability(i, j) = abs(power_law.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::HyperbolicSecant) {
                    probability(i, j) = abs(hyperbolic_secant.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::Logistic) {
                    probability(i, j) = abs(logistic.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::Gamma) {
                    probability(i, j) = abs(gamma.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::ExponentialPower) {
                    probability(i, j) = abs(exponential_power.pdf(distance_to_center));
                }
                sum += probability(i, j);
            }
        }
        // normalize based on the sum of all probabilities in the raster
        probability /= sum;
    }

    /*! Generates a new position for the spread.
     *
     *  Creates a copy of the probability matrix to mark where dispersers are
     * assigned. New window created any time a new cell is selected from
     * simulation.disperse
     *
     *  Selects next row/col value based on the cell with the highest probability
     *  in the window.
     *
     */
    template<class Generator>
    std::tuple<int, int> operator()(Generator& generator, int row, int col)
    {
        UNUSED(generator);  // Deterministic does not need random numbers.
        if (kernel_type_ != DispersalKernelType::Cauchy
            && kernel_type_ != DispersalKernelType::Exponential
            && kernel_type_ != DispersalKernelType::Weibull
            && kernel_type_ != DispersalKernelType::Normal
            && kernel_type_ != DispersalKernelType::LogNormal
            && kernel_type_ != DispersalKernelType::HyperbolicSecant
            && kernel_type_ != DispersalKernelType::PowerLaw
            && kernel_type_ != DispersalKernelType::Logistic
            && kernel_type_ != DispersalKernelType::Gamma
            && kernel_type_ != DispersalKernelType::ExponentialPower) {
            throw std::invalid_argument(
                "DeterministicDispersalKernel: Unsupported dispersal kernel type");
        }
        // reset the window if considering a new cell
        if (row != prev_row || col != prev_col) {
            proportion_of_dispersers = 1.0 / (double)dispersers_(row, col);
            probability_copy = probability;
        }

        int row_movement = 0;
        int col_movement = 0;

        double max = (double)-std::numeric_limits<int>::max();
        int max_prob_row = 0;
        int max_prob_col = 0;

        // find cell with highest probability
        for (int i = 0; i < number_of_rows; i++) {
            for (int j = 0; j < number_of_columns; j++) {
                if (probability_copy(i, j) > max) {
                    max = probability_copy(i, j);
                    max_prob_row = i;
                    max_prob_col = j;
                    row_movement = i - mid_row;
                    col_movement = j - mid_col;
                }
            }
        }

        // subtracting 1/number of dispersers ensures we always move the same
        // proportion of the individuals to each cell no matter how many are
        // dispersing
        probability_copy(max_prob_row, max_prob_col) -= proportion_of_dispersers;
        prev_row = row;
        prev_col = col;

        // return values in terms of actual location
        return std::make_tuple(row + row_movement, col + col_movement);
    }

    /*! Returns true if the kernel class support a given kernel type
     *
     * \warning This function is experimental and may be removed or
     * changed at any time.
     */
    static bool supports_kernel(const DispersalKernelType type)
    {
        static const std::array<DispersalKernelType, 10> supports = {
            DispersalKernelType::Cauchy,
            DispersalKernelType::Exponential,
            DispersalKernelType::Weibull,
            DispersalKernelType::Normal,
            DispersalKernelType::LogNormal,
            DispersalKernelType::PowerLaw,
            DispersalKernelType::HyperbolicSecant,
            DispersalKernelType::Gamma,
            DispersalKernelType::ExponentialPower,
            DispersalKernelType::Logistic};
        auto it = std::find(supports.cbegin(), supports.cend(), type);
        return it != supports.cend();
    }
};

}  // namespace pops

#endif  // POPS_DETERMINISTIC_KERNEL_HPP
