/*
 * PoPS model - natural dispersal kernel
 *
 * Copyright (C) 2019-2021 by the authors.
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

#ifndef POPS_NATURAL_KERNEL_HPP
#define POPS_NATURAL_KERNEL_HPP

#include "radial_kernel.hpp"
#include "deterministic_kernel.hpp"
#include "uniform_kernel.hpp"
#include "neighbor_kernel.hpp"
#include "kernel_types.hpp"
#include "kernel_base.hpp"
#include "config.hpp"

#include <memory>

namespace pops {

/**
 * @brief Create natural kernel from configuration
 *
 * Kernel parameters are taken from the configuration.
 *
 * @param config Configuration for the kernel
 * @param dispersers The disperser raster (reference, for deterministic kernel)
 *
 * @return Created kernel
 */
template<typename Generator, typename IntegerRaster, typename RasterIndex>
std::unique_ptr<KernelInterface<Generator>>
create_natural_kernel(const Config& config, const IntegerRaster& dispersers)
{
    auto natural_kernel = kernel_type_from_string(config.natural_kernel_type);
    if (natural_kernel == DispersalKernelType::Uniform) {
        using Kernel = DynamicWrapperKernel<UniformDispersalKernel, Generator>;
        // This can be std::make_unique in C++14.
        return std::unique_ptr<Kernel>(new Kernel(config.rows, config.cols));
    }
    else if (natural_kernel == DispersalKernelType::DeterministicNeighbor) {
        using Kernel =
            DynamicWrapperKernel<DeterministicNeighborDispersalKernel, Generator>;
        return std::unique_ptr<Kernel>(
            new Kernel(direction_from_string(config.natural_direction)));
    }
    else if (!config.dispersal_stochasticity) {
        using Kernel = DynamicWrapperKernel<
            DeterministicDispersalKernel<IntegerRaster>,
            Generator>;
        return std::unique_ptr<Kernel>(new Kernel(
            natural_kernel,
            dispersers,
            config.dispersal_percentage,
            config.ew_res,
            config.ns_res,
            config.natural_scale,
            config.shape));
    }
    else {
        using Kernel =
            DynamicWrapperKernel<RadialDispersalKernel<IntegerRaster>, Generator>;
        return std::unique_ptr<Kernel>(new Kernel(
            config.ew_res,
            config.ns_res,
            natural_kernel,
            config.natural_scale,
            direction_from_string(config.natural_direction),
            config.natural_kappa,
            config.shape));
    }
}

}  // namespace pops

#endif  // POPS_NATURAL_KERNEL_HPP
