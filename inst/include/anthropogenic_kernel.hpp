/*
 * PoPS model - anthropogenic dispersal kernel
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

#ifndef POPS_ANTHROPOGENIC_KERNEL_HPP
#define POPS_ANTHROPOGENIC_KERNEL_HPP

#include "radial_kernel.hpp"
#include "deterministic_kernel.hpp"
#include "uniform_kernel.hpp"
#include "neighbor_kernel.hpp"
#include "network_kernel.hpp"
#include "kernel_types.hpp"
#include "kernel_base.hpp"
#include "config.hpp"

#include <memory>

namespace pops {

/**
 * @brief Create anthropogenic kernel from configuration
 *
 * Same structure as the natural kernel, but the parameters are for anthropogenic
 * kernel when available.
 *
 * @param dispersers The disperser raster (reference, for deterministic kernel)
 * @param network Network (initialized or not)
 * @return Created kernel
 */
template<typename Generator, typename IntegerRaster, typename RasterIndex>
std::unique_ptr<KernelInterface<Generator>> create_anthro_kernel(
    const Config& config,
    const IntegerRaster& dispersers,
    const Network<RasterIndex>& network)
{
    auto anthro_kernel = kernel_type_from_string(config.anthro_kernel_type);
    if (anthro_kernel == DispersalKernelType::Uniform) {
        using Kernel = DynamicWrapperKernel<UniformDispersalKernel, Generator>;
        return std::unique_ptr<Kernel>(new Kernel(config.rows, config.cols));
    }
    else if (anthro_kernel == DispersalKernelType::DeterministicNeighbor) {
        using Kernel =
            DynamicWrapperKernel<DeterministicNeighborDispersalKernel, Generator>;
        return std::unique_ptr<Kernel>(
            new Kernel(direction_from_string(config.anthro_direction)));
    }
    else if (anthro_kernel == DispersalKernelType::Network) {
        using Kernel =
            DynamicWrapperKernel<NetworkDispersalKernel<RasterIndex>, Generator>;
        if (config.network_movement == "teleport")
            return std::unique_ptr<Kernel>(new Kernel(network));
        bool jump = config.network_movement == "jump" ? true : false;
        return std::unique_ptr<Kernel>(new Kernel(
            network, config.network_min_distance, config.network_max_distance, jump));
    }
    else if (!config.dispersal_stochasticity) {
        using Kernel = DynamicWrapperKernel<
            DeterministicDispersalKernel<IntegerRaster>,
            Generator>;
        return std::unique_ptr<Kernel>(new Kernel(
            anthro_kernel,
            dispersers,
            config.dispersal_percentage,
            config.ew_res,
            config.ns_res,
            config.anthro_scale,
            config.shape));
    }
    else {
        using Kernel =
            DynamicWrapperKernel<RadialDispersalKernel<IntegerRaster>, Generator>;
        return std::unique_ptr<Kernel>(new Kernel(
            config.ew_res,
            config.ns_res,
            anthro_kernel,
            config.anthro_scale,
            direction_from_string(config.anthro_direction),
            config.anthro_kappa,
            config.shape));
    }
}

}  // namespace pops

#endif  // POPS_ANTHROPOGENIC_KERNEL_HPP
