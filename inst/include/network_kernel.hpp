/*
 * PoPS model - network dispersal kernel
 *
 * Copyright (C) 2020-2021 by the authors.
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

#ifndef POPS_NETWORK_KERNEL_HPP
#define POPS_NETWORK_KERNEL_HPP

#include "kernel_types.hpp"
#include "network.hpp"

namespace pops {

/*!
 * @brief Dispersal kernel for dispersal over a network.
 *
 * Network node must be present in the cell to start traveling, so is_cell_eligible()
 * needs to be called first to see if the kernel can be used with the given row and
 * column.
 */
template<typename RasterIndex>
class NetworkDispersalKernel
{
public:
    /**
     * @brief Create kernel which travels through a network.
     *
     * The kernel assumes that the *network* is already initialized. It does not modify
     * the network.
     *
     * The *min_distance* and *max_distance* parameters are used as a range for uniform
     * real distribution which determines the travel distance (cost) through the network
     * for one trip if the network movement is walking (and not teleporting).
     *
     * @param network Existing network
     * @param min_distance Minimum travel distance (cost)
     * @param max_distance Maximum travel distance (cost)
     * @param jump End always on a node (snaps result to the closest node)
     */
    NetworkDispersalKernel(
        const Network<RasterIndex>& network,
        double min_distance,
        double max_distance,
        bool jump = false)
        : network_(network),
          distance_distribution_(min_distance, max_distance),
          jump_(jump)
    {}

    /**
     * @brief Create kernel which teleports from one node to another.
     *
     * The kernel assumes that the *network* is already initialized. It does not modify
     * the network.
     *
     * @param network Existing network
     */
    NetworkDispersalKernel(const Network<RasterIndex>& network)
        : network_(network), teleport_{true}
    {}

    /*! \copybrief RadialDispersalKernel::operator()()
     *
     * @throws std::invalid_argument if there is no network node at a given *row*
     * and *col* (can be checked with is_cell_eligible() beforehand)
     */
    template<typename Generator>
    std::tuple<int, int> operator()(Generator& generator, int row, int col)
    {
        if (teleport_) {
            return network_.teleport(row, col, generator);
        }
        double distance = distance_distribution_(generator);
        std::tie(row, col) = network_.walk(row, col, distance, generator);

        return std::make_tuple(row, col);
    }

    /**
     * @brief Test if cell is eligible to be used with the kernel.
     *
     * @param row Row to be used with the kernel
     * @param col Column to be used with the kernel
     * @return true if cell can be used, false otherwise
     */
    bool is_cell_eligible(int row, int col)
    {
        return network_.has_node_at(row, col);
    }

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
    static bool supports_kernel(const DispersalKernelType type)
    {
        return type == DispersalKernelType::Network;
    }

protected:
    /** Reference to the network */
    const Network<RasterIndex>& network_;
    /** Travel distance (cost) distribution */
    std::uniform_real_distribution<double> distance_distribution_;
    /** Step through network instead of traveling between nodes */
    bool teleport_{false};
    /** Snap to nodes when traveling between nodes */
    bool jump_{false};
};

}  // namespace pops

#endif  // POPS_NETWORK_KERNEL_HPP
