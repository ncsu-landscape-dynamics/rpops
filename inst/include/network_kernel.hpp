/*
 * PoPS model - network dispersal kernel
 *
 * Copyright (C) 2020-2025 by the authors.
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
#include "multi_network.hpp"

namespace pops {

/*!
 * @brief Dispersal kernel for dispersal over a network.
 *
 * Network node must be present in the cell to start traveling, so is_cell_eligible()
 * needs to be called first to see if the kernel can be used with the given row and
 * column.
 */
template<typename NetworkType>
class NetworkDispersalKernel
{
public:
    /**
     * @brief Create kernel which travels through a network.
     *
     * The kernel assumes that the *network* is already initialized. It does not modify
     * the network.
     *
     * @param network Existing network
     */
    NetworkDispersalKernel(const NetworkType& network) : network_(network) {}

    /*! \copybrief RadialDispersalKernel::operator()()
     *
     * @throws std::invalid_argument if there is no network node at a given *row*
     * and *col* (can be checked with is_cell_eligible() beforehand)
     */
    template<typename Generator>
    std::tuple<int, int> operator()(Generator& generator, int row, int col)
    {
        return network_.move(row, col, generator);
    }

    /**
     * @brief Test if cell is eligible to be used with the kernel.
     *
     * Tests if node is present at a cell in the underlying network.
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
    const NetworkType& network_;
};

}  // namespace pops

#endif  // POPS_NETWORK_KERNEL_HPP
