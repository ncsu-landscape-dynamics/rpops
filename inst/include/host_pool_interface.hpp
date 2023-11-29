/*
 * PoPS model - pest or pathogen spread simulation
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

#ifndef POPS_HOST_POOL_INTERFACE_HPP
#define POPS_HOST_POOL_INTERFACE_HPP

namespace pops {

/**
 * Interface declaration for a host pool.
 *
 * Currently, the interface is providing only total number of hosts because that's the
 * only function needed throughout the code.
 *
 * @note The usefulness of this interface will be evaluated later on, e.g., if we need
 * functionally different hosts or if we have dependencies which the interface would
 * address.
 */
template<typename RasterIndexType>
class HostPoolInterface
{
public:
    using RasterIndex = RasterIndexType;

    virtual ~HostPoolInterface() {}

    /**
     * @brief Get total number of hosts for a cell
     *
     * @param row Row index of the cell
     * @param col Column index the cell
     *
     * @return Number of hosts
     */
    virtual int total_hosts_at(RasterIndex row, RasterIndex col) const = 0;
};

}  // namespace pops

#endif  // POPS_HOST_POOL_INTERFACE_HPP
