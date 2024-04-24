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

#ifndef POPS_PEST_POOL_HPP
#define POPS_PEST_POOL_HPP

#include <tuple>
#include <vector>

namespace pops {

/**
 * A class to track pests and manage pest which are not on hosts.
 *
 * This class can be viewed as a pests pool where pests can be
 * when they are not on hosts. However, it also tracks other information
 * about pests which is not captured by hosts.
 */
template<typename IntegerRaster, typename FloatRaster, typename RasterIndex>
class PestPool
{
public:
    /**
     * @brief Create an object with linked data
     * @param dispersers Generated dispersers
     * @param established_dispersers Established dispersers from a cell
     * @param outside_dispersers Dispersers which left the area
     *
     * The object provides method-based access to the underlying rasters.
     *
     * The object does not copy or take ownership of the objects passed in the
     * constructor.
     */
    PestPool(
        IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        std::vector<std::tuple<int, int>>& outside_dispersers)
        : dispersers_(dispersers),
          established_dispersers_(established_dispersers),
          outside_dispersers_(outside_dispersers)
    {}
    /**
     * @brief Set number of dispersers
     * @param row Row number
     * @param col Column number
     * @param value The new value
     */
    void set_dispersers_at(RasterIndex row, RasterIndex col, int value)
    {
        dispersers_(row, col) = value;
    }
    /**
     * @brief Return number of dispersers
     * @param row Row number
     * @param col Column number
     * @return The current value
     */
    int dispersers_at(RasterIndex row, RasterIndex col) const
    {
        return dispersers_(row, col);
    }
    /**
     * @brief Get raster of dispersers
     *
     * @note Most interaction should be done by dispersers_at().
     *
     * @return Read-only reference to dispersers raster
     */
    const IntegerRaster& dispersers()
    {
        return dispersers_;
    }
    /**
     * @brief Set number of established dispersers
     *
     * Established are dispersers which left cell (row, col) and established themselves
     * elsewhere, i.e., origin of the established dispersers is tracked.
     *
     * @param row Row number
     * @param col Column number
     * @param value The new value
     */
    void set_established_dispersers_at(RasterIndex row, RasterIndex col, int value)
    {
        established_dispersers_(row, col) = value;
    }
    // TODO: The following function should not be necessary because pests can't
    // un-establish. It exists just because it mirrors how the raster was handled in the
    // original Simulation code.
    /**
     * @brief Remove established dispersers
     * @param row Row number
     * @param col Column number
     * @param count How many dispers to remove
     */
    void remove_established_dispersers_at(RasterIndex row, RasterIndex col, int count)
    {
        established_dispersers_(row, col) -= count;
    }
    /**
     * @brief Add a disperser which left the study area
     * @param row Row number
     * @param col Column number
     */
    void add_outside_disperser_at(RasterIndex row, RasterIndex col)
    {
        outside_dispersers_.emplace_back(row, col);
    }
    /**
     * @brief Add a dispersers which left the study area
     * @param row Row number
     * @param col Column number
     * @param count Number of dispersers which left
     */
    void add_outside_dispersers_at(RasterIndex row, RasterIndex col, int count)
    {
        outside_dispersers_.reserve(outside_dispersers_.size() + count);
        for (int pest = 0; pest < count; ++pest)
            outside_dispersers_.emplace_back(row, col);
    }

private:
    /// Generated dispersers
    IntegerRaster& dispersers_;
    /// Origins of established dispersers
    IntegerRaster& established_dispersers_;
    /// Destinations of dispersers which left
    std::vector<std::tuple<int, int>>& outside_dispersers_;
};

}  // namespace pops

#endif  // POPS_PEST_POOL_HPP
