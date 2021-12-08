/*
 * PoPS model - disperal kernels
 *
 * Copyright (C) 2021 by the authors.
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

#ifndef KERNEL_BASE_HPP
#define KERNEL_BASE_HPP

#include "kernel_types.hpp"

namespace pops {

/**
 * Interface which all dynamically used kernels need to inherit from.
 */
template<typename Generator>
class KernelInterface
{
public:
    /*! \copydoc RadialDispersalKernel::operator()()
     */
    virtual std::tuple<int, int> operator()(Generator& generator, int row, int col) = 0;

    virtual bool is_cell_eligible(int row, int col) = 0;

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
    virtual bool supports_kernel(const DispersalKernelType type) = 0;

    virtual ~KernelInterface() = default;
};

/**
 * A class which turns any kernel into a kernel usable with KernelInterface.
 *
 * Kernel needs to implement operator() and is_cell_eligible function.
 */
template<typename ActualKernel, typename Generator>
class DynamicWrapperKernel : public KernelInterface<Generator>
{
public:
    /** Creates an interal copy of the provided kernel */
    DynamicWrapperKernel(const ActualKernel& kernel) : kernel_(kernel) {}

    /** Creates a new internal kernel instance from the provided parameters */
    template<typename... Args>
    DynamicWrapperKernel(Args&&... args) : kernel_(std::forward<Args>(args)...)
    {}

    /*! \copydoc RadialDispersalKernel::operator()()
     */
    std::tuple<int, int> operator()(Generator& generator, int row, int col) override
    {
        return kernel_.operator()(generator, row, col);
    }

    /*! \copydoc RadialDispersalKernel::is_cell_eligible()
     */
    bool is_cell_eligible(int row, int col) override
    {
        return DynamicWrapperKernel<ActualKernel, Generator>::kernel_.is_cell_eligible(
            row, col);
    }

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
    bool supports_kernel(const DispersalKernelType type) override
    {
        return ActualKernel::supports_kernel(type);
    }

protected:
    ActualKernel kernel_;
};

}  // namespace pops

#endif  // KERNEL_BASE_HPP
