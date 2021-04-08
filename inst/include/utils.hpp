#ifndef POPS_UTILS_HPP
#define POPS_UTILS_HPP

/*
 * PoPS model - general utility functions (unrelated to the model)
 *
 * Copyright (C) 2020 by the authors.
 *
 * Authors: Vaclav Petras <wenzeslaus gmail com>
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

/*!
 * Macro to mark unused variables (including parameters) and silence the warning
 * while documenting that it is intentionally unused.
 *
 * It is recommended to also document why the variable is unused, but left in the code.
 *
 * Usage:
 *
 * ```
 * UNUSED(variable_name);  // Parameter needed for backwards compatibility.
 * ```
 *
 * To be replaced by `[[maybe_unused]]` once we migrate to C++17 or higher.
 */
#define UNUSED(expr) (void)(expr)

#define M_PI 3.14159265358979323846
#define PI M_PI

#include <array>

typedef std::tuple<int, int, int, int> BBoxInt;
typedef std::tuple<double, double, double, double> BBoxFloat;
typedef std::tuple<bool, bool, bool, bool> BBoxBool;

/*! Spread direction
 *
 * Spread, typically wind, direction.
 * Values are in degrees and are used in computations.
 * `None` means that there is no wind.
 */
enum class Direction
{
    N = 0,  //!< North
    NE = 45,  //!< Northeast
    E = 90,  //!< East
    SE = 135,  //!< Southeast
    S = 180,  //!< South
    SW = 225,  //!< Southwest
    W = 270,  //!< West
    NW = 315,  //!< Northwest
    None  //!< No direction (non-directional)
};

template<int... Indices>
struct indices
{
    using next = indices<Indices..., sizeof...(Indices)>;
};

template<int Size>
struct build_indices
{
    using type = typename build_indices<Size - 1>::type::next;
};

template<>
struct build_indices<0>
{
    using type = indices<>;
};

template<typename T>
using Bare = typename std::remove_cv<typename std::remove_reference<T>::type>::type;

template<typename Tuple>
constexpr typename build_indices<std::tuple_size<Bare<Tuple>>::value>::type
make_indices()
{
    return {};
}

template<typename Tuple, int... Indices>
std::array<
    typename std::tuple_element<0, Bare<Tuple>>::type,
    std::tuple_size<Bare<Tuple>>::value>
to_array(Tuple&& tuple, indices<Indices...>)
{
    using std::get;
    return {{get<Indices>(std::forward<Tuple>(tuple))...}};
}

template<typename Tuple>
auto to_array(Tuple&& tuple)
    -> decltype(to_array(std::declval<Tuple>(), make_indices<Tuple>()))
{
    return to_array(std::forward<Tuple>(tuple), make_indices<Tuple>());
}

std::string quarantine_enum_to_string(Direction type)
{
    switch (type) {
    case Direction::N:
        return "N";
    case Direction::S:
        return "S";
    case Direction::E:
        return "E";
    case Direction::W:
        return "W";
    case Direction::None:
        return "None";
    default:
        return "Invalid direction";
    }
}

#endif  // POPS_UTILS_HPP
