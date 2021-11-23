#ifndef POPS_UTILS_HPP
#define POPS_UTILS_HPP

/*
 * PoPS model - general utility functions (unrelated to the model)
 *
 * Copyright (C) 2020-2021 by the authors.
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

#include <algorithm>
#include <array>
#include <vector>

/**
 * Return true if _container_ contains _value_.
 */
template<typename Container, typename Value>
bool container_contains(const Container& container, const Value& value)
{
    return container.find(value) != container.end();
}

/**
 * Return true if _string_ contains _value_ (character or substring).
 */
template<typename String, typename Value>
bool string_contains(const String& string, const Value& value)
{
    return string.find(value) != String::npos;
}

// Replace by ranges::shuffle in C++20.
/**
 * Reorder items in container.
 */
template<typename Container, typename Generator>
void shuffle_container(Container& container, Generator& generator)
{
    std::shuffle(container.begin(), container.end(), generator);
}

/**
 * Select a random item from a container.
 *
 * May be slow if it takes a long time to increase the iterator (e.g., for std::set) and
 * the container is large.
 */
template<typename Container, typename Generator>
typename Container::value_type
pick_random_item(const Container& container, Generator& generator)
{
    // Replace .size() call by std::size in C++17.
    std::uniform_int_distribution<size_t> distribution(0, container.size() - 1);
    auto index = distribution(generator);
    // For small containers, this is expected to be fast for both sets and vectors.
    return *std::next(container.begin(), index);
}

/**
 * \brief A const iterator which encapsulates either forward or reverse iterator.
 *
 * Uses the iterator was passed in the constructor, but it can contain either
 * forward or reverse iterator, so it allows writing direction agnostic code
 * where the direction is determined in the runtime.
 */
template<typename Container>
class ConstEitherWayIterator
{
public:
    ConstEitherWayIterator(typename Container::const_iterator it)
        : is_forward_(true), forward_it_(it)
    {}
    ConstEitherWayIterator(typename Container::const_reverse_iterator it)
        : is_forward_(false), reverse_it_(it)
    {}
    ConstEitherWayIterator& operator++()
    {
        if (is_forward_)
            ++forward_it_;
        else
            ++reverse_it_;
        return *this;
    }
    const typename Container::value_type& operator*() const
    {
        if (is_forward_)
            return *forward_it_;
        return *reverse_it_;
    }
    bool operator!=(const ConstEitherWayIterator& other)
    {
        if (is_forward_)
            return forward_it_ != other.forward_it_;
        return reverse_it_ != other.reverse_it_;
    }

private:
    // Can be improved by std::optional or any in C++17.
    bool is_forward_;
    typename Container::const_iterator forward_it_;
    typename Container::const_reverse_iterator reverse_it_;
};

/**
 * \brief A const view of a container which can be iterated.
 *
 * Depending on which iterators are passed in the constructor, the view
 * is in the forward direction of the container or in reverse. The view can also be
 * a subset. The view supports only const iterators.
 */
template<typename Container>
class ContainerView
{
public:
    ContainerView(
        typename Container::const_iterator first,
        typename Container::const_iterator last)
        : begin_(first), end_(last)
    {}
    ContainerView(
        typename Container::const_reverse_iterator first,
        typename Container::const_reverse_iterator last)
        : begin_(first), end_(last)
    {}
    ConstEitherWayIterator<Container> begin() const
    {
        return begin_;
    }
    ConstEitherWayIterator<Container> end() const
    {
        return end_;
    }

private:
    ConstEitherWayIterator<Container> begin_;
    ConstEitherWayIterator<Container> end_;
};

typedef std::tuple<int, int, int, int> BBoxInt;
typedef std::tuple<double, double, double, double> BBoxFloat;
typedef std::tuple<bool, bool, bool, bool> BBoxBool;

// C++20 NumericType
template<typename Number>
struct BBox
{
    Number north;
    Number south;
    Number east;
    Number west;

    BBox() : north(0), south(0), east(0), west(0) {}
};

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

/**
 * Create a list of suitable cells in from a host raster.
 *
 * Suitable cell is defined as cell with value higher than zero, i.e., there is at least
 * one host.
 *
 * Suitable cells datastructure is vector of vectors where the nested vector always has
 * size equal to two and contains row and column index. The type was chosen to work well
 * with Rcpp.
 *
 * First template parameter is the index type for the resulting sutibale cell indices.
 * Second template parameter is deduced automatically from the function parameter.
 */
template<typename RasterIndex, typename RasterType>
std::vector<std::vector<RasterIndex>> find_suitable_cells(const RasterType& raster)
{
    std::vector<std::vector<RasterIndex>> cells;
    // The assumption is that the raster is sparse (otherwise we would not be doing
    // this), so we have no number for the reserve method.
    for (RasterIndex row = 0; row < raster.rows(); ++row) {
        for (RasterIndex col = 0; col < raster.cols(); ++col) {
            if (raster(row, col) > 0) {
                cells.push_back({row, col});
            }
        }
    }
    return cells;
}

#endif  // POPS_UTILS_HPP
