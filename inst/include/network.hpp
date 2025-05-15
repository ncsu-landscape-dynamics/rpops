/*
 * PoPS model - network dispersal kernel
 *
 * Copyright (C) 2020-2022 by the authors.
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

#ifndef POPS_NETWORK_HPP
#define POPS_NETWORK_HPP

#include "raster.hpp"
#include "utils.hpp"

#include <algorithm>
#include <set>
#include <random>
#include <string>
#include <sstream>
#include <map>
#include <cmath>
#include <vector>

namespace pops {

/** Edge geometry, i.e., cells connecting two nodes (aka segment) */
template<typename RasterCell>
class EdgeGeometry : public std::vector<RasterCell>
{
public:
    using typename std::vector<RasterCell>::const_reference;
    using typename std::vector<RasterCell>::size_type;

    /** Get index from cost */
    size_type index_from_cost(double cost) const
    {
        // This is like
        // index = (max_index / total_cost) * cost
        // but cost_per_cell == total_cost / max_index, so
        // index = 1 / (total_cost / max_index) * cost.
        return std::lround(cost / this->cost_per_cell());
    }

    /** Get cell by cost instead of an index */
    const_reference cell_by_cost(double cost) const
    {
        auto index = this->index_from_cost(cost);
        return this->operator[](index);
    }

    /** Get cost of the whole segment (edge). */
    double cost() const
    {
        if (static_cast<bool>(total_cost_))
            return total_cost_;
        // This is short for ((size - 2) + (2 * 1/2)) * cost per cell.
        return (this->size() - 1) * cost_per_cell_;
    }

    void set_total_cost(double value)
    {
        total_cost_ = value;
    }

    /** Get cost per cell for the segment (edge). */
    double cost_per_cell() const
    {
        if (static_cast<bool>(total_cost_))
            return total_cost_ / (this->size() - 1);
        return cost_per_cell_;
    }

    void set_cost_per_cell(double value)
    {
        cost_per_cell_ = value;
    }

    /** Get probability of an edge */
    double probability() const
    {
        return probability_;
    }

    /** Set probability of an edge */
    void set_probability(double value)
    {
        probability_ = value;
    }

private:
    double cost_per_cell_ = 0;
    double total_cost_ = 0;
    double probability_ = 0;
};

/** Constant view of a edge geometry (to iterate a segment in either direction)
 *
 * Notably, the view uses iterators to flip the direction of the geometry,
 * but the total cost is still for the whole geometry, i.e., this is not view of
 * a part of the geometry, but view of a potentially reversed geometry.
 */
template<typename EdgeGeometryType>
class EdgeGeometryView : public ContainerView<EdgeGeometryType>
{
public:
    EdgeGeometryView(
        typename EdgeGeometryType::const_iterator first,
        typename EdgeGeometryType::const_iterator last,
        const EdgeGeometryType& segment)
        : ContainerView<EdgeGeometryType>(first, last), segment_(segment)
    {}
    EdgeGeometryView(
        typename EdgeGeometryType::const_reverse_iterator first,
        typename EdgeGeometryType::const_reverse_iterator last,
        const EdgeGeometryType& segment)
        : ContainerView<EdgeGeometryType>(first, last), segment_(segment)
    {}

    /** Get cell by cost instead of an index */
    typename EdgeGeometryType::const_reference cell_by_cost(double cost) const
    {
        // The index is without a direction, so we can use it in reverse too.
        auto index = segment_.index_from_cost(cost);
        return this->operator[](index);
    }

    /** Get cost of the whole underlying segment (edge).
     *
     * Notably, this function assumes that the view represents the whole segment,
     * not just part of it.
     */
    double cost() const
    {
        return segment_.cost();
    }

    /** Get cost per cell for the underlying segment (edge). */
    double cost_per_cell() const
    {
        return segment_.cost_per_cell();
    }

private:
    const EdgeGeometryType& segment_;
};

/**
 * Network structure and algorithms
 *
 * The network consists of nodes connected by edges. Edges are spatially represented
 * as segments. Edges themselves don't carry cost and have meaning only as indicators of
 * existence of a connection between nodes. Cost for each edge is a travel distance
 * (cost) determined by advancing through all cells in a segment.
 *
 * Nodes are the hop-on locations for dispersers. The dispersers can hop-off anywhere.
 *
 * The general workflow is constructing the object (with the constructor) and loading
 * the data (with the load() function). Then the network is ready to be used for
 * simulating trips over the network (with the walk() or teleport() functions).
 *
 * When the walk() or teleport() functions are used from a kernel, user of the network
 * directly calls only the setup functions.
 *
 * The class exposes number of functions as public which are meant for testing or other
 * special workflows.
 *
 * RasterIndex template parameter is an integral type used for indexing in the grid,
 * specifically used as a type for rows and columns.
 *
 * Only non-const (write) functions for network are those loading the network, so an
 * existing network can be shared by multiple kernel or model instances.
 */
template<typename RasterIndex>
class Network
{
public:
    using NodeId = int;  ///< Type for node IDs
    using Statistics = std::map<std::string, int>;  ///< Type for summary statistics

    /**
     * @brief Construct empty network.
     *
     * @param bbox Bounding box of the raster grid (in real world coordinates)
     * @param ew_res East-west resolution of the raster grid
     * @param ns_res North-south resolution of the raster grid
     */
    Network(BBox<double> bbox, double ew_res, double ns_res)
        : bbox_(bbox),
          ew_res_(ew_res),
          ns_res_(ns_res),
          max_row_(0),
          max_col_(0),
          distance_per_cell_((ew_res + ns_res) / 2)
    {
        std::tie(max_row_, max_col_) = xy_to_row_col(bbox_.east, bbox_.south);
    }

    /**
     * @brief Create an empty network not meant for futher use.
     *
     * This is useful when a network object is needed to construct a kernel, but it will
     * not be used in runtime.
     *
     * @return New Network object
     */
    static Network null_network()
    {
        return Network(BBox<double>(), 0, 0);
    }

    /**
     * @brief Convert real world coordinates to row and column in the raster grid.
     *
     * @param x The X coordinate (east-west or column direction)
     * @param y The Y coordinate (north-south or row direction)
     * @return Row and column
     */
    std::pair<RasterIndex, RasterIndex> xy_to_row_col(double x, double y) const
    {
        double col = (x - bbox_.west) / ew_res_;
        double row = (bbox_.north - y) / ns_res_;
        return {std::floor(row), std::floor(col)};
    }

    /**
     * Coverts coordinates as xy_to_row_col(double, double) but converts from strings.
     */
    std::pair<RasterIndex, RasterIndex>
    xy_to_row_col(const std::string& x, const std::string& y) const
    {
        try {
            return xy_to_row_col(std::stod(x), std::stod(y));
        }
        catch (const std::invalid_argument&) {
            if (string_contains(x, '"') || string_contains(x, '\'')) {
                throw std::invalid_argument(
                    std::string("Text for coordinate X cannot contain quotes: ") + x);
            }
            if (string_contains(y, '"') || string_contains(y, '\'')) {
                throw std::invalid_argument(
                    std::string("Text for coordinate Y cannot contain quotes: ") + y);
            }
            if (x.empty()) {
                throw std::invalid_argument(
                    "Text for coordinate X cannot be an empty string, Y is " + y);
            }
            if (y.empty()) {
                throw std::invalid_argument(
                    "Text for coordinate Y cannot be an empty string, X is " + y);
            }
            throw std::invalid_argument(
                std::string("Text cannot be converted to X and Y coordinates: X: ") + x
                + " Y: " + y);
        }
        catch (const std::out_of_range&) {
            throw std::out_of_range(
                std::string("Numerical value too large for X or Y coordinate: X: ") + x
                + " Y: " + y);
        }
    }

    /**
     * @brief Test if coordinates XY are in the bounding box.
     * @param x Real world coordinate X
     * @param y Real world coordinate Y
     * @return True if XY is in the bounding box, false otherwise.
     */
    bool xy_out_of_bbox(double x, double y) const
    {
        return x > bbox_.east || x < bbox_.west || y > bbox_.north || y < bbox_.south;
    }

    /**
     * @brief Test if a cell is in the bounding box.
     * @param row Row in the raster grid
     * @param col Column in the raster grid
     * @return True if cells is in the bounding box, false otherwise.
     */
    bool row_col_out_of_bbox(RasterIndex row, RasterIndex col) const
    {
        return row > max_row_ || row < 0 || col > max_col_ || col < 0;
    }

    /**
     * @brief Test if a cell is in the bounding box.
     * @param cell Pair consisting of a row and column indices
     * @return True if cell is in the bounding box, false otherwise.
     */
    bool cell_out_of_bbox(const std::pair<RasterIndex, RasterIndex>& cell) const
    {
        return cell.first > max_row_ || cell.first < 0 || cell.second > max_col_
               || cell.second < 0;
    }

    /**
     * @brief Load network from an input stream.
     *
     * @param stream Input stream containing text records for network
     * @param allow_empty True if the loaded network can be empty
     *
     * Network stream is text in a custom format resembling CSV geared towards
     * parsing, not readability, with rows `node_id_1,node_id_2,X1;Y1;X2;Y2;X3;Y3;...`.
     * Each line represents an edge (pair of connected nodes) and the spatial
     * representation of the edge (here called a segment). Node coordinates are taken
     * from the first and last coordinate pair in the segment.
     *
     * Optionally, a header can be provided in the CSV which can specify additional
     * input format variations. The default format described above has either no header
     * or `node_1,node_2,geometry`. To specify cost for each edge, an additional cost
     * column is needed and needs to be placed after the node columns and before the
     * geometry column, i.e., the header will look like `node_1,node_2,cost,geometry`.
     *
     * The function can take any input stream which behaves like std::ifstream or
     * std::istringstream. These are also the two expected streams to be used
     * (for reading from a file and reading from a string, respectively).
     *
     * A simple code for reading from files without any checking may look like this:
     *
     * ```
     * Network network(...);
     * std::string filename = "network.txt"
     * std::ifstream stream{filename};
     * network.load(stream);
     * ```
     *
     * The input network is a list of segments where each segment consists of a pairs of
     * nodes and list of coordinates. In the terminology of the this class, the pair of
     * nodes is called an edge and the edge with coordinates (or rows and columns) is
     * called a segment.
     *
     * The input is used as-is. The network is agnostic
     * towards how the segments look like in terms of crossing each other or other
     * nodes.
     *
     * All XY coordinates of segment points including nodes coordinates are converted
     * to row and column using bounding box and resolution.
     *
     * Coordinates for nodes are taken from the beginning and the end of each segment,
     * so input segment coordinates should include the node coordinates.
     *
     * Only edges (segments) with both nodes in the bounding box
     * are considered, i.e., input network is reduced to the bounding box during
     * reading, but clipping happens on the level of whole edges, not in the middle of
     * a segment.
     *
     * Standalone nodes (without any connection) are theoretically allowed in the
     * internal representation of the network and handled
     * as no movement from the source cell, but the input always needs to contain an
     * edge.
     */
    template<typename InputStream>
    void load(InputStream& stream, bool allow_empty = false)
    {
        load_segments(stream);
        if (node_matrix_.empty()) {
            if (allow_empty)
                return;
            else
                throw std::runtime_error("Network: No nodes within the extend");
        }
    }

    /**
     * @brief Get a list of nodes at a given cell
     *
     * Returns a const reference to an internally stored set. If there are no nodes
     * at a given cell, a reference to an empty set is returned.
     *
     * @param row Row in the raster grid
     * @param col Column in the raster grid
     * @return List of node IDs as a set
     */
    const std::set<NodeId>& get_nodes_at(RasterIndex row, RasterIndex col) const
    {
        auto it = nodes_by_row_col_.find(std::make_pair(row, col));
        if (it != nodes_by_row_col_.end())
            return it->second;
        // If not found, return an empty set (const reference to local static).
        static std::set<NodeId> empty;
        return empty;
    }

    /**
     * @brief Get a randomly selected node at a given cell
     *
     * There can be more than one node at a given cell. This function selects one
     * of these nodes randomly.
     *
     * @param row Row in the raster grid
     * @param col Column in the raster grid
     * @param generator Random number generator
     * @return Node ID
     * @throws std::invalid_argument if there is no node at a given cell
     */
    template<typename Generator>
    NodeId
    get_random_node_at(RasterIndex row, RasterIndex col, Generator& generator) const
    {
        const auto& nodes = get_nodes_at(row, col);
        auto num_nodes = nodes.size();
        if (num_nodes == 1) {
            return *nodes.begin();
        }
        else if (num_nodes > 1) {
            return pick_random_item(nodes, generator);
        }
        else {
            throw std::invalid_argument(
                "No nodes at a given row and column: " + std::to_string(row) + ", "
                + std::to_string(col));
        }
    }

    /**
     * @brief Test if the network has a node at a given row and column
     * @param row Row
     * @param col Column
     * @return True if there is at least one node, false otherwise
     */
    bool has_node_at(RasterIndex row, RasterIndex col) const
    {
        // Replace by contains in C++20.
        // return nodes_by_row_col_.count(std::make_pair(row, col)) > 0;
        auto it = nodes_by_row_col_.find(std::make_pair(row, col));
        if (it == nodes_by_row_col_.end())
            return false;
        return true;
    }

    /**
     * @brief Get row and column for a node.
     * @param node Node id to get the coordinates for
     * @return Row and column pair
     */
    std::pair<RasterIndex, RasterIndex> get_node_row_col(NodeId node) const
    {
        for (const auto& item : nodes_by_row_col_) {
            if (container_contains(item.second, node))
                return item.first;
        }
        throw std::invalid_argument("No node with a given id");
    }

    /**
     * Walk a given distance (cost) in the network from given row and column.
     *
     * All previously visited nodes are tracked and, if possible, excluded
     * from further traveling.
     *
     * The function assumes there is a node at the given *row* and *column*, i.e., that
     * the decision to call this function was based on the caller knowing there is a
     * node. If there is no node, an std::invalid_argument exception is thrown.
     * If there is more than one node at the given *row* and *column*, a random node is
     * picked and used as a next walking destination.
     *
     * If *jump* is true, then results are snapped to the closest node, otherwise
     * result can be anywhere in between the nodes based on the edge geomerty (segment).
     *
     * @returns Final row and column pair
     */
    template<typename Generator>
    std::tuple<int, int> walk(
        RasterIndex row,
        RasterIndex col,
        double distance,
        Generator& generator,
        bool jump = false) const
    {
        auto node_id = get_random_node_at(row, col, generator);
        std::set<NodeId> visited_nodes;
        while (distance >= 0) {
            auto next_node_id = next_node(node_id, visited_nodes, generator);
            // We have visited the current node (initial start node or end node from
            // last iteration. (There is no need to tell next_node that the current node
            // is visited, but we need to tell it the next time because it won't be
            // current anymore.)
            visited_nodes.insert(node_id);
            // If there is no segment from the node, return the start cell.
            if (next_node_id == node_id)
                return std::make_tuple(row, col);
            auto segment = get_segment(node_id, next_node_id);
            // Set node ID for the next iteration.
            node_id = next_node_id;

            if (distance > segment.cost()) {
                // Go over the whole segment.
                distance -= segment.cost();
                continue;
            }
            if (jump) {
                if (distance < segment.cost() / 2) {
                    // Less than half snaps to the start node.
                    return segment.front();
                }
                // Half or more snaps to the end node.
                return segment.back();
            }
            // No jumping (snapping), advance in a segment.
            // This includes the special cases when distance is 0 or total segment cost.
            return segment.cell_by_cost(distance);
        }
        throw std::invalid_argument("Distance must be greater than or equal to zero");
    }

    /**
     * Teleport to a different node in the network from given row and column.
     *
     * Returns any node of the nodes connected to the start node possibly based on the
     * edge probability if probability was assigned to the edges without considering
     * cost to travel from one node to the next one.
     *
     * If *num_steps* is greater than 1, multiple steps are performed and the last node
     * is returned. In each node, the probability of picking a specific connection is
     * either determined by the provided edge probabilities or is equal among the
     * connections. Consequently, previously visited nodes can be visited again. In
     * other words, for highly probable connections, most likely next step is a step
     * back to the starting node (or generally previous node for `num_steps >= 3`).
     *
     * The function assumes a node is at the *row*, *col* coordinates, i.e., that this
     * was either checked beforehand or otherwise ensured. If there is no node, an
     * std::invalid_argument exception is thrown.
     * If there is more than one node at the given *row* and *column*, a random node is
     * picked and used.
     *
     * @returns Destination row and column pair
     */
    template<typename Generator>
    std::tuple<int, int> teleport(
        RasterIndex row, RasterIndex col, Generator& generator, int num_steps = 1) const
    {
        auto node_id = get_random_node_at(row, col, generator);
        for (int i = 0; i < num_steps; ++i) {
            node_id = next_probable_node(node_id, generator);
        }
        return get_node_row_col(node_id);
    }

    /**
     * @brief Get all nodes as vector of all ids with their row and column.
     *
     * If there is more than one node at a given cell (row and column),
     * each of these nodes is returned as a separate item in the list.
     *
     * This translates the internal representation and returns a new object.
     *
     * @return A vector of pairs of id and row and col pair.
     */
    std::vector<std::pair<NodeId, std::pair<RasterIndex, RasterIndex>>>
    get_all_nodes() const
    {
        std::vector<std::pair<NodeId, std::pair<RasterIndex, RasterIndex>>> nodes;
        nodes.reserve(nodes_by_row_col_.size());  // It will be at least this big.
        for (const auto& item : nodes_by_row_col_) {
            RasterIndex row = item.first.first;
            RasterIndex col = item.first.second;
            for (const auto node_id : item.second) {
                nodes.emplace_back(node_id, std::make_pair(row, col));
            }
        }
        return nodes;
    }

    /**
     * @brief Collect statistics about the network
     *
     * Computes statistics such as number of nodes and segments.
     *
     * @return Associative array with statistics
     */
    Statistics collect_statistics() const
    {
        std::map<std::string, int> stats;
        std::set<NodeId> node_ids;
        for (const auto& item : nodes_by_row_col_) {
            for (const auto node_id : item.second) {
                node_ids.insert(node_id);
            }
        }
        stats["num_nodes"] = node_ids.size();
        // We store segments in both directions, so each segment is stored twice.
        stats["num_segments"] = node_matrix_.size() / 2;
        std::set<NodeId> nodes_with_segments;
        for (const auto& item : node_matrix_) {
            nodes_with_segments.insert(item.first);
        }
        stats["num_nodes_with_segments"] = nodes_with_segments.size();
        // TODO: Only count the nodes here. Output them in the dump network output.
        int num_standalone_nodes = 0;
        for (NodeId node_id : node_ids) {
            if (nodes_with_segments.find(node_id) == nodes_with_segments.end()) {
                stats["standalone_node_" + std::to_string(++num_standalone_nodes)] =
                    node_id;
            }
        }
        stats["num_standalone_nodes"] = num_standalone_nodes;
        // TODO: This can be faster. The list is a ordered set and we iterate over
        // elements above.
        const auto minmax_node_id = minmax_element(node_ids.begin(), node_ids.end());
        stats["min_node_id"] = *minmax_node_id.first;
        stats["max_node_id"] = *minmax_node_id.second;
        return stats;
    }

    /**
     * @brief Output network to a stream.
     *
     * Writes statistics and different views of the network as a YAML file.
     */
    template<typename OutputStream>
    void dump_yaml(OutputStream& stream) const
    {
        stream << "network:\n";
        stream << "  statistics:\n";
        auto stats = collect_statistics();
        for (const auto& item : stats) {
            stream << "    " << item.first << ": " << item.second << "\n";
        }
        stream << "  extent:\n";
        stream << "    north: " << bbox_.north << "\n";
        stream << "    south: " << bbox_.south << "\n";
        stream << "    east: " << bbox_.east << "\n";
        stream << "    west: " << bbox_.west << "\n";
        stream << "  resolution:\n";
        stream << "    ns: " << ns_res_ << "\n";
        stream << "    ew: " << ew_res_ << "\n";
        stream << "  raster:\n";
        RasterIndex min_row;
        RasterIndex min_col;
        RasterIndex max_row;
        RasterIndex max_col;
        std::tie(min_row, min_col) =
            xy_to_row_col(bbox_.west + ew_res_ / 2, bbox_.north - ns_res_ / 2);
        std::tie(max_row, max_col) =
            xy_to_row_col(bbox_.east - ew_res_ / 2, bbox_.south + ew_res_ / 2);
        stream << "    min_row: " << min_row << "\n";
        stream << "    min_col: " << min_col << "\n";
        stream << "    max_row: " << max_row << "\n";
        stream << "    max_col: " << max_col << "\n";
        stream << "    num_rows: " << max_row - min_row + 1 << "\n";
        stream << "    num_cols: " << max_col - min_col + 1 << "\n";
        stream << "  cost:\n";
        stream << "    distance_per_cell: " << distance_per_cell_ << "\n";
        std::set<std::pair<NodeId, NodeId>> edges;
        for (const auto& item : node_matrix_) {
            for (const auto& node_id : item.second.second) {
                edges.emplace(item.first, node_id);
            }
        }
        stream << "  edges:\n";
        for (const auto& item : edges) {
            stream << "    - [" << item.first << ", " << item.second << "]\n";
        }
        stream << "  nodes:\n";
        for (const auto& item : nodes_by_row_col_) {
            auto row = item.first.first;
            auto col = item.first.second;
            for (const auto& node_id : item.second) {
                stream << "    - id: " << node_id << "\n";
                stream << "      row: " << row << "\n";
                stream << "      col: " << col << "\n";
            }
        }
        stream << "  segments:\n";
        for (const auto& item : segments_by_nodes_) {
            stream << "    - start_node: " << item.first.first << "\n";
            stream << "      end_node: " << item.first.second << "\n";
            stream << "      cells: [";
            for (const auto& cell : item.second) {
                stream << "[" << cell.first << ", " << cell.second << "], ";
            }
            stream << "]\n";
        }
    }

protected:
    /** Node connections (edges)
     *
     * The matrix is sparse. One node can be connected to multiple nodes.
     * Connection in both directions is stored explicitly (see load_segments() source
     * code).
     */
    using NodeMatrix =
        std::map<NodeId, std::pair<std::vector<double>, std::vector<NodeId>>>;

    /** Smallest component of edge geometry */
    using Cell = std::pair<RasterIndex, RasterIndex>;

    /** Cells connecting two nodes (segment between nodes) */
    using Segment = EdgeGeometry<Cell>;

    /** Constant view of a segment (to iterate a segment in either direction) */
    using SegmentView = EdgeGeometryView<Segment>;

    /** Segments by nodes (edges) */
    using SegmentsByNodes = std::map<std::pair<NodeId, NodeId>, Segment>;

    /**
     * @brief Convert string to node ID
     *
     * @param text String with node ID
     * @return Node ID
     */
    static NodeId node_id_from_text(const std::string& text)
    {
        try {
            return std::stoi(text);
        }
        catch (const std::invalid_argument& err) {
            if (string_contains(text, '"') || string_contains(text, '\'')) {
                throw std::invalid_argument(
                    std::string("Text for a node ID cannot contain quotes "
                                "(only digits are allowed): ")
                    + text);
            }
            if (text.empty()) {
                throw std::invalid_argument(
                    "Text for a node ID cannot be an empty string");
            }
            else {
                throw std::invalid_argument(
                    std::string("Text cannot be converted to a node ID "
                                "(only digits are allowed): ")
                    + text);
            }
        }
        catch (const std::out_of_range& err) {
            throw std::out_of_range(
                std::string("Numerical value too large for a node ID: ") + text);
        }
    }

    /**
     * @brief Convert string to probability
     *
     * @param text String with probability
     * @return Probability as number
     */
    static double probability_from_text(const std::string& text)
    {
        double value{0};
        try {
            value = std::stod(text);
        }
        catch (const std::invalid_argument& err) {
            if (string_contains(text, '"') || string_contains(text, '\'')) {
                throw std::invalid_argument(
                    std::string("Text for connection probabilty cannot contain quotes "
                                "(only digits are allowed): ")
                    + text);
            }
            if (text.empty()) {
                throw std::invalid_argument(
                    "Text for connection probabilty cannot be an empty string");
            }
            else {
                throw std::invalid_argument(
                    std::string("Text cannot be converted to connection probabilty "
                                "(only digits are allowed): ")
                    + text);
            }
        }
        catch (const std::out_of_range& err) {
            throw std::out_of_range(
                std::string("Numerical value too large for connection probabilty: ")
                + text);
        }
        if (value < 0) {
            throw std::invalid_argument(
                std::string("Probability needs to be >=0, not ") + text);
        }
        return value;
    }

    /**
     * @brief Convert string to cost
     *
     * @param text String with cost
     * @return Cost as number
     */
    static double cost_from_text(const std::string& text)
    {
        try {
            return std::stod(text);
        }
        catch (const std::invalid_argument& err) {
            if (string_contains(text, '"') || string_contains(text, '\'')) {
                throw std::invalid_argument(
                    std::string("Text for cost cannot contain quotes "
                                "(only digits are allowed): ")
                    + text);
            }
            if (text.empty()) {
                throw std::invalid_argument("Text for cost cannot be an empty string");
            }
            else {
                throw std::invalid_argument(
                    std::string("Text cannot be converted to cost "
                                "(only digits are allowed): ")
                    + text);
            }
        }
        catch (const std::out_of_range& err) {
            throw std::out_of_range(
                std::string("Numerical value too large for cost: ") + text);
        }
    }

    template<typename InputStream>
    static std::pair<bool, bool> stream_has_columns(InputStream& stream, char delimeter)
    {
        bool has_cost{false};
        bool has_probability{false};
        // Get header to determine what is included.
        auto starting_position = stream.tellg();
        std::string line;
        std::getline(stream, line);
        std::istringstream line_stream{line};
        std::string label;
        int column_number = 0;
        while (std::getline(line_stream, label, delimeter)) {
            column_number++;
            if (column_number == 1 && label != "node_1") {
                // The right label is not there. Assuming it is without a header.
                stream.seekg(starting_position);
                break;
            }
            if (label == "probability") {
                if (has_cost) {
                    // Detailed check to give a more relevant error message.
                    throw std::runtime_error(
                        "The cost column must be after the probability column");
                }
                if (column_number != 3) {
                    throw std::runtime_error(
                        "The probability column must be the third column");
                }
                has_probability = true;
                continue;
            }
            if (label == "cost") {
                if (!(column_number == 3 || column_number == 4)) {
                    throw std::runtime_error(
                        "The cost column must be the third or fourth column");
                }
                has_cost = true;
                continue;
            }
        }
        return {has_cost, has_probability};
    }

    /**
     * @brief Read segments from a stream.
     *
     * @param stream Input stream with segments
     */
    template<typename InputStream>
    void load_segments(InputStream& stream)
    {
        char delimeter{','};
        bool has_cost{false};
        bool has_probability{false};
        std::tie(has_cost, has_probability) = stream_has_columns(stream, delimeter);

        std::string line;
        while (std::getline(stream, line)) {
            std::istringstream line_stream{line};
            std::string node_1_text;
            std::getline(line_stream, node_1_text, delimeter);
            std::string node_2_text;
            std::getline(line_stream, node_2_text, delimeter);
            auto node_1_id = node_id_from_text(node_1_text);
            auto node_2_id = node_id_from_text(node_2_text);
            if (node_1_id < 1 || node_2_id < 1) {
                throw std::runtime_error(
                    std::string("Node ID must be greater than zero (node 1, node 2): ")
                    + node_1_text + ", " + node_2_text + ", line: " + line);
            }
            Segment segment;

            if (has_probability) {
                std::string probability_text;
                std::getline(line_stream, probability_text, delimeter);
                double connection_probability = probability_from_text(probability_text);
                segment.set_probability(connection_probability);
            }
            if (has_cost) {
                std::string cost_text;
                std::getline(line_stream, cost_text, delimeter);
                double cost = cost_from_text(cost_text);
                segment.set_total_cost(cost);
            }
            else {
                segment.set_cost_per_cell(distance_per_cell_);
            }

            std::string segment_text;
            std::getline(line_stream, segment_text, delimeter);
            std::istringstream segment_stream{segment_text};
            char in_cell_delimeter{';'};
            std::string x_coord_text;
            std::string y_coord_text;

            long int loaded_coord_pairs = 0;
            while (std::getline(segment_stream, x_coord_text, in_cell_delimeter)
                   && std::getline(segment_stream, y_coord_text, in_cell_delimeter)) {
                // The same cell is possibly repeated if raster resolution is lower than
                // the detail ("resolution") of each segment, so we check if the last
                // added coordinate pair converted to same cell.
                auto new_point = xy_to_row_col(x_coord_text, y_coord_text);
                if (segment.empty() || segment.back() != new_point)
                    segment.emplace_back(new_point);
                ++loaded_coord_pairs;
            }

            if (segment.empty()) {
                throw std::runtime_error(
                    std::string("Row for an edge between nodes ") + node_1_text
                    + " and " + node_2_text + " does not have any node coordinates");
            }
            if (loaded_coord_pairs < 2) {
                throw std::runtime_error(
                    std::string("Row for an edge between nodes ") + node_1_text
                    + " and " + node_2_text + " has only 1 coordinate pair "
                    + "(at least two are needed, one for each node), "
                    + "the one coordinate pair was: " + x_coord_text + ", "
                    + y_coord_text);
            }
            // If the two nodes has the same coordinates (e.g., after computing the cell
            // from the real-world coordinates, the segment geometry list contains only
            // one pair, so size equal to one is considered correct.
            // However, we need put back the missing coordinates for the second node.
            if (segment.size() == 1)
                segment.push_back(typename Segment::value_type(segment.back()));

            // If either node of the segment is not in the extent, skip the segment.
            // This means that a segment is ignored even if one of the nodes and
            // significant portion of the segment is in the area of iterest.
            // TODO: Part of the segment may still be out, so that needs to be checked.
            // Cut to extend
            if (cell_out_of_bbox(segment.front()) || cell_out_of_bbox(segment.back()))
                continue;
            // We are done with the segment data, so we move them to the attribute.
            segments_by_nodes_.emplace(
                std::make_pair(node_1_id, node_2_id), std::move(segment));
        }

        for (const auto& node_segment : segments_by_nodes_) {
            const auto& start_node_id{node_segment.first.first};
            const auto& end_node_id{node_segment.first.second};
            const auto& segment{node_segment.second};
            nodes_by_row_col_[segment.front()].insert(start_node_id);
            nodes_by_row_col_[segment.back()].insert(end_node_id);
            if (has_probability) {
                node_matrix_[start_node_id].first.push_back(segment.probability());
                node_matrix_[end_node_id].first.push_back(segment.probability());
            }
            node_matrix_[start_node_id].second.push_back(end_node_id);
            node_matrix_[end_node_id].second.push_back(start_node_id);
        }
    }

    /**
     * @brief Get a segment between the two given nodes
     *
     * For a given pair of nodes, the function finds the connecting segment
     * and returns a view of the segment oriented according to the order of the nodes.
     *
     * @param start Start node id
     * @param end End node id
     * @returns View of the segment oriented from start to end
     * @throws std::invalid_argument if there is no segment connecting the nodes
     */
    SegmentView get_segment(NodeId start, NodeId end) const
    {
        auto it = segments_by_nodes_.find(std::make_pair(start, end));
        if (it != segments_by_nodes_.end()) {
            return SegmentView(it->second.cbegin(), it->second.cend(), it->second);
        }
        it = segments_by_nodes_.find(std::make_pair(end, start));
        if (it != segments_by_nodes_.end()) {
            return SegmentView(it->second.crbegin(), it->second.crend(), it->second);
        }
        throw std::invalid_argument(std::string(
            "No segment for given nodes: " + std::to_string(start) + " "
            + std::to_string(end)));
    }

    /**
     * @brief Get nodes connected by an edge to a given node.
     *
     * @param node Node to get connections from
     * @return List of connected nodes
     */
    const std::vector<NodeId>& nodes_connected_to(NodeId node) const
    {
        return node_matrix_.at(node).second;
    }

    /**
     * @brief Pick a next node from the given node.
     *
     * If there is more than one edge leading from the given node, a random node is
     * picked. If there are no edges leading from the given node, the node is returned.
     *
     * The random node is picked from candidate nodes which are nodes connected to
     * a given node. The candidate nodes which are in the *ignore* list are excluded
     * from the random selection. If all candidate nodes are in the *ignore* list,
     * the *ignore* list is ignored and all candidate nodes are used.
     *
     * The function always returns a node to go to even if it means going to an ignored
     * node or returning the value of the *node* parameter.
     */
    template<typename Generator>
    NodeId
    next_node(NodeId node, const std::set<NodeId>& ignore, Generator& generator) const
    {
        // Get all candidate nodes.
        const auto& all_nodes = nodes_connected_to(node);

        // Resolve disconnected node and dead end cases.
        auto num_nodes = all_nodes.size();
        if (!num_nodes)
            return node;
        else if (num_nodes == 1)
            return all_nodes[0];

        // Filter out the ignored nodes.
        std::vector<int> nodes;
        std::back_insert_iterator<std::vector<int>> back_it(nodes);
        std::remove_copy_if(
            all_nodes.begin(), all_nodes.end(), back_it, [&ignore](NodeId id) {
                return container_contains(ignore, id);
            });

        // Pick a random node. Fallback to all candidate nodes if all are in ignore.
        num_nodes = nodes.size();
        if (!num_nodes)
            return pick_random_item(all_nodes, generator);
        else if (num_nodes == 1)
            return nodes[0];
        return pick_random_item(nodes, generator);
    }

    /**
     * @brief Pick a probable node from the given node.
     *
     * If there is more than one edge leading from the given node, a random node is
     * picked based on edge probabilities. If there are no probabilities assigned, a
     * random node is picked. If there are no edges leading from the given node, the
     * node itself is returned.
     */
    template<typename Generator>
    NodeId next_probable_node(NodeId node, Generator& generator) const
    {
        // Get all candidate nodes.
        const auto& record{node_matrix_.at(node)};
        const auto& probabilities{record.first};
        const auto& nodes{record.second};

        // Resolve disconnected node and dead end cases.
        auto num_nodes = nodes.size();
        if (!num_nodes)
            return node;
        else if (num_nodes == 1)
            return nodes[0];

        // Pick nodes based on edge probabilities if they are available.
        if (!probabilities.empty()) {
            std::discrete_distribution<int> dd{
                probabilities.begin(), probabilities.end()};
            return nodes.at(dd(generator));
        }
        // Pick a connected node with equal edge probabilities.
        return pick_random_item(nodes, generator);
    }

    BBox<double> bbox_;  ///< Bounding box of the network grid in real world coordinates
    double ew_res_;  ///< East-west resolution of the grid
    double ns_res_;  ///< North-south resolution of the grid
    RasterIndex max_row_;  ///< Maximum row index in the grid
    RasterIndex max_col_;  ///< Maximum column index in the grid
    double distance_per_cell_;  ///< Distance (cost) to walk through one cell
    /** Node IDs stored by row and column (multiple nodes per cell) */
    std::map<std::pair<RasterIndex, RasterIndex>, std::set<NodeId>> nodes_by_row_col_;
    NodeMatrix node_matrix_;  ///< List of node neighbors by node ID (edges)
    SegmentsByNodes segments_by_nodes_;  ///< Lists of cells connecting nodes
};

}  // namespace pops

#endif  // POPS_NETWORK_HPP
