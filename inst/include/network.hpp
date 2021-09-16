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

/**
 * Network structure and algorithms
 *
 * The network consists of nodes connected by edges. Edges are spatially represented
 * as segments. Edges themselves don't carry cost and have meaning only as indicators of
 * existence of a connection between nodes. Cost is travel time determined by advancing
 * through the segment.
 *
 * Nodes are the hop-on locations for dispersers. The dispersers can hop-off anywhere.
 *
 * The general workflow is contructing the object (with the constructor) and loading the
 * data (with the load() function). Then the network is ready to be used for simulating
 * trips over the network (with the travel() function).
 *
 * When the travel() function is used from a kernel, user of the network directly calls
 * only the setup functions.
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
     * @param speed Speed of travel in the same units as coordinates
     */
    Network(BBox<double> bbox, double ew_res, double ns_res, double speed)
        : bbox_(bbox),
          ew_res_(ew_res),
          ns_res_(ns_res),
          cell_travel_time_(((ew_res + ns_res) / 2) / speed)
    {}

    /**
     * @brief Create an empty network not meant for futher use.
     *
     * This is useful when a network object is needed to contruct a kernel, but it will
     * not be used in runtime.
     *
     * @return New Network object
     */
    static Network null_network()
    {
        return Network(BBox<double>(), 0, 0, 0);
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
        return xy_to_row_col(std::stod(x), std::stod(y));
    }

    /**
     * @brief Test if coordinates XY are in the bounding box.
     * @param x Real world coordinate X
     * @param y Real world coordinate Y
     * @return True if XY is in the bounding box, false otherwise.
     */
    bool out_of_bbox(double x, double y) const
    {
        return x > bbox_.east || x < bbox_.west || y > bbox_.north || y < bbox_.south;
    }

    /**
     * @brief Load nodes and segments from input streams.
     *
     * @param node_stream Input stream with records for nodes
     * @param segment_stream Input stream with records for segments
     * @param allow_empty True if the loaded network can be empty (no nodes)
     *
     * The nodes stream is a CSV with rows `node_id,X,Y`.
     * Segment stream is custom text format resembling CSV geared towards parsing, not
     * readability, with rows `node_id_1,node_id_2,X1;Y1;X2;Y2;X3;Y3;...`.
     *
     * The function can take any input stream which behaves like std::ifstream or
     * std::istringstream. These are also the two expected streams to be used
     * (for reading from a file and reading from a string, respectively).
     *
     * A simple code for reading from files without any checking may look like this:
     *
     * ```
     * Network network(...);
     * std::string node_filename = "nodes.txt";
     * std::string segment_filename = "segments.txt"
     * std::ifstream node_stream{node_filename};
     * std::ifstream segment_stream{segment_filename};
     * network.load(node_stream, segment_stream);
     * ```
     *
     * The input network is a list of nodes and their coordinates and a list of edges
     * and the corresponding segments. The input is used as-is. The network is agnostic
     * towards how the segments look like in terms of crossing each other or other
     * nodes.
     *
     * All XY coordinates of nodes and segment points are converted to row and column
     * using bounding box and resolution.
     *
     * Coordinates for nodes and segment are treated separately, so input segment
     * coordinates should include the node coordinates.
     *
     * Standalone nodes (without any connection) are allowed and handled as no movement
     * from the source cell. Only edges (segments) with both nodes in the bounding box
     * are considered, i.e., input network is reduced to the bounding box during
     * reading, but clipping happens on the level of whole edges, not in the middle of
     * a segment.
     */
    template<typename InputStream>
    void load(
        InputStream& node_stream, InputStream& segment_stream, bool allow_empty = false)
    {
        std::set<NodeId> node_ids;
        load_nodes(node_stream, node_ids);
        if (node_ids.empty()) {
            if (allow_empty)
                return;
            else
                throw std::runtime_error("Network: No nodes within the extend");
        }
        load_segments(segment_stream, node_ids);
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
            throw std::invalid_argument("No nodes at a given row and column");
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
     * Travel in the network from given row and column for a given time.
     *
     * All previously visited nodes are tracked and, if possible, excluded
     * from further traveling.
     *
     * @returns Final row and column pair
     */
    template<typename Generator>
    std::tuple<int, int>
    travel(RasterIndex row, RasterIndex col, double time, Generator& generator) const
    {
        // We assume there is a node here, i.e., that we are made decision
        // to use this kernel knowing there is a node.
        auto node_id = get_random_node_at(row, col, generator);
        std::set<NodeId> visited_nodes;
        while (time >= 0) {
            auto next_node_id = next_node(node_id, visited_nodes, generator);
            visited_nodes.insert(node_id);
            // If there is no segment from the node, return the start cell.
            if (next_node_id == node_id)
                return std::make_tuple(row, col);
            auto segment = get_segment(node_id, next_node_id);
            // nodes may need special handling
            for (const auto& cell : segment) {
                time -= cell_travel_time_;
                if (time <= 0) {
                    return cell;
                    // Given the while condition, this subsequently ends the while loop
                    // as well.
                    // break;
                }
            }
            node_id = next_node_id;
        }
        throw std::invalid_argument("Time must be greater than or equal to zero");
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
        // We store segements in both directions, so each segment is stored twice.
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
        stream << "    cell_travel_time: " << cell_travel_time_ << "\n";
        std::set<std::pair<NodeId, NodeId>> edges;
        for (const auto& item : node_matrix_) {
            for (const auto& node_id : item.second) {
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
     * Connection is both directions is stored explicitly (see load_segments() source
     * code).
     */
    using NodeMatrix = std::map<NodeId, std::vector<NodeId>>;
    /** Cells connecting two nodes (segment between nodes) */
    using Segment = std::vector<std::pair<RasterIndex, RasterIndex>>;
    /** Constant view of a segment (to iterate a segment in either direction) */
    using SegmentView = ContainerView<Segment>;
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
        return std::stoi(text);
    }

    /**
     * @brief Read nodes from a stream.
     *
     * @param stream Input stream with nodes and their coordinates
     * @param[out] node_ids Container to store node IDs in
     */
    template<typename InputStream>
    void load_nodes(InputStream& stream, std::set<NodeId>& node_ids)
    {
        std::string line;
        while (std::getline(stream, line)) {
            std::istringstream line_stream{line};
            char delimeter{','};
            std::string node_text;
            std::getline(line_stream, node_text, delimeter);
            std::string x_coord_text;
            std::getline(line_stream, x_coord_text, delimeter);
            std::string y_coord_text;
            std::getline(line_stream, y_coord_text, delimeter);
            RasterIndex row;
            RasterIndex col;
            double x = std::stod(x_coord_text);
            double y = std::stod(y_coord_text);
            // Cut to extend
            if (out_of_bbox(x, y))
                continue;
            std::tie(row, col) = xy_to_row_col(x_coord_text, y_coord_text);
            NodeId node_id = node_id_from_text(node_text);
            if (node_id < 1) {
                std::runtime_error("Node ID must be greater than zero");
            }
            nodes_by_row_col_[std::make_pair(row, col)].insert(node_id);
            node_ids.insert(node_id);
        }
    }

    /**
     * @brief Read segments from a stream.
     *
     * @param stream Input stream with segments
     * @param node_ids Container with node IDs to use
     */
    template<typename InputStream>
    void load_segments(InputStream& stream, const std::set<NodeId>& node_ids)
    {
        std::set<std::pair<NodeId, NodeId>> node_pairs;
        std::string line;
        while (std::getline(stream, line)) {
            std::istringstream line_stream{line};
            char delimeter{','};
            std::string node_1_text;
            std::getline(line_stream, node_1_text, delimeter);
            std::string node_2_text;
            std::getline(line_stream, node_2_text, delimeter);
            auto node_1_id = node_id_from_text(node_1_text);
            auto node_2_id = node_id_from_text(node_2_text);
            // If either end nodes of the segment is not in the extent, skip it.
            // This means that a segment is ignored even if one of the nodes and
            // significant portion of the segment is in the area of iterest.
            // TODO: Part of the segment may still be out, so that needs to be checked.
            // Replace find by contains for C++20.
            if (node_ids.find(node_1_id) == node_ids.end()
                || node_ids.find(node_2_id) == node_ids.end())
                continue;
            if (node_1_id == node_2_id) {
                std::runtime_error(
                    std::string("Segment cannot begin and end with the same node: ")
                    + node_1_text + " " + node_2_text);
            }
            std::string segment_text;
            std::getline(line_stream, segment_text, delimeter);
            // We don't know which way the nodes are ordered, so instead of checking
            // the order, we create a symmetric matrix since we allocated the memory
            // anyway.
            node_pairs.emplace(node_1_id, node_2_id);
            std::istringstream segment_stream{segment_text};
            char in_cell_delimeter{';'};
            std::string x_coord_text;
            std::string y_coord_text;
            Segment segment;
            while (std::getline(segment_stream, x_coord_text, in_cell_delimeter)
                   && std::getline(segment_stream, y_coord_text, in_cell_delimeter)) {
                // The same cell is possibly repeated if raster resolution is lower than
                // the detail ("resolution") of each segment, so we check if the last
                // added coordinate pair converted to same cell.
                auto new_point = xy_to_row_col(x_coord_text, y_coord_text);
                if (segment.empty() || segment.back() != new_point)
                    segment.emplace_back(new_point);
            }
            segments_by_nodes_.emplace(
                std::make_pair(node_1_id, node_2_id), std::move(segment));
        }

        for (auto node_id : node_ids) {
            std::vector<NodeId> nodes;
            for (const auto& item : node_pairs) {
                if (item.first == node_id)
                    nodes.push_back(item.second);
                else if (item.second == node_id)
                    nodes.push_back(item.first);
            }
            node_matrix_[node_id] = std::move(nodes);
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
        for (const auto& item : segments_by_nodes_) {
            const auto& key{item.first};
            if (key.first == start && key.second == end)
                return SegmentView(item.second.cbegin(), item.second.cend());
            if (key.second == start && key.first == end) {
                return SegmentView(item.second.crbegin(), item.second.crend());
            }
        }
        throw std::invalid_argument("No segment for given nodes");
    }

    /**
     * @brief Get nodes connected by an edge to a given node.
     *
     * @param node Node to get connections from
     * @return List of connected nodes
     */
    const std::vector<NodeId>& nodes_connected_to(NodeId node) const
    {
        return node_matrix_.at(node);
    }

    /**
     * @brief Pick a next node from the given node.
     *
     * If there is more than one edge leading from the given node, a random node is
     * picked. If there are no edges leading from the given node, the node is returned.
     *
     * The random node is picked from candidate nodes which are nodes connected to
     * a given node. The candidate nodes which are in the *ignore* list are excluded
     * from the random selection. If all candiate nodes are in the *ignore* list,
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

    BBox<double> bbox_;  ///< Bounding box of the network grid in real world coordinates
    double ew_res_;  ///< East-west resolution of the grid
    double ns_res_;  ///< North-south resolution of the grid
    double cell_travel_time_;  ///< Time to travel through one cell
    /** Node IDs stored by row and column (multiple nodes per cell) */
    std::map<std::pair<RasterIndex, RasterIndex>, std::set<NodeId>> nodes_by_row_col_;
    NodeMatrix node_matrix_;  ///< List of node neighbors by node ID (edges)
    SegmentsByNodes segments_by_nodes_;  ///< Lists of cells connecting nodes
};

}  // namespace pops

#endif  // POPS_NETWORK_HPP
