/*
 * PoPS model - Competency table for hosts and pest
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

#ifndef POPS_COMPETENCY_TABLE_HPP
#define POPS_COMPETENCY_TABLE_HPP

#include <vector>
#include <stdexcept>
#include <string>

#include "config.hpp"

namespace pops {

/**
 * Competency table holding combinations of host presences and absences and the
 * corresponding competency score.
 */
template<typename HostPool>
class CompetencyTable
{
public:
    using RasterIndex = typename HostPool::RasterIndex;
    using Environment = typename HostPool::Environment;

    /**
     * @brief Create an empty competency table
     *
     * @param environment Reference to the environment
     */
    CompetencyTable(const Environment& environment) : environment_(environment) {}

    /**
     * @brief Create a competency table using values in config
     *
     * @param config Configuration with competency table data
     * @param environment Reference to the environment
     */
    CompetencyTable(const Config& config, const Environment& environment)
        : environment_(environment)
    {
        for (const auto& row : config.competency_table_data()) {
            competency_table_.emplace_back(row.presence_absence, row.competency);
        }
    }

    /**
     * @brief Add competencies for a combination of host presences and absences
     *
     * Order of presence-absence data needs to be the same as host order in the
     * environment.
     *
     * Order of calls does not matter.
     *
     * @param presence_absence Presence (true) and absence (false) for each host.
     * @param competency Competency for a given presence-absence combination.
     */
    void
    add_host_competencies(const std::vector<bool>& presence_absence, double competency)
    {
        competency_table_.emplace_back(presence_absence, competency);
    }

    /**
     * @brief Get competency at a given cell
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param host Pointer to a specific host pool which is asking about competency
     *
     * @return Competency score
     *
     * @see find_competency()
     */
    double competency_at(RasterIndex row, RasterIndex col, const HostPool* host) const
    {
        auto presence_absence = environment_.host_presence_at(row, col);
        auto host_index = environment_.host_index(host);
        return find_competency(presence_absence, host_index);
    }

private:
    /**
     * @brief Find competency score in the table
     *
     * The *host_index* parameter is used to skip records which don't include the given
     * host.
     *
     * @param presence_absence Actual presence-absence data at a given place (cell)
     * @param host_index Index of the specific pool host originating the request
     *
     * @return Competency score
     */
    double
    find_competency(const std::vector<bool>& presence_absence, size_t host_index) const
    {
        // Go over all the rows and find the highest competency which fulfilled the
        // presence criteria.
        double competency = 0;
        for (const auto& table_row : competency_table_) {
            // probably faster if we just wait for the iteration over the whole thing
            if (!table_row.presence_absence[host_index])
                continue;
            if (presence_absence.size() != table_row.presence_absence.size()) {
                throw std::invalid_argument(
                    "Number of hosts in the environment is not the same as "
                    "the number of hosts in the competency table ("
                    + std::to_string(presence_absence.size())
                    + " != " + std::to_string(table_row.presence_absence.size()) + ")");
            }
            // Pick always the highest competency. Don't test row which has lower one.
            if (table_row.competency <= competency)
                continue;
            // Row is considered only if all required hosts are present.
            bool condition_fulfilled = true;
            for (size_t i = 0; i < presence_absence.size(); ++i) {
                if (!table_row.presence_absence[i]) {
                    continue;
                }
                if (!presence_absence[i]) {
                    condition_fulfilled = false;
                    break;
                }
            }
            if (condition_fulfilled)
                competency = table_row.competency;
        }
        return competency;
    }

    /**
     * @brief One row of the comptenecy table
     */
    struct TableRow
    {
        std::vector<bool> presence_absence;
        double competency;

        // Constructor can be removed with C++20, GCC 10, clang 16.
        TableRow(const std::vector<bool>& presence_absence, double competency)
            : presence_absence(presence_absence), competency(competency)
        {}
    };
    std::vector<TableRow> competency_table_;  ///< Internal table for lookup
    const Environment& environment_;  ///< Environment (for host index)
};

}  // namespace pops

#endif  // POPS_COMPETENCY_TABLE_HPP
