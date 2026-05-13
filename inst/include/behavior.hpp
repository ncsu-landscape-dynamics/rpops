/*
 * PoPS model - grower/manager behavior module
 *
 * Copyright (C) 2024 by the authors.
 *
 * Authors: Rachel Seibel
 *
 * This file is part of PoPS.
 *
 * PoPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * PoPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PoPS. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef POPS_BEHAVIOR_HPP
#define POPS_BEHAVIOR_HPP

#include "date.hpp"
#include "scheduling.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace pops {

/**
 * @brief Per-grower decision parameters.
 *
 * These are the global behavioral parameters that govern how each grower
 * responds to perceived disease pressure. They are stored inline in the main
 * YAML config under `behavior_params` and loaded into a lookup table keyed
 * by grower_id (the integer value in the grower_id raster).
 *
 * All probability fields must be in [0, 1].
 */
struct GrowerParams {
    /** Probability that an infected cell within the management unit is detected. */
    double detection_prob{1.0};
    /**
     * Probability of deciding to treat given that perceived prevalence
     * exceeds decision_threshold.
     */
    double willingness_to_treat{0.5};
    /** Fractional reduction in infection applied by the treatment. */
    double treatment_efficacy{1.0};
    /**
     * Perceived prevalence (fraction of host cells that appear infected)
     * above which the grower considers treating.
     */
    double decision_threshold{0.0};
};

/**
 * @brief Grower/manager behavior module for PoPS.
 *
 * Activated by `use_behavior_module: true` in the YAML config; disabled by
 * default. When active, the module intercepts the model loop at each
 * configured decision date, computes unit-level prevalence, applies a
 * stochastic perception filter and decision rule, and injects a treatment
 * map back into the host pool.
 *
 * ## Data model
 * - **grower_id raster** (`grower_id_`): integer raster (same extent/resolution
 *   as the host raster) with one integer ID per management unit cell; 0 = no
 *   unit (background).  Produced by `delineate_management_units()` in R and
 *   saved as a GeoTIFF, then read into PoPS alongside the host raster.
 * - **parameter table** (`params_`): `std::map<int, GrowerParams>` keyed by
 *   grower ID.  Populated at construction from the YAML `behavior_params`
 *   list via `pops.cpp`.
 *
 * ## Decision loop (called from `pops.cpp` model step loop)
 * At each step the model checks `is_decision_step(step)`.  When true,
 * `generate_treatment_map()` is called, which:
 *   1. Aggregates infected and host cell counts per management unit.
 *   2. Applies binomial thinning with `detection_prob` to obtain perceived
 *      prevalence.
 *   3. Compares to `decision_threshold`; if exceeded, treats with probability
 *      `willingness_to_treat`.
 *   4. Returns a FloatRaster of treatment efficacy values (0 or
 *      `treatment_efficacy` per cell) to be fed into `Treatments::apply()`.
 *
 * @tparam IntegerRaster Raster type for integer data (infected, host counts,
 *   grower IDs).
 * @tparam FloatRaster   Raster type for floating-point data (treatment maps).
 * @tparam Generator     Random number generator type (e.g.
 *   `std::default_random_engine`).
 */
template<typename IntegerRaster, typename FloatRaster, typename Generator>
class BehaviorModule
{
public:
    /**
     * @brief Construct the behavior module.
     *
     * @param grower_id_raster Integer raster mapping cells to management unit
     *   IDs (0 = background). Must have the same dimensions as the host raster.
     * @param params Lookup table of per-grower decision parameters, keyed by
     *   the integer IDs present in @p grower_id_raster.
     * @param decision_dates Sorted vector of simulation dates on which growers
     *   make treatment decisions.  Dates outside [start_date, end_date] are
     *   silently ignored.
     */
    BehaviorModule(
        const IntegerRaster& grower_id_raster,
        const std::map<int, GrowerParams>& params,
        const std::vector<Date>& decision_dates)
        : grower_id_(grower_id_raster),
          params_(params),
          decision_dates_(decision_dates)
    {}

    /**
     * @brief Returns true if @p step contains a configured decision date.
     *
     * Called from the main simulation loop in `pops.cpp` to gate the
     * (relatively expensive) prevalence aggregation.
     */
    bool is_decision_step(const Step& step) const
    {
        for (const auto& d : decision_dates_) {
            if (d >= step.start_date() && d <= step.end_date())
                return true;
        }
        return false;
    }

    /**
     * @brief Compute per-unit prevalence and generate a treatment map.
     *
     * @param infected Current infected-host counts (IntegerRaster).
     * @param total_hosts Current total-host counts (IntegerRaster).
     * @param generator Random number generator for stochastic decisions.
     * @return FloatRaster with treatment efficacy values:
     *   - 0.0 where no treatment is applied (background or untreated units)
     *   - `GrowerParams::treatment_efficacy` where the unit decides to treat
     */
    FloatRaster generate_treatment_map(
        const IntegerRaster& infected,
        const IntegerRaster& total_hosts,
        Generator& generator) const
    {
        const int rows = infected.rows();
        const int cols = infected.cols();
        FloatRaster tmap(rows, cols, 0.0);

        // --- Step 1: aggregate infected and host cell counts per unit --------
        // Use std::map so we only iterate over unique IDs found in the raster.
        std::map<int, long long> unit_infected;
        std::map<int, long long> unit_hosts;

        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int uid = grower_id_(r, c);
                if (uid == 0) continue;  // background cell
                unit_infected[uid] += static_cast<long long>(infected(r, c));
                unit_hosts[uid]    += static_cast<long long>(total_hosts(r, c));
            }
        }

        // --- Step 2: decision per unit ----------------------------------------
        std::uniform_real_distribution<double> uniform(0.0, 1.0);

        for (auto& kv : unit_hosts) {
            int uid = kv.first;
            long long n_hosts    = kv.second;
            long long n_infected = unit_infected.count(uid) ? unit_infected[uid] : 0LL;

            if (n_hosts == 0) continue;

            auto it = params_.find(uid);
            if (it == params_.end()) continue;  // no params for this unit
            const GrowerParams& p = it->second;

            // Apply detection probability (binomial thinning)
            long long observed_infected = n_infected;
            if (p.detection_prob < 1.0 && n_infected > 0) {
                std::binomial_distribution<long long> binom(n_infected, p.detection_prob);
                observed_infected = binom(generator);
            }
            double perceived_prevalence =
                static_cast<double>(observed_infected) / static_cast<double>(n_hosts);

            // Decision rule
            bool treats = false;
            if (perceived_prevalence > p.decision_threshold) {
                treats = (uniform(generator) < p.willingness_to_treat);
            }

            if (!treats) continue;

            // --- Step 3: write efficacy to all cells in this unit -------------
            double eff = p.treatment_efficacy;
            for (int r = 0; r < rows; ++r) {
                for (int c = 0; c < cols; ++c) {
                    if (grower_id_(r, c) == uid && total_hosts(r, c) > 0)
                        tmap(r, c) = static_cast<typename FloatRaster::value_type>(eff);
                }
            }
        }

        return tmap;
    }

private:
    IntegerRaster grower_id_;
    std::map<int, GrowerParams> params_;
    std::vector<Date> decision_dates_;
};

}  // namespace pops

#endif  // POPS_BEHAVIOR_HPP
