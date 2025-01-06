/*
 * PoPS model - treatments
 *
 * Copyright (C) 2015-2020 by the authors.
 *
 * Authors: Anna Petrasova <akratoc gmail com>
 *          Vaclav Petras <wenzeslaus gmail com>
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#ifndef POPS_TREATMENTS_HPP
#define POPS_TREATMENTS_HPP

#include "raster.hpp"
#include "date.hpp"
#include "scheduling.hpp"
#include "utils.hpp"
#include "host_pool.hpp"

#include <map>
#include <vector>
#include <string>
#include <functional>
#include <stdexcept>

// only temporarily for direct host pool creation
#include <random>

namespace pops {

/**
 * @brief The enum to decide how treatment is applied
 */
enum class TreatmentApplication
{
    Ratio,  ///< A ratio is applied to all treated rasters
    AllInfectedInCell  ///< All infected individuals are removed, rest by ratio
};

/**
 * @brief Get treatment application enum from string

 * Throws an std::invalid_argument exception if the value
 * is not supported.
 */
inline TreatmentApplication treatment_app_enum_from_string(const std::string& text)
{
    std::map<std::string, TreatmentApplication> mapping{
        {"ratio_to_all", TreatmentApplication::Ratio},
        {"ratio", TreatmentApplication::Ratio},
        {"all_infected_in_cell", TreatmentApplication::AllInfectedInCell},
        {"all infected", TreatmentApplication::AllInfectedInCell}};
    try {
        return mapping.at(text);
    }
    catch (const std::out_of_range&) {
        throw std::invalid_argument(
            "treatment_application_enum_from_string:"
            " Invalid value '"
            + text + "' provided");
    }
}

/*!
 * Abstract interface for treatment classes
 *
 * The class is meant for better internal code
 * layout and, at this point, it is not meant
 * as a universal matured interface for treatments.
 * Functions apply_treatment and end_treatment
 * are examples where we account for the current
 * concrete classes and will introduce more
 * general set of parameters only when
 * needed for additional classes.
 */
template<typename HostPool, typename FloatRaster>
class AbstractTreatment
{
public:
    virtual unsigned get_start() = 0;
    virtual unsigned get_end() = 0;
    virtual bool should_start(unsigned step) = 0;
    virtual bool should_end(unsigned step) = 0;
    virtual void apply_treatment(HostPool& host_pool) = 0;
    virtual void end_treatment(HostPool& host_pool) = 0;
    virtual ~AbstractTreatment() {}
};

/*!
 * Base treatment class.
 * Holds functions common between all treatment classes.
 */
template<typename HostPool, typename FloatRaster>
class BaseTreatment : public AbstractTreatment<HostPool, FloatRaster>
{
protected:
    unsigned start_step_;
    unsigned end_step_;
    FloatRaster map_;
    TreatmentApplication application_;

public:
    BaseTreatment(
        const FloatRaster& map,
        unsigned start,
        TreatmentApplication treatment_application)
        : start_step_(start),
          end_step_(start),
          map_(map),
          application_(treatment_application)
    {}
    unsigned get_start() override
    {
        return start_step_;
    }
    unsigned get_end() override
    {
        return end_step_;
    }

    // returning double allows identical results with the previous version
    int get_treated(int i, int j, int count)
    {
        return get_treated(i, j, count, this->application_);
    }

    int get_treated(int i, int j, int count, TreatmentApplication application)
    {
        if (application == TreatmentApplication::Ratio) {
            return std::lround(count * this->map_(i, j));
        }
        else if (application == TreatmentApplication::AllInfectedInCell) {
            return static_cast<bool>(this->map_(i, j)) ? count : 0;
        }
        throw std::runtime_error(
            "BaseTreatment::get_treated: unknown TreatmentApplication");
    }
};

/*!
 * Simple treatment class.
 * Removes percentage (given by treatment efficiency)
 * of infected and susceptible host (e.g. cut down trees).
 */
template<typename HostPool, typename FloatRaster>
class SimpleTreatment : public BaseTreatment<HostPool, FloatRaster>
{
public:
    SimpleTreatment(
        const FloatRaster& map,
        unsigned start,
        TreatmentApplication treatment_application)
        : BaseTreatment<HostPool, FloatRaster>(map, start, treatment_application)
    {}
    bool should_start(unsigned step) override
    {
        if (this->start_step_ == step)
            return true;
        return false;
    }
    bool should_end(unsigned) override
    {
        return false;
    }
    void apply_treatment(HostPool& host_pool) override
    {
        for (auto indices : host_pool.suitable_cells()) {
            int i = indices[0];
            int j = indices[1];
            int remove_susceptible = this->get_treated(
                i, j, host_pool.susceptible_at(i, j), TreatmentApplication::Ratio);
            // Treated infected are computed as a sum of treated in mortality groups.
            int remove_infected = 0;
            std::vector<int> remove_mortality;
            for (int count : host_pool.mortality_by_group_at(i, j)) {
                int remove = this->get_treated(i, j, count);
                remove_mortality.push_back(remove);
                remove_infected += remove;
            }
            // Will need to use infected directly if not mortality.

            std::vector<int> remove_exposed;
            for (int count : host_pool.exposed_by_group_at(i, j)) {
                remove_exposed.push_back(this->get_treated(i, j, count));
            }
            host_pool.completely_remove_hosts_at(
                i,
                j,
                remove_susceptible,
                remove_exposed,
                remove_infected,
                remove_mortality);
        }
    }
    void end_treatment(HostPool& host_pool) override
    {
        UNUSED(host_pool);
        return;
    }
};

/*!
 * Pesticide treatment class.
 * Removes percentage (given by treatment efficiency)
 * of infected and susceptible to resistant pool
 * and after certain number of days back to susceptible.
 */
template<typename HostPool, typename FloatRaster>
class PesticideTreatment : public BaseTreatment<HostPool, FloatRaster>
{
public:
    PesticideTreatment(
        const FloatRaster& map,
        unsigned start,
        unsigned end,
        TreatmentApplication treatment_application)
        : BaseTreatment<HostPool, FloatRaster>(map, start, treatment_application)
    {
        this->end_step_ = end;
    }
    bool should_start(unsigned step) override
    {
        if (this->start_step_ == step)
            return true;
        return false;
    }
    bool should_end(unsigned step) override
    {
        if (this->end_step_ == step)
            return true;
        return false;
    }
    void apply_treatment(HostPool& host_pool) override
    {
        for (auto indices : host_pool.suitable_cells()) {
            int i = indices[0];
            int j = indices[1];
            int susceptible_resistant = this->get_treated(
                i, j, host_pool.susceptible_at(i, j), TreatmentApplication::Ratio);
            std::vector<int> resistant_exposed_list;
            for (const auto& number : host_pool.exposed_by_group_at(i, j)) {
                resistant_exposed_list.push_back(this->get_treated(i, j, number));
            }
            int infected = 0;
            std::vector<int> resistant_mortality_list;
            for (const auto& number : host_pool.mortality_by_group_at(i, j)) {
                int remove = this->get_treated(i, j, number);
                resistant_mortality_list.push_back(remove);
                infected += remove;
            }
            host_pool.make_resistant_at(
                i,
                j,
                susceptible_resistant,
                resistant_exposed_list,
                infected,
                resistant_mortality_list);
        }
    }
    void end_treatment(HostPool& host_pool) override
    {
        for (auto indices : host_pool.suitable_cells()) {
            int i = indices[0];
            int j = indices[1];
            if (this->map_(i, j) > 0) {
                host_pool.remove_resistance_at(i, j);
            }
        }
    }
};

/*!
 * Treatments class manages all treatments.
 * Treatments can be simple (host removal)
 * and using pesticide (temporarily removed).
 * Each treatment can have unique date, type (simple, pesticide),
 * length (in case of pesticide), and treatment application.
 *
 * Pesticide treatments should not overlap spatially AND temporally
 * because of single resistance raster. In that case all resistant
 * populations that overlap both spatially and temporally
 * are returned to the susceptible pool when the first treatment ends.
 */
template<typename HostPool, typename FloatRaster>
class Treatments
{
private:
    std::vector<AbstractTreatment<HostPool, FloatRaster>*> treatments;
    Scheduler scheduler_;

public:
    Treatments(const Scheduler& scheduler) : scheduler_(scheduler) {}
    ~Treatments()
    {
        for (auto item : treatments) {
            delete item;
        }
    }
    /*!
     * \brief Add treatment, based on parameters it is distinguished
     * which treatment it will be.
     *
     * This works internally like a factory function
     * separating the user from all treatment classes.
     *
     * \param map treatment raster
     * \param start_date date when treatment is applied
     * \param num_days for simple treatments should be 0, otherwise number of days host
     *        is resistant
     * \param treatment_application if efficiency < 100% how should it be
     *        applied to infected/susceptible
     */
    void add_treatment(
        const FloatRaster& map,
        const Date& start_date,
        int num_days,
        TreatmentApplication treatment_application)
    {
        unsigned start = scheduler_.schedule_action_date(start_date);
        if (num_days == 0)
            treatments.push_back(new SimpleTreatment<HostPool, FloatRaster>(
                map, start, treatment_application));
        else {
            Date end_date(start_date);
            end_date.add_days(num_days);
            unsigned end = scheduler_.schedule_action_date(end_date);
            treatments.push_back(new PesticideTreatment<HostPool, FloatRaster>(
                map, start, end, treatment_application));
        }
    }
    /*!
     * \brief Do management if needed.
     * Should be called before every simulation step.
     * Decides internally whether any treatment needs to be
     * activated/deactivated.
     *
     * \param current simulation step
     * \param host_pool Host
     *
     * \return true if any management action was necessary
     */
    bool manage(unsigned current, HostPool& host_pool)
    {
        bool changed = false;
        for (unsigned i = 0; i < treatments.size(); i++) {
            if (treatments[i]->should_start(current)) {
                treatments[i]->apply_treatment(host_pool);
                changed = true;
            }
            else if (treatments[i]->should_end(current)) {
                treatments[i]->end_treatment(host_pool);
                changed = true;
            }
        }
        return changed;
    }
    /*!
     * \brief Used to remove treatments after certain step.
     * Needed for computational steering.
     * \param step simulation step
     */
    void clear_after_step(unsigned step)
    {
        for (auto& treatment : treatments) {
            if (treatment->get_start() > step) {
                delete treatment;
                treatment = nullptr;
            }
        }
        treatments.erase(
            std::remove(treatments.begin(), treatments.end(), nullptr),
            treatments.end());
    }
};

}  // namespace pops
#endif  // POPS_TREATMENTS_HPP
