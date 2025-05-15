/*
 * PoPS model - pest or pathogen spread simulation
 *
 * Copyright (C) 2015-2023 by the authors.
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

#ifndef POPS_SIMULATION_HPP
#define POPS_SIMULATION_HPP

#include <cmath>
#include <memory>
#include <tuple>
#include <vector>
#include <random>
#include <string>
#include <stdexcept>

#include "actions.hpp"
#include "utils.hpp"
#include "soils.hpp"
#include "model_type.hpp"
#include "pest_pool.hpp"
#include "host_pool.hpp"

namespace pops {

/*! A class to control the spread simulation.
 *
 * @note
 * The class is deprecated for external use in favor of individual action classes and a
 * higher-level Model. The class corresponding to the original Simulation class before
 * too much code accumulated in Simulation is SpreadAction. The class is now used only
 * in tests.
 *
 * The Simulation class handles the mechanics of the model, but the
 * timing of the events or steps should be handled outside of this
 * class unless noted otherwise. The notable exceptions are exposed
 * hosts in the SEI model type and mortality.
 *
 * The template parameters IntegerRaster and FloatRaster are raster
 * image or matrix types. Any 2D numerical array should work as long as
 * it uses function call operator to access the values, i.e. it provides
 * indexing for reading and writing values using `()`. In other words,
 * the operations such as the two following ones should be possible:
 *
 * ```
 * a(i, j) = 1;
 * a(i, j) == 1;
 * ```
 *
 * The PoPS library offers a Raster template class to fill this role,
 * but other classes can be used as well.
 *
 * Template parameter RasterIndex is type used for maximum indices of
 * the used rasters and should be the same as what the actual raster
 * types are using. However, at the same time, comparison with signed
 * type are performed and a signed type might be required in the future.
 * A default is provided, but it can be changed in the future.
 */
template<
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex = int,
    typename Generator = DefaultSingleGeneratorProvider>
class Simulation
{
private:
    RasterIndex rows_;
    RasterIndex cols_;
    bool dispersers_stochasticity_;
    bool establishment_stochasticity_;
    bool movement_stochasticity_;
    ModelType model_type_;
    unsigned latency_period_;
    /// Non-owning pointer to environment for weather
    // can be const in simulation, except for need to add host now
    Environment<IntegerRaster, FloatRaster, RasterIndex, Generator>* environment_{
        nullptr};
    /**
     * Optional soil pool
     */
    std::shared_ptr<SoilPool<IntegerRaster, FloatRaster, RasterIndex, Generator>>
        soil_pool_{nullptr};
    /**
     * Percentage (0-1 ratio) of disperers to be send to soil
     */
    double to_soil_percentage_{0};

public:
    // Host pool has the provider from model, but in test, it gets plain engine.
    using StandardHostPool =
        HostPool<IntegerRaster, FloatRaster, RasterIndex, Generator>;
    using StandardPestPool = PestPool<IntegerRaster, FloatRaster, RasterIndex>;

    /**
     * Creates simulation object with the values which are fixed during the simulation.
     *
     * The number or rows and columns needs to be the same as the size
     * of rasters used with the Simulation object
     * (potentially, it can be also smaller).
     *
     * @param model_type Type of the model (SI or SEI)
     * @param latency_period Length of the latency period in steps
     * @param rows Number of rows
     * @param cols Number of columns
     * @param dispersers_stochasticity Enable stochasticity in generating of dispersers
     * @param establishment_stochasticity Enable stochasticity in establishment step
     * @param movement_stochasticity Enable stochasticity in movement of hosts
     */
    Simulation(
        RasterIndex rows,
        RasterIndex cols,
        ModelType model_type = ModelType::SusceptibleInfected,
        unsigned latency_period = 0,
        bool dispersers_stochasticity = true,
        bool establishment_stochasticity = true,
        bool movement_stochasticity = true)
        : rows_(rows),
          cols_(cols),
          dispersers_stochasticity_(dispersers_stochasticity),
          establishment_stochasticity_(establishment_stochasticity),
          movement_stochasticity_(movement_stochasticity),
          model_type_(model_type),
          latency_period_(latency_period)
    {}

    Simulation() = delete;

    /**
     * @brief Set environment used for weather to provided environment
     * @param environment Pointer to an existing environment
     *
     * The simulation object does not take ownership of the environment.
     */
    // parameter and attribute should be const
    void set_environment(
        Environment<IntegerRaster, FloatRaster, RasterIndex, Generator>* environment)
    {
        this->environment_ = environment;
    }

    /**
     * @brief Get environment used in the simulation
     *
     * @param allow_empty if true, empty (non-functional) environment is returned
     * @return Const pointer to the environment
     * @throw std::logic_error when environment is not set
     */
    Environment<IntegerRaster, FloatRaster, RasterIndex, Generator>*
    environment(bool allow_empty = false)
    {
        static Environment<IntegerRaster, FloatRaster, RasterIndex, Generator> empty;
        if (!this->environment_) {
            if (allow_empty)
                return &empty;
            throw std::logic_error("Environment used in Simulation, but not provided");
        }
        return this->environment_;
    }

    /** Remove infected based on temperature tolerance
     *
     * @param infected Currently infected hosts
     * @param susceptible Currently susceptible hosts
     * @param exposed Currently exposed hosts
     * @param total_exposed Total exposed in all exposed cohorts
     * @param mortality_tracker_vector Mortality tracker
     * @param lethal_temperature temperature at which lethal conditions occur
     * @param suitable_cells used to run model only where host are known to occur
     * @param generator Provider of random number generators
     */
    template<typename GeneratorProvider>
    void remove(
        IntegerRaster& infected,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& total_exposed,
        std::vector<IntegerRaster>& mortality_tracker_vector,
        double lethal_temperature,
        std::vector<std::vector<int>>& suitable_cells,
        GeneratorProvider& generator)
    {
        IntegerRaster empty;
        StandardHostPool hosts(
            model_type_,
            susceptible,
            exposed,
            0,
            infected,
            total_exposed,
            empty,
            mortality_tracker_vector,
            empty,
            empty,
            *environment(true),
            false,
            0,
            false,
            0,
            0,
            0,
            suitable_cells);
        RemoveByTemperature<
            StandardHostPool,
            IntegerRaster,
            FloatRaster,
            RasterIndex,
            GeneratorProvider>
            remove(*environment(false), lethal_temperature);
        remove.action(hosts, generator);
        this->environment(true)->remove_hosts();
    }

    /** Removes percentage of exposed and infected
     *
     * @param infected Currently infected hosts
     * @param susceptible Currently susceptible hosts
     * @param mortality_tracker_vector Hosts that are infected at a specific time step
     * @param exposed Exposed hosts per cohort
     * @param total_exposed Total exposed in all exposed cohorts
     * @param survival_rate Raster between 0 and 1 representing pest survival rate
     * @param suitable_cells used to run model only where host are known to occur
     * @param generator Provider of random number generators
     */
    template<typename GeneratorProvider>
    void remove_percentage(
        IntegerRaster& infected,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& mortality_tracker_vector,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& total_exposed,
        const FloatRaster& survival_rate,
        std::vector<std::vector<int>>& suitable_cells,
        GeneratorProvider& generator)
    {
        IntegerRaster empty;
        StandardHostPool hosts(
            model_type_,
            susceptible,
            exposed,
            0,
            infected,
            total_exposed,
            empty,
            mortality_tracker_vector,
            empty,
            empty,
            *environment(true),
            false,
            0,
            false,
            0,
            0,
            0,
            suitable_cells);
        SurvivalRateAction<StandardHostPool, IntegerRaster, FloatRaster> survival(
            survival_rate);
        survival.action(hosts, generator);
        this->environment(true)->remove_hosts();
    }

    /** kills infected hosts based on mortality rate and timing. In the last year
     * of mortality tracking the first index all remaining tracked infected hosts
     * are removed. In indexes that are in the mortality_time_lag no mortality occurs.
     * In all other indexes the number of tracked individuals is multiplied by the
     * mortality rate to calculate the number of hosts that die that time step. The
     * mortality_tracker_vector has a minimum size of mortality_time_lag + 1.
     *
     * @param infected Currently infected hosts
     * @param total_hosts All hosts
     * @param mortality_rate percent of infected hosts that die each time period
     * @param mortality_time_lag time lag prior to mortality beginning
     * @param died dead hosts during time step
     * @param mortality_tracker_vector vector of matrices for tracking infected
     * host infection over time. Expectation is that mortality tracker is of
     * length (1/mortality_rate + mortality_time_lag)
     * @param suitable_cells used to run model only where host are known to occur
     */
    void mortality(
        IntegerRaster& infected,
        IntegerRaster& total_hosts,
        double mortality_rate,
        int mortality_time_lag,
        IntegerRaster& died,
        std::vector<IntegerRaster>& mortality_tracker_vector,
        std::vector<std::vector<int>>& suitable_cells)
    {
        IntegerRaster empty;
        std::vector<IntegerRaster> empty_vector;
        StandardHostPool hosts{
            model_type_,
            empty,
            empty_vector,
            0,
            infected,
            empty,
            empty,
            mortality_tracker_vector,
            died,
            total_hosts,
            *environment(true),
            false,
            0,
            false,
            0,
            0,
            0,
            suitable_cells};
        Mortality<StandardHostPool, IntegerRaster, FloatRaster> mortality(
            mortality_rate, mortality_time_lag);
        mortality.action(hosts);
        this->environment(true)->remove_hosts();
    }

    /** Moves hosts from one location to another
     *
     * @note Note that unlike the other functions, here, *total_hosts*,
     * i.e., number of hosts is required, not number of all hosts
     * and non-host individuals.
     *
     * @param infected Currently infected hosts
     * @param susceptible Currently susceptible hosts
     * @param mortality_tracker_vector Hosts that are infected at a specific time step
     * @param total_hosts All host individuals in the area. Is equal to
     *        infected + exposed + susceptible in the cell.
     * @param total_exposed Total exposed in all exposed cohorts
     * @param exposed Exposed hosts per cohort
     * @param resistant Resistant hosts
     * @param step the current step of the simulation
     * @param last_index the last index to not be used from movements
     * @param movements a vector of ints with row_from, col_from, row_to, col_to, and
     *        num_hosts
     * @param movement_schedule a vector matching movements with the step at which the
     *        movement from movements are applied
     * @param suitable_cells List of indices of cells with hosts
     * @param generator Provider of random number generators
     *
     * @note Mortality and non-host individuals are not supported in movements.
     */
    unsigned movement(
        IntegerRaster& infected,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& mortality_tracker_vector,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& resistant,
        IntegerRaster& total_hosts,
        IntegerRaster& total_exposed,
        unsigned step,
        unsigned last_index,
        const std::vector<std::vector<int>>& movements,
        std::vector<unsigned> movement_schedule,
        std::vector<std::vector<int>>& suitable_cells,
        Generator& generator)
    {
        HostMovement<StandardHostPool, IntegerRaster, FloatRaster, RasterIndex>
            host_movement{step, last_index, movements, movement_schedule};
        IntegerRaster empty;
        StandardHostPool hosts{
            model_type_,
            susceptible,
            exposed,
            0,
            infected,
            total_exposed,
            resistant,
            mortality_tracker_vector,
            empty,
            total_hosts,
            *environment(true),
            false,
            0,
            false,
            0,
            0,
            0,
            suitable_cells};
        auto ret = host_movement.action(hosts, generator);
        this->environment(true)->remove_hosts();
        return ret;
    }

    /** Generates dispersers based on infected
     *
     * @param[out] dispersers  (existing values are ignored)
     * @param[out] established_dispersers Dispersers for a cell established in
     * another cell (later modified to final form in disperse())
     * @param infected Currently infected hosts
     * @param weather Whether to use the weather coefficient
     * @param reproductive_rate reproductive rate (used unmodified when weather
     *        coefficient is not used)
     * @param[in] suitable_cells List of indices of cells with hosts
     * @param generator Provider of random number generators
     */
    void generate(
        IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        const IntegerRaster& infected,
        bool weather,
        double reproductive_rate,
        const std::vector<std::vector<int>>& suitable_cells,
        Generator& generator)
    {
        IntegerRaster empty;
        std::vector<IntegerRaster> empty_vector;
        StandardHostPool host_pool{
            model_type_,
            empty,
            empty_vector,
            0,
            const_cast<IntegerRaster&>(infected),
            empty,
            empty,
            empty_vector,
            empty,
            empty,
            *environment(!weather),
            dispersers_stochasticity_,
            reproductive_rate,
            false,
            0,
            0,
            0,
            const_cast<std::vector<std::vector<int>>&>(suitable_cells)};
        std::vector<std::tuple<int, int>> empty_outside_dispersers;
        StandardPestPool pests{
            dispersers, established_dispersers, empty_outside_dispersers};
        std::default_random_engine unused_kernel;
        SpreadAction<
            StandardHostPool,
            StandardPestPool,
            IntegerRaster,
            FloatRaster,
            RasterIndex,
            std::default_random_engine,
            Generator>
            spread_action{unused_kernel};
        spread_action.activate_soils(soil_pool_, to_soil_percentage_);
        spread_action.generate(host_pool, pests, generator);
        this->environment(!weather)->remove_hosts();
    }

    /** Creates dispersal locations for the dispersing individuals
     *
     * Depending on what data is provided as the *exposed_or_infected*
     * parameter, this function can be part of the S to E step or the
     * S to I step.
     *
     * Typically, the generate() function is called beforehand to
     * create dispersers. In SEI model, the infect_exposed() function is
     * typically called afterwards.
     *
     * DispersalKernel is callable object or function with one parameter
     * which is the random number engine (generator). The return value
     * is row and column in the raster (or outside of it). The current
     * position is passed as parameters. The return value is in the
     * form of a tuple with row and column so that std::tie() is usable
     * on the result, i.e. function returning
     * `std::make_tuple(row, column)` fulfills this requirement.
     *
     * The *total_populations* can be total number of hosts in the basic case
     * or it can be the total size of population of all relevant species
     * both host and non-host if dilution effect should be applied.
     *
     * If establishment stochasticity is disabled,
     * *establishment_probability* is used to decide whether or not
     * a disperser is established in a cell. Value 1 means that all
     * dispersers will establish and value 0 means that no dispersers
     * will establish.
     *
     * @param[in] dispersers Dispersing individuals ready to be dispersed
     * @param[in,out] established_dispersers Dispersers for a cell established in
     * another cell
     * @param[in,out] susceptible Susceptible hosts
     * @param[in,out] exposed Exposed hosts
     * @param[in,out] infected Infected hosts
     * @param[in,out] mortality_tracker Newly infected hosts (if applicable)
     * @param[in,out] total_exposed Total exposed in all exposed cohorts
     * @param[in] total_populations All host and non-host individuals in the area
     * @param[in,out] outside_dispersers Dispersers escaping the raster
     * @param weather Whether or not weather coefficients should be used
     * @param dispersal_kernel Dispersal kernel to move dispersers
     * @param establishment_probability Probability of establishment with no
     *        stochasticity
     * @param[in] suitable_cells List of indices of cells with hosts
     * @param generator Provider of random number generators
     *
     * @note If the parameters or their default values don't correspond
     * with the disperse_and_infect() function, it is a bug.
     */
    template<typename DispersalKernel>
    void disperse(
        IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& infected,
        std::vector<IntegerRaster>& mortality_tracker,
        const IntegerRaster& total_populations,
        IntegerRaster& total_exposed,
        std::vector<std::tuple<int, int>>& outside_dispersers,
        bool weather,
        DispersalKernel& dispersal_kernel,
        std::vector<std::vector<int>>& suitable_cells,
        double establishment_probability,
        Generator& generator)
    {
        // The interaction does not happen over the member variables yet, use empty
        // variables. This requires SI/SEI to be fully resolved in host and not in
        // disperse_and_infect.
        IntegerRaster empty;
        StandardHostPool host_pool{
            model_type_,
            susceptible,
            exposed,
            0,
            infected,
            total_exposed,
            empty,
            mortality_tracker,
            empty,
            empty,
            *environment(!weather),
            false,
            0,
            establishment_stochasticity_,
            establishment_probability,
            rows_,
            cols_,
            suitable_cells};
        // This would be part of the main initialization process.
        if (environment_) {
            environment_->set_total_population(&total_populations);
        }
        StandardPestPool pests{dispersers, established_dispersers, outside_dispersers};
        SpreadAction<
            StandardHostPool,
            StandardPestPool,
            IntegerRaster,
            FloatRaster,
            RasterIndex,
            DispersalKernel,
            Generator>
            spread_action{dispersal_kernel};
        spread_action.activate_soils(soil_pool_, to_soil_percentage_);
        spread_action.disperse(host_pool, pests, generator);
        this->environment(!weather)->remove_hosts();
    }

    // For backwards compatibility for tests (without exposed and mortality)
    template<typename DispersalKernel>
    void disperse(
        const IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        IntegerRaster& susceptible,
        IntegerRaster& infected,
        IntegerRaster& mortality_tracker,
        const IntegerRaster& total_populations,
        IntegerRaster& total_exposed,
        std::vector<std::tuple<int, int>>& outside_dispersers,
        bool weather,
        DispersalKernel& dispersal_kernel,
        std::vector<std::vector<int>>& suitable_cells,
        double establishment_probability = 0.5)
    {
        std::vector<IntegerRaster> tmp;
        tmp.push_back(mortality_tracker);
        std::vector<IntegerRaster> empty_vector;
        disperse(
            dispersers,
            established_dispersers,
            susceptible,
            empty_vector,
            infected,
            tmp,  // mortality
            total_populations,
            total_exposed,
            outside_dispersers,
            weather,
            dispersal_kernel,
            suitable_cells,
            establishment_probability);
        mortality_tracker = tmp.back();
    }

    // For backwards compatibility for tests (without exposed and mortality)
    template<typename DispersalKernel>
    void disperse(
        IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        IntegerRaster& susceptible,
        IntegerRaster& infected,
        IntegerRaster& mortality_tracker,
        const IntegerRaster& total_populations,
        IntegerRaster& total_exposed,
        std::vector<std::tuple<int, int>>& outside_dispersers,
        bool weather,
        DispersalKernel& dispersal_kernel,
        std::vector<std::vector<int>>& suitable_cells,
        double establishment_probability,
        Generator& generator)
    {
        std::vector<IntegerRaster> tmp;
        tmp.push_back(mortality_tracker);
        std::vector<IntegerRaster> empty_vector;
        disperse(
            dispersers,
            established_dispersers,
            susceptible,
            empty_vector,
            infected,
            tmp,  // mortality
            total_populations,
            total_exposed,
            outside_dispersers,
            weather,
            dispersal_kernel,
            suitable_cells,
            establishment_probability,
            generator);
        mortality_tracker = tmp.back();
    }

    /** Move overflowing pest population to other hosts.
     *
     * When the number of pests (pest population) is too high, part of them moves
     * to a different location. Number of infected/infested hosts is considered to be
     * the number of pests (groups of pest) in a raster cell.
     *
     * The movement happens in two stages. First, all the leaving pests are identified
     * and removed from the source cells. Second, the move to the target cells is
     * performed. This means that even if the resulting number of pests in the target
     * cell is considered too high, it is left as is and the move is performed the next
     * time this function is called.
     *
     * If the pests (pest population) cannot be accommodated in the target cell due to
     * the insufficient number of susceptible hosts, the excessive pests die.
     *
     * @param[in,out] susceptible Susceptible hosts
     * @param[in,out] infected Infected hosts
     * @param total_hosts All host individuals in the area. Is equal to
     * infected + exposed + susceptible in the cell.
     * @param[in,out] outside_dispersers Dispersers escaping the rasters
     * @param dispersal_kernel Dispersal kernel to move dispersers (pests)
     * @param[in] suitable_cells List of indices of cells with hosts
     * @param overpopulation_percentage Percentage of occupied hosts when the cell is
     *        considered to be overpopulated
     * @param leaving_percentage Percentage pests leaving an overpopulated cell
     * @param generator Provider of random number generators
     *
     * @note Exposed hosts do not count towards total number of pest,
     *       i.e., *total_host* is assumed to be S + E in SEI model.
     * @note Mortality is not supported by this function, i.e., the mortality rasters
     *       are not modified while the infected are.
     */
    template<typename DispersalKernel>
    void move_overpopulated_pests(
        IntegerRaster& susceptible,
        IntegerRaster& infected,
        const IntegerRaster& total_hosts,
        std::vector<std::tuple<int, int>>& outside_dispersers,
        DispersalKernel& dispersal_kernel,
        std::vector<std::vector<int>>& suitable_cells,
        double overpopulation_percentage,
        double leaving_percentage,
        Generator& generator)
    {
        UNUSED(total_hosts);  // Total hosts is computed now.
        IntegerRaster empty;
        std::vector<IntegerRaster> empty_vector;
        StandardHostPool hosts{
            model_type_,
            susceptible,
            empty_vector,
            0,
            infected,
            empty,
            empty,
            empty_vector,
            empty,
            empty,
            *environment(true),
            false,
            0,
            false,
            0,
            0,
            0,
            suitable_cells};
        StandardPestPool pests{empty, empty, outside_dispersers};
        MoveOverpopulatedPests<
            StandardHostPool,
            StandardPestPool,
            IntegerRaster,
            FloatRaster,
            RasterIndex,
            DispersalKernel>
            move_pest{
                dispersal_kernel,
                overpopulation_percentage,
                leaving_percentage,
                rows_,
                cols_};
        move_pest.action(hosts, pests, generator);
        this->environment(true)->remove_hosts();
    }

    /** Disperse, expose, and infect based on dispersers
     *
     * This function wraps disperse() and infect_exposed() for use in SI
     * and SEI models.
     *
     * In case of SEI model, before calling this function, last item in
     * the exposed vector needs to be ready to be used for exposure,
     * i.e., typically, it should be empty in the sense that there are
     * no hosts in the raster. This is normally taken care of by a
     * previous call to this function. The initial state of the exposed
     * vector should be such that size is latency period in steps plus 1
     * and each raster is empty, i.e., does not contain any hosts
     * (all values set to zero).
     *
     * See the infect_exposed() function for the details about exposed
     * vector, its size, and its items.
     *
     * See disperse() and infect_exposed() for a detailed list of
     * parameters and behavior. The disperse() parameter documentation
     * can be applied as is except that disperse() function's parameter
     * *exposed_or_infested* is expected to change based on the context
     * while this function's parameter *infected* is always the infected
     * individuals. Besides parameters from disperse(), this function
     * has parameter *exposed* which is the same as the one from the
     * infect_exposed() function.
     */
    template<typename DispersalKernel>
    void disperse_and_infect(
        unsigned step,
        IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& infected,
        std::vector<IntegerRaster>& mortality_tracker,
        const IntegerRaster& total_populations,
        IntegerRaster& total_exposed,
        std::vector<std::tuple<int, int>>& outside_dispersers,
        bool weather,
        DispersalKernel& dispersal_kernel,
        std::vector<std::vector<int>>& suitable_cells,
        double establishment_probability,
        Generator& generator)
    {
        this->disperse(
            dispersers,
            established_dispersers,
            susceptible,
            exposed,
            infected,
            mortality_tracker,
            total_populations,
            total_exposed,
            outside_dispersers,
            weather,
            dispersal_kernel,
            suitable_cells,
            establishment_probability,
            generator);
        IntegerRaster empty;
        StandardHostPool host_pool{
            model_type_,
            susceptible,
            exposed,
            latency_period_,
            infected,
            total_exposed,
            empty,
            mortality_tracker,
            empty,
            empty,
            *environment(!weather),
            false,
            0,
            establishment_stochasticity_,
            establishment_probability,
            rows_,
            cols_,
            suitable_cells};
        host_pool.step_forward(step);
        this->environment(!weather)->remove_hosts();
    }

    template<typename DispersalKernel>
    void disperse_and_infect(
        unsigned step,
        IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& infected,
        IntegerRaster& mortality_tracker,
        const IntegerRaster& total_populations,
        IntegerRaster& total_exposed,
        std::vector<std::tuple<int, int>>& outside_dispersers,
        bool weather,
        DispersalKernel& dispersal_kernel,
        std::vector<std::vector<int>>& suitable_cells,
        double establishment_probability,
        Generator& generator)
    {
        std::vector<IntegerRaster> tmp;
        tmp.push_back(mortality_tracker);
        disperse_and_infect(
            step,
            dispersers,
            established_dispersers,
            susceptible,
            exposed,
            infected,
            tmp,
            total_populations,
            total_exposed,
            outside_dispersers,
            weather,
            dispersal_kernel,
            suitable_cells,
            establishment_probability,
            generator);
        mortality_tracker = tmp.back();
    }

    /**
     * @brief Activate storage of dispersers in soil
     *
     * Calling this function activates the soils. By default, the soil pool is not used.
     * The parameters are soil pool used to store the dispersers and
     * a percentage (0-1 ratio) of dispersers which will be send to the soil (and may
     * establish or not depending on the soil pool object).
     *
     * Soil pool is optional and implemented in more general (but experimental) way.
     * This function needs to be called separately some time after the object is created
     * to active the soil part of the simulation. This avoids the need for many
     * constructors or many optional parameters which need default values.
     *
     * @param soil_pool Soils pool object to use for storage
     * @param dispersers_percentage Percentage of dispersers moving to the soil
     */
    void activate_soils(
        std::shared_ptr<SoilPool<IntegerRaster, FloatRaster, RasterIndex, Generator>>
            soil_pool,
        double dispersers_percentage)
    {
        this->soil_pool_ = soil_pool;
        this->to_soil_percentage_ = dispersers_percentage;
    }
};

}  // namespace pops

#endif  // POPS_SIMULATION_HPP
