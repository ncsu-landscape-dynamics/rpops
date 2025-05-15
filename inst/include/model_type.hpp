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

#ifndef POPS_MODEL_TYPE_HPP
#define POPS_MODEL_TYPE_HPP

namespace pops {

/** The type of a epidemiological model (SI or SEI)
 */
enum class ModelType
{
    SusceptibleInfected,  ///< SI (susceptible - infected)
    SusceptibleExposedInfected  ///< SEI (susceptible - exposed - infected)
};

/*! Get a corresponding enum value for a string which is a model type name.
 *
 * Throws an std::invalid_argument exception if the value was not
 * found or is not supported (which is the same thing).
 */
inline ModelType model_type_from_string(const std::string& text)
{
    if (text == "SI" || text == "SusceptibleInfected" || text == "susceptible-infected"
        || text == "susceptible_infected")
        return ModelType::SusceptibleInfected;
    else if (
        text == "SEI" || text == "SusceptibleExposedInfected"
        || text == "susceptible-exposed-infected"
        || text == "susceptible_exposed_infected")
        return ModelType::SusceptibleExposedInfected;
    else
        throw std::invalid_argument(
            "model_type_from_string: Invalid"
            " value '"
            + text + "' provided");
}

/*! Overload which allows to pass C-style string which is nullptr (NULL)
 *
 * @see model_type_from_string(const std::string& text)
 */
inline ModelType model_type_from_string(const char* text)
{
    // call the string version
    return model_type_from_string(text ? std::string(text) : std::string());
}

}  // namespace pops

#endif  // POPS_MODEL_TYPE_HPP
