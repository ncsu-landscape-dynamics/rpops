/*
  * PoPS model - mortality
  *
  * Copyright (C) 2015-2020 by the authors.
  *
  * Authors: Chris Jones >cjones1688 gmail com
  *          Anna Petrasova <akratoc gmail com>
  *          Vaclav Petras <wenzeslaus gmail com>
  *
  * The code contained herein is licensed under the GNU General Public
  * License. You may obtain a copy of the GNU General Public License
  * Version 2 or later at the following locations:
  *
  * http://www.opensource.org/licenses/gpl-license.html
  * http://www.gnu.org/copyleft/gpl.html
*/
  
#ifndef POPS_MORTALITY_HPP
#define POPS_MORTALITY_HPP
  
#include "raster.hpp"
#include "date.hpp"
#include "scheduling.hpp"
  
#include <map>
#include <vector>
#include <string>
#include <functional>

namespace pops {
  

void mortality(
  IntegerRaster& infected,
  double mortality_rate,
  int current_year,
  int first_mortality_year,
  IntegerRaster& mortality,
  std::vector<IntegerRaster>& mortality_tracker_vector)
{
  if (current_step >= (first_mortality_year)) {
    
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        for (int index = 0; index <= max_index;
             index++) {
          int mortality_in_index = 0;
          if (mortality_tracker_vector[index](i, j) > 0) {
            mortality_in_index =
              mortality_rate
            * mortality_tracker_vector[index](i, j);
            mortality_tracker_vector[index](i, j) -=
              mortality_in_index;
            mortality(i, j) += mortality_in_index;
            mortality_current_year += mortality_in_index;
            if (infected(i, j) > 0) {
              infected(i, j) -= mortality_in_index;
            }
          }
        }
      }
    }
  }
}

if (model_type_ == ModelType::SusceptibleExposedInfected) {
  if (step >= latency_period_) {
    // Oldest item needs to be in the front
    auto& oldest = exposed.front();
    // Move hosts to infected raster
    infected += oldest;
    mortality_tracker += oldest;
    // Reset the raster
    // (hosts moved from the raster)
    oldest.fill(0);
  }
  // Age the items and the used one to the back
  // elements go one position to the left
  // new oldest goes to the front
  // old oldest goes to the back
  rotate_left_by_one(exposed);
}