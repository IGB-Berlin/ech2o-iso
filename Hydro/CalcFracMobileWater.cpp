/*******************************************************************************
 * Ech2o, a spatially-distributed, ecohydrologic simulator
 * Copyright (c) 2016 Marco Maneta <marco.maneta@umontana.edu>
 *
 *     This file is part of ech2o, a hydrologic model developed at the 
 *     University of Montana.
 *
 *     Ech2o is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     Ech2o is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with Ech2o.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contributors:
 *    Marco Maneta, Sylvain Kuppel
 *******************************************************************************/
/*
 * CalcFracMobileWater.cpp
 *
 *  Created on: May 3, 2018
 *      Author: Sylvain Kuppel
 */

#include"Basin.h"

int Basin::CalcFracMobileWater() {

  UINT4 r, c;
  REAL8 theta_MW1 = 0; // tightly-bound to mobile transition moisture
  REAL8 theta_MW2 = 0; // tightly-bound to mobile transition moisture
  
  REAL8 theta1 = 0; // water content first layer
  REAL8 theta2 = 0; // water content 2nd layer

#pragma omp parallel default(shared)\
  private(r,c, theta1, theta2, theta_MW1, theta_MW2)
  {
#pragma omp for nowait
    
    for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
      
      r = _vSortedGrid.cells[j].row;
      c = _vSortedGrid.cells[j].col;
      
      theta1 = _soilmoist1->matrix[r][c];
      theta2 = _soilmoist2->matrix[r][c];

      theta_MW1 = _moist_MW1->matrix[r][c];
      theta_MW2 = _moist_MW2->matrix[r][c];
            
      _fracMW1->matrix[r][c] = 1 - std::min<double>(1, theta_MW1 / theta1);
      _fracMW2->matrix[r][c] = 1 - std::min<double>(1, theta_MW2 / theta2);
            
    }
  }   
  return EXIT_SUCCESS;
}
