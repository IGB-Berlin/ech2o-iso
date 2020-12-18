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
 * CalcTPDMoisture.cpp
 * --> translation of tightly-bound domain tension in mobile water moisture
 *
 *  Created on: Apr 6, 2018
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"

int Basin::CalcTPDMoisture(Control &ctrl){

  UINT4 r, c;
  
  UINT4 length = _vSortedGrid.cells.size();
  
#pragma omp parallel default(none)		\
  private(  r,c), shared(length, ctrl)
  for (UINT4 j = 0; j < length ; j++)
    {
      r = _vSortedGrid.cells[j].row;
      c = _vSortedGrid.cells[j].col;
      
      // Lower limit is field capacity for now (more often much higher).
      // It's necessary since GW is here part of MW only.
      // It could be improved / made more flexible
      _moist_MW1->matrix[r][c] = 
	std::min<double>(_fieldcapL1->matrix[r][c],
			 powl(_psi_ae->matrix[r][c] / _psi_MW->matrix[r][c], 
			      1/_BClambda->matrix[r][c])
			 * (_porosityL1->matrix[r][c] - _theta_rL1->matrix[r][c]) + _theta_rL1->matrix[r][c]);

      _moist_MW2->matrix[r][c] = 
	std::min<double>(_fieldcapL2->matrix[r][c],
			 powl(_psi_ae->matrix[r][c] / _psi_MW->matrix[r][c], 
			      1/_BClambda->matrix[r][c])
			 * (_porosityL2->matrix[r][c] - _theta_rL2->matrix[r][c]) + _theta_rL2->matrix[r][c]);
      
    }
  
  return EXIT_SUCCESS;
}
