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
 * CalculateSatArea.cpp
 *
 *  Created on: Nov 7, 2016
 *      Author: Sylvain Kuppel
 */

#include"Basin.h"

int Basin::CalculateSatArea(Atmosphere &atm, Control &ctrl) {

  UINT4 r, c;
  REAL8 poros; //porosity

  REAL8 IsSaturated = 0; //
  REAL8 theta1 = 0; // water content first layer

#pragma omp parallel default(none)\
  private(r,c, theta1, poros, IsSaturated) \
  shared(atm, ctrl)
  {
    #pragma omp for nowait
   
    for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
      
      r = _vSortedGrid.cells[j].row;
      c = _vSortedGrid.cells[j].col;
      IsSaturated = 0;
      
      //surface routing stuff
      theta1 = _soilmoist1->matrix[r][c];
      poros = _porosityL1->matrix[r][c];
      
      if (fabs(poros - theta1) < RNDOFFERR)
	IsSaturated = 1;
      
      _IsSaturated->matrix[r][c] = IsSaturated;
      
    }
  }   
  return EXIT_SUCCESS;
}
