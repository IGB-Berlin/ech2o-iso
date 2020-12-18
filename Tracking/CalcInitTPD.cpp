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
 * CalcInitTPD.cpp
 * --> derives initial tightly-bound and mobile water maps of tracers, from input soil maps
 *
 *  Created on: Apr 25, 2018
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"
#include "Tracking.h"

int Tracking::CalcInitTPD(Basin &bsn, Control &ctrl){

  UINT4 r, c;
  UINT4 length = bsn.getSortedGrid().cells.size();

#pragma omp parallel default(none)		\
  private(  r,c), shared(length, bsn, ctrl)
  for (UINT4 j = 0; j < length ; j++)
    {
      r = bsn.getSortedGrid().cells[j].row;
      c = bsn.getSortedGrid().cells[j].col;
      
      // d2H
      if(ctrl.sw_2H){
	_d2H_TB1->matrix[r][c] = _d2Hsoil1->matrix[r][c];
	_d2H_MW1->matrix[r][c] = _d2Hsoil1->matrix[r][c];
	_d2H_TB2->matrix[r][c] = _d2Hsoil2->matrix[r][c];
	_d2H_MW2->matrix[r][c] = _d2Hsoil2->matrix[r][c];
      }
      // d18O
      if(ctrl.sw_18O){
	_d18O_TB1->matrix[r][c] = _d18Osoil1->matrix[r][c];
	_d18O_MW1->matrix[r][c] = _d18Osoil1->matrix[r][c];
	_d18O_TB2->matrix[r][c] = _d18Osoil2->matrix[r][c];
	_d18O_MW2->matrix[r][c] = _d18Osoil2->matrix[r][c];
	}
      // Age
      if(ctrl.sw_Age){
	_Age_TB1->matrix[r][c] = _Agesoil1->matrix[r][c];
	_Age_MW1->matrix[r][c] = _Agesoil1->matrix[r][c];
	_Age_TB2->matrix[r][c] = _Agesoil2->matrix[r][c];
	_Age_MW2->matrix[r][c] = _Agesoil2->matrix[r][c];
      }      
    }
  
  return EXIT_SUCCESS;
}
