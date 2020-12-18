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
 * CalcTrck_L1L2.cpp
 *
 *  Created on: Jun 27, 2017
 *      Author: Sylvain Kuppel
 */

#include "Tracking.h"
#include "Grid.h"

// Calculates isotopes / age weighted average over the top two soil layers
int Tracking::Calcd2Hsoil_12(Basin &bsn){
  
  double d1, d2;
  int r, c;
#pragma omp parallel default(shared) private(d1, d2, r,c)
  {
#pragma omp for nowait
    for (unsigned int j = 0; j < bsn.getSortedGrid().cells.size(); j++) {
      r = bsn.getSortedGrid().cells[j].row;
      c = bsn.getSortedGrid().cells[j].col;
      d1 = bsn.getSoilDepth1()->matrix[r][c];
      d2 = bsn.getSoilDepth2()->matrix[r][c];
      _d2Hsoil_12->matrix[r][c] = (_d2Hsoil1->matrix[r][c] * d1
				+ _d2Hsoil2->matrix[r][c] * d2) / (d1+d2);
    }
  }
  return EXIT_SUCCESS;
}

int Tracking::Calcd18Osoil_12(Basin &bsn){
  
  double d1, d2;
  int r, c;
#pragma omp parallel default(shared) private(d1, d2, r,c)
  {
#pragma omp for nowait
    for (unsigned int j = 0; j < bsn.getSortedGrid().cells.size(); j++) {
      r = bsn.getSortedGrid().cells[j].row;
      c = bsn.getSortedGrid().cells[j].col;
      d1 = bsn.getSoilDepth1()->matrix[r][c];
      d2 = bsn.getSoilDepth2()->matrix[r][c];
      _d18Osoil_12->matrix[r][c] = (_d18Osoil1->matrix[r][c] * d1
				    + _d18Osoil2->matrix[r][c] * d2) / (d1+d2);
    }
  }
  return EXIT_SUCCESS;
  }

int Tracking::CalcAgesoil_12(Basin &bsn){
  
  double d1, d2;
  int r, c;
#pragma omp parallel default(shared) private(d1, d2, r,c)
  {
#pragma omp for nowait
    for (unsigned int j = 0; j < bsn.getSortedGrid().cells.size(); j++) {
      r = bsn.getSortedGrid().cells[j].row;
      c = bsn.getSortedGrid().cells[j].col;
      d1 = bsn.getSoilDepth1()->matrix[r][c];
      d2 = bsn.getSoilDepth2()->matrix[r][c];
      _Agesoil_12->matrix[r][c] = (_Agesoil1->matrix[r][c] * d1
				   + _Agesoil2->matrix[r][c] * d2) / (d1+d2);
    }
  }
  return EXIT_SUCCESS;
}
