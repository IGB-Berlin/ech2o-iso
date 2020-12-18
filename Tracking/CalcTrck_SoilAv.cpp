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
 * CalcTrck_SoilAv.cpp
 *
 *  Created on: Jan 31, 2018
 *      Author: Sylvain Kuppel
 */

#include "Tracking.h"

// Calculates isotopes / age weighted average over the soil layers

int Tracking::Calcd2Hsoil_Av(Basin &bsn){
  
  double depth, d1, d2, d3;
  //double fc; 
  double theta1, theta2, theta3;
  int r, c;
#pragma omp parallel default(shared) private(r,c,depth, d1, d2, d3, \
					     theta1, theta2, theta3)
  {
#pragma omp for nowait
    for (unsigned int j = 0; j < bsn.getSortedGrid().cells.size(); j++) {
      r = bsn.getSortedGrid().cells[j].row;
      c = bsn.getSortedGrid().cells[j].col;
      depth = bsn.getSoilDepth()->matrix[r][c];
      d1 = bsn.getSoilDepth1()->matrix[r][c];
      d2 = bsn.getSoilDepth2()->matrix[r][c];
      d3 = depth - d1 - d2;
      //fc = bsn.getFieldCapacity()->matrix[r][c];
      theta1 = bsn.getSoilMoist1()->matrix[r][c];
      theta2 = bsn.getSoilMoist2()->matrix[r][c];
      theta3 = bsn.getSoilMoist3()->matrix[r][c];
      
      _d2HsoilAv->matrix[r][c] = (_d2Hsoil1->matrix[r][c] * d1 * theta1
				 + _d2Hsoil2->matrix[r][c] * d2 * theta2
				 + _d2Hsoil3->matrix[r][c] * d3 * theta3 )/ 
	(d1*theta1+d2*theta2+d3*theta3);
      //std::min<double>(fc,theta3)
      //+ _d2Hgroundwater->matrix[r][c] * d3 * std::max<double>(0.0,theta3-fc))
    }
  }
  return EXIT_SUCCESS;
}

int Tracking::Calcd18Osoil_Av(Basin &bsn){
  
  double depth, d1, d2, d3;
  //double fc; 
  double theta1, theta2, theta3;
  int r, c;
#pragma omp parallel default(shared) private(r,c,depth, d1, d2, d3, \
					     theta1, theta2, theta3)
  {
#pragma omp for nowait
    for (unsigned int j = 0; j < bsn.getSortedGrid().cells.size(); j++) {
      r = bsn.getSortedGrid().cells[j].row;
      c = bsn.getSortedGrid().cells[j].col;
      depth = bsn.getSoilDepth()->matrix[r][c];
      d1 = bsn.getSoilDepth1()->matrix[r][c];
      d2 = bsn.getSoilDepth2()->matrix[r][c];
      d3 = depth - d1 - d2;
      //fc = bsn.getFieldCapacity()->matrix[r][c];
      theta1 = bsn.getSoilMoist1()->matrix[r][c];
      theta2 = bsn.getSoilMoist2()->matrix[r][c];
      theta3 = bsn.getSoilMoist3()->matrix[r][c];
      
      _d18OsoilAv->matrix[r][c] = (_d18Osoil1->matrix[r][c] * d1 * theta1
				 + _d18Osoil2->matrix[r][c] * d2 * theta2
				 + _d18Osoil3->matrix[r][c] * d3 * theta3 )/ 
	(d1*theta1+d2*theta2+d3*theta3);
      //std::min<double>(fc,theta3)
      //+ _d18Ogroundwater->matrix[r][c] * d3 * std::max<double>(0.0,theta3-fc)) 

    }
  }
  return EXIT_SUCCESS;
}

int Tracking::CalcAgesoil_Av(Basin &bsn){
  
  double depth, d1, d2, d3;
  //double fc; 
  double theta1, theta2, theta3;
  int r, c;
#pragma omp parallel default(shared) private(r,c,depth, d1, d2, d3, \
					     theta1, theta2, theta3)
  {
#pragma omp for nowait
    for (unsigned int j = 0; j < bsn.getSortedGrid().cells.size(); j++) {
      r = bsn.getSortedGrid().cells[j].row;
      c = bsn.getSortedGrid().cells[j].col;
      depth = bsn.getSoilDepth()->matrix[r][c];
      d1 = bsn.getSoilDepth1()->matrix[r][c];
      d2 = bsn.getSoilDepth2()->matrix[r][c];
      d3 = depth - d1 - d2;
      //fc = bsn.getFieldCapacity()->matrix[r][c];
      theta1 = bsn.getSoilMoist1()->matrix[r][c];
      theta2 = bsn.getSoilMoist2()->matrix[r][c];
      theta3 = bsn.getSoilMoist3()->matrix[r][c];
      
      _AgesoilAv->matrix[r][c] = (_Agesoil1->matrix[r][c] * d1 * theta1
				 + _Agesoil2->matrix[r][c] * d2 * theta2
				  + _Agesoil3->matrix[r][c] * d3 * theta3 )/ 
	(d1*theta1+d2*theta2+d3*theta3);
      //std::min<double>(fc,theta3)
      //+ _Agegroundwater->matrix[r][c] * d3 * std::max<double>(0.0,theta3-fc)) 

    }
  }
  return EXIT_SUCCESS;
}
