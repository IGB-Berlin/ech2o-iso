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
 * CalcTPDtoLayers.cpp
 * --> derives layer-averaged values from tightly-ound and mobile domain signatures
 *
 *  Created on: Apr 25, 2018
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"
#include "Tracking.h"

int Tracking::CalcTPDtoLayers(Basin &bsn, Control &ctrl){

  UINT4 r, c;
  UINT4 length = bsn.getSortedGrid().cells.size();
  double theta1, theta2;
  double theta_MW1, theta_MW2;
  double d1, d2;

#pragma omp parallel default(none)		\
  private(  r,c, theta1, theta2, theta_MW1, theta_MW2, d1, d2),	\
  shared(length, bsn, ctrl, cout)
  for (UINT4 j = 0; j < length ; j++)
    {
      r = bsn.getSortedGrid().cells[j].row;
      c = bsn.getSortedGrid().cells[j].col;
      theta1 = bsn.getSoilMoist1()->matrix[r][c];
      theta2 = bsn.getSoilMoist2()->matrix[r][c];
      theta_MW1 = bsn.getMoistureMW1()->matrix[r][c];
      theta_MW2 = bsn.getMoistureMW2()->matrix[r][c];
      d1 = bsn.getSoilDepth1()->matrix[r][c];
      d2 = bsn.getSoilDepth2()->matrix[r][c];

      // d2H
      if(ctrl.sw_2H){
	_d2Hsoil1->matrix[r][c] = (std::min<double>(theta_MW1,theta1)*_d2H_TB1->matrix[r][c] +
				   std::max<double>(0,theta1-theta_MW1)*_d2H_MW1->matrix[r][c])
	  / theta1;

	_d2Hsoil2->matrix[r][c] = (std::min<double>(theta_MW2,theta2)*_d2H_TB2->matrix[r][c] +
				   std::max<double>(0,theta1-theta_MW2)*_d2H_MW2->matrix[r][c])
	  / theta2;
      }
      // d18O
      if(ctrl.sw_18O){
	_d18Osoil1->matrix[r][c] = (std::min<double>(theta_MW1,theta1)*_d18O_TB1->matrix[r][c] +
				   std::max<double>(0,theta1-theta_MW1)*_d18O_MW1->matrix[r][c])
	  / theta1;
	_d18Osoil2->matrix[r][c] = (std::min<double>(theta_MW2,theta2)*_d18O_TB2->matrix[r][c] +
				   std::max<double>(0,theta1-theta_MW2)*_d18O_MW2->matrix[r][c])
	  / theta2;
	}
      // Age
      if(ctrl.sw_Age){
	_Agesoil1->matrix[r][c] = (std::min<double>(theta_MW1,theta1)*_Age_TB1->matrix[r][c] +
				   std::max<double>(0,theta1-theta_MW1)*_Age_MW1->matrix[r][c])
	  / theta1;
	_Agesoil2->matrix[r][c] = (std::min<double>(theta_MW2,theta2)*_Age_TB2->matrix[r][c] +
				   std::max<double>(0,theta1-theta_MW2)*_Age_MW2->matrix[r][c])
	  / theta2;

	// Vertical average
	if(ctrl.Rep_Age_TBUp or ctrl.RepTs_Age_TBUp)
	  _Age_TB12->matrix[r][c] = (_Age_TB1->matrix[r][c]*d1*std::min<double>(theta_MW1,theta1) +
	     _Age_TB2->matrix[r][c]*d2*std::min<double>(theta_MW2,theta2)) / 
	    (std::min<double>(theta_MW1,theta1)*d1 + std::min<double>(theta_MW2,theta2)*d2);

	if(ctrl.Rep_Age_MWUp or ctrl.RepTs_Age_MWUp)
	  _Age_MW12->matrix[r][c] = std::max<double>(0,theta1-theta_MW1) + 
	    std::max<double>(0,theta2-theta_MW2) > RNDOFFERR ?
	    (_Age_MW1->matrix[r][c]*d1*std::max<double>(0,theta1-theta_MW1) +
	     _Age_MW2->matrix[r][c]*d2*std::max<double>(0,theta2-theta_MW2)) /
	    (std::max<double>(0,theta1-theta_MW1)*d1 + std::max<double>(0,theta2-theta_MW2)*d2) : 
	    0.0;
      }      
    }
  
  return EXIT_SUCCESS;
}
