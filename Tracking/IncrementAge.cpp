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
 * IncrementAge.cpp
 *
 *  Created on: Nov 15, 2016
 *      Author: Sylvain Kuppel
 */

#include "Tracking.h"

int Tracking::IncrementAge(Basin &bsn, Control &ctrl){

  UINT4 r, c, s;
  REAL8 dt = ctrl.dt / 86400 ; // units: days
  UINT4 nsp = bsn.getNumSpecies();

#pragma omp parallel default(shared) private(r,c, s)
  {

#pragma omp for

    for (unsigned int j = 0; j < bsn.getSortedGrid().cells.size(); j++) {

      r = bsn.getSortedGrid().cells[j].row;
      c = bsn.getSortedGrid().cells[j].col;
      
      _Agecanopy_sum->matrix[r][c] += dt; // Canopy interception storage
      for (s = 0; s < nsp - 1; s++)
	bsn.setAgecanopy(s, r, c, bsn.getAgecanopy(s)->matrix[r][c] + dt);
      
      _Agesnowpack->matrix[r][c] = bsn.getSnowWaterEquiv()->matrix[r][c] > RNDOFFERR ?
	                                 _Agesnowpack->matrix[r][c] + dt : 0.0 ; 	// Snowpack
      _Agesurface->matrix[r][c] = bsn.getPondingWater()->matrix[r][c] > RNDOFFERR ?
	                                 _Agesurface->matrix[r][c] + dt : 0.0 ; 	// Ponding
      _Agegroundwater->matrix[r][c] = bsn.getGrndWater()->matrix[r][c] > RNDOFFERR ?
	                                 _Agegroundwater->matrix[r][c] + dt : 0.0; 	// Groundwater
      if(ctrl.sw_deepGW)
	_AgedeepGW->matrix[r][c] = bsn.getDeepGW()->matrix[r][c] > RNDOFFERR ? 
					_AgedeepGW->matrix[r][c] + dt : 0.0; 		//Extra GW

      if(ctrl.sw_TPD){
	_Age_TB1->matrix[r][c] += dt; // Mobile water layer 1
	_Age_MW1->matrix[r][c] += dt; // Mobile water layer 1
	_Age_TB2->matrix[r][c] += dt; // Tightly-bound layer 2
	_Age_MW2->matrix[r][c] += dt; // Mobile water layer 2
      }
      _Agesoil1->matrix[r][c] += dt; // Vadose layer 1
      _Agesoil2->matrix[r][c] += dt; // Vadose layer 2
      _Agesoil3->matrix[r][c] += dt; // Vadose layer 3

      // The outputs as well, for mass balance consistency (notably)
      _AgeevapS_sum->matrix[r][c] += dt;
      _AgeevapI_sum->matrix[r][c] += dt;
      _AgeevapT_sum->matrix[r][c] += dt;
      _Ageleakage->matrix[r][c] += dt;
      for (s = 0; s < nsp - 1; s++){
	bsn.setAgeevapS(s, r, c, bsn.getAgeevapS(s)->matrix[r][c] + dt);
	bsn.setAgeevapI(s, r, c, bsn.getAgeevapI(s)->matrix[r][c] + dt);
	bsn.setAgeevapT(s, r, c, bsn.getAgeevapT(s)->matrix[r][c] + dt);
      }
    }
    
#pragma omp for
    for (UINT4 i = 0; i < _AgeOvlndOutput.cells.size(); i++){ //b->getSortedGrid().cells.size();; i++){
      _AgeOvlndOutput.cells[i].val += dt ;
      _AgeGwtrOutput.cells[i].val += dt ;
      if(ctrl.sw_deepGW)
	_AgeDeepGwtrOutput.cells[i].val += dt;
      }
  }
  
  return EXIT_SUCCESS;
}
