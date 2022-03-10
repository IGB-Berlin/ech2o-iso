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
 *    Marco Maneta
 *******************************************************************************/
/*
 * Intercep_Rutter.cpp
 *
 *  Created on: Oct 12, 2009
 *      Author: Marco Maneta
 */

#include "Forest.h"

int Forest::CanopyInterception(Atmosphere &atm, Control &ctrl, REAL8 &DelCanStor, REAL8 &D, UINT4 s, UINT4 r, UINT4 c)
{

    REAL8 S=0, C=0, Pp=0;
    REAL8 dC=0;
    REAL8 dt = ctrl.dt;
    REAL8 p = 0; //throughfall coefficient

    REAL8 ratio = 0;
    REAL8 dC_net=0;

    Pp = atm.getPrecipitation()->matrix[r][c];

    UINT4 nsp = getNumSpecies();

    if(s == nsp -1) //for bare soil, water reaching the ground is pp times its proportion of the cell
      D = Pp;
    else{
      S = getMaxCanopyStorage(s, r, c); //returns max canopy storage as calculated with current LAI
      C = this->getIntercWater(s, r, c); //holds current water storage in the canopy

      if(fabs(C) < RNDOFFERR){
	_species[s]._WaterStorage->matrix[r][c] = 0.0;
	C = 0.0;
      }
      // For checking the fluxes to tracking
      dC = Pp * dt;                       //change in canopy in m
      p = _species[s].Throughfall_coeff;

      if(ctrl.sw_intercept==1){
	ratio = C / S;                    // ratio of canopy storage to maximum canopy storage [ref]
	dC_net = dC * (1 - p);            // rainfall that is potentially intercepted by the canopy
	D = dC_net * ratio + dC * p;      // proportion of direct throughfall + exponential decrease in canopy storage with fullness
	D += ( C+(dC_net*(1-ratio)) <= S ? 0.0 : (C+ (dC_net*(1-ratio)) - S) );
	
      }	else {
	
	D = dC * p;                       // proportion of precipitation that is direct throughfall 
	D += ( C+(dC*(1-p)) <= S ? 0.0 : (C+ (dC*(1-p)) - S) );   //estimate additional throughfall from overflow of canopy storage
	//D = min<REAL8>(D, C+dC);

      }
      
      DelCanStor = Pp * dt - D;           //Change in canopy storage
      
      D /= dt; //make D into flux ms-1

      if(fabs(DelCanStor) < RNDOFFERR)
	DelCanStor = 0;

      _species[s]._WaterStorage->matrix[r][c] += DelCanStor;
     
    }
    return EXIT_SUCCESS;
}
