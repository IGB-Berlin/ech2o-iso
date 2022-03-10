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
 * SnowOutputPhase.cpp
 *
 *  Created on: Nov 25, 2009
 *      Author: marco.maneta
 */

#include"Basin.h"

double Basin::SnowOutput(Atmosphere &atm, Control &ctrl, Tracking &trck,
			 const double &meltheat, int row, int col) {

	double h = 0; // depth of snow water equivalent
	double h0 = 0; // depth of snow water equivalent previous time-step
	double dh = 0; // depth of snow output - decrease in snow water equivalent depth

	h = _snow->matrix[row][col] < RNDOFFERR ? 0.0 : _snow->matrix[row][col];

	h0 = _snow_old->matrix[row][col] < RNDOFFERR ? 0.0 : _snow_old->matrix[row][col];

	_snow->matrix[row][col] = h;
	
	//if there is snow pack and latent heat of melt < 0
	if (h > RNDOFFERR && meltheat < RNDOFFERR){
	  dh = -meltheat * ctrl.dt / (lat_heat_fus*rho_w); 	//transform LE of melt into snowmelt depth
	  if (dh > h)						//if this water energy equivalent > SWE
	    dh = h; 						// output = water in snowpack
	  _snow->matrix[row][col] -= dh;

	} else
	  h = 0.0;
    
	// Flux tracking after snowmelt
	if(ctrl.sw_trck){
	  _FluxSnowtoSrf->matrix[row][col] = dh;
	  trck.MixingV_snow(atm, *this, ctrl, h, h0, dh, row, col);
	}

	_snow_old->matrix[row][col] = _snow->matrix[row][col]; // update the snow for the next time-step

	return dh;
	
}
