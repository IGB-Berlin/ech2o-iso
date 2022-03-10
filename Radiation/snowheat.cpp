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
 * snwwheat.cpp
 *
 *  Created on: Nov 20, 2009
 *      Author: marco.maneta
 */

#include"Basin.h"

double Basin::SnowHeat(Atmosphere &atm, Control &ctrl, const double &Ts, int row, int col){

	double Ts_old = 0; // temperature of surface at time step t-1
	double h = 0; // depth of snow water equivalent
	double R = 0; // energy advected by rain Wm-2
	float dt = ctrl.dt;

	Ts_old = _Temp_s_old->matrix[row][col];
	h = _snow->matrix[row][col];

	return
		-spec_heat_ice * rho_w * h * (Ts - Ts_old) * (1 / dt) - R; // TODO: solve the energy advected by rain
									   //TODO: change the energy differential to include temperature of snowpack

}
