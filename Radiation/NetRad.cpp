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
 * NetRad.cpp
 *
 *  Created on: Nov 18, 2009
 *      Author: marco.maneta
 */

#include"Basin.h"

double Basin::NetRad(Atmosphere &atm, const double &Ts, REAL8 Kbeers, REAL8 lai, REAL8 ec, REAL8 Tc, INT4 iswater, int row, int col){

    	double es, albedo;

	if(iswater==1){
	  es = 0.97;
	  albedo = 0.06;
	} else {
	  es = _emiss_surf->matrix[row][col];
	  //weights albedo with soil and concrete using impervious
	  albedo = _albedo->matrix[row][col] * (1-_fImperv->matrix[row][col]) + conc_albedo * _fImperv->matrix[row][col];
	  albedo = _snow->matrix[row][col] > RNDOFFERR ? max_snow_albedo : _albedo->matrix[row][col]; //TODO: include albedo decay with time and with covered area (covered area a function of snowdepth?)
	}

	return	atm.getIncomingShortWave()->matrix[row][col] * (1 - albedo) * ( expl(-1*Kbeers * lai) )
 			+ es * (1 - ec) * atm.getIncomingLongWave()->matrix[row][col]
			+ es * ec * stefboltz * powl(Tc + 273.2, 4)
			- es * stefboltz * powl(Ts + 273.2, 4);

}
