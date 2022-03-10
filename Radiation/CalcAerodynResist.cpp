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
 * CalcAerodynResist.cpp
 *
 *  Created on: Nov 19, 2009
 *      Author: marco.maneta
 */

#include"Basin.h"

REAL8 Basin::CalcAerodynResist(REAL8 u_za, REAL8 z_a, REAL8 z_0u, REAL8 z_du, REAL8 z_0o, REAL8 z_do, REAL8 Ht, REAL8 LAI,
		REAL8 Ts, REAL8 Ta, INT4 option, bool surface) {

	REAL8 rlog_u = 0;
	REAL8 rlog_o = 0;
	REAL8 rexp = 0;
	REAL8 vt = 0;
	REAL8 a = 0;
	REAL8 lm = 2;
	REAL8 zt = (z_do + z_0o) * 0.8; // for the overstory, exponential resistance is for only half the canopy

	if(surface == true){
		rlog_u = rlog(u_za, z_a, z_du, z_0u, option);
		zt = z_du + z_0u;
	}

	if(LAI > 0){ //if there is overstory, calculate the resistance associated with the exponential profile

		lm = sqrt(4*0.05*Ht/(PI*LAI)); //average leaf separation assuming an average leave with of 5 cm (0.05 m)
		vt = u_za * vonkarman * vonkarman * (Ht - z_do) / log( (z_a - z_do)/z_0o ); // conductance at the top of the tree, which is the bottom of the overstory wind log profile
		a = min<double>(max<double>(1.0,0.2 * LAI * Ht / lm), 4); //wind extinction factor as per foken 2008, Campbell and Norman 1998

		rexp = (1 / (a * vt) ) * Ht * expl(a) * ( expl(-a * ( zt ) / Ht) - expl(-a * (z_do + z_0o) / Ht));
		rlog_o = rlog(u_za, z_a, z_do, z_0o, option); // add the resistance due to the logarithmic profile over the

		//rlog_u=0; //for the time being the exponential profile goes all the way down to 0 velocity (zt = z_du+z_0u) for surface

	}
	//otherwise resistance is simply the aerodynamic resistance associated with a logarithmic wind profile

	return (rlog_u + rexp + rlog_o);


}
