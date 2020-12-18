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
 * NetRad_water.cpp
 *
 *  Created on: Nov 18, 2009
 *      Author: marco.maneta
 */

#include"Basin.h"

double Basin::NetRad_water(Atmosphere &atm, const double &Tw, int row, int col){

		double es = 0.97; //emisivity of surface
		double ea = AirEmissivity(atm.getTemperature()->matrix[row][col]);
		double albedo = 0.06;

		return	atm.getIncomingShortWave()->matrix[row][col] * (1 - albedo)
		  + es * (1 - ea) * atm.getIncomingLongWave()->matrix[row][col]
		  + es * ea * stefboltz * powl(atm.getTemperature()->matrix[row][col] + 273.2 , 4)
   		  - es * stefboltz * powl(Tw + 273.2, 4);
		
}
