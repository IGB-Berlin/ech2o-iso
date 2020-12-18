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
 * checkForestDatabase.cpp
 *
 *  Created on: Mar 29, 2012
 *      Author: marco
 */

#include "Forest.h"

void Forest::checkForestDatabase()
{
	UINT4 i = 0;
	try{
	for(i = 0; i < _Nsp-1; i++){ // substract 1 because _Nsp includes bare soil

		/*
		 * PARAMETER FILE CHECKS
		 */

		//Checks if tmax is larger than tmin and that topt is in between
		if(_species[i].TempMax < _species[i].TempMin)
			throw("Parameter TMin is larger than TMax for species ");
		if( (_species[i].TempOpt > _species[i].TempMax) || (_species[i].TempOpt < _species[i].TempMin) )
			throw("Optimal temperature parameter is out of the Tmin - Tmax range given for species ");


		/*
		 * STATE VARIABLE CHECKS
		 *
		 */

	}
	}catch(std::string &e){

		cerr << "error: " << e << i;
		throw;

	}


}
