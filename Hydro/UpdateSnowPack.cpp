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
 * UpdateSnowPack.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: marco.maneta
 */

#include"Basin.h"

int Basin::UpdateSnowPack(Atmosphere &atm, Control &ctrl){

	int r, c;
	double dt = ctrl.dt; //secs

	for (unsigned int j = 0; j < _vSortedGrid.cells.size() ; j++)	{
		r = _vSortedGrid.cells[j].row;
		c = _vSortedGrid.cells[j].col;

		if(atm.getTemperature()->matrix[r][c] < 0)
			_snow->matrix[r][c] += atm.getPrecipitation()->matrix[r][c] * dt;
		else
			_incident_water_depth->matrix[r][c] += atm.getPrecipitation()->matrix[r][c] * dt;


	}
	return EXIT_SUCCESS;
}
