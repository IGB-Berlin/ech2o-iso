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
 * GrowRoots.cpp
 *
 *  Created on: Jul 1, 2010
 *      Author: Marco.Maneta
 */

#include "Forest.h"

int Forest::GrowRoots(UINT4 spec, UINT4 row, UINT4 col, REAL8 dt){

	_species[spec]._RootMass->matrix[row][col] += max<REAL8>(-0.95*_species[spec]._RootMass->matrix[row][col], 
							_species[spec]._Del_RootMass->matrix[row][col] - 
							_species[spec]._RootMass->matrix[row][col] * _species[spec].RootTurnover * dt);

  	return EXIT_SUCCESS;

}
