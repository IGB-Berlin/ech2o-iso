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
 * NetRadCanopy.cpp
 *
 *  Created on: Jul 9, 2010
 *      Author: Marco.Maneta
 */

#include"Forest.h"

double Forest::NetRadCanopy(Atmosphere &atm, const double &Ts, REAL8 emiss, REAL8 albedo, REAL8 Kbeers, REAL8 lai, int row, int col){

	return	atm.getIncomingShortWave()->matrix[row][col] * (1 - albedo) * ( 1 - expl(-Kbeers * lai) ) +
				emiss*atm.getIncomingLongWave()->matrix[row][col] -
				emiss * stefboltz * powl(Ts + 273.2, 4);

}
