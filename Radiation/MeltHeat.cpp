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
 * MeltHeat.cpp
 *
 *  Created on: Oct 1, 2010
 *      Author: Marco.Maneta
 */

#include"Basin.h"

double Basin::MeltHeat(Atmosphere &atm, Control &ctrl, const double &Ts, const double &swe, const double &M, int row, int col){

	double L1 = 0;
	double L2 = 0;

	if (Ts < 0)
		return 0.0;

	L1 = rho_w * lat_heat_fus * swe / ctrl.dt;
	L2 = rho_w * lat_heat_fus * M * Ts;

	if (L1 < L2)
		return -L1;
	else
		return -L2;

}
