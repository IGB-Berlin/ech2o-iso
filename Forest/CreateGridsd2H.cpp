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
 * CreateGrovesd2H.cpp
 *
 *  Created on: Mar 17, 2017
 *      Author: Sylvain Kuppel
 */

#include "Grove.h"

int Grove::CreateGridsd2H(grid *base){

	try{
		_d2Hcanopy = new grid (*base);
		_d2Hthroughfall = new grid (*base);
		_d2HevapI = new grid (*base);
		_d2HevapT = new grid (*base);
		_d2HevapS = new grid (*base);

		_d2HevapI_Vap = new grid (*base);
		_d2HevapT_Vap = new grid (*base);

	}catch(const exception& e){

		cerr << "Failed allocate memory for Grove grid (2H) object \n" << e.what() << endl;
		exit(EXIT_FAILURE);
	}

	return EXIT_SUCCESS;
}
