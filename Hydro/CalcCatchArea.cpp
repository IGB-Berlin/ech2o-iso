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
 * CalcCatchArea.cpp
 *
 *  Created on: Jul 29, 2010
 *      Author: Marco.Maneta
 */

#include "Basin.h"


int Basin::CalcCatchArea(){

	UINT4 r, c, d, length;
	REAL8 area;

	length = _vSortedGrid.cells.size();

	for (UINT4 j = 0; j < length ; j++) {
		r = _vSortedGrid.cells[j].row;
		c = _vSortedGrid.cells[j].col;
		d = _vSortedGrid.cells[j].dir;

		area = _catcharea->matrix[r][c] + _dx * _dx;

		//add the previously calculated discharge to the downstream cell
 	        switch (d) {
			case 1:   _catcharea->matrix[r+1][c-1] = area; break;
			case 2:   _catcharea->matrix[r+1][c]= area; break;
			case 3:   _catcharea->matrix[r+1][c+1]= area; break;
			case 4:   _catcharea->matrix[r][c-1]= area; break;
			case 5: _catcharea->matrix[r][c] = area; break;//adds the area of the outlet cell
			case 6:   _catcharea->matrix[r][c+1]= area; break;
			case 7:   _catcharea->matrix[r-1][c-1]= area; break;
			case 8:   _catcharea->matrix[r-1][c]= area; break;
			case 9:   _catcharea->matrix[r-1][c+1]+= area; break;
			default: return -1;
		}
	}
	return EXIT_SUCCESS;
}
