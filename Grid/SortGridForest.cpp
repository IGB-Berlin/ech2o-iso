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
 * SortGridForest.cpp
 *
 *  Created on: May 20, 2010
 *      Author: Marco.Maneta
 */

#include "Forest.h"

vectCells Forest::SortGrid(){

    REAL4 value = 0;

    UINT4 r,c;
    UINT4 nr = _patches->r;
    UINT4 nc = _patches->c;

    vectCells NotSortedArray; //array to hold the non-mv cells as read left to right and top to bottom
/*
    #pragma omp parallel for\
     private(r,c, value)
*/
    for(r=1; r < nr-1; r++)
    {
     for(c=1; c < nc-1; c++)
     {
      value = (int)_patches->matrix[r][c];
      if (value != _patches->nodata)
          NotSortedArray.cells.push_back(cell(r, c, (int)value));
     }
    }

   return NotSortedArray;

}
