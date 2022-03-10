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
 *    Marco Maneta, Sylvain Kuppel
 *******************************************************************************/
/*
 * AccountArea.cpp
 *
 *  Created on: Nov 7, 2016
 *      Author: Sylvain Kuppel
 */
#include "Budget.h"

double Budget::AccountRelArea(const grid *map1, const grid *map2, const Basin *b)
{

  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 result2= 0;
  REAL8 dx = b->getCellSize();
#pragma omp parallel for			\
  default(shared) private(r,c)			\
  reduction (+:result)

  for (UINT4 i = 0; i< length; i++){
    
    r = b->getSortedGrid().cells[i].row;
    c = b->getSortedGrid().cells[i].col;
    //result += (map->matrix[r][c]/length);
    result2+= dx*dx*map2->matrix[r][c];
    result += (map1->matrix[r][c]*dx*dx*map2->matrix[r][c]);
  }
  
  return result/result2;
}
