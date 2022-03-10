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
 * CSF2Grid.cpp
 *
 *  Created on: Oct 13, 2009
 *      Author: Marco Maneta
 */

#include <fstream>
#include <cmath>
#include "Grid.h"
#include "ConstAndFuncs.h"

int CSF2Grid(std::string fname, grid *pnt){

  UINT4 cells = 0; REAL8 data;
  MAP *pcrMap = Mopen(fname.c_str(), M_READ);

  if (Merrno != 0)
    {cerr << MstrError() << ": " << fname << endl; exit(-1);}

  RuseAs(pcrMap, CR_REAL8);
  
  pnt->c = RgetNrCols(pcrMap);
  pnt->r = RgetNrRows(pcrMap);
  pnt->dx = RgetCellSize(pcrMap);
  pnt->nodata = -9999;
  pnt->north = RgetYUL(pcrMap);
  pnt->west = RgetXUL(pcrMap);
  pnt->east = pnt->west + (pnt->dx * pnt->c);

  if (MgetProjection(pcrMap) == PT_YINCT2B)
    pnt->south = pnt->north + (pnt->dx * pnt->r);
  if (MgetProjection(pcrMap) == PT_YDECT2B)
    pnt->south = pnt->north - (pnt->dx * pnt->r);

  //initializes the grid matrix
  try{
    pnt->matrix = new REAL8*[pnt->r];
    for (UINT4 i=0; i < pnt->r; i++)
      pnt->matrix[i] = new REAL8[pnt->c];
  }catch(std::bad_alloc&){cerr << "unable to allocate memory...";
    cin.get(); exit(-1);}
  
  for (UINT4 i=0; i<pnt->r; i++)
    {
      for (UINT4 j=0; j<pnt->c; j++)
	{
	  cells +=RgetCell(pcrMap, i, j, &data);
	  if(IsMV(pcrMap, &data)) data = pnt->nodata;
	  pnt->matrix[i][j] = data;
	}
    }
  Mclose(pcrMap);

  if(cells!=(pnt->r*pnt->c)) return -1;
  return cells;
}
