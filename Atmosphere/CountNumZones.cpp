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
 * CountNumZones.cpp
 *
 *  Created on: Oct 15, 2009
 *      Author: Marco Maneta
 */

#include "Atmosphere.h"

void Atmosphere::CountNumZones(){

  unsigned int r,c,i;
  bool found = false;
  _nzones = 0;
  //makes space for potentially as many zones as number of cells in the map
  int *zone = new int[_NRows*_NCols]; 

  //initialization to prevent valgrind form complaining on uninitialized variables
  for (unsigned int i = 0; i<_NRows*_NCols; i++) 
    zone[i] = 0;
  
  // check nr of zones
  zone[0] =(int)_zones->nodata;
  for ( r = 0; r < _NRows; r++)
    for ( c = 0; c < _NCols; c++)
      if (_zones->matrix[r][c] != _zones->nodata){
	found = false;
	for( i = 0; i <= _nzones; i++)
	  if (zone[i] == _zones->matrix[r][c]) found = true;
	
	if(!found && i < _NRows*_NCols){
	  zone[_nzones] = (int)_zones->matrix[r][c];
	  _vSortedGrid.push_back(vectCells(zone[_nzones]));
	  _nzones++;
	}
      }

  _zoneId = new UINT4[_nzones]; //will be destroyed in Atmosphere destructor

  delete[] zone;
}
