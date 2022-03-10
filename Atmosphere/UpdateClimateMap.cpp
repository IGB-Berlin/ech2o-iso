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
 * UpdateClimateMap.cpp
 *
 *  Created on: Oct 17, 2009
 *      Author: Marco Maneta
 */

#include "Atmosphere.h"
#include "ConstAndFuncs.h"

int Atmosphere::UpdateClimateMap(ifstream &ifHandle, grid &ClimMap){

  float *data;
  int data_written = 0;

  data = new float[_NZns]; //creates the array to hold the data

  ifHandle.read((char *)data, sizeof(float)*_NZns); //reads data for all zones

  int r, c;

#pragma omp parallel default(none) private(r,c) shared(data_written, data, ClimMap)
  {
#pragma omp for reduction(+:data_written)
  for (unsigned int i = 0; i < _vSortedGrid.size() ; i++)
    for (unsigned int j = 0; j < _vSortedGrid[i].cells.size() ; j++)
      {
	r = _vSortedGrid[i].cells[j].row;
	c = _vSortedGrid[i].cells[j].col;

	ClimMap.matrix[r][c] = data[_zoneId[i]];
	data_written++;
      }
  }


  delete[] data;

  return data_written;

}
