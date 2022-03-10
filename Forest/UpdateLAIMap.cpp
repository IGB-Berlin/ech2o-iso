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
 * UpdateLAIMap.cpp
 *
 *  Created on: Dec 7, 2018
 *      Author: Sylvain Kuppel
 */

#include "Forest.h"
#include "ConstAndFuncs.h"

int Forest::UpdateLAIMap(ifstream &ifHandle, grid &DataMap){

  float *data;
  int data_written = 0;

  data = new float[1]; //creates the array to hold the data

  ifHandle.read((char *)data, sizeof(1)); //reads data for all zones
  
  int r, c;

  for (unsigned int k = 0; k < _vSortedGrid.cells.size() ; k++) {
      r = _vSortedGrid.cells[k].row;
      c = _vSortedGrid.cells[k].col;
      
      DataMap.matrix[r][c] = data[0];
      data_written++;

  }

  delete[] data;

  return data_written;

}
