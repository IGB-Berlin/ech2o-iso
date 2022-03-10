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
 *    Marco Maneta, Sylvain Kuppel, Aaron Smith
 *******************************************************************************/
/*
 * UpdateBCMap_Iso.cpp
 *
 *  Created on: Nov 12, 2020
 *      Author: Aaron Smith
 */

#include "Basin.h"
#include "ConstAndFuncs.h"

int Tracking::UpdateBCMap_Iso(ifstream &ifHandle, grid &BCMap, Atmosphere &atm){

  float *data;
  int data_written = 0;

  data = new float[atm.getNZns()]; //creates the array to hold the data

  ifHandle.read((char *)data, sizeof(float)*atm.getNZns()); //reads data for all zones

  int r, c;

  for (unsigned int a = 0; a < atm.getSortedGrid().size(); a++ ){ //loops only over the number of zones in the climate zone map, not in the climate dataset
    for (unsigned int k = 0; k < atm.getSortedGrid()[a].cells.size() ; k++)
      {
	r = atm.getSortedGrid()[a].cells[k].row;
	c = atm.getSortedGrid()[a].cells[k].col;
      
	BCMap.matrix[r][c] = data[atm.getzoneId()[a]];
	data_written++;

      }
  }


  delete[] data;

  return data_written;

}
