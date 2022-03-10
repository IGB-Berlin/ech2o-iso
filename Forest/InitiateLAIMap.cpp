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
 *  Created on: Oct 16, 2009
 *      Author: Marco Maneta
 */

#include "Forest.h"
#include "ConstAndFuncs.h"

UINT4 Forest::InitiateLAIMap(ifstream &ifHandle, grid &DataMap){

  char comment[256];
  UINT4 nZns; //number of zones as read from the climatic data file
  UINT4  nTs; //number of time steps
  UINT4 *Zns = NULL; //array to hold the zone ids
  float *TS = NULL;
  float *data = NULL;
  UINT4 data_written = 0;
  UINT4 nzones = 1;
  
  try{
    //reads 256 bytes of comments
    ifHandle.read(comment, sizeof(comment));
    //reads one int with the number of time steps
    ifHandle.read((char *)&nTs, sizeof(int));
    TS = new float[nTs];
    //reads nTs floats into the TS array
    ifHandle.read((char *)TS, sizeof(float)*nTs);
    //reads one int with the number of zones
    ifHandle.read((char *)&nZns, sizeof(int));
    
    // reads nZns ints with zone id information into array Zns
    Zns = new UINT4[nZns];
    ifHandle.read((char *)Zns, sizeof(UINT4)*nZns); //read as many zones as there are zones in teh map

    //    cout << "nZns: " << nZns << " |nzones: "<< nzones<< endl;
    
    if(nZns- nzones < 0){
      cout << "FATAL ERROR: there are " << nzones << " counted in the map " << "and only " << nZns << "in the climate dataset" << endl;
      throw(EXIT_FAILURE);
    }

    data = new float[nZns]; //creates the array to hold the data
    ifHandle.read((char *)data, sizeof(1)); //reads data for all zones

    //DEBUG
    //cout << *data << endl;
    
    int r, c;
    
    //store the index of the climate data array that corresponds to the clim zones map
    // so when the map is updated we do not have to cycle anymore though all the climate data
    //that is not needed
    for (unsigned int k = 0; k < _vSortedGrid.cells.size() ; k++)
      {
	r = _vSortedGrid.cells[k].row;
	c = _vSortedGrid.cells[k].col;
	
	DataMap.matrix[r][c] = data[0];
	data_written++;
	
      }
    
    delete[] TS;
    delete[] data;
    delete[] Zns;
    
  }catch(int &i){//TODO: clear this crap and implement a decent exception handling system
    
    if(TS)
      delete[] TS;
    if(data)
      delete[] data;
    if(Zns)
      delete[] Zns;
    
    throw;
    
  }
  return data_written;
  
}
