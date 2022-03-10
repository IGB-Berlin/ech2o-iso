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

#include "Atmosphere.h"
#include "ConstAndFuncs.h"

UINT4 Atmosphere::InitiateClimateMap(ifstream &ifHandle, grid &ClimMap){

  char comment[256];
  UINT4 nZns; //number of zones as read from the climatic data file
  UINT4  nTs; //number of time steps
  UINT4 *Zns = NULL; //array to hold the zone ids
  float *TS = NULL;
  float *data = NULL;
  UINT4 data_written = 0;
  
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

    //    cout << "nZns: " << nZns << " |nzones " << _nzones << endl;
    
    if(_NZns!=0) //if this is not the first time a climatic map is read and the number of zone is being set
      if(_NZns != (UINT4)nZns){ //make sure the number of zones in all climate time series are the same
	cout << "The number of zones in two climatic time series differ. Number of zones in all climatic times series must be equal" << endl;
	throw(EXIT_FAILURE);
      }
    _NZns = nZns; //stores the number of zones of the climate time series in the object. It is rewritten for each map initialization but since there are all the same it doesnt matter

    // reads nZns ints with zone id information into array Zns
    Zns = new UINT4[nZns];
    ifHandle.read((char *)Zns, sizeof(UINT4)*nZns); //read as many zones as there are zones in teh map
        
    if(nZns-_nzones < 0){
      cout << "FATAL ERROR: there are " << _nzones << " counted in the map " << "and only " << nZns << "in the climate dataset" << endl;
      throw(EXIT_FAILURE);
    }
    //ifHandle.seekg(sizeof(int)*(nZns-_nzones), ios::cur); //skip the rest of the zones

    data = new float[nZns]; //creates the array to hold the data
    
    ifHandle.read((char *)data, sizeof(float)*nZns); //reads data for all zones

    int r, c;
#pragma omp parallel for default(none) private(r,c) shared(data,Zns, ClimMap, nZns) reduction(+:data_written)        
    for (UINT4 a = 0; a < nZns; a++ ){
      for (UINT4 i = 0; i < _vSortedGrid.size(); i++ ){
	if(_vSortedGrid[i].zone == Zns[a])
	  {
	    //store the index of the climate data array that corresponds to the clim zones map
	    // so when the map is updated we do not have to cycle anymore though all the climate data
	    //that is not needed
	    _zoneId[i] = a;
	    for (unsigned int j = 0; j < _vSortedGrid[i].cells.size() ; j++)
	      {
		r = _vSortedGrid[i].cells[j].row;
		c = _vSortedGrid[i].cells[j].col;
		
		ClimMap.matrix[r][c] = data[a];

		data_written++;
		
	      }
	  }
      }
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
