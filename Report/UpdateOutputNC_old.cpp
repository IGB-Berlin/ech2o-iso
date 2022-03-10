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
 * UpdateOutputNC.cpp
 *
 *  Created on: Feb 10, 2021
 *      Author: Xiaoqiang Yang
 *      This script includes two functions creating (ts = 1) and updating the nc files
 *      for water fluxes (Water_Fluxes_States.nc), deutrium tracer (D_Fluxes_States.nc),
 *      oxygen-18 tracer (O18_Fluxes_States.nc), water age (Age_Fluxes_States.nc) and
 *      vegtation dynamics (VegDyn_Fluxes_States.nc).
 *  Technicial notes:
 *      netCDF C++4 library should be installed (e.g., in UFZ eve-cluster: module load netCDF-C++4) 
 *      and configured in the compiling folder (e.g., add "-lnetcdf_c++4" in ./Release-Linux/objects.mk)
 */
#include <fstream>
#include <netcdf>
#include "Report.h"
#include "Basin.h"
#include "Sativa.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;



// Names of things. 
#define LAT_NAME "latitude"
#define LON_NAME "longitude"
#define REC_NAME "time"

string  UNITS = "units";
string  DEGREES_EAST =  "degrees_east";
string  DEGREES_NORTH = "degrees_north";

// For the units attributes. 
string LAT_UNITS = "degrees_north";
string LON_UNITS = "degrees_east";

// Return this code to the OS in case of failure.
#define NC_ERR 2

int Report::UpdateOutputNC(const grid *input, string varname, string outTP){
  
  UINT4 NLAT, NLON;
  int TScount;
  string filename;
  //from Basin
  NLAT = oBasin->getNumRows();
  NLON = oBasin->getNumCols();  
  TScount = oControl->current_ts_count;
  REAL8 outdata[NLAT][NLON];   
  //get the outdata from input map
  for(unsigned i = 0; i < NLAT; i++){
    for(unsigned j = 0; j < NLON; j++){
      outdata[i][j] = input->matrix[i][j];
    }
  }   
  //indexing the record location
  //c++ starts with index 0 AND the first nodata record should be replaced so
  TScount -= 1; 
  
  try
    {
      //open the right nc file according to the output type "outTP"
      if(outTP == "W") {
	filename =  oControl->path_ResultsFolder + "Water_Fluxes_States.nc";
      }else if(outTP == "TD"){
	filename =  oControl->path_ResultsFolder + "D_Fluxes_States.nc";
      }else if(outTP == "TO"){
	filename =  oControl->path_ResultsFolder + "O18_Fluxes_States.nc";	
      }else if(outTP == "TA"){
	filename =  oControl->path_ResultsFolder + "Age_Fluxes_States.nc";	
      }else if(outTP == "VG"){
	filename =  oControl->path_ResultsFolder + "VegDyn_Fluxes_States.nc";
      }
      //open the file
      NcFile file1(filename, NcFile::write);
      if(file1.isNull()){
	cout << "ERROR*** file " << filename << " doesn't exist!" << endl;
	return NC_ERR;
      }
      /* 	// Get the  latitude and longitude coordinate variables 
		NcVar latVar, lonVar;
		latVar = dataFile.getVar("latitude");
		if(latVar.isNull()) return NC_ERR;
		lonVar = dataFile.getVar("longitude");
		if(lonVar.isNull()) return NC_ERR; */  
      //Get the record variable
      NcVar recVar;
      recVar = file1.getVar(varname.c_str());
      if(recVar.isNull()){
	cout << "ERROR*** variable " << varname.c_str() << " doesn't exist!" << endl;
	return NC_ERR;
      }

      //Get using the current_ts_count for the start index of new record  
      //cout << varname.c_str() << "existing number of records: " << TScount+1 << endl;
      vector<size_t> startp,countp;
      startp.push_back(TScount);
      startp.push_back(0);
      startp.push_back(0);
      countp.push_back(1);
      countp.push_back(NLAT);
      countp.push_back(NLON);
      recVar.putVar(startp,countp,outdata);
 
      // close of nc takes place in destructor
      return 0;
    }
  catch(NcException& e)
    {
      e.what();
      cout<<"FAILURE**You've probably done something wrong!!"<<endl;	  
      return NC_ERR;
    }	
}
