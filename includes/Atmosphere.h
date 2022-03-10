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
 * Atmosphere.h
 *
 *  Created on: Oct 14, 2009
 *      Author: Marco Maneta
 */

#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include <fstream>
#include "Grid.h"
#include "InitConf.h"
#include "SortGrid.h"
//#include "exceptions_echo.h"

using namespace std;

class Atmosphere{

  UINT4 _NRows;
  UINT4 _NCols;
  REAL8 _dx;
  grid *_zones; //map with the zone extents and identifiers
  UINT4 _nzones; // number of zones in map
  UINT4 _NZns; //number of zones in binary data time series
  UINT4 *_zoneId; // array with zone identifiers from the climate dataset
  vector<vectCells> _vSortedGrid;
  grid *_isohyet; //map with rainfall multipliers to spatially distribute precipitation
  UINT4 _vSsortedGridTotalCellNumber;
  
  REAL8 _rain_snow_temp; //scalar with air temperature threshold for rain-snow partition
  
  grid *_Ldown; //longwave downward radiation Wm-2
  grid *_Sdown; //shortwave downward radiation Wm-2
  grid *_Tp; //Average air temperature C
  grid *_MaxTp; //Maximum air temperature C
  grid *_MinTp; //Minimum air temperature C
  grid *_Precip; // Precipitation intensity ms-1
  grid *_Rel_humid; //relative humidity of air [0-1]
  grid *_Wind_speed; //windspeed in ms-1
  grid *_Pa; //atmospheric pressure in Pa
  grid *_Anthr_Heat; //anthropogenic heat [W.m-2]
  grid *_d2Hprecip; //Isotopic signature of precipitation (2H per mil)
  grid *_d18Oprecip; //Isotopic signature of precipitation (18O per mil)
  
  void CountNumZones();
  vectCells SortGrid(int zoneId);
  UINT4 InitiateClimateMap(ifstream & ifHandle, grid & ClimMap);
  // internal function that updates a climate map
  int UpdateClimateMap(ifstream & ifHandle, grid & ClimMap);
  int AdjustPrecip();
  
  //climate data file handles
  ifstream ifLdown;
  ifstream ifSdown;
  ifstream ifTp;
  ifstream ifMaxTp;
  ifstream ifMinTp;
  ifstream ifPrecip;
  ifstream ifRelHumid;
  ifstream ifWindSpeed;
  ifstream ifPa;
  ifstream ifAnthrHeat;
  ifstream ifd2Hprecip;
  ifstream ifd18Oprecip;
 
 public:
  
  Atmosphere();
  Atmosphere(Control &ctrl);
  ~Atmosphere();
  
  int AdvanceClimateMaps(Control &ctrl); //external interface that updates all climate maps by calling UpdateClimateMap
  
  /*Getters and setters*/
  //get methods (inline)
  
  REAL8 getCellSize() const {
    return _dx;
  }

  UINT4 getNZns() const {
    return _NZns;
  }

  UINT4 getnzones() const {
    return _nzones;
  }

  UINT4 *getzoneId() const {
    return _zoneId;
  }

  UINT4 getSsortedGridTotalCellNumber() const{
    return _vSsortedGridTotalCellNumber;
  }

  grid *getzones() const {
    return _zones;
  }
  
  const vector<vectCells> &getSortedGrid() const {
    return _vSortedGrid;
  }
  
  REAL8 getRainSnowTempThreshold() const{
    return _rain_snow_temp;
  }
  grid *getIncomingLongWave() const
  {
    return _Ldown;
  }
  
  grid *getIncomingShortWave() const
  {
    return _Sdown;
  }
  
  grid *getTemperature() const
  {
    return _Tp;
  }
  grid *getMaxTemperature() const
  {
    return _MaxTp;
  }
  grid *getMinTemperature() const
  {
    return _MinTp;
  }
  grid *getPrecipitation() const
  {
    return _Precip;
  }
  
  grid *getRelativeHumidty() const
  {
    return _Rel_humid;
  }
  
  grid *getWindSpeed() const
  {
    return _Wind_speed;
  }
  
  grid *getPressure() const
  {
    return _Pa;
  }

  grid *getAnthrHeat() const
  {
    return _Anthr_Heat;
  }

  // Isotope tracking
  grid *getd2Hprecip() const
  {
    return _d2Hprecip;
  }
  grid *getd18Oprecip() const
  {
    return _d18Oprecip;
  }
  
  //setter
  void setPrecip(UINT4 row, UINT4 col, REAL8 value)
  {
    this->_Precip->matrix[row][col] = value;
  }
};
#endif /* ATMOSPHERE_H_ */
