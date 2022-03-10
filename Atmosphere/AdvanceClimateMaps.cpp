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
 * AdvanceClimateMaps.cpp
 *
 *  Created on: Oct 18, 2009
 *      Author: Marco Maneta
 */

#include "Atmosphere.h"

int Atmosphere::AdvanceClimateMaps(Control &ctrl){

  int grid = (int) _vSsortedGridTotalCellNumber;
  
  if(UpdateClimateMap(ifLdown, *_Ldown)!=grid)
    throw string("Error advancing long-wave time-step");
  if(UpdateClimateMap(ifSdown, *_Sdown)!=grid)
    throw string("Error advancing short wave time-step");
  if(UpdateClimateMap(ifTp, *_Tp)!=grid)
    throw string("Error advancing avg air temp time-step");
  if(UpdateClimateMap(ifMaxTp, *_MaxTp)!=grid)
    throw string("Error advancing max air temp time-step");
  if(UpdateClimateMap(ifMinTp, *_MinTp)!=grid)
    throw string("Error advancing min air temp time-step");
  if(UpdateClimateMap(ifPrecip, *_Precip)!=grid)
    throw string("Error advancing precipitation time-step");
  AdjustPrecip(); // adjust precipitation with the isohyet map
  if(UpdateClimateMap(ifRelHumid, *_Rel_humid)!=grid)
    throw string("Error advancing RH time-step");
  if(UpdateClimateMap(ifWindSpeed, *_Wind_speed)!=grid)
    throw string("Error advancing wind speed time-step");
  if(UpdateClimateMap(ifPa, *_Pa)!=grid)
    throw string("Error advancing avg air pressure time-step");
  // --------------------------------------
  //Urban ech2o
  // --------------------------------------
  if(ctrl.sw_anthr_heat)
    if(UpdateClimateMap(ifAnthrHeat, *_Anthr_Heat)!=grid)
      throw string("Error advancing anthropogenic heat time-step");  
  // --------------------------------------
  //Tracking
  // --------------------------------------
  if(ctrl.sw_trck && ctrl.sw_2H){
    if(UpdateClimateMap(ifd2Hprecip, *_d2Hprecip)!=grid)
      throw string("Error advancing d2H precipitation time-step");}
  if(ctrl.sw_trck && ctrl.sw_18O){
    if(UpdateClimateMap(ifd18Oprecip, *_d18Oprecip)!=grid)
      throw string("Error advancing d18O precipitation time-step");}
  return EXIT_SUCCESS;
  
}
