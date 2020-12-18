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

  if(UpdateClimateMap(ifLdown, *_Ldown)!=_vSsortedGridTotalCellNumber)
    throw;
  if(UpdateClimateMap(ifSdown, *_Sdown)!=_vSsortedGridTotalCellNumber)
    throw;
  if(UpdateClimateMap(ifTp, *_Tp)!=_vSsortedGridTotalCellNumber)
    throw;
  if(UpdateClimateMap(ifMaxTp, *_MaxTp)!=_vSsortedGridTotalCellNumber)
    throw;
  if(UpdateClimateMap(ifMinTp, *_MinTp)!=_vSsortedGridTotalCellNumber)
    throw;
  if(UpdateClimateMap(ifPrecip, *_Precip)!=_vSsortedGridTotalCellNumber)
    throw;
  AdjustPrecip(); // adjust precipitation with the isohyet map
  if(UpdateClimateMap(ifRelHumid, *_Rel_humid)!=_vSsortedGridTotalCellNumber)
    throw;
  if(UpdateClimateMap(ifWindSpeed, *_Wind_speed)!=_vSsortedGridTotalCellNumber)
    throw;
  
  // Tracking
  if(ctrl.sw_trck && ctrl.sw_2H){
    if(UpdateClimateMap(ifd2Hprecip, *_d2Hprecip)!=_vSsortedGridTotalCellNumber)
      throw;}
  if(ctrl.sw_trck && ctrl.sw_18O){
    if(UpdateClimateMap(ifd18Oprecip, *_d18Oprecip)!=_vSsortedGridTotalCellNumber)
      throw;}
  
  return EXIT_SUCCESS;
  
}
