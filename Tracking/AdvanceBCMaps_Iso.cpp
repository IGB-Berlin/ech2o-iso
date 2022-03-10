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
 * AdvanceBCMaps_Iso.cpp
 *
 *  Created on: Nov 12, 2020
 *      Author: Aaron Smith
 */

#include "Sativa.h"
#include "Basin.h"


int Tracking::AdvanceBCMaps_Iso(Atmosphere &atm, Control &ctrl, int iso) {

  if(iso == 0){ // deuterium
    // Update for surface water
    UpdateBCMap_Iso(ifd2HsurfaceBC, *_d2HsurfaceBC, atm);
    // Update for groundwater
    UpdateBCMap_Iso(ifd2HgroundwaterBC, *_d2HgroundwaterBC, atm);
    if(ctrl.sw_deepGW) // Update deep groundwater
      UpdateBCMap_Iso(ifd2HdeepGWBC, *_d2HdeepGWBC, atm);
  }
  if(iso == 1){ //oxygen-18
    // Update for surface water
    UpdateBCMap_Iso(ifd18OsurfaceBC, *_d18OsurfaceBC, atm);
    // Update for groundwater
    UpdateBCMap_Iso(ifd18OgroundwaterBC, *_d18OgroundwaterBC, atm);
    if(ctrl.sw_deepGW) // Update deep groundwater
      UpdateBCMap_Iso(ifd18OdeepGWBC, *_d18OdeepGWBC, atm);    
  }
  if(iso == 2){ //age
    // Update for surface water
    UpdateBCMap_Iso(ifAgesurfaceBC, *_AgesurfaceBC, atm);
    // Update for groundwater
    UpdateBCMap_Iso(ifAgegroundwaterBC, *_AgegroundwaterBC, atm);
    if(ctrl.sw_deepGW) // Update deep groundwater
      UpdateBCMap_Iso(ifAgedeepGWBC, *_AgedeepGWBC, atm);        
  }


  return EXIT_SUCCESS;
}
