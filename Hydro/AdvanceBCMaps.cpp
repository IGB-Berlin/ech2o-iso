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
 * AdvanceBCMaps.cpp
 *
 *  Created on: Nov 12, 2020
 *      Author: Aaron Smith
 */

#include "Sativa.h"
#include "Basin.h"


int Basin::AdvanceBCMaps(Atmosphere &atm, Control &ctrl) {

  //  cout << "testing update BCmaps" << endl;
  // Update for surface water
  UpdateBCMap(ifBCsurface, *_BCsurface, atm);
  // Update for layer 1 water
  //  UpdateBCMap(ifBClayer1, *_BClayer1, atm);
  // Update for layer 2 water
  //  UpdateBCMap(ifBClayer2, *_BClayer2, atm);
  // Update for groundwater
  UpdateBCMap(ifBCgroundwater, *_BCgroundwater, atm);
  //Update for deep groundwater
  if(ctrl.sw_deepGW)
    UpdateBCMap(ifBCdeepgwtr, *_BCdeepgwtr, atm);
  

  return EXIT_SUCCESS;
}
