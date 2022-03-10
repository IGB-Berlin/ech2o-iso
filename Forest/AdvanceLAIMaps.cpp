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
 * AdvanceLAIMaps.cpp
 *
 *  Created on: Dec 7, 2018
 *      Author: Sylvain Kuppel
 */

#include "Forest.h"

int Forest::AdvanceLAIMaps(){

  for(UINT4 i = 0; i < _Nsp - 1; i++){
    UpdateLAIMap(_species[i].ifLAI, *_species[i]._LAI);
    // Now also update the height map
    UpdateLAIMap(_species[i].ifhgt, *_species[i]._Height);    
  }
   
  return EXIT_SUCCESS;
  
}
