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
 * MassBalanceError.cpp
 *
 *  Created on: Mar 18, 2010
 *      Author: Marco Maneta
 */

#include "Budget.h"

void Budget::MassBalanceError()
{
  double inputs = 0.0;
  double outputs = 0.0;
  double ds = 0.0;

  inputs = precipitation + initsnowpack + initponding + initchan + initL1 + initL2 + initL3 + initGW;
	
  outputs = evaporationS + evaporationI + transpiration + ovlndflow + gwtrflow + leakage;

  ds = canopy + snowpack + ponding + chan_store + soilL1 + soilL2 + soilL3 + grndwater;

  if(inputs>0) 
    MBErr = 100/inputs*(inputs-outputs - ds);
  else 
    MBErr = 0;

}
