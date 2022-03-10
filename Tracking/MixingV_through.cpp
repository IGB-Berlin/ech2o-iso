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
 * MixingV_through.cpp
 *
 *  Created on: Jun 21, 2017
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"

void Tracking::MixingV_through(Atmosphere &atm, Basin &bsn, Control &ctrl,
			       double &rain, double &p, int s, int r, int c) //time step
{
  
  double dt = ctrl.dt;
  double pond_old = bsn.getPondingWater()->matrix[r][c];
  //**********************************************************************************
  //--- Deuterium --------------------------------------------------------------------
  //**********************************************************************************
  if(ctrl.sw_2H)
    _d2Hsurface->matrix[r][c] = InputMix(pond_old, _d2Hsurface->matrix[r][c],
					rain*p*dt,bsn.getd2Hthroughfall(s)->matrix[r][c]);

  //**********************************************************************************
  //--- Oxygen 18 --------------------------------------------------------------------
  //**********************************************************************************
  if(ctrl.sw_18O)
    _d18Osurface->matrix[r][c] = InputMix(pond_old, _d18Osurface->matrix[r][c],
					  rain*p*dt,bsn.getd18Othroughfall(s)->matrix[r][c]);

  //**********************************************************************************  
  // Water age
  //**********************************************************************************  
  if(ctrl.sw_Age)
    _Agesurface->matrix[r][c] = InputMix(pond_old, _Agesurface->matrix[r][c],
					 rain*p*dt,bsn.getAgethroughfall(s)->matrix[r][c]);  

}


