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
 *    Marco Maneta, Aaron Smith
 *******************************************************************************/
/*
 * BCupstreamMixing.cpp
 *
 *  Created on: Nov 24, 2020
 *      Author: Aaron Smith
 */

#include "Basin.h"

void Tracking::BCupstreamMixing(Basin &bsn, Control &ctrl, double BCsrf, double BCGW, 
				double BCdeepGW, double dx, double dt, int r, int c)
{
  // Deuterium
  if(ctrl.sw_2H){
    _Fd2HLattoChn->matrix[r][c] += (BCsrf*dt) / (dx * dx) * _d2HsurfaceBC->matrix[r][c];    
    _Fd2HLattoGW->matrix[r][c] += (BCGW * dt/dx) * _d2HgroundwaterBC->matrix[r][c];
    if(ctrl.sw_deepGW)
      _Fd2HLattoDeepGW->matrix[r][c] += (BCdeepGW * dt/dx) * _d2HdeepGWBC->matrix[r][c];
  }
  
  // Oxygen 18
  if(ctrl.sw_18O){
    _Fd18OLattoChn->matrix[r][c] += (BCsrf*dt) / (dx * dx) * _d18OsurfaceBC->matrix[r][c];
    _Fd18OLattoGW->matrix[r][c] += (BCGW * dt/dx) * _d18OgroundwaterBC->matrix[r][c];
    if(ctrl.sw_deepGW)
      _Fd18OLattoDeepGW->matrix[r][c] += (BCdeepGW * dt/dx) * _d18OdeepGWBC->matrix[r][c];    
  }
  
  // Water age
  if(ctrl.sw_Age){
    _FAgeLattoChn->matrix[r][c] += (BCsrf*dt) / (dx * dx) * _AgesurfaceBC->matrix[r][c];
    _FAgeLattoGW->matrix[r][c] += (BCGW * dt/dx) * _AgegroundwaterBC->matrix[r][c];
    if(ctrl.sw_deepGW)
      _FAgeLattoDeepGW->matrix[r][c] += (BCdeepGW * dt/dx) * _AgedeepGWBC->matrix[r][c];        
  }
  
}


