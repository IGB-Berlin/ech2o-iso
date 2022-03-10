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
 * UpdateSoilMoist.cpp
 *
 *  Created on: Mar 8, 2010
 *      Author: Marco Maneta
 */
#include "Basin.h"

void Basin::SoilEvapotranspiration(REAL8 &LE, //input latent heat
				   REAL8 Ts, //input surface temperature
				   REAL8 lambda,
				   REAL8 rs,// input the potential exfiltration capacity
				   REAL8 &etp, //output soil evaporation
				   REAL8 &theta,//output updated soil moisture
				   REAL8 dt, //time step
				   UINT4 r,
				   UINT4 c)
{
  REAL8 ETP;
  REAL8 theta_r = _theta_rL1->matrix[r][c];
  REAL8 sd = _depth_layer1->matrix[r][c];
  REAL8 le = lambda; //Ts < 0 ?  lat_heat_vap + lat_heat_fus : lat_heat_vap;
  REAL8 SW = std::max<REAL8>(0.0,(theta - theta_r));
  
  if(LE<0){
    etp = 0;
    return;
  }
  etp = min<REAL8>(1/rs , LE/(rho_w*le));
  ETP = etp * dt;

  if (SW * sd < ETP){ //makes sure we are not evaporating more that it is available above residual theta
    ETP = SW * sd;
    etp = ETP / dt;
  }

  // Correct the LE for the maximum water available
  if(LE>0){
    LE = etp * (rho_w*le);
  }

  theta -= ETP / sd;
}
