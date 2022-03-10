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
 * GrndHeat.cpp
 *
 *  Created on: Nov 20, 2009
 *      Author: marco.maneta
 */

#include"Basin.h"

double Basin::GrndHeat(Atmosphere &atm, Control &ctrl, const double &theta, const double &Ts, const double &Td, int row, int col){


  double Ts_old = 0; 	// temperature of surface at time step t-1
  double d = 0; 	// damping depth m
  double C = 0; 	// soil heat capacity Jm-3C-1
  double K = 0; 	// soil thermal conductivity Wm-1C-1
  double n = 0; 	// porosity
  float dt = ctrl.dt;
  double P = dt < 86400 ? 86400*14 : 31536000; //period is biweekly if less than one day (same as 1day:1year ratio) and yearly if timestep is daily or larger

  Ts_old = _Temp_s_old->matrix[row][col];

  n = _porosityL1->matrix[row][col];

  C = SoilHeatCapacity(_soil_dry_heatcap->matrix[row][col],
		       n, theta, Ts,atm.getPressure()->matrix[row][col]);
  K = SoilHeatConductivity(_soil_dry_thermcond->matrix[row][col],
			   n, theta);

  d = sqrt( (K/C) / ( 2 * ( 2 * PI / P) ) );

  return   C * d * (Ts_old - Ts)/dt + C * d * 2 * PI * (Td - Ts) / P;

}
