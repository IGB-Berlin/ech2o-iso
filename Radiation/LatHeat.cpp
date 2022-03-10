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
 * LatHeat.cpp
 *
 *  Created on: Nov 19, 2009
 *      Author: marco.maneta
 */

#include"Basin.h"

double Basin::LatHeat(Atmosphere &atm, double beta, double ra, double rs, double rc, const double &Ts, int row, int col){

  double airdens = 0; //density of air in kgm-3
  double es = 0; //saturated vapor pressure at temperature Ts
  double ea = 0; //vapor pressure at air temperature
  double gamma = 0;
  double z = 0;
  
  airdens = AirDensity(atm.getTemperature()->matrix[row][col], atm.getPressure()->matrix[row][col]); //kgm-3
  es = SatVaporPressure(Ts); //vapor pressure at temp Ts in Pa
  ea = SatVaporPressure(atm.getTemperature()->matrix[row][col]) * atm.getRelativeHumidty()->matrix[row][col]; //vapor pressure at air temp in Pa
  z = _DEM->matrix[row][col];
  gamma = PsychrometricConst(atm.getPressure()->matrix[row][col], z); //todo: replace with Pressure when implemented
  
  return
    (1/gamma) * (1/(ra + rs) + rc) * airdens * spec_heat_air * beta *(ea - es);
  
}
