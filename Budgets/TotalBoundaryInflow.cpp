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
 * TotalBoundaryInflow.cpp
 *
 *  Created on: Mar 8, 2010
 *      Author: Marco Maneta, Sylvain Kuppel, Aaron Smith
 */

#include "Budget.h"

void Budget::TotalBoundaryInflow(const grid* map1, const grid* map2,
				 const grid* map3, const grid* map4, const Basin *b, const Control *ctrl)
{
  ovlndinflow += AccountBCFluxesQ(map1,map4,b);
  grndinflow += AccountBCFluxes(map2,map4,b);
  if(ctrl->sw_deepGW)
    deepgrndinflow += AccountBCFluxes(map3,map4,b);
}

void Budget::TotalBoundaryInflow_d2H(const grid* map1, const grid* map2,
				     const grid* map3,
				     const grid* map4, const grid* map5,
				     const grid* map6,
				     const grid* map7,
				     const Basin *b, const Control *ctrl)
{
  // map 1-3 - surface - GW    - deepGW
  // map 4-6 - d2Hsurf - d2HGW - d2HdeepGW
  // map 7   - area fraction
  ovlndinflow_d2H += AccountTrckBCFluxesQ(map1,map4,map7,b);
  grndinflow_d2H += AccountTrckBCFluxes(map2,map5,map7,b);
  if(ctrl->sw_deepGW)
    deepgrndinflow_d2H += AccountTrckBCFluxes(map3,map6,map7,b);
}

void Budget::TotalBoundaryInflow_d18O(const grid* map1, const grid* map2,
				     const grid* map3,
				     const grid* map4, const grid* map5,
				     const grid* map6,
				     const grid* map7,
				     const Basin *b, const Control *ctrl)
{
  // map 1-3 - surface - GW    - deepGW
  // map 4-6 - d2Hsurf - d2HGW - d2HdeepGW
  // map 7   - area fraction
  ovlndinflow_d18O += AccountTrckBCFluxesQ(map1,map4,map7,b);
  grndinflow_d18O += AccountTrckBCFluxes(map2,map5,map7,b);
  if(ctrl->sw_deepGW)
    deepgrndinflow_d18O += AccountTrckBCFluxes(map3,map6,map7,b);
}

void Budget::TotalBoundaryInflow_Age(const grid* map1, const grid* map2,
				     const grid* map3,
				     const grid* map4, const grid* map5,
				     const grid* map6,
				     const grid* map7,
				     const Basin *b, const Control *ctrl)
{
  // map 1-3 - surface - GW    - deepGW
  // map 4-6 - d2Hsurf - d2HGW - d2HdeepGW
  // map 7   - area fraction
  ovlndinflow_Age += AccountTrckBCFluxesQ(map1,map4,map7,b);
  grndinflow_Age += AccountTrckBCFluxes(map2,map5,map7,b);
  if(ctrl->sw_deepGW)
    deepgrndinflow_Age += AccountTrckBCFluxes(map3,map6,map7,b);
}
/*
void Budget::TotalPrecipitation_d2H(const grid* map1, const grid* map2,
				    const grid* map3, const Atmosphere *atm)
{
  precipitation_d2H += AccountTrckFluxes(map1, map2, map3, atm);
  //precipitation_d2H = AccountTrckFluxes(map1, map2, atm);
}

void Budget::TotalPrecipitation_d18O(const grid* map1, const grid* map2,
				     const grid* map3, const Atmosphere *atm)
{
  precipitation_d18O += AccountTrckFluxes(map1, map2, map3, atm);
  //precipitation_d18O = AccountTrckFluxes(map1, map2, atm);
}

// the water that already entered is kept in the balance and "aging" as well
void Budget::TotalPrecipitation_Age()
{
  precipitation_Age += precipitation * dt / 86400;
  //precipitation_Age = 0;
}

// For Basin*Summary.txt
void Budget::InstPrecipitation_d2H(const grid* map1, const grid* map2, const Basin *b)
{
  d2Hprecip = AccountTrckFluxes2(map1, map2, b);
}
void Budget::InstPrecipitation_d18O(const grid* map1, const grid* map2, const Basin *b)
{
  d18Oprecip = AccountTrckFluxes2(map1, map2, b);
}
*/
