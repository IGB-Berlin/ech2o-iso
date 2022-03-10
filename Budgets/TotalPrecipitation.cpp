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
 * TotalPrecipitation.cpp
 *
 *  Created on: Mar 8, 2010
 *      Author: Marco Maneta, Sylvain Kuppel
 */

#include "Budget.h"

void Budget::TotalPrecipitation(const grid* map1, const grid* map2, const Atmosphere *b)
{
  precipitation += AccountFluxes(map1, map2, b);
}

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
