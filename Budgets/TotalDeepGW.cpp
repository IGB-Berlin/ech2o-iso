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
 *    Marco Maneta, Sylvain Kuppl
 *******************************************************************************/
/*
 * TotalDeepGWtoChn.cpp
 *
 *  Created on: May 14, 2020
 *      Author: Xiaoqiang Yang
 *  Following TotalGWtoChn.cpp
 */

#include "Budget.h"

void Budget::TotalDeepGW(const grid* map1, const grid* map2, const Basin *b)
{
  deep_grndwater = AccountStorages(map1, map2, b);
}
/*
void Budget::InstDeepGW_d2H(const grid* map1, const grid* map2, const Basin *b)
{
  d2Hdeepgwtochn = AccountTrckFluxes2(map1, map2, b);
}

void Budget::InstDeepGW_d18O(const grid* map1, const grid* map2, const Basin *b)
{
  d18Odeepgwtochn = AccountTrckFluxes2(map1, map2, b);
}

// the water that already left is kept in the balance and "aging" as well
void Budget::InstDeepGW_Age(const grid* map1, const grid* map2, const Basin *b)
{
  Agedeepgwtochn = AccountTrckFluxes2(map1, map2, b);

}
*/
