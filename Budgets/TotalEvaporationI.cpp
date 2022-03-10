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
 * TotalEvaporationI.cpp
 *
 *  Created on: Jan 17, 2017
 *      Author: Sylvain Kuppel
 */

#include "Budget.h"

void Budget::TotalEvaporationI(const grid* map1, const grid* map2, const Basin *b)
{
        evaporationI += AccountFluxes(map1, map2, b);
}

void Budget::TotalEvaporationI_d2H(const grid* map1, const grid* map2,
				   const grid* map3, const Basin *b)
{
  evaporationI_d2H += AccountTrckFluxes(map1, map2, map3, b);
  //evaporationI_d2H = AccountTrckFluxes(map1, map2, b);
}

void Budget::TotalEvaporationI_d18O(const grid* map1, const grid* map2,
				    const grid* map3, const Basin *b)
{
  evaporationI_d18O += AccountTrckFluxes(map1, map2, map3, b);
  //evaporationI_d18O = AccountTrckFluxes(map1, map2, b);
}

// the water that already left is kept in the balance and "aging" as well
void Budget::TotalEvaporationI_Age(const grid* map1, const grid* map2,
				   const grid* map3, const Basin *b)
{
  evaporationI_Age += evaporationI * dt /86400 + AccountTrckFluxes(map1, map2, map3, b);
  //evaporationI_Age = AccountTrckFluxes(map1, map2, b);
}

// Instantaneous tracer reporting
void Budget::InstEvaporationI_d2H(const grid* map1, const grid* map2, const Basin *b)
{
  d2HevapI = AccountTrckFluxes2(map1, map2, b);
}
void Budget::InstEvaporationI_d18O(const grid* map1, const grid* map2, const Basin *b)
{
  d18OevapI = AccountTrckFluxes2(map1, map2, b);
}
void Budget::InstEvaporationI_Age(const grid* map1, const grid* map2, const Basin *b)
{
  AgeevapI = AccountTrckFluxes2(map1, map2, b);
}
