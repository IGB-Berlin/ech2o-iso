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
 * TotalDeepGrndFlow.cpp
 *
 *  Created on: May 14, 2020
 *      Author: Xiaoqiang Yang
 *  Following TotalGrndFlow.cpp
 */

#include "Budget.h"

void Budget::TotalDeepGrndFlow(const vectCells *timeseries, const Basin *b)
{
  deepgwtrflow += AccountFluxes(timeseries, b); 
}

void Budget::TotalDeepGrndFlow_d2H(const vectCells* timeseries1, const vectCells* timeseries2)
{
  deepgwtrflow_d2H += AccountTrckFluxes(timeseries1, timeseries2);
}

void Budget::TotalDeepGrndFlow_d18O(const vectCells* timeseries1, const vectCells* timeseries2)
{
  deepgwtrflow_d18O += AccountTrckFluxes(timeseries1, timeseries2);
}

// the water that already left is kept in the balance and "aging" as well
void Budget::TotalDeepGrndFlow_Age(const vectCells* timeseries1, const vectCells* timeseries2)
{
  deepgwtrflow_Age += deepgwtrflow * dt / 86400 + AccountTrckFluxes(timeseries1, timeseries2); 
}

// Instantaneous d2H reporting
void Budget::InstDeepGrndFlow_d2H(const vectCells* timeseries1, const vectCells* timeseries2)
{
  d2HdeepGWOut = AccountTrckFluxes2(timeseries1, timeseries2);
}
// Instantaneous d18O reporting
void Budget::InstDeepGrndFlow_d18O(const vectCells* timeseries1, const vectCells* timeseries2)
{
  d18OdeepGWOut = AccountTrckFluxes2(timeseries1, timeseries2);
}

// Instantaneous age reporting
void Budget::InstDeepGrndFlow_Age(const vectCells* timeseries1, const vectCells* timeseries2)
{
  AgedeepGWOut = AccountTrckFluxes2(timeseries1, timeseries2);
}
