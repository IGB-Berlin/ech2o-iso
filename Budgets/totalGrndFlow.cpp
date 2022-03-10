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
 * TotalGrndFlow.cpp
 *
 *  Created on: Dec 10, 2010
 *      Author: Marco.Maneta
 */

#include "Budget.h"

void Budget::TotalGrndFlow(const vectCells *timeseries, const Basin *b)
{
  gwtrflow += AccountFluxes(timeseries, b);
}

void Budget::TotalGrndFlow_d2H(const vectCells* timeseries1, const vectCells* timeseries2)
{
  gwtrflow_d2H += AccountTrckFluxes(timeseries1, timeseries2);
}

void Budget::TotalGrndFlow_d18O(const vectCells* timeseries1, const vectCells* timeseries2)
{
  gwtrflow_d18O += AccountTrckFluxes(timeseries1, timeseries2);
}

// the water that already left is kept in the balance and "aging" as well
void Budget::TotalGrndFlow_Age(const vectCells* timeseries1, const vectCells* timeseries2)
{
  gwtrflow_Age += gwtrflow * dt / 86400 + AccountTrckFluxes(timeseries1, timeseries2);
}

// Instantaneous d2H reporting
void Budget::InstGrndFlow_d2H(const vectCells* timeseries1, const vectCells* timeseries2)
{
  d2HGWOut = AccountTrckFluxes2(timeseries1, timeseries2);
}
// Instantaneous d18O reporting
void Budget::InstGrndFlow_d18O(const vectCells* timeseries1, const vectCells* timeseries2)
{
  d18OGWOut = AccountTrckFluxes2(timeseries1, timeseries2);
}
// Instantaneous Age reporting
void Budget::InstGrndFlow_Age(const vectCells* timeseries1, const vectCells* timeseries2)
{
  AgeGWOut = AccountTrckFluxes2(timeseries1, timeseries2);
}
