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
 * TotalOvlndFlow.cpp
 *
 *  Created on: Jul 29, 2010
 *      Author: Marco Maneta, Sylvain Kuppel
 */

#include "Budget.h"

void Budget::TotalOvlndFlow(const vectCells *timeseries, const Basin *b)
{
	ovlndflow += AccountFluxes(timeseries, b);
}

void Budget::TotalOvlndFlow_d2H(const vectCells* timeseries1, const vectCells* timeseries2)
{
  ovlndflow_d2H += AccountTrckFluxes(timeseries1, timeseries2);
  //ovlndflow_d2H = AccountTrckFluxes(timeseries1, timeseries2);
}

void Budget::TotalOvlndFlow_d18O(const vectCells* timeseries1, const vectCells* timeseries2)
{
  ovlndflow_d18O += AccountTrckFluxes(timeseries1, timeseries2);
  //ovlndflow_d18O = AccountTrckFluxes(timeseries1, timeseries2);
}

// the water that already left is kept in the balance and "aging" as well
void Budget::TotalOvlndFlow_Age(const vectCells* timeseries1, const vectCells* timeseries2)
{
  ovlndflow_Age += ovlndflow * dt / 86400 + AccountTrckFluxes(timeseries1, timeseries2);
  //ovlndflow_Age = AccountTrckFluxes(timeseries1, timeseries2);
}

// Instantaneous d2H reporting
void Budget::InstOvlndFlow_d2H(const vectCells* timeseries1, const vectCells* timeseries2)
{
  d2HOvlndOut = AccountTrckFluxes2(timeseries1, timeseries2);
}
// Instantaneous d18O reporting
void Budget::InstOvlndFlow_d18O(const vectCells* timeseries1, const vectCells* timeseries2)
{
  d18OOvlndOut = AccountTrckFluxes2(timeseries1, timeseries2);
}
// Instantaneous age reporting
void Budget::InstOvlndFlow_Age(const vectCells* timeseries1, const vectCells* timeseries2)
{
  AgeOvlndOut = AccountTrckFluxes2(timeseries1, timeseries2);
}
