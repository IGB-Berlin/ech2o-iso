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
 * TotalEvaporation.cpp
 *
 *  Created on: Mar 17, 2010
 *      Author: Marco Maneta
 */

#include "Budget.h"

void Budget::TotalEvaporation(const grid* map1, const grid* map2, const Basin *b)
{
        evaporation += AccountFluxes(map1, map2, b);
}

// For Basind2HSummary.txt --------------------------------------------------------------
void Budget::InstOut_d2H(const grid* evapS, const grid* CevapS,
			 const grid* evapI, const grid* CevapI,
			 const grid* evapT, const grid* CevapT,
			 const grid* leakage, const grid* Cleakage,
			 const vectCells *OvlndOut, const vectCells *COvlndOut,
			 const vectCells *GWOut, const vectCells *CGWOut,
			 const grid* ttarea,			 
			 const Basin *b)
{
  d2HOut = AccountTrckOut(evapS, CevapS,evapI, CevapI, evapT, CevapT,
			  leakage, Cleakage, OvlndOut, COvlndOut, GWOut, CGWOut, ttarea, b);
}

void Budget::InstEvaporation_d2H(const grid* evapS, const grid* CevapS,
				 const grid* evapI, const grid* CevapI,
				 const grid* evapT, const grid* CevapT,
				 const Basin *b)
{
  d2HET = AccountTrckET(evapS, CevapS,evapI, CevapI, evapT, CevapT, b);
}

// For Basind18OSummary.txt -------------------------------------------------------------
void Budget::InstOut_d18O(const grid* evapS, const grid* CevapS,
			 const grid* evapI, const grid* CevapI,
			 const grid* evapT, const grid* CevapT,
			 const grid* leakage, const grid* Cleakage,
			 const vectCells *OvlndOut, const vectCells *COvlndOut,
			 const vectCells *GWOut, const vectCells *CGWOut,
			 const grid* ttarea,			  
			 const Basin *b)
{
  d18OOut = AccountTrckOut(evapS, CevapS,evapI, CevapI, evapT, CevapT,
			   leakage, Cleakage, OvlndOut, COvlndOut, GWOut, CGWOut, ttarea, b);
}

void Budget::InstEvaporation_d18O(const grid* evapS, const grid* CevapS,
				 const grid* evapI, const grid* CevapI,
				 const grid* evapT, const grid* CevapT,
				 const Basin *b)
{
  d18OET = AccountTrckET(evapS, CevapS,evapI, CevapI, evapT, CevapT, b);
}

// For BasinAgeSummary.txt -------------------------------------------------------------
void Budget::InstOut_Age(const grid* evapS, const grid* CevapS,
			 const grid* evapI, const grid* CevapI,
			 const grid* evapT, const grid* CevapT,
			 const grid* leakage, const grid* Cleakage,
			 const vectCells *OvlndOut, const vectCells *COvlndOut,
			 const vectCells *GWOut, const vectCells *CGWOut,
			 const grid* ttarea,			 
			 const Basin *b)
{
  // In days
  AgeOut = AccountTrckOut(evapS, CevapS,evapI, CevapI, evapT, CevapT,
			  leakage, Cleakage, OvlndOut, COvlndOut, GWOut, CGWOut, ttarea, b);
}

void Budget::InstEvaporation_Age(const grid* evapS, const grid* CevapS,
				 const grid* evapI, const grid* CevapI,
				 const grid* evapT, const grid* CevapT,
				 const Basin *b)
{
  AgeET = AccountTrckET(evapS, CevapS,evapI, CevapI, evapT, CevapT, b);
}
