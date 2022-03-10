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
 * AccountInputFluxes.cpp
 *
 *  Created on: Mar 8, 2010
 *      Author: Marco Maneta
 */
#include "Budget.h"

double Budget::AccountFluxes(const grid *map1, const grid *map2, const Basin *b) {

  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = b->getCellSize();

#pragma omp parallel for			\
  default(shared) private(r,c)			\
  reduction (+:result)

  for (UINT4 i = 0; i < length; i++) {

    r = b->getSortedGrid().cells[i].row;
    c = b->getSortedGrid().cells[i].col;

    result += (map1->matrix[r][c] * dx * dx * dt * map2->matrix[r][c]);
    //map2 is the proportion of cell within the catchment
  }

  return result;
}

double Budget::AccountFluxes(const grid *map1, const grid *map2, const Atmosphere *b) {

  UINT4 zones = b->getSortedGrid().size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = b->getCellSize();

#pragma omp parallel for			\
  default(shared) private(r,c)			\
  reduction (+:result)

  for (UINT4 i = 0; i < zones; i++)
    for (UINT4 j = 0; j < b->getSortedGrid()[i].cells.size(); j++) {

      r = b->getSortedGrid()[i].cells[j].row;
      c = b->getSortedGrid()[i].cells[j].col;

      result += (map1->matrix[r][c] * dx * dx * dt * map2->matrix[r][c]);
      //map2 is the proportion of cell within the catchment
    }

  return result;
}

double Budget::AccountFluxes(const vectCells *timeseries, const Basin *b) {

  UINT4 length = timeseries->cells.size(); //b->getSortedGrid().cells.size();

  REAL8 result = 0;

#pragma omp parallel for			\
  reduction (+:result)

  for (UINT4 i = 0; i < length; i++) {

    result += timeseries->cells[i].val * dt;
  }

  return result;
}

// --- Tracking : uses two maps (or vectors) to multiply -----------------------------------

double Budget::AccountTrckFluxes(const grid *map1, const grid *map2, const grid *map3, const Basin *b) {
  
  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = b->getCellSize();
  
#pragma omp parallel for			\
  default(shared) private(r,c)			\
  reduction (+:result)
  
  for (UINT4 i = 0; i < length; i++) {
    
    r = b->getSortedGrid().cells[i].row;
    c = b->getSortedGrid().cells[i].col;
    
    result +=  (map1->matrix[r][c] * map2->matrix[r][c] *dx*dx*dt*map3->matrix[r][c]);
    //map3 is the proportion of cell within the catchment
  }
  
  return result;
}

// Precip isotopes
//map3:ttarea
double Budget::AccountTrckFluxes(const grid *map1, const grid *map2, const grid *map3, const Atmosphere *a) {
  
  UINT4 zones = a->getSortedGrid().size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = a->getCellSize();
  
#pragma omp parallel for			\
  default(shared) private(r,c)			\
  reduction (+:result)
  
  for (UINT4 i = 0; i < zones; i++)
    for (UINT4 j = 0; j < a->getSortedGrid()[i].cells.size(); j++) {
      
      r = a->getSortedGrid()[i].cells[j].row;
      c = a->getSortedGrid()[i].cells[j].col;
      
      result += (map1->matrix[r][c] * map2->matrix[r][c] * dx * dx * dt * map3->matrix[r][c]);
    }
  
  return result;
}

// Precip age: it's 0 at entry, but the budgets must account for aging previously-input precip!
//map2 = ttarea
double Budget::AccountTrckFluxes(const grid *map1, const grid *map2, const Atmosphere *a){
  
  UINT4 zones = a->getSortedGrid().size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = a->getCellSize();
  
#pragma omp parallel for			\
  default(shared) private(r,c)			\
  reduction (+:result)
  
  for (UINT4 i = 0; i < zones; i++)
    for (UINT4 j = 0; j < a->getSortedGrid()[i].cells.size(); j++) {
      
      r = a->getSortedGrid()[i].cells[j].row;
      c = a->getSortedGrid()[i].cells[j].col;
      
      result += map1->matrix[r][c] * map2->matrix[r][c]* dx * dx * dt * dt / 86400; 
    }
  
  return result;
}

double Budget::AccountTrckFluxes(const vectCells *timeseries1, const vectCells *timeseries2) {
  
  UINT4 length = timeseries1->cells.size(); //b->getSortedGrid().cells.size();
  REAL8 result = 0;
  
#pragma omp parallel for default(shared)	\
  reduction (+:result)
  
  for (UINT4 i = 0; i < length; i++)
    result +=  (timeseries1->cells[i].val * timeseries2->cells[i].val * dt);
  
  return result;
}

// == AgeReporting stuff ------------------------------------------------------------------------------

double Budget::AccountTrckFluxes2(const grid *map1, const grid *map2, const Basin *b) {
  
  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 numer = 0;
  REAL8 denom = 0;
  
#pragma omp parallel default(shared) private(r,c)	\
  
  {
#pragma omp for reduction (+:numer, denom)
  
    for (UINT4 i = 0; i < length; i++) {
    
      r = b->getSortedGrid().cells[i].row;
      c = b->getSortedGrid().cells[i].col;
    
      numer += map1->matrix[r][c] * map2->matrix[r][c];
      denom += map1->matrix[r][c];
    
    }
  }
  result = numer / denom;
  
  return result;
}

double Budget::AccountTrckFluxes2(const vectCells *timeseries1, const vectCells *timeseries2) {
  
  UINT4 length = timeseries1->cells.size(); //b->getSortedGrid().cells.size();
  REAL8 result = 0;
  REAL8 numer = 0;
  REAL8 denom = 0;
  
  //cout << length << endl;
  
#pragma omp parallel default(shared)	
  {

#pragma omp for reduction (+:numer, denom)
    
    for (UINT4 i = 0; i < length; i++){
      numer +=  timeseries1->cells[i].val * timeseries2->cells[i].val; 
      //cout << timeseries1->cells[i].val << endl;
      //cout << timeseries2->cells[i].val << endl;
      denom +=  timeseries1->cells[i].val;
    }
  }
  result = numer / denom;
  return result;
}

// -- Flux-weighted average across evaporative losses
double Budget::AccountTrckET(const grid* evapS, const grid* CevapS,
			     const grid* evapI, const grid* CevapI,
			     const grid* evapT, const grid* CevapT,
			     const Basin *b)
{
  
  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 numer = 0;
  REAL8 denom = 0;
  
#pragma omp parallel default(shared) private(r,c) 
  
  {
#pragma omp for reduction (+:numer, denom)
    for (UINT4 i = 0; i< length; i++){
      
      r = b->getSortedGrid().cells[i].row;
      c = b->getSortedGrid().cells[i].col;
      
      numer += evapS->matrix[r][c]* CevapS->matrix[r][c] +
	evapI->matrix[r][c]* CevapI->matrix[r][c] +
	evapT->matrix[r][c]* CevapT->matrix[r][c] ;
      
      denom += evapS->matrix[r][c] + evapI->matrix[r][c] + evapT->matrix[r][c] ;
    }
  }   
  
  result = numer / denom ;
  return result;
}

// -- Flux-weighted average across outputs
double Budget::AccountTrckOut(const grid* evapS, const grid* CevapS,
			      const grid* evapI, const grid* CevapI,
			      const grid* evapT, const grid* CevapT,
			      const grid* leakage, const grid* Cleakage,
			      const vectCells *OvlndOut, const vectCells *COvlndOut,
			      const vectCells *GWOut, const vectCells *CGWOut,
                              const grid* ttarea,
			      const Basin *b)
{
  
  UINT4 length1 = b->getSortedGrid().cells.size();
  UINT4 length2 = OvlndOut->cells.size(); //b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 numer1 = 0;
  REAL8 numer2 = 0;
  REAL8 denom1 = 0;
  REAL8 denom2 = 0;
  REAL8 dx = b->getCellSize();
  
#pragma omp parallel default(shared) private(r,c) 
  
  {
#pragma omp for reduction (+:numer1, denom1)

    for (UINT4 i = 0; i< length1; i++){
	
      r = b->getSortedGrid().cells[i].row;
      c = b->getSortedGrid().cells[i].col;
      
      numer1 += (evapS->matrix[r][c]* CevapS->matrix[r][c] +
		 evapI->matrix[r][c]* CevapI->matrix[r][c] +
		 evapT->matrix[r][c]* CevapT->matrix[r][c] +
		 leakage->matrix[r][c]* Cleakage->matrix[r][c]) *dx*dx * ttarea->matrix[r][c];
      
      denom1 += (evapS->matrix[r][c] + evapI->matrix[r][c] + evapT->matrix[r][c] +
		 leakage->matrix[r][c]) *dx*dx  * ttarea->matrix[r][c];
    }
 
    //cout << numer1 << " " << denom1 << endl;
    //cout << numer2 << " " << denom2 << endl;
  
#pragma omp for reduction (+:numer2, denom2)

    for (UINT4 j = 0; j < length2; j++) {
      numer2 +=  OvlndOut->cells[j].val * COvlndOut->cells[j].val +
	GWOut->cells[j].val * CGWOut->cells[j].val;
      
      denom2 +=  OvlndOut->cells[j].val + GWOut->cells[j].val;
      
    }
  }
  //cout << numer1 << " " << denom1 << endl;
  //cout << numer2 << " " << denom2 << endl;
  
  result = (numer1 + numer2) / (denom1 + denom2);

  //cout << result << endl;

  return result;
}
