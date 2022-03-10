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
 *    Sylvain Kuppel, Xiaoqiang Yang
 *******************************************************************************/
/*
 * AccountTrckFluxes.cpp
 *
 *  Created on: Feb 28, 2018
 *      Author: Sylvain Kuppel, Aaron Smith
 */
#include "Budget.h"

double Budget::AccountTrckFluxes(const grid *map1, const grid *map2, const grid *map3, const Basin *b) {
  
  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  //  REAL8 dx = b->getCellSize();
  
#pragma omp parallel for		     \
  default(shared) private(r,c)		     \
  reduction (+:result)
  
  for (UINT4 i = 0; i < length; i++) {
    
    r = b->getSortedGrid().cells[i].row;
    c = b->getSortedGrid().cells[i].col;
    
    result +=  (map1->matrix[r][c] * map2->matrix[r][c] *dx*dx*dt*map3->matrix[r][c] );
    
  }
  
  return result;
}

// incoming surface water
double Budget::AccountTrckBCFluxesQ(const grid *map1, const grid *map2, const grid *map3, const Basin *b) {
  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = b->getCellSize();
  #pragma omp parallel for		     \
  default(shared) private(r,c)		     \
  reduction (+:result)
  
  for (UINT4 i = 0; i < length; i++) {
    
    r = b->getSortedGrid().cells[i].row;
    c = b->getSortedGrid().cells[i].col;
    // map 1 is already in cms
    result +=  (map1->matrix[r][c] * map2->matrix[r][c] *dt*map3->matrix[r][c] );
    
  }
  
  return result;
}

// incoming surface water
double Budget::AccountTrckBCFluxes(const grid *map1, const grid *map2, const grid *map3, const Basin *b) {
  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = b->getCellSize();
  #pragma omp parallel for		     \
  default(shared) private(r,c)		     \
  reduction (+:result)
  
  for (UINT4 i = 0; i < length; i++) {
    
    r = b->getSortedGrid().cells[i].row;
    c = b->getSortedGrid().cells[i].col;
    // map 1 is already in m2/s
    result +=  (map1->matrix[r][c] * map2->matrix[r][c] *dx*dt*map3->matrix[r][c] );
    
  }
  
  return result;
}

// Precip isotopes
double Budget::AccountTrckFluxes(const grid *map1, const grid *map2, const grid *map3, const Atmosphere *b) {
  
  UINT4 zones = b->getSortedGrid().size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = b->getCellSize();
  
#pragma omp parallel for		     \
  default(shared) private(r,c)		     \
  reduction (+:result)
  
  for (UINT4 i = 0; i < zones; i++)
    for (UINT4 j = 0; j < b->getSortedGrid()[i].cells.size(); j++) {
      
      r = b->getSortedGrid()[i].cells[j].row;
      c = b->getSortedGrid()[i].cells[j].col;
      
      result += (map1->matrix[r][c] * map2->matrix[r][c] * dx * dx * dt * map3->matrix[r][c]);
      //result += (map1->matrix[r][c] * Delta2Ratio(map2->matrix[r][c], iso) * dx * dx * dt);
    }
  
  return result;
}

// Precip age: it's 0 at entry, but the budgets must account for aging previously-input precip!
double Budget::AccountTrckFluxes(const grid *map1, const grid *map2, const Atmosphere *b){//, const Control *ctrl) {
  
  UINT4 zones = b->getSortedGrid().size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = b->getCellSize();
  
#pragma omp parallel for		     \
  default(shared) private(r,c)		     \
  reduction (+:result)
  
  for (UINT4 i = 0; i < zones; i++)
    for (UINT4 j = 0; j < b->getSortedGrid()[i].cells.size(); j++) {
      
      r = b->getSortedGrid()[i].cells[j].row;
      c = b->getSortedGrid()[i].cells[j].col;
      
      result += map1->matrix[r][c] * map2->matrix[r][c] * dx * dx * dt * dt / 86400; 
      // new input: has age 1 day at end-of-step
    }
  
  return result;
}


double Budget::AccountTrckFluxes(const vectCells *timeseries1, const vectCells *timeseries2) {
  
  UINT4 length = timeseries1->cells.size(); //b->getSortedGrid().cells.size();
  
  REAL8 result = 0;
  
#pragma omp parallel for			\
  reduction (+:result)
  
  for (UINT4 i = 0; i < length; i++) {
    
    //if(iso==0 or iso==1)
      //  result += timeseries1->cells[i].val * Delta2Ratio(timeseries2->cells[i].val, iso) * dt;
    //else if(iso ==2)
      //result +=  (timeseries1->cells[i].val * (timeseries2->cells[i].val+dt/86400) * dt);
      result +=  (timeseries1->cells[i].val * timeseries2->cells[i].val * dt);
  }
  
  return result;
}
