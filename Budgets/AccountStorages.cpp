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
 * AccountStorages.cpp
 *
 *  Created on: Mar 22, 2010
 *      Author: Marco Maneta
 */
#include "Budget.h"

double Budget::AccountStorages(const grid *map, const Basin *b)
{

  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = b->getCellSize();

#pragma omp parallel for			\
  default(shared) private(r,c)			\
  reduction (+:result)
  for (UINT4 i = 0; i< length; i++){

    r = b->getSortedGrid().cells[i].row;
    c = b->getSortedGrid().cells[i].col;

    result += (map->matrix[r][c]*dx*dx);
  }

  return result;
}

double Budget::AccountStorages(const grid *map1, const grid *map2, const Basin *b)
{

  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = b->getCellSize();

#pragma omp parallel for			\
  default(shared) private(r,c)			\
  reduction (+:result)
  for (UINT4 i = 0; i< length; i++){

    r = b->getSortedGrid().cells[i].row;
    c = b->getSortedGrid().cells[i].col;

    result += (map1->matrix[r][c]*map2->matrix[r][c]*dx*dx);
  }

  return result;
}

double Budget::AccountStorages(const grid *map1, const grid *map2, const grid *map3, const Basin *b)
{

  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = b->getCellSize();

#pragma omp parallel for			\
  default(shared) private(r,c)			\
  reduction (+:result)
  for (UINT4 i = 0; i< length; i++){

    r = b->getSortedGrid().cells[i].row;
    c = b->getSortedGrid().cells[i].col;

    result += (map1->matrix[r][c]*map2->matrix[r][c]*dx*dx*map3->matrix[r][c]);
  }

  return result;
}

// --- Tracking ------------------------------------------------------------------------

// -- Sum tracer * storage
/*
double Budget::AccountTrckStorages(const grid *map1, const grid *map2, const grid *map3,const Basin *b)
{
  
  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 result = 0;
  REAL8 dx = b->getCellSize();
  
#pragma omp parallel for			\
  default(shared) private(r,c)			\
  reduction (+:result)
  for (UINT4 i = 0; i< length; i++){
    
    r = b->getSortedGrid().cells[i].row;
    c = b->getSortedGrid().cells[i].col;

    result += (map1->matrix[r][c]* map2->matrix[r][c] *dx*dx*map3->matrix[r][c]);

  }
  
  return result;
}
*/
// -- Sum tracer * storage / storage
double Budget::AccountTrckStorages2(const grid *map1, const grid *map2, const Basin *b)
{
  
  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 numer = 0;
  REAL8 denom = 0;
  REAL8 result = 0;
  
#pragma omp parallel default(shared) private(r,c)
  {
    
#pragma omp for reduction (+:numer, denom)
    for (UINT4 i = 0; i< length; i++){
      
      r = b->getSortedGrid().cells[i].row;
      c = b->getSortedGrid().cells[i].col;
    
      numer += map1->matrix[r][c]* map2->matrix[r][c];
      denom += map1->matrix[r][c];
    
    }
  }
  result = numer / denom;
  
  return result;
}


// -- Volume-weighted average across the entire domain
double Budget::AccountTrckDomain(const grid *mapCanopy, const grid *mapCCanopy, 
				 const grid *mapSnow, const grid *mapCSnow, 
				 const grid *mapSurface, const grid *mapCSurface, 
				 const grid *mapL1, const grid *mapCL1, 
				 const grid *mapL2, const grid *mapCL2, 
				 const grid *mapL3, const grid *mapCL3, 
				 const grid *mapGW, const grid *mapCGW, 
				 const Basin *b)
{
  
  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 numer = 0;
  REAL8 denom = 0;
  REAL8 result = 0;
  
  
#pragma omp parallel default(shared) private(r,c)
{
#pragma omp for reduction (+:numer, denom)
    for (UINT4 i = 0; i< length; i++){
    
    r = b->getSortedGrid().cells[i].row;
    c = b->getSortedGrid().cells[i].col;

    numer += mapCanopy->matrix[r][c]* mapCCanopy->matrix[r][c] +
      mapSnow->matrix[r][c]* mapCSnow->matrix[r][c] +
      mapSurface->matrix[r][c]* mapCSurface->matrix[r][c] +
      mapL1->matrix[r][c]* mapCL1->matrix[r][c] +
      mapL2->matrix[r][c]* mapCL2->matrix[r][c] +
      mapL3->matrix[r][c]* mapCL3->matrix[r][c] +
      mapGW->matrix[r][c]* mapCGW->matrix[r][c] ;
    
    denom += mapCanopy->matrix[r][c] + mapSnow->matrix[r][c] + mapSurface->matrix[r][c] + 
      mapL1->matrix[r][c] + mapL2->matrix[r][c] + 
      mapL3->matrix[r][c] + mapGW->matrix[r][c] ;
    
  }
    }
  //}
  result = numer / denom ;

  return result;
}

// -- Volume-weighted average across the soil domain (including GW)
double Budget::AccountTrckRootZone(const grid *mapL1, const grid *mapCL1, const grid *mappL1,
				   const grid *mapL2, const grid *mapCL2, const grid *mappL2,
				   const grid *mapL3, const grid *mapCL3, const grid *mappL3,
				   const grid *mapGW, const grid *mapCGW, 
				   const Basin *b)
{
  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 numer = 0;
  REAL8 denom = 0;
  REAL8 result = 0;
  
#pragma omp parallel default(shared) private(r,c)
  {
    
#pragma omp for reduction (+:numer, denom)
    for (UINT4 i = 0; i< length; i++){
    
      r = b->getSortedGrid().cells[i].row;
      c = b->getSortedGrid().cells[i].col;

      numer += mapL1->matrix[r][c] * mapCL1->matrix[r][c] * mappL1->matrix[r][c] +
	mapL2->matrix[r][c]* mapCL2->matrix[r][c] * mappL2->matrix[r][c] +
	(mapL3->matrix[r][c]* mapCL3->matrix[r][c] +
	 mapGW->matrix[r][c]* mapCGW->matrix[r][c]) * mappL3->matrix[r][c] ;
      
      denom += mapL1->matrix[r][c] * mappL1->matrix[r][c] + 
	mapL2->matrix[r][c]* mappL2->matrix[r][c] + 
	(mapL3->matrix[r][c] + mapGW->matrix[r][c])* mappL3->matrix[r][c] ;
    }
  }
  result = numer / denom ;
  
  return result;
}


// -- Volume-weighted average across the soil domain (including GW)
double Budget::AccountTrckVadose(const grid *mapL1, const grid *mapCL1, 
				 const grid *mapL2, const grid *mapCL2, 
				 const grid *mapL3, const grid *mapCL3, 
				 const grid *mapGW, const grid *mapCGW, 
				 const Basin *b)
{
  UINT4 length = b->getSortedGrid().cells.size();
  UINT4 r, c;
  REAL8 numer = 0;
  REAL8 denom = 0;
  REAL8 result = 0;
  
#pragma omp parallel default(shared) private(r,c)
  {
    
#pragma omp for reduction (+:numer, denom)
    for (UINT4 i = 0; i< length; i++){
    
      r = b->getSortedGrid().cells[i].row;
      c = b->getSortedGrid().cells[i].col;

      numer += mapL1->matrix[r][c]* mapCL1->matrix[r][c] +
	mapL2->matrix[r][c]* mapCL2->matrix[r][c] +
	mapL3->matrix[r][c]* mapCL3->matrix[r][c] +
	mapGW->matrix[r][c]* mapCGW->matrix[r][c] ;
    
      denom += mapL1->matrix[r][c] + mapL2->matrix[r][c] + 
	mapL3->matrix[r][c] + mapGW->matrix[r][c] ;
    }
  }
  result = numer / denom ;
  
  return result;
}
