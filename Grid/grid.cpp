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
 * grid.cpp
 *
 *  Created on: May 21, 2009
 *      Author: Marco Maneta
 *
 *      grid.cpp, grid.h
 *      implementation of grid constructor, destructor and operator overloads
 */

#include <fstream>
#include "Grid.h"

int Grass2Grid(std::string fname, grid *pnt);
int CSF2Grid(std::string fname, grid *pnt);
int Table2Grid(std::string fname, grid *pnt);

grid::grid(UINT4 row, UINT4 cols): r(row), c(cols)
{
  try{
    matrix = new REAL8*[row];
    for (UINT4 i=0; i < row; i++)
      matrix[i] = new REAL8[cols];
  }catch(std::bad_alloc &){cerr << "unable to allocate memory...";
    cin.get(); exit(-1);}

#pragma omp parallel for
  for (UINT4 i = 0; i<r; i++)
    for (UINT4 j = 0; j<c; j++)
      matrix[i][j] = 0;

  north = south = east = west = nodata = dx = 0;

}

grid::grid(std::string fname, std::string type){

  if(!type.compare("grass"))
    Grass2Grid(fname, this);
  else if (!type.compare("csf"))
    CSF2Grid(fname, this);
  else if (!type.compare("table"))
    Table2Grid(fname, this);
  else{
    cerr << "Illegal type " << type << endl;
    exit(EXIT_FAILURE);
  }


}

grid::grid(const grid &m)
{
  r = m.r;
  c = m.c;
  try{
    matrix = new REAL8*[r];
    for (UINT4 i=0; i < r; i++)
      matrix[i] = new REAL8[c];
  }catch(std::bad_alloc &){
    cerr << "unable to allocate memory...";
    cin.get();
    exit(-1);
  }
  for (UINT4 i = 0; i<r; i++)
    for (UINT4 j = 0; j<c; j++)
      matrix[i][j] = 0.0;
  //matrix[i][j] = m.matrix[i][j];

  north = m.north;
  south = m.south;
  east = m.east;
  west = m.west;
  nodata = m.nodata;
  dx = m.dx;

}


void grid::reset(){

  for (UINT4 i=0; i<r; i++)
    for (UINT4 j=0; j<c; j++)
      this->matrix[i][j]=0;
}

/*
void grid::resetIso(){

  for (UINT4 i=0; i<r; i++)
    for (UINT4 j=0; j<c; j++)
      this->matrix[i][j]= -1000;
}
*/

grid& grid::operator=(const grid &m)
{
  if(this == &m) return *this;

  for (UINT4 i = 0; i < r; i++)
    delete[] matrix[i];
  delete[] matrix;
  r = m.r;
  c = m.c;
  north = m.north;
  south = m.south;
  east = m.east;
  west = m.west;
  nodata = m.nodata;
  dx = m.dx;

  try{
    matrix = new REAL8*[r];
    for (UINT4 i=0; i < r; i++)
      matrix[i] = new REAL8[c];
  }catch(std::bad_alloc &){
    std::cerr << "unable to allocate memory...";
    cin.get();
    exit(-1);
  }
  for (UINT4 i = 0; i<r; i++)
    for (UINT4 j = 0; j<c; j++)
      matrix[i][j] = m.matrix[i][j];
  return *this;
}


grid grid::operator+(const grid &m)  //check for memory leaks here
{
  UINT4 r = m.r;
  UINT4 c = m.c;
  struct grid result(r,c) ;
  for (UINT4 i=0; i<r; i++)
    for (UINT4 j=0; j<c; j++)
      result.matrix[i][j] =
	this->matrix[i][j] +
	m.matrix[i][j];
  return result;
}


grid grid::operator-(const grid &m)
{
  UINT4 r = m.r;
  UINT4 c = m.c;
  struct grid result(r,c) ;
  for (UINT4 i=0; i<r; i++)
    for (UINT4 j=0; j<c; j++)
      result.matrix[i][j] =
	this->matrix[i][j] -
	m.matrix[i][j];
  return result;
}

grid& grid::operator+=(const grid &m)  //and here
{
  UINT4 r = m.r;
  UINT4 c = m.c;
  for (UINT4 i=0; i<r; i++)
    for (UINT4 j=0; j<c; j++)
      matrix[i][j] +=m.matrix[i][j];
  return *this;
}

grid& grid::operator-=(const grid &m)
{
  UINT4 r = m.r;
  UINT4 c = m.c;
  for (UINT4 i=0; i<r; i++)
    for (UINT4 j=0; j<c; j++)
      matrix[i][j] -=m.matrix[i][j];
  return *this;
}

bool grid::operator==(const grid&b) const
{

  if(r!=b.r || c!=b.c)
    return false;

  for (UINT4 i=0; i<r; i++)
    for (UINT4 j=0; j<c; j++){
      if(matrix[i][j]!=b.matrix[i][j])
	return false;
    }
  return true;
}


bool grid::operator!=(const grid&b) const{
  return !(*this == b);

}

grid::~grid()
{
  for (UINT4 i = 0; i < r; i++)
    delete[] matrix[i];
  delete[] matrix;
}
