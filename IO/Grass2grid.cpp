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
 * Grass2grid.cpp
 *
 *  Created on: Oct 13, 2009
 *      Author: Marco Maneta
 */

#include <fstream>
#include <cmath>
#include "Grid.h"

int Grass2Grid(std::string fname, grid *pnt){

  ifstream input;
  string tags;
  REAL8 values;

  try{
    input.open(fname.c_str());
    if(!input.is_open()){
      cerr << "Couldn't open file " << fname << endl;
      exit(EXIT_FAILURE);
    }
  }catch(const exception& e){
    cerr << "Couldn't open file " << fname << " with message " << e.what() << endl;
    exit(EXIT_FAILURE);
  }

  for(int i=0; i < 6; i++)
    {
      input >> tags;
      input >> values;
      
      if(tags.compare("north:")==0)
	pnt->north = values;
      else if(tags.compare("south:")==0)
	pnt->south = values;
      else if(tags.compare("east:")==0)
	pnt->east = values;
      else if(tags.compare("west:")==0)
	pnt->west = values;
      else if(tags.compare("cols:")==0)
	pnt->c = (int)values;
      else if(tags.compare("rows:")==0)
	pnt->r = (int)values;
      else{
	cerr << "Check the format of file " << fname << endl;
	input.close();
	exit(EXIT_FAILURE);
      }
    }
  pnt->dx = (abs(pnt->east - pnt->west))/pnt->c;
	  
  //initializes the grid matrix
  try{
    pnt->matrix = new REAL8*[pnt->r];
    for (UINT4 i=0; i < pnt->r; i++)
      pnt->matrix[i] = new REAL8[pnt->c];
  }catch(std::bad_alloc&){cerr << "unable to allocate memory...";
    cin.get(); exit(-1);}


  //reads from the file and fills-in the matrix
  for(UINT4 i = 0; i < pnt->r; i++)
    {
      for(UINT4 j = 0; j < pnt->c; j++)
	{
	  input >> values;
	  if(input.eof())
	    {
	      cerr << "Unexpected end of file " << fname << endl;
	      exit(EXIT_FAILURE);
	    }
	  pnt->matrix[i][j] = values;
	}
    }
  
  try{
    if(input.is_open())
      input.close();

  }catch(const exception& e){
    cerr << "Couldn't close file with message " << e.what() << endl;
    exit(EXIT_FAILURE);
  }
  return EXIT_SUCCESS;
}
