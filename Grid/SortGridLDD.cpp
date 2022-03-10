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
 * SortGridLDD.cpp
 *
 *  Created on: Oct 17, 2009
 *      Author: Marco Maneta
 */

#include "Basin.h"
#include "ConstAndFuncs.h"

vectCells Basin::SortGridLDD(){

		REAL4 value = 0;
	    REAL4 value1 = 0;
	    REAL4 value2 = 0;
	    REAL4 value3 = 0;
	    REAL4 value4 = 0;
	    REAL4 value6 = 0;
	    REAL4 value7 = 0;
	    REAL4 value8 = 0;
	    REAL4 value9 = 0;

	    //	    REAL4 flag;
	    REAL8 checked = _DEM->nodata;

	    grid *temp;
	    UINT4 r,c;
	    UINT4 nr = _DEM->r;
	    UINT4 nc = _DEM->c;

	    temp = new grid(nr, nc); //init grid

	    vectCells map2array; //instantiate the array to hold the sorted map

	    ////////////////////////////////////
	    //////makes a temporary copy of ldd.
	    /////////////////////////////////////
	    *temp = *_ldd;
	    ///////////////////////////////

	    //calculates how many cells are in the ldd area. Needed for stop condition
	       UINT4 counter = 0;

           #pragma omp parallel for\
           default(shared) private(r,c, value) \
           reduction (+:counter)

	       for(r=1; r < nr-1; r++)
	       {
	        for(c=1; c < nc-1; c++)
	        {
	         value = (int)temp->matrix[r][c];
	         if (value > 0 && value <=9)
	             counter ++;
	        }
	       }

	    //this loop sorts the map filling in the array with the right order of the
	    //cells and creates a colored map with the incremental order of the cells.
	    //It stops when the array has the same number of ldd cells than the map
	    UINT4 x = 0;
	    UINT4 i = 0;
	    std::cout << "Sorting DEM drainage network... " << std::endl;

	    do
	    {
	    x++;
	      for(r=1; r < nr-1; r++)
	      {
	          for(c=1; c < nc-1; c++)
	          {
	          value = (int)temp->matrix[r][c];
	              if (value == checked)continue;
	          value7 = temp->matrix[r-1][c-1];
	          value8 = temp->matrix[r-1][c];
	          value9 = temp->matrix[r-1][c+1];
	          value4 = temp->matrix[r][c-1];
	          value6 = temp->matrix[r][c+1];
	          value1 = temp->matrix[r+1][c-1];
	          value2 = temp->matrix[r+1][c];
	          value3 = temp->matrix[r+1][c+1];
	              if (value7 != 3 &&
	                  value8 != 2 &&
	                  value9 != 1 &&
	                  value4 != 6 &&
	                  value6 != 4 &&
	                  value1 != 9 &&
	                  value2 != 8 &&
	                  value3 != 7 &&
	                  value != _DEM->nodata)
	              {
			//	               flag = x;

	               map2array.cells.push_back(cell(r,c,(int)value));
	              }
	           }
	      }
	      for(; i < map2array.cells.size(); i++)
	      {
	          r = map2array.cells[i].row;
	          c = map2array.cells[i].col;
	          temp->matrix[r][c] = (REAL4)checked;
	      }

	     printProgBar(map2array.cells.size() * 100 / counter);
	    }while(map2array.cells.size()< counter);

	    delete temp;
	    return map2array;
	    std::cout << endl;

}

void loadSortedGrid(vectCells &v, const char *filename)
{

	std::ifstream ifs(filename);
	boost::archive::text_iarchive ia(ifs);
	ia >> v;
	ifs.close();
}

void saveSortedGrid(vectCells &v, const char *filename){
	std::ofstream ofs(filename);
	boost::archive::text_oarchive oa(ofs);
	oa << v;
	ofs.close();
}

