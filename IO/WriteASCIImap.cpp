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
 * GridIO.cpp
 *
 *  Created on: May 23, 2009
 *      Author: Marco Maneta
 *
 *
 *      writeASCIIMap writes the grid file into an ASCII file with grass format
 *
 *
 */

#include <fstream>
#include "Grid.h"



int grid::writeASCIIMap(std::string fname)
{


	fstream output;

	try{

		  output.open(fname.c_str(), ios_base::out);
		        if(!output.is_open()){
		         cerr << "Couldn't open file " << fname << endl;
		         exit(EXIT_FAILURE);
		        }
		  }catch(const exception& e){

		    cerr << "Couldn't open file " << fname << " with message " << e.what() << endl;
		    exit(EXIT_FAILURE);
		 }


		  try{

	        output << "north: " << north << endl;
	        output << "south:  " << south << endl;
	        output << "east:  " << east << endl;
	        output << "west:  " << west << endl;
	        output << "rows: " << r << endl;
	        output << "cols: " << c << endl;

	        for(UINT4 i=0; i < r; i++)
	        {
	          for(UINT4 j=0; j < c; j++)
	          {
	            output << matrix[i][j] << " ";
	          }
	          output << endl;
	        }

	        output.close();

	     //   cout << "File "<< fname << " succesfully written." << endl;

	        }catch(ios_base::failure& e){

	         cerr << "Exception occured while reading file " << fname<<
	         " with message " << e.what() << endl;

	         exit(EXIT_FAILURE);
	        }

	  return EXIT_SUCCESS;
}
