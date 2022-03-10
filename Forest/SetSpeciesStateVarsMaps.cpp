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
 * SetSpeciesStateVarsMaps.cpp
 *
 *  Created on: Sep 20, 2011
 *      Author: Marco.Maneta
 */

#include "Forest.h"

void Forest::SetStateVarsMaps(Control &ctrl){

  UINT4 r,c;
  stringstream fn;

  //for PCRASTER maps

  try{
    if (!ctrl.MapType.compare("csf")) {
      for (unsigned int j = 0; j < _Nsp - 1; j++) { //loads all species expect
	//loads the maps to their respective grids
	fn.str(""); fn << "p_" << j << ".map";
	if(_species[j]._fraction->PCRMap2grid(ctrl.path_BasinFolder + fn.str())==-1)
	  throw fn.str();
	
	if(ctrl.toggle_veg_dyn!=2){
	  fn.str(""); fn << "lai_" << j << ".map";
	  if(_species[j]._LAI->PCRMap2grid(ctrl.path_BasinFolder + fn.str())==-1)
	    throw fn.str();
	  fn.str(""); fn << "hgt_" << j << ".map";
	  if(_species[j]._Height->PCRMap2grid(ctrl.path_BasinFolder + fn.str())==-1)
	    throw fn.str();	  
	}

	*_species[j]._grassLAI_g = *_species[j]._LAI; //needed to initialize lai for grasses
	
	fn.str(""); fn << "age_" << j << ".map";
	if(_species[j]._AGE->PCRMap2grid(ctrl.path_BasinFolder + fn.str())==-1)
	  throw fn.str();
	fn.str(""); fn << "bas_" << j << ".map";
	if(_species[j]._BasalArea->PCRMap2grid(ctrl.path_BasinFolder + fn.str())==-1)
	  throw fn.str();
	fn.str(""); fn << "ntr_" << j << ".map";
	if(_species[j]._StemDensity->PCRMap2grid(ctrl.path_BasinFolder + fn.str())==-1)
	  throw fn.str();
	fn.str(""); fn << "root_" << j << ".map";
	if(_species[j]._RootMass->PCRMap2grid(ctrl.path_BasinFolder + fn.str())==-1)
	  throw fn.str();

      }
    }
    else{

      cerr << "Reading forest state vars from maps not implemented yet for ascii maps. Pleas, use csf (PCRASTER) maps" << endl;
      exit(EXIT_FAILURE);
    }
  }catch(string &e){
    cerr << "Initial forest map " << e << " not found" << endl;
    exit(EXIT_FAILURE);
  }

  //calculate proportion of of bare soil

  for (UINT4 k = 0; k < _vSortedGrid.cells.size() ; k++)
    {
      r = _vSortedGrid.cells[k].row;
      c = _vSortedGrid.cells[k].col;
      _species[_Nsp-1]._fraction->matrix[r][c] = 1;

      for(UINT4 j = 0; j < _Nsp - 1; j++)
	_species[_Nsp-1]._fraction->matrix[r][c] -= _species[j]._fraction->matrix[r][c];

      if(_species[_Nsp-1]._fraction->matrix[r][c] < 0){
	cerr << "Proportion of species is larger than 1 in cell row: " << r << " col: " << c << endl;
	exit(EXIT_FAILURE);
      }

    }



}
