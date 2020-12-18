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
 * OutletVals.cpp
 *
 *  Created on: Mar 1, 2018
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"

void Tracking::OutletVals(Control &ctrl, int mode, int r, int c)
{
    
  if (mode == 0){ // Initialization
    if(ctrl.sw_2H){
      _d2HOvlndOutput.cells.clear();
      _d2HL1wtrOutput.cells.clear();
      _d2HL2wtrOutput.cells.clear();
      _d2HGwtrOutput.cells.clear();
    }
    if(ctrl.sw_18O){
      _d18OOvlndOutput.cells.clear();
      _d18OL1wtrOutput.cells.clear();
      _d18OL2wtrOutput.cells.clear();
      _d18OGwtrOutput.cells.clear();
    }
    if(ctrl.sw_Age){
      _AgeOvlndOutput.cells.clear();
      _AgeL1wtrOutput.cells.clear();
      _AgeL2wtrOutput.cells.clear();
      _AgeGwtrOutput.cells.clear();
    }
  }

  else if (mode == 1){ // Active mode

    // Deuterium
    if(ctrl.sw_2H){
      _d2HL1wtrOutput.cells.push_back(cell(r, c, _d2Hsoil1_Out->matrix[r][c]));
      _d2HL2wtrOutput.cells.push_back(cell(r, c, _d2Hsoil2_Out->matrix[r][c]));
      _d2HGwtrOutput.cells.push_back(cell(r, c, _d2Hgroundwater->matrix[r][c]));
      _d2HOvlndOutput.cells.push_back(cell(r, c, _d2Hsurface->matrix[r][c])); 
    }
    // Oxygen 18
    if(ctrl.sw_18O){
      _d18OL1wtrOutput.cells.push_back(cell(r, c, _d18Osoil1_Out->matrix[r][c]));
      _d18OL2wtrOutput.cells.push_back(cell(r, c, _d18Osoil2_Out->matrix[r][c]));
      _d18OGwtrOutput.cells.push_back(cell(r, c, _d18Ogroundwater->matrix[r][c]));
      _d18OOvlndOutput.cells.push_back(cell(r, c, _d18Osurface->matrix[r][c])); 
    }
    // Age
    if(ctrl.sw_Age){
      _AgeL1wtrOutput.cells.push_back(cell(r, c, _Agesoil1_Out->matrix[r][c]));
      _AgeL2wtrOutput.cells.push_back(cell(r, c, _Agesoil2_Out->matrix[r][c]));
      _AgeGwtrOutput.cells.push_back(cell(r, c, _Agegroundwater->matrix[r][c]));
      _AgeOvlndOutput.cells.push_back(cell(r, c, _Agesurface->matrix[r][c])); 
    }
  }
}


