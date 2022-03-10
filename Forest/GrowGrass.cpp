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
 * GrowGrass.cpp
 *
 *  Created on: Nov 11, 2012
 *      Author: marco
 *
 *      This is the grass growth module as per
 *      Istanbulloglu, Wang and Wedin (2012). Evaluation of
 *      ecohydrologic model parsimonry at local and regional scales in a
 *      semiarid grassland ecosystem. Ecohydrology, 5
 *
 *      We only use the carbon allocation and senescence expressions.
 *      The calculation of GPP and NPP is the same as for trees
 *
 */


#include "Forest.h"

int Forest::GrowGrass(UINT4 j, UINT4 r, UINT4 c, REAL8 dt){

  REAL8 lai_g, lai_d, lai_max;
  REAL8 nr, nf;

  lai_g = _species[j]._grassLAI_g->matrix[r][c];
  lai_d = _species[j]._grassLAI_d->matrix[r][c];
  lai_max = 5; //assumes lai_max is 5 as per Istanbulloglu (2012) eq 15a

  nf = 1 - ( lai_g / (lai_max - lai_d) );
  nr = 1 - nf;

  //Increase root Mass
  _species[j]._Del_RootMass->matrix[r][c] = _species[j]._NPP->matrix[r][c] * nr;
  //IncraseFoliageMass
  _species[j]._Del_FoliageMass->matrix[r][c] = _species[j]._NPP->matrix[r][c] * nf;

  GrowGrassLAI(j, r, c, dt); //Increase LAI of species j in r,c cell
  GrowRoots(j, r, c, dt); //Increase root density of species j in r,c cell

  return EXIT_SUCCESS;
}


