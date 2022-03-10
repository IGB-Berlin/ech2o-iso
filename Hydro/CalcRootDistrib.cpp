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
 * CalcRootDistrib.cpp
 *
 *  Created on: Feb 6th, 2017
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"


int Basin::CalcRootDistrib(){

  UINT4 r, c;
  UINT4 s, nsp;
  REAL8 frac1, frac2;
  REAL8 k, d1, d2, d, d95;
  REAL8 p_veg = 0;
  REAL8 p = 0;
  
  nsp = fForest->getNumSpecies();

#pragma omp parallel default(none)		\
  private(s,r,c,k,d,d1,d2,d95,frac1,frac2,p,p_veg)	\
  shared(nsp)
  { 
#pragma omp for nowait

    for (UINT4 j = 0; j < _vSortedGrid.cells.size(); j++) {
	r = _vSortedGrid.cells[j].row;
	c = _vSortedGrid.cells[j].col;

	d = _soildepth->matrix[r][c];
	d1 = _depth_layer1->matrix[r][c];
	d2 = _depth_layer2->matrix[r][c];
	d95 = 0; 

	p_veg = 0;
	
	for (s = 0; s < nsp; s++) {
	  if (s == nsp - 1) { //if this is bare ground set fracs to 0
	    frac1 = 0;
	    frac2 = 0;
	  } 
	  else {
	    // use exponential profile
	    k = fForest->getKRoot(s);
	    frac1 = (1 - expl(-k*d1))/(1-expl(-k*d));
	    frac2 = (expl(-k*d1) - expl(-k*(d1+d2)))/(1-expl(-k*d));

	    // Contribution of each layer to rootzone: use depth at which 95% of the roots
	    // are found
	    d95 = std::min<double>(d,log2l(0.05+0.95*expl(-k*d))/(-k));
	    // average over species is made using the vegetated fraction sum (p_veg<=1)
	    p = fForest->getPropSpecies(s, r, c);

	    _ProotzoneL1->matrix[r][c] += d95 >= d1 ? p : p*d95/d1 ;
	    _ProotzoneL2->matrix[r][c] += d95 >= d1+d2 ? p : std::max<double>(0.0,d95-d1)*p/d2 ;
	    _ProotzoneL3->matrix[r][c] += d95 >= d ? p : std::max<double>(0.0,d95-d1-d2)*p/(d-d1-d2);

	    // Average over species fraction
	    _Zroot95->matrix[r][c] += p*d95 ;
	    
	    p_veg += p;

	  }

	  fForest->setRootFrac1Species(s, r, c, frac1);
	  fForest->setRootFrac2Species(s, r, c, frac2);
	}

	// Average contribution of each layer to root zone storage
	// (using p_veg allows to discard 100% bare soil pixels in catchment budgets)
	_ProotzoneL1->matrix[r][c] = p_veg > RNDOFFERR ? _ProotzoneL1->matrix[r][c] / p_veg : 0;
	_ProotzoneL2->matrix[r][c] = p_veg > RNDOFFERR ? _ProotzoneL2->matrix[r][c] / p_veg : 0;
	_ProotzoneL3->matrix[r][c] = p_veg > RNDOFFERR ? _ProotzoneL3->matrix[r][c] / p_veg : 0;
	_Zroot95->matrix[r][c] = p_veg > RNDOFFERR ? _Zroot95->matrix[r][c] / p_veg : 0;

      } // for
  } //end omp parallel block
  return EXIT_SUCCESS;

}
