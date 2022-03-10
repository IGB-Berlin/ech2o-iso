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
 * CheckMaps.cpp
 *
 *  Created on: Feb 10, 2011
 *      Author: Marco.Maneta
 */

#include "Basin.h"

void Basin::CheckMaps(Control &ctrl) {

  UINT4 r, c;
  UINT4 j = 0;
  bool excep_thrown = false; //  poor man way  to rethrow exception outside omp pragma
  UINT4 length = _vSortedGrid.cells.size();

#pragma omp parallel for			\
  default(shared) private(r,c,j)
  for (j = 0; j < length; j++) {
    r = _vSortedGrid.cells[j].row;
    c = _vSortedGrid.cells[j].col;
    try {

      _soildepth->matrix[r][c] += _depth_layer1->matrix[r][c] + _depth_layer2->matrix[r][c];

      if (_vSortedGrid.cells[j].dir < 1
	  || _vSortedGrid.cells[j].dir > 9) {
	string e("The drain direction map is not a valid ldd map...\n");
	throw e;
      }

      if (_slope->matrix[r][c] == _slope->nodata) {
	string e(
		 "Slope map contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_slope->matrix[r][c] == 0)
	_slope->matrix[r][c] = MIN_SLOPE;
      if (_slope->matrix[r][c] <= 0) {
	string e(
		 "Slope map contains negative or zero values inside the valid domain...\n");
	throw e;
      }

      if (_Ksat0->matrix[r][c] == _Ksat0->nodata) {
	string e(
		 "Ksat0 map contains no data values inside the valid domain...\n");
	throw e;
      }
      if(ctrl.sw_expKsat){
	if (_kKsat->matrix[r][c] == _kKsat->nodata) {
	  string e(
		   "Ksat profile map contains no data values inside the valid domain...\n");
	  throw e;
	}
      }

      if(_fImperv->matrix[r][c] < 0)
	throw string("Fraction of impervious area is negative inside the valid domain ...\n");
      if(_fImperv->matrix[r][c] > 1)
	throw string("Fraction of impervious area is larger than 1 inside the valid domain ...\n");
      
      if (_slope->matrix[r][c] <= 0) {
	string e(
		 "Slope map contains negative or zero values inside the valid domain...\n");
	throw e;
      }

      if (_porosity0->matrix[r][c] == _porosity0->nodata) {
	string e(
		 "Base porosity map contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_porosity0->matrix[r][c] <= 0) {
	string e(
		 "Base porosity map contains negative or zero values inside the valid domain...\n");
	throw e;
      }
      if(ctrl.sw_expPoros){
	if (_kporos->matrix[r][c] == _kporos->nodata) {
	  string e(
		   "Porosity profile map contains no data values inside the valid domain...\n");
	  throw e;
	}
	if (_kporos->matrix[r][c] <= 0) {
	  string e(
		   "Porosity profile map contains negative or zero values inside the valid domain...\n");
	  throw e;
	}
      }

      if (_psi_ae->matrix[r][c] == _psi_ae->nodata) {
	string e(
		 "Air entry pressure map contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_psi_ae->matrix[r][c] <= 0) {
	string e(
		 "Air entry pressure map contains negative or zero values inside the valid domain...\n");
	throw e;
      }

      if (_BClambda->matrix[r][c] == _BClambda->nodata) {
	string e(
		 "Brooks and Corey lambda map contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_BClambda->matrix[r][c] <= 1) {
	string e(
		 "WARNING: Brooks and Corey lambda map is too small, switching to minimum value of 2...\n");
	cout << e;
	_BClambda->matrix[r][c] = 1;
	//throw e;
      }

      // Initial check of theta_r using L1 (it is the same for all 3 layers,
      // at least before comparing to _porosityL*) 
      if (_theta_rL1->matrix[r][c] == _theta_rL1->nodata) {
	string e(
		 "residual moisture map contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_theta_rL1->matrix[r][c] <= 0) {
	string e(
		 "residual moisture map contains negative or zero values inside the valid domain...\n");
	throw e;
      }
      if (_theta_rL1->matrix[r][c] > _porosity0->matrix[r][c]) {
	string e(
		 "Input residual soil moisture map is larger than top-of-profile porosity inside the valid domain...\n");
	throw e;
      }
      
      if (_soildepth->matrix[r][c] == _soildepth->nodata) {
	string e(
		 "soil depth map contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_soildepth->matrix[r][c] < 0) {
	string e(
		 "soil depth map contains negative values inside the valid domain...\n");
	throw e;
      }
      /*
	if (_Kroot->matrix[r][c] == _Kroot->nodata) {
	string e(
	"root profile map contains no data inside the valid domain...\n");
	throw e;
	}
      */
      if (_paramWc->matrix[r][c] == _paramWc->nodata) {
	string e(
		 "Soil parameter Wc map contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_paramWc->matrix[r][c] <= 0) {
	string e(
		 "Soil parameter Wc map contains negative or zero values inside the valid domain...\n");
	throw e;
      }

      if (_paramWp->matrix[r][c] == _paramWp->nodata) {
	string e(
		 "Soil parameter Wp map contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_paramWp->matrix[r][c] <= 0) {
	string e(
		 "Soil parameter Wp map contains negative or zero values inside the valid domain...\n");
	throw e;
      }

      if (_meltCoeff->matrix[r][c] == _meltCoeff->nodata) {
	string e(
		 "Snowmelt coefficient map contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_meltCoeff->matrix[r][c] <= 0) {
	string e(
		 "Snowmelt coefficient map contains negative or zero values inside the valid domain...\n");
	throw e;
      }

      if (_snow->matrix[r][c] == _snow->nodata) {
	string e(
		 "Initial SWE map contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_snow->matrix[r][c] < 0) {
	string e(
		 "Initial SWE map contains negative values inside the valid domain...\n");
	throw e;
      }

      if (_albedo->matrix[r][c] == _albedo->nodata) {
	string e(
		 "Albedo map contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_albedo->matrix[r][c] <= 0) {
	string e(
		 "Albedo map contains negative or zero values inside the valid domain...\n");
	throw e;
      }

      if (_emiss_surf->matrix[r][c] == _emiss_surf->nodata) {
	string e(
		 "Surface emissivity map contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_emiss_surf->matrix[r][c] <= 0) {
	string e(
		 "Surface emissivity map contains negative or zero values inside the valid domain...\n");
	throw e;
      }

      if (_soil_dry_heatcap->matrix[r][c] == _soil_dry_heatcap->nodata) {
	string e(
		 "Soil dry heat capacity map contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_soil_dry_heatcap->matrix[r][c] <= 0) {
	string e(
		 "Soil dry heat capacity map contains negative or zero values inside the valid domain...\n");
	throw e;
      }

      if (_soil_dry_thermcond->matrix[r][c]
	  == _soil_dry_thermcond->nodata) {
	string e(
		 "Dry soil thermal conductivity map contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_soil_dry_thermcond->matrix[r][c] <= 0) {
	string e(
		 "Dry soil thermal conductivity map contains negative or zero values inside the valid domain...\n");
	throw e;
      }

      if (_dampdepth->matrix[r][c] == _dampdepth->nodata) {
	string e(
		 "Thermal soil damping depth map contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_dampdepth->matrix[r][c] <= 0) {
	string e(
		 "Thermal soil damping depth map contains negative or zero values inside the valid domain...\n");
	throw e;
      }

      if (_Temp_d->matrix[r][c] == _Temp_d->nodata) {
	string e(
		 "Soil temp at damp depth map contains no data values inside the valid domain...\n");
	throw e;
      }

      // Compare theta_rL* and porosityL* (to avoid inconsistencies when exp profiel is activated)
      if (_theta_rL1->matrix[r][c] > 0.25 * _porosityL1->matrix[r][c]) {
	string e("WARNING: Topsoil residual soil moisture is > 0.25 * topsoil porosity, let's tone it down...\n");
	cout << e;
	_theta_rL1->matrix[r][c] = 0.25 * _porosityL1->matrix[r][c];
	// throw e;
      }
      if (_theta_rL2->matrix[r][c] > 0.25 * _porosityL2->matrix[r][c]) {
	string e("WARNING: Residual soil moisture in layer 2 is > 0.25*porosity, let's tone it down...\n");
	cout << e;
	_theta_rL2->matrix[r][c] = 0.25 * _porosityL2->matrix[r][c];
	// throw e;
      }
      if (_theta_rL3->matrix[r][c] > 0.25 * _porosityL3->matrix[r][c]) {
	string e("WARNING: Residual soil moisture in layer 3 is > 0.25 * porosity, let's tone it down...\n");
	cout << e;
	_theta_rL3->matrix[r][c] = 0.25 * _porosityL3->matrix[r][c];
	// throw e;
      }

      // Check soil moisture itself
      if (_soilmoist1->matrix[r][c] == _soilmoist1->nodata) {
	string e(
		 "Initial soil moisture map contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_soilmoist2->matrix[r][c] == _soilmoist2->nodata) {
	string e(
		 "Initial soil moisture map in layer 2 contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_soilmoist3->matrix[r][c] == _soilmoist3->nodata) {
	string e(
		 "Initial soil moisture map in layer 3 contains no data values inside the valid domain...\n");
	throw e;
      }
      if (_soilmoist1->matrix[r][c] < _theta_rL1->matrix[r][c]) {
	string e("WARNING: Topsoil Initial soil moisture map is lower than residual moisture, let's jump it up...\n");
	cout << e;
	_soilmoist1->matrix[r][c] = (_theta_rL1->matrix[r][c] + 3*_porosityL1->matrix[r][c]) / 4;
	// throw e;
      }
      if (_soilmoist2->matrix[r][c] < _theta_rL2->matrix[r][c]) {
	string e("WARNING: Initial soil moisture map in layer 2 is lower than residual moisture, let's jump it up...\n");
	cout << e;
	_soilmoist2->matrix[r][c] = (_theta_rL2->matrix[r][c] + 3*_porosityL2->matrix[r][c]) / 4;
	// throw e ;
      }
      if (_soilmoist3->matrix[r][c] < _theta_rL3->matrix[r][c]) {
	string e("WARNING: Initial soil moisture map in layer 3 is lower than residual moisture, let's jump it up...\n");
	cout << e;
	_soilmoist3->matrix[r][c] = (_theta_rL3->matrix[r][c] + 3*_porosityL3->matrix[r][c]) / 4;
	//throw e;
      }
      if (_soilmoist1->matrix[r][c] > _porosityL1->matrix[r][c]) {
	string e("WARNING: Topsoil Initial soil moisture map is larger than porosity, let's bring it down...\n");
	cout << e;
	_soilmoist1->matrix[r][c] = (_theta_rL1->matrix[r][c] + 3*_porosityL1->matrix[r][c]) / 4;
	//throw e;
      }
      if (_soilmoist2->matrix[r][c] > _porosityL2->matrix[r][c]) {
	string e("WARNING: Initial soil moisture in layer 2 is larger than porosity, let's bring it down...\n");
	cout << e;
	_soilmoist2->matrix[r][c] = (_theta_rL2->matrix[r][c] + 3*_porosityL2->matrix[r][c]) / 4;
	//throw e;
      }
      if (_soilmoist3->matrix[r][c] > _porosityL3->matrix[r][c]) {
	string e("WARNING: Initial soil moisture in layer 3 is larger than porosity, let's bring it down...\n");
	cout << e;
	_soilmoist3->matrix[r][c] = (_theta_rL3->matrix[r][c] + 3*_porosityL3->matrix[r][c]) / 4;
	//throw e;
      }
      if (_channelwidth->matrix[r][c] < 0) {
	string e("The channel width map contains negative values\n");
	throw e;
      }
      if (_channellength->matrix[r][c] < 0) {
	string e("The channel length map contains negative values\n");
	throw e;
      }      
      if(ctrl.sw_deepGW){
	if (_Hydrofrac_DeepGW->matrix[r][c] > 0.5* _DeepGW->matrix[r][c]) {
	  string e(
		   "The hydrologically active fraction of deep groundwater is larger than 0.5...\n");
	  throw e;
	}
      }
    } catch (string &e) {
      cout << e;
      cout << "In row " << r << " col " << c << endl;
    }

  }

  if (excep_thrown)
    throw;
}
