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
 * CalculateCanopyConduct.cpp
 *
 *  Created on: Jun 23, 2010
 *      Author: Marco.Maneta
 */

#include "Forest.h"

int Forest::CalculateCanopyConduct(const Basin &bas, const Atmosphere &atm,
		const Control &ctrl, const double &lwp, const double &lwp_max, const double &lwp_min,
		double &dgsdlwp, UINT4 j, UINT4 r, UINT4 c) {

	REAL8 sw_rad;
	REAL8 gsmax;
	REAL8 f_light_coeff, f_vpd_coeff;
	REAL8 lai;
	REAL8 airTemp, optTemp, maxTemp, minTemp;
	REAL8 es, ea, vpd;
	REAL8 f_light, f_temp, f_vpd, f_psi;
	REAL8 shelter_factor;
	REAL8 gs;
	REAL8 gsmin = 0.00000000000001; //arbitrary number to set as minimum canopy conductance so its inverse (resistance) doesn't explode when divided by zero conductance

	shelter_factor = 1; //account to the shading effect of leaves

	sw_rad = atm.getIncomingShortWave()->matrix[r][c];
	gsmax = _species[j].gsmax;
	lai = _species[j]._LAI->matrix[r][c];
	airTemp = atm.getTemperature()->matrix[r][c];
	optTemp = _species[j].TempOpt;
	maxTemp = _species[j].TempMax;
	minTemp = _species[j].TempMin;

	f_light_coeff = _species[j].gs_light_coeff;
	f_vpd_coeff = _species[j].gs_vpd_coeff;
	es = SatVaporPressure(airTemp); //TODO urgent: calculate vapor pressure deficit with respect to soil
	ea = es * atm.getRelativeHumidty()->matrix[r][c];
	vpd = es - ea;

	f_light = Calculate_gs_light(sw_rad, f_light_coeff); //TODO: implement light extinction factor when other species in the same cell is taller and has large LAI
	f_temp = Calculate_ft(airTemp, maxTemp, minTemp, optTemp);
	// Commented out because it is redundant: vpd limitation is already taken into account in LET calculation
	// EDIT: it is left here but the parameter value can be set as very small to effectively have no dependence
	f_vpd = Calculate_gs_vpd(vpd, f_vpd_coeff);

	// Nonlinear or linear stomatal model
	if(ctrl.toggle_sm == 0){
	  f_psi = Calculate_gs_lwp_nonlinear(lwp, lwp_max, lwp_min);
	} else {
	  f_psi = Calculate_gs_lwp_linear(lwp, lwp_max, lwp_min);	
	}

	gs = gsmax * lai * shelter_factor * f_light * f_temp * f_vpd * f_psi;

	_species[j]._CanopyConductance->matrix[r][c] = gs < gsmin ? gsmin : gs;

	dgsdlwp =  gsmax * lai* shelter_factor * f_light * f_temp * f_vpd;

	return EXIT_SUCCESS;
}
