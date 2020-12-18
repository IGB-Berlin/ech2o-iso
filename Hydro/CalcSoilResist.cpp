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
 * CalcSoilResist.cpp
 *
 *  Created on: Mar 5, 2010
 *      Author: Marco Maneta
 */

#include"Basin.h"

REAL8 Basin::CalcSoilResist(double &theta, int row, int col, UINT4 option){


  // Passerat de Silans et al. (1989)
	const double rmin = 3.8113e4; //minimum soil resistance in s.m-1
	const double Cs = 13.515; //empirical parameter

	// Values for CLM model formulation as per Sakaguchi and Zeng (2009):
	//"Effects of soil wetness, plant litter, and under-canopy atmospheric
	// stability on ground evaporation in the Community Land Model", JGR
	const double Do = 2.2e-5; //Molecular Diffusion coefficient of water vapor in the air (m2s-1)
	const double w = 5; //empirical shape parameter
	const double e = 2.718; //empirical shape parameter
	double D = 0;
	double L = 0;

	double d1 = 0;

	double thetas = _porosityL1->matrix[row][col];
	double thetafc = _fieldcapL1->matrix[row][col];
	double thetar = _theta_rL1->matrix[row][col];


	double S = (theta - thetar) / (thetafc - thetar);

	if (theta < RNDOFFERR || thetafc < RNDOFFERR)
		return 0;

	if (option == 0)
		return 0;
	// Passerat de Silans et al. (1989)
	if (option == 1)
		return rmin * expl(-Cs*(S));
	// Sellers et al. (1992)
	else if (option == 2)
		return ( expl(8.206-4.255*S) );
	// Sakaguchi and Zeng (2009)
	else if (option == 3){
		 d1 = _depth_layer1->matrix[row][col];
		 L = d1 * ( ( expl(powl(1 - (theta/thetas),w) ) - 1 ) / (e - 1) );
	     D = Do*thetas*thetas* powl(1 - (thetar/thetas), 2+3 * _BClambdaL1->matrix[row][col]);
	     return L/D;
	}
	else{
		cout << "Wrong option in the Soil_resistance_opt toggle switch. Please verify the configuration file" << endl;
		exit(EXIT_FAILURE);
	}



}
