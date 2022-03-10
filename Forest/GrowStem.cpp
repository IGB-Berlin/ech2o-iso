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
 * GrowStem.cpp
 *
 *  Created on: Jun 30, 2010
 *      Author: Marco.Maneta
 */

#include "Forest.h"

int Forest::GrowStem(UINT4 spec, UINT4 row, UINT4 col){

	REAL8 delD;

	REAL8 numTrees = _dx * _dx * _species[spec]._fraction->matrix[row][col] *
	  			_species[spec]._StemDensity->matrix[row][col];
	REAL8 H = _species[spec]._Height->matrix[row][col];
	REAL8 D = 2 * powl(_species[spec]._BasalArea->matrix[row][col] / PI, 0.5);
	REAL8 G = _species[spec]._Del_StemMass->matrix[row][col] / _species[spec]._StemDensity->matrix[row][col]; //average mass increment per individual tree (grams per tree)
	REAL8 gamma = _species[spec].WoodDensity;
	REAL8 phi = _species[spec].TreeShapePar;
	REAL8 lambda = _species[spec].Crown2Stem;
	REAL8 maxAge = _species[spec].MaxAge;
	REAL8 Age = _species[spec]._AGE->matrix[row][col];
	REAL8 Fhdmin = _species[spec].Fhdmin;
	REAL8 Fhdmax = _species[spec].Fhdmax;
	REAL8 Stc; //relative coverage of the stand (m2 m-2) to calculate competition
	REAL8 Fhd = 0;

	Stc = (PI*0.25) * (lambda * D) * (lambda * D) * numTrees / (_dx*_dx);

	if (Stc < 1 && H/D >= Fhdmin)
		Fhd = Fhdmin;
	else if (Stc >= 1 && H/D <= Fhdmax)
		Fhd = Age > 0.5*maxAge ? Fhdmax : Fhdmin;
	else if( (H/D) < Fhdmin )
		Fhd = Fhdmax;
	else if( (H/D) > Fhdmax )
		Fhd = 0.5 * Fhdmin;

	if (Age > 0.7*maxAge)
		Fhd = 0;

	delD = 4*G /( gamma*PI*phi*D*D*(2*(H/D)+Fhd) );

	D += delD;

	_species[spec]._BasalArea->matrix[row][col] = D * D * 0.25 * 3.14159;

	_species[spec]._Height->matrix[row][col] += Fhd * delD;

	return EXIT_SUCCESS;
}
