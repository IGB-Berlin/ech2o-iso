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
 * GrowGrassLAI.cpp
 *
 *  Created on: Nov 11, 2012
 *      Author: marco
 */

#include "Forest.h"

int Forest::GrowGrassLAI(UINT4 spec, UINT4 row, UINT4 col, REAL8 dt){

	REAL8 SLA_g = _species[spec].SLA;

	REAL8 DryLeafTurnoverRate = _species[spec].Fprn;
	REAL8 DryLeafTurnoverRateAdjustmentParam = _species[spec].Fpra;
	REAL8 DryDecayAdjustment = min<REAL8>(max<REAL8>(_species[spec]._Temp_c->matrix[row][col] / DryLeafTurnoverRateAdjustmentParam,0), 1);
	REAL8 GreenLaiDecay;
	REAL8 DryLaiDecay;

	GreenLaiDecay = _species[spec].LeafTurnover * dt * _species[spec]._grassLAI_g->matrix[row][col];

	DryLaiDecay = DryLeafTurnoverRate * DryDecayAdjustment * dt * _species[spec]._grassLAI_d->matrix[row][col];

	_species[spec]._grassLAI_g->matrix[row][col] += max<REAL8>(-0.95*_species[spec]._grassLAI_g->matrix[row][col], 
							_species[spec]._Del_FoliageMass->matrix[row][col] * SLA_g - GreenLaiDecay);

	_species[spec]._grassLAI_d->matrix[row][col] += GreenLaiDecay - DryLaiDecay;

	_species[spec]._LAI->matrix[row][col] = _species[spec]._grassLAI_g->matrix[row][col] +
  						_species[spec]._grassLAI_d->matrix[row][col];

	return EXIT_SUCCESS;
}


