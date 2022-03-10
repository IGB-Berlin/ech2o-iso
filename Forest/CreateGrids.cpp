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
 * CreateGroves.cpp
 *
 *  Created on: Jun 18, 2010
 *      Author: marco
 */

#include "Grove.h"

int Grove::CreateGrids(grid *base){

	try{
		_fraction = new grid (*base);
		_StemDensity = new grid (*base);
		_LAI = new grid (*base);
		_grassLAI_g = new grid (*base);
		_grassLAI_d = new grid (*base);
		_AGE = new grid (*base);
		_CanopyConductance = new grid (*base);
		_GPP = new grid (*base);
		_NPP = new grid (*base);
		_BasalArea = new grid (*base);
		_Height = new grid (*base);
		_RootMass = new grid (*base);
		_Del_FoliageMass = new grid (*base);
		_Del_StemMass = new grid (*base);
		_Del_RootMass = new grid (*base);
		_Temp_c = new grid(*base);
		_RUptakeL1 = new grid (*base);
		_RUptakeL2 = new grid (*base);
		_RUptakeL3 = new grid (*base);		
		_NetR_Can = new grid (*base);
		_LatHeat_CanE = new grid (*base);
		_LatHeat_CanT = new grid (*base);
		_SensHeat_Can = new grid (*base);
		_WaterStorage = new grid (*base);
		_ET = new grid (*base);
		_Einterception = new grid (*base);
		_Transpiration = new grid (*base);
		_TranspirationFlux = new grid (*base);
		_Esoil = new grid (*base);
		_SoilWatPot = new grid (*base);		
		_LeafWatPot = new grid (*base);
		_SapVelocity = new grid (*base);
		_rootfrac1 = new grid (*base);
		_rootfrac2 = new grid (*base);

	}catch(const exception& e){

		cerr << "Failed allocate memory for Grove grid object \n" << e.what() << endl;
		exit(EXIT_FAILURE);
	}

	return EXIT_SUCCESS;
}
