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
 * GroveDestruct.cpp
 *
 *  Created on: Jun 17, 2010
 *      Author: marco
 */

#include "Grove.h"

Grove::~Grove(){

	if(_fraction)
		delete _fraction;
	if(_StemDensity)
		delete _StemDensity;
	if(_LAI)
		delete _LAI;
	if(_grassLAI_g)
		delete _grassLAI_g;
	if(_grassLAI_d)
		delete _grassLAI_d;
	if(_AGE)
		delete _AGE;
	if(_CanopyConductance)
		delete _CanopyConductance;
	if(_GPP)
		delete _GPP;
	if(_NPP)
		delete _NPP;
	if(_BasalArea)
		delete _BasalArea;
	if(_Height)
		delete _Height;
	if(_RootMass)
		delete _RootMass;
	if(_Del_FoliageMass)
		delete _Del_FoliageMass;
	if(_Del_StemMass)
		delete _Del_StemMass;
	if(_Del_RootMass)
		delete _Del_RootMass;
	if(_Temp_c)
		delete _Temp_c;
        if(_RUptakeL1)
		delete _RUptakeL1;
	if(_RUptakeL2)
		delete _RUptakeL2;
	if(_RUptakeL3)
		delete _RUptakeL3;	
	if(_NetR_Can)
		delete _NetR_Can;
	if(_LatHeat_CanE)
		delete _LatHeat_CanE;
	if(_LatHeat_CanT)
		delete _LatHeat_CanT;
	if(_SensHeat_Can)
		delete _SensHeat_Can;
	if(_WaterStorage)
		delete _WaterStorage;
	if(_ET)
		delete _ET;
	if(_Transpiration)
		delete _Transpiration;
	if(_Einterception)
		delete _Einterception;
	if(_Esoil)
			delete _Esoil;
	if(_SoilWatPot)
		delete _SoilWatPot;	
	if(_LeafWatPot)
		delete _LeafWatPot;
	if(_SapVelocity)
	  delete _SapVelocity;
	if(_rootfrac1)
		delete _rootfrac1;
	if(_rootfrac2)
		delete _rootfrac2;

	// Tracking
	if(_d2Hcanopy)
		delete _d2Hcanopy;
	if(_d18Ocanopy)
		delete _d18Ocanopy;
	if(_Agecanopy)
		delete _Agecanopy;

	if(_d2Hthroughfall)
		delete _d2Hthroughfall;
	if(_d18Othroughfall)
		delete _d18Othroughfall;
	if(_Agethroughfall)
		delete _Agethroughfall;

	if(_d2HevapT)
		delete _d2HevapT;
	if(_d18OevapT)
		delete _d18OevapT;
	if(_AgeevapT)
		delete _AgeevapT;

	if(_d2HevapI)
		delete _d2HevapI;
	if(_d18OevapI)
		delete _d18OevapI;
	if(_AgeevapI)
		delete _AgeevapI;

	if(_d2HevapS)
		delete _d2HevapS;
	if(_d18OevapS)
		delete _d18OevapS;
	if(_AgeevapS)
		delete _AgeevapS;

	if(ifLAI.is_open())
	  ifLAI.close();
}
