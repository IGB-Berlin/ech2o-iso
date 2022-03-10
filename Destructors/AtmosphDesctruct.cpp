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
 * AtmosphDesctruct.cpp
 *
 *  Created on: Oct 15, 2009
 *      Author: Marco Maneta
 */

#include "Atmosphere.h"

Atmosphere::~Atmosphere(){

	if(_zones)
		delete _zones;
	if(_Ldown)
		delete _Ldown;
	if(_Sdown)
		delete _Sdown;
	if(_Tp)
		delete _Tp;
	if(_MaxTp)
		delete _MaxTp;
	if(_MinTp)
		delete _MinTp;
	if(_Precip)
		delete _Precip;
	if(_Rel_humid)
		delete _Rel_humid;
	if(_Wind_speed)
		delete _Wind_speed;
        if(_Pa)
                delete _Pa;
	if(_Anthr_Heat)
		delete _Anthr_Heat;	
	if(_zoneId)
		delete[] _zoneId;
	if(_isohyet)
		delete _isohyet;
	if(_d2Hprecip)
		delete _d2Hprecip;
	if(_d18Oprecip)
		delete _d18Oprecip;


	if(ifLdown.is_open())
		ifLdown.close();
	if(ifSdown.is_open())
		ifSdown.close();
	if(ifTp.is_open())
		ifTp.close();
	if(ifMaxTp.is_open())
		ifMaxTp.close();
	if(ifMinTp.is_open())
		ifMinTp.close();
	if(ifPrecip.is_open())
		ifPrecip.close();
	if(ifRelHumid.is_open())
		ifRelHumid.close();
	if(ifWindSpeed.is_open())
		ifWindSpeed.close();
	if(ifAnthrHeat.is_open())
		ifAnthrHeat.close();	
	if(ifd2Hprecip.is_open())
		ifd2Hprecip.close();
	if(ifd18Oprecip.is_open())
		ifd18Oprecip.close();
}
