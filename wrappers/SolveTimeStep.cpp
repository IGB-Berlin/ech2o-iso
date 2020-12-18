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
 * SolveTimeStep.cpp
 *
 *  Created on: Aug 2, 2010
 *      Author: Marco.Maneta
 */
#include <iostream>
#include "Sativa.h"


int SolveTimeStep(){

  oBasin->SolveCanopyFluxes(*oAtmosphere, *oControl, *oTracking);
  oBasin->SolveSurfaceFluxes(*oAtmosphere, *oControl, *oTracking);
  oBasin->CalculateGrowForest(*oAtmosphere, *oControl);
  oBasin->DailyGWRouting(*oAtmosphere, *oControl, *oTracking);
  oBasin->CalculateSatArea(*oAtmosphere, *oControl);

  // If tracking...	  
  if(oControl->sw_trck){
    // If Two-pore...
    if(oControl->sw_TPD){
      // Calculate the soil-averaged tracer values
      oTracking->CalcTPDtoLayers(*oBasin, *oControl);
      // Calculate the relative fraction of tightly-bound domain
      oBasin->CalcFracMobileWater();
    }

    if(oControl->sw_2H){ 
      if(oControl->Rep_d2HsoilUp || oControl->RepTs_d2HsoilUp)
	oTracking->Calcd2Hsoil_12(*oBasin);
      if(oControl->Rep_d2HsoilAv || oControl->RepTs_d2HsoilAv)
	oTracking->Calcd2Hsoil_Av(*oBasin);
    }
    
    if(oControl->sw_18O){
      if(oControl->Rep_d18OsoilUp || oControl->RepTs_d18OsoilUp)
	oTracking->Calcd18Osoil_12(*oBasin);
      if(oControl->Rep_d18OsoilAv || oControl->RepTs_d18OsoilAv)
	oTracking->Calcd18Osoil_Av(*oBasin);
    }
	  
    if(oControl->sw_Age){
      // Increment age by one time step duration
      oTracking->IncrementAge(*oBasin, *oControl);
      // Reported quantities
      if(oControl->Rep_AgesoilUp || oControl->RepTs_AgesoilUp)
	oTracking->CalcAgesoil_12(*oBasin);
      if(oControl->Rep_AgesoilAv || oControl->RepTs_AgesoilAv)
	oTracking->CalcAgesoil_Av(*oBasin);

    }
  }

  return EXIT_SUCCESS;
}
