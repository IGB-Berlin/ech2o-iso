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
 * Report2screen.cpp
 *
 *  Created on: Aug 2, 2010
 *      Author: Marco.Maneta
 */

#include "Sativa.h"

int Report2Screen(){

  UINT4 ncell = oBasin->getSortedGrid().cells.size();
  REAL8 dx = oBasin->getCellSize();
  REAL8 area = ncell * dx * dx;

  /*
  // A few tracking reports here, so that 
  // the age mass balance check uses beg-of-time-step values, more simple!
  if(oControl->sw_trck){
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
    }*/

  // ==== BasinSummary.txt --------------------------------------------------------
  // -----------------------------------------------------

  // For the splash screen, using mm for flux and stores
  // makes things more obvious than m3 (to check quickly)!
  // m3 is kept for BasinSummary.txt
  
  //printf("\nTotal Precipitation (m3):  %.2f \t", oBudget->precipitation);
  printf("\nTotal Precipitation (mm):  %.2f \t", 1000 * oBudget->precipitation / area);
  ofSummary << oBudget->precipitation << "\t";

  //printf("SWE (m3): %.2f \n", oBudget->snowpack);
  printf("SWE (mm): %.2f \n", 1000 * oBudget->snowpack / area);
  ofSummary << oBudget->snowpack << "\t";

  //printf("Canopy Storage (m3): %.2f \t", oBudget->canopy);
  printf("Canopy Storage (mm): %.2f \t", 1000 * oBudget->canopy /area);
  ofSummary << oBudget->canopy << "\t";

  //printf("Ponding (m3): %.2f \n", oBudget->ponding);
  printf("Ponding (mm): %.2f \n", 1000 * oBudget->ponding / area);
  ofSummary << oBudget->ponding << "\t";

  //printf("Channel (m3): %.2f \n", oBudget->chan_store);
  printf("Channel (mm): %.2f \n", 1000 * oBudget->chan_store / area);
  ofSummary << oBudget->chan_store << "\t";

  //printf("Soil water (m3): %.2f \t", oBudget->vadose);
  printf("Soil water (mm): %.2f \t", 1000 * oBudget->vadose / area);
  ofSummary << oBudget->vadose << "\t";

  // Only for budget Summary: each layer
  ofSummary << oBudget->soilL1 << "\t";
  ofSummary << oBudget->soilL2 << "\t";
  ofSummary << oBudget->soilL3 << "\t";
  
  //printf("of which Groundwater (m3): %.2f \n", oBudget->grndwater);
  printf("of which RootZone water (mm): %.2f \t", 1000 * oBudget->rootzone / area);
  ofSummary << oBudget->rootzone << "\t";
  
  //printf("of which Groundwater (m3): %.2f \n", oBudget->grndwater);
  printf(", and Groundwater (mm): %.2f \n", 1000 * oBudget->grndwater / area);
  ofSummary << oBudget->grndwater << "\t";

  //printf("Total Evapotranspiration (m3): %.2f \t", oBudget->evaporation);
  printf("Total Evapotranspiration (mm): %.2f \t", 1000 * oBudget->evaporation / area);
  ofSummary << oBudget->evaporation << "\t";

  //printf("Total Soil Evaporation (m3): %.2f \n", oBudget->evaporationS);
  printf("Total Soil Evaporation (mm): %.2f \n", 1000 * oBudget->evaporationS /area);
  ofSummary << oBudget->evaporationS << "\t";

  printf("Total Channel Evaporation (mm): %.2f \n", 1000 * oBudget->evaporationC /area);
  ofSummary << oBudget->evaporationC << "\t";

  //printf("Total Canopy Evaporation (m3): %.2f \t", oBudget->evaporationI);
  printf("Total Canopy Evaporation (mm): %.2f \t", 1000 * oBudget->evaporationI / area);
  ofSummary << oBudget->evaporationI << "\t";

  //printf("Total Transpiration (m3): %.2f \n", oBudget->transpiration);
  printf("Total Transpiration (mm): %.2f \n", 1000 * oBudget->transpiration / area);
  ofSummary << oBudget->transpiration << "\t";

  //printf("Bedrock Leak (m3): %.2f \n", oBudget->leakage);
  printf("Bedrock Leak (mm): %.2f \n", 1000 * oBudget->leakage /area);
  ofSummary << oBudget->leakage << "\t";

  //printf("Total OvlndFlow output (m3): %.2f \t", oBudget->ovlndflow);
  printf("Total OvlndFlow output (mm): %.2f \t", 1000 * oBudget->ovlndflow / area);
  ofSummary << oBudget->ovlndflow << "\t";

  //printf("Total GWFlow output (m3): %.2f \n", oBudget->gwtrflow);
  printf("Total GWFlow output (mm): %.2f \n", 1000 * oBudget->gwtrflow / area);
  ofSummary << oBudget->gwtrflow << "\t";

  //printf("Run-off to channel (m3): %.2f \t", oBudget->srftochn);
  printf("Run-off to channel (m3): %.2f \t", 1000 * oBudget->srftochn / area);
  ofSummary << oBudget->srftochn << "\t";

  //printf("GW to channel (m3): %.2f \n", oBudget->gwtochn);
  printf("GW to channel (m3): %.2f \n", 1000 * oBudget->gwtochn / area);
  ofSummary << oBudget->gwtochn << "\t";

  //printf("GW recharge (m3): %.2f \t", oBudget->recharge);
  printf("GW recharge (mm): %.2f \t", 1000 * oBudget->recharge / area);
  ofSummary << oBudget->recharge << "\t";

  // Saturated area (% of the catchment)
  printf("Saturated area fraction: %.2f \n", oBudget->satarea);
  ofSummary << oBudget->satarea << "\t";

  printf("Mass Balance Error (%): %e \n", oBudget->MBErr);
  ofSummary << oBudget->MBErr ;//<< "\t";

  if(oControl->sw_trck and oControl->sw_2H){
    printf("Deuterium Mass Balance Error (%): %e \n", oBudget->MBErr_d2H);
    ofSummary << "\t" << oBudget->MBErr_d2H ;
  }

  if(oControl->sw_trck and oControl->sw_18O){
    printf("Oxygen 18 Mass Balance Error (%): %e \n", oBudget->MBErr_d18O);
    ofSummary << "\t" << oBudget->MBErr_d18O ;
  }

  if(oControl->sw_trck and oControl->sw_Age){
    printf("Age Mass Balance Error (%): %e \n", oBudget->MBErr_Age);
    ofSummary << "\t" << oBudget->MBErr_Age ;
  }

  ofSummary << "\n";

  // ==== Basind2HSummary.txt --------------------------------------------------------
  // -----------------------------------------------------
  if(oControl->sw_trck and oControl->sw_2H){
    ofd2HSummary << oBudget->d2HTot << "\t";
    ofd2HSummary << oBudget->d2Hsnowpack << "\t";
    ofd2HSummary << oBudget->d2Hcanopy << "\t";
    ofd2HSummary << oBudget->d2Hponding << "\t";
    ofd2HSummary << oBudget->d2Hchannel << "\t";
    ofd2HSummary << oBudget->d2Hvadose << "\t";
    ofd2HSummary << oBudget->d2HsoilL1 << "\t";
    ofd2HSummary << oBudget->d2HsoilL2 << "\t";
    ofd2HSummary << oBudget->d2HsoilL3 << "\t";
    ofd2HSummary << oBudget->d2Hrootzone << "\t";
    ofd2HSummary << oBudget->d2Hgrndwater << "\t";
    ofd2HSummary << oBudget->d2HET << "\t";
    ofd2HSummary << oBudget->d2HevapS << "\t";
    ofd2HSummary << oBudget->d2HevapC << "\t";
    ofd2HSummary << oBudget->d2HevapI << "\t";
    ofd2HSummary << oBudget->d2HevapT << "\t";
    ofd2HSummary << oBudget->d2Hleakage << "\t";
    ofd2HSummary << oBudget->d2HOvlndOut << "\t";
    ofd2HSummary << oBudget->d2HGWOut << "\t";
    ofd2HSummary << oBudget->d2HOut << "\t";
    ofd2HSummary << oBudget->d2Hsrftochn << "\t";
    ofd2HSummary << oBudget->d2Hgwtochn << "\t";   
    ofd2HSummary << oBudget->d2Hrecharge << "\t";
    ofd2HSummary << oBudget->d2Hprecip << "\n";

  }

  // ==== Basind18OSummary.txt --------------------------------------------------------
  // -----------------------------------------------------
  if(oControl->sw_trck and oControl->sw_18O){
    ofd18OSummary << oBudget->d18OTot << "\t";
    ofd18OSummary << oBudget->d18Osnowpack << "\t";
    ofd18OSummary << oBudget->d18Ocanopy << "\t";
    ofd18OSummary << oBudget->d18Oponding << "\t";
    ofd18OSummary << oBudget->d18Ochannel << "\t";
    ofd18OSummary << oBudget->d18Ovadose << "\t";
    ofd18OSummary << oBudget->d18OsoilL1 << "\t";
    ofd18OSummary << oBudget->d18OsoilL2 << "\t";
    ofd18OSummary << oBudget->d18OsoilL3 << "\t";
    ofd18OSummary << oBudget->d18Orootzone << "\t";
    ofd18OSummary << oBudget->d18Ogrndwater << "\t";
    ofd18OSummary << oBudget->d18OET << "\t";
    ofd18OSummary << oBudget->d18OevapS << "\t";
    ofd18OSummary << oBudget->d18OevapC << "\t";
    ofd18OSummary << oBudget->d18OevapI << "\t";
    ofd18OSummary << oBudget->d18OevapT << "\t";
    ofd18OSummary << oBudget->d18Oleakage << "\t";
    ofd18OSummary << oBudget->d18OOvlndOut << "\t";
    ofd18OSummary << oBudget->d18OGWOut << "\t";
    ofd18OSummary << oBudget->d18OOut << "\t";
    ofd18OSummary << oBudget->d18Osrftochn << "\t";
    ofd18OSummary << oBudget->d18Ogwtochn << "\t";   
    ofd18OSummary << oBudget->d18Orecharge << "\t";
    ofd18OSummary << oBudget->d18Oprecip << "\n";

  }

  // ==== BasinAgeSummary.txt --------------------------------------------------------
  // -----------------------------------------------------
  if(oControl->sw_trck and oControl->sw_Age){
    ofAgeSummary << oBudget->AgeTot << "\t";
    ofAgeSummary << oBudget->Agesnowpack << "\t";
    ofAgeSummary << oBudget->Agecanopy << "\t";
    ofAgeSummary << oBudget->Ageponding << "\t";
    ofAgeSummary << oBudget->Agechannel << "\t";
    ofAgeSummary << oBudget->Agevadose << "\t";
    ofAgeSummary << oBudget->AgesoilL1 << "\t";
    ofAgeSummary << oBudget->AgesoilL2 << "\t";
    ofAgeSummary << oBudget->AgesoilL3 << "\t";
    ofAgeSummary << oBudget->Agerootzone << "\t";
    ofAgeSummary << oBudget->Agegrndwater << "\t";
    ofAgeSummary << oBudget->AgeET << "\t";
    ofAgeSummary << oBudget->AgeevapS << "\t";
    ofAgeSummary << oBudget->AgeevapI << "\t";
    ofAgeSummary << oBudget->AgeevapT << "\t";
    ofAgeSummary << oBudget->Ageleakage << "\t";
    ofAgeSummary << oBudget->AgeOvlndOut << "\t";
    ofAgeSummary << oBudget->AgeGWOut << "\t";
    ofAgeSummary << oBudget->AgeOut << "\t";
    ofAgeSummary << oBudget->Agesrftochn << "\t";
    ofAgeSummary << oBudget->Agegwtochn << "\t";   
    ofAgeSummary << oBudget->Agerecharge << "\n";
  }

  return EXIT_SUCCESS;
}
