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

  // ==== SpashScreen -------------------------------------------------------------
  // For the splash screen, using mm for flux and stores
  // makes things more obvious than m3 (to check quickly)!
  printf("\n%-31s%-15.2f","Total Precipitation (mm): ", 1000 * oBudget->precipitation / area);
  printf("%-31s%-15.2f","SWE (mm): ", 1000 * oBudget->snowpack / area);
  if(oControl->sw_BC){
    printf("\n%-31s%-15.2f","Total Surface BC Inflow (mm): ",1000 * oBudget->ovlndinflow / area);
    printf("%-31s%-15.2f","Total GW BC Inflow (mm): ",1000 * oBudget->grndinflow / area);
    if(oControl->sw_deepGW)
      printf("\n%-31s%-15.2f","Total DeepGW BC Inflow (mm): ",1000 * oBudget->deepgrndinflow / area);
  }
  printf("\n%-31s%-15.2f","Canopy Storage (mm): ", 1000 * oBudget->canopy /area);
  printf("%-31s%-15.2f","Ponding (mm): ", 1000 * oBudget->ponding / area);
  printf("\n%-31s%-15.2f","Channel Storage (mm): ", 1000 * oBudget->chan_store / area);
  printf("%-31s%-15.2f","Saturated area fraction: ", oBudget->satarea);
  printf("\n%-31s%-15.2f","Soil water (mm): ", 1000 * oBudget->vadose / area);
  printf("%-31s%-15.2f","of which RootZone water (mm): ", 1000 * oBudget->rootzone / area);
  printf("\n%-31s%-15s%-31s%-15.2f"," "," ","and Groundwater (mm): ", 1000 * oBudget->grndwater / area);
  printf("\n%-31s%-15.2f","Total Evapotranspiration (mm): ", 1000 * oBudget->evaporation / area);
  printf("%-31s%-15.2f","Total Transpiration (mm): ", 1000 * oBudget->transpiration / area);
  printf("\n%-31s%-15.2f","Total Canopy Evaporation (mm): ", 1000 * oBudget->evaporationI / area);
  printf("%-31s%-15.2f","Total Soil Evaporation (mm): ", 1000 * oBudget->evaporationS /area);
  printf("\n%-31s%-15.2f","Total Channel Evaporation (mm): ", 1000 * oBudget->evaporationC /area);
  printf("\n%-31s%-15.2f","Bedrock Leak (mm): ", 1000 * oBudget->leakage /area);  
  if(oControl->sw_deepGW)
    printf("%-31s%-15.2f","Deep GW (mm): ", 1000 * oBudget->deep_grndwater / area);
  printf("\n%-31s%-15.2f","Total OvlndFlow output (mm): ", 1000 * oBudget->ovlndflow / area);  
  printf("%-31s%-15.2f","Total GWFlow output (mm): ", 1000 * oBudget->gwtrflow / area);
  printf("\n%-31s%-15.2f","Run-off to channel (mm): ", 1000 * oBudget->srftochn / area);  
  printf("%-31s%-15.2f","GW to channel (mm): ", 1000 * oBudget->gwtochn / area);  
  if(oControl->sw_deepGW){
    printf("\n%-31s%-15.2f","Total Deep GWFlow output (mm): ", 1000 * oBudget->deepgwtrflow / area);
    printf("%-31s%-15.2f","Deep GW to channel (mm): ", 1000 * oBudget->deepgwtochn / area);
  }
  printf("\n%-31s%-15.2f","GW recharge (mm): ", 1000 * oBudget->recharge / area);
  printf("\n%-31s%-15e","Mass Balance Error (%): ", oBudget->MBErr);
  if(oControl->sw_trck){
    if(oControl->sw_2H)    
      printf("\n%-31s%-15e","Deuterium Mass Balance Error (%%): %e \n", oBudget->MBErr_d2H);
    if(oControl->sw_18O)
      printf("\n%-31s%-15e","Oxygen 18 Mass Balance Error (%%): %e \n", oBudget->MBErr_d18O);
    if(oControl->sw_Age)
      printf("\n%-31s%-15e","Age Mass Balance Error (%%): %e \n", oBudget->MBErr_Age);
  }
  printf("\n");
  // ==== BasinSummary.txt --------------------------------------------------------
  // -----------------------------------------------------
  // m3 is kept for BasinSummary.txt

  ofSummary << oBudget->precipitation << "\t";
  ofSummary << oBudget->snowpack << "\t";
  if(oControl->sw_BC){
    ofSummary << oBudget->ovlndinflow << "\t";
    ofSummary << oBudget->grndinflow << "\t";
    if(oControl->sw_deepGW)
      ofSummary << oBudget->deepgrndinflow << "\t";
  }
  ofSummary << oBudget->canopy << "\t";
  ofSummary << oBudget->ponding << "\t";
  ofSummary << oBudget->chan_store << "\t";
  ofSummary << oBudget->vadose << "\t";         //all soil water
  ofSummary << oBudget->soilL1 << "\t";         //layer 1 soil
  ofSummary << oBudget->soilL2 << "\t";         //layer 2 soil
  ofSummary << oBudget->soilL3 << "\t";
  ofSummary << oBudget->rootzone << "\t";
  ofSummary << oBudget->grndwater << "\t";
  ofSummary << oBudget->evaporation << "\t";    //total ET
  ofSummary << oBudget->evaporationS << "\t";   //total Es
  ofSummary << oBudget->evaporationC << "\t";   //total Echannel
  ofSummary << oBudget->evaporationI << "\t";   //total Einter
  ofSummary << oBudget->transpiration << "\t";  //total transp
  ofSummary << oBudget->leakage << "\t";        //Bedrock leakage
  if(oControl->sw_deepGW)    
    ofSummary << oBudget->deep_grndwater << "\t"; //Deep GW
  ofSummary << oBudget->ovlndflow << "\t";      //overland flow
  ofSummary << oBudget->gwtrflow << "\t";       //GW flow
  if(oControl->sw_deepGW)    
    ofSummary << oBudget->deepgwtrflow << "\t";   //Deep GW flow
  ofSummary << oBudget->srftochn << "\t";       //surface to channel
  ofSummary << oBudget->gwtochn << "\t";        //GW to channel
  if(oControl->sw_deepGW)  
    ofSummary << oBudget->deepgwtochn << "\t";    //Deep GW to channel
  ofSummary << oBudget->recharge << "\t";
  ofSummary << oBudget->satarea << "\t";
  ofSummary << oBudget->MBErr ;
  if(oControl->sw_trck){
    if(oControl->sw_2H)    
      ofSummary << "\t" << oBudget->MBErr_d2H ;
    if(oControl->sw_18O)
      ofSummary << "\t" << oBudget->MBErr_d18O ;
    if(oControl->sw_Age)
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
