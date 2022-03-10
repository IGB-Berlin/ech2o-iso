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
 * CreateWorld.cpp
 *
 *  Created on: Jul 30, 2010
 *      Author: Marco.Maneta
 */

#include "Sativa.h"

int CreateWorld(char* argv[]){

  oControl = new Control;
  cout << "Control created ok... " << "\n";
  oControl->ReadConfigFile(argv[1]);
  cout << "Config.ini read ok... " << "\n";

  if(oControl->toggle_hydrologic_engine>0 and oControl->sw_trck)
    throw std::ios::failure("Error: Richard's not available with Tracking right now\n");
  
  oAtmosphere = new Atmosphere(*oControl);
  cout << "Atmosphere created ok... " << "\n";
  
  oBasin = new Basin(*oControl, *oAtmosphere);
  cout << "Basin created ok... " << "\n";

  oReport = new Report(*oControl);
  cout << "Report created ok... " << "\n";

  oTracking = new Tracking(*oControl, *oBasin, *oAtmosphere);
  if(oControl->sw_trck)
    cout << "Isotope module created ok... " << "\n";

  oBudget = new Budget(oBasin, oControl, oTracking);
  cout << "Budget created ok... " << "\n";

  // == Basin Summary ==========================================================
  // ---
  try{
    ofSummary.open((oControl->path_ResultsFolder + "BasinSummary.txt").c_str());
    if(!ofSummary)
      throw std::ios::failure("Error opening BasinSummary.txt buffer\n");

  }catch(const std::exception &e){
    cout << e.what() << endl;
    throw;
  }
  // Headers for BasinSummary
  ofSummary << "Precip(m3)\t";
  ofSummary << "SWE(m3)\t";
  if(oControl->sw_BC){
    ofSummary << "BCSurf(m3)\t";
    ofSummary << "BCGW(m3)\t";
    if(oControl->sw_deepGW)
      ofSummary << "DeepBCGW(m3)\t";
  }
  ofSummary << "Canopy(m3)\t";
  ofSummary << "Surface(m3)\t";
  ofSummary << "Channel(m3)\t";
  ofSummary << "SoilW(m3)\t";
  ofSummary << "SoilL1(m3)\t";
  ofSummary << "SoilL2(m3)\t";
  ofSummary << "SoilL3(m3)\t";
  ofSummary << "RZW(m3)\t";
  ofSummary << "GW(m3)\t";
  ofSummary << "ET(m3)\t";
  ofSummary << "EvapS(m3)\t";
  ofSummary << "EvapC(m3)\t";  
  ofSummary << "EvapI(m3)\t";
  ofSummary << "EvapT(m3)\t";
  ofSummary << "Leakage(m3)\t";
  if(oControl->sw_deepGW)  
    ofSummary << "DeepGW(m3)\t";
  ofSummary << "SrfOut(m3)\t";
  ofSummary << "GWOut(m3)\t";
  if(oControl->sw_deepGW)    
    ofSummary << "DeepGWOut(m3)\t";  
  ofSummary << "SrftoCh\t";
  ofSummary << "GWtoCh\t";
  if(oControl->sw_deepGW)    
    ofSummary << "DeepGWtoCh\t";  
  ofSummary << "Rchrge\t";
  ofSummary << "SatExt\t";
  ofSummary << "MBErr";
  if(oControl->sw_trck and oControl->sw_2H)
    ofSummary << "\t2H_MBE";
  if(oControl->sw_trck and oControl->sw_18O)
    ofSummary << "\t18O_MBE";
  if(oControl->sw_trck and oControl->sw_Age)
    ofSummary << "\tAge_MBE";
  ofSummary << "\n";

  // == d2H Summary ==========================================================
  // ---  
  if(oControl->sw_trck and oControl->sw_2H){
    try{
      ofd2HSummary.open((oControl->path_ResultsFolder + "Basind2HSummary.txt").c_str());
      if(!ofd2HSummary)
	throw std::ios::failure("Error opening Basind2HSummary.txt buffer\n");
      
    }catch(const std::exception &e){
      cout << e.what() << endl;
      throw;
    }
    // Headers for BasinSummary
    ofd2HSummary << "StorAll\t";
    ofd2HSummary << "Snow\t";
    ofd2HSummary << "Intercp\t";
    ofd2HSummary << "Surface\t";
    ofd2HSummary << "Soil\t";
    ofd2HSummary << "SoilL1\t";
    ofd2HSummary << "SoilL2\t";
    ofd2HSummary << "SoilL3\t";
    ofd2HSummary << "RZW\t";
    ofd2HSummary << "GW\t";
    ofd2HSummary << "ET\t";
    ofd2HSummary << "EvapS\t";
    ofd2HSummary << "EvapI\t";
    ofd2HSummary << "EvapT\t";
    ofd2HSummary << "Leakage\t";
    ofd2HSummary << "SrfOut\t";
    ofd2HSummary << "GWOut\t";
    ofd2HSummary << "AllOut\t";
    ofd2HSummary << "SrftoCh\t";
    ofd2HSummary << "GWtoCh\t";
    ofd2HSummary << "Rchrge\t";
    ofd2HSummary << "Precip\n";
  }

  // == d18O Summary ==========================================================
  // ---  
  if(oControl->sw_trck and oControl->sw_18O){
    try{
      ofd18OSummary.open((oControl->path_ResultsFolder + "Basind18OSummary.txt").c_str());
      if(!ofd18OSummary)
	throw std::ios::failure("Error opening Basind18OSummary.txt buffer\n");
      
    }catch(const std::exception &e){
      cout << e.what() << endl;
      throw;
    }
    // Headers for BasinSummary
    ofd18OSummary << "StorAll\t";
    ofd18OSummary << "Snow\t";
    ofd18OSummary << "Intercp\t";
    ofd18OSummary << "Surface\t";
    ofd18OSummary << "Soil\t";
    ofd18OSummary << "SoilL1\t";
    ofd18OSummary << "SoilL2\t";
    ofd18OSummary << "SoilL3\t";
    ofd18OSummary << "RZW\t";
    ofd18OSummary << "GW\t";
    ofd18OSummary << "ET\t";
    ofd18OSummary << "EvapS\t";
    ofd18OSummary << "EvapI\t";
    ofd18OSummary << "EvapT\t";
    ofd18OSummary << "Leakage\t";
    ofd18OSummary << "SrfOut\t";
    ofd18OSummary << "GWOut\t";
    ofd18OSummary << "AllOut\t";
    ofd18OSummary << "SrftoCh\t";
    ofd18OSummary << "GWtoCh\t";
    ofd18OSummary << "Rchrge\t";
    ofd18OSummary << "Precip\n";
  }

  // == Age Summary ==========================================================
  // ---  
  if(oControl->sw_trck and oControl->sw_Age){
    try{
      ofAgeSummary.open((oControl->path_ResultsFolder + "BasinAgeSummary.txt").c_str());
      if(!ofAgeSummary)
	throw std::ios::failure("Error opening BasinAgeSummary.txt buffer\n");
      
    }catch(const std::exception &e){
      cout << e.what() << endl;
      throw;
    }
    // Headers for BasinSummary
    ofAgeSummary << "StorAll\t";
    ofAgeSummary << "Snow\t";
    ofAgeSummary << "Intercp\t";
    ofAgeSummary << "Surface\t";
    ofAgeSummary << "Soil\t";
    ofAgeSummary << "SoilL1\t";
    ofAgeSummary << "SoilL2\t";
    ofAgeSummary << "SoilL3\t";
    ofAgeSummary << "RZW\t";
    ofAgeSummary << "GW\t";
    ofAgeSummary << "ET\t";
    ofAgeSummary << "EvapS\t";
    ofAgeSummary << "EvapI\t";
    ofAgeSummary << "EvapT\t";
    ofAgeSummary << "Leakage\t";
    ofAgeSummary << "SrfOut\t";
    ofAgeSummary << "GWOut\t";
    ofAgeSummary << "AllOut\t";
    ofAgeSummary << "SrftoCh\t";
    ofAgeSummary << "GWtoCh\t";
    ofAgeSummary << "Rchrge\n";

  }
  
  return EXIT_SUCCESS;
}
