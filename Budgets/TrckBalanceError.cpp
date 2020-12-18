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
 * MassBalanceError.cpp
 *
 *  Created on: Mar 18, 2010
 *      Author: Marco Maneta
 */

#include "Budget.h"

void Budget::TrckBalanceError(const Control *ctrl)
{
	double inputs = 0.0;
	double outputs = 0.0;
	double Snew = 0.0;
	double Sold = 0.0;

	// Tracking -------------
	double inputs_d2H = 0.0;
	double outputs_d2H = 0.0;
	double Sold_d2H = 0.0;
	double Snew_d2H = 0.0;
	double inputs_d18O = 0.0;
	double outputs_d18O = 0.0;
	double Sold_d18O = 0.0;
	double Snew_d18O = 0.0;
	double inputs_Age = 0.0;
	double outputs_Age = 0.0;
	double Sold_Age = 0.0;
	double Snew_Age = 0.0;

	/* inputs = precipitation + initsnowpack + initponding + initvadose;
	// + initgravwater + initgrndwater; -> obsolete: gravity and groundwater embedded in vadose...	
	outputs = evaporation + ovlndflow + gwtrflow + leakage;
	
	ds = canopy + snowpack + ponding + vadose;*/
	//+ gravwater + grndwater;  -> obsolete: gravity and groundwater embedded in vadose...

	inputs = precipitation;
	


	Sold = initsnowpack + initponding + initchan + initL1 + initL2 + initL3 + initGW;

	outputs = evaporationS + evaporationC + evaporationI + transpiration + ovlndflow + gwtrflow + leakage;

	Snew = canopy + snowpack + ponding + chan_store + soilL1 + soilL2 + soilL3 + grndwater;
	/*
	cout << "inputs: " << inputs << endl;	
	cout << "Sold: " << Sold << endl;
	cout << "Sold--- " << "Snowpack:" << initsnowpack << " | Ponding:" << initponding <<endl; 
	cout << "Sold--- " << "Layer1  :" << initL1       << " | Layer2 :" << initL2 << endl;
	cout << "Sold--- " << "Layer3  :" << initL3       << " | GW     :" << initGW << endl << endl;
	cout << "output: " << outputs << endl;
    	cout << "output--- " << "EvapS :" << evaporationS << " | EvapI :" << evaporationI <<endl; 
	cout << "output--- " << "Trans :" << transpiration<< " | Ovrlnd:" << ovlndflow << endl;
	cout << "output--- " << "GWflo :" << gwtrflow     << " | Leak  :" << leakage << endl;
	cout << "output--- " << "ChanE :" << evaporationC << endl << endl;	
	cout << "Snew: " << Snew << endl;
  	cout << "Snew--- " << "Snowpack:" << snowpack << " | Canopy :" << canopy <<endl; 
	cout << "Snew--- " << "Ponding :" << ponding  << " | Layer1 :" << soilL1 << endl;
	cout << "Snew--- " << "Layer2  :" << soilL2   << " | Layer3 :" << soilL3 << " | GW     :" << grndwater << endl << endl;	  
	*/
		
	if(inputs+Sold>0) 
	  MBErr = 100/(inputs+Sold)*(inputs-outputs +Sold-Snew);
	else 
	  MBErr = 0;

	//Deuterium mass balance
	if(ctrl->sw_2H){
	  inputs_d2H = precipitation_d2H ;

	  Sold_d2H = initsnowpack_d2H + initponding_d2H + initchan_d2H + initL1_d2H + initL2_d2H + initL3_d2H + initGW_d2H;
	  
	  outputs_d2H = evaporationS_d2H + evaporationC_d2H + evaporationI_d2H + transpiration_d2H + ovlndflow_d2H + gwtrflow_d2H + leakage_d2H;
	  
	  Snew_d2H = canopy_d2H + snowpack_d2H + ponding_d2H + chan_d2H + soilL1_d2H + soilL2_d2H + soilL3_d2H + grndwater_d2H;
	  /*
	  cout << "Input d2H: " << inputs_d2H << endl;
	  cout << "Sold: " << Sold_d2H << endl;
	  cout << "Sold d2H--- " << "Snowpack(d2H):" << initsnowpack_d2H << " | Ponding(d2H):" << initponding_d2H <<endl; 
	  cout << "Sold d2H--- " << "Layer1(d2H)  :" << initL1_d2H       << " | Layer2(d2H) :" << initL2_d2H << endl;
	  cout << "Sold d2H--- " << "Layer3(d2H)  :" << initL3_d2H       << " | GW(d2H)     :" << initGW_d2H << endl;
	  cout << "Sold d2H--- " << "Chan(d2H)  :" << initchan_d2H << endl << endl;
	  cout << "output d2H: " << outputs_d2H << endl;
	  cout << "output d2H--- " << "EvapS(d2H) :" << evaporationS_d2H << " | EvapI(d2H) :" << evaporationI_d2H <<endl; 
	  cout << "output d2H--- " << "Trans(d2H) :" << transpiration_d2H<< " | Ovrlnd(d2H):" << ovlndflow_d2H << endl;
	  cout << "output d2H--- " << "GWflo(d2H) :" << gwtrflow_d2H     << " | Leak(d2H)  :" << leakage_d2H << endl;
	  cout << "output d2H--- " << "ChanE      :" << evaporationC_d2H << endl << endl;
	  cout << "Snew: " << Snew_d2H << endl;
	  cout << "Snew d2H--- " << "Snowpack(d2H):" << snowpack_d2H << " | Canopy(d2H) :" << canopy_d2H <<endl; 
	  cout << "Snew d2H--- " << "Ponding(d2H) :" << ponding_d2H  << " | Layer1(d2H) :" << soilL1_d2H << endl;
	  cout << "Snew d2H--- " << "Layer2(d2H)  :" << soilL2_d2H   << " | Layer3(d2H) :" << soilL3_d2H << " | GW(d2H)     :" << grndwater_d2H << endl;	  
	  cout << "Sold d2H--- " << "Chan(d2H)  :" << chan_d2H << endl << endl;
	  */
	  
	  if(inputs+Sold>0) 
	    MBErr_d2H = 100*(inputs_d2H-outputs_d2H + Sold_d2H-Snew_d2H)/(inputs_d2H+Sold_d2H);
	  else 
	    MBErr_d2H = 0;
	  
	}
	
	// Oxygen 18 mass balance
	if(ctrl->sw_18O){
	  inputs_d18O = precipitation_d18O ;

	  Sold_d18O = initsnowpack_d18O + initponding_d18O + initchan_d18O + initL1_d18O + initL2_d18O + initL3_d18O + initGW_d18O;
	  
	  outputs_d18O = evaporationS_d18O + evaporationC_d18O + evaporationI_d18O + transpiration_d18O + ovlndflow_d18O + gwtrflow_d18O + leakage_d18O;
	  
	  Snew_d18O = canopy_d18O + snowpack_d18O + ponding_d18O + chan_d18O + soilL1_d18O + soilL2_d18O + soilL3_d18O + grndwater_d18O;

	  cout << "Input d18O: " << inputs_d18O << endl;	  
  	  cout << "Sold d18O--- " << "Snowpack(d18O):" << initsnowpack_d18O << " | Ponding(d18O):" << initponding_d18O <<endl; 
	  cout << "Sold d18O--- " << "Layer1(d18O)  :" << initL1_d18O       << " | Layer2(d18O) :" << initL2_d18O << endl;
	  cout << "Sold d18O--- " << "Layer3(d18O)  :" << initL3_d18O       << " | GW(d18O)     :" << initGW_d18O << endl;
	  cout << "Sold d18O--- " << "Chan(d18O)  :" << initchan_d18O << endl << endl;
	  cout << "output d18O: " << outputs_d18O << endl;
	  cout << "output d18O--- " << "EvapS(d18O) :" << evaporationS_d18O << " | EvapI(d18O) :" << evaporationI_d18O <<endl; 
	  cout << "output d18O--- " << "Trans(d18O) :" << transpiration_d18O<< " | Ovrlnd(d18O):" << ovlndflow_d18O << endl;
	  cout << "output d18O--- " << "GWflo(d18O) :" << gwtrflow_d18O     << " | Leak(d18O)  :" << leakage_d18O << endl;
	  cout << "output d18O--- " << "ChanE      :" << evaporationC_d18O << endl << endl;
	  cout << "Snew: " << Snew_d18O << endl;
	  cout << "Snew d18O--- " << "Snowpack(d18O):" << snowpack_d18O << " | Canopy(d18O) :" << canopy_d18O <<endl; 
	  cout << "Snew d18O--- " << "Ponding(d18O) :" << ponding_d18O  << " | Layer1(d18O) :" << soilL1_d18O << endl;
	  cout << "Snew d18O--- " << "Layer2(d18O)  :" << soilL2_d18O   << " | Layer3(d18O) :" << soilL3_d18O << " | GW(d18O)     :" << grndwater_d18O << endl;	  
	  cout << "Sold d18O--- " << "Chan(d18O)  :" << chan_d18O << endl << endl;

	  if(inputs+Sold>0) 
	    MBErr_d18O = 100/(inputs_d18O+Sold_d18O)*(inputs_d18O-outputs_d18O + Sold_d18O-Snew_d18O);
	  else 
	    MBErr_d18O = 0;
	}
	
	// Age mass balance
	if(ctrl->sw_Age){
	  inputs_Age = precipitation_Age ; // gradually aging precip  

	  // Ageing initial storage
	  initsnowpack_Age += initsnowpack * ctrl->dt / 86400 ;
	  initponding_Age += initponding * ctrl->dt / 86400;
	  initchan_Age += initchan * ctrl->dt / 86400;
	  initL1_Age += initL1 * ctrl->dt / 86400 ;
	  initL2_Age += initL2 * ctrl->dt / 86400 ;
	  initL3_Age += initL3 * ctrl->dt / 86400 ;
	  initGW_Age += initGW * ctrl->dt / 86400;	  

	  Sold_Age = initsnowpack_Age + initponding_Age + initchan_Age + 
	    initL1_Age + initL2_Age + initL3_Age + initGW_Age;

	  outputs_Age = evaporationS_Age +evaporationC_Age + evaporationI_Age + transpiration_Age +
	    ovlndflow_Age + gwtrflow_Age + leakage_Age;
	  
	  Snew_Age = canopy_Age + snowpack_Age + ponding_Age + chan_Age + 
	    soilL1_Age + soilL2_Age + soilL3_Age + grndwater_Age;
	  
	  cout << "Input Age: " << inputs_Age << endl;	  
	  cout << "Sold_Age   " << Sold_Age << endl;
  	  cout << "Sold Age--- " << "Snowpack(Age):" << initsnowpack_Age << " | Ponding(Age):" << initponding_Age <<endl; 
	  cout << "Sold Age--- " << "Layer1(Age)  :" << initL1_Age       << " | Layer2(Age) :" << initL2_Age << endl;
	  cout << "Sold Age--- " << "Layer3(Age)  :" << initL3_Age       << " | GW(Age)     :" << initGW_Age << endl;
	  cout << "Sold Age--- " << "Chan(Age)  :" << initchan_Age << endl << endl;
	  cout << "output Age: " << outputs_Age << endl;
	  cout << "output Age--- " << "EvapS(Age) :" << evaporationS_Age << " | EvapI(Age) :" << evaporationI_Age <<endl; 
	  cout << "output Age--- " << "Trans(Age) :" << transpiration_Age<< " | Ovrlnd(Age):" << ovlndflow_Age << endl;
	  cout << "output Age--- " << "GWflo(Age) :" << gwtrflow_Age     << " | Leak(Age)  :" << leakage_Age << endl;
	  cout << "output Age--- " << "ChanE      :" << evaporationC_Age << endl << endl;
	  cout << "Snew: " << Snew_Age << endl;
	  cout << "Snew Age--- " << "Snowpack(Age):" << snowpack_Age << " | Canopy(Age) :" << canopy_Age <<endl; 
	  cout << "Snew Age--- " << "Ponding(Age) :" << ponding_Age  << " | Layer1(Age) :" << soilL1_Age << endl;
	  cout << "Snew Age--- " << "Layer2(Age)  :" << soilL2_Age   << " | Layer3(Age) :" << soilL3_Age << " | GW(Age)     :" << grndwater_Age << endl;	  
	  cout << "Sold Age--- " << "Chan(Age)  :" << chan_Age << endl << endl;

	  if(inputs_Age+Sold_Age > RNDOFFERR) 
	    MBErr_Age = 100/(inputs_Age+Sold_Age)*(inputs_Age-outputs_Age+Sold_Age-Snew_Age);
	  else 
	    MBErr_Age = 0;
	}

	/*	
	// Debugging verbose
	cout << "== Fluxes/Stores ===---------------------------------------------" << endl;
	cout << endl << "precip: " << precipitation << 
	  ", sum(Isnow): " << initsnowpack << 
	  ", sum(Isurf): " << initponding << 
	  //", sum(Ivadose): " << initvadose << endl;
	  ", sum(IsoilL1): " << initL1 << 
	  ", sum(IsoilL2): " << initL2 << 	
	  ", sum(IsoilL3): " << initL3 << 
	  ", sum(IGW): " << initGW << endl;
	
	cout << "sum(evapS): " << evaporationS << 
	  ", sum(evapI): " << evaporationI << 
	  ", sum(evapT): " << transpiration << 
	  ", sum(OutSurf): " << ovlndflow << 
	  ", sum(OutGW): " << gwtrflow << 
	  ", sum(OutLeak): " << leakage << endl;
	
	cout << "sum(canopy): " << canopy << 
	  ", sum(snowpack): " << snowpack << 
	  ", sum(ponding): " << ponding << 
	  //", sum(vadose): " << vadose << endl << endl;
	  ", sum(soilL1): " << soilL1 << 
	  ", sum(soilL2): " << soilL2 << 	
	  ", sum(soilL3): " << soilL3 << 
	  ", sum(GW): " << grndwater << endl;
	
	cout << "sum(inputs): " << inputs << 
	  ", sum(outputs): " << outputs <<
	  ", sum(S) : " << Snew <<endl;

	if(ctrl->sw_2H){
	  cout << "== d2H ===------------------------------------------------------" << endl;
	  // cout << "precip*d2H: " << precipitation_d2H << 
	  //   ", Icanopy*d2H: " << initcanopy_d2H << 
	  //   ", Isnow*d2H: " << initsnowpack_d2H << 
	  //   ", Isurf*d2H: " << initponding_d2H << 
	  //   ", IsoilL1*d2H: " << initL1_d2H << 
	  //   ", IsoilL2*d2H: " << initL2_d2H << 	
	  //   ", IsoilL3*d2H: " << initL3_d2H << 
	  //   ", IGW*d2H: " << initGW_d2H << endl;	  
	  // cout << "evapS*d2H: " << evaporationS_d2H << 
	  //   ", evapI*d2H: " << evaporationI_d2H << 
	  //   ", evapT*d2H: " << transpiration_d2H << 
	  //   ", OutSurf*d2H: " << ovlndflow_d2H << 
	  //   ", OutGW*d2H: " << gwtrflow_d2H << 
	  //   ", OutLeak*d2H: " << leakage_d2H << endl;
	  // cout << "canopy*d2H: " << canopy_d2H << 
	  //   ", snowpack*d2H: " << snowpack_d2H << 
	  //   ", ponding*d2H: " << ponding_d2H << 
	  //   ", soilL1*d2H: " << soilL1_d2H << 
	  //   ", soilL2*d2H: " << soilL2_d2H << 
	  //   ", soilL3*d2H: " << soilL3_d2H << 
	  //   ", GW*d2H: " << grndwater_d2H << endl;
	  
	  cout << "precip: " << precipitation_d2H / precipitation << endl ;
	  cout << //"Icanopy: " << initcanopy_d2H / canopy << 
	    "Isnow: " << initsnowpack_d2H / initsnowpack << 
	    ", Isurf: " << initponding_d2H / initponding << 
	    ", IsoilL1: " << initL1_d2H / initL1 << 
	    ", IsoilL2: " << initL2_d2H / initL2 << 	
	    ", IsoilL3: " << initL3_d2H / initL3 << 
	    ", IGW: " << initGW_d2H / initGW << endl;	  
	  cout << "evapS: " << evaporationS_d2H / evaporationS << 
	    ", evapI: " << evaporationI_d2H / evaporationI << 
	    ", evapT: " << transpiration_d2H / transpiration << 
	    ", OutSurf: " << ovlndflow_d2H / ovlndflow << 
	    ", OutGW: " << gwtrflow_d2H / gwtrflow << 
	    ", OutLeak: " << leakage_d2H / leakage << endl;
	  cout << "canopy: " << canopy_d2H /canopy << 
	    ", snowpack: " << snowpack_d2H / snowpack << 
	    ", ponding: " << ponding_d2H / ponding << 
	    ", soilL1: " << soilL1_d2H / soilL1 << 
	    ", soilL2: " << soilL2_d2H / soilL2 << 
	    ", soilL3: " << soilL3_d2H / soilL3 << 
	    ", GW: " << grndwater_d2H / grndwater << endl ;
	  
	  cout << "inputs: " << inputs_d2H/inputs << ", outputs: " << outputs_d2H/outputs << 
	    ", Sold: " << Sold_d2H/Sold << ", Snew: " << Snew_d2H/Snew <<endl;
	  cout << "inputs*d2H: " << inputs_d2H << ", outputs*d2H: " << 
	    outputs_d2H << ", Sold*d2H: " << Sold_d2H << ", Snew*d2H: " << Snew_d2H <<endl;	  
	}
	
	if(ctrl->sw_18O){
	  cout << "== d18O ===-----------------------------------------------------" << endl;
	  // cout << "precip*d18O: " << precipitation_d18O << 
	  //   ", Icanopy*d18O: " << initcanopy_d18O << 
	  //   ", Isnow*d18O: " << initsnowpack_d18O << 
	  //   ", Isurf*d18O: " << initponding_d18O << 
	  //   ", IsoilL1*d18O: " << initL1_d18O << 
	  //   ", IsoilL2*d18O: " << initL2_d18O << 	
	  //   ", IsoilL3*d18O: " << initL3_d18O << 
	  //   ", IGW*d18O: " << initGW_d18O << endl;
       	  // cout << "evapS*d18O: " << evaporationS_d18O << 
	  //   ", evapI*d18O: " << evaporationI_d18O << 
	  //   ", evapT*d18O: " << transpiration_d18O << 
	  //   ", OutSurf*d18O: " << ovlndflow_d18O << 
	  //   ", OutGW*d18O: " << gwtrflow_d18O << 
	  //   ", OutLeak*d18O: " << leakage_d18O << endl;
	  // cout << "canopy*d18O: " << canopy_d18O << 
	  //   ", snowpack*d18O: " << snowpack_d18O << 
	  //   ", ponding*d18O: " << ponding_d18O << 
	  //   ", soilL1*d18O: " << soilL1_d18O << 
	  //   ", soilL2*d18O: " << soilL2_d18O << 
	  //   ", soilL3*d18O: " << soilL3_d18O << 
	  //   ", GW*d18O: " << grndwater_d18O << endl;
	  cout << "precip: " << precipitation_d18O / precipitation << endl ;
	  cout << //", Icanopy: " << initcanopy_d18O / canopy << 
	    "Isnow: " << initsnowpack_d18O / initsnowpack << 
	    ", Isurf: " << initponding_d18O / initponding << 
	    ", IsoilL1: " << initL1_d18O / initL1 << 
	    ", IsoilL2: " << initL2_d18O / initL2 << 	
	    ", IsoilL3: " << initL3_d18O / initL3 << 
	    ", IGW: " << initGW_d18O / initGW << endl;
       	  cout << "evapS: " << evaporationS_d18O / evaporationS << 
	    ", evapI: " << evaporationI_d18O / evaporationI << 
	    ", evapT: " << transpiration_d18O / transpiration << 
	    ", OutSurf: " << ovlndflow_d18O / ovlndflow << 
	    ", OutGW: " << gwtrflow_d18O / gwtrflow << 
	    ", OutLeak: " << leakage_d18O / leakage << endl;
	  cout << "canopy: " << canopy_d18O / canopy << 
	    ", snowpack: " << snowpack_d18O /snowpack << 
	    ", ponding: " << ponding_d18O / ponding << 
	    ", soilL1: " << soilL1_d18O / soilL1 << 
	    ", soilL2: " << soilL2_d18O / soilL2 << 
	    ", soilL3: " << soilL3_d18O / soilL3 << 
	    ", GW: " << grndwater_d18O / grndwater << endl;

	  cout << "inputs: " << inputs_d18O/inputs << ", outputs: " << outputs_d18O/outputs << 
	    ", Sold: " << Sold_d18O/Sold << ", Snew: " << Snew_d18O/Snew <<endl;
	  cout << "inputs*d18O: " << inputs_d18O << ", outputs*d18O: " << 
	    outputs_d18O << ", Sold*d18O: " << Sold_d18O << ", Snew*d18O: " << Snew_d18O <<endl;	  
	}
   
	if(ctrl->sw_Age){

	  cout << "== AGE ===-----------------------------------------------------" << endl;
	  cout << "precip: " << precipitation_Age/precipitation << endl;
	  cout << ", Isnow*Age: " << snowpack_Age / initsnowpack << 
	    ", Isurf*Age: " << initponding_Age / initponding << 
	    ", IsoilL1*Age: " << initL1_Age / initL1 << 
	    ", IsoilL2*Age: " << initL2_Age / initL2 << 	
	    ", IsoilL3*Age: " << initL3_Age / initL3 << 
	    ", IGW*Age: " << initGW_Age / initGW << endl;
	  
	  // cout << "evapS*Age: " << evaporationS_Age << 
	  //   ", evapI*Age: " << evaporationI_Age << 
	  //   ", evapT*Age: " << transpiration_Age << 
	  //   ", OutSurf*Age: " << ovlndflow_Age << 
	  //   ", OutGW*Age: " << gwtrflow_Age << 
	  //   ", OutLeak*Age: " << leakage_Age << endl;
	  // cout << "canopy*Age: " << canopy_Age << 
	  //   ", snowpack*Age: " << snowpack_Age << 
	  //   ", pondingT*Age: " << ponding_Age << 
	  //   ", soilL1*Age: " << soilL1_Age << 
	  //   ", soilL2*Age: " << soilL2_Age << 
	  //   ", soilL3*Age: " << soilL3_Age << 
	  //   ", GW*Age: " << grndwater_Age << endl;

	  cout << "evapS: " << evaporationS_Age / evaporationS << 
	    ", evapI: " << evaporationI_Age / evaporationI << 
	    ", evapT: " << transpiration_Age / transpiration << 
	    ", OutSurf: " << ovlndflow_Age / ovlndflow << 
	    ", OutGW: " << gwtrflow_Age / gwtrflow_Age << 
	    ", OutLeak: " << leakage_Age / leakage << endl;
	  cout << "canopy: " << canopy_Age/canopy << 
	    ", snowpack: " << snowpack_Age / snowpack << 
	    ", ponding: " << ponding_Age / ponding << 
	    ", soilL1: " << soilL1_Age / soilL1 << 
	    ", soilL2: " << soilL2_Age / soilL2 << 
	    ", soilL3: " << soilL3_Age / soilL3 << 
	    ", GW: " << grndwater_Age / grndwater << endl;
	  	  
	  cout << "inputs: " << inputs_Age/inputs << 
	    ", outputs: " << outputs_Age/outputs <<
	    ", Sold: " << Sold_Age/Sold <<
	    ", Snew: " << Snew_Age/Snew <<endl;
	  cout << "inputs*Age: " << inputs_Age << 
	    ", outputs*Age: " << outputs_Age <<
	    ", Sold*Age: " << Sold_Age <<
	    ", Snew*Age: " << Snew_Age <<endl;
	}
	cout << "-----------------------------------------------------" << endl;	
	*/
	// === Save last time step storage --------------------------------------------------
	/*
	if(ctrl->sw_2H){
	  initcanopy_d2H = canopy_d2H;
	  initsnowpack_d2H = snowpack_d2H;
	  initponding_d2H = ponding_d2H;
	  initL1_d2H = soilL1_d2H;
	  initL2_d2H = soilL2_d2H;
	  initL3_d2H = soilL3_d2H;
	  initGW_d2H = grndwater_d2H;
	}
	if(ctrl->sw_18O){
	  initcanopy_d18O = canopy_d18O;
	  initsnowpack_d18O = snowpack_d18O;
	  initponding_d18O = ponding_d18O;
	  initL1_d18O = soilL1_d18O;
	  initL2_d18O = soilL2_d18O;
	  initL3_d18O = soilL3_d18O;
	  initGW_d18O = grndwater_d18O;
	}
	
	if(ctrl->sw_Age){
	  initcanopy_Age = canopy_Age;
	  initsnowpack_Age = snowpack_Age;
	  initponding_Age = ponding_Age;
	  initL1_Age = soilL1_Age;
	  initL2_Age = soilL2_Age;
	  initL3_Age = soilL3_Age;
	  initGW_Age = grndwater_Age;	  
	}
	*/
	
}
