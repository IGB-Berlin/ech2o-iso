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
  double inputsBC = 0.0;
  double outputs = 0.0;
  double Snew = 0.0;
  double Sold = 0.0;

  // Tracking -------------
  double inputs_d2H = 0.0;
  double inputsBC_d2H = 0.0;
  double outputs_d2H = 0.0;
  double Sold_d2H = 0.0;
  double Snew_d2H = 0.0;
  double inputs_d18O = 0.0;
  double inputsBC_d18O = 0.0;  
  double outputs_d18O = 0.0;
  double Sold_d18O = 0.0;
  double Snew_d18O = 0.0;
  double inputs_Age = 0.0;
  double inputsBC_Age = 0.0;  
  double outputs_Age = 0.0;
  double Sold_Age = 0.0;
  double Snew_Age = 0.0;

  if(ctrl->sw_BC){
    inputsBC = ovlndinflow + grndinflow;
    if(ctrl->sw_deepGW)
      inputsBC += deepgrndinflow;
  }
  
  inputs = precipitation + inputsBC;

  Sold = initsnowpack + initponding + initchan + initL1 + initL2 + initL3 +
              initGW + initDeepGW;
  outputs = evaporationS + evaporationC + evaporationI + transpiration +
    ovlndflow + gwtrflow + (ctrl->sw_deepGW ? 0 : leakage) + deepgwtrflow;

  Snew = canopy + snowpack + ponding + chan_store + soilL1 + soilL2 + soilL3 + grndwater + deep_grndwater;
	
  if(inputs+Sold>0) 
    MBErr = 100/(inputs+Sold)*(inputs-outputs +Sold-Snew);
  else 
    MBErr = 0;

  //-----------------------------------------------------------------------
  //Deuterium mass balance
  //-----------------------------------------------------------------------  
  if(ctrl->sw_2H){

    if(ctrl->sw_BC){
      inputsBC_d2H = ovlndinflow_d2H + grndinflow_d2H;
      if(ctrl->sw_deepGW)
	inputsBC_d2H += deepgrndinflow_d2H;
    }
    inputs_d2H = precipitation_d2H  + inputsBC_d2H;
    
    Sold_d2H = initsnowpack_d2H + initponding_d2H + initchan_d2H + initL1_d2H +
                    initL2_d2H + initL3_d2H + initGW_d2H;
    
    outputs_d2H = evaporationS_d2H + evaporationC_d2H + evaporationI_d2H + transpiration_d2H +
                    ovlndflow_d2H + gwtrflow_d2H + leakage_d2H;
	  
    Snew_d2H = canopy_d2H + snowpack_d2H + ponding_d2H + chan_d2H + soilL1_d2H +
                    soilL2_d2H + soilL3_d2H + grndwater_d2H;
    
    if(inputs+Sold>0) 
      MBErr_d2H = 100*(inputs_d2H-outputs_d2H + Sold_d2H-Snew_d2H)/(inputs_d2H+Sold_d2H);
    else 
      MBErr_d2H = 0;
  }

  //-----------------------------------------------------------------------
  // Oxygen 18 mass balance
  //-----------------------------------------------------------------------  
  if(ctrl->sw_18O){
    if(ctrl->sw_BC){
      inputsBC_d18O = ovlndinflow_d18O + grndinflow_d18O;
      if(ctrl->sw_deepGW)
	inputsBC_d18O += deepgrndinflow_d18O;
    }

    inputs_d18O = precipitation_d18O + inputsBC_d18O;
    
    Sold_d18O = initsnowpack_d18O + initponding_d18O + initchan_d18O + initL1_d18O +
                     initL2_d18O + initL3_d18O + initGW_d18O;
    
    outputs_d18O = evaporationS_d18O + evaporationC_d18O + evaporationI_d18O + transpiration_d18O +
                     ovlndflow_d18O + gwtrflow_d18O + leakage_d18O;
    
    Snew_d18O = canopy_d18O + snowpack_d18O + ponding_d18O + chan_d18O + soilL1_d18O +
                     soilL2_d18O + soilL3_d18O + grndwater_d18O;

    if(inputs+Sold>0) 
      MBErr_d18O = 100/(inputs_d18O+Sold_d18O)*(inputs_d18O-outputs_d18O + Sold_d18O-Snew_d18O);
    else 
      MBErr_d18O = 0;
  }

  //-----------------------------------------------------------------------  
  // Age mass balance
  //-----------------------------------------------------------------------  
  if(ctrl->sw_Age){
    if(ctrl->sw_BC){
      inputsBC_Age = ovlndinflow_Age + grndinflow_Age;
      if(ctrl->sw_deepGW)
	inputsBC_Age += deepgrndinflow_Age;
    }

    inputs_Age = precipitation_Age + inputsBC_Age; // gradually aging precip  
    
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

    if(inputs_Age+Sold_Age > RNDOFFERR) 
      MBErr_Age = 100/(inputs_Age+Sold_Age)*(inputs_Age-outputs_Age+Sold_Age-Snew_Age);
    else 
      MBErr_Age = 0;
  }


  /*
    cout << "inputs: " << inputs << endl;	
    cout << "Sold: " << Sold << endl;
    cout << "Sold--- " << "Snowpack:" << initsnowpack << " | Ponding:" << initponding <<endl; 
    cout << "Sold--- " << "Layer1  :" << initL1       << " | Layer2 :" << initL2 << endl;
    cout << "Sold--- " << "Layer3  :" << initL3       << " | GW     :" << initGW << endl << endl;
    cout << "output: " << outputs << endl;
    cout << "output--- " << "EvapS :" << evaporationS << " | EvapI :" << evaporationI <<endl;
    cout << "output--- " << "EvapC :" << evaporationC << endl;
    cout << "output--- " << "Trans :" << transpiration<< " | Ovrlnd:" << ovlndflow << endl;
    cout << "output--- " << "GWflo :" << gwtrflow     << " | Leak  :" << leakage << endl;
    cout << "output--- " << "Deepgw:" << deepgwtrflow << endl << endl;	
    cout << "Snew: " << Snew << endl;
    cout << "Snew--- " << "Snowpack:" << snowpack << " | Canopy :" << canopy    <<endl; 
    cout << "Snew--- " << "Ponding :" << ponding  << " | Chan   :" << chan_store<< " | Layer1 :" << soilL1 << endl;
    cout << "Snew--- " << "Layer2  :" << soilL2   << " | Layer3 :" << soilL3    << " | GW     :" << grndwater << endl << endl;	  
  */
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

    /*
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
    */

    /*
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
    */
  
}
