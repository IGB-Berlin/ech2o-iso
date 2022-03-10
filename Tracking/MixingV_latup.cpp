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
 * MixingV_latup.cpp
 *
 *  Created on: Apr 25, 2018
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"

void Tracking::MixingV_latup(Basin &bsn, Control &ctrl, 
			     double &d1, double &d2, double &d3, double &fc,
			     double &Qk1, double &dtdx, double &dx, int r, int c)
{
  int mixmod = ctrl.toggle_mix;

  //****************************************************************************
  // Soil state before routing
  //****************************************************************************
  double pond_old = bsn.getPondingWater()->matrix[r][c];
  double theta1_old = bsn.getSoilMoist1()->matrix[r][c];
  double theta2_old = bsn.getSoilMoist2()->matrix[r][c];
  double theta3_old = bsn.getSoilMoist3()->matrix[r][c];

  double S1 = std::max<double>((theta1_old - bsn.getSoilMoistR1()->matrix[r][c])/
			       (bsn.getPorosityL1()->matrix[r][c] - bsn.getSoilMoistR1()->matrix[r][c]),0);
  double S2 = std::max<double>((theta2_old - bsn.getSoilMoistR2()->matrix[r][c])/
			       (bsn.getPorosityL2()->matrix[r][c] - bsn.getSoilMoistR2()->matrix[r][c]),0);
  double S3 = std::max<double>((theta3_old - bsn.getSoilMoistR3()->matrix[r][c])/
			       (bsn.getPorosityL3()->matrix[r][c] - bsn.getSoilMoistR3()->matrix[r][c]),0);

  //****************************************************************************
  // Vertical/Lateral fluxes
  //****************************************************************************
  double L1toSrf = bsn.getFluxExfilt()->matrix[r][c]; 			// L1 to surface
  double L2toL1 = bsn.getFluxL2toL1()->matrix[r][c]; 			// L2 to L1
  double L3toL2 = bsn.getFluxL3toL2()->matrix[r][c]; 			// L3 to L2
  // Lateral out
  double GWtoLat = bsn.getFluxGWtoLat()->matrix[r][c]; 			// GW to next cell
  double ChntoLat = Qk1*dtdx/dx; 					// Channel to next cell
  // Lateral in
  double LattoGW = bsn.getFluxLattoGW()->matrix[r][c];  		// GW in from u/s cell
  double LattoChn = bsn.getFluxLattoChn()->matrix[r][c];  		// Chan in from u/s cell
  double GWtoChn = bsn.getFluxGWtoChn()->matrix[r][c]; 			// GW to channel

  //----------------------------------------------------------------------------
  // Deep groundwater additions
  //----------------------------------------------------------------------------
  double DeepGWtoChn,DeepGWtoLat,LattoDeepGW,DeepGW_old,Leakage;
  double DeepGW_d2H,DeepGW_d18O,DeepGW_Age;
  DeepGWtoChn=DeepGWtoLat=LattoDeepGW=DeepGW_old=Leakage=0;
  DeepGW_d2H=DeepGW_d18O=DeepGW_Age=0;
  if(ctrl.sw_deepGW){
    DeepGWtoChn = bsn.getFluxDeepGWtoChn()->matrix[r][c]; 		// DeepGW to channel
    DeepGWtoLat = bsn.getFluxDeepGWtoLat()->matrix[r][c]; 		// DeepGW to next cell
    LattoDeepGW = bsn.getFluxLattoDeepGW()->matrix[r][c]; 		// DeepGW in from u/s cell
    DeepGW_old = bsn.getDeepGW()->matrix[r][c];
    Leakage = bsn.getBedrockLeakage()->matrix[r][c] * ctrl.dt; 
  }
  //----------------------------------------------------------------------------

  //****************************************************************************
  // For GW and surface (pond+channel), equivalent lateral inputs values
  // Lateral ponded water inflow and snowmelt in MixingV_down
  //****************************************************************************
  double FinSrf = L1toSrf + pond_old;
  double d2Hsurf,d18Osurf,Agesurf;
  double d2Hin,d18Oin,Agein;
  d2Hin=d18Oin=Agein=d2Hsurf=d18Osurf=Agesurf=0;
  
  // for channel mixing

  double chan_store_old = bsn.getChanStoreOld()->matrix[r][c];
  double FinSrfChn = (L1toSrf + pond_old) + GWtoChn + LattoChn + DeepGWtoChn; // all channel inflow

  // Two-pore stuff
  double theta_MW1,theta_MW2,L3toTB2,L3toMW2,MW2toTB1,MW2toMW1;
  theta_MW1=theta_MW2=L3toTB2=L3toMW2=MW2toTB1=MW2toMW1=0;
  
  // save the outflow tracer values of the fluxes
  double L2toL1_d2,L2toL1_o18,L2toL1_age;
  double L1toSrf_d2,L1toSrf_o18,L1toSrf_age;
  double Srfout_d2,Srfout_o18,Srfout_age;
  double Chnout_d2,Chnout_o18,Chnout_age;
  L2toL1_d2=L2toL1_o18=L2toL1_age=-1000;
  L1toSrf_d2=L1toSrf_o18=L1toSrf_age=-1000;
  Srfout_d2=Srfout_o18=Srfout_age=-1000;
  Chnout_d2=Chnout_o18=Chnout_age=-1000;
  
  if(ctrl.sw_TPD){
    theta_MW1 = bsn.getMoistureMW1()->matrix[r][c];
    theta_MW2 = bsn.getMoistureMW2()->matrix[r][c];
    // Return flow to L2: weighted between TB2 (if there's deficit there) and MW2
    L3toTB2 = std::min<double>(L3toL2,//*(theta_MW-theta_r)/(porosity-theta_r),
				std::max<double>(0,d2*(theta_MW2-theta2_old)));
    L3toMW2 = std::max<double>(0, L3toL2 - L3toTB2);
    // Return flow to L1: : weighted between TB1 (if there's deficit there) and MW1
    MW2toTB1 = std::min<double>(L2toL1,//*(theta_MW-theta_r)/(porosity-theta_r),
				std::max<double>(0,d1*(theta_MW1-theta1_old)));
    MW2toMW1 = std::max<double>(0,L2toL1 - MW2toTB1);
  }

  //*******************************************************************************************
  //--- Deep groundwater ----------------------------------------------------------------------
  //*******************************************************************************************
  if(ctrl.sw_deepGW){
    if(LattoDeepGW>RNDOFFERR){
      if(ctrl.sw_2H){
        d2Hin = (_d2Hleakage->matrix[r][c]*Leakage +  _Fd2HLattoDeepGW->matrix[r][c]) / (LattoDeepGW + Leakage);
        TracerMixing(bsn,ctrl,DeepGW_old,_d2HdeepGW->matrix[r][c],_d2HdeepGW->matrix[r][c],
	  	     LattoDeepGW,d2Hin,DeepGWtoLat,_d2HdeepGWQ->matrix[r][c],1.0,mixmod,r,c);
	DeepGW_d2H = _d2HdeepGWQ->matrix[r][c]; 	//for channel mixing
      }
      if(ctrl.sw_18O){
        d18Oin = (_d18Oleakage->matrix[r][c]*Leakage + _Fd18OLattoDeepGW->matrix[r][c]) / (LattoDeepGW + Leakage);
        TracerMixing(bsn,ctrl,DeepGW_old,_d18OdeepGW->matrix[r][c],_d18OdeepGW->matrix[r][c],
	  	     LattoDeepGW,d18Oin,DeepGWtoLat,_d18OdeepGWQ->matrix[r][c],1.0,mixmod,r,c);
	DeepGW_d18O = _d18OdeepGWQ->matrix[r][c]; 	//for channel mixing
      }
      if(ctrl.sw_Age){
        Agein = (_Ageleakage->matrix[r][c]*Leakage + _FAgeLattoDeepGW->matrix[r][c]) / (LattoDeepGW + Leakage);
        TracerMixing(bsn,ctrl,DeepGW_old,_AgedeepGW->matrix[r][c],_AgedeepGW->matrix[r][c],
	  	     LattoDeepGW,Agein,DeepGWtoLat,_AgedeepGWQ->matrix[r][c],1.0,mixmod,r,c);
	DeepGW_Age = _AgedeepGWQ->matrix[r][c]; 	//for channel mixing
      }
    } // end LattoDeepGW
  } // end deepGW mixing

  //*******************************************************************************************
  //--- Layer 3 (GW included) -----------------------------------------------------------------
  //--- Equivalent input signature for GW -----------------------------------------------------
  //*******************************************************************************************
  if(LattoGW > RNDOFFERR){  
    if(ctrl.sw_2H){
      d2Hin = _Fd2HLattoGW->matrix[r][c] / LattoGW ; 		//Lateral GW d2H in
      TracerMixing(bsn,ctrl,theta3_old*d3,_d2Hsoil3->matrix[r][c],_d2Hsoil3->matrix[r][c],
		   LattoGW,d2Hin,(L3toL2+GWtoLat),_d2Hgroundwater->matrix[r][c],S3,mixmod,r,c);
      _d2HGWtoChn->matrix[r][c] = _d2Hgroundwater->matrix[r][c];
    } // end deuterium
    if(ctrl.sw_18O){
      d18Oin = _Fd18OLattoGW->matrix[r][c] / LattoGW ;
      TracerMixing(bsn,ctrl,theta3_old*d3,_d18Osoil3->matrix[r][c],_d18Osoil3->matrix[r][c],
		   LattoGW,d18Oin,(L3toL2+GWtoLat),_d18Ogroundwater->matrix[r][c],S3,mixmod,r,c);
      _d18OGWtoChn->matrix[r][c] = _d18Ogroundwater->matrix[r][c];
    } // end oyxgen-18
    if(ctrl.sw_Age){
      Agein = _FAgeLattoGW->matrix[r][c] / LattoGW ;
      TracerMixing(bsn,ctrl,theta3_old*d3,_Agesoil3->matrix[r][c],_Agesoil3->matrix[r][c],
		   LattoGW,Agein,(L3toL2+GWtoLat),_Agegroundwater->matrix[r][c],S3,mixmod,r,c);
      _AgeGWtoChn->matrix[r][c] = _Agegroundwater->matrix[r][c];
    } // end age
  } // end gw mixing

  //*******************************************************************************************
  //--- Layer 2 -------------------------------------------------------------------------------
  //*******************************************************************************************
  //+++ If two-pore domain activated ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if(ctrl.sw_TPD and L3toL2 > RNDOFFERR){
    // Tightly-bound
    if(ctrl.sw_2H){
      _d2H_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2,_d2H_TB2->matrix[r][c],
					L3toTB2, _d2Hgroundwater->matrix[r][c]);
      if(std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR){
	TracerMixing(bsn,ctrl,theta2_old*d2,_d2H_MW2->matrix[r][c],_d2H_MW2->matrix[r][c],
		     L3toMW2,_d2Hgroundwater->matrix[r][c],L2toL1,L2toL1_d2,S2,mixmod,r,c);
      } else {
	_d2H_MW2->matrix[r][c] = _d2Hgroundwater->matrix[r][c];
      }
    } //end d2H
    if(ctrl.sw_18O){
      _d18O_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2,_d18O_TB2->matrix[r][c],
					 L3toTB2, _d18Ogroundwater->matrix[r][c]);
      if(std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR){
	TracerMixing(bsn,ctrl,theta2_old*d2,_d18O_MW2->matrix[r][c],_d18O_MW2->matrix[r][c],
		     L3toMW2,_d18Ogroundwater->matrix[r][c],L2toL1,L2toL1_o18,S2,mixmod,r,c);
      } else {
	_d18O_MW2->matrix[r][c] = _d18Ogroundwater->matrix[r][c];
      }
    } //end d18O
    if(ctrl.sw_Age){
      _Age_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2,_Age_TB2->matrix[r][c],
					L3toTB2, _Agegroundwater->matrix[r][c]);
      if(std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR){
	TracerMixing(bsn,ctrl,theta2_old*d2,_Age_MW2->matrix[r][c],_Age_MW2->matrix[r][c],
		     L3toMW2,_Agegroundwater->matrix[r][c],L2toL1,L2toL1_age,S2,mixmod,r,c);
      } else {
	_Age_MW2->matrix[r][c] = _Agegroundwater->matrix[r][c];
      }
    }//end age
  //+++ Soil-averaged values +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  } else if (L3toL2 > RNDOFFERR) { // Soil-averaged values
    if(ctrl.sw_2H)
      TracerMixing(bsn,ctrl,theta2_old*d2,_d2Hsoil2->matrix[r][c],_d2Hsoil2->matrix[r][c],
		   L3toL2,_d2Hgroundwater->matrix[r][c],L2toL1,L2toL1_d2,S2,mixmod,r,c);
    if(ctrl.sw_18O)
      TracerMixing(bsn,ctrl,theta2_old*d2,_d18Osoil2->matrix[r][c],_d18Osoil2->matrix[r][c],
		   L3toL2,_d18Ogroundwater->matrix[r][c],L2toL1,L2toL1_o18,S2,mixmod,r,c);
    if(ctrl.sw_Age)
      TracerMixing(bsn,ctrl,theta2_old*d2,_Agesoil2->matrix[r][c],_Agesoil2->matrix[r][c],
		   L3toL2,_Agegroundwater->matrix[r][c],L2toL1,L2toL1_age,S2,mixmod,r,c);
  } // end soil averaged

  //*******************************************************************************************
  //--- Layer 1 -------------------------------------------------------------------------------
  //*******************************************************************************************
  //+++ If two-pore domain activated ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if(ctrl.sw_TPD and L2toL1 > RNDOFFERR){
    // Tightly-bound
    if(ctrl.sw_2H){
      _d2H_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1,_d2H_TB1->matrix[r][c],	
					MW2toTB1, _d2H_MW2->matrix[r][c]);
      if(std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR){
	L2toL1_d2 = L2toL1_d2 == -1000 ? _d2Hsoil2->matrix[r][c] : L2toL1_d2;
	TracerMixing(bsn,ctrl,theta1_old*d1,_d2H_MW1->matrix[r][c],_d2H_MW1->matrix[r][c],
		     MW2toMW1,L2toL1_d2,L1toSrf,L1toSrf_d2,S1,mixmod,r,c);
      } else {
	_d2H_MW1->matrix[r][c] = _d2H_MW1->matrix[r][c];
      }
    } // end d2H
    if(ctrl.sw_18O){
      _d18O_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1,_d18O_TB1->matrix[r][c],
					 MW2toTB1, _d18O_MW2->matrix[r][c]);
      if(std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR){
	L2toL1_o18 = L2toL1_o18 == -1000 ? _d18Osoil2->matrix[r][c] : L2toL1_o18;
	TracerMixing(bsn,ctrl,theta1_old*d1,_d18O_MW1->matrix[r][c],_d18O_MW1->matrix[r][c],
		     MW2toMW1,L2toL1_o18,L1toSrf,L1toSrf_o18,S1,mixmod,r,c);
      } else {
	_d18O_MW1->matrix[r][c] = _d18O_MW1->matrix[r][c];
      }
    } // end d18O
    if(ctrl.sw_Age){
      _Age_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1,_Age_TB1->matrix[r][c],
					MW2toTB1, _Age_MW2->matrix[r][c]);
      if(std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR){
	L2toL1_age = L2toL1_age == -1000 ? _Agesoil2->matrix[r][c] : L2toL1_age;
	TracerMixing(bsn,ctrl,theta1_old*d1,_Age_MW1->matrix[r][c],_Age_MW1->matrix[r][c],
		     MW2toMW1,L2toL1_age,L1toSrf,L1toSrf_age,S1,mixmod,r,c);
      } else {
	_Age_MW1->matrix[r][c] = _Age_MW1->matrix[r][c];
      }
    } // end age
    //+++ Soil-averaged values +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  } else if (L2toL1 > RNDOFFERR) { // Soil-averaged values    
    if(ctrl.sw_2H){
      L2toL1_d2 = L2toL1_d2 == -1000 ? _d2Hsoil2->matrix[r][c] : L2toL1_d2;
      TracerMixing(bsn,ctrl,theta1_old*d1,_d2Hsoil1->matrix[r][c],_d2Hsoil1->matrix[r][c],
		   L2toL1,L2toL1_d2,L1toSrf,L1toSrf_d2,S1,mixmod,r,c);
    }
    if(ctrl.sw_18O){
      L2toL1_o18 = L2toL1_o18 == -1000 ? _d18Osoil2->matrix[r][c] : L2toL1_o18;
      TracerMixing(bsn,ctrl,theta1_old*d1,_d18Osoil1->matrix[r][c],_d18Osoil1->matrix[r][c],
		   L2toL1,L2toL1_o18,L1toSrf,L1toSrf_o18,S1,mixmod,r,c);
    }
    if(ctrl.sw_Age){
      L2toL1_age = L2toL1_age == -1000 ? _Agesoil2->matrix[r][c] : L2toL1_age;
      TracerMixing(bsn,ctrl,theta1_old*d1,_Agesoil1->matrix[r][c],_Agesoil1->matrix[r][c],
		   L2toL1,L2toL1_age,L1toSrf,L1toSrf_age,S1,mixmod,r,c);
    }
  }

  //*******************************************************************************************
  //--- Surface -------------------------------------------------------------------------------
  // If non-channel, return flow + runoff (all water) (runon + snow done in MixingV_Down).
  // If channel cell add input discharge, seepage, river outflow.
  // If surface storage initially null, propagated signature is that of inputs
  // Order:
  //    1)      Get exfiltration tracer
  //    2)      Mix exfiltrated and ponded water (ponded already includes lateral ponded additions)
  //    3)      Get tracer of the input to be mixed (pond + GW + DeepGW + Chan_in)
  //    4)      Mix new input with existing channel storage to get d2H for Qk1
  //******************************************************************************************* 
  if(FinSrfChn > RNDOFFERR) { //if there is water incoming to channel
    if(ctrl.sw_2H){
      d2Hsurf = L1toSrf_d2==-1000 ? (ctrl.sw_TPD ? _d2H_MW1->matrix[r][c] : _d2Hsoil1->matrix[r][c]) : L1toSrf_d2;
      TracerMixing(bsn,ctrl,pond_old,_d2Hsurface->matrix[r][c],_d2Hsurface->matrix[r][c],
 			 L1toSrf,d2Hsurf,FinSrf,Srfout_d2,1.0,0,r,c);
      if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){	  
	_d2HSrftoChn->matrix[r][c] = _d2Hsurface->matrix[r][c];
	_d2HGWtoChn->matrix[r][c] = _d2Hgroundwater->matrix[r][c];
        d2Hin= (FinSrf*_d2Hsurface->matrix[r][c] + GWtoChn*_d2Hgroundwater->matrix[r][c] + 
		         DeepGWtoChn*DeepGW_d2H + _Fd2HLattoChn->matrix[r][c]) / FinSrfChn;
        TracerMixing(bsn,ctrl,chan_store_old,_d2Hchan->matrix[r][c],_d2Hchan->matrix[r][c],
			 FinSrfChn,d2Hin,ChntoLat,Chnout_d2,1.0,0,r,c);
      }// close if channel mixing
    } // close d2H if statement
    if(ctrl.sw_18O){
      d18Osurf = L1toSrf_o18==-1000 ? (ctrl.sw_TPD ? _d18O_MW1->matrix[r][c] : _d18Osoil1->matrix[r][c]): L1toSrf_o18;
      TracerMixing(bsn,ctrl,pond_old,_d18Osurface->matrix[r][c],_d18Osurface->matrix[r][c],
			L1toSrf,d18Osurf,FinSrf,Srfout_o18,1.0,0,r,c);
      if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){	
        _d18OSrftoChn->matrix[r][c] = _d18Osurface->matrix[r][c];
	_d18OGWtoChn->matrix[r][c] = _d18Ogroundwater->matrix[r][c];
	d18Oin= (FinSrf*_d18Osurface->matrix[r][c] + GWtoChn*_d18Ogroundwater->matrix[r][c] + 
		   	DeepGWtoChn*DeepGW_d18O + _Fd18OLattoChn->matrix[r][c]) / FinSrfChn;
        TracerMixing(bsn,ctrl,chan_store_old,_d18Ochan->matrix[r][c],_d18Ochan->matrix[r][c],
			FinSrfChn,d18Oin,ChntoLat,Chnout_o18,1.0,0,r,c);
      }// close if channel mixing
    } // close d18O if statement
    if(ctrl.sw_Age){
      Agesurf = L1toSrf_age == -1000 ? (ctrl.sw_TPD ? _Age_MW1->matrix[r][c] : _Agesoil1->matrix[r][c]): L1toSrf_age;
      TracerMixing(bsn,ctrl,pond_old,_Agesurface->matrix[r][c],_Agesurface->matrix[r][c],
			 L1toSrf,Agesurf,FinSrf,Srfout_age,1.0,0,r,c);
      if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){	
	  _AgeSrftoChn->matrix[r][c] = _Agesurface->matrix[r][c];
	  _AgeGWtoChn->matrix[r][c] = _Agegroundwater->matrix[r][c];
	  Agein= (FinSrf*_Agesurface->matrix[r][c] + GWtoChn*_Agegroundwater->matrix[r][c] + 
		         DeepGWtoChn*DeepGW_Age + _FAgeLattoChn->matrix[r][c]) / FinSrfChn;
          TracerMixing(bsn,ctrl,chan_store_old,_Agechan->matrix[r][c],_Agechan->matrix[r][c],
			 FinSrfChn,Agein,ChntoLat,Chnout_age,1.0,0,r,c);
      }// close if channel mixing
    } // close Age if statement
  } // close FinSrfChn incoming water statement

} // end MixingV_latup

