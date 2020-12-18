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

  // Soil state before routing
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

  // Vertical fluxes
  double L1toSrf = bsn.getFluxExfilt()->matrix[r][c];
  double L2toL1 = bsn.getFluxL2toL1()->matrix[r][c];
  double L3toL2 = bsn.getFluxL3toL2()->matrix[r][c];
  double GWtoChn = bsn.getFluxGWtoChn()->matrix[r][c];
  double SnowtoSrf = bsn.getFluxSnowtoSrf()->matrix[r][c];
  // Lateral out
  double GWtoLat = bsn.getFluxGWtoLat()->matrix[r][c];
  double ChntoLat = Qk1*dtdx/dx;
  double SrftoLat = bsn.getFluxSrftoLat()->matrix[r][c];
  // Lateral in
  double LattoGW = bsn.getFluxLattoGW()->matrix[r][c]; 
  double LattoChn = bsn.getFluxLattoChn()->matrix[r][c]; 
  double LattoSrf = bsn.getFluxLattoSrf()->matrix[r][c]; 

  // For GW and surface (pond+channel), equivalent lateral inputs values
  double FinSrf = L1toSrf + GWtoChn + LattoChn;
  double FinSrf2 = L1toSrf + SnowtoSrf + LattoSrf;
  double d2Hin = 0;
  double d18Oin = 0;
  double Agein= 0;

  // for channel mixing
  double d2Hsurf;
  double d18Osurf;
  double Agesurf;
  double chan_store_old = bsn.getChanStoreOld()->matrix[r][c];
  double FinSrfChn = GWtoChn + LattoChn;

  // Two-pore stuff
  double theta_MW1 = 0;
  double theta_MW2 = 0;
  //double theta_r = 0;
  //double porosity = 0;
  double L3toTB2 = 0;
  double L3toMW2 = 0;
  double MW2toTB1 = 0;
  double MW2toMW1 = 0;

  // save the outflow tracer values of the fluxes
  double L2toL1_d2   = -1000;
  double L2toL1_o18  = -1000;
  double L2toL1_age  = -1000;
  double L1toSrf_d2  = -1000;
  double L1toSrf_o18 = -1000;
  double L1toSrf_age = -1000;
  double Srfout_d2   = -1000;
  double Srfout_o18  = -1000;
  double Srfout_age  = -1000;
  double Chnout_d2   = -1000;
  double Chnout_o18  = -1000;
  double Chnout_age  = -1000;

  if(ctrl.sw_TPD){
    //theta_r = bsn.getSoilMoistR()->matrix[r][c];
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

  // Layer 3 (GW included) --------------------------------------------------------------------

  if(LattoGW > RNDOFFERR){  
    if(ctrl.sw_2H){
      // Equivalent input signature for GW
      d2Hin = _Fd2HLattoGW->matrix[r][c] / LattoGW ;
      TracerMixing(bsn,ctrl,theta3_old*d3,_d2Hsoil3->matrix[r][c],_d2Hsoil3->matrix[r][c],
		   LattoGW,d2Hin,(L3toL2+GWtoLat),_d2Hgroundwater->matrix[r][c],S3,mixmod,r,c);
    } // end deuterium
    if(ctrl.sw_18O){
      d18Oin = _Fd18OLattoGW->matrix[r][c] / LattoGW ;
      TracerMixing(bsn,ctrl,theta3_old*d3,_d18Osoil3->matrix[r][c],_d18Osoil3->matrix[r][c],
		   LattoGW,d18Oin,(L3toL2+GWtoLat),_d18Ogroundwater->matrix[r][c],S3,mixmod,r,c);
    } // end oyxgen-18
    if(ctrl.sw_Age){
      Agein = _FAgeLattoGW->matrix[r][c] / LattoGW ;
      TracerMixing(bsn,ctrl,theta3_old*d3,_Agesoil3->matrix[r][c],_Agesoil3->matrix[r][c],
		   LattoGW,Agein,(L3toL2+GWtoLat),_Agegroundwater->matrix[r][c],S3,mixmod,r,c);
    } // end age
  } // end gw mixing

  // Layer 2 ------------------------------------------------------------------------

  // If two-pore domain activated
  if(ctrl.sw_TPD and L3toL2 > RNDOFFERR){
    // Tightly-bound
    if(ctrl.sw_2H)
      _d2H_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2,_d2H_TB2->matrix[r][c],
					L3toTB2, _d2Hgroundwater->matrix[r][c]);
    if(ctrl.sw_18O)
      _d18O_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2,_d18O_TB2->matrix[r][c],
					 L3toTB2, _d18Ogroundwater->matrix[r][c]);
    if(ctrl.sw_Age)
      _Age_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2,_Age_TB2->matrix[r][c],
					L3toTB2, _Agegroundwater->matrix[r][c]);
    
    // Mobile water
    if(ctrl.sw_2H){
      if(std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR){
	TracerMixing(bsn,ctrl,theta2_old*d2,_d2H_MW2->matrix[r][c],_d2H_MW2->matrix[r][c],
		     L3toMW2,_d2Hgroundwater->matrix[r][c],L2toL1,L2toL1_d2,S2,mixmod,r,c);
      } else {
	_d2H_MW2->matrix[r][c] = _d2Hgroundwater->matrix[r][c];
      }
    }// end deuterium

    if(ctrl.sw_18O){
      if(std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR){
	TracerMixing(bsn,ctrl,theta2_old*d2,_d18O_MW2->matrix[r][c],_d18O_MW2->matrix[r][c],
		     L3toMW2,_d18Ogroundwater->matrix[r][c],L2toL1,L2toL1_o18,S2,mixmod,r,c);
      } else {
	_d18O_MW2->matrix[r][c] = _d18Ogroundwater->matrix[r][c];
      }
    }// end oxygen-18

    if(ctrl.sw_Age){
      if(std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR){
	TracerMixing(bsn,ctrl,theta2_old*d2,_Age_MW2->matrix[r][c],_Age_MW2->matrix[r][c],
		     L3toMW2,_Agegroundwater->matrix[r][c],L2toL1,L2toL1_age,S2,mixmod,r,c);
      } else {
	_Age_MW2->matrix[r][c] = _Agegroundwater->matrix[r][c];
      }
    }//end age

  } else if (L3toL2 > RNDOFFERR) { // Soil-averaged values
    if(ctrl.sw_2H){
      TracerMixing(bsn,ctrl,theta2_old*d2,_d2Hsoil2->matrix[r][c],_d2Hsoil2->matrix[r][c],
		   L3toL2,_d2Hgroundwater->matrix[r][c],L2toL1,L2toL1_d2,S2,mixmod,r,c);
    } //end deuterium

    if(ctrl.sw_18O){
      TracerMixing(bsn,ctrl,theta2_old*d2,_d18Osoil2->matrix[r][c],_d18Osoil2->matrix[r][c],
		   L3toL2,_d18Ogroundwater->matrix[r][c],L2toL1,L2toL1_o18,S2,mixmod,r,c);
    } // end oxygen-18

    if(ctrl.sw_Age){
      TracerMixing(bsn,ctrl,theta2_old*d2,_Agesoil2->matrix[r][c],_Agesoil2->matrix[r][c],
		   L3toL2,_Agegroundwater->matrix[r][c],L2toL1,L2toL1_age,S2,mixmod,r,c);
    } // end water age
  }

  // Layer 1 ------------------------------------------------------------------------

  // If two-pore domain activated: return flow only from MW2
  if(ctrl.sw_TPD and L2toL1 > RNDOFFERR){
    // Tightly-bound
    if(ctrl.sw_2H){
      _d2H_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1,_d2H_TB1->matrix[r][c],	
					MW2toTB1, _d2H_MW2->matrix[r][c]);
    }

    if(ctrl.sw_18O){
      _d18O_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1,_d18O_TB1->matrix[r][c],
					 MW2toTB1, _d18O_MW2->matrix[r][c]);
    }

    if(ctrl.sw_Age)
      _Age_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1,_Age_TB1->matrix[r][c],
					MW2toTB1, _Age_MW2->matrix[r][c]);

    // Mobile water
    if(ctrl.sw_2H){
      if(std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR){
	L2toL1_d2 = L2toL1_d2 == -1000 ? _d2Hsoil2->matrix[r][c] : L2toL1_d2;
	TracerMixing(bsn,ctrl,theta1_old*d1,_d2H_MW1->matrix[r][c],_d2H_MW1->matrix[r][c],
		     MW2toMW1,L2toL1_d2,L1toSrf,L1toSrf_d2,S1,mixmod,r,c);
      } else {
	_d2H_MW1->matrix[r][c] = _d2H_MW1->matrix[r][c];
      }
    } // end deuterium

    if(ctrl.sw_18O){
      if(std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR){
	L2toL1_o18 = L2toL1_o18 == -1000 ? _d18Osoil2->matrix[r][c] : L2toL1_o18;
	TracerMixing(bsn,ctrl,theta1_old*d1,_d18O_MW1->matrix[r][c],_d18O_MW1->matrix[r][c],
		     MW2toMW1,L2toL1_o18,L1toSrf,L1toSrf_o18,S1,mixmod,r,c);
      } else {
	_d18O_MW1->matrix[r][c] = _d18O_MW1->matrix[r][c];
      }
    }

    if(ctrl.sw_Age){
      if(std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR){
	L2toL1_age = L2toL1_age == -1000 ? _Agesoil2->matrix[r][c] : L2toL1_age;
	TracerMixing(bsn,ctrl,theta1_old*d1,_Age_MW1->matrix[r][c],_Age_MW1->matrix[r][c],
		     MW2toMW1,L2toL1_age,L1toSrf,L1toSrf_age,S1,mixmod,r,c);
      } else {
	_Age_MW1->matrix[r][c] = _Age_MW1->matrix[r][c];
      }
    } // end age

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

  // Surface --------------------------------------------------------------------------------
  
  // If non-channel, return flow + runon + runoff.
  // If channel cell add input discharge, seepage, river outflow.
  // If surface storage initially null, propagated signature is that of inputs

  if(FinSrf > RNDOFFERR) {
    // If two-pore domain activated: return flow only from MW1
    if(ctrl.sw_TPD){
      if(ctrl.sw_2H){
	// If it's a channel cell, there's no reinfiltration, so that all lateral+snow+exfilt
	// inputs ultimately feed the channel
	if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){	  
	  // 1)      Mix surface pond with runoff, snowmelt, and exfiltration
	  d2Hsurf = L1toSrf_d2 == -1000 ? _d2H_MW1->matrix[r][c] : L1toSrf_d2;                                   // All return flow is from the mobile water
	  // 2)      Mix the runoff with exfiltrated and ponded water
	  if(pond_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,pond_old,_d2Hsurface->matrix[r][c],_d2Hsurface->matrix[r][c],
			 L1toSrf,d2Hsurf,(pond_old + L1toSrf),Srfout_d2,1.0,0,r,c);
	  } else {
	    _d2Hsurface->matrix[r][c] = d2Hsurf;
	  }
	  // 3)      Get d2H of the input to be mixed
	  d2Hin= ((pond_old + L1toSrf) * _d2Hsurface->matrix[r][c] + GWtoChn * _d2Hgroundwater->matrix[r][c] + 
		  _Fd2HLattoChn->matrix[r][c]) / (FinSrfChn + pond_old + L1toSrf);
	  // 4)      Mix new input with existing channel storage to get d2H for Qk1
	  if(chan_store_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,chan_store_old,_d2Hchan->matrix[r][c],_d2Hchan->matrix[r][c],
			 (FinSrfChn + pond_old + L1toSrf),d2Hin,(ChntoLat + SrftoLat),Chnout_d2,1.0,0,r,c);
	  } else {
	    _d2Hchan->matrix[r][c] = d2Hin;
	  }
	  _d2HGWtoChn->matrix[r][c] = _d2Hgroundwater->matrix[r][c];
	  _d2HSrftoChn->matrix[r][c] = FinSrf2 > RNDOFFERR ? (L1toSrf*d2Hsurf + _Fd2HLattoSrf->matrix[r][c] + 
							      SnowtoSrf*_d2Hsnowmelt->matrix[r][c]) / FinSrf2 : d2Hsurf;
	} else { 
	  // Equivalent input signature for surface inputs
	  d2Hsurf = L1toSrf_d2 == -1000 ? _d2H_MW1->matrix[r][c] : L1toSrf_d2;                                   // All return flow is from the mobile water
	  d2Hin = (L1toSrf*d2Hsurf + GWtoChn *_d2Hgroundwater->matrix[r][c] + _Fd2HLattoChn->matrix[r][c]) / FinSrf ;
	  if(pond_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,pond_old,_d2Hsurface->matrix[r][c],_d2Hsurface->matrix[r][c],
			 FinSrf,d2Hin,(ChntoLat+SrftoLat),Srfout_d2,1.0,0,r,c);
	    _d2Hchan->matrix[r][c] = _d2Hsurface->matrix[r][c];
	  } else {
	    _d2Hsurface->matrix[r][c] = d2Hin;
	    _d2Hchan->matrix[r][c] = _d2Hsurface->matrix[r][c];
	  }
	}// close channel mixing
      } // close d2H if statement

      if(ctrl.sw_18O){
	// If it's a channel cell, there's no reinfiltration, so that all lateral+snow+exfilt
	// inputs ultimately feed the channel
	if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){	  
	  // 1)      Mix surface pond with runoff, snowmelt, and exfiltration
	  d18Osurf = L1toSrf_o18 == -1000 ? _d18O_MW1->matrix[r][c] : L1toSrf_o18;                           // All return flow is from the mobile water
	  // 2)      Mix the runoff with exfiltrated and ponded water
	  if(pond_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,pond_old,_d18Osurface->matrix[r][c],_d18Osurface->matrix[r][c],
			 L1toSrf,d18Osurf,(pond_old + L1toSrf),Srfout_o18,1.0,0,r,c);
	  } else {
	    _d18Osurface->matrix[r][c] = d18Osurf;
	  }
	  // 3)      Get d18O of the input to be mixed
	  d18Oin= ((pond_old + L1toSrf) * _d18Osurface->matrix[r][c] + GWtoChn * _d18Ogroundwater->matrix[r][c] + 
		   _Fd18OLattoChn->matrix[r][c]) / (FinSrfChn + pond_old + L1toSrf);
	  // 4)      Mix new input with existing channel storage to get d18O for Qk1
	  if(chan_store_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,chan_store_old,_d18Ochan->matrix[r][c],_d18Ochan->matrix[r][c],
			 (FinSrfChn + pond_old + L1toSrf),d18Oin,(ChntoLat + SrftoLat),Chnout_o18,1.0,0,r,c);
	  } else {
	    _d18Ochan->matrix[r][c] = d18Oin;
	  }	  
	  _d18OGWtoChn->matrix[r][c] = _d18Ogroundwater->matrix[r][c];
	  _d18OSrftoChn->matrix[r][c] = FinSrf2 > RNDOFFERR ? (L1toSrf*d18Osurf + _Fd18OLattoSrf->matrix[r][c] + 
							       SnowtoSrf*_d18Osnowmelt->matrix[r][c]) / FinSrf2 : d18Osurf;
	} else {
	  // Equivalent input signature for surface inputs
	  d18Osurf = L1toSrf_o18 == -1000 ? _d18O_MW1->matrix[r][c] : L1toSrf_o18;                                   // All return flow is from the mobile water
	  d18Oin = (L1toSrf*d18Osurf + GWtoChn *_d18Ogroundwater->matrix[r][c] + _Fd18OLattoChn->matrix[r][c]) / FinSrf ;
	  if(pond_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,pond_old,_d18Osurface->matrix[r][c],_d18Osurface->matrix[r][c],
			 FinSrf,d18Oin,(ChntoLat+SrftoLat),Srfout_o18,1.0,0,r,c);
	    _d18Ochan->matrix[r][c] = _d18Osurface->matrix[r][c];
	  } else {
	    _d18Osurface->matrix[r][c] = d18Oin;
	    _d18Ochan->matrix[r][c] = _d18Osurface->matrix[r][c];
	  }
	} // close channel mixing
      } // close d18o if statement

      if(ctrl.sw_Age){
	// If it's a channel cell, there's no reinfiltration, so that all lateral+snow+exfilt
	// inputs ultimately feed the channel
	if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){	  
	  // 1)      Mix surface pond with runoff, snowmelt, and exfiltration
	  Agesurf = L1toSrf_age == -1000 ? _Age_MW1->matrix[r][c] : L1toSrf_age;                           // All return flow is from the mobile water
	  // 2)      Mix the runoff with exfiltrated and ponded water
	  if(pond_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,pond_old,_Agesurface->matrix[r][c],_Agesurface->matrix[r][c],
			 L1toSrf,Agesurf,(pond_old + L1toSrf),Srfout_age,1.0,0,r,c);
	  } else {
	    _Agesurface->matrix[r][c] = Agesurf;
	  }
	  // 3)      Get Age of the input to be mixed
	  Agein= ((pond_old + L1toSrf) * _Agesurface->matrix[r][c] + GWtoChn * _Agegroundwater->matrix[r][c] + 
		  _FAgeLattoChn->matrix[r][c]) / (FinSrfChn + pond_old + L1toSrf);
	  // 4)      Mix new input with existing channel storage to get Age for Qk1
	  if(chan_store_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,chan_store_old,_Agechan->matrix[r][c],_Agechan->matrix[r][c],
			 (FinSrfChn + pond_old + L1toSrf),Agein,(ChntoLat + SrftoLat),Chnout_age,1.0,0,r,c);
	  } else {
	    _Agechan->matrix[r][c] = d18Oin;
	  }	  
	  _AgeGWtoChn->matrix[r][c] = _Agegroundwater->matrix[r][c];
	  _AgeSrftoChn->matrix[r][c] = FinSrf2 > RNDOFFERR ? (L1toSrf*Agesurf + _FAgeLattoSrf->matrix[r][c] + 
							      SnowtoSrf*_Agesnowmelt->matrix[r][c]) / FinSrf2 : 0.0;
	} else { 
	  // Equivalent input signature for surface inputs
	  Agesurf = L1toSrf_age == -1000 ? _Age_MW1->matrix[r][c] : L1toSrf_age;                           // All return flow is from the mobile water
	  Agein = (L1toSrf*Agesurf + GWtoChn *_Agegroundwater->matrix[r][c] + _FAgeLattoChn->matrix[r][c]) / FinSrf ;
	  if(pond_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,pond_old,_Agesurface->matrix[r][c],_Agesurface->matrix[r][c],
			 FinSrf,Agein,(ChntoLat+SrftoLat),Srfout_age,1.0,0,r,c);
	    _Agechan->matrix[r][c] = _Agesurface->matrix[r][c];
	  } else {
	    _Agesurface->matrix[r][c] = Agein;
	    _Agechan->matrix[r][c] = _Agesurface->matrix[r][c];
	  }
	} //close channel mixing
      } // close age if statement
      
    } else {      // soil averaged water
      if(ctrl.sw_2H){
	// If it's a channel cell, there's no reinfiltration, so that all lateral+snow+exfilt
	// inputs ultimately feed the channel
	if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){	  
	  // 1)      Mix surface pond with runoff, snowmelt, and exfiltration
	  d2Hsurf = L1toSrf_d2 == -1000 ? _d2Hsoil1->matrix[r][c] : L1toSrf_d2;                                   // All return flow is from the mobile water
	  // 2)      Mix the runoff with exfiltrated and ponded water
	  if(pond_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,pond_old,_d2Hsurface->matrix[r][c],_d2Hsurface->matrix[r][c],
			 L1toSrf,d2Hsurf,(pond_old + L1toSrf),Srfout_d2,1.0,0,r,c);
	  } else {
	    _d2Hsurface->matrix[r][c] = d2Hsurf;
	  }
	  // 3)      Get d2H of the input to be mixed
	  d2Hin= ((pond_old + L1toSrf) * _d2Hsurface->matrix[r][c] + GWtoChn * _d2Hgroundwater->matrix[r][c] + 
		  _Fd2HLattoChn->matrix[r][c]) / (FinSrfChn + pond_old + L1toSrf);
	  // 4)      Mix new input with existing channel storage to get d2H for Qk1
	  if(chan_store_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,chan_store_old,_d2Hchan->matrix[r][c],_d2Hchan->matrix[r][c],
			 (FinSrfChn + pond_old + L1toSrf),d2Hin,(ChntoLat + SrftoLat),Chnout_d2,1.0,0,r,c);
	  } else {
	    _d2Hchan->matrix[r][c] = d2Hin;
	  }
	  _d2HGWtoChn->matrix[r][c] = _d2Hgroundwater->matrix[r][c];
	  _d2HSrftoChn->matrix[r][c] = FinSrf2 > RNDOFFERR ? (L1toSrf*d2Hsurf + _Fd2HLattoSrf->matrix[r][c] +  
							      SnowtoSrf*_d2Hsnowmelt->matrix[r][c]) / FinSrf2 : d2Hsurf;
	} else {// close channel mixing
	  // Equivalent input signature for surface inputs
	  d2Hsurf = L1toSrf_d2 == -1000 ? _d2Hsoil1->matrix[r][c] : L1toSrf_d2;                                   // All return flow is from the mobile water
	  d2Hin = (L1toSrf*d2Hsurf + GWtoChn *_d2Hgroundwater->matrix[r][c] + _Fd2HLattoChn->matrix[r][c]) / FinSrf ;
	  if(pond_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,pond_old,_d2Hsurface->matrix[r][c],_d2Hsurface->matrix[r][c],
			 FinSrf,d2Hin,(ChntoLat+SrftoLat),Srfout_d2,1.0,0,r,c);
	    _d2Hchan->matrix[r][c] = _d2Hsurface->matrix[r][c];
	  } else {
	    _d2Hsurface->matrix[r][c] = d2Hin;
	    _d2Hchan->matrix[r][c] = _d2Hsurface->matrix[r][c];
	  }
	} //close channel mixing
      } // close d2H if statement
      
      if(ctrl.sw_18O){
	// If it's a channel cell, there's no reinfiltration, so that all lateral+snow+exfilt
	// inputs ultimately feed the channel
	if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){	  
	  // 1)      Mix surface pond with runoff, snowmelt, and exfiltration
	  d18Osurf = L1toSrf_o18 == -1000 ? _d18Osoil1->matrix[r][c] : L1toSrf_o18;                                 // All return flow is from the mobile water
	  // 2)      Mix the runoff with exfiltrated and ponded water
	  if(pond_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,pond_old,_d18Osurface->matrix[r][c],_d18Osurface->matrix[r][c],
			 L1toSrf,d18Osurf,(pond_old + L1toSrf),Srfout_o18,1.0,0,r,c);
	  } else {
	    _d18Osurface->matrix[r][c] = d18Osurf;
	  }
	  // 3)      Get d2H of the input to be mixed
	  d18Oin= ((pond_old + L1toSrf) * _d18Osurface->matrix[r][c] + GWtoChn * _d18Ogroundwater->matrix[r][c] + 
		  _Fd18OLattoChn->matrix[r][c]) / (FinSrfChn + pond_old + L1toSrf);
	  // 4)      Mix new input with existing channel storage to get d2H for Qk1
	  if(chan_store_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,chan_store_old,_d18Ochan->matrix[r][c],_d18Ochan->matrix[r][c],
			 (FinSrfChn + pond_old + L1toSrf),d18Oin,(ChntoLat + SrftoLat),Chnout_o18,1.0,0,r,c);
	  } else {
	    _d2Hchan->matrix[r][c] = d2Hin;
	  }
	  _d18OGWtoChn->matrix[r][c] = _d18Ogroundwater->matrix[r][c];
	  _d18OSrftoChn->matrix[r][c] = FinSrf2 > RNDOFFERR ? (L1toSrf*d18Osurf + _Fd18OLattoSrf->matrix[r][c] +  
							      SnowtoSrf*_d18Osnowmelt->matrix[r][c]) / FinSrf2 : d18Osurf;
	} else {// close channel mixing
	  // Equivalent input signature for surface inputs
	  d18Osurf = L1toSrf_o18 == -1000 ? _d18Osoil1->matrix[r][c] : L1toSrf_o18;                                   // All return flow is from the mobile water
	  d18Oin = (L1toSrf*d18Osurf + GWtoChn *_d18Ogroundwater->matrix[r][c] + _Fd18OLattoChn->matrix[r][c]) / FinSrf ;
	  if(pond_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,pond_old,_d18Osurface->matrix[r][c],_d18Osurface->matrix[r][c],
			 FinSrf,d18Oin,(ChntoLat+SrftoLat),Srfout_o18,1.0,0,r,c);
	    _d18Ochan->matrix[r][c] = _d18Osurface->matrix[r][c];
	  } else {
	    _d18Osurface->matrix[r][c] = d18Oin;
	    _d18Ochan->matrix[r][c] = _d18Osurface->matrix[r][c];
	  }
	} //close channel mixing
      } // close d18o if statement
      
      if(ctrl.sw_Age) {

	// If it's a channel cell, there's no reinfiltration, so that all lateral+snow+exfilt
	// inputs ultimately feed the channel
	if(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0){	  
	  // 1)      Mix surface pond with runoff, snowmelt, and exfiltration
	  Agesurf = L1toSrf_age == -1000 ? _Agesoil1->matrix[r][c] : L1toSrf_age;                                 // All return flow is from the mobile water 
	  // 2)      Mix the runoff with exfiltrated and ponded water
	  if(pond_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,pond_old,_Agesurface->matrix[r][c],_Agesurface->matrix[r][c],
			 L1toSrf,Agesurf,(pond_old + L1toSrf),Srfout_age,1.0,0,r,c);
	  } else {
	    _Agesurface->matrix[r][c] = Agesurf;
	  }
	  // 3)      Get Age of the input to be mixed
	  Agein= ((pond_old + L1toSrf) * _Agesurface->matrix[r][c] + GWtoChn * _Agegroundwater->matrix[r][c] +
		  _FAgeLattoChn->matrix[r][c]) / (FinSrfChn + pond_old + L1toSrf);
	  // 4)      Mix new input with existing channel storage to get Age for Qk1
	  if(chan_store_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,chan_store_old,_Agechan->matrix[r][c],_Agechan->matrix[r][c],
			 (FinSrfChn + pond_old + L1toSrf),Agein,(ChntoLat + SrftoLat),Chnout_age,1.0,0,r,c);
	  } else {
	    _Agechan->matrix[r][c] = Agein;
	  }
	  _AgeGWtoChn->matrix[r][c] = _Agegroundwater->matrix[r][c];
	  _AgeSrftoChn->matrix[r][c] = FinSrf2 > RNDOFFERR ? (L1toSrf*Agesurf + _FAgeLattoSrf->matrix[r][c] + 
							      SnowtoSrf*_Agesnowmelt->matrix[r][c]) / FinSrf2 : Agesurf;
	} else {
	  Agesurf = L1toSrf_age == -1000 ? _Agesoil1->matrix[r][c] : L1toSrf_age;
	  Agein = (L1toSrf*Agesurf + GWtoChn *_Agegroundwater->matrix[r][c] + _FAgeLattoChn->matrix[r][c]) / FinSrf ;
	  if(pond_old > RNDOFFERR){
	    TracerMixing(bsn,ctrl,pond_old,_Agesurface->matrix[r][c],_Agesurface->matrix[r][c],
		       FinSrf,Agein,(ChntoLat+SrftoLat),Srfout_age,1.0,0,r,c);
	    _Agechan->matrix[r][c] = _Agesurface->matrix[r][c];
	  } else {
	    _Agesurface->matrix[r][c] = Agein;
	    _Agechan->matrix[r][c] = _Agesurface->matrix[r][c];
	  }
	} //close channel mixing
      } // close age if statement
    }
  }
  // -------------------------------------------------------------------------------
}

