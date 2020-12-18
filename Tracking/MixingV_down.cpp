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
 * MixingV_down.cpp
 *
 *  Created on: Jun 21, 2017
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"

void Tracking::MixingV_down(Basin &bsn, Control &ctrl,
			    double &d1, double &d2, double &d3, double &fc,
			    int r, int c, 
			    int step) // 0 for SolveSurface, 1 for GWrouting
{

  double theta1_old = bsn.getSoilMoist1()->matrix[r][c];
  double theta2_old = bsn.getSoilMoist2()->matrix[r][c];
  double theta3_old = bsn.getSoilMoist3()->matrix[r][c];

  double S1 = std::max<double>((theta1_old - bsn.getSoilMoistR1()->matrix[r][c])/
			       (bsn.getPorosityL1()->matrix[r][c] - bsn.getSoilMoistR1()->matrix[r][c]),0);
  double S2 = std::max<double>((theta2_old - bsn.getSoilMoistR2()->matrix[r][c])/
			       (bsn.getPorosityL2()->matrix[r][c] - bsn.getSoilMoistR2()->matrix[r][c]),0);
  double S3 = std::max<double>((theta3_old - bsn.getSoilMoistR3()->matrix[r][c])/
			       (bsn.getPorosityL3()->matrix[r][c] - bsn.getSoilMoistR3()->matrix[r][c]),0);

  double LattoSrf = 0;
  double SnowtoSrf = 0;
  double FinSrf = 0;
  double pond_old = 0;
  double d2Hin = 0 ;
  double d18Oin = 0;
  double Agein = 0; 

  double SrftoL1 = bsn.getFluxSrftoL1()->matrix[r][c];
  double L1toL2 = bsn.getFluxL1toL2()->matrix[r][c];
  double L2toL3 = bsn.getFluxL2toL3()->matrix[r][c];
  double leak = step == 0 ? bsn.getFluxLeak()->matrix[r][c] : 0;

  int mixmod = ctrl.toggle_mix;

  double theta_MW1 = 0;
  double theta_MW2 = 0;
  //double theta_r = 0;
  //double porosity = 0;
  double SrftoTB1 = 0;
  double SrftoMW1 = 0;
  double MW1toTB2 = 0;
  double MW1toMW2 = 0;

  // save the outflow tracer values of the fluxes
  double SrftoL1_d2  = -1000;
  double SrftoL1_o18 = -1000;
  double SrftoL1_age = -1000;
  double L1toL2_d2   = -1000;
  double L1toL2_o18  = -1000;
  double L1toL2_age  = -1000;
  double L2toL3_d2   = -1000;
  double L2toL3_o18  = -1000;
  double L2toL3_age  = -1000;
  double L3out_d2    = -1000;
  double L3out_o18   = -1000;
  double L3out_age   = -1000;
 
  if(ctrl.sw_2H and step==0){
    _d2HGWtoChn->reset();
    _d2HSrftoChn->reset();
    _d2HRecharge->reset();
  }
  if(ctrl.sw_18O and step==0){
    _d18OGWtoChn->reset();
    _d18OSrftoChn->reset();
    _d18ORecharge->reset();
  }
  if(ctrl.sw_Age and step==0){
    _AgeGWtoChn->reset();
    _AgeSrftoChn->reset();
    _AgeRecharge->reset();
  }
    
  if(ctrl.sw_TPD){
    //theta_r = bsn.getSoilMoistR()->matrix[r][c];
    theta_MW1 = bsn.getMoistureMW1()->matrix[r][c];
    theta_MW2 = bsn.getMoistureMW2()->matrix[r][c];
    //porosity = bsn.getPorosity()->matrix[r][c];
    // Routing to tightly-bound domain: fraction from relative volume of "inactive"-pore domain, 
    // limited by the "available space" in tightly-bound domain max(0,d1*(theta_MW-theta1_old))
    SrftoTB1 = std::min<double>(SrftoL1,//*(theta_MW-theta_r)/(porosity-theta_r),
				std::max<double>(0,d1*(theta_MW1-theta1_old)));
    // Remainder directly goes to mobile water
    SrftoMW1 = std::max<double>(0,SrftoL1-SrftoTB1);
    // L2: Percolation only from mobile water in L1
    MW1toTB2 = std::min<double>(L1toL2,//*(theta_MW-theta_r)/(porosity-theta_r),
				std::max<double>(0,d2*(theta_MW2-theta2_old)));
    MW1toMW2 = std::max<double>(0,L1toL2-MW1toTB2);
    // L3: Percolation only from mobile water water in L2 (MW2toL3 = L2toL3). 
  }

  // Surface (if reinfiltrating) ----------------------------------------------------
  if(step == 1){

    // Fluxes in
    SnowtoSrf = bsn.getFluxSnowtoSrf()->matrix[r][c];
    LattoSrf = bsn.getFluxLattoSrf()->matrix[r][c]; 
    FinSrf = SnowtoSrf + LattoSrf;
    
    if(FinSrf > RNDOFFERR) {
      
      pond_old = bsn.getPondingWater()->matrix[r][c] - FinSrf;

      if(ctrl.sw_2H){    
	d2Hin = (_Fd2HLattoSrf->matrix[r][c] + SnowtoSrf*_d2Hsnowmelt->matrix[r][c])/ FinSrf ; // input to ponded water mixing
	if(pond_old > RNDOFFERR){
	  TracerMixing(bsn,ctrl,pond_old,_d2Hsurface->matrix[r][c],_d2Hsurface->matrix[r][c],
		       FinSrf,d2Hin,SrftoL1,SrftoL1_d2,1.0,ctrl.toggle_mix,r,c);
	} else {
	  _d2Hsurface->matrix[r][c] = d2Hin;
	  SrftoL1_d2 = d2Hin;
	}
      } // end deuterium

      if(ctrl.sw_18O){
	// Update surface
	d18Oin = (_Fd18OLattoSrf->matrix[r][c] + SnowtoSrf*_d18Osnowmelt->matrix[r][c])/ FinSrf ;
	if(pond_old > RNDOFFERR){
	  TracerMixing(bsn,ctrl,pond_old,_d18Osurface->matrix[r][c],_d18Osurface->matrix[r][c],
		       FinSrf,d18Oin,SrftoL1,SrftoL1_o18,1.0,ctrl.toggle_mix,r,c);
	} else {
	  _d18Osurface->matrix[r][c] = d18Oin;
	  SrftoL1_o18 = d18Oin;
	}
      } // end oxygen-18

      if(ctrl.sw_Age){
	Agein = (_FAgeLattoSrf->matrix[r][c] + SnowtoSrf*_Agesnowmelt->matrix[r][c])/ FinSrf ;
	if(pond_old > RNDOFFERR){
	  TracerMixing(bsn,ctrl,pond_old,_Agesurface->matrix[r][c],_Agesurface->matrix[r][c],
		       FinSrf,Agein,SrftoL1,SrftoL1_age,1.0,ctrl.toggle_mix,r,c);
	} else {
	  _Agesurface->matrix[r][c] = Agein;
	  SrftoL1_age = Agein;
	}
      } // end water age
    } // end the surface water      
  } // end if step == 1




  // === Only if initial infitlration or reinfilt + non-channel
  if(step == 0 or (step==1 and ctrl.sw_reinfilt and 
		   !(ctrl.sw_channel and bsn.getChannelWidth()->matrix[r][c] > 0))){

    // Layer 1 ------------------------------------------------------------------------
  
    // If two-pore domain activated, and different groundwater origin (MW2 instead of soil2)
    if(ctrl.sw_TPD and SrftoL1>RNDOFFERR){
      //TODO update the saturation for each TB and MW
      if(ctrl.sw_2H){
	SrftoL1_d2 = SrftoL1_d2 == -1000 ? _d2Hsurface->matrix[r][c] : SrftoL1_d2;

	_d2H_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1, _d2H_TB1->matrix[r][c],
					  SrftoTB1, SrftoL1_d2);
     
	if(std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR){
	  TracerMixing(bsn,ctrl,std::max<double>(0,theta1_old-theta_MW1)*d1,_d2H_MW1->matrix[r][c],
		       _d2H_MW1->matrix[r][c],SrftoMW1,SrftoL1_d2,L1toL2,L1toL2_d2,S1,mixmod,r,c);
	} else {
	  _d2H_MW1->matrix[r][c] = SrftoL1_d2;
	}
      }// end deuterium
    
      if(ctrl.sw_18O){
	SrftoL1_o18 = SrftoL1_o18 == -1000 ? _d18Osurface->matrix[r][c] : SrftoL1_o18;

	_d18O_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1, _d18O_TB1->matrix[r][c],
					   SrftoTB1, SrftoL1_o18);
      
	if(std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR){
	  TracerMixing(bsn,ctrl,std::max<double>(0,theta1_old-theta_MW1)*d1,_d18O_MW1->matrix[r][c],
		       _d18O_MW1->matrix[r][c],SrftoMW1,SrftoL1_o18,L1toL2,L1toL2_o18,S1,mixmod,r,c);
	} else {
	  _d18O_MW1->matrix[r][c] = SrftoL1_o18;
	}
      }// end oxygen-18
    
      if(ctrl.sw_Age){
	SrftoL1_age = SrftoL1_age == -1000 ? _Agesurface->matrix[r][c] : SrftoL1_age;

	_Age_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1, _Age_TB1->matrix[r][c], 
					  SrftoTB1, SrftoL1_age);

	if(std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR){
	  TracerMixing(bsn,ctrl,std::max<double>(0,theta1_old-theta_MW1)*d1,_Age_MW1->matrix[r][c],
		       _Age_MW1->matrix[r][c],SrftoMW1,SrftoL1_age,L1toL2,L1toL2_age,S1,mixmod,r,c);
	} else {
	  _Age_MW1->matrix[r][c] = SrftoL1_age;
	}
      }// end water age

    } else if (SrftoL1>RNDOFFERR){ // Soil-averaged
      if(ctrl.sw_2H){
	SrftoL1_d2 = SrftoL1_d2 == -1000 ? _d2Hsurface->matrix[r][c] : SrftoL1_d2;
	TracerMixing(bsn,ctrl,theta1_old*d1,_d2Hsoil1->matrix[r][c],_d2Hsoil1->matrix[r][c],
		     SrftoL1,SrftoL1_d2,L1toL2,L1toL2_d2,S1,ctrl.toggle_mix,r,c);
      } // end deuterium

      if(ctrl.sw_18O){
	SrftoL1_o18 = SrftoL1_o18 == -1000 ? _d18Osurface->matrix[r][c] : SrftoL1_o18;
	TracerMixing(bsn,ctrl,theta1_old*d1,_d18Osoil1->matrix[r][c],_d18Osoil1->matrix[r][c],
		     SrftoL1,SrftoL1_o18,L1toL2,L1toL2_o18,S1,ctrl.toggle_mix,r,c);
      } // end oxygen-18

      if(ctrl.sw_Age){
	SrftoL1_age = SrftoL1_age == -1000 ? _Agesurface->matrix[r][c] : SrftoL1_age;
	TracerMixing(bsn,ctrl,theta1_old*d1,_Agesoil1->matrix[r][c],_Agesoil1->matrix[r][c],
		     SrftoL1,SrftoL1_age,L1toL2,L1toL2_age,S1,ctrl.toggle_mix,r,c);
      } // end water age
    }
  
    // Layer 2 ------------------------------------------------------------------------
  
    // If two-pore domain activated: only MW1 percolates
    if(ctrl.sw_TPD and L1toL2>RNDOFFERR){
      //TODO update the saturation for each TB and MW
      if(ctrl.sw_2H){
	_d2H_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2,_d2H_TB2->matrix[r][c],
					  MW1toTB2, L1toL2_d2);

	if(std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR){
	  TracerMixing(bsn,ctrl,std::max<double>(0,theta2_old-theta_MW2)*d2,_d2H_MW2->matrix[r][c],
		       _d2H_MW2->matrix[r][c],MW1toMW2,L1toL2_d2,L2toL3,L2toL3_d2,S2,mixmod,r,c);
	} else {
	  _d2H_MW2->matrix[r][c] = _d2H_MW1->matrix[r][c];
	}
      }// end deuterium

      if(ctrl.sw_18O){
	_d18O_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2,_d18O_TB2->matrix[r][c], 
					   MW1toTB2, L1toL2_o18);

	if(std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR){
	  TracerMixing(bsn,ctrl,std::max<double>(0,theta2_old-theta_MW2)*d2,_d18O_MW2->matrix[r][c],
		       _d18O_MW2->matrix[r][c],MW1toMW2,L1toL2_o18,L2toL3,L2toL3_o18,S2,mixmod,r,c);
	} else {
	  _d18O_MW2->matrix[r][c] = _d18O_MW1->matrix[r][c];
	}
      
      } // end oxygen-18
      if(ctrl.sw_Age){
	_Age_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2,_Age_TB2->matrix[r][c], 
					  MW1toTB2, L1toL2_age);

	if(std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR){
	  TracerMixing(bsn,ctrl,std::max<double>(0,theta2_old-theta_MW2)*d2,_Age_MW2->matrix[r][c],
		       _Age_MW2->matrix[r][c],MW1toMW2,L1toL2_age,L2toL3,L2toL3_age,S2,mixmod,r,c);
	} else {
	  _Age_MW2->matrix[r][c] = _Age_MW1->matrix[r][c];
	}
      }//end water age

    } else if (L1toL2 > RNDOFFERR){ // Soil-averaged
      if(ctrl.sw_2H){
	L1toL2_d2 = L1toL2_d2 == -1000 ? _d2Hsoil1->matrix[r][c] : L1toL2_d2;
    	TracerMixing(bsn,ctrl,theta2_old*d2,_d2Hsoil2->matrix[r][c],_d2Hsoil2->matrix[r][c],
		     L1toL2,L1toL2_d2,L2toL3,L2toL3_d2,S2,ctrl.toggle_mix,r,c);
      } // end deuterium
      if(ctrl.sw_18O){
	L1toL2_o18 = L1toL2_o18 == -1000 ? _d18Osoil1->matrix[r][c] : L1toL2_o18;
	TracerMixing(bsn,ctrl,theta2_old*d2,_d18Osoil2->matrix[r][c],_d18Osoil2->matrix[r][c],
		     L1toL2,L1toL2_o18,L2toL3,L2toL3_o18,S2,ctrl.toggle_mix,r,c);
      } // end oxygen-18
      if(ctrl.sw_Age){
	L1toL2_age = L1toL2_age == -1000 ? _Agesoil1->matrix[r][c] : L1toL2_age;
    	TracerMixing(bsn,ctrl,theta2_old*d2,_Agesoil2->matrix[r][c],_Agesoil2->matrix[r][c],
		     L1toL2,L1toL2_age,L2toL3,L2toL3_age,S2,ctrl.toggle_mix,r,c);
      }
    }
  
    // Layer 3 ------------------------------------------------------------------------
    
    // If two-pore domain activated: only MW2 percolates
    if(ctrl.sw_TPD and L2toL3>RNDOFFERR){
      if(ctrl.sw_2H){
    	TracerMixing(bsn,ctrl,theta3_old*d3,_d2Hsoil3->matrix[r][c],_d2Hsoil3->matrix[r][c],
		     L2toL3,L2toL3_d2,leak,L3out_d2,S3,ctrl.toggle_mix,r,c);

	_d2HRecharge->matrix[r][c] = step == 0 ? L2toL3_d2 :
	  InputMix(std::max<double>(0.0,bsn.getFluxRecharge()->matrix[r][c] - L2toL3), 
		   _d2HRecharge->matrix[r][c], L2toL3, L2toL3_d2);
      }

      if(ctrl.sw_18O){
    	TracerMixing(bsn,ctrl,theta3_old*d3,_d18Osoil3->matrix[r][c],_d18Osoil3->matrix[r][c],
		     L2toL3,L2toL3_o18,leak,L3out_o18,S3,ctrl.toggle_mix,r,c);

	_d18ORecharge->matrix[r][c] = step == 0 ? L2toL3_o18 :
	  InputMix(std::max<double>(0.0,bsn.getFluxRecharge()->matrix[r][c] - L2toL3), 
		   _d18ORecharge->matrix[r][c], L2toL3, L2toL3_o18);
      }

      if(ctrl.sw_Age){
    	TracerMixing(bsn,ctrl,theta3_old*d3,_Agesoil3->matrix[r][c],_Agesoil3->matrix[r][c],
		     L2toL3,L2toL3_age,leak,L3out_age,S3,ctrl.toggle_mix,r,c);

	_AgeRecharge->matrix[r][c] = step == 0 ? L2toL3_age :
	  InputMix(std::max<double>(0.0,bsn.getFluxRecharge()->matrix[r][c] - L2toL3), 
		   _AgeRecharge->matrix[r][c], L2toL3, L2toL3_age);
      }
      
    } else if (L2toL3>RNDOFFERR){// Soil-averaged
      
      if(ctrl.sw_2H){
	L2toL3_d2 = L2toL3_d2 == -1000 ? _d2Hsoil2->matrix[r][c] : L2toL3_d2;
    	TracerMixing(bsn,ctrl,theta3_old*d3,_d2Hsoil3->matrix[r][c],_d2Hsoil3->matrix[r][c],
		     L2toL3,L2toL3_d2,leak,L3out_d2,S3,ctrl.toggle_mix,r,c);

	_d2HRecharge->matrix[r][c] = step == 0 ? L2toL3_d2 :
	       InputMix(std::max<double>(0.0,bsn.getFluxRecharge()->matrix[r][c] - L2toL3), 
		   _d2HRecharge->matrix[r][c], L2toL3, L2toL3_d2);

      }
      if(ctrl.sw_18O){
	L2toL3_o18 = L2toL3_o18 == -1000 ? _d18Osoil2->matrix[r][c] : L2toL3_o18;
    	TracerMixing(bsn,ctrl,theta3_old*d3,_d18Osoil3->matrix[r][c],_d18Osoil3->matrix[r][c],
		     L2toL3,L2toL3_o18,leak,L3out_o18,S3,ctrl.toggle_mix,r,c);

	_d18ORecharge->matrix[r][c] = step == 0 ? L2toL3_o18 :
	       InputMix(std::max<double>(0.0,bsn.getFluxRecharge()->matrix[r][c] - L2toL3), 
		   _d18ORecharge->matrix[r][c], L2toL3, L2toL3_o18);
      }
      if(ctrl.sw_Age){
	L2toL3_age = L2toL3_age == -1000 ? _Agesoil2->matrix[r][c] : L2toL3_age;
    	TracerMixing(bsn,ctrl,theta3_old*d3,_Agesoil3->matrix[r][c],_Agesoil3->matrix[r][c],
		     L2toL3,L2toL3_age,leak,L3out_age,S3,ctrl.toggle_mix,r,c);

	_AgeRecharge->matrix[r][c] = step == 0 ? L2toL3_age :
	       InputMix(std::max<double>(0.0,bsn.getFluxRecharge()->matrix[r][c] - L2toL3), 
		   _AgeRecharge->matrix[r][c], L2toL3, L2toL3_age);
      }
    }

  // Groundwater ------------------------------------------------------------------------
  
    if(ctrl.sw_2H)
      _d2Hgroundwater->matrix[r][c] = L3out_d2 == -1000 ? _d2Hsoil3->matrix[r][c] : L3out_d2;
	//_d2Hgroundwater->matrix[r][c] = _d2Hsoil3->matrix[r][c];
      
    if(ctrl.sw_18O)
      _d18Ogroundwater->matrix[r][c] = L3out_o18 == -1000 ? _d18Osoil3->matrix[r][c] : L3out_o18;
        //_d18Ogroundwater->matrix[r][c] = _d18Osoil3->matrix[r][c];

    if(ctrl.sw_Age)
      _Agegroundwater->matrix[r][c] = L3out_age == -1000 ? _Agesoil3->matrix[r][c] : L3out_age;
        //_Agegroundwater->matrix[r][c] = _Agesoil3->matrix[r][c];

  }

  // -- Leakage : only if first infiltration (in SolveSurfaceFluxes)
  if(step == 0){
    //    cout << "L3out_d2: " << L3out_d2 << " |d2hGW: " << _d2Hgroundwater->matrix[r][c] << endl;
    // Backup pre-routing GW d2H for use in MixingRouting2.cpp
    if(ctrl.sw_2H){
      //_d2Hleakage->matrix[r][c] = _d2Hgroundwater->matrix[r][c];
      _d2Hleakage->matrix[r][c] = L3out_d2 == -1000 ? _d2Hgroundwater->matrix[r][c] : L3out_d2;
    }
    if(ctrl.sw_18O){
      //_d18Oleakage->matrix[r][c] = _d18Ogroundwater->matrix[r][c]; 
      _d18Oleakage->matrix[r][c] = L3out_o18 == -1000 ? _d18Ogroundwater->matrix[r][c] : L3out_o18;
    }
    if(ctrl.sw_Age){
      //_Ageleakage->matrix[r][c] = _Agegroundwater->matrix[r][c]; 
      _Ageleakage->matrix[r][c] = L3out_age == -1000 ? _Agegroundwater->matrix[r][c] : L3out_age;
    }
  }  
}


