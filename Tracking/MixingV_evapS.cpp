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
 * MixingV_evapS.cpp
 *
 *  Created on: Jun 21, 2017
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"

void Tracking::MixingV_evapS(Atmosphere &atm, Basin &bsn, Control &ctrl, 
			     double &d1, double &theta_new,
			     double Ts, double &etp, double &beta,
			     double &d2H_evap, double &d18O_evap, double &Age_evap,
			     int r, int c)
{
  // p-weighted soil evap
  double evapS   = etp * ctrl.dt; // etp > RNDOFFERR ? etp * ctrl.dt : 0.0 ;
  //vegetation-summed Es (until this point)
  double etotS   = bsn.getEvaporationS_all()->matrix[r][c] * ctrl.dt; 
  //double SrftoL1 = bsn.getFluxSrftoL1()->matrix[r][c];

  // useful effective volume
  double V_new = theta_new*d1;
  double V_old = V_new + evapS;
  double V_TBnew = 0;
  double V_MWnew = 0;
  double V_TBold = 0;
  double V_MWold = 0;
  double evapS_TB = 0;
  double evapS_MW = 0;

  double d2H_new, d18O_new;//, Age_new;
  double d2H_TBevap, d18O_TBevap;//, Age_TBevap;
  double d2H_MWevap, d18O_MWevap;//, Age_MWevap;
  
  double theta_MW1 = 0;
  double poros = 0;

  double d_old = 0;

  if(ctrl.sw_TPD){
    theta_MW1 = bsn.getMoistureMW1()->matrix[r][c];
    poros = bsn.getPorosityL1()->matrix[r][c];
    V_TBold = min<double>(V_old, theta_MW1*d1);
    V_MWold = max<double>(0, V_old - V_TBold);
    //evapS_MW = min<double>(evapS*(poros - theta_MW1)/(poros - bsn.getSoilMoistR()->matrix[r][c]),
    //			   V_MWold);
    // Evaporation takes from both domains proportionnally to their content (and not of 
    // relavitve pore volume, as before)
    evapS_MW = evapS * V_MWold / V_old;
    evapS_TB = evapS - evapS_MW;
    V_MWnew = V_MWold - evapS_MW;
    V_TBnew = V_TBold - evapS_TB;
  }
  
  // Two-pore domain ----------------------------------------------------------------------------
  
  if(ctrl.sw_TPD){
    // Safeguard with minimum E value: avoids diverging d2H_evap values
    if(ctrl.sw_frac and V_new/V_old <0.999) {
	if(ctrl.sw_2H){
	  // Tightly-bound
	  Frac_Esoil(atm, bsn, ctrl, V_TBold, V_TBnew, beta, 
		     _d2H_TB1->matrix[r][c], d2H_new, d2H_TBevap, Ts, r, c, 0);
	  if(abs(d2H_new)>100)
	    cout << r << " " << c << "| d2H " << 
	      "| dTB1_new:" << d2H_new << "| dTB1_old:" << _d2H_TB1->matrix[r][c] << 
	      "| VTBold:" << V_TBold << "| VTBnew:" << V_TBnew << "| evapS_TB: " << evapS_TB << endl;
	  _d2H_TB1->matrix[r][c] = d2H_new;
	  // Mobile water
	  if(V_MWold > RNDOFFERR){
	    Frac_Esoil(atm, bsn, ctrl, V_MWold, V_MWnew, beta, 
		       _d2H_MW1->matrix[r][c], d2H_new, d2H_MWevap, Ts, r, c, 0);
	    if(abs(d2H_new)>100)
	      cout << r << " " << c << "| d2H " << 
		"| dMW1_new:" << d2H_new << "| dMW1_old:" << _d2H_MW1->matrix[r][c] << 
		"| VMWold:" << V_MWold << "| VMWnew:" << V_MWnew <<  "| evapS_MW: " << evapS_MW << endl;
	    _d2H_MW1->matrix[r][c] = d2H_new;
	  }
	}
	
	
	if(ctrl.sw_18O){
	  // Tightly-bound
	  Frac_Esoil(atm, bsn, ctrl, V_TBold, V_TBnew, beta, 
		     _d18O_TB1->matrix[r][c], d18O_new, d18O_TBevap, Ts, r, c, 1);
	  _d18O_TB1->matrix[r][c] = d18O_new;
	  // Mobile water
	  if(V_MWold > RNDOFFERR){
	    Frac_Esoil(atm, bsn, ctrl, V_MWold, V_MWnew, beta, 
		       _d18O_MW1->matrix[r][c], d18O_new, d18O_MWevap, Ts, r, c, 1);
	    _d18O_MW1->matrix[r][c] = d18O_new;
	  }
	}
	
    } else { // if no fractionation or very small E, or Age: only update E, weighted means
      if(ctrl.sw_2H)
	d2H_evap = evapS > RNDOFFERR ? 
	  (evapS_TB*_d2H_TB1->matrix[r][c] + evapS_MW*_d2H_MW1->matrix[r][c]) / evapS : -1000;
      if(ctrl.sw_18O)
	d18O_evap = evapS > RNDOFFERR ? 
	  (evapS_TB*_d18O_TB1->matrix[r][c] + evapS_MW*_d18O_MW1->matrix[r][c]) / evapS : -1000;
    }

    // Age: no fractionation, wighted sum of both domains
    if(ctrl.sw_Age)
      Age_evap = evapS > RNDOFFERR ? 
	(evapS_TB*_Age_TB1->matrix[r][c] + evapS_MW*_Age_MW1->matrix[r][c]) / evapS : 0;

    // --------------------------------------------------------------------------------
  } else { // L1-averaged values
    // Safeguard with minimum E value: avoids diverging d2H_evap values
    if(ctrl.sw_frac and V_new/V_old <0.999) {
      if(ctrl.sw_2H){
	Frac_Esoil(atm, bsn, ctrl, V_old, V_new, beta, 
		   _d2Hsoil1->matrix[r][c], d2H_new, d2H_evap, Ts, r, c, 0);
	_d2Hsoil1->matrix[r][c] = d2H_new;
      }
      if(ctrl.sw_18O){
	Frac_Esoil(atm, bsn, ctrl, V_old, V_new, beta, 
		   _d18Osoil1->matrix[r][c], d18O_new, d18O_evap, Ts, r, c, 1);
	_d18Osoil1->matrix[r][c] = d18O_new; 
      }
    } else { // If no fractionation (or negligible E), only update E
      if(ctrl.sw_2H)
	d2H_evap = _d2Hsoil1->matrix[r][c]; //evapS > RNDOFFERR ? _d2Hsoil1->matrix[r][c] : -1000;
      if(ctrl.sw_18O)
	d18O_evap = _d18Osoil1->matrix[r][c]; //evapS > RNDOFFERR ? _d18Osoil1->matrix[r][c] : -1000;
    }

    if(ctrl.sw_Age)
	Age_evap = _Agesoil1->matrix[r][c]; //evapS > RNDOFFERR ? _Agesoil1->matrix[r][c] : 0;	      
    
  }
    
  // Finally, update summed evap values (soil values already done)
  // flow from mobile to tightly-bound (if two--pore domain activated)
  if(ctrl.sw_2H){
    _d2HevapS_sum->matrix[r][c] = InputMix(etotS, _d2HevapS_sum->matrix[r][c], evapS, d2H_evap);
      
    if(ctrl.sw_TPD){
      d_old = _d2H_TB1->matrix[r][c];
      _d2H_TB1->matrix[r][c] = InputMix(V_TBnew, _d2H_TB1->matrix[r][c],
					min<double>(theta_MW1*d1 - V_TBnew, V_MWnew),
					_d2H_MW1->matrix[r][c]) ;
      if(abs(_d2H_TB1->matrix[r][c])>100)
	cout << r << " " << c << "| d2H " << 
	  "| dTB1_new:" << _d2H_TB1->matrix[r][c] << "| dTB1_old:" << d_old << 
	  "| dMW1:" << _d2H_MW1->matrix[r][c] << 
	  "| VTBold:" << V_TBold << "| VTBnew:" << V_TBnew <<  "| E_TB:" << evapS_TB << 
	  "| VMWold:" << V_MWold << "| VTBnew:" << V_MWnew <<  "| E_MW:" << evapS_MW << 
	  "| V_old:" << V_old << "| V_new:" << V_new <<  "| E:" << evapS <<
	  "| VMW:" << theta_MW1*d1 << "| theta_MW1:" << theta_MW1 <<
	  "| MWtoTB:" << min<double>(theta_MW1*d1 - V_TBnew, V_MWnew) << endl;
    }
  }
  if(ctrl.sw_18O){
    _d18OevapS_sum->matrix[r][c] = InputMix(etotS, _d18OevapS_sum->matrix[r][c], evapS, d18O_evap);
      
    if(ctrl.sw_TPD)
      _d18O_TB1->matrix[r][c] = InputMix(V_TBnew, _d18O_TB1->matrix[r][c],
					 min<double>(theta_MW1*d1 - V_TBnew, V_MWnew),
					 _d18O_MW1->matrix[r][c]) ;
  }
    
  if(ctrl.sw_Age){
    _AgeevapS_sum->matrix[r][c] =  InputMix(etotS, _AgeevapS_sum->matrix[r][c], evapS, Age_evap);
      
    if(ctrl.sw_TPD)
      _Age_TB1->matrix[r][c] = InputMix(V_TBnew, _Age_TB1->matrix[r][c],
					min<double>(theta_MW1*d1 - V_TBnew, V_MWnew),
					_Age_MW1->matrix[r][c]) ;
      
  }
// -----------------------------------------------------------------------------------
}

