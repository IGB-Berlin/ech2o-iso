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
 *    Marco Maneta, Sylvain Kuppel, Aaron Smith
 *******************************************************************************/
/*
 * Fractionation_Estorage.cpp
 *
 *  Created on: Apr 5, 2017
 *      Author: Sylvain Kuppel
 */

#include "Tracking.h"

int Tracking::Frac_Estorage(Atmosphere &atm, Basin &bsn, Control &ctrl,
			    REAL8 V_old, REAL8 V_new, REAL8 &di_old,
			    REAL8 &di_new, REAL8 &di_evap, REAL8 &beta,REAL8 &Ts, 
			    int issoil, int r, int c, int iso){

  //*****************************************************************************************
  //This is a general fractionation for all open water/storage
  //*****************************************************************************************
  REAL8 Ta = atm.getTemperature()->matrix[r][c] + 273.15 ; // Air temperature (K)
  REAL8 ha = atm.getRelativeHumidty()->matrix[r][c]; // Atmospheric relative humidity (fraction)
  REAL8 es_s, ea_s;
  REAL8 bclambda,psiae,th_r,th_s,d1;
  es_s = ea_s = bclambda = psiae = th_r = th_s = d1 = 0;
  if(issoil==1){
    es_s = SatVaporPressure(Ts-273.15) ; //saturated vapor pressure in the atmosphere (Pa)
    ea_s = SatVaporPressure(Ta-273.15) ; // saturated vapor pressure at the surface (Pa)
  }
  
  REAL8 ha_p;			// Corrected relatvie air humidity above the surface (fraction)
  REAL8 hs; 			// Soil vapor saturation at the surface (fraction)
  REAL8 f; 			// Water loss fraction after evaporation (fraction)
  REAL8 alpha_p = 0; 		// equilibrium isotope fractionation factor (fraction)
  REAL8 eps_p; 			// equilibrium isotope fractionation factor (per mil)
  REAL8 eps_k = 0; 		// kinetic isotope fractionation factor (per mil)
  REAL8 eps; 			// total isotope fractionation factor (per mil)
  REAL8 di_atm = 0; 		// Isotopic signatures (permil)
  REAL8 di_s; 			// Limiting isotopic composition (per mil)
  REAL8 m; 			// Calculation factor (-)
  REAL8 n;			// Parameter translating dominant water transport mode (-)
  REAL8 DiD;
  REAL8 u;  

  ha_p = (issoil == 1) ? ( ha*ea_s/es_s > 1.0 ? 1.0 : ha*ea_s/es_s) : (ha + 1)/2;  
  //*****************************************************************************************
  //--- Humidity in the soil (if this is soil) ---------------------------------------------
  //*****************************************************************************************
  if(issoil){
    if(ctrl.toggle_hs == 0) {             // Just 1
      hs = 1;
    } else if (ctrl.toggle_hs == 1) {     // Corrected h for surface temperature and beta
      hs = beta + (1-beta)*ha_p;
    } else if (ctrl.toggle_hs == 2) {     // Sorderberg et al. (2012), orginially for dry soils (departs from EcH2O's evaporation conceptualization of evap!)
      th_r = bsn.getSoilMoistR1()->matrix[r][c];
      th_s = bsn.getPorosityL1()->matrix[r][c];
      psiae = bsn.getPsiAEL1()->matrix[r][c];
      bclambda = bsn.getBClambdaL1()->matrix[r][c];
      d1 = bsn.getSoilDepth1()->matrix[r][c];
      hs = expl(psiae*powl((V_old/d1-th_r)/(th_s-th_r),-bclambda)*18.0145/(8.3145*Ts*1000)); // waterpot*molarmass/(R*Ts*waterdensity)
    } else {
      std::cout << "Wrong option in the surface humidity toggle switch for fractionation. Please verify the configuration file" << std::endl;
      exit(EXIT_FAILURE);  
    }
    hs = hs < ha_p ? ha_p : hs;          // soil humidity cannot be higher than atmosphere (math errors)
  }
  //*****************************************************************************************
  //--- Fractionation variables -------------------------------------------------------------
  //*****************************************************************************************
  // ****  Horita and Wesolowski (1994) **** 
  if(iso == 0) // deuterium
    alpha_p = expl((1158.8*powl(Ta,3)*1e-9 - 1620.1*powl(Ta,2)*1e-6 + 794.84*Ta*1e-3 - \
		    161.04 + 2.9992*1e9/powl(Ta,3))/1000);
  else if (iso == 1) // oxygen 18
    alpha_p = expl((-7.685 + 6.7123*1000/Ta - 1.6664*1e6/powl(Ta,2) +
		    0.35041*1e9/powl(Ta,3))/1000);
  
  // **** Equilibrium enrichment (Skrzypek et al., 2015) ****
  eps_p = (alpha_p - 1) *1000;
  
  // **** Atmospheric vapour estimation (Gat, 1995) + (Gibson and Reid, 2014) ****
  if(iso == 0) // deuterium
    di_atm = (atm.getd2Hprecip()->matrix[r][c] - eps_p)/ alpha_p;
  else if(iso ==1) // oxygen 18
    di_atm = (atm.getd18Oprecip()->matrix[r][c] - eps_p)/ alpha_p;
  
  // **** Water transport mode: from diffusive (=1, dry soil) to turbulent (=0.5, water body) ****
  if(issoil == 1){
    if(ctrl.toggle_n == 0) {
      n = 1;
    } else if (ctrl.toggle_n == 1) {
      // Mathieu and Bariac (1996) + Braud et al. (2005)
      n = 1 - 0.5*(V_old/d1-bsn.getSoilMoistR1()->matrix[r][c])/
        (bsn.getPorosityL1()->matrix[r][c]-bsn.getSoilMoistR1()->matrix[r][c]);
    } else {
      std::cout << "Wrong option in the soil fractionation n toggle switch. Please verify the configuration file" << std::endl;
      exit(EXIT_FAILURE);  
    }
  } else { //Open water body
    n = 0.5;
  }

  // **** Kinetic fractionation factor epsilon_k ****
  if(ctrl.toggle_ek == 0) {                  		// Value of Di/D from Merlivat (1978)
    DiD   = (iso==0) ? 0.9757 : 0.9727;      		// Deuterium / Oxygen-18
    eps_k = (1-ha_p)*(1-DiD)*1000*n;
  } else if (ctrl.toggle_ek == 1) {          		// Value of Di/D from Vogt (1976)
    DiD   = (iso==0) ? 0.9877 : 0.9859;      		// Deuterium / Oxygen-18
    eps_k = (1-ha_p)*(1-DiD)*1000*n;
  } else if (ctrl.toggle_ek == 2) {          		// From Merlivat and Jouzel (1979) adapted by Haese et al. (2013)
    u = atm.getWindSpeed()->matrix[r][c];    		// wind speed (m.s-1)
    if(iso==0) // deuterium
      eps_k = u > 7.0 ? (1-ha_p)*0.88*(0.285*u+0.82) : (1-ha_p)*0.88*6;
    else if (iso==1) // oxygen 18
      eps_k = u > 7.0 ? (1-ha_p)*(0.285*u+0.82) : (1-ha_p)*6;
  }

  // **** Generalized following Good et al. (2014) ****
  eps = eps_p + eps_k;					// Gibson and Reid (2010)
  di_s = (ha_p*di_atm + eps) / (ha_p - eps/1000);  	// (Gat and Levy, 1978) + (Gat, 1981)
  m = (1 - ha_p + eps_k/1000) < RNDOFFERR ? 0 : 
	(ha_p - eps/1000) / (1 - ha_p + eps_k/1000);    // (Welhan and Fritz, 1977) + (Allison and Leaney, 1982)

  //*****************************************************************************************
  //--- Isotopic Mixing ---------------------------------------------------------------------
  //*****************************************************************************************  
  REAL8 evap = (V_old - V_new);
  REAL8 Vavg = (V_old + V_new)/2;
  REAL8 Vdiff = (V_new - V_old);
  REAL8 corr = 0; 					// Correction factor (mass-balance closure)
  if(V_old > RNDOFFERR){
    f = V_new/V_old;					// Evaporative loss fraction
    // New isotopic composition
    di_new = di_s - (di_s - di_old) * (powl(f,m) > 1 ? 1 : powl(f,m));
    // Isotopic signature of evaporated water
    di_evap = (1 - ha_p + eps_k/1000) < RNDOFFERR ? di_new : (di_new - ha_p*di_atm - eps)/ (1 - ha_p + eps_k/1000);
    if(abs(di_evap - di_new) > RNDOFFERR){
      corr = -((Vdiff*di_old/2 - Vavg*di_old - evap*(ha_p*di_atm+eps)/(1-ha_p+eps_k/1000)) + 
	      di_new*(Vavg +Vdiff/2 + evap/(1-ha_p+eps_k/1000)))/
     	      (Vavg + Vdiff/2 + evap/(1-ha_p+eps_k/1000));
      di_new = di_new + corr;
      di_evap = (1 - ha_p + eps_k/1000) < RNDOFFERR ? di_new :(di_new - ha_p*di_atm - eps)/ (1 - ha_p + eps_k/1000);
      if(di_evap > di_new){      // Check here if di_evap > di_new (not possible)
	di_new = di_old;
	di_evap = di_old;
      }
    } //close (di_evap-di_new)
  }

  return EXIT_SUCCESS;
  
}
