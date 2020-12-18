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
 * Fractionation_Echan.cpp
 *
 *  Created on: Apr 5, 2017
 *      Author: Sylvain Kuppel
 */

#include "Tracking.h"

int Tracking::Frac_Echan(Atmosphere &atm, Basin &bsn, Control &ctrl,
			 REAL8 V_old, REAL8 V_new,
			 REAL8 &di_old, REAL8 &di_new, REAL8 &di_evap,
			 REAL8 &Ts, int r, int c, int iso){
  
  REAL8 Ta = atm.getTemperature()->matrix[r][c] + 273.15 ; // Air temperature (K)
  REAL8 ha = atm.getRelativeHumidty()->matrix[r][c]; // Atmospheric relative humidity (fraction)
  //  REAL8 es_s = SatVaporPressure(Ts-273.15) ; //saturated vapor pressure in the atmosphere (Pa)
  //  REAL8 ea_s = SatVaporPressure(Ta-273.15) ; // saturated vapor pressure at the surface (Pa)
  
  REAL8 ha_p; // Corrected relatvie air humidity above the surface (fraction)
  REAL8 f; // Water loss fraction after evaporation (fraction)
  REAL8 alpha_p = 0; // equilibrium isotope fractionation factor (fraction)
  REAL8 eps_p; // equilibrium isotope fractionation factor (per mil)
  REAL8 eps_k = 0; // kinetic isotope fractionation factor (per mil)
  REAL8 eps; // total isotope fractionation factor (per mil)
  REAL8 di_atm = 0; // Isotopic signatures (permil)
  REAL8 di_s; // Limiting isotopic composition (per mil)
  REAL8 m; // Calculation factor (-)
  REAL8 n; // Parameter translating dominant water transport mode (-)
  
  // Corrected relative humidity at the surface
  ha_p = (ha + 1)/2;
  //ha_p = std::min<double>(ha*ea_s/es_s,1.0);

  n = 0.5;
  
  // Horita and Wesolowski (1994)
  if(iso == 0) // deuterium
    alpha_p = expl((1158.8*powl(Ta,3)*1e-9 - 1620.1*powl(Ta,2)*1e-6 + 794.84*Ta*1e-3 - \
		    161.04 + 2.9992*1e9/powl(Ta,3))/1000);
  else if (iso == 1) // oxygen 18
    alpha_p = expl((-7.685 + 6.7123*1000/Ta - 1.6664*1e6/powl(Ta,2) +
		    0.35041*1e9/powl(Ta,3))/1000);
  
  // Skrzypek et al. (2015)
  eps_p = (1 - 1/alpha_p)*1000;
  
  // (Gat, 1995) + (Gibson and Reid, 2014)
  if(iso == 0) // deuterium
    di_atm = (atm.getd2Hprecip()->matrix[r][c] - eps_p)/ alpha_p;
  else if(iso ==1) // oxygen 18
    di_atm = (atm.getd18Oprecip()->matrix[r][c] - eps_p)/ alpha_p;

  
  // Kinetic fractionation factor epsilon_k
  if(ctrl.toggle_ek == 0) {
    // Value of Di/D from Merlivat (1978)
    if(iso==0) 
      eps_k = (1-ha_p)*(1-0.9757)*1000*n;
    else if(iso==1) 
      eps_k = (1-ha_p)*(1-0.9727)*1000*n;

  } else if (ctrl.toggle_ek == 1) {
    // Value of Di/D from Vogt (1976)
    if(iso==0) // deuterium
      eps_k = (1-ha_p)*(1-0.9877)*1000*n; 
    else if(iso==1) // oxygen 18
      eps_k = (1-ha_p)*(1-0.9859)*1000*n; 

  } else if (ctrl.toggle_ek == 2) {
    REAL8 u = atm.getWindSpeed()->matrix[r][c]; // wind speed (m.s-1)
    // From Merlivat and Jouzel (1979) adapted by Haese et al. (2013)
    if(iso==0) // deuterium
      eps_k = u > 7.0 ? (1-ha_p)*0.88*(0.285*u+0.82) : (1-ha_p)*0.88*6;
    else if (iso==1) // oxygen 18
      eps_k = u > 7.0 ? (1-ha_p)*(0.285*u+0.82) : (1-ha_p)*6;
  }

  // --- Generalized following Good et al. (2014) -------------
  // Gibson and Reid (2010)
  eps = eps_p + eps_k;
  // (Gat and Levy, 1978) + (Gat, 1981)
  di_s = (ha_p*di_atm + eps) / (ha_p - eps/1000);
  // (Welhan and Fritz, 1977) + (Allison and Leaney, 1982)
  m = (1 - ha_p + eps_k/1000) < RNDOFFERR ? 0 : (ha_p - eps/1000) / (1 - ha_p + eps_k/1000);
  // ----------------------------------------------------------
  
  // Evaporative loss fraction
  f = V_new/V_old;
  
  // (Hamilton et al., 2005)
  // New isotopic signature in topsoil
  di_new = di_s - (di_s - di_old) * (powl(f,m) > 1 ? 1 : powl(f,m));

  REAL8 evap = (V_old - V_new);
  REAL8 Vavg = (V_old + V_new)/2;
  REAL8 Vdiff = (V_new - V_old);
  
  // Isotopic signature of evaporated water
  di_evap = (1 - ha_p + eps_k/1000) < RNDOFFERR ? di_new : (di_new - ha_p*di_atm - eps)/ (1 - ha_p + eps_k/1000);

  REAL8 corr = 0;
  
  if(abs(di_evap - di_new) > RNDOFFERR){
     corr = -((Vdiff*di_old/2 - Vavg*di_old - evap*(ha_p*di_atm+eps)/(1-ha_p+eps_k/1000)) + 
	di_new*(Vavg +Vdiff/2 + evap/(1-ha_p+eps_k/1000)))/
	(Vavg + Vdiff/2 + evap/(1-ha_p+eps_k/1000));
     di_new = di_new + corr;
     di_evap = (1 - ha_p + eps_k/1000) < RNDOFFERR ? di_new :(di_new - ha_p*di_atm - eps)/ (1 - ha_p + eps_k/1000);
  }
  
  return EXIT_SUCCESS;
  
}
