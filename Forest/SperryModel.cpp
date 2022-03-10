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
 * SperryModel.cpp
 *
 * Created on: Apr 16, 2019
 *  Author: Aaron Smith
 */
#define ARMA_NO_DEBUG //disables armadillo bound checks for speed optimization
#include <armadillo>
#include "Forest.h"

using namespace arma;

int Forest::SperryModel(Basin &bas, Atmosphere &atm, Control &ctrl,
			REAL8 rootdepth,REAL8 Sold,REAL8 Keff, REAL8 psiae, REAL8 bclambda,
			REAL8 airTp,REAL8 airRH, REAL8 rho_a, REAL8 gamma, REAL8 ra,
			REAL8 poros, REAL8 thetar,REAL8 thetawp, REAL8 evap_a,REAL8 fA,
			REAL8 fB, REAL8 fC, REAL8 leavesurfRH,REAL8 leafRH,REAL8 &LET,
			REAL8 &LE, REAL8 &H,REAL8 &temp0,REAL8 &temp1,
		        REAL8 &temp2,REAL8 &temp3, UINT4 s, UINT4 r, UINT4 c){


  const REAL8 grav = 9.8;
  const REAL8 Vw = 18e-6; // partial molal volume of water m3 mol-1
  const REAL8 Rbar = 8.31446; // Universal gas constant in J K-1 mol-1
  REAL8 dt = ctrl.dt;
  REAL8 desdTs;
  REAL8 ra_t;
  REAL8 dleafRHdT = 0;
  REAL8 dleafRHdpsi_l = 0;
  REAL8 dgcdfgspsi = 0;
  REAL8 dEdlwp = 0;
  REAL8 dEdT = 0;
  REAL8 dLETdlwp = 0;
  REAL8 dLETdT = 0;
  REAL8 gc = 0;
  REAL8 lwp_min, lwp_max;
  REAL8 albedo,emissivity,BeerK,LAI;
  REAL8 E = 0;
  REAL8 ea, es;
  REAL8 maxAv;
  REAL8 lambda;
  REAL8 temp = 0;
  REAL8 LEmax;

  REAL8 sperry_c = _species[s].sperry_c; //Weibull parameter (exponent)
  REAL8 sperry_d = _species[s].sperry_d; //Weibull parameter (denominator)
  REAL8 sperry_ks = _species[s].sperry_Kp;//maximum vegetation conductance
  REAL8 Plnt_Dyn= 0; // fail criteria - switchs to 1 if solution fails to solve

  REAL8 dgs_psidpsi;
  REAL8 RAI = _species[s]._RootMass->matrix[r][c] * _species[s].SRA;
  REAL8 gp = 0;
  REAL8 gsr = 0;
  REAL8 gsrp = 0;
  REAL8 dgp = 0;
  REAL8 dgsr = 0;
  REAL8 denfac = 0;
  /***
  * state variables:
  * x[0]: S - (degree of saturation at time t+1)
  * x[1]: psi_s - soil water potential
  * x[2]: psi_l - leaf water potential
  * x[3]: Ts - Leaf temperature
  ***/

  albedo = _species[s].albedo;
  emissivity = _species[s].emissivity;
  BeerK = _species[s].KBeers;
  LAI = _species[s]._LAI->matrix[r][c];

  ea = SatVaporPressure(airTp) * airRH;

  if(ctrl.toggle_sm == 0){ 
    lwp_min = _species[s].lwp_min;
    lwp_max = _species[s].lwp_max;
  } else {
    lwp_min = psiae * 0.0098 / powl(Saturation(thetawp,thetar,poros),bclambda) * _species[s].lwp_min;
    lwp_max = psiae * 0.0098 / powl(Saturation(thetawp,thetar,poros),bclambda) * _species[s].lwp_max;
  }

  int k = 0;

  colvec x(4);
  colvec deltax(4);
  colvec F(4);
  mat J = zeros<mat>(4,4);
  //-------- SOLVE THE SOILS IN LAYER 1 ---------------------------------
  //provide initial guess  for loop
  x[0] = Sold;
  x[1] = psiae * 0.0098 / powl(x[0], bclambda);
  x[2] = x[1];
  x[3] = airTp;

  maxAv = (poros - thetar);

  CalculateCanopyConduct(bas, atm, ctrl, x[2],lwp_max,lwp_min,dgcdfgspsi, s, r, c);
  do {

      if(ctrl.toggle_plant_hydr == 0 or Plnt_Dyn == 1)
      x[2] = x[1];

    lambda = x[3] < 0 ? lat_heat_vap + lat_heat_fus : lat_heat_vap;

    //lwp_min - larger value	  //lwp_max - smaller value
    if(ctrl.toggle_sm == 0){ // nonlinear model //TODO CLEAN UP IN CONSTANDFUNCS
      gc = dgcdfgspsi * Calculate_gs_lwp_nonlinear(x[2],lwp_max,lwp_min);
      dgs_psidpsi = - dgcdfgspsi * lwp_max * powl(x[2]/lwp_min, lwp_max) /
	(x[2] * ( powl(x[2]/lwp_min, lwp_max) + 1) * ( powl(x[2]/lwp_min, lwp_max) + 1));
    } else { // linear model
      gc = dgcdfgspsi*Calculate_gs_lwp_linear(x[2],lwp_max,lwp_min);
      dgs_psidpsi = dgcdfgspsi * ((x[2] > lwp_min) || (x[2] < lwp_max) ? 0 : -1/(lwp_min - lwp_max));
    }
      
    if (gc < 1e-13)
      gc = 1e-13;
      
    ra_t = ra + (1 / gc);

    x[0] = x[0] < RNDOFFERR ? 0.01 : x[0]; //check saturation

    temp = -x[2] * rho_w * grav * Vw / (Rbar*(x[3]+273.15));
    leafRH = (temp > -708.4) ? std::min<REAL8>(1,expl(temp)) : 0;	

    if(abs(evap_a)<RNDOFFERR){
      LE = 0.0;
    } else if(LatHeatCanopy(bas, atm, leavesurfRH, ra, x[3], r, c) > 0.0){
      LE = 0.0;
    } else {
      LEmax = evap_a * (rho_w * lambda); 
      LE = LatHeatCanopy(bas, atm, leavesurfRH, ra, x[3], r, c) < 
	LEmax ? LEmax : LatHeatCanopy(bas, atm, leavesurfRH, ra, x[3], r, c);
    }
    LET = LatHeatCanopy(bas, atm, leafRH, ra_t, x[3], r, c);
    H = SensHeatCanopy(atm, ra, x[3], r, c);
    E = std::max<REAL8>(0.0,-LET / (rho_w * lambda));

    F[0] = (x[0] - Sold) * maxAv * rootdepth / dt + E;
    F[1] = psiae * 0.0098 / powl(x[0], bclambda) - x[1];
    F[2] = 0;
    F[3] = NetRadCanopy(atm, x[3], emissivity, albedo, BeerK, LAI, r, c) + LE + H + LET;
    // --------------------------------------------------------------
    // Plant Hydraulics - Sperry Stuff (only if Plant_Hyd is on)
    //    the conversion coefficients are: 
    //          1e6 in the numerator to convert Keff from m/s to micrometers/s
    //          1000000 in the numerator to convert Pa (rootdepth*grav*rho_w) to MPa in the denominator
    // --------------------------------------------------------------	
    if(ctrl.toggle_plant_hydr and Plnt_Dyn == 0){
      gsr = Keff * 1e6 * powl(x[0], 2 * bclambda + 3) * sqrt(RAI)
	/ ( PI * rootdepth * grav * rho_w) * 1000000;
	
      gp = (x[2] < 0) ? sperry_ks : sperry_ks * expl(-powl(x[2] / sperry_d, sperry_c));
	
      denfac = gsr + gp * LAI;

      gsrp = denfac < 1e-12 ? 0 : LAI * gsr * gp / denfac;    //prevent division by zero

      F[2] = gc < 1e-12 ? 0 : gsrp * std::max<REAL8>(0,(x[2] - x[1])) - E * 1e6; //if no demand, shutdown supply
    }

    //*************************************************
    // Stuff for Jacobian calculation of F
    //*************************************************
    dleafRHdT = leafRH *  x[2] * rho_w * grav * Vw / (Rbar*powl((x[3]+273.15),2));

    dleafRHdpsi_l = leafRH == 1 ? 0 : - rho_w * grav * Vw * leafRH  /(Rbar*(x[3]+273.15));

    es = SatVaporPressure(x[3]);
    desdTs = es * ((17.3 / (x[3] + 237.3)) - 17.3*x[3] / (powl(x[3] + 237.3, 2)));      

    if ((dgcdfgspsi == 0) || (dgs_psidpsi == 0)){
      dLETdlwp = 0;
      if(gc < 1e-12)
	if(ctrl.toggle_sm == 1 && x[2] > lwp_min)
	  LET = E = 0;
    } else{
      dLETdlwp = spec_heat_air * rho_a * (ea * (dgs_psidpsi/powl(gc,2) ) -
					  es * (dleafRHdpsi_l * ra_t + leafRH * dgs_psidpsi/powl(gc,2)) ) /  (gamma * powl(ra_t,2));
    }
      
    dLETdT = - rho_a * spec_heat_air / (ra_t * gamma) * (desdTs*leafRH + es*dleafRHdT);
      
    dEdlwp = - dLETdlwp / (rho_w * lambda);
	
    dEdT = - dLETdT / (rho_w * lambda);
    
    //------------------------------------------------
    // Fill out the jacobian
    //------------------------------------------------
    J(0, 0) = rootdepth * maxAv / dt;
    J(0, 1) = (ctrl.toggle_plant_hydr == 0 or Plnt_Dyn == 1) ? dEdlwp : 0; //toggle between hydr on and off
    J(0, 2) = (ctrl.toggle_plant_hydr == 1 and Plnt_Dyn == 0) ? dEdlwp : 0; //toggle between hydr on and off	
    J(0, 3) = (E==0) ? 0 : dEdT;
    
    J(1,0) = -bclambda * psiae * 0.0098 * powl(x[0], -(bclambda + 1));
    J(1,1) = -1;

    if (gc < 1e-12 || (ctrl.toggle_sm == 1 && x[2] > lwp_min) || ctrl.toggle_plant_hydr==0 || Plnt_Dyn == 1) {
      J(2, 0) = J(2, 1) = J(2, 3) = 0;
      J(2, 2) = 1;
    } else { 
      dgsr = Keff * 1e6 * powl(x[0], 2 * bclambda + 2) * sqrt(RAI)* (2 * bclambda + 3) / ( PI * rootdepth * grav * rho_w)	* 1000000;
      dgp = -gp * sperry_c * powl(x[2] / sperry_d, sperry_c) / x[2];
      J(2, 0) = (gp * LAI * (dgsr * denfac - gsr * dgsr) / powl(denfac , 2)) * (x[2] - x[1]); //quotient rule dF[2]/dx[0]
      J(2, 1) = -gsrp;
      J(2, 2) = ((gsr * LAI * dgp * denfac - gsr * powl(LAI,2) * dgp * gp) /
		 powl(denfac,2)) * (x[2] - x[1]) + gsrp - dEdlwp * 1e6;  //quotient rule dF[2]/dx[2]
      J(2, 3) = (-dEdT);
    }
    
    J(3, 1) = (ctrl.toggle_plant_hydr == 0 or Plnt_Dyn == 1) ? dLETdlwp : 0;
    J(3, 2) = (ctrl.toggle_plant_hydr == 1 and Plnt_Dyn == 0) ? dLETdlwp : 0;	
    J(3, 3) = fA * powl(x[3] + 273.2, 3) + fB * desdTs * leavesurfRH + fC + dLETdT;

    //------------------------------------------------
    // End of the jacobian
    //------------------------------------------------
    // solve system
    if (!solve(deltax, J, -F)) {
      cout << "Singular Jacobian found in Newton solver for canopy balance.\n";
      if(Plnt_Dyn==1){
	k = MAX_ITER;
      } else {
	cout << "Plant hydraulics failed - run with hydrualics off. \n";
	Plnt_Dyn = 1;
	F[0] = 1; F[1] = 1; F[2] = 0 ;F[3] = 1;
	x[0] = Sold;
	x[1] = psiae * 0.0098 / powl(x[0], bclambda); // 0.0098 convert from m to MPa
	x[2] = x[1];
	x[3] = airTp;	  
	k = 0;
      }
    } else {
      x += deltax;
      k++;
      
      x[2] = x[2] < x[1] ? x[1] : x[2]; //no negative flow (leaf to root)
	
      if ((x[2] < 0) & (x[1] < 0) & (Plnt_Dyn == 0)){
	// Resets the solver and runs with plant dynamics off
	cout << "Plant hydraulics failed - run with hydrualics off. \n";
	Plnt_Dyn = 1;
	F[0] = 1; F[1] = 1; F[2] = 0 ;F[3] = 1;
	x[0] = Sold;
	x[1] = psiae * 0.0098 / powl(x[0], bclambda); // 0.0098 convert from m to MPa
	x[2] = x[1];
	x[3] = airTp;	  
	k = 0;
      } else if((x[2] < 0) & (x[1] < 0) & (Plnt_Dyn == 1)){
	    k = MAX_ITER;
      }
    } //else - singular jacobian
  } while (norm(F, "inf") > 0.000001 && k < MAX_ITER);

  if (k >= MAX_ITER)
    std::cout
      << "WARNING: non-convergence in canopy energy balance at cell row: "
      << r << " col: " << c << " closure err: " << norm(deltax, 2)
      << endl;

  LET = LatHeatCanopy(bas, atm, leafRH, ra_t, x[3], r, c);
  if(abs(evap_a)<RNDOFFERR){
    LE = 0.0;
  } else {
    LEmax = evap_a * (rho_w * lambda); 
    LE = LatHeatCanopy(bas, atm, leavesurfRH, ra, x[3], r, c) < LEmax ? LEmax : LatHeatCanopy(bas, atm, leavesurfRH, ra, x[3], r, c);
  }
  H = SensHeatCanopy(atm, ra, x[3], r, c);

  temp0 = x[0];
  temp1 = x[1];
  temp2 = x[2];
  temp3 = x[3];

  return EXIT_SUCCESS;
}

