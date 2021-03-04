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
 * SolveEnergyBalance.cpp
 *
 *  Created on: Jul 9, 2010
 *      Author: Marco.Maneta
 */
#define ARMA_NO_DEBUG //disables armadillo bound checks for speed optimization
#include <armadillo>
#include "Forest.h"

using namespace arma;

UINT4 Forest::SolveCanopyEnergyBalance(Basin &bas, Atmosphere &atm, Control &ctrl,
				       REAL8 theta_r1, REAL8 theta_r2, REAL8 theta_r3,
				       REAL8 rootdepth,
				       REAL8 psiae, REAL8 bclambda, REAL8 ra, 
				       REAL8 &DelCanStor, REAL8 &evap_a, 
				       REAL8 &transp_a, REAL8 &netR_a, REAL8 &le_a,
				       UINT4 s, UINT4 r, UINT4 c) {

  //some constants
  const REAL8 grav = 9.8;
  const REAL8 Vw = 18e-6; // partial molal volume of water m3 mol-1
  const REAL8 Rbar = 8.31446; // Universal gas constant in J K-1 mol-1
  //energy balance parameters
  REAL8 dt = ctrl.dt;
  REAL8 fA, fB, fC; //pooling factors
  REAL8 temp = 0;
  REAL8 rho_a; //density of air
  REAL8 airRH; //air humidity
  REAL8 airTp; // air temperature
  REAL8 es, ea; // saturated vapor pressure
  REAL8 desdTs; // derivative of saturation vapor pressure function with respect to Ts
  REAL8 emissivity; //canopy emissivity
  REAL8 albedo; //canopy albedo
  REAL8 LAI;
  REAL8 BeerK; //Beers-Lambert coefficient
  REAL8 LE, LEmax, H;
  REAL8 LE1, LE2, LE3;
  REAL8 LET1, LET2, LET3;
  REAL8 H1, H2, H3;
  REAL8 z; //terrain height
  REAL8 gamma; //psychrometric constant
  REAL8 f1, f2, f3; //root fractions

  REAL8 lambda = lat_heat_vap;
  REAL8 ra_t; //resistance to transpiration ra_t = ra + 1/gc
  REAL8 LET; //latent heat of transpiration
  REAL8 CanStor = 0;
  REAL8 MaxCanStor = 0;

  REAL8 leafRH; //soil relative humidty use in teh calculation of soil vapor pressure for latent heat exchanges
  REAL8 leavesurfRH; //relative humidity of the leave surface. 1 when leave is saturated with intercepted water, airRH when no water
  REAL8 dleafRHdT = 0;
  REAL8 dleafRHdpsi_l = 0;
  //REAL8 dgcdlwp = 0;
  REAL8 dEdlwp = 0;
  REAL8 dEdT = 0;
  REAL8 dLETdlwp =0;
  REAL8 dLETdT =0;

  // variables for Sperry's model
  REAL8 Sold = 0;  // Soil Saturation at beginning of t
  REAL8 E = 0; // Temporary calculatio of transpiration
  REAL8 maxAv = 0; // Equivalent (weighted) maximum available water 
  
  REAL8 gc = 0;
  REAL8 lwp_min, lwp_max; // lower and higher boundaries for linear lwp stomatal model
  REAL8 dgcdfgspsi = 0;

  /* ----------------------------------------- */
  /* IF DEPTH DEPENDENT SOIL */
  REAL8 bclambda1 = bclambda;  
  REAL8 bclambda2 = bas.getBClambdaL2()->matrix[r][c];
  REAL8 bclambda3 = bas.getBClambdaL3()->matrix[r][c];
  REAL8 psiae1 = psiae;
  REAL8 psiae2 = bas.getPsiAEL2()->matrix[r][c];
  REAL8 psiae3 = bas.getPsiAEL3()->matrix[r][c];
  REAL8 d1 = bas.getSoilDepth1()->matrix[r][c];
  REAL8 d2 = bas.getSoilDepth2()->matrix[r][c];
  REAL8 d3 = bas.getSoilDepth()->matrix[r][c] - d1 - d2;

  REAL8 Sold1 = 0;
  REAL8 Sold2 = 0;
  REAL8 Sold3 = 0;
  REAL8 x1,x10,x11,x12,x2,x20,x21,x22;
  /* ----------------------------------------- */

  //  cout << "bclambda (1-2-3): " << bclambda << " | " << bclambda1 << " | " << bclambda2 << " | " << bclambda3 << endl;

  REAL8 theta_wp = 0; //wilting point
  REAL8 theta1 = bas.getSoilMoist1()->matrix[r][c];
  REAL8 theta2 = bas.getSoilMoist2()->matrix[r][c];
  REAL8 theta3 = bas.getSoilMoist3()->matrix[r][c];
  //REAL8 theta = 0;
  REAL8 poros1 = bas.getPorosityL1()->matrix[r][c];
  REAL8 poros2 = bas.getPorosityL2()->matrix[r][c];
  REAL8 poros3 = bas.getPorosityL3()->matrix[r][c];

  UINT4 nsp = getNumSpecies();

  if (s == nsp - 1) //for bare soil, water reaching the ground is pp times its proportion of the cell
    evap_a = transp_a = netR_a = le_a = 0;
  else {

    airTp = atm.getTemperature()->matrix[r][c];
    rho_a = AirDensity(airTp); //kgm-3
    z = bas.getDEM()->matrix[r][c];
    gamma = PsychrometricConst(101325, z);
    airRH = atm.getRelativeHumidty()->matrix[r][c];

    ea = SatVaporPressure(airTp) * airRH;


    albedo = _species[s].albedo;
    emissivity = _species[s].emissivity;
    BeerK = _species[s].KBeers;
    LAI = _species[s]._LAI->matrix[r][c];

    lwp_min = _species[s].lwp_min;
    lwp_max = _species[s].lwp_max;

    // Root fractiona
    f1 = _species[s]._rootfrac1->matrix[r][c];
    f2 = _species[s]._rootfrac2->matrix[r][c];
    f3 = 1 - f1 - f2;

    CanStor = getIntercWater(s, r, c);
    MaxCanStor = getMaxCanopyStorage(s, r, c);

    leafRH = 1; //min<REAL8>(1.0,Calculate_gs_theta(theta, fc, _species[s].WiltingPoint, 2.0)); 
    //calculates soil pore relative humidity
    leavesurfRH = airRH + ((1 - airRH) / MaxCanStor) * CanStor;

    fA =  -4 * emissivity * stefboltz;//pools together net radiation factors
    fB = (-1 / (ra * gamma)) * rho_a * spec_heat_air; // pools together the latent heat factors
    fC = (-1 / (ra)) * rho_a * spec_heat_air; // pools together the sensible heat factors

    evap_a = CanStor < MaxCanStor ? -CanStor / ctrl.dt * powl((CanStor / MaxCanStor), 0.6) : -CanStor / ctrl.dt;

    // Calculation of conductances necessary for the implementation of Sperry's model
    // to calculate leaf water potential.
    
    // Soil to root conductance. 
    // Adapted from Rodriguez-Iturbe and Porporato (eq 6.4, page 181) for units of hydraulic head
    // Here using a weighted average because porosity changes with depth
    //Sold = f1* (theta1 - theta_r1) / (poros1 - theta_r1) + f2* (theta2 - theta_r2) / (poros2 - theta_r2)+ f3* (theta3 - theta_r3) / (poros3 - theta_r3);
    // Maximum available water
    //maxAv = (poros1 - theta_r1) * (poros2 - theta_r2) * (poros3 - theta_r3) /(f1*(poros2-theta_r2)*(poros3-theta_r3)+f2*(poros1-theta_r1)*(poros3-theta_r3)+f3*(poros1-theta_r1)*(poros2-theta_r2));

    // EDIT: this formulation makes things much more simple
    maxAv = f1*(poros1 - theta_wp) + f2*(poros2-theta_wp)+f3*(poros3-theta_wp);
    Sold = (f1*std::max(0.0,(theta1-theta_wp))+f2*std::max(0.0,(theta2-theta_wp))+f3*std::max(0.0,(theta3-theta_wp))) / maxAv ;
		
    int k = 0;
    
    /***
     * state variables:
     * x[0]: S - (degree of saturation at time t+1)
     * x[1]: psi_s - soil water potential
     * x[3]: Ts - Leaf temperature
     ***/

    colvec x(3);
    colvec deltax(3);
    colvec F(3);
    mat J = zeros<mat>(3,3);

    if(ctrl.sw_ddSoilPar){
      //      cout << "Depth dependent Sperry is running" << endl;
      Sold1 = std::min<REAL8>(1.0,std::max<REAL8>(0.0,(theta1 - theta_wp) / (poros1 - theta_wp)));
      Sold2 = std::min<REAL8>(1.0,std::max<REAL8>(0.0,(theta2 - theta_wp) / (poros2 - theta_wp)));
      Sold3 = std::min<REAL8>(1.0,std::max<REAL8>(0.0,(theta3 - theta_wp) / (poros3 - theta_wp)));

      SperryModel(bas, atm, ctrl,d1,Sold1,psiae1,bclambda1,airTp,airRH,rho_a,gamma,ra,poros1,theta_wp,evap_a,fA,fB,fC,
		  leavesurfRH,leafRH,LET1,LE1,H1,x10,x20,s,r,c);

      SperryModel(bas, atm, ctrl,d2,Sold2,psiae2,bclambda2,airTp,airRH,rho_a,gamma,ra,poros3,theta_wp,evap_a,fA,fB,fC,
		  leavesurfRH,leafRH,LET2,LE2,H2,x11,x21,s,r,c);

      SperryModel(bas, atm, ctrl,d3,Sold3,psiae3,bclambda3,airTp,airRH,rho_a,gamma,ra,poros2,theta_wp,evap_a,fA,fB,fC,
		  leavesurfRH,leafRH,LET3,LE3,H3,x12,x22,s,r,c);

      x1 = f1*x10 + f2*x11 + f3*x12;
      x2 = f1*x20 + f2*x21 + f3*x22;
      lambda = x2 < 0 ? lat_heat_vap + lat_heat_fus : lat_heat_vap;

      LE = LE1*f1 + LE2*f2 + LE3*f3;
      LET = LET1*f1 + LET2*f2 + LET3*f3;
      H = H1*f1 + H2*f2 + H3*f3;
      // ---- SAVE THE FINAL VALUES FOR THE CANOPY BALANCE
      _species[s]._Temp_c->matrix[r][c] = x2;
      netR_a = NetRadCanopy(atm, x2, emissivity, albedo, BeerK, LAI, r, c);
      CalculateCanopyConduct(bas, atm, ctrl, x1, dgcdfgspsi, s, r, c);       //Updates canopy conductance with final values of soil potential
    } else {

      //provide initial guess  for loop
      x[0] = Sold;
      x[1] = psiae * 0.0098 / powl(x[0], bclambda); // 0.0098 convert from m to MPa
      x[2] = airTp;

      //used to calculate the gc factors other than f_lwp
      // this information is contained in dgcdfgspsi and is used to calculate gc in the solution scheme below
      CalculateCanopyConduct(bas, atm, ctrl, x[1], dgcdfgspsi, s, r, c);

      do {

	lambda = x[2] < 0 ? lat_heat_vap + lat_heat_fus : lat_heat_vap;

	gc = dgcdfgspsi*std::max<REAL8>(0,std::min<REAL8>(1,(lwp_min - x[1]) /
							  (lwp_min - lwp_max)));
	//gc = dgcdfgspsi * 1 / (1 + powl(x[1]/lwp_den, lwp_c));

	if (gc < 1e-13)
	  gc = 1e-13;

	ra_t = ra + (1 / gc);

	if(x[0]<0)
	  x[0] = 0.01;
	temp = -x[1] * rho_w * grav * Vw / (Rbar*(x[2]+273.15));
	if (temp >-708.4)
	  leafRH = std::min<REAL8>(1,expl(temp));
	else
	  leafRH = 0;

	//leafRH = 1;
	if(abs(evap_a)<RNDOFFERR){
	  LE = 0.0;
	} else {
	  LEmax = evap_a * (rho_w * lambda); // maximum latent heat based on water available
	  LE = LatHeatCanopy(bas, atm, leavesurfRH, ra, x[2], r, c) < LEmax ? LEmax : LatHeatCanopy(bas, atm, leavesurfRH, ra, x[2], r, c);
	}
	LET = LatHeatCanopy(bas, atm, leafRH, ra_t, x[2], r, c);
	H = SensHeatCanopy(atm, ra, x[2], r, c);

	E = -LET / (rho_w * lambda);
	E= std::max<REAL8>(0.0,E);


	F[0] = (x[0] - Sold) * maxAv * rootdepth / dt + E;
	F[1] = psiae * 0.0098 / powl(x[0], bclambda) - x[1];
	F[2] = NetRadCanopy(atm, x[2], emissivity, albedo, BeerK, LAI, r, c)
	  + LE + H + LET;


	dleafRHdT = leafRH *  x[1] * rho_w * grav * Vw / (Rbar*(x[2]+273.15)*(x[2]+273.15));

	if(leafRH == 1)
	  dleafRHdpsi_l = 0;
	else
	  dleafRHdpsi_l = - rho_w * grav * Vw * leafRH  / (Rbar*(x[2]+273.15));

	es = SatVaporPressure(x[2]);
	desdTs = es * ((17.3 / (x[2] + 237.3)) - 17.3 * x[2] / (powl(x[2] + 237.3, 2)));

	//dgcdlwp = gc == 1e-13 ? 0 : - dgcdfgspsi * lwp_c * powl(x[1]/lwp_den, lwp_c) / (x[1] * ( powl(x[1]/lwp_den, lwp_c) + 1) * ( powl(x[1]/lwp_den, lwp_c) + 1));
	if (gc < 1e-12 || x[1] > lwp_min)
	  dLETdlwp = LET = E = 0;
	else if (x[1] < lwp_max)
	  dLETdlwp = -spec_heat_air * rho_a * (ea - es) * (lwp_min - lwp_max)
	    / (gamma * dgcdfgspsi * (lwp_min - lwp_max) * (lwp_min - lwp_max) * ra_t * ra_t);
	else
	  dLETdlwp = -spec_heat_air * rho_a * (ea - es) * (lwp_min - lwp_max)
	    / (gamma * dgcdfgspsi * (lwp_min - x[1])  * (lwp_min - x[1]) * ra_t * ra_t);

	dLETdT = - rho_a * spec_heat_air / (ra_t * gamma) * (desdTs*leafRH + es*dleafRHdT);

	dEdlwp = - dLETdlwp / (rho_w * lambda);
	dEdT = - dLETdT / (rho_w * lambda);

	// Fill out the jacobian
	J(0,0) = rootdepth * maxAv / dt;
	J(0,1) = E==0 ?  0 : dEdlwp;
	J(0,2) = E==0 ? 0 : dEdT;

	J(1,0) = -bclambda * psiae * 0.0098 * powl(x[0], -(bclambda + 1));
	J(1,1) = -1;

	J(2,1) = dLETdlwp;
	J(2,2) = fA * powl(x[2] + 273.2, 3) + fB * desdTs * leavesurfRH + fC + dLETdT;

	// solve system
	if (!solve(deltax, J, -F)) {
	  cout << "Singular Jacobian found in Newton solver for canopy balance.\n";
	}

	x += deltax;

	k++;

      } while (norm(F, "inf") > 0.000001 && k < MAX_ITER);

      if (k >= MAX_ITER)
	std::cout
	  << "WARNING: non-convergence in canopy energy balance at cell row: "
	  << r << " col: " << c << " closure err: " << norm(deltax, 2)
	  << endl;

      if (x[2] < atm.getTemperature()->matrix[r][c]) { //if the calculated canopy temperature is lower than air temperature make it air temperature
	LET = LatHeatCanopy(bas, atm, leafRH, ra_t, x[2], r, c);
	H = SensHeatCanopy(atm, ra, x[2], r, c);
	netR_a = NetRadCanopy(atm, x[2], emissivity,	albedo, BeerK, LAI, r, c);
	if(abs(evap_a)<RNDOFFERR){
	  LE = 0.0;
	} else {
	  LEmax = evap_a * (rho_w * lambda); 
	  LE = LatHeatCanopy(bas, atm, leavesurfRH, ra, x[2], r, c) < LEmax ? LEmax : LatHeatCanopy(bas, atm, leavesurfRH, ra, x[2], r, c);

	} //if
	x[2] = atm.getTemperature()->matrix[r][c];
      } else {
	netR_a = NetRadCanopy(atm, x[2], emissivity,	albedo, BeerK, LAI, r, c);
	H = SensHeatCanopy(atm, ra, x[2], r, c);
      }
      // ---- SAVE THE FINAL VALUES FOR THE CANOPY BALANCE
      _species[s]._Temp_c->matrix[r][c] = x[2];
      CalculateCanopyConduct(bas, atm, ctrl, x[1], dgcdfgspsi, s, r, c);       //Updates canopy conductance with final values of soil potential

    }

    evap_a = std::max<REAL8>(0.0,std::min<REAL8>(-evap_a, fabs(-LE / (rho_w * lambda)))); //swap sign since outgoing evaporation is negative and we accumulate it as positive. Also checks for negative evap
    transp_a = std::max<REAL8>(0.0, -LET / (rho_w * lambda)); //swap sign since outgoing evaporation is negative and we accumulate it as positive. Also checks for negative evap

    //////////////////////////////////////////
    // This chunk of code is to make sure we are not transpiring below residual moisture content
    // Probably not needed anymore since mass balance is enforced in the system of eqs.
    // solved in this function
    REAL8 Tp = transp_a * ctrl.dt;
    REAL8 WaterAv = f1* (std::max<REAL8>(0.0,(theta1 - theta_wp))) + 
                    f2* (std::max<REAL8>(0.0,(theta2 - theta_wp))) + 
                    f3* (std::max<REAL8>(0.0,(theta3 - theta_wp)));
        
    ///TODO: change to wilting point (not residual water content)
    if (WaterAv * rootdepth < Tp) {
      Tp = WaterAv * rootdepth;
      transp_a = Tp / ctrl.dt;
    }

    // There is a need to updated the latent heat based on limiations on evap
    if (evap_a > RNDOFFERR){
      LE = -evap_a * (rho_w * lambda);
    }

    // There is a need to updated the latent heat based on limiations of transpiration
    if (transp_a > RNDOFFERR){
      LET = -transp_a *(rho_w * lambda) ;
    }

    DelCanStor -= evap_a * ctrl.dt;

    le_a = LE + LET;

    _species[s]._NetR_Can->matrix[r][c] = netR_a ; //Net radiation
    _species[s]._LatHeat_CanE->matrix[r][c] = LE ; // Latent heat of canopy evap
    _species[s]._LatHeat_CanT->matrix[r][c] = LET; // Latent heat of transpiration
    _species[s]._SensHeat_Can->matrix[r][c] = H; // Sensible heat


    /* THIS PREVIOUS UPDATE CAN CAUSE PROBLEMS, NOT USED FOR NOW
    // Here we also check that theta is larger than both wp and thetar, otherwise
    // Tp could end up being negative if e.g. wp > thetar
    if ((theta - std::max<double>(_species[s].WiltingPoint,thetar)) * rootdepth < Tp) {
      theta = (theta - std::max<double>(_species[s].WiltingPoint,thetar)) < RNDOFFERR ? 0 :
	(theta - std::max<double>(_species[s].WiltingPoint,thetar));
      Tp = theta * rootdepth;
      transp_a = Tp / ctrl.dt;
    }
    */
    //////////////////////////////////////////

    _species[s]._Einterception->matrix[r][c] = evap_a;
    _species[s]._Transpiration->matrix[r][c] = transp_a;
    _species[s]._ET->matrix[r][c] = transp_a + evap_a; // then Es is added in surface routines

    _species[s]._WaterStorage->matrix[r][c] -= evap_a * ctrl.dt;

  }

  return EXIT_SUCCESS;
}

