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
				       REAL8 &transp_a, REAL8 &netR_a, REAL8 &le_a,REAL8 &h_a,
				       UINT4 s, UINT4 r, UINT4 c) {

  //energy balance parameters
  REAL8 fA, fB, fC; //pooling factors
  REAL8 rho_a; //density of air
  REAL8 airRH; //air humidity
  REAL8 airTp; // air temperature
  REAL8 airPr; // air pressure
  REAL8 emissivity; //canopy emissivity
  REAL8 albedo; //canopy albedo
  REAL8 LAI;
  REAL8 BeerK; //Beers-Lambert coefficient
  REAL8 LE, H;
  REAL8 LE1, LE2, LE3;
  REAL8 LET1, LET2, LET3;
  REAL8 H1, H2, H3;
  REAL8 z; //terrain height
  REAL8 gamma; //psychrometric constant
  REAL8 f1, f2, f3; //root fractions

  REAL8 lambda = lat_heat_vap;
  REAL8 LET; //latent heat of transpiration
  REAL8 CanStor = 0;
  REAL8 MaxCanStor = 0;

  REAL8 leafRH; //soil relative humidty use in teh calculation of soil vapor pressure for latent heat exchanges
  REAL8 leavesurfRH; //relative humidity of the leave surface. 1 when leave is saturated with intercepted water, airRH when no water

  // variables for Sperry's model
  REAL8 Sold = 0;  // Soil Saturation at beginning of t
  REAL8 maxAv = 0; // Equivalent (weighted) maximum available water 
  REAL8 lwp_min, lwp_max; // lower and higher boundaries for linear lwp stomatal model
  REAL8 dgcdfgspsi = 0;

  // Sapflow estimation
  REAL8 stem_d; //diameter of stem
  REAL8 xylem_d;//diamter of xylem
  REAL8 sapArea;//area of sapflow [m2]
  REAL8 sap_vel = 0;

  // harmonic mean of soil conductivity
  REAL8 Keff1 = 2 * bas.getKsatL1()->matrix[r][c]*bas.getKvKsL1()->matrix[r][c] /
    (1 + bas.getKvKsL1()->matrix[r][c]);
  REAL8 Keff2 = 2 * bas.getKsatL2()->matrix[r][c]*bas.getKvKsL2()->matrix[r][c] /
    (1 + bas.getKvKsL2()->matrix[r][c]);
  REAL8 Keff3 = 2 * bas.getKsatL3()->matrix[r][c]*bas.getKvKsL3()->matrix[r][c] /
    (1 + bas.getKvKsL3()->matrix[r][c]);  
  REAL8 Keff = 0;

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
  REAL8 x00,x01,x02,x10,x11,x12;
  REAL8 x20,x21,x22,x30,x31,x32;
  /* ----------------------------------------- */

  REAL8 theta_wp = 0; //wilting point
  REAL8 theta1 = bas.getSoilMoist1()->matrix[r][c];
  REAL8 theta2 = bas.getSoilMoist2()->matrix[r][c];
  REAL8 theta3 = bas.getSoilMoist3()->matrix[r][c];
  REAL8 poros1 = bas.getPorosityL1()->matrix[r][c];
  REAL8 poros2 = bas.getPorosityL2()->matrix[r][c];
  REAL8 poros3 = bas.getPorosityL3()->matrix[r][c];

  UINT4 nsp = getNumSpecies();

  if (s == nsp - 1) //for bare soil, water reaching the ground is pp times its proportion of the cell
    evap_a = transp_a = netR_a = le_a = 0;
  else {
    airTp = atm.getTemperature()->matrix[r][c];
    airPr = atm.getPressure()->matrix[r][c];
    rho_a = AirDensity(airTp,airPr); //kgm-3
    z = bas.getDEM()->matrix[r][c];
    gamma = PsychrometricConst(airPr, z);
    airRH = atm.getRelativeHumidty()->matrix[r][c];

    albedo = _species[s].albedo;
    emissivity = _species[s].emissivity;
    BeerK = _species[s].KBeers;
    LAI = _species[s]._LAI->matrix[r][c];
    theta_wp = _species[s].WiltingPoint;

    stem_d = powl(_species[s]._BasalArea->matrix[r][c]/PI,0.5) * 2 * 100; //diameter in cm
    //This is estimated for self-supporting angiosperms Olson&Rossel 2012
    xylem_d = powl(10,log10(stem_d) * 0.36 + 1.485) * 1e-6;     //xylem diameter in the stem in meters
    sapArea = _species[s]._BasalArea->matrix[r][c] - (PI * powl((stem_d/100 - 2 * xylem_d) /2 , 2));

    // Linear model has automatic "cutoff" at wilting point
    if(ctrl.toggle_sm == 0){ // nonlinear model //TODO CLEAN UP IN CONSTANDFUNCS
      lwp_min = _species[s].lwp_min;
      lwp_max = _species[s].lwp_max;
    } else {
      lwp_min = psiae * 0.0098 / powl(Saturation(theta_wp,theta_r1,poros1),bclambda) * _species[s].lwp_min;
      lwp_max = psiae * 0.0098 / powl(Saturation(theta_wp,theta_r1,poros1),bclambda) * _species[s].lwp_max;
    }

    // Root fractiona
    f1 = _species[s]._rootfrac1->matrix[r][c];
    f2 = _species[s]._rootfrac2->matrix[r][c];
    f3 = 1 - f1 - f2;

    CanStor = getIntercWater(s, r, c);
    MaxCanStor = getMaxCanopyStorage(s, r, c);

    // Calculate Humidity on the leaf and soil pore relative humidity
    leafRH = 1; 
    leavesurfRH = airRH + ((1 - airRH) / MaxCanStor) * CanStor;

    fA =  -4 * emissivity * stefboltz;//pools together net radiation factors
    fB = (-1 / (ra * gamma)) * rho_a * spec_heat_air; // pools together the latent heat factors
    fC = (-1 / (ra)) * rho_a * spec_heat_air; // pools together the sensible heat factors

    evap_a = CanStor < MaxCanStor ? -CanStor / ctrl.dt * powl((CanStor / MaxCanStor), 0.6) : -CanStor / ctrl.dt;
    //********************************************************************************************
    // Calculation of conductances necessary for the implementation of Sperry's model to calculate leaf water potential.
    // Soil to root conductance. 
    // Adapted from Rodriguez-Iturbe and Porporato (eq 6.4, page 181) for units of hydraulic head
    // Here using a weighted average because porosity changes with depth
    // EDIT: this formulation makes things much more simple
    maxAv = f1*(poros1 - theta_r1) + f2*(poros2-theta_r2)+f3*(poros3-theta_r3);
    Keff = 1 / (f1 / Keff1 + f2 / Keff2 + f3 / Keff3);
    Sold = (f1*std::max(0.0,(theta1-theta_r1))+f2*std::max(0.0,(theta2-theta_r2))+f3*std::max(0.0,(theta3-theta_r3))) / maxAv ;

    /***
     * state variables:
     * x[0]: S - (degree of saturation at time t+1)
     * x[1]: psi_s - soil water potential
     * x[2]: psi_l - leaf water potential
     * x[3]: Ts - Leaf temperature
     ***/

    colvec x(4);
    colvec deltax(4);
    colvec F(4);
    mat J = zeros<mat>(4,4);

    if(ctrl.toggle_soil_prop > 0){ //Each soil layer has different soil properties
      Sold1 = Saturation(theta1,theta_r1,poros1);
      Sold2 = Saturation(theta2,theta_r2,poros2);
      Sold3 = Saturation(theta3,theta_r3,poros3);

      SperryModel(bas, atm, ctrl,d1,Sold1,Keff1,psiae1,bclambda1,airTp,airRH,rho_a,gamma,ra,poros1,theta_r1,theta_wp,evap_a,fA,fB,fC,
		  leavesurfRH,leafRH,LET1,LE1,H1,x00,x10,x20,x30,s,r,c);
      SperryModel(bas, atm, ctrl,d2,Sold2,Keff2,psiae2,bclambda2,airTp,airRH,rho_a,gamma,ra,poros2,theta_r2,theta_wp,evap_a,fA,fB,fC,
		  leavesurfRH,leafRH,LET2,LE2,H2,x01,x11,x21,x31,s,r,c);
      SperryModel(bas, atm, ctrl,d3,Sold3,Keff3,psiae3,bclambda3,airTp,airRH,rho_a,gamma,ra,poros3,theta_r3,theta_wp,evap_a,fA,fB,fC,
		  leavesurfRH,leafRH,LET3,LE3,H3,x02,x12,x22,x32,s,r,c);

      //Combine all estimated fluxes from each soil layer -- update vector [x]
      x[0] = f1*x00 + f2*x01 + f3*x02;
      x[1] = f1*x10 + f2*x11 + f3*x12;
      x[2] = f1*x20 + f2*x21 + f3*x22;
      x[3] = f1*x30 + f2*x31 + f3*x32;
      netR_a = NetRadCanopy(atm, x[3], emissivity, albedo, BeerK, LAI, r, c);
      if (x[3] < atm.getTemperature()->matrix[r][c])
	x[3] = atm.getTemperature()->matrix[r][c];
      _species[s]._Temp_c->matrix[r][c] = x[3];
      lambda = x[3] < 0 ? lat_heat_vap + lat_heat_fus : lat_heat_vap;

      LE = LE1*f1 + LE2*f2 + LE3*f3;
      LET = LET1*f1 + LET2*f2 + LET3*f3;
      H = H1*f1 + H2*f2 + H3*f3;
      CalculateCanopyConduct(bas, atm, ctrl, x[1], lwp_max,lwp_min,dgcdfgspsi, s, r, c);       //Updates canopy conductance with final values of soil potential
    } else {
      SperryModel(bas, atm, ctrl,d1,Sold,Keff,psiae,bclambda,airTp,airRH,rho_a,gamma,ra,
		  poros1,theta_r1,theta_wp,evap_a,fA,fB,fC,
		  leavesurfRH,leafRH,LET,LE,H,x00,x10,x20,x30,s,r,c);
      x[0] = x00;
      x[1] = x10;
      x[2] = x20;
      x[3] = x30;
      lambda = x[3] < 0 ? lat_heat_vap + lat_heat_fus : lat_heat_vap;      
      netR_a = NetRadCanopy(atm, x[3], emissivity, albedo, BeerK, LAI, r, c);
      if (x[3] < atm.getTemperature()->matrix[r][c])
	x[3] = atm.getTemperature()->matrix[r][c];
      _species[s]._Temp_c->matrix[r][c] = x[3];
      CalculateCanopyConduct(bas, atm, ctrl, x[1], lwp_max,lwp_min,dgcdfgspsi, s, r, c);
      //Updates canopy conductance with final values of soil potential
    }

    // Evap and Transp - swap sign since outgoing evaporation is negative and we accumulate it as positive. Also checks for negative evap
    evap_a = std::max<REAL8>(0.0,-LE / (rho_w * lambda)); 
    transp_a = std::max<REAL8>(0.0, -LET / (rho_w * lambda)); 
    
    //////////////////////////////////////////
    // This chunk of code is to make sure we are not transpiring below wilting point
    // Probably not needed anymore since mass balance is enforced in the system of eqs.
    // solved in this function
    REAL8 Tp = transp_a * ctrl.dt;
    REAL8 WaterAv = f1* (std::max<REAL8>(0.0,(theta1 - theta_wp))) + 
                    f2* (std::max<REAL8>(0.0,(theta2 - theta_wp))) + 
                    f3* (std::max<REAL8>(0.0,(theta3 - theta_wp)));
        
    if (WaterAv * rootdepth < Tp) {
      Tp = WaterAv * rootdepth;
      transp_a = Tp / ctrl.dt;
    }

    LE  = (evap_a > RNDOFFERR)   ? (-evap_a * (rho_w * lambda))  : 0.0 ;
    LET = (transp_a > RNDOFFERR) ? (-transp_a *(rho_w * lambda)) : 0.0 ;

    DelCanStor -= evap_a * ctrl.dt;

    le_a = LE + LET;
    h_a = H;

    sap_vel = transp_a * powl(_dx,2) *_species[s]._StemDensity->matrix[r][c] / (_species[s]._LAI->matrix[r][c] * sapArea);

    _species[s]._RUptakeL1->matrix[r][c] = (std::max(0.0,(theta1-theta_wp))*f1) / WaterAv;
    _species[s]._RUptakeL2->matrix[r][c] = (std::max(0.0,(theta2-theta_wp))*f2) / WaterAv;
    _species[s]._RUptakeL3->matrix[r][c] = (std::max(0.0,(theta3-theta_wp))*f3) / WaterAv;

    _species[s]._NetR_Can->matrix[r][c] = netR_a ;    //Net radiation
    _species[s]._LatHeat_CanE->matrix[r][c] = LE ;    // Latent heat of canopy evap
    _species[s]._LatHeat_CanT->matrix[r][c] = LET;    // Latent heat of transpiration
    _species[s]._SensHeat_Can->matrix[r][c] = H;      // Sensible heat
    _species[s]._SoilWatPot->matrix[r][c] = -x[1];    //Soil water potential
    _species[s]._LeafWatPot->matrix[r][c] = -x[2];    //Leaf water potential
    _species[s]._SapVelocity->matrix[r][c] = sap_vel; //sap velocity [m/s]
    _species[s]._Einterception->matrix[r][c] = evap_a;
    _species[s]._Transpiration->matrix[r][c] = transp_a;
    _species[s]._ET->matrix[r][c] = transp_a + evap_a; // then Es is added in surface routines
    _species[s]._WaterStorage->matrix[r][c] -= evap_a * ctrl.dt;

  }

  return EXIT_SUCCESS;
}

