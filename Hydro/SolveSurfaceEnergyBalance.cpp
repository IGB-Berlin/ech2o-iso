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
 * SolveSurfaceEnergyBalance.cpp
 *
 *  Created on: Jul 10, 2010
 *      Author: Marco.Maneta
 */

#include "Basin.h"

int Basin::SolveSurfaceEnergyBalance(Atmosphere &atm,
				     Control &ctrl,
				     Tracking &trck,
				     REAL8 ra,
				     REAL8 rs,
				     REAL8 rc,
				     REAL8 Kbeers,
				     REAL8 lai,
				     REAL8 emis_can,
				     REAL8 Temp_can,
				     REAL8 &nrad,
				     REAL8 &latheat,
				     REAL8 &sensheat,
				     REAL8 &grndheat,
				     REAL8 &snowheat,
				     REAL8 &meltheat,
				     REAL8 &Tsold,
				     REAL8 &etp,
				     REAL8 &pond,
				     REAL8 &theta,
				     REAL8 &Ts1,
				     REAL8 &Tdold,
				     REAL8 &TdL1,
				     REAL8 &TdL2,
				     REAL8 &TdL3,
				     REAL8 p,
				     UINT4 bsoil,
				     UINT4 r,
				     UINT4 c,
				     UINT4 s){
  
  float dt = ctrl.dt; //time step
  REAL8 fA, fB, fC, fD, fG, fH, fHa; //pooling factors
  REAL8 C; // soil heat capacity
  REAL8 K; // soil thermal heat conductivity
  REAL8 Pe = dt < 86400 ? 86400*14 : 31536000; //period is daily if time step is less than a day adn yearly if time step is daily or larger
  REAL8 Pres = atm.getPressure()->matrix[r][c];
  REAL8 Ts = _Temp_s->matrix[r][c];
  REAL8 Td = _Temp_d->matrix[r][c];
  REAL8 Td1= _Temp_L1->matrix[r][c];
  REAL8 Td2= _Temp_L2->matrix[r][c];
  REAL8 Td3= _Temp_L3->matrix[r][c];
  REAL8 fTs;
  REAL8 dfTs;
  REAL8 desdTs;
  REAL8 d; // temperature fluctuation damping depth
  REAL8 d0; //bottom depth of bottom thermal layer
  REAL8 z;
  REAL8 gamma;
  REAL8 LE = 0;
  REAL8 H = 0;
  REAL8 G = 0;
  REAL8 S = 0;
  REAL8 LM = 0;
  REAL8 R = 0; // the last two variables are the latent heat of melt and the heat advected by rain
  REAL8 RainIntensity;
  REAL8 W; //anthropogenic heat
  REAL8 LEp = 0; // the latent heat fraction of vegetation proportion
  REAL8 MeltFac; //snowmelt factor
  REAL8 h; //snow water equivalent
  REAL8 n; //porosity
  REAL8 fc; //field capacity 
  //REAL8 thetar; //residual moisture content
  REAL8 beta; // adjustment of soil relative humidity to account for pores (hs=beta+(1-beta)*ha)
  //REAL8 ea; //emissivity of air
  REAL8 nrad_a; // net radiation
  REAL8 rho_a; //density of air
  //REAL8 RainIntensity; //ms-1
  //REAL8 exfilt;
  REAL8 lambda = Ts1 < 0 ?  lat_heat_vap + lat_heat_fus : lat_heat_vap;
  REAL8 d1; // depth of the first layer
  REAL8 d2;
  REAL8 d3;
  REAL8 fImperv = _fImperv->matrix[r][c]; //impervious fraction
  double d2H_evap, d18O_evap, Age_evap;

  z = _DEM->matrix[r][c];
  gamma =PsychrometricConst(Pres, z);
  d0 = _dampdepth->matrix[r][c];

  d1 = _depth_layer1->matrix[r][c];
  d2 = d1 + _depth_layer2->matrix[r][c];
  d3 = _soildepth->matrix[r][c] - d1 - d2;

  MeltFac = _meltCoeff->matrix[r][c];

  h = _snow->matrix[r][c];
  n = _porosityL1->matrix[r][c];
  rho_a = AirDensity(atm.getTemperature()->matrix[r][c],Pres); //kgm-3
  
  C = SoilHeatCapacity(_soil_dry_heatcap->matrix[r][c],  n, theta, Ts1,Pres);
  K = SoilHeatConductivity(_soil_dry_thermcond->matrix[r][c], n,theta);
  
  d = sqrt( (K/C) / ( 2 * ( 2 * PI / Pe) ) );
  
  fc = _fieldcapL1->matrix[r][c];
  if (h>0.005){
    beta = 1; //RH in snow pores is assumed to be saturated. Switch when there is at least 1 cms of snow
    rs = 0; //no extra resistance to evaporation
  }  else {
    beta = Calculate_beta(theta, fc);
  }
  
  //energy balance factors that do not need updating in the N-R loop
  fA = -4*_emiss_surf->matrix[r][c] * stefboltz;	//pools together net radiation factors
  fB = (-1/gamma) * (1/(ra + rs) + rc) * rho_a * spec_heat_air; //(-1/((ra + rs + rc) * gamma)) * rho_a * spec_heat_air; // pools together the latent heat factors
  fC = (-1/(ra)) * rho_a * spec_heat_air; // pools together the sensible heat factors
  fD =  -(d * C /(2*dt)) - PI * d * C / Pe ; //same for ground heat flux (both terms). Assumes C does not depend on Ts (tiny dependency does nto affect derivative)
  fG = -spec_heat_ice * rho_w * h * (1 / dt); //and heat fluxes into the snow
  fH = -1*lat_heat_fus * rho_w * MeltFac; // last value is M factor
  
  RainIntensity = pond * p / dt;
  R = RainHeat(atm, RainIntensity, r, c); //heat advected by rain only if there is snowpack present

  W = 0;
  if(ctrl.sw_anthr_heat)
    W = atm.getAnthrHeat()->matrix[r][c];
  
  int k = 0;
  
  do{
    Ts = _Temp_s->matrix[r][c];
    Td = _Temp_d->matrix[r][c];
    desdTs = 611 * ( (17.3/( Ts + 237.7)) - ((17.3 * Ts)/(powl(Ts + 237.2 , 2))) )
      * exp(17.3 * Ts /( Ts + 237.7));
    
    Td = -( ((d/d0) * 2 * PI * (Td - Ts) / Pe) * dt ) + Td;
    
    if (h > 0.005){
      LE = fB = 0;
      G = 0;
    }
    else{ // is snowpack is thicker than 5 mm, insulate the soil and shutoff ground heat and evaporation
      G = GrndHeat(atm, ctrl, theta, Ts, Td, r, c);
      LE = LatHeat(atm, beta, ra, rs, rc, Ts, r, c);// * temp;
    }
    // for impervious surfaces reduce latent heat
    if(bsoil == 1 and fImperv > 0){ //no vegetation and impervious surfaces
      LE = LE * (1 - fImperv);
      G  = G * (1- fImperv);
    }
      
    H = SensHeat(atm, ra, Ts, r, c);
    S = SnowHeat(atm, ctrl, Ts, r, c);
    LM = MeltHeat(atm, ctrl, Ts, h, MeltFac, r, c);
    if( Ts < 0 || h==0 )
      fHa = 0;
    else
      fHa = fH;
    
    fTs = NetRad(atm, Ts, Kbeers, lai, emis_can, Temp_can,  0, r, c) + LE + H + G + S + LM + R + W;
    dfTs = fA*powl(Ts + 273.2, 3) + fB * desdTs * beta + fC + ((G==0)? 0 : fD) + fG + fHa;
    
    Ts1 = Ts - (fTs/dfTs);
    _Temp_s->matrix[r][c] = Ts1;
    k++;
  }while(fabs(Ts1 - Ts) > 0.00001 && k < MAX_ITER);
  
  if (k>=MAX_ITER)
    std::cout << "WARNING: non-convergence in surface energy balance at cell row: " << r << " col: " << c << " closure err: " << (Ts1 - Ts) << endl;

  nrad_a = NetRad(atm, Ts1, Kbeers, lai, emis_can, Temp_can, 0, r, c);
  
  if(s == fForest->getNumSpecies()-1)
    _Rn_sum->matrix[r][c] += nrad_a * p;
 
  LEp = -LE * p;
  SoilEvapotranspiration(LEp, Ts1, lambda, p*rs, etp, theta, dt, r, c);

  // Estimate the soil temperature at different depths
  Td1 = -( ((d/d1) * 2 * PI * (Td1 - Ts) / Pe) * dt ) + Td1;
  Td2 = -( ((d/d2) * 2 * PI * (Td2 - Ts) / Pe) * dt ) + Td2;
  Td3 = -( ((d/d3) * 2 * PI * (Td3 - Ts) / Pe) * dt ) + Td3;  

  // Fraction of energy by the vegetation proportion
  nrad += nrad_a * p;
  latheat += -LEp;
  sensheat += SensHeat(atm, ra, Ts1, r, c) * p;
  grndheat += GrndHeat(atm, ctrl, theta, Ts1, Td, r, c) * p;
  snowheat += SnowHeat(atm, ctrl, Ts1, r, c) * p;
  meltheat += MeltHeat(atm, ctrl, Ts1, h, MeltFac, r, c) * p;
  Tsold += Ts1 * p;
  Tdold += Td * p;
  TdL1 += Td1 * p;
  TdL2 += Td2 * p;
  TdL3 += Td3 * p;
  
  // Flux tracking after evap
  if(ctrl.sw_trck){

    // Convert Ts to Kelvins for fractionation  
    trck.MixingV_evapS(atm, *this, ctrl, d1, theta, Ts+273.15, etp, beta,
    		       d2H_evap, d18O_evap, Age_evap,r, c);   

    // Vegetation-dependent values are assigned here because fForest is private (to be fixed?)
    if(ctrl.sw_2H){
      fForest->setd2HevapS(s,r,c, d2H_evap);
    }
    if(ctrl.sw_18O){
      fForest->setd18OevapS(s,r,c, d18O_evap);
    }
    if(ctrl.sw_Age)
      fForest->setAgeevapS(s,r,c, Age_evap);

  }
  return EXIT_SUCCESS;
}
