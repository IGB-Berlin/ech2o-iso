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
 *    Aaron Smith
 *******************************************************************************/
/*
 * SolveChannelEnergyBalance.cpp
 *
 *  Created on: Feb 2, 2020
 *      Author: Aaron Smith
 */

#include "Basin.h"

int Basin::SolveChannelEnergyBalance(Atmosphere &atm,
				     Control &ctrl,
				     Tracking &trck,
				     REAL8 chan_prec,
				     REAL8 qout,
				     UINT4 r,
				     UINT4 c){
  
  float dt = ctrl.dt; //time step
  REAL8 fA, fB, fC, fD; //pooling factors
  //Temperature inputs
  REAL8 Ts = _Temp_s->matrix[r][c]; // temperature of the surface
  REAL8 Td = _Temp_d->matrix[r][c]; // temperature of the groundwater
  REAL8 Tw = _Temp_w->matrix[r][c]; // temperature of the water in the channel
  REAL8 Tw_old = _Temp_w->matrix[r][c]; // temperature of the water in the channel
  REAL8 Tw_init;
  REAL8 Tw1;
  REAL8 Ta = atm.getTemperature()->matrix[r][c];
  REAL8 fTw;
  REAL8 dfTw;
  REAL8 desdTw;
  //Energy balance compoents
  REAL8 z;
  REAL8 gamma;
  REAL8 LE = 0;
  REAL8 H = 0;
  REAL8 Qx = 0;
  REAL8 beta = 1; // adjustment of soil relative humidity to account for pores (hs=beta+(1-beta)*ha)
  REAL8 rho_a; //density of air
  REAL8 lambda = Tw < 0 ?  lat_heat_vap + lat_heat_fus : lat_heat_vap;
  REAL8 ra_H, ra_LE;    // aero resistance
  REAL8 wind = atm.getWindSpeed()->matrix[r][c];
  REAL8 chan_rough = powl(10,_chan_roughness->matrix[r][c]);
  REAL8 rough_LE = min<REAL8>(_random_roughness->matrix[r][c],_channelwidth->matrix[r][c]/chan_rough);
  REAL8 za_LE = rough_LE + 2 ;
  REAL8 z0u_LE = rough_LE * 0.1;
  REAL8 zdu_LE = rough_LE * 0.7;

  REAL8 rough_H = _random_roughness->matrix[r][c];
  REAL8 za_H = rough_H + 2 ;
  REAL8 z0u_H = max<REAL8>(0.000005,rough_H * 0.1);
  REAL8 zdu_H = rough_H * 0.7;
  REAL8 evap = 0; // initialize the channel evaporation

  //For the temperature mixing
  REAL8 GWtoChn = _FluxGWtoChn->matrix[r][c];
  REAL8 LattoChn = _FluxLattoChn->matrix[r][c];
  REAL8 T_LattoChn = _FTemp_w->matrix[r][c];
  REAL8 L1toSrf = _FluxExfilt->matrix[r][c];
  REAL8 chan_store = _chan_store->matrix[r][c];
  REAL8 chan_store_old = _chan_store_old->matrix[r][c];
  UINT4 mixmod = 0;//ctrl.toggle_Tmix;
  REAL8 qin = GWtoChn + LattoChn + L1toSrf + chan_prec;
  REAL8 Tin = qin < RNDOFFERR ? Ts : (GWtoChn*Td + T_LattoChn + L1toSrf * Ts + chan_prec * Ta) / qin;
  REAL8 rs = 0; //no extra resistance to evaporation
  // Canopy influence
  REAL8 rc = 0; //extra resistance due to overhanging canopy
  UINT4 nsp = fForest->getNumSpecies();
  UINT4 s;
  REAL8 za = 0;
  REAL8 p;
  REAL8 z0o;
  REAL8 zdo;
  REAL8 treeheight;

  z = _DEM->matrix[r][c];

  gamma =PsychrometricConst(101325, z);
 
  rho_a = AirDensity(Ta); //kgm-3
  
  ra_LE = CalcAerodynResist(wind , za_LE , z0u_LE , zdu_LE , 0 , 0 , 0 , 0 ,Tw,Ta,ctrl.toggle_ra,true);
  //                     wind , za , z0u , zdu ,z0u,zdo,hgt,LAI,Tw,Ta, toggle, true

  ra_H  = CalcAerodynResist(wind , za_H , z0u_H , zdu_H , 0 , 0 , 0 , 0 ,Tw,Ta,ctrl.toggle_ra,true);
  //                     wind , za , z0u , zdu ,z0u,zdo,hgt,LAI,Tw,Ta, toggle, true

  // Resistance due to overhanging vegetation
  for (s=0; s < nsp; s++) {
    p = fForest->getPropSpecies(s,r,c);
    if (p == 0)
      continue;

    if (s == nsp - 1){
      rc += 0 * p;
    } else {
      treeheight = max<REAL8>(0.01, fForest->getTreeHeight(s, r, c));
      za = treeheight + 2;
      z0o = powl(treeheight, 1.19) * 0.057544;
      zdo = powl(treeheight, 0.98) * 0.707946;
      rc += CalcAerodynResist(wind, za, 0, 0, z0o, zdo, treeheight,
				 fForest->getLAISpecies(s, r, c),
				 getCanopyTemp(s)->matrix[r][c],
				 atm.getTemperature()->matrix[r][c], ctrl.toggle_ra,
				 false) * p;
    }
  }

  //energy balance factors that do not need updating in the N-R loop
  fA = -4* 0.995 * stefboltz;	        //pools together net radiation factors ( emissivity of water is 0.995)
  fB = (-1/gamma) * (1/(ra_LE + rs + rc)) * rho_a * spec_heat_air; // pools together the latent heat factors
  fC = (-1/(ra_H)) * rho_a * spec_heat_air; // pools together the sensible heat factors

  fD = (-spec_heat_water * rho_w / ( dt * _channelwidth->matrix[r][c] ) * _channelwidth->matrix[r][c] * chan_store);
  //  fD = (-spec_heat_water * rho_w / ( _dx * _channelwidth->matrix[r][c] ) * ( qout ));
  
  Tw_init = InOutMix_Temperature(chan_store_old,Tw_old,qin,Tin,qout,mixmod);

  int k = 0;
  
  do{

    Tw = _Temp_w->matrix[r][c];

    desdTw = 611 * ( (17.3/( Tw + 237.7)) - ((17.3 * Tw)/(powl(Tw + 237.2 , 2))) ) * exp(17.3 * Tw /( Tw + 237.7));

    //This is the water advected energy accounting needed for input/output of evaporation 
    //    Qv = spec_heat_water * rho_w / ( _dx * _channelwidth->matrix[r][c] ) * ( qin * Tin - qout * Tw);

    //This is the latent heat of evaporation
    LE = LatHeat(atm, beta, ra_LE, rc, 0.0, Tw, r, c); // latent heat flux

    if (LE < 0){
      evap = std::min<REAL8>(LE/(rho_w*lambda) * dt,chan_store);
    } else {
      evap = 0;
    }
    //This is the water stored energy
    Qx = spec_heat_water * rho_w / (dt * _channelwidth->matrix[r][c]) * (_channelwidth->matrix[r][c] * (Tw_init * chan_store - Tw * (chan_store - evap < RNDOFFERR ? 0.0 : chan_store - evap)) );

    //This is the sensible heat
    H = SensHeat(atm, ra_H, Tw, r, c); // sensible heat flux

    fTw = NetRad_water(atm, Tw,  r, c) + LE + H - Qx;  // total energy balance

    dfTw = fA*powl(Tw + 273.2, 3) + fB * desdTw + fC - fD;

    Tw1 = Tw - (fTw/dfTw);
    
    _Temp_w->matrix[r][c] = Tw1;

    k++;

  }while(fabs(Tw1 - Tw) > 0.00001 && k < MAX_ITER);
  
  if (k>=MAX_ITER){
    std::cout << "WARNING: non-convergence in channel energy balance at cell row: " << r << " col: " << c << " closure err: " << (Tw1 - Tw) << endl;
  }

  if (Tw1 > 35){ //correct the water temperature to a max/min value
    Tw1 = Tw_init;
    _Temp_w->matrix[r][c] = Tw_init;
    LE = LatHeat(atm, beta, ra_LE, rc, 0.0, Tw1, r, c); // latent heat flux
  }

  if (Tw1 < 0){ //correct the water temperature to a max/min value
    Tw1 = 0;
    _Temp_w->matrix[r][c] = 0.0;
    LE = LatHeat(atm, beta, ra_LE, rc, 0.0, Tw1, r, c); // latent heat flux
  }

  // If there is no evaporation, channel water temperature is a mixture of what is input and what is stored
  //  if (LE > 0){
  //  Tw = InOutMix_Temperature(chan_store_old,Tw_old,qin,Tin,qout,mixmod);
  //  _Temp_w->matrix[r][c] = Tw;
  //}

  evap = 0;

  ChannelEvaporation(-LE,lambda,ra_LE,evap, chan_store,dt,r,c);

   // Flux tracking after evap
  if(ctrl.sw_trck){
    // Convert Tw to Kelvins for fractionation  
    trck.MixingV_evapW(atm, *this, ctrl, _Temp_w->matrix[r][c]+273.15, evap, chan_store, r, c);   
  }
  return EXIT_SUCCESS;
}
