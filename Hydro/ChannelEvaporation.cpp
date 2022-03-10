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
 * ChannelEvaporation.cpp
 *
 *  Created on: Feb 2, 2020
 *      Author: Aaron Smith
 */

#include "Basin.h"

int Basin::ChannelEvaporation(Atmosphere &atm,
				     Control &ctrl,
				     Tracking &trck,
				     REAL8 chan_prec,
				     REAL8 Qk1,
				     bool lat_ok,	
				     UINT4 r,
				     UINT4 c,
				     UINT4 rind,
				     UINT4 cind){
  
  float dt = ctrl.dt; //time step
  REAL8 fA, fB, fC, fD; //pooling factors
  //Temperature inputs
  REAL8 Ts = _Temp_s->matrix[r][c]; // temperature of the surface
  REAL8 Td = _Temp_d->matrix[r][c]; // temperature of the groundwater
  REAL8 Tw,Tw_old; // temperature of the water in the channel
  Tw = Tw_old  = 0;
  REAL8 Tw_init, Tw1, fTw, dfTw, desdTw;
  REAL8 Pa = atm.getPressure()->matrix[r][c];
  REAL8 Ta = atm.getTemperature()->matrix[r][c];
  //Energy balance compoents
  REAL8 z, gamma, rho_a, ra_H, ra_LE, es, ea;
  REAL8 evap,LE,H,Qx;
  REAL8 beta = 1; // adjustment of soil relative humidity to account for pores (hs=beta+(1-beta)*ha)
  REAL8 lambda = Tw < 0 ?  lat_heat_vap + lat_heat_fus : lat_heat_vap;
  REAL8 wind = atm.getWindSpeed()->matrix[r][c];
  REAL8 chan_rough,rough_LE,za_LE,z0u_LE,zdu_LE,rough_H;
  REAL8 za_H,z0u_H,zdu_H;

  //mass-transfer parameters
  REAL8 b0,b1;
  
  //For the temperature mixing
  REAL8 GWtoChn = _FluxGWtoChn->matrix[r][c];
  REAL8 LattoChn = _FluxLattoChn->matrix[r][c];
  REAL8 T_LattoChn =0;
  REAL8 L1toSrf = _FluxExfilt->matrix[r][c];
  REAL8 chan_store = _chan_store->matrix[r][c];
  UINT4 mixmod = 0;//ctrl.toggle_Tmix;
  REAL8 qin = GWtoChn + LattoChn + L1toSrf + chan_prec;
  REAL8 Tin = 0;
  REAL8 rs = 0; //no extra resistance to evaporation
  // Canopy influence
  REAL8 rc,za; //extra resistance due to overhanging canopy
  UINT4 nsp = fForest->getNumSpecies();
  UINT4 s;
  REAL8 p, z0o, zdo, treeheight;

  evap = LE = H = Qx = rc = za = 0; //Initialize
  z = _DEM->matrix[r][c];
  gamma =PsychrometricConst(Pa, z);
  rho_a = AirDensity(Ta,Pa); //kgm-3
  // ****************************************************************************
  // Estimate evaporation
  // ****************************************************************************
  if(ctrl.toggle_chan_evap == 0){ // no channel evaporation
    LE = 0; 							// LE is automatically 0 - no evaporation
    Tw1 = Ta; 							// only for tracer mixing
  }  else if(ctrl.toggle_chan_evap == 1){ // energy balance method
    T_LattoChn = _FTemp_w->matrix[r][c];
    Tin = qin < RNDOFFERR ? Ts : (GWtoChn*Td + T_LattoChn + L1toSrf * Ts + chan_prec * Ta) / qin;
    Tw = Tw_old  = _Temp_w->matrix[r][c];
    chan_rough = powl(10,_chan_roughness->matrix[r][c]);
    rough_LE = min<REAL8>(_random_roughness->matrix[r][c],_channelwidth->matrix[r][c]/chan_rough);
    za_LE = rough_LE + 2 ;
    z0u_LE = rough_LE * 0.1;
    zdu_LE = rough_LE * 0.7;
    rough_H = _random_roughness->matrix[r][c];
    za_H = rough_H + 2 ;
    z0u_H = max<REAL8>(0.000005,rough_H * 0.1);
    zdu_H = rough_H * 0.7;
  
    ra_LE = CalcAerodynResist(wind , za_LE , z0u_LE , zdu_LE , 0 , 0 , 0 , 0 ,Tw,Ta,ctrl.toggle_ra,true);
    ra_H  = CalcAerodynResist(wind , za_H , z0u_H , zdu_H , 0 , 0 , 0 , 0 ,Tw,Ta,ctrl.toggle_ra,true);
    // ----------------------------------------------------------------------------
    // Resistance due to overhanging vegetation
    // ----------------------------------------------------------------------------
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
				 Ta, ctrl.toggle_ra,false) * p;
      }
    }
    // ----------------------------------------------------------------------------
    // Energy balance factors that do not need updating in the N-R loop
    // ----------------------------------------------------------------------------
    fA = -4* 0.995 * stefboltz;	        //pools together net radiation factors ( emissivity of water is 0.995)
    fB = (-1/gamma) * (1/(ra_LE + rs + rc)) * rho_a * spec_heat_air; // pools together the latent heat factors
    fC = (-1/(ra_H)) * rho_a * spec_heat_air; // pools together the sensible heat factors
    fD = (-spec_heat_water*rho_w / (dt) * chan_store);

    trck.TracerMixing(*this,ctrl,chan_store,Tw_old,Tw_init,qin,Tin,Qk1,Tw,0,mixmod,r,c);
    // ----------------------------------------------------------------------------
    // Energy balance loop
    // ----------------------------------------------------------------------------
    int k = 0;
    do{
      Tw = _Temp_w->matrix[r][c];
      desdTw = 611 * ( (17.3/( Tw + 237.7)) - ((17.3 * Tw)/(powl(Tw + 237.2 , 2))) ) * exp(17.3 * Tw /( Tw + 237.7));

      LE = LatHeat(atm, beta, ra_LE, rc, 0.0, Tw, r, c); 			// latent heat flux
      evap = (LE < 0) ? std::min<REAL8>(LE/(rho_w*lambda) * dt,chan_store) : 0;
      Qx = spec_heat_water * rho_w / (dt * _channelwidth->matrix[r][c]) * 
                (_channelwidth->matrix[r][c] * (Tw_init * chan_store - Tw * 
                (chan_store - evap < RNDOFFERR ? 0.0 : chan_store - evap)) );   // Water stored energy
      H = SensHeat(atm, ra_H, Tw, r, c); 					// sensible heat flux
      fTw = NetRad(atm, Tw,  0, 0, AirEmissivity(Ta), Ta, 1 , r, c) + LE + H - Qx;  			// total energy balance
      dfTw = fA*powl(Tw + 273.2, 3) + fB * desdTw + fC - fD;
      Tw1 = Tw - (fTw/dfTw); 							// Update temperature
      _Temp_w->matrix[r][c] = Tw1;
      k++;
    }while(fabs(Tw1 - Tw) > 0.00001 && k < MAX_ITER);
    if (k>=MAX_ITER){
      std::cout << "WARNING: non-convergence in channel energy balance at cell row: " << r << " col: " << c << " closure err: " << (Tw1 - Tw) << endl;
    }
    // ----------------------------------------------------------------------------
    // Temperature correction and update
    // ----------------------------------------------------------------------------
    if (Tw1 > 35 || Tw1 < 0){ //correct the water temperature to a max/min value
      Tw1 = (Tw>35) ? Tw_init : 0;
      _Temp_w->matrix[r][c] = (Tw>35) ? Tw_init : 0;
      LE = LatHeat(atm, beta, ra_LE, rc, 0.0, Tw1, r, c); // latent heat flux
    } // end correction of temperature
    if(lat_ok) //update channel temperature
      _FTemp_w->matrix[rind][cind] += Qk1*(ctrl.dt/powl(_dx,2)) * _Temp_w->matrix[r][c];

  } else if(ctrl.toggle_chan_evap == 2){ //mass transfer approach
    es = SatVaporPressure(Ts);                     				// in Pa
    ea = SatVaporPressure(Ta) *atm.getRelativeHumidty()->matrix[r][c];          // in Pa

    b0 = 0; 									//generally set to value of 0
    b1 = powl(vonkarman,2) / powl( log((10 / 2.3e-3)) , 2) * molec_water_vap * 
 	   rho_a / (Pa * rho_w);   						// in Pa-1
    LE = -((b0 + b1 * wind) * (es - ea) * rho_w * lambda); 			//negative b/c outgoing
    ra_LE = 1;                                                                  //no additional resistance
    Tw1 = Ta; 									// only for tracer mixing
  }

  // ****************************************************************************
  // Make sure we don't evaporation too much
  // ****************************************************************************
  if( (-LE)>0 ) { //only evaporation if LE is negative
    evap = std::min<REAL8>(1/ra_LE , -LE/(rho_w*lambda));
    evap = (chan_store < evap * dt) ? chan_store/dt : evap;
  }
  chan_store -= evap * dt;
  _chan_store->matrix[r][c] = chan_store;
  _chan_evap->matrix[r][c] = evap;

  // ****************************************************************************
  // Flux tracking after evap
  // ****************************************************************************
  if(ctrl.sw_trck)    // Convert Tw to Kelvins for fractionation 
     trck.MixingV_evapW(atm, *this, ctrl, Tw1 , evap, chan_store, r, c);   

  return EXIT_SUCCESS;
}
