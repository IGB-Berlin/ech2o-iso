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
 * GWrouting.cpp
 *
 *  Created on: Dec 2, 2010
 *      Author: Marco.Maneta
 */

#include"Basin.h"

int Basin::DailyGWRouting(Atmosphere &atm, Control &ctrl, Tracking &trck) {

  int r, c, d;
  int rind, cind;
  bool lat_ok;
  REAL8 dtdx;
  REAL8 alpha;
  REAL8 qj1i,R,hji1,hj1i1; //groundwater flow and head
  REAL8 Deep_qj1i,Deep_R,Deep_hj1i1;
  REAL8 dt = ctrl.dt;
  REAL8 poros1, poros2, poros3; //porosity
  REAL8 soildepth, d1, d2, d3; //guess
  REAL8 fc; //field capacity
  REAL8 deficit; //soil water deficit to reach field capacity in m

  //surface routing parameters
  REAL8 contArea = 0; //contributing area
  REAL8 ponding = 0;
  REAL8 chan_store = 0;
  REAL8 dcdx = 0; //ratio of channel length to cell size
  REAL8 theta1,theta2,theta3;
  theta1 = theta2 = theta3 = 0;
  REAL8 DeepGW_all = 0; //Deep groundwater volume
  REAL8 DeepGW;
  REAL8 Fhydro = 0; //hydrologically active proportion of deep water
  REAL8 f = 0;
  REAL8 F = 0;
  REAL8 gw = 0; //gravitational water
  REAL8 returnflow = 0; //flow from gw in excess of the available soil storage
  REAL8 qc = 0; // water transfered from the subsurface system to the channel
  REAL8 Deep_qc = 0; //water transfered from deep GW to channel
  REAL8 qall = 0; //lateral inflows to channel
  REAL8 Qij1 = 0; //new discharge from the upstream boundary
  REAL8 Qk1 = 0; //new discharge out of the cell
  REAL8 Si1j1 = 0; //storage in channel at the end of time step
  
  REAL8 leak = 0;

  dtdx = dt / _dx;

  // Reinitialize to zero the fluxes modified earlier / in the previous time step
  _FluxExfilt->reset();
  //  resetLateralfluxes(ctrl,1);
  _FluxLattoSrf->reset();
  _FluxLattoChn->reset();      
  _FluxLattoGW->reset();
  if(ctrl.sw_trck){
    resetTrckfluxes(ctrl,1);
    trck.resetFTracerLat(ctrl);
  }

  // --------------------------------------------------------------------------------------
  for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
    r = _vSortedGrid.cells[j].row;
    c = _vSortedGrid.cells[j].col;
    d = _vSortedGrid.cells[j].dir;

    //surface routing stuff
    returnflow = 0;
    Qij1 = _Disch_upstreamBC->matrix[r][c]; 						//[m3/s]
    qall = 0; 										//[m2/s]
    ponding = _ponding->matrix[r][c]; 							//[m]
    chan_store = _chan_store->matrix[r][c];
    theta1 = _soilmoist1->matrix[r][c];
    theta2 = _soilmoist2->matrix[r][c];
    theta3 = _soilmoist3->matrix[r][c];
    gw = _GrndWater->matrix[r][c]; 							//[m]
	  
    contArea = powl(_dx,2) * _ttarea->matrix[r][c]; 					//[m2] - Contributing area of pixel

    fc = _fieldcapL3->matrix[r][c]; 
    soildepth = _soildepth->matrix[r][c]; 						//[m]
    d1 = _depth_layer1->matrix[r][c]; 							//[m]
    d2 = _depth_layer2->matrix[r][c]; 							//[m]
    d3 = soildepth - d1 - d2; 								//[m]

    // ****************************************************************************
    // Infiltrate Water
    // if reinfiltration switch is on
    // ****************************************************************************  
    if (ctrl.sw_reinfilt)
      Infilt_GreenAmpt(ctrl, f, F, theta1, theta2, theta3, ponding, gw, dt, r, c);
    
    // ****************************************************************************
    // Always mixes when tracking is on 
    //   snowmelt and later input (pond) and vertical when reinfilt is on
    // ****************************************************************************
    if(ctrl.sw_trck)	      // Mixing across the profile, accounting for snowmelt and lateral input
      trck.MixingV_down(*this, ctrl, d1, d2, d3, fc, r, c, 1);

    // Back up soil moisture before vertical redistrib
    _soilmoist1->matrix[r][c] = theta1; 						//soil moisture at t=t+1
    _soilmoist2->matrix[r][c] = theta2;							//soil moisture at t=t+1
    _soilmoist3->matrix[r][c] = theta3;							//soil moisture at t=t+1
    _ponding->matrix[r][c] = ponding; 

    // ****************************************************************************
    // GW to the Channel (if there is a channel)
    //   For the rest of the routine, theta3 is only the content of the non-saturated part of L3
    //   Calculates the maximum gravitational water that can go
    // ****************************************************************************
    gw = (theta3 > fc) ? (theta3 - fc) * d3 : 0 ;					//[m]
    theta3 = (theta3 > fc) ? fc : theta3;

    if (ctrl.sw_channel && _channelwidth->matrix[r][c] > 0) {
      dcdx = _channellength->matrix[r][c]/_dx;                                          //[m/m]
      qc = _KsatL3->matrix[r][c] * gw * (1 - expl(-_chGWparam->matrix[r][c] * gw)); 	//[m2/s]
      gw -= qc * dtdx * dcdx; 				                                //[m]
    }
    _GravityWater->matrix[r][c] = gw; 	//Update groundwater
	  
    // ****************************************************************************
    // Solve groundwater flow (j is timestep, i is spacestep)
    // ****************************************************************************
    poros1 = _porosityL1->matrix[r][c];
    poros2 = _porosityL2->matrix[r][c];
    poros3 = _porosityL3->matrix[r][c];
    alpha = _KsatL3->matrix[r][c] * sin(atan(_slope->matrix[r][c])); 			//[m/s]
	  
    deficit = (fabs(fc - theta3) > RNDOFFERR) ? (fc - theta3) * d3 : 0; 		//[m]

    qj1i = _GWupstreamBC->matrix[r][c];	//upstream discharge/unit width at t+1
    hji1 = 0;    			//Not used - local GW head is embedded in the theta3 portion
    R = _GravityWater->matrix[r][c];    //recharge to the groundwater system at t+1 [m]
    _GravityWater->matrix[r][c] = 0;    //gravity water becomes groundwater
	  
    // Solution of the kinematic wave (hj1i1 = head "ready to go" downstream)
    hj1i1 = (dtdx * qj1i + hji1 + R - returnflow - deficit)
      / (1 + alpha * dtdx); //R is in meters so no need to multiply times dt here
   
    // If there's deficit and negative head -> capillary flow to L3, no GW outflow
    if (deficit > 0 && hj1i1 < 0) {
      theta3 += (dtdx * qj1i + hji1 + R - returnflow) / d3;
      hj1i1 = 0;
    // If there's deficit and positive head -> capillary flow to L3
    } else if (deficit > 0 && hj1i1 >= 0){
      theta3 += (dtdx * qj1i + hji1 + R - returnflow - hj1i1 * (1 + alpha * dtdx)) / d3;
    }

    // ****************************************************************************
    // Solve for Returnflow
    // If the new amount of water in the cell is larger than the soil storage:
    // ****************************************************************************
    if (((poros3 - theta3) * d3) < hj1i1) {
      returnflow = -(poros3 - theta3) * d3 * (1 + alpha * dtdx)	\
	+ (dtdx * qj1i) + hji1 + R - deficit;
      theta2 += returnflow / d2;
      _ReturnL2->matrix[r][c] = returnflow;
      if (theta2 > poros2) {
	theta1 += (theta2 - poros2) * d2 / d1;
	_ReturnL1->matrix[r][c] = (theta2 - poros2) * d2 ;
	theta2 = poros2;
      }
      if (theta1 > poros1) {
	ponding += (theta1 - poros1) * d1;
	_FluxExfilt->matrix[r][c] = (theta1 - poros1) * d1 ;	// Tracking (always on)
	theta1 = poros1;
      }
      // Update head
      hj1i1 = (poros3 - theta3) * d3; 							//[m]
    }
	
    // ****************************************************************************  
    // Deep GW terrestrial dynamics
    // ****************************************************************************
    Deep_hj1i1 = 0;
    if (ctrl.sw_deepGW){
      DeepGW_all = _DeepGW->matrix[r][c];
      Fhydro = _Hydrofrac_DeepGW->matrix[r][c];
      leak = _BedrockLeakageFlux->matrix[r][c] * dt;            				//[m]
      DeepGW = DeepGW_all * Fhydro + leak;							//[m]
      if (ctrl.sw_channel && _channelwidth->matrix[r][c] > 0) { 
        Deep_qc = _KsatL3->matrix[r][c] * DeepGW * 
			(1 - expl(-_chDeepGWparam->matrix[r][c] * DeepGW));			//[m2/s]
        DeepGW -= Deep_qc * dtdx * dcdx; 				                        //[m]
      }
      Deep_qj1i = _DeepGWupstreamBC->matrix[r][c]; 						//[m2/s]
      Deep_R = DeepGW; 										//[m]
      Deep_hj1i1 = (dtdx*Deep_qj1i + Deep_R)/(1+alpha*dtdx); 					//[m]
    }
    // ****************************************************************************  
    // Channel routing
    // ****************************************************************************
    if (ctrl.sw_channel && _channelwidth->matrix[r][c] > 0) {
      chan_store = _chan_store->matrix[r][c] * _dx/dt;           				//[m2/s]
      chan_store += ponding * contArea / (_dx*dt);                                              //[m2/s]
      // ****************************
      //chan_store += _ttarea->matrix[r][c] * ((qc + Deep_qc) * dcdx) ;                           //[m2/s]
      // added m in channel --> (qc * dtdx * dcdx)
      // convert to m2/s    --> (qc * dtdx * dcdx) * contArea / (_dx*dt)
      chan_store += ((qc +Deep_qc)* dtdx * dcdx) * contArea / (_dx*dt);
      // ****************************      
      qall = chan_store ;					                                //[m2/s]
      KinematicWave(Qk1, Si1j1, Qij1, qall, dt, r, c); 						//[m3/s]

      // Not all channel water get routed...the amount of water routed AFTER all GW and U/S Q 
      _FluxGWtoChn->matrix[r][c] = qc*dtdx *(dcdx); 		                                //[m]
      _FluxSrftoChn->matrix[r][c] = ponding;                                                    //[m]
      
      _AccGWtoChn->matrix[r][c] += _FluxGWtoChn->matrix[r][c]; 					//[m]
      _AccSrftoChn->matrix[r][c] += _FluxSrftoChn->matrix[r][c]; 				//[m]	  

      if (ctrl.sw_deepGW){
        _FluxDeepGWtoChn->matrix[r][c] = Deep_qc*dtdx*dcdx;
        _AccDeepGWtoChn->matrix[r][c] += _FluxDeepGWtoChn->matrix[r][c];
      }
      ponding = 0; 		//all ponded water moved to channel
      chan_store = 0; 		//channel storage updated later
    }	  

    // ****************************************************************************
    // Locate downstream cell (if it exists) and move water
    // ****************************************************************************
    lat_ok = 0;
    rind = r + rr(d);
    cind = c+ cc(d);
    if(d==5){
      _dailyGwtrOutput.cells.push_back(cell(r, c, (alpha * hj1i1 * _dx))); 			//[m3/s]
      _dailyOvlndOutput.cells.push_back(cell(r, c, Qk1+ponding * contArea / dt));  		//[m3/s]
      if(ctrl.sw_deepGW)
	_dailyDeepGwtrOutput.cells.push_back(cell(r, c, (alpha * Deep_hj1i1 * _dx))); 		//[m3/s]	
    } else {
      lat_ok = 1;
      _FluxLattoSrf->matrix[rind][cind] += ponding ; 						//[m]
      _FluxLattoChn->matrix[rind][cind] += Qk1*dtdx/_dx; 					//[m]
      _FluxLattoGW->matrix[rind][cind] += hj1i1 * alpha * dtdx; 				//[m]
      // Accumulated fluxes
      _AccLattoSrf->matrix[rind][cind] += ponding ; 						//[m]
      _AccLattoChn->matrix[rind][cind] += Qk1*dtdx/_dx ; 					//[m]
      _AccLattoGW->matrix[rind][cind] += hj1i1 * alpha * dtdx; 					//[m]
      if(ctrl.sw_deepGW){
       _FluxLattoDeepGW->matrix[rind][cind] += Deep_hj1i1 * alpha * dtdx; 			//[m]
       _AccLattoDeepGW->matrix[rind][cind] += Deep_hj1i1 * alpha * dtdx; 			//[m]
       _DeepGWupstreamBC->matrix[rind][cind] += Deep_hj1i1 * alpha; 				//[m2/s]
      }
      // Add the previously calculated *discharge* (not elevation) to the downstream cell
      _GWupstreamBC->matrix[rind][cind] += hj1i1 * alpha; 					//[m2/s]
      _Disch_upstreamBC->matrix[rind][cind] += Qk1; 						//[m3/s]
      _ponding->matrix[rind][cind] += ponding;	 						//[m]
      if(ctrl.sw_channel && _channelwidth->matrix[rind][cind]){
	_chan_store->matrix[rind][cind] += chan_store;                                          //[m]
      }
    }

    // Outgoing water (outside of lat_ok because can be 0)
    _FluxSrftoLat->matrix[r][c] = lat_ok != 0 ? ponding : 0.0; 					//[m]
    _FluxGWtoLat->matrix[r][c] = lat_ok != 0 ? hj1i1 * alpha * dtdx : 0.0; 			//[m]
    _FluxChntoLat->matrix[r][c] = lat_ok != 0 ? Qk1*dtdx / _dx : 0.0; 				//[m]
    // Accumulated fluxes
    _AccSrftoLat->matrix[r][c] += _FluxSrftoLat->matrix[r][c]; 					//[m]
    _AccGWtoLat->matrix[r][c] += _FluxGWtoLat->matrix[r][c]; 					//[m]
    _AccChntoLat->matrix[r][c] += _FluxChntoLat->matrix[r][c]; 					//[m]
    if (ctrl.sw_deepGW){
      _FluxDeepGWtoLat->matrix[r][c] = lat_ok != 0 ? Deep_hj1i1 * alpha * dtdx : 0.0; 		//[m]		
      _AccLattoDeepGW->matrix[r][c] += _FluxDeepGWtoLat->matrix[r][c]; 				//[m]
    } 

    // ****************************************************************************
    // Tracking of lateral in/out + return + seepage + channel evaporation
    // ****************************************************************************
    if(ctrl.sw_trck){
      trck.MixingV_latup(*this, ctrl, d1, d2, d3, fc,Qk1, dtdx, _dx, r, c);
      if(lat_ok == 1)      	// Tracking lateral inputs to the downstream cell
	trck.FCdownstream(*this, ctrl, Qk1, dtdx, _dx, r, c, rind, cind);
      else			// Catchment outlets' values
	trck.OutletVals(ctrl, 1, r, c);		
    }

    // ****************************************************************************
    // Update ponding, channel, and evaporation
    // ****************************************************************************
    _ponding->matrix[r][c] = 0.0;
    if (ctrl.sw_channel && _channelwidth->matrix[r][c] > 0){
      _chan_store->matrix[r][c] = Si1j1 / (_dx * _dx) ;						//[m]
      ChannelEvaporation(atm,ctrl,trck,0,Qk1,lat_ok,r,c,rind,cind);
    } else {
      _chan_store->matrix[r][c] = 0.0; 								//[m]
      _chan_evap->matrix[r][c] = 0.0; 								//[m]
      if (ctrl.toggle_chan_evap == 1) // only if energy balance approach is used
	_Temp_w->matrix[r][c] = 0.0; 								//[oC]
    }

    // ****************************************************************************
    // Update moisture and accumulated fluxes
    // ****************************************************************************
    _soilmoist1->matrix[r][c] = theta1;
    _soilmoist2->matrix[r][c] = theta2;
    _soilmoist3->matrix[r][c] = theta3 + hj1i1 / d3;
    _GrndWater->matrix[r][c] = hj1i1;
    if(ctrl.sw_deepGW)
      _DeepGW->matrix[r][c] = Deep_hj1i1 + DeepGW_all*(1-Fhydro);

    _Disch_old->matrix[r][c] = Qk1; 		        // Save river discharge
    Qk1 = 0; 						// reset discharge

    _AccInfilt->matrix[r][c] += _FluxInfilt->matrix[r][c];
    _AccExfilt->matrix[r][c] += _FluxExfilt->matrix[r][c];
    _AccPercolL2->matrix[r][c] += _FluxPercolL2->matrix[r][c];
    _AccL2toL1->matrix[r][c] += _ReturnL1->matrix[r][c];
    _AccPercolL3->matrix[r][c] += _FluxPercolL3->matrix[r][c];
    _AccL3toL2->matrix[r][c] += _ReturnL2->matrix[r][c];
    _AccRecharge->matrix[r][c] += _FluxRecharge->matrix[r][c];
    _AccEvaporationS->matrix[r][c] += _EvaporationS_all->matrix[r][c] *dt;

  } // end for loop
	
  *_chan_store_old = *_chan_store;

  return EXIT_SUCCESS;
}
