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
  int rr, cc;
  bool lat_ok;
  REAL8 dtdx;
  REAL8 alpha;
  REAL8 qj1i;
  REAL8 hji1;
  REAL8 hj1i1;
  REAL8 R;
  REAL8 dt = ctrl.dt;
  REAL8 poros1, poros2, poros3; //porosity
  REAL8 soildepth, d1, d2, d3; //guess
  REAL8 fc; //field capacity
  REAL8 deficit; //soil water deficit to reach field capacity in m

  //surface routing parameters
  REAL8 ponding = 0;
  REAL8 chan_store = 0;
  REAL8 theta1 = 0;
  REAL8 theta2 = 0;
  REAL8 theta3 = 0;
  REAL8 f = 0;
  REAL8 F = 0;
  REAL8 gw = 0; //gravitational water
  REAL8 returnflow = 0; //flow from gw in excess of the available soil storage
  REAL8 qc = 0; // water transfered from the subsurface system to the channel
  REAL8 qall = 0; //lateral inflows to channel
  REAL8 Qij1 = 0; //new discharge from the upstream boundary
  REAL8 Qk1 = 0; //new discharge out of the cell
  REAL8 Si1j1 = 0; //storage in channel at the end of time step

  REAL8 leak = 0;
  UINT4 reinf = 1;

  dtdx = dt / _dx;

  // Reinitialize to zero the fluxes modified earlier / in the previous time step
  _FluxExfilt->reset();
  _FluxL2toL1->reset();
  _FluxL3toL2->reset();
  if(ctrl.sw_trck){
    _FluxSrftoL1->reset();
    _FluxL1toL2->reset();
    _FluxL2toL3->reset();
    if(ctrl.sw_2H)
      trck.resetFd2HLat();
    if(ctrl.sw_18O)
      trck.resetFd18OLat();
    if(ctrl.sw_Age)
      trck.resetFAgeLat();
  }
  // Initialize others
  _FluxLattoSrf->reset();
  _FluxLattoChn->reset();
  _FluxLattoGW->reset();

  if(ctrl.sw_chan_evap)
    _FTemp_w->reset();

  // --------------------------------------------------------------------------------------
  for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {
    r = _vSortedGrid.cells[j].row;
    c = _vSortedGrid.cells[j].col;
    d = _vSortedGrid.cells[j].dir;

    //surface routing stuff
    returnflow = 0;
    Qij1 = _Disch_upstreamBC->matrix[r][c];
    qall = 0;
    ponding = _ponding->matrix[r][c];
    theta1 = _soilmoist1->matrix[r][c];
    theta2 = _soilmoist2->matrix[r][c];
    theta3 = _soilmoist3->matrix[r][c];
    //ca = _catcharea->matrix[r][c];
    gw = _GrndWater->matrix[r][c];
	  
    fc = _fieldcapL3->matrix[r][c];
    soildepth = _soildepth->matrix[r][c];
    d1 = _depth_layer1->matrix[r][c];
    d2 = _depth_layer2->matrix[r][c];
    d3 = soildepth - d1 - d2;
	  
    //if reinfiltration switch is on is not a channel cell or the channel switch is off     
    //if (ctrl.sw_reinfilt && !(ctrl.sw_channel && _channelwidth->matrix[r][c] > 0))
    if (ctrl.sw_reinfilt){

      Infilt_GreenAmpt(ctrl, f, F, theta1, theta2, theta3, ponding, gw, dt, r, c);
      // Tracking
      if(ctrl.sw_trck){
        // Mixing across the profile, accounting for snowmelt and lateral input
        trck.MixingV_down(*this, ctrl, d1, d2, d3, fc, r, c, 1);
        // Back up soil moisture before vertical redistrib
        _soilmoist1->matrix[r][c] = theta1; //soil moisture at t=t+1
        _soilmoist2->matrix[r][c] = theta2;
        _soilmoist3->matrix[r][c] = theta3;
        // Back up ponding before GW seepage (for rivers)
        _ponding->matrix[r][c] = ponding;
      } // end tracking
    }

    // For the rest of the routine, theta3 is only the content of the non-saturated
    // part of L3
    if (theta3 > fc) {
      gw = (theta3 - fc) * d3;
      theta3 = fc;
    } else
      gw = 0;

    if (ctrl.sw_channel && _channelwidth->matrix[r][c] > 0) { 
      //if this is a channel cell and channels are activated
      //maxR = ( _porosity->matrix[r][c] - fc ) * soildepth; 
      //calculates the maximum gravitational water that can go

      chan_store = _chan_store->matrix[r][c];

      qc = _KsatL3->matrix[r][c] * gw
	* (1 - expl(-_chGWparam->matrix[r][c] * gw));
      gw -= qc * dtdx;

    }

    _GravityWater->matrix[r][c] = gw;
	  
    ///enter groundwater
    poros1 = _porosityL1->matrix[r][c];
    poros2 = _porosityL2->matrix[r][c];
    poros3 = _porosityL3->matrix[r][c];
    alpha = _KsatL3->matrix[r][c] * sin(atan(_slope->matrix[r][c]));
	  
    deficit = 0;
    if (fabs(fc - theta3) > RNDOFFERR) {
      deficit = (fc - theta3) * d3;
    }
	  
    // discharge (j is timestep) so j1i is total upstream discharge per unit width at t+1
    qj1i = _GWupstreamBC->matrix[r][c];	
    //Not used since local GW head is embedded in the updated theta3 portion, becoming R	 
    hji1 = 0;
    //recharge to the groundwater system at the end of the time step in meters
    R = _GravityWater->matrix[r][c]; 
    //gravity water becomes groundwater
    _GravityWater->matrix[r][c] = 0; 
	  
    // Solution of the kinematic wave (hj1i1 = head "ready to go" downstream)
    hj1i1 = (dtdx * qj1i + hji1 + R - returnflow - deficit)
      / (1 + alpha * dtdx); //R is in meters so no need to multiply times dt here
	  
	  
    // If there's deficit and negative head -> capillary flow to L3, no GW outflow
    if (deficit > 0 && hj1i1 < 0) {
      theta3 += (dtdx * qj1i + hji1 + R - returnflow) / d3;
      hj1i1 = 0;
    // If there's deficit and positive head -> capillary flow to L3
    } else if (deficit > 0 && hj1i1 >= 0){
      theta3 += (dtdx * qj1i + hji1 + R - returnflow
		 - hj1i1 * (1 + alpha * dtdx)) / d3;
    }
	  
    // If the new amount of water in the cell is larger than the soil storage:
    // --> Solve for return flow
    if (((poros3 - theta3) * d3) < hj1i1) {
      returnflow = -(poros3 - theta3) * d3 * (1 + alpha * dtdx)	\
	+ (dtdx * qj1i) + hji1 + R - deficit;
      theta2 += returnflow / d2;
	    
      // Tracking (always on)
      _FluxL3toL2->matrix[r][c] = returnflow ;
	    
      if (theta2 > poros2) {
	theta1 += (theta2 - poros2) * d2 / d1;

	// Tracking (always on)
	_FluxL2toL1->matrix[r][c] = (theta2 - poros2) * d2 ;

	theta2 = poros2;
      }
      if (theta1 > poros1) {
	ponding += (theta1 - poros1) * d1;

	// Tracking (always on)
	_FluxExfilt->matrix[r][c] = (theta1 - poros1) * d1 ;

	theta1 = poros1;
      }
	    
      //ponding += returnflow;
      // Update head
      hj1i1 = (poros3 - theta3) * d3;
    }
	  
    //channel routing
    if (ctrl.sw_channel && _channelwidth->matrix[r][c] > 0) {

      chan_store += ponding + qc * dtdx;

      qall = chan_store * _dx / dt;
	    
      KinematicWave(Qk1, Si1j1, Qij1, qall, dt, r, c);
      //Qk1 = ponding * _dx*_dx/dt  + Qij1 ; // oooold (before kinematic wave)

      // Not all ponding water get routed once in the channel..
      // what is the actual amount of run-off that contributes to streamflow then?
      // For lack of a better solution, it is the amount of ponding that is effectively routed
      // AFTER all groundwater and upstream streamflow has been used
      // (since groundwater effectively enters the stream)
      _FluxGWtoChn->matrix[r][c] = qc*dtdx ;	    
      _FluxSrftoChn->matrix[r][c] = std::max<double>(0.0,(Qk1 - (Qij1 + qc*_dx))*dtdx/_dx);
      // Accumulated fluxes
      _AccGWtoChn->matrix[r][c] += _FluxGWtoChn->matrix[r][c];
      _AccSrftoChn->matrix[r][c] += _FluxSrftoChn->matrix[r][c];	  

      ponding = 0;
	    
      chan_store = 0;

    }	  

    // Locate downstream cell (if it exists)
    lat_ok = 0;
    switch (d) 
      {
      case 1:
	rr = r+1;
	cc = c-1;
	lat_ok = 1;
	break;
      case 2:
	rr = r+1;
	cc = c;
	lat_ok = 1;
	break;
      case 3:
	rr = r+1;
	cc = c+1;
	lat_ok = 1;
	break;
      case 4:
	rr = r;
	cc = c-1;
	lat_ok = 1;
	break;
      case 5: //if it is an outlet store the outflow m3s-1
	_dailyGwtrOutput.cells.push_back(cell(r, c, (alpha * hj1i1 * _dx)));
	_dailyOvlndOutput.cells.push_back(cell(r, c, Qk1+ponding * _dx *_dx / dt)); 
	//second term needed to account for outer at outlets with no channel	      
	break;
      case 6:
	rr = r;
	cc = c+1;
	lat_ok = 1;
	break;
      case 7:
	rr = r-1;
	cc = c-1;
	lat_ok = 1;
	break;
      case 8:
	rr = r-1;
	cc = c;
	lat_ok = 1;
	break;
      case 9:
	rr = r-1;
	cc = c+1;
	lat_ok = 1;
	break;
      default:
	return -1;
      }
	  
    // Check there is downstream cell
    if(lat_ok){

      // Input water for downstream cells (additive)
      _FluxLattoSrf->matrix[rr][cc] += ponding ;
      _FluxLattoChn->matrix[rr][cc] += Qk1*dtdx/_dx;
      _FluxLattoGW->matrix[rr][cc] += hj1i1 * alpha * dtdx;
      // Accumulated fluxes
      _AccLattoSrf->matrix[rr][cc] += ponding ;
      _AccLattoChn->matrix[rr][cc] += Qk1*dtdx/_dx ;
      _AccLattoGW->matrix[rr][cc] += hj1i1 * alpha * dtdx;
    
      // Add the previously calculated *discharge* (not elevation) to the downstream cell
      _GWupstreamBC->matrix[rr][cc] += hj1i1 * alpha;
      _Disch_upstreamBC->matrix[rr][cc] += Qk1;
      _ponding->matrix[rr][cc] += ponding;	
      if(ctrl.sw_channel && _channelwidth->matrix[rr][cc])
	_chan_store->matrix[rr][cc] += chan_store;
    }

    // Outgoing water (outside of lat_ok because can be 0)
    //if(ctrl.sw_trck){
    _FluxSrftoLat->matrix[r][c] = lat_ok != 0 ? ponding : 0.0;
    _FluxGWtoLat->matrix[r][c] = lat_ok != 0 ? hj1i1 * alpha * dtdx : 0.0;
    _FluxChntoLat->matrix[r][c] = lat_ok != 0 ? Qk1*dtdx / _dx : 0.0;
    // Accumulated fluxes
    _AccSrftoLat->matrix[r][c] += _FluxSrftoLat->matrix[r][c];
    _AccGWtoLat->matrix[r][c] += _FluxGWtoLat->matrix[r][c];
    _AccChntoLat->matrix[r][c] += _FluxChntoLat->matrix[r][c];

    // Tracking of lateral in/out + return + seepage
    if(ctrl.sw_trck){
      trck.MixingV_latup(*this, ctrl, d1, d2, d3, fc, 
			 Qk1, dtdx, _dx, r, c);
      // Tracking lateral inputs to the downstream cell
      // Summed tracking contribution downstream cells (for mixing)
      if(lat_ok == 1)
	trck.FCdownstream(*this, ctrl, Qk1, dtdx, _dx, r, c, rr, cc);
      else
	// Catchment outlets' values
	trck.OutletVals(ctrl, 1, r, c);
    }

    // Update ponding and water contents
    if (ctrl.sw_channel && _channelwidth->matrix[r][c] > 0){
      _chan_store->matrix[r][c] = Si1j1 / (_dx * _dx);
      _ponding->matrix[r][c] = 0.0;
      if (ctrl.sw_chan_evap){
	SolveChannelEnergyBalance(atm,ctrl,trck,0,Qk1,r,c);
	if (lat_ok){
	  _FTemp_w->matrix[rr][cc] += Qk1*dtdx/_dx * _Temp_w->matrix[r][c];
	}
      } else {
	_chan_evap->matrix[r][c] = 0.0; // no channel evaporation
	if (ctrl.sw_trck)
	  trck.MixingV_evapW(atm,*this,ctrl,atm.getTemperature()->matrix[r][c],
			     _chan_evap->matrix[r][c],_chan_store->matrix[r][c],r,c);//to get a value for evap isotopes (no frac)
      }
    } else {
      _ponding->matrix[r][c] = 0.0;
      _chan_store->matrix[r][c] = 0.0;
      _chan_evap->matrix[r][c] = 0.0;
      if (ctrl.sw_chan_evap){
	_Temp_w->matrix[r][c] = 0.0;
      }
    }

    _soilmoist1->matrix[r][c] = theta1;
    _soilmoist2->matrix[r][c] = theta2;
    _soilmoist3->matrix[r][c] = theta3 + hj1i1 / d3;
    _GrndWater->matrix[r][c] = hj1i1;

    // Save river discharge
    _Disch_old->matrix[r][c] = Qk1;
    Qk1 = 0;

    // Accumulated fluxes
    _AccInfilt->matrix[r][c] += _FluxInfilt->matrix[r][c];
    _AccExfilt->matrix[r][c] += _FluxExfilt->matrix[r][c];
    _AccPercolL2->matrix[r][c] += _FluxPercolL2->matrix[r][c];
    _AccL2toL1->matrix[r][c] += _FluxL2toL1->matrix[r][c];
    _AccPercolL3->matrix[r][c] += _FluxPercolL3->matrix[r][c];
    _AccL3toL2->matrix[r][c] += _FluxL3toL2->matrix[r][c];
    _AccRecharge->matrix[r][c] += _FluxRecharge->matrix[r][c];
    _AccEvaporationS->matrix[r][c] += _EvaporationS_all->matrix[r][c] *dt;
  }
	
  // Save previous GW and surface state
  *_GrndWater_old = *_GrndWater;
  *_ponding_old = *_ponding;
  *_chan_store_old = *_chan_store;

  return EXIT_SUCCESS;
}
