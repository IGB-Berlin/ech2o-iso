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
 * SolveEnergyBalance.cpp
 *
 *  Created on: Nov 20, 2009
 *      Author: marco.maneta
 */

#include"Basin.h"

int Basin::SolveSurfaceFluxes(Atmosphere &atm, Control &ctrl, Tracking &trck) {

  int r, c, d;
  float dt = ctrl.dt; //time step
  
  //energy balance parameters
  
  REAL8 ra; 			//soil aerodynamic resistance
  REAL8 rs; 			//bare soil resistance (a function of soil moisture)
  REAL8 Ta = 0; 		//air temperature
  REAL8 Ts = 0; 		//soil surface temperature
  REAL8 Tsold = 0; 		//old surface temperature
  REAL8 TdL1  = 0; 		//temperature of bottom of layer 1
  REAL8 TdL2  = 0; 		//temperature of bottom of layer 2
  REAL8 TdL3  = 0; 		//temperature of bottom of layer 3
  REAL8 Tdold = 0; 		//temperature of lower soil thermal layer
  
  REAL8 LAI = 0;
  REAL8 BeersK = 0;
  REAL8 Temp_can = 0; 		//temperature of the canopy
  REAL8 emis_can = 0; 		//emissivity of the canopy
  
  REAL8 evap = 0; 		//evaporation
  
  //infiltration parameters
  REAL8 infcap = 0;
  REAL8 accinf = 0;
  REAL8 theta = 0; 		//soil moisture for entire soil profile or for first soil layer
  REAL8 theta2 = 0; 		//for second and third soil moisture
  REAL8 theta3 = 0; 		//layers in case Richard's equation is chosen
  REAL8 ponding = 0;
  REAL8 gw = 0; 		//gravitational water

  double d1, d2, d3; 		// soil layers' depths
  double fc; 			//field capacity in third layer (for tracking)
  
  //aerodynamic resistance parameters
  REAL8 za; 			//height of wind speed measurements
  REAL8 z0u; 			// roughness length for understory
  REAL8 zdu; 			//zero plane displacement for understory
  REAL8 z0o; 			// roughness length for overrstory
  REAL8 zdo; 			//zero plane displacement for overstory
  
  REAL8 wind; 			//wind speed
  REAL8 treeheight;
  
  REAL8 nr, le, sens, grndh, snowh, mltht, dh_snow, etc;
  
  UINT4 nsp;
  REAL8 p;			//fraction of species s
  int bsoil;                    //switch for if bare soil
  REAL8 BC_deepGW;              //deep groundwater boundary condition
  
  //needed in the water routing routines
  resetBCfluxes(ctrl);
  resetVerticalfluxes();
  resetLateralfluxes(ctrl,0);
  resetTrckfluxes(ctrl,0);
  
   if(ctrl.sw_trck){
    _FluxSnowtoSrf->reset();
    resetTrckfluxes(ctrl,0);
    trck.resetFTracerLat(ctrl);
    trck.OutletVals(ctrl, 0, 0, 0);    
  }

#pragma omp parallel default(shared) \
  private(r, c, ra, rs, Ta, Ts, Tsold, Tdold,TdL1,TdL2,TdL3, LAI, BeersK, Temp_can, emis_can, \
	  evap, infcap, accinf, theta, theta2, theta3,ponding,gw, za, z0u, zdu, z0o, \
	  zdo, wind, treeheight,nr, le, sens, grndh, snowh, mltht, dh_snow, p, etc,d1, d2, d3,\
	  fc, d, bsoil)
  {
#pragma omp for nowait
    for (unsigned int j = 0; j < _vSortedGrid.cells.size() ; j++)
      {
	r = _vSortedGrid.cells[j].row;
	c = _vSortedGrid.cells[j].col;
	d = _vSortedGrid.cells[j].dir;	

        //*****************************************************************************************
	//--- Set the boundary conditions ---------------------------------------------------------
        //*****************************************************************************************
	BC_deepGW = (ctrl.sw_deepGW and ctrl.sw_BC) ? _BCdeepgwtr->matrix[r][c] : 0;
	if(ctrl.sw_BC == 1){
	  if(ctrl.sw_trck){ // if there is tracking
	    trck.BCupstreamMixing(*this,ctrl,_BCsurface->matrix[r][c],
				  _BCgroundwater->matrix[r][c], BC_deepGW,_dx,dt,r,c); //TODO update tracking isotopes & deepGW
	  }
	  _GWupstreamBC->matrix[r][c]     += _BCgroundwater->matrix[r][c];   // [m2.s-1] per unit width
	  // Surface water is added directly to ponded/channel water rather than flux
	  if(ctrl.sw_channel && _channelwidth->matrix[r][c]) {
	    _chan_store->matrix[r][c]     += _BCsurface->matrix[r][c] * dt / (_dx * _dx); // [m]
	  } else {
	    _ponding->matrix[r][c]       += _BCsurface->matrix[r][c] * dt / (_dx * _dx); // [m]
	  }
	  _FluxGWtoLat->matrix[r][c]      += _BCgroundwater->matrix[r][c]*dt/_dx;     // [m]
	  if(ctrl.sw_deepGW){
	    _DeepGWupstreamBC->matrix[r][c] += _BCdeepgwtr->matrix[r][c];    // [m2.s-1]
	    _FluxDeepGWtoLat->matrix[r][c] += _BCdeepgwtr->matrix[r][c] * dt/_dx;// [m]
	  }
	}

        //*****************************************************************************************	
	//--- Initialization of conditions --------------------------------------------------------
        //*****************************************************************************************
	wind = atm.getWindSpeed()->matrix[r][c];
	
	theta = _soilmoist1->matrix[r][c]; //soil moisture at time t
	theta2 = _soilmoist2->matrix[r][c];
	theta3 = _soilmoist3->matrix[r][c];
	ponding = _ponding->matrix[r][c]; //surface ponding at time t
	gw = _GrndWater->matrix[r][c]; //gravity water at time t

	fc = _fieldcapL3->matrix[r][c];
	d1 = _depth_layer1->matrix[r][c];
	d2 = _depth_layer2->matrix[r][c];
	d3 = _soildepth->matrix[r][c] - d1 - d2;
	
	nr = le = sens = grndh = snowh = mltht = 0;
	Ts = _Temp_s->matrix[r][c];
	Ta = atm.getTemperature()->matrix[r][c];
	Tsold = Tdold = dh_snow = BeersK = TdL1 = TdL2 = TdL3 = 0;

	nsp = fForest->getNumSpecies();
	treeheight = 0;

        //*****************************************************************************************
     	//--- Solve the energy balance ------------------------------------------------------------
        //*****************************************************************************************
	for(UINT4 s = 0; s < nsp ; s++) {
	    p = fForest->getPropSpecies(s, r, c);
	    if(p == 0)
	      continue;//if no species j present, continue
	    if(s == nsp -1){
	      //for bare soil, water reaching the ground is Pp times its proportion of the cell
	      // Impervious fractions affect bare soil
	      LAI = emis_can = Temp_can = 0;
	      za = _random_roughness->matrix[r][c] + 2;
	      z0u = max<REAL8>(0.000005,_random_roughness->matrix[r][c] * 0.1);
	      zdu = _random_roughness->matrix[r][c] * 0.7;
	      z0o = 0; //no overstory
	      zdo = 0;
	      bsoil = 1;
	    } else {
	      treeheight = max<REAL8>(0.01,fForest->getTreeHeight(s, r, c)); //equations only apply to 40% of the tree as per Campbell and Norman 1998
	      LAI = fForest->getLAISpecies(s, r, c);
	      BeersK = fForest->getBeersCoeff(s, r, c);
	      Temp_can = fForest->getCanopyTemp(s, r, c);
	      emis_can = fForest->getCanopyEmissivity(s, r, c);
	      za = treeheight + 2;
	      z0o = powl(treeheight,1.19)*0.057544;     //treeheight > 1 ? 0.1 : treeheight * 0.1;
	      zdo = powl(treeheight,0.98)*0.707946; //treeheight > 1 ? 0.1 : treeheight * 0.7;
 	      zdu = min<double>(_random_roughness->matrix[r][c], zdo * 0.1);
	      z0u = 0.1*zdu/0.7;
	      bsoil = 0;
	    }
	    
	    ra = CalcAerodynResist(wind, za, z0u, zdu, z0o, zdo, treeheight, LAI, Ts, 
				  Ta, ctrl.toggle_ra,true);
	    rs = CalcSoilResist(theta, r, c, ctrl.toggle_rs);
	    
	    SolveSurfaceEnergyBalance(atm, ctrl, trck, ra, rs, 0.0, BeersK, LAI,
				      emis_can, Temp_can, nr, le, sens, grndh, snowh,
				      mltht, Tsold, evap, ponding, theta, Ts,
				      Tdold, TdL1, TdL2, TdL3, p, bsoil, r, c, s);

	    _Evaporation->matrix[r][c] += evap; 		//evaporation at t=t+1 (w by p)
	    _EvaporationS_all->matrix[r][c] += evap; 		//soil evaporation at t=t+1 (w by p)
	    if(s != nsp -1){
	      fForest->setEsoilSpecies(s, r, c, evap/p);	// Esoil
	      etc = fForest->getEvapoTransp(s, r, c);           // ET-soil
	      fForest->setETSpecies(s, r, c, etc + evap/p);     // ET total
	    }
	}//for over the species

	_incident_water_depth->matrix[r][c] += SnowOutput(atm, ctrl, trck, mltht, r, c);
	
        //*****************************************************************************************
	//--- Infiltration + percolation if exceeds porosity --------------------------------------
        //*****************************************************************************************
	if(ctrl.toggle_hydrologic_engine > 0){
	  RichardsEquation(ctrl,infcap,accinf, theta, theta2, theta3, ponding, dt, r,c,d);
	} else{
	  ponding += _incident_water_depth->matrix[r][c];
	  _incident_water_depth->matrix[r][c] = 0;
          Infilt_GreenAmpt(ctrl, infcap, accinf, theta, theta2, theta3, ponding, gw, dt, r, c);
          // Percolation if exceeding field capacity (L1 and L2), 
          // goes to GW in L3 (and bedrock leakage if activated)
          SoilWaterRedistribution(ctrl, accinf, theta, theta2, theta3, ponding, gw, dt, r, c);
          if(ctrl.sw_trck)
            trck.MixingV_down(*this, ctrl, d1, d2, d3, fc, r, c, 0); //TODO update tracking here
        }

	// Calculates the soil moisture profile to derive equivalent water table depth
	if(ctrl.Rep_WaterTableDepth == 1 || ctrl.RepTs_WaterTableDepth == 1)
	  CalcSoilMoistureProfile(atm, ctrl, getSoilMoist_av()->matrix[r][c], r,c);

        //*****************************************************************************************
	//--- Update global objects ---------------------------------------------------------------
        //*****************************************************************************************
	_soilmoist1->matrix[r][c] = theta; //soil moisture at t=t+1
	_soilmoist2->matrix[r][c] = theta2;
	_soilmoist3->matrix[r][c] = theta3;

	_Rn->matrix[r][c] = nr;
	_latheat->matrix[r][c] = le;
	_latheat_sum->matrix[r][c] += le;
	_sensheat->matrix[r][c] = sens;
	_sensheat_sum->matrix[r][c] += sens;
	_grndheat->matrix[r][c] = grndh;
	_snwheat->matrix[r][c] = snowh;
	_Temp_s_old->matrix[r][c] = Tsold;
	_Temp_s->matrix[r][c] = Tsold; //
	
	_Temp_d->matrix[r][c] = Tdold;
	_Temp_L1->matrix[r][c] = TdL1;
	_Temp_L2->matrix[r][c] = TdL2;
	_Temp_L3->matrix[r][c] = TdL3;	
	
      }//for all cells
  }//end omp parallel block
  return EXIT_SUCCESS;
}
