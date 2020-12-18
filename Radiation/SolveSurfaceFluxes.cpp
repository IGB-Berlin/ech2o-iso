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
  
  REAL8 ra; //soil aerodynamic resistance
  REAL8 rs; //bare soil resistance (a function of soil moisture)
  REAL8 Ts = 0; //
  REAL8 Tsold = 0; //old surface temperature
  REAL8 Tdold = 0; //temperature of lower soil thermal layer
  
  REAL8 LAI = 0;
  REAL8 BeersK = 0;
  REAL8 Temp_can = 0; //temperature of the canopy
  REAL8 emis_can = 0; //emissivity of the canopy
  
  REAL8 evap = 0; //evaporation
  
  //infiltration parameters
  REAL8 infcap = 0;
  REAL8 accinf = 0;
  REAL8 theta = 0; //soil moisture for entire soil profile or for first soil layer
  REAL8 theta2 = 0; //for second and third soil moisture
  REAL8 theta3 = 0; //layers in case Richard's equation is chosen
  REAL8 ponding = 0;
  REAL8 gw = 0; //gravitational water
  REAL8 leak = 0; //bedrock leakage flux;

  UINT4 reinf;
  
  double d1, d2, d3; // soil layers' depths
  double fc; //field capacity in third layer (for tracking)
  
  //aerodynamic resistance parameters
  REAL8 za; //height of wind speed measurements
  REAL8 z0u; // roughness length for understory
  REAL8 zdu; //zero plane displacement for understory
  REAL8 z0o; // roughness length for overrstory
  REAL8 zdo; //zero plane displacement for overstory
  
  REAL8 wind; //wind speed
  REAL8 treeheight;
  
  REAL8 nr, le, sens, grndh, snowh, mltht, dh_snow, etc;
  
  UINT4 nsp;
  REAL8 p;//fraction of species s

  //needed in the water routing routines
  _dailyOvlndOutput.cells.clear();
  _dailyL1wtrOutput.cells.clear();
  _dailyL2wtrOutput.cells.clear();
  _dailyGwtrOutput.cells.clear();
  _Layer1Outlet->reset();
  _Layer2Outlet->reset();
  _GWOutlet->reset();
  _Layer1upstreamBC->reset();
  _Layer2upstreamBC->reset();
  _GWupstreamBC->reset();
  _Disch_upstreamBC->reset();
  // Infilt, Percol and Recharge
  _FluxInfilt->reset();
  _FluxExfilt->reset();
  _FluxPercolL2->reset();
  _FluxPercolL3->reset();
  _FluxRecharge->reset();
  _BedrockLeakageFlux->reset();
  _FluxL1toLat->reset();
  _FluxL2toLat->reset();
  _FluxGWtoLat->reset();
  _FluxL1toChn->reset();
  _FluxL2toChn->reset();
  _FluxGWtoChn->reset();
  _FluxLattoL1->reset();
  _FluxLattoL2->reset();
  _FluxLattoGW->reset();
  if(ctrl.sw_trck){
    if(ctrl.sw_2H)
      trck.resetFd2HLat();
    if(ctrl.sw_18O)
      trck.resetFd18OLat();
    if(ctrl.sw_Age)
      trck.resetFAgeLat();
  }
  // Set EvapS to zero before looping over baresoil/understory
  _EvaporationS_all->reset();

  if(ctrl.sw_trck)
    trck.OutletVals(ctrl, 0, 0, 0);

#pragma omp parallel default(shared) \
  private(r, c, ra, rs, Ts, Tsold, Tdold, LAI, BeersK, Temp_can, emis_can, \
	  evap, infcap, accinf, theta, theta2, theta3, ponding,leak,  \
	  gw, za, z0u, zdu, z0o, zdo, wind, treeheight,			\
	  nr, le, sens, grndh, snowh, mltht, dh_snow, p, etc,		\
	  d1, d2, d3, fc, reinf, d)
  {
    //thre = omp_get_num_threads();
    //#pragma omp single
    //printf("\nnum threads %d: ", thre);
#pragma omp for nowait
    for (unsigned int j = 0; j < _vSortedGrid.cells.size() ; j++)
      {
	r = _vSortedGrid.cells[j].row;
	c = _vSortedGrid.cells[j].col;
	d = _vSortedGrid.cells[j].dir;	
	
	wind = atm.getWindSpeed()->matrix[r][c];
	
	theta = _soilmoist1->matrix[r][c]; //soil moisture at time t
	theta2 = _soilmoist2->matrix[r][c];
	theta3 = _soilmoist3->matrix[r][c];
	ponding = _ponding->matrix[r][c]; //surface ponding at time t
	gw = _GrndWater->matrix[r][c]; //gravity water at time t
	leak = 0;
	reinf = 0;

	fc = _fieldcapL3->matrix[r][c];
	d1 = _depth_layer1->matrix[r][c];
	d2 = _depth_layer2->matrix[r][c];
	d3 = _soildepth->matrix[r][c] - d1 - d2;
	
	nr = 0;
	le = 0;
	sens = 0;
	grndh = 0;
	snowh = 0;
	mltht = 0;
	Ts = _Temp_s->matrix[r][c];
	Tsold = 0;
	Tdold = 0;
	dh_snow = 0;
	BeersK = 0;

	// Infiltration + percolation if exceeds porosity

        // Tracking
        if (ctrl.sw_trck){
          _FluxSrftoL1->matrix[r][c] = 0;
          _FluxL1toL2->matrix[r][c] = 0;
          _FluxL2toL3->matrix[r][c] = 0;
        }
        Infilt_GreenAmpt(ctrl, infcap, accinf, theta, theta2, theta3, ponding, gw, dt, r, c);
        // Percolation if exceeding field capacity (L1 and L2), 
        // goes to GW in L3 (and bedrock leakage if activated)
        SoilWaterRedistribution(ctrl, accinf, theta, theta2, theta3, ponding, gw, leak, dt, r, c);
        // Tracking
        if(ctrl.sw_trck)
          trck.MixingV_down(*this, ctrl, d1, d2, d3, fc, r, c, 0);
        // Update global objects
        _ponding->matrix[r][c] = ponding;                       // [m]
        _GravityWater->matrix[r][c] = gw;                       // [m]
        _GrndWater->matrix[r][c] = gw;                          // [m]
        _BedrockLeakageFlux->matrix[r][c] = leak;               // [m.s-1]

	// Calculates the soil moisture profile to derive equivalent water table depth
	if(ctrl.Rep_WaterTableDepth == 1 || ctrl.RepTs_WaterTableDepth == 1)
	  CalcSoilMoistureProfile(atm, ctrl, getSoilMoist_av()->matrix[r][c], r,c);

	// Tracking
	if(ctrl.sw_trck)
	  _soilmoist1->matrix[r][c] = theta; 

	nsp = fForest->getNumSpecies();
	treeheight = 0;
	
	for(UINT4 s = 0; s < nsp ; s++)
	  {
	    p = fForest->getPropSpecies(s, r, c);
	    if(p == 0)
	      continue;//if no species j present, continue
	    
	    if(s == nsp -1){ //for bare soil, water reaching the ground is pp times its proportion of the cell
	      LAI = 0;
	      emis_can = 0;
	      Temp_can = 0;
	      za = _random_roughness->matrix[r][c] + 2;
	      z0u = max<REAL8>(0.000005,_random_roughness->matrix[r][c] * 0.1);
	      zdu = _random_roughness->matrix[r][c] * 0.7;
	      z0o = 0; //no overstory
	      zdo = 0;
	    }
	    else
	      {
		treeheight = max<REAL8>(0.01,fForest->getTreeHeight(s, r, c)); //equations only apply to 40% of the tree as per Campbell and Norman 1998
		LAI = fForest->getLAISpecies(s, r, c);
		BeersK = fForest->getBeersCoeff(s, r, c);
		Temp_can = fForest->getCanopyTemp(s, r, c);
		emis_can = fForest->getCanopyEmissivity(s, r, c);
		za = treeheight + 2;
		z0o = powl(treeheight,1.19)*0.057544;     //treeheight > 1 ? 0.1 : treeheight * 0.1;
		zdo = powl(treeheight,0.98)*0.707946; //treeheight > 1 ? 0.1 : treeheight * 0.7;
		zdu = min<double>(_random_roughness->matrix[r][c], zdo * 0.1);//min<double>(treeheight * 0.1, zdo * 0.1);
		z0u = 0.1*zdu/0.7;
		
	      }
	    
	    
	    ra = CalcAerodynResist(wind, za, z0u, zdu, z0o, zdo, treeheight,
				   LAI, Ts, atm.getTemperature()->matrix[r][c], ctrl.toggle_ra,
				   true);
	    
	    //rs = CalcSoilResist(_soilmoist1->matrix[r][c], r, c, ctrl.toggle_rs);
	    rs = CalcSoilResist(theta, r, c, ctrl.toggle_rs);
	    //rs =  1/max<double>( 0.0000000000001, ExfiltrationCapacity(theta, dt, r, c) );
	    
	    //theta_old = theta;

	    SolveSurfaceEnergyBalance(atm, ctrl, trck, ra, rs, 0.0, BeersK, LAI,
				      emis_can, Temp_can, nr, le, sens, grndh, snowh, mltht,
				      Tsold, evap, ponding, theta, Ts, Tdold, p, r, c, s);

	    
	    _Evaporation->matrix[r][c] += evap; //evaporation at t=t+1 (weighted by p)
	    _EvaporationS_all->matrix[r][c] += evap; //soil evaporation at t=t+1 (weighted by p)
	    // individual component of Esoil and ET (below vegetation only, de-weighted!)
	    if(s != nsp -1){
	      fForest->setEsoilSpecies(s, r, c, evap/p);
	      etc = fForest->getEvapoTransp(s, r, c);
	      fForest->setETSpecies(s, r, c, etc + evap/p);
	    }
	    
	    
	  }//for over the species

	// Update soil moisture objects
	_soilmoist1->matrix[r][c] = theta; //soil moisture at t=t+1
	_soilmoist2->matrix[r][c] = theta2;
	_soilmoist3->matrix[r][c] = theta3;
	
	_Rn->matrix[r][c] = nr;
	_latheat->matrix[r][c] = le;
	_latheat_sum->matrix[r][c] += le;
	_sensheat->matrix[r][c] = sens;
	_grndheat->matrix[r][c] = grndh;
	_snwheat->matrix[r][c] = snowh;
	_Temp_s_old->matrix[r][c] = Tsold;
	_Temp_s->matrix[r][c] = Tsold; //
	
	_Temp_d->matrix[r][c] = Tdold;
	       
	// Update surface pool
	_ponding->matrix[r][c] += SnowOutput(atm, ctrl, trck, mltht, r, c);
	// Back up before routing
	_GrndWater_old->matrix[r][c] = _GrndWater->matrix[r][c];
	
      }//for
  }//end omp parallel block
  
  return EXIT_SUCCESS;
}
