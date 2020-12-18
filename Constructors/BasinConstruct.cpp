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
 * BasinConstruct.cpp
 *
 *  Created on: Oct 9, 2009
 *      Author: Marco Maneta
 */
#include <unistd.h>
#include  <fstream>
#include "Basin.h"
#include <errno.h>
#include <stdio.h>
#include <string.h>


Basin::Basin(Control &ctrl)
{

  try{

    //Read the base map and writes the dimensions of the grid
    _DEM = new grid(ctrl.path_BasinFolder + ctrl.fn_dem, ctrl.MapType);
    _NRows = _DEM->r;
    _NCols = _DEM->c;
    _dx = _DEM->dx;

    _ldd = new grid(ctrl.path_BasinFolder + ctrl.fn_ldd, ctrl.MapType);

    printf("Checking if file %s exists...\n", (ctrl.path_BasinFolder + ctrl.fn_dem + ".serialized.svf").c_str());
    /*
     * Checks if there is a _vSordtedGrid object with the correct name in the spatial folder
     */
    if (access((ctrl.path_BasinFolder + ctrl.fn_dem + ".serialized.svf").c_str(), F_OK) != -1) {
      printf("File Found!. Loading object...\n");
      loadSortedGrid(_vSortedGrid, (ctrl.path_BasinFolder + ctrl.fn_dem + ".serialized.svf").c_str());
    }
    else{
      printf("File not found!. Initializing and sorting grid...\n");
      printf("WARNING -- if the progress bar stalls for too long, please consider checking your DEM map: it should\n");
      printf("contain a buffer of at least 1 cell of no-data (mv) around the valid domain (see documentation) --\n");
      /*sorts the basin with data cells according
       * to the ldd after _DEM and _ldd have been created*/
      _vSortedGrid = Basin::SortGridLDD();

      printf("Sorting done. Saving serialized sorted grid object for subsequent runs...\n");

      saveSortedGrid(_vSortedGrid, (ctrl.path_BasinFolder + ctrl.fn_dem + ".serialized.svf").c_str());
    }

    fForest = new Forest(ctrl); //constructs the Forest object

    /*basin parameters and properties*/
    _slope = new grid(ctrl.path_BasinFolder + ctrl.fn_slope, ctrl.MapType);

    _Ksat0 = new grid(ctrl.path_BasinFolder + ctrl.fn_Ksat0, ctrl.MapType);
    _kKsat = NULL;
    if(ctrl.sw_expKsat)
      _kKsat = new grid(ctrl.path_BasinFolder + ctrl.fn_kKsat, ctrl.MapType);
    
    _KvKs = new grid(ctrl.path_BasinFolder + ctrl.fn_kvkh, ctrl.MapType);
    _random_roughness = new grid(ctrl.path_BasinFolder + ctrl.fn_randrough, ctrl.MapType);

    _porosity0 = new grid(ctrl.path_BasinFolder + ctrl.fn_poros0, ctrl.MapType);
    _kporos = NULL;
    if(ctrl.sw_expPoros)
      _kporos = new grid(ctrl.path_BasinFolder + ctrl.fn_kporos, ctrl.MapType);
    
    _psi_ae  = new grid(ctrl.path_BasinFolder + ctrl.fn_psi_ae, ctrl.MapType);
    _BClambda = new grid(ctrl.path_BasinFolder + ctrl.fn_BClambda, ctrl.MapType);
    _theta_rL1 = new grid(ctrl.path_BasinFolder + ctrl.fn_theta_r, ctrl.MapType);
    _theta_rL2 = new grid(ctrl.path_BasinFolder + ctrl.fn_theta_r, ctrl.MapType);
    _theta_rL3 = new grid(ctrl.path_BasinFolder + ctrl.fn_theta_r, ctrl.MapType);
    _soildepth = new grid(ctrl.path_BasinFolder + ctrl.fn_soildepth, ctrl.MapType);
    _depth_layer1 = new grid(ctrl.path_BasinFolder + ctrl.fn_depth_layer1, ctrl.MapType);
    _depth_layer2 = new grid(ctrl.path_BasinFolder + ctrl.fn_depth_layer2, ctrl.MapType);
    _bedrock_leak = new grid(ctrl.path_BasinFolder + ctrl.fn_bedrock_leak, ctrl.MapType);
    _paramWc = new grid(ctrl.path_BasinFolder + ctrl.fn_paramWc, ctrl.MapType);
    _paramWp = new grid(ctrl.path_BasinFolder + ctrl.fn_paramWp, ctrl.MapType);
    _meltCoeff = new grid(ctrl.path_BasinFolder + ctrl.fn_snowCf, ctrl.MapType);

    /*state variables for depth-dependent soil parameters*/
    if(ctrl.sw_ddSoilPar) {
      _KsatL1 = new grid(ctrl.path_BasinFolder + ctrl.fn_Ksat0, ctrl.MapType);
      _KsatL2 = new grid(ctrl.path_BasinFolder + ctrl.fn_Ksat2, ctrl.MapType);
      _KsatL3 = new grid(ctrl.path_BasinFolder + ctrl.fn_Ksat3, ctrl.MapType);
      _KvKsL1 = new grid(ctrl.path_BasinFolder + ctrl.fn_kvkh, ctrl.MapType);
      _KvKsL2 = new grid(ctrl.path_BasinFolder + ctrl.fn_kvkh2, ctrl.MapType);
      _KvKsL3 = new grid(ctrl.path_BasinFolder + ctrl.fn_kvkh3, ctrl.MapType);
      _porosityL1 = new grid(ctrl.path_BasinFolder + ctrl.fn_poros0, ctrl.MapType);
      _porosityL2 = new grid(ctrl.path_BasinFolder + ctrl.fn_poros2, ctrl.MapType);
      _porosityL3 = new grid(ctrl.path_BasinFolder + ctrl.fn_poros3, ctrl.MapType);
      _psi_aeL1= new grid(ctrl.path_BasinFolder + ctrl.fn_psi_ae, ctrl.MapType);
      _psi_aeL2 = new grid(ctrl.path_BasinFolder + ctrl.fn_psi_ae2, ctrl.MapType);
      _psi_aeL3 = new grid(ctrl.path_BasinFolder + ctrl.fn_psi_ae3, ctrl.MapType);
      _BClambdaL1 = new grid(ctrl.path_BasinFolder + ctrl.fn_BClambda, ctrl.MapType);
      _BClambdaL2 = new grid(ctrl.path_BasinFolder + ctrl.fn_BClambda2, ctrl.MapType);
      _BClambdaL3 = new grid(ctrl.path_BasinFolder + ctrl.fn_BClambda3, ctrl.MapType);
    } else {
      _KsatL1 = new grid(*_DEM);
      _KsatL2 = new grid(*_DEM);
      _KsatL3 = new grid(*_DEM);
      _KvKsL1 = new grid(*_DEM);
      _KvKsL2 = new grid(*_DEM);
      _KvKsL3 = new grid(*_DEM);
      _porosityL1 = new grid(*_DEM);
      _porosityL2 = new grid(*_DEM);
      _porosityL3 = new grid(*_DEM);
      _psi_aeL1= new grid(*_DEM);
      _psi_aeL2 = new grid(*_DEM);
      _psi_aeL3 = new grid(*_DEM);      
      _BClambdaL1 = new grid(*_DEM);
      _BClambdaL2 = new grid(*_DEM);
      _BClambdaL3 = new grid(*_DEM);
    }

    //_rootfrac1 = new grid(ctrl.path_BasinFolder + ctrl.fn_root_fraction_lay1, ctrl.MapType);
    //_rootfrac2 = new grid(ctrl.path_BasinFolder + ctrl.fn_root_fraction_lay2, ctrl.MapType);
    //_Kroot = new grid(ctrl.path_BasinFolder + ctrl.fn_Kroot, ctrl.MapType);

    _snow = new grid(ctrl.path_BasinFolder + ctrl.fn_swe, ctrl.MapType);
    _snow_old = new grid(ctrl.path_BasinFolder + ctrl.fn_swe, ctrl.MapType);
    _albedo = new grid(ctrl.path_BasinFolder + ctrl.fn_albedo, ctrl.MapType);
    _emiss_surf = new grid(ctrl.path_BasinFolder + ctrl.fn_emiss, ctrl.MapType);
    _soil_dry_heatcap = new grid(ctrl.path_BasinFolder + ctrl.fn_soilheatcap, ctrl.MapType);
    _soil_dry_thermcond = new grid(ctrl.path_BasinFolder + ctrl.fn_soilthermcond, ctrl.MapType);
    _dampdepth = new grid(ctrl.path_BasinFolder + ctrl.fn_dampdepth, ctrl.MapType);
    _Temp_d = new grid(ctrl.path_BasinFolder + ctrl.fn_tempdamp, ctrl.MapType);

    _channelwidth = new grid(ctrl.path_BasinFolder + ctrl.fn_chwidth, ctrl.MapType);
    _chGWparam = new grid(ctrl.path_BasinFolder + ctrl.fn_chgwparam, ctrl.MapType);
    _Manningn = new grid(ctrl.path_BasinFolder + ctrl.fn_chmanningn, ctrl.MapType);

    /*state variables initialized with user map*/
    _soilmoist1 = new grid(ctrl.path_BasinFolder + ctrl.fn_soilmoist, ctrl.MapType);  //soil moisture volumetric
    _soilmoist2 = new grid(ctrl.path_BasinFolder + ctrl.fn_soilmoist2, ctrl.MapType);  //soil moisture 2nd layer volumetric
    _soilmoist3 = new grid(ctrl.path_BasinFolder + ctrl.fn_soilmoist3, ctrl.MapType);  //soil moisture 3rd layer volumetric
    _Temp_s_old = new grid(ctrl.path_BasinFolder + ctrl.fn_soiltemp, ctrl.MapType);  //initial soil temperature C
    _Disch_old =  new grid(ctrl.path_BasinFolder + ctrl.fn_streamflow, ctrl.MapType);  //initial soil temperature C

    /*state variables initialized with the base map*/
    _catcharea = new grid(*_DEM);
    _fieldcapL1 = new grid(*_DEM);
    _fieldcapL2 = new grid(*_DEM);
    _fieldcapL3 = new grid(*_DEM);
    _Rn = new grid(*_DEM);
    _Rn_sum = new grid(*_DEM);
    _latheat = new grid(*_DEM);
    _latheat_sum = new grid(*_DEM);
    _sensheat = new grid(*_DEM);
    _grndheat = new grid(*_DEM);
    _snwheat = new grid(*_DEM);
    _Temp_s = new grid(*_DEM);
    //_Temp_s_old = new grid(*_DEM);
    _CanopyStorage = new grid(*_DEM);
    //_rootfrac1 = new grid(*_DEM);
    //_rootfrac2 = new grid(*_DEM);

    _chan_store = new grid(*_DEM);
    _chan_store_old = new grid(*_DEM);
    
    _Layer1Outlet = new grid(*_DEM);
    _Layer2Outlet = new grid(*_DEM);
    _GWOutlet = new grid(*_DEM);

    // Depth at which 95% of roots are found, and 
    // contribution of each layer to the root zone
    // (both averaged over pixel-sharing species)
    _Zroot95 = new grid(*_DEM);
    _ProotzoneL1 = new grid(*_DEM);
    _ProotzoneL2 = new grid(*_DEM);
    _ProotzoneL3 = new grid(*_DEM);

    _IsSaturated = new grid(*_DEM);  //saturation map
    _soilmoist_av = new grid(*_DEM); //average volumetric soil moisture of the first 10 cm of the soil as calculated using a hydrstatic equilibrium moisture profile
    _soilmoist_12 = new grid(*_DEM); //average volumetric soil moisture of the upper two layers
    _ponding = new grid(*_DEM);
    _infilt_cap = new grid(*_DEM); //infilt cap m h-1
    _Evaporation = new grid(*_DEM); //actual evaporation in m s-1
    _BedrockLeakageFlux = new grid(*_DEM); //bedrock leakage flux in m s-1
    _SoilWaterDepth = new grid(*_DEM); //soil moisture depth m
    _SoilWaterDepthL1 = new grid(*_DEM); //soil moisture depth L1 m
    _SoilWaterDepthL2 = new grid(*_DEM); //soil moisture depth L2 m
    _SoilWaterDepthL3 = new grid(*_DEM); //soil moisture depth L3 (vadose) m
    _WaterTableDepth = new grid(*_DEM); //reconstructed WTD in meters
    _SoilSatDeficit = new grid(*_DEM); //soil moisture including water below and above field capacity
    _GravityWater = new grid(*_DEM); //soil water storage beyond
    _ponding_old = new grid(*_DEM); //surface storage at teh beginning of the time step
    _GrndWater_old = new grid(*_DEM); //groundwater storage at teh beginning of the time step
    _GrndWater = new grid(*_DEM); //groundwater storage at the end of the time step
    _Layer1upstreamBC = new grid(*_DEM); //Layer1 flux upstream boundary condition (ms-1)
    _Layer2upstreamBC = new grid(*_DEM); //Layer2 flux upstream boundary condition (ms-1)
    _GWupstreamBC = new grid(*_DEM); //gw flux upstream boundary condition (ms-1)
    _Disch_upstreamBC = new grid(*_DEM);
    _EvaporationS_all = new grid(*_DEM); //actual soil evaporation in m s-1
    _EvaporationI_all = new grid(*_DEM); //actual evaporation from summed interception in m s-1
    _Transpiration_all = new grid(*_DEM); //transpiration from summed in m s-1

    _FluxInfilt = new grid(*_DEM); // surface to L1 (summed over timestep)
    _FluxExfilt = new grid(*_DEM); // L1 to surface (return flow)
    _FluxPercolL2 = new grid(*_DEM); // L1 to L2 (summed over timestep)
    _FluxL2toL1 = new grid(*_DEM); // capillary + return flow, L2 to L1
    _FluxPercolL3 = new grid(*_DEM); // L2 to L3 (summed over timestep)
    _FluxRecharge = new grid(*_DEM); // Recharge to GW
    _FluxL3toL2 = new grid(*_DEM); // return from L3 to L2
    _FluxLattoChn = new grid(*_DEM); // channel inflow
    _FluxLattoSrf = new grid(*_DEM); // surface run-on (excluding streamflow)
    _FluxLattoL1 = new grid(*_DEM); // layer 1 lateral outflow
    _FluxLattoL2 = new grid(*_DEM); // layer 2 lateral outflow
    _FluxLattoGW = new grid(*_DEM); // groundwater lateral outflow
    _FluxChntoLat = new grid(*_DEM); // surface run-off (only streamflow)
    _FluxSrftoLat = new grid(*_DEM); // surface run-off (excluding streamflow)
    _FluxL1toLat = new grid(*_DEM); // layer 1 lateral inflow	
    _FluxL2toLat = new grid(*_DEM); // layer 2 lateral inflow	
    _FluxGWtoLat = new grid(*_DEM); // groundwater lateral inflow	
    _FluxL1toChn = new grid(*_DEM); // intra-cell groundwater to channel
    _FluxL2toChn = new grid(*_DEM); // intra-cell groundwater to channel
    _FluxGWtoChn = new grid(*_DEM); // intra-cell groundwater to channel
    _FluxSrftoChn = new grid(*_DEM); // intra-cell ponding to channel

    // Accumulated fluxes
    _AccInfilt = new grid(*_DEM);
    _AccExfilt = new grid(*_DEM);
    _AccPercolL2 = new grid(*_DEM); // L1 to L2 (summed over timestep)
    _AccL2toL1 = new grid(*_DEM); // capillary + return flow, L2 to L1
    _AccPercolL3 = new grid(*_DEM); // L2 to L3 (summed over timestep)
    _AccRecharge = new grid(*_DEM); // Recharge to GW
    _AccLeak = new grid(*_DEM); //leakance from L3
    _AccEvaporationS = new grid(*_DEM);
    _AccTranspiL1 = new grid(*_DEM);
    _AccTranspiL2 = new grid(*_DEM);
    _AccTranspiL3 = new grid(*_DEM);
    _AccL3toL2 = new grid(*_DEM); // return from L3 to L2
    _AccLattoChn = new grid(*_DEM);
    _AccLattoSrf = new grid(*_DEM);
    _AccLattoGW = new grid(*_DEM);
    _AccChntoLat = new grid(*_DEM);
    _AccSrftoLat = new grid(*_DEM);
    _AccL1toLat = new grid(*_DEM);
    _AccL2toLat = new grid(*_DEM);
    _AccGWtoLat = new grid(*_DEM);
    _AccL1toChn = new grid(*_DEM); // groundwater to channel (accumlated)
    _AccL2toChn = new grid(*_DEM); // groundwater to channel (accumlated)
    _AccGWtoChn = new grid(*_DEM); // groundwater to channel (accumlated)
    _AccSrftoChn = new grid(*_DEM); // ponding to channel (accumulated)

    // Maps specifically for channel evaporation 
    _Temp_w = NULL;
    _FTemp_w = NULL;
    _chan_evap = new grid(*_DEM);// channel evaporation
    _chan_roughness = NULL;
    if(ctrl.sw_chan_evap) {
      _Temp_w = new grid(ctrl.path_BasinFolder + ctrl.fn_temp_w, ctrl.MapType);  // channel water temperature
      _FTemp_w = new grid(*_DEM);
      _chan_roughness = new grid(ctrl.path_BasinFolder + ctrl.fn_chanrough, ctrl.MapType);
    }


    // Tracking fluxes specifically for isotopes
    _psi_MW = NULL;
    _moist_MW1 = NULL ;
    _moist_MW2 = NULL ;
    _fracMW1 = NULL ;
    _fracMW2 = NULL;
    _fracMW12 = NULL;
    _FluxCnptoSrf = NULL;
    _FluxCnptoSnow = NULL;
    _FluxSnowtoSrf = NULL;
    _FluxSrftoL1 = NULL;
    _FluxL1toL2 = NULL;
    _FluxL2toL3 = NULL;
    //_FluxL2toGW = NULL;
    //_FluxL3toGW = NULL;
    //_FluxGWtoL2 = NULL;
    //_FluxGWtoL3 = NULL;
    _FluxLeak = NULL;
    _IncompAlpha = NULL;
		
    if(ctrl.sw_trck){
      _FluxCnptoSrf = new grid(*_DEM); // canopy/sky to surface
      _FluxCnptoSnow = new grid(*_DEM); // canopy/sky to snowpac
      _FluxSnowtoSrf = new grid(*_DEM); // snowpack to surface
      _FluxSrftoL1 = new grid(*_DEM); // surface to first layer
      _FluxL1toL2 = new grid(*_DEM); // percolation L1 to L2
      _FluxL2toL3 = new grid(*_DEM); // percolation L2 to L3
      //_FluxL2toGW = new grid(*_DEM); // recharge L3 to groundwater
      //_FluxL3toGW = new grid(*_DEM); // recharge L3 to groundwater
      //_FluxGWtoL2 = new grid(*_DEM); // return flow, groundwater to L2
      //_FluxGWtoL3 = new grid(*_DEM); // discharge, groundwater to L3
      _FluxLeak = new grid(*_DEM);  // groundwater to leakance
      if(ctrl.toggle_mix == 3)
	_IncompAlpha = new grid(ctrl.path_BasinFolder+ctrl.fn_IncompMix, ctrl.MapType);
      
      if(ctrl.sw_TPD){
	_psi_MW = new grid(ctrl.path_BasinFolder+ctrl.fn_psi_MW, ctrl.MapType);
	_moist_MW1 = new grid(*_DEM);
	_moist_MW2 = new grid(*_DEM);
	_fracMW1 = new grid(*_DEM);
	_fracMW2 = new grid(*_DEM);
	_fracMW12 = new grid(*_DEM);
      }
    }

    //    cout << errno << endl ;
    
    // reset errno just to make sure that any constructor error here are not
    // "leftovers" from earlier constrcutors' failures 
    errno = 0;

    try{
      //calculate the value of Ksat for each layer (integrated from expo profile)
      CalcKsatLayers(ctrl); 
      if(errno!=0){
	cout << "Error creating Ksat for each layer: " << endl;
	throw string("Ksat, soil depths, kKsat");
      }
      //calculate the value of porosity for each layer (integrated from expo profile)
      CalcPorosLayers(ctrl);
      if(errno!=0){
	cout << "Error creating Kporos for each layer: " << endl;
	throw string("porosm, kporos, soil depths");
      }
      //Partial check of maps mainly to make sure no no data is written within the valid domain
      CheckMaps(ctrl);
      if(errno!=0){
	cout << "Error creating maps: " << endl;
	throw string(" ");
      }
      //Fills-in the _catcharea map with the upstream catchment area of each cell
      CalcCatchArea();
      if(errno!=0){
	cout << "Error creating catchment area: " << strerror(errno) << endl;
	throw string(" ");
      }
      // Calculate the value of field capacity using the Brooks and Corey Formula
      CalcFieldCapacity(ctrl);	
      if(errno!=0){
	cout << "Error creating field capacity: " << endl;
	throw string("psi_ae, BClambda, poros, theta_r ");
      }
	
      //Reads initial streamflow map and populate _ponding variable with initial storage in stream
      CalcInitialStreamStorage();
      if(errno!=0){
	cout << "Error creating stream storage maps: " << endl;
	throw string("slope, chanmask, chanmann, chanwidth ");
      }

      // From the species root profile and layer depth --> root frac in each layer
      CalcRootDistrib();
      if(errno!=0){
	cout << "Error creating root distributions: " << endl;
	throw string("kroot, soil depths  ");
      }

      // From the specified tightly-bound / mobile tension threshold, the threshold moisture is calculated using the Brooks-Corey formula 
      if(ctrl.sw_trck and ctrl.sw_TPD){
	CalcTPDMoisture(ctrl);
	if(errno!=0){
	  cout << "Error creating the TPD moisture threshold: " << endl;
	  throw string(" psi_ae, psi_MW, bclambda, poros, theta_r  ");
	}
      }

    } catch (string e){
      cout << "Check the  " << e << " maps, error " << strerror(errno) << endl;
      throw;
    }

  }catch (std::bad_alloc &)
    { cerr << " Cleaning basin objects..." << "\n";
      if(_DEM)
	delete _DEM;
      if(_ldd)
	delete _ldd;
      if(_snow)
	delete _snow;
      if(_snow_old)
	delete _snow_old;
      if(_Rn)
	delete _Rn;
      if(_Rn_sum)
	delete _Rn_sum;
      if(_latheat)
	delete _latheat;
      if(_latheat_sum)
	delete _latheat_sum;
      if(_sensheat)
	delete _sensheat;
      if(_grndheat)
	delete _grndheat;
      if(_snwheat)
	delete _snwheat;
      if(_Temp_s)
	delete _Temp_s;
      if(_Temp_s_old)
	delete _Temp_s_old;
      if(_albedo)
	delete _albedo;
      if(_emiss_surf)
	delete _emiss_surf;
      if(_soil_dry_heatcap)
	delete _soil_dry_heatcap;
      if(_soil_dry_thermcond)
	delete _soil_dry_thermcond;
      if(_dampdepth)
	delete _dampdepth;
      if(_Temp_d)
	delete _Temp_d;
      if(_ponding)
	delete _ponding;
      if(_Ksat0)
	delete _Ksat0;
      if(_kKsat)
	delete _kKsat;
      if(_KsatL1)
	delete _KsatL1;
      if(_KsatL2)
	delete _KsatL2;
      if(_KsatL3)
	delete _KsatL3;
      if(_KvKs)
	delete _KvKs;
      if(_KvKsL1)
	delete _KvKsL1;
      if(_KvKsL2)
	delete _KvKsL2;
      if(_KvKsL3)
	delete _KvKsL3;
      if(_random_roughness)
	delete _random_roughness;
      if(_slope)
	delete _slope;
      if(_porosity0)
	delete _porosity0;
      if(_kporos)
	delete _kporos;
      if(_porosityL1)
	delete _porosityL1;
      if(_porosityL2)
	delete _porosityL2;
      if(_porosityL3)
	delete _porosityL3;
      if(_psi_ae)
	delete _psi_ae;
      if(_psi_aeL1)
	delete _psi_aeL1;
      if(_psi_aeL2)
	delete _psi_aeL2;
      if(_psi_aeL3)
	delete _psi_aeL3;
      if(_BClambda)
	delete _BClambda;
      if(_BClambdaL1)
	delete _BClambdaL1;
      if(_BClambdaL2)
	delete _BClambdaL2;
      if(_BClambdaL3)
	delete _BClambdaL3;
      if(_theta_rL1)
	delete _theta_rL1;
      if(_theta_rL2)
	delete _theta_rL2;
      if(_theta_rL3)
	delete _theta_rL3;
      if(_infilt_cap)
	delete _infilt_cap;
      if(_soilmoist1)
	delete _soilmoist1;
      if(_soildepth)
	delete _soildepth;
      if(_depth_layer1)
	delete _depth_layer1;
      if(_depth_layer2)
	delete _depth_layer2;

      if(_chan_store)
	delete _chan_store;
      if(_chan_store_old)
	delete _chan_store_old;

      if(_Layer1Outlet)
	delete _Layer1Outlet;
      if(_Layer2Outlet)
	delete _Layer2Outlet;
      if(_GWOutlet)
	delete _GWOutlet;

      if(_Temp_w)
	delete _Temp_w;
      if(_FTemp_w)
	delete _FTemp_w;
      if(_chan_evap)
	delete _chan_evap;
      if(_chan_roughness)
	delete _chan_roughness;

      if(_Zroot95)
	delete _Zroot95;
      if(_ProotzoneL1)
	delete _ProotzoneL1;
      if(_ProotzoneL2)
	delete _ProotzoneL2;
      if(_ProotzoneL3)
	delete _ProotzoneL3;
      if(_fieldcapL1)
	delete _fieldcapL1;
      if(_fieldcapL2)
	delete _fieldcapL2;
      if(_fieldcapL3)
	delete _fieldcapL3;
      if(_paramWc)
	delete _paramWc;
      if(_paramWp)
	delete _paramWp;
      if(_meltCoeff)
	delete _meltCoeff;
      if(_Evaporation)
	delete _Evaporation;
      if(_BedrockLeakageFlux)
	delete _BedrockLeakageFlux;
      if(_SoilWaterDepth)
	delete _SoilWaterDepth;
      if(_SoilWaterDepthL1)
	delete _SoilWaterDepthL1;
      if(_SoilWaterDepthL2)
	delete _SoilWaterDepthL2;
      if(_SoilWaterDepthL3)
	delete _SoilWaterDepthL3;
      if(_WaterTableDepth)
	delete _WaterTableDepth;
      if(_SoilSatDeficit)
	delete _SoilSatDeficit;
      if(_CanopyStorage)
	delete _CanopyStorage;
      if(_Disch_old)
	delete _Disch_old;
      if(_Disch_upstreamBC)
	delete _Disch_upstreamBC;
      if(_catcharea)
	delete _catcharea;
      if(_GravityWater)
	delete _GravityWater;
      if(_GrndWater_old)
	delete _GrndWater_old;
      if(_GrndWater)
	delete _GrndWater;
      if(_Layer1upstreamBC)
	delete _Layer1upstreamBC;
      if(_Layer2upstreamBC)
	delete _Layer2upstreamBC;
      if(_GWupstreamBC)
	delete _GWupstreamBC;
      if(_channelwidth)
	delete _channelwidth;
      if(_chGWparam)
	delete _chGWparam;
      if(_Manningn)
	delete _Manningn;
      if(_soilmoist_av)
	delete _soilmoist_av;
      if(_soilmoist_12)
	delete _soilmoist_12;
      if(_soilmoist2)
	delete _soilmoist2;
      if(_soilmoist3)
	delete _soilmoist3;
      if(_bedrock_leak)
	delete _bedrock_leak;
      if(_IsSaturated)
	delete _IsSaturated;
      if(_EvaporationS_all)
	delete _EvaporationS_all;
      if(_EvaporationI_all)
	delete _EvaporationI_all;
      if(_Transpiration_all)
	delete _Transpiration_all;


      if(fForest)
	delete fForest;

      if(_ponding_old)
	delete _ponding_old;
      if(_FluxCnptoSrf)
	delete _FluxCnptoSrf;
      if(_FluxCnptoSnow)
	delete _FluxCnptoSnow;
      if(_FluxSnowtoSrf)
	delete _FluxSnowtoSrf;
      if(_FluxSrftoL1)
	delete _FluxSrftoL1;
      if(_FluxLeak)
	delete _FluxLeak;
      if(_FluxInfilt)
	delete _FluxInfilt;
      if(_FluxExfilt)
	delete _FluxExfilt;
      if(_FluxPercolL2)
	delete _FluxPercolL2;
      if(_FluxL2toL1)
	delete _FluxL2toL1;
      if(_FluxPercolL3)
	delete _FluxPercolL3;
      if(_FluxRecharge)
	delete _FluxRecharge;
      if(_FluxL3toL2)
	delete _FluxL3toL2;
      if(_FluxLattoSrf)
	delete _FluxLattoSrf;
      if(_FluxLattoL1)
	delete _FluxLattoL1;
      if(_FluxLattoL2)
	delete _FluxLattoL2;
      if(_FluxLattoGW)
	delete _FluxLattoGW;
      if(_FluxLattoChn)
	delete _FluxLattoChn;
      if(_FluxChntoLat)
	delete _FluxChntoLat;
      if(_FluxSrftoLat)
	delete _FluxSrftoLat;
      if(_FluxL1toLat)
	delete _FluxL1toLat;
      if(_FluxL2toLat)
	delete _FluxL2toLat;
      if(_FluxGWtoLat)
	delete _FluxGWtoLat;
      if(_FluxL1toChn)
	delete _FluxL1toChn;
      if(_FluxL2toChn)
	delete _FluxL2toChn;
      if(_FluxGWtoChn)
	delete _FluxGWtoChn;
      if(_FluxSrftoChn)
	delete _FluxSrftoChn;

      if(_AccInfilt)
	delete _AccInfilt;
      if(_AccExfilt)
	delete _AccExfilt;
      if(_AccPercolL2)
	delete _AccPercolL2;
      if(_AccL2toL1)
	delete _AccL2toL1;
      if(_AccPercolL3)
	delete _AccPercolL3;
      if(_AccRecharge)
	delete _AccRecharge;
      if(_AccLeak)
	delete _AccLeak;
      if(_AccL3toL2)
	delete _AccL3toL2;
      if(_AccEvaporationS)
	delete _AccEvaporationS;
      if(_AccTranspiL1)
	delete _AccTranspiL1;
      if(_AccTranspiL2)
	delete _AccTranspiL2;
      if(_AccTranspiL3)
	delete _AccTranspiL3;
      if(_AccLattoGW)
	delete _AccLattoGW;
      if(_AccLattoSrf)
	delete _AccLattoSrf;
      if(_AccLattoChn)
	delete _AccLattoChn;
      if(_AccChntoLat)
	delete _AccChntoLat;
      if(_AccSrftoLat)
	delete _AccSrftoLat;
      if(_AccL1toLat)
	delete _AccL1toLat;
      if(_AccL2toLat)
	delete _AccL2toLat;
      if(_AccGWtoLat)
	delete _AccGWtoLat;
      if(_AccL1toChn)
	delete _AccL1toChn;
      if(_AccL2toChn)
	delete _AccL2toChn;
      if(_AccGWtoChn)
	delete _AccGWtoChn;
      if(_AccSrftoChn)
	delete _AccSrftoChn;

      // Tracking
      if(_psi_MW)
	delete _psi_MW;
      if(_moist_MW1)
	delete _moist_MW1;
      if(_moist_MW2)
	delete _moist_MW2;
      if(_fracMW1)
	delete _fracMW1;
      if(_fracMW2)
	delete _fracMW2;
      if(_fracMW12)
	delete _fracMW12;
      if(_FluxL1toL2)
	delete _FluxL1toL2;
      if(_FluxL2toL3)
	delete _FluxL2toL3;
      if(_IncompAlpha)
	delete _IncompAlpha;
      //if(_FluxL2toGW)
      //		delete _FluxL2toGW;
      //if(_FluxL3toGW)
      //	delete _FluxL3toGW;
      //if(_FluxGWtoL2)
      //	delete _FluxGWtoL2;
      //if(_FluxGWtoL3)
      //delete _FluxGWtoL3;


      throw;
    }

  /*	for(unsigned int r = 0; r <  _NRows; r++)//remove!!!!!
	for(unsigned int c = 0; c <  _NCols; c++)
	{
	_Temp_s->matrix[r][c] = 2;
	_Temp_s_old->matrix[r][c] = 1;
	}*/
}
