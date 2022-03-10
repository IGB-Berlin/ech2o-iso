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
#include <armadillo>
#include "Basin.h"
#include <errno.h>
#include <stdio.h>
#include <string.h>


Basin::Basin(Control &ctrl, Atmosphere &atm)
{

  try{

    // -------------------------------------------------------------------
    // Read the base map and writes the dimensions of the grid
    // -------------------------------------------------------------------
    _DEM = new grid(ctrl.path_BasinFolder + ctrl.fn_dem, ctrl.MapType);
    _NRows = _DEM->r;
    _NCols = _DEM->c;
    _dx = _DEM->dx;
    //geo coordinations for nc files
    _north = _DEM->north;
    _south = _DEM->south;
    _west = _DEM->west;
    _east = _DEM->east;
    _nodata = _DEM->nodata;

    // -------------------------------------------------------------------
    // If boundary conditions are to be est - read them in
    // -------------------------------------------------------------------
    _BCsurface     = NULL;
    _BCgroundwater = NULL;
    _BCdeepgwtr    = NULL;
    if(ctrl.sw_BC){
      errno = 0; //reset the error
      _BCsurface     = new grid(*_DEM); //incoming surface water
      _BCgroundwater = new grid(*_DEM); //incoming groundwater water
      if(ctrl.sw_deepGW)
	_BCdeepgwtr = new grid(*_DEM); //incoming groundwater water

      *_BCsurface    = *_DEM;
      *_BCgroundwater= *_DEM;
      if(ctrl.sw_deepGW)
	*_BCdeepgwtr = *_DEM;
	
      try{
	ifBCsurface.open((ctrl.path_ClimMapsFolder + ctrl.fn_BCsurface).c_str(),ios::binary);
	if(errno!=0) throw ctrl.fn_BCsurface;
	ifBCgroundwater.open((ctrl.path_ClimMapsFolder + ctrl.fn_BCgroundwater).c_str(),ios::binary);
	if(errno!=0) throw ctrl.fn_BCgroundwater;
	if(ctrl.sw_deepGW){
	  ifBCdeepgwtr.open((ctrl.path_ClimMapsFolder + ctrl.fn_BCdeepgwtr).c_str(),ios::binary);
	  if(errno!=0) throw ctrl.fn_BCdeepgwtr;
	}
      } catch (string e){
	cout << "Dang!!: cannot find/read the " << e << " file: error " << strerror(errno) << endl;
	throw;
      }

      try{
	if(InitiateBCMap(ifBCsurface, *_BCsurface, atm) != atm.getSsortedGridTotalCellNumber())
	  throw string("BC surface");
	if(InitiateBCMap(ifBCgroundwater, *_BCgroundwater, atm) != atm.getSsortedGridTotalCellNumber())
	  throw string("BC groundwater");
	if(ctrl.sw_deepGW){
	  if(InitiateBCMap(ifBCdeepgwtr, *_BCdeepgwtr, atm) != atm.getSsortedGridTotalCellNumber())
	    throw string("BC deep groundwater");
	}
      } catch (string e) {
	cout << "Error: some sections of the domain were not filled with " << e << " data." << endl;
	cout << "Please verify that all the boundary zones in the map are presented in the binary" << endl;
	cout << "and that the n boundary zones present are the first n zones in the file" << endl;
      }    } //toggle BC
    // -------------------------------------------------------------------    


    // ---------------------------------------------------------------------------------------------------------------
    /*Catchment map check/creation of sorted grid*/
    // ---------------------------------------------------------------------------------------------------------------        
    _ldd = new grid(ctrl.path_BasinFolder + ctrl.fn_ldd, ctrl.MapType);

    printf("Checking if file %s exists...\n", (ctrl.path_BasinFolder + ctrl.fn_dem + ".serialized.svf").c_str());
    /* Checks if there is a _vSordtedGrid object with the correct name in the spatial folder */
    if (access((ctrl.path_BasinFolder + ctrl.fn_dem + ".serialized.svf").c_str(), F_OK) != -1) {
      printf("File Found!. Loading object...\n");
      loadSortedGrid(_vSortedGrid, (ctrl.path_BasinFolder + ctrl.fn_dem + ".serialized.svf").c_str());
    }
    else{
      printf("File not found!. Initializing and sorting grid...\n");
      printf("WARNING -- if the progress bar stalls for too long, please consider checking your DEM map: it should\n");
      printf("contain a buffer of at least 1 cell of no-data (mv) around the valid domain (see documentation) --\n");
      /*sorts the basin with data cells according to the ldd after _DEM and _ldd have been created*/
      _vSortedGrid = Basin::SortGridLDD();
      printf("Sorting done. Saving serialized sorted grid object for subsequent runs...\n");
      saveSortedGrid(_vSortedGrid, (ctrl.path_BasinFolder + ctrl.fn_dem + ".serialized.svf").c_str());
    }

    fForest = new Forest(ctrl); //constructs the Forest object

    // ---------------------------------------------------------------------------------------------------------------
    /*Basin parameters and properties*/
    // ---------------------------------------------------------------------------------------------------------------    
    _slope = new grid(ctrl.path_BasinFolder + ctrl.fn_slope, ctrl.MapType);           //Basin Slope
    _ttarea = new grid(ctrl.path_BasinFolder + ctrl.fn_ttarea, ctrl.MapType);         //Cell % within catchment
    _fImperv = new grid(ctrl.path_BasinFolder + ctrl.fn_fimperv, ctrl.MapType);       //Percent impervious area of cell
    
    _Ksat0 = new grid(ctrl.path_BasinFolder + ctrl.fn_Ksat0, ctrl.MapType);           //Horizontal conductivity near surface
    _KsatTopSoil = new grid(ctrl.path_BasinFolder + ctrl.fn_ksat_skin, ctrl.MapType); //Skin horizontal conductivity    
    _kKsat = NULL;                                                                    //Exp parameter dec. conductivity with depth
    if(ctrl.sw_expKsat)
      _kKsat = new grid(ctrl.path_BasinFolder + ctrl.fn_kKsat, ctrl.MapType);
    
    _KvKs = new grid(ctrl.path_BasinFolder + ctrl.fn_kvkh, ctrl.MapType);             //Anisotropy ratio for conductivity
    _random_roughness = new grid(ctrl.path_BasinFolder + ctrl.fn_randrough, ctrl.MapType);//Random surface roughness
    _porosity0 = new grid(ctrl.path_BasinFolder + ctrl.fn_poros0, ctrl.MapType);      //Soil porosity near surface
    _kporos = NULL;                                                                   //Exp parameter dec. porosity with depth
    if(ctrl.sw_expPoros)
      _kporos = new grid(ctrl.path_BasinFolder + ctrl.fn_kporos, ctrl.MapType);
    
    _psi_ae  = new grid(ctrl.path_BasinFolder + ctrl.fn_psi_ae, ctrl.MapType);        //Air-entry pressure [m]
    _BClambda = new grid(ctrl.path_BasinFolder + ctrl.fn_BClambda, ctrl.MapType);     //Brooks-Corey lambda
    _theta_rL1 = new grid(ctrl.path_BasinFolder + ctrl.fn_theta_r, ctrl.MapType);     //Residual moisture layer 1
    _theta_rL2 = new grid(ctrl.path_BasinFolder + ctrl.fn_theta_r, ctrl.MapType);     //Residual moisture layer 2
    _theta_rL3 = new grid(ctrl.path_BasinFolder + ctrl.fn_theta_r, ctrl.MapType);     //Residual moisture layer 3
    _soildepth = new grid(ctrl.path_BasinFolder + ctrl.fn_soildepth, ctrl.MapType);   //Soil depth (depth of layer 3) [m]
    _depth_layer1 = new grid(ctrl.path_BasinFolder + ctrl.fn_depth_layer1, ctrl.MapType);//Depth of soil layer 1 [m]
    _depth_layer2 = new grid(ctrl.path_BasinFolder + ctrl.fn_depth_layer2, ctrl.MapType);//Depth of soil layer 2 [m]
    _bedrock_leak = new grid(ctrl.path_BasinFolder + ctrl.fn_bedrock_leak, ctrl.MapType);//Bedrock leakance (frac of Ksat in L3)
    _paramWc = new grid(ctrl.path_BasinFolder + ctrl.fn_paramWc, ctrl.MapType);       //Veg water use parameter 
    _paramWp = new grid(ctrl.path_BasinFolder + ctrl.fn_paramWp, ctrl.MapType);       //Veg water use parameter
    // ---------------------------------------------------------------------------------------------------------------    
    /*state variables for depth-dependent soil parameters*/
    // ---------------------------------------------------------------------------------------------------------------    
    if(ctrl.toggle_soil_prop == 2) {
      _KsatL1 = new grid(ctrl.path_BasinFolder + ctrl.fn_Ksat0, ctrl.MapType);        //Horizontal conductivity L1
      _KsatL2 = new grid(ctrl.path_BasinFolder + ctrl.fn_Ksat2, ctrl.MapType);        //Horizontal conductivity L2
      _KsatL3 = new grid(ctrl.path_BasinFolder + ctrl.fn_Ksat3, ctrl.MapType);        //Horizontal conductivity L3
      _KvKsL1 = new grid(ctrl.path_BasinFolder + ctrl.fn_kvkh,  ctrl.MapType);        //Anisotropy ratio L1
      _KvKsL2 = new grid(ctrl.path_BasinFolder + ctrl.fn_kvkh2, ctrl.MapType);        //Anisotropy ratio L2
      _KvKsL3 = new grid(ctrl.path_BasinFolder + ctrl.fn_kvkh3, ctrl.MapType);        //Anisotropy ratio L3
      _porosityL1 = new grid(ctrl.path_BasinFolder + ctrl.fn_poros0, ctrl.MapType);   //Porosity L1
      _porosityL2 = new grid(ctrl.path_BasinFolder + ctrl.fn_poros2, ctrl.MapType);   //Porosity L2
      _porosityL3 = new grid(ctrl.path_BasinFolder + ctrl.fn_poros3, ctrl.MapType);   //Porosity L3
      _psi_aeL1 = new grid(ctrl.path_BasinFolder + ctrl.fn_psi_ae,  ctrl.MapType);    //Air-entry pressure L1
      _psi_aeL2 = new grid(ctrl.path_BasinFolder + ctrl.fn_psi_ae2, ctrl.MapType);    //Air-entry pressure L2
      _psi_aeL3 = new grid(ctrl.path_BasinFolder + ctrl.fn_psi_ae3, ctrl.MapType);    //Air-entry pressure L3
      _BClambdaL1 = new grid(ctrl.path_BasinFolder + ctrl.fn_BClambda,  ctrl.MapType);//Brooks-Corey lambda L1
      _BClambdaL2 = new grid(ctrl.path_BasinFolder + ctrl.fn_BClambda2, ctrl.MapType);//Brooks-Corey lambda L2
      _BClambdaL3 = new grid(ctrl.path_BasinFolder + ctrl.fn_BClambda3, ctrl.MapType);//Brooks-Corey lambda L3
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

    // ---------------------------------------------------------------------------------------------------------------    
    /*state variables for deep groundwater*/
    // ---------------------------------------------------------------------------------------------------------------    
    _DeepGW = NULL;                                                                  //Baseflow to chan water transfer parameter [m]
    _Hydrofrac_DeepGW = NULL;                                                        //Baseflow to chan water transfer parameter [0-1]
    _chDeepGWparam = NULL;                                                           //Deep GW storage for each grid cell
    _FluxDeepGWtoChn = NULL;                                                         //Deep baseflow from the deep GW storge
    _FluxLattoDeepGW = NULL;                                                         //Lateral deep GW input
    _FluxDeepGWtoLat = NULL;                                                         //Lateral deep GW output
    _AccDeepGWtoChn = NULL;                                                          //Accumulated baseflow export to channel
    _AccLattoDeepGW = NULL;                                                          //Accumulated baseflow input
    _AccDeepGWtoLat = NULL;                                                          //Accumulated baseflow output
    _DeepGWupstreamBC = NULL;
    if(ctrl.sw_deepGW){
      _DeepGW = new grid(ctrl.path_BasinFolder + ctrl.fn_deepGW, ctrl.MapType);
      _Hydrofrac_DeepGW = new grid(ctrl.path_BasinFolder + ctrl.fn_hydro_deepGW, ctrl.MapType);
      _chDeepGWparam = new grid(ctrl.path_BasinFolder + ctrl.fn_chdeepgwparam, ctrl.MapType);
      _FluxDeepGWtoChn = new grid(*_DEM);
      _FluxLattoDeepGW = new grid(*_DEM);
      _FluxDeepGWtoLat = new grid(*_DEM);
      _AccDeepGWtoChn = new grid(*_DEM); 
      _AccLattoDeepGW = new grid(*_DEM); 
      _AccDeepGWtoLat  = new grid(*_DEM);
      _DeepGWupstreamBC = new grid(*_DEM);
    }
    
    // ---------------------------------------------------------------------------------------------------------------    
    /*state variables for energy balance*/
    // ---------------------------------------------------------------------------------------------------------------    
    _meltCoeff = new grid(ctrl.path_BasinFolder + ctrl.fn_snowCf, ctrl.MapType);      //Degree-day snowmelt coefficient
    _albedo = new grid(ctrl.path_BasinFolder + ctrl.fn_albedo, ctrl.MapType);         //Soil albedo
    _emiss_surf = new grid(ctrl.path_BasinFolder + ctrl.fn_emiss, ctrl.MapType);      //Soil emissivity
    _soil_dry_heatcap = new grid(ctrl.path_BasinFolder + ctrl.fn_soilheatcap, ctrl.MapType);//Soil heat capacity
    _soil_dry_thermcond = new grid(ctrl.path_BasinFolder + ctrl.fn_soilthermcond, ctrl.MapType);//Soil heat conductivity
    _dampdepth = new grid(ctrl.path_BasinFolder + ctrl.fn_dampdepth, ctrl.MapType);   //Damping depth
    _Temp_d = new grid(ctrl.path_BasinFolder + ctrl.fn_tempdamp, ctrl.MapType);       //Temperature at damping depth

    _Temp_L1 = new grid(ctrl.path_BasinFolder + ctrl.fn_tempdamp, ctrl.MapType);       //Temperature at damping depth
    _Temp_L2 = new grid(ctrl.path_BasinFolder + ctrl.fn_tempdamp, ctrl.MapType);       //Temperature at damping depth
    _Temp_L3 = new grid(ctrl.path_BasinFolder + ctrl.fn_tempdamp, ctrl.MapType);       //Temperature at damping depth

    // ---------------------------------------------------------------------------------------------------------------    
    /*state variables for channels*/
    // ---------------------------------------------------------------------------------------------------------------    
    _channelwidth = new grid(ctrl.path_BasinFolder + ctrl.fn_chwidth, ctrl.MapType);  //channel width [m]
    _channellength = new grid(ctrl.path_BasinFolder + ctrl.fn_chlength,ctrl.MapType); //channel length [m]
    _chGWparam = new grid(ctrl.path_BasinFolder + ctrl.fn_chgwparam, ctrl.MapType);   //GW to channel parameter
    _Manningn = new grid(ctrl.path_BasinFolder + ctrl.fn_chmanningn, ctrl.MapType);   //Mannings n
    _Temp_w = NULL;                                                                   //Channel water temp [oC]
    _FTemp_w = NULL;                                                                  //Flux*Water Temp
    _chan_roughness = NULL;                                                           //Channel roughness
    if(ctrl.toggle_chan_evap==1) {
      _Temp_w = new grid(ctrl.path_BasinFolder + ctrl.fn_temp_w, ctrl.MapType);       //Channel water temperature
      _FTemp_w = new grid(*_DEM);                                                     //Flux*Water Temp
      _chan_roughness = new grid(ctrl.path_BasinFolder + ctrl.fn_chanrough, ctrl.MapType);//Channel roughness
    }

    // ---------------------------------------------------------------------------------------------------------------        
    /*state variables initialized with user map*/
    // ---------------------------------------------------------------------------------------------------------------        
    _soilmoist1 = new grid(ctrl.path_BasinFolder + ctrl.fn_soilmoist, ctrl.MapType);  //soil moisture volumetric L1
    _soilmoist2 = new grid(ctrl.path_BasinFolder + ctrl.fn_soilmoist2, ctrl.MapType); //soil moisture volumetric L2
    _soilmoist3 = new grid(ctrl.path_BasinFolder + ctrl.fn_soilmoist3, ctrl.MapType); //soil moisture volumetric L3
    _snow = new grid(ctrl.path_BasinFolder + ctrl.fn_swe, ctrl.MapType);              //Current SWE [m]
    _snow_old = new grid(ctrl.path_BasinFolder + ctrl.fn_swe, ctrl.MapType);          //Previous SWE [m]
    _Temp_s_old = new grid(ctrl.path_BasinFolder + ctrl.fn_soiltemp, ctrl.MapType);   //Previous surface temperature [oC]
    _Disch_old =  new grid(ctrl.path_BasinFolder + ctrl.fn_streamflow, ctrl.MapType); //Previous discharge 

    // ---------------------------------------------------------------------------------------------------------------
    /*state variables initialized with the base map*/
    // ---------------------------------------------------------------------------------------------------------------
    // Energy balance and soil properties
    _catcharea = new grid(*_DEM);                                                     //catchment area
    _fieldcapL1 = new grid(*_DEM);                                                    //Field capacity in L1
    _fieldcapL2 = new grid(*_DEM);                                                    //Field capacity in L2
    _fieldcapL3 = new grid(*_DEM);                                                    //Field capacity in L3
    _infilt_cap = new grid(*_DEM);                                                    //Infilt capacity [m.h-1]
    _Rn = new grid(*_DEM);                                                            //Surface Net Radiation [W.m-2]
    _Rn_sum = new grid(*_DEM);                                                        //Total Net Radiation [W.m-2]
    _latheat = new grid(*_DEM);                                                       //Surface Latent Heat [W.m-2]
    _latheat_sum = new grid(*_DEM);                                                   //Total Latent Heat [W.m-2]
    _sensheat = new grid(*_DEM);                                                      //Surface Sensible Heat [W.m-2]
    _sensheat_sum = new grid(*_DEM);                                                      //Surface Sensible Heat [W.m-2]    
    _grndheat = new grid(*_DEM);                                                      //Ground Heat [W.m-2]
    _snwheat = new grid(*_DEM);                                                       //Snow Heat [W.m-2]
    _Temp_s = new grid(*_DEM);                                                        //Surface temperature [oC]
    _CanopyStorage = new grid(*_DEM);                                                 //Canopy Storage [m]
    // Channel storage
    _chan_store = new grid(*_DEM);                                                    //Channel Storage [m]
    _chan_store_old = new grid(*_DEM);                                                //Previous channel storage [m]
    _chan_evap = new grid(*_DEM);                                                     //Channel evaporation
    // Vegetation root properties
    _Zroot95 = new grid(*_DEM);                                                       //Depth at with 95% of roots are found
    _ProotzoneL1 = new grid(*_DEM);                                                   //Root contribution L1
    _ProotzoneL2 = new grid(*_DEM);                                                   //Root contribution L2
    _ProotzoneL3 = new grid(*_DEM);                                                   //Root contribution L3
    // Soil storage
    _IsSaturated = new grid(*_DEM);                                                   //Saturation map
    _incident_water_depth = new grid(*_DEM);  
    _soilmoist_av = new grid(*_DEM);                                                  //Avg VSM of 10 cm (hydrstatic equil)
    _soilmoist_12 = new grid(*_DEM);                                                  //Age VSM L1+L2
    _ponding = new grid(*_DEM);                                                       //Ponding Depth [m]
    _SoilWaterDepth = new grid(*_DEM);                                                //Soil water depth [m]
    _SoilWaterDepthL1 = new grid(*_DEM);                                              //Soil water depth L1 [m]
    _SoilWaterDepthL2 = new grid(*_DEM);                                              //Soil water depth L2 [m]
    _SoilWaterDepthL3 = new grid(*_DEM);                                              //Soil water depth L3 (vadose) [m]
    _WaterTableDepth = new grid(*_DEM);                                               //Reconstructed WTD [m]
    _SoilSatDeficit = new grid(*_DEM);                                                //SM inc. water below and above field capacity
    _GravityWater = new grid(*_DEM);                                                  //Soil water storage above FC
    _GrndWater = new grid(*_DEM);                                                     //GW storage [m]
    // Lateral cell fluxes
    _Layer1upstreamBC = new grid(*_DEM);                                              //L1 flux upstream boundary condition (m.s-1)
    _Layer2upstreamBC = new grid(*_DEM);                                              //L2 flux upstream boundary condition (m.s-1)
    _GWupstreamBC = new grid(*_DEM);                                                  //GW flux upstream boundary condition (m.s-1)
    _Disch_upstreamBC = new grid(*_DEM);                                              //Disch upstream boundary condition (m.s-1)
    _matGWupstreamBC = arma::cube(_NRows, _NCols, 4, arma::fill::zeros);               //Matrix of upstream groundwater
    _FluxLattoChn = new grid(*_DEM);                                                  //Channel inflow
    _FluxLattoSrf = new grid(*_DEM);                                                  //Surface run-on (excluding streamflow)
    _FluxLattoGW = new grid(*_DEM);                                                   //GW lateral outflow
    _FluxChntoLat = new grid(*_DEM);                                                  //Surface run-off (only streamflow)
    _FluxSrftoLat = new grid(*_DEM);                                                  //Surface run-off (excluding streamflow)
    _FluxL1toLat = new grid(*_DEM);                                                   //L1 lateral inflow	
    _FluxL2toLat = new grid(*_DEM);                                                   //L2 lateral inflow	
    _FluxGWtoLat = new grid(*_DEM);                                                   //GW lateral inflow	
    _FluxL1toChn = new grid(*_DEM);                                                   //Intra-cell L1 to channel
    _FluxL2toChn = new grid(*_DEM);                                                   //Intra-cell L2 to channel
    _FluxGWtoChn = new grid(*_DEM);                                                   //Intra-cell GW to channel
    _FluxSrftoChn = new grid(*_DEM);                                                  //Intra-cell ponding to channel
    _Layer1Outlet = new grid(*_DEM);                                                  //Soil outflow @ outlet from L1
    _Layer2Outlet = new grid(*_DEM);                                                  //Soil outflow @ outlet from L2
    _GWOutlet = new grid(*_DEM);                                                      //Soil outflow @ outlet from L3
    // Vertical cell fluxes
    _BedrockLeakageFlux = new grid(*_DEM);                                            //Bedrock leakage flux in [m.s-1]    
    _Evaporation = new grid(*_DEM);                                                   //Actual total evaporation in [m.s-1]
    _EvaporationS_all = new grid(*_DEM);                                              //Actual total soil evaporation in [m.s-1]
    _EvaporationI_all = new grid(*_DEM);                                              //Actual total interc. evaporation [m.s-1]
    _Transpiration_all = new grid(*_DEM);                                             //Actural total transpiration [m.s-1]
    _FluxInfilt = new grid(*_DEM);                                                    //Surface to L1 (summed over timestep)
    _FluxExfilt = new grid(*_DEM);                                                    //L1 to surface (return flow)
    _FluxPercolL2 = new grid(*_DEM);                                                  //L1 to L2 (summed over timestep)
    _FluxPercolL3 = new grid(*_DEM);                                                  //L2 to L3 (summed over timestep)
    _FluxRecharge = new grid(*_DEM);                                                  //Recharge to GW
    _ReturnL1 = new grid(*_DEM);
    _ReturnL2 = new grid(*_DEM);

    // ---------------------------------------------------------------------------------------------------------------
    /* Accumulated fluxes */
    // ---------------------------------------------------------------------------------------------------------------
    _AccInfilt = new grid(*_DEM);                                                     //Infiltration
    _AccExfilt = new grid(*_DEM);                                                     //Exfiltration
    _AccPercolL2 = new grid(*_DEM);                                                   //L1 to L2 (summed over timestep)
    _AccL2toL1 = new grid(*_DEM);                                                     //Capillary + return flow, L2 to L1
    _AccPercolL3 = new grid(*_DEM);                                                   //L2 to L3 (summed over timestep)
    _AccRecharge = new grid(*_DEM);                                                   //Recharge to GW
    _AccLeak = new grid(*_DEM);                                                       //Leakance from L3
    _AccEvaporationS = new grid(*_DEM);                                               //Soil Evaporation
    _AccL3toL2 = new grid(*_DEM);                                                     //Return from L3 to L2
    _AccLattoChn = new grid(*_DEM);                                                   //Channel inflow
    _AccLattoSrf = new grid(*_DEM);                                                   //Surface run-on
    _AccLattoGW = new grid(*_DEM);                                                    //GW lateral outflow
    _AccChntoLat = new grid(*_DEM);                                                   //Surface runoff (channel only)
    _AccSrftoLat = new grid(*_DEM);                                                   //Surface runoff (excluding channel)
    _AccL1toLat = new grid(*_DEM);                                                    //L1 outflow
    _AccL2toLat = new grid(*_DEM);                                                    //L2 outflow
    _AccGWtoLat = new grid(*_DEM);                                                    //GW outflow
    _AccL1toChn = new grid(*_DEM);                                                    //L1 to channel (accumlated)
    _AccL2toChn = new grid(*_DEM);                                                    //L2 to channel (accumlated)
    _AccGWtoChn = new grid(*_DEM);                                                    //GW to channel (accumlated)
    _AccSrftoChn = new grid(*_DEM);                                                   //Ponding to channel (accumulated)

    // ---------------------------------------------------------------------------------------------------------------        
    /* Tracking fluxes specifically for isotopes */
    // ---------------------------------------------------------------------------------------------------------------            
    _psi_MW = NULL;                                                                   //Press separating 2 water pools
    _moist_MW1 = NULL ;                                                               //Moisture of mobile water in L1
    _moist_MW2 = NULL ;                                                               //Moisture of mobile water in L2
    _fracMW1 = NULL ;                                                                 //Frac of mobile water L1
    _fracMW2 = NULL;                                                                  //Frac of mobile water L2
    _fracMW12 = NULL;                                                                 //Frac of mobile water L1+L2
    //Downward
    _FluxCnptoSrf = NULL;                                                             //Canopy water to surface [m]
    _FluxCnptoSnow = NULL;                                                            //Canopy water to snow [m]
    _FluxSnowtoSrf = NULL;                                                            //Snowmelt water to surface [m]
    _FluxSrftoL1 = NULL;                                                              //Surface water to layer 1 [m]
    _FluxL1toL2 = NULL;                                                               //Layer 1 to layer 2 [m]
    _FluxL2toL3 = NULL;                                                               //Layer 2 to layer 3 [m]
    _FluxLeak = NULL;
    //Upward
    _FluxL1toSrf = NULL;                                                              //Layer 1 to surface [m]
    //Mixing
    _IncompAlpha = NULL;                                                              //Frac of water used in mixing (alt. to MW)
		
    if(ctrl.sw_trck){
      _FluxCnptoSrf = new grid(*_DEM);                                                //Canopy/sky to surface
      _FluxCnptoSnow = new grid(*_DEM);                                               //Canopy/sky to snowpac
      _FluxSnowtoSrf = new grid(*_DEM);                                               //Snowpack to surface
      _FluxSrftoL1 = new grid(*_DEM);                                                 //Surface to L1
      _FluxL1toSrf = new grid(*_DEM);                                                 //L1 to surface
      _FluxL1toL2 = new grid(*_DEM);                                                  //Percolation L1 to L2
      _FluxL2toL3 = new grid(*_DEM);                                                  //Percolation L2 to L3
      _FluxLeak = new grid(*_DEM);                                                    //GW to leakance
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

    // ---------------------------------------------------------------------------------------------------------------
    /* Error checking of maps  */
    // ---------------------------------------------------------------------------------------------------------------        
    errno = 0;    // reset errno just to make sure that any constructor error here are not "leftovers"

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
      /* Boundary conditions */
      if(_BCsurface)
	delete _BCsurface;
      if(_BCgroundwater)
	delete _BCgroundwater;
      if(_BCdeepgwtr)
	delete _BCdeepgwtr;      
      /* General basin properties */
      if(fForest)
	delete fForest;
      if(_DEM)
	delete _DEM;
      if(_ldd)
	delete _ldd;
      if(_slope)
	delete _slope;
      if(_fImperv)
	delete _fImperv;
      /* General soil properties */
      if(_Ksat0)
	delete _Ksat0;
      if(_KsatTopSoil)
	delete _KsatTopSoil;
      if(_kKsat)
	delete _kKsat;
      if(_KvKs)
	delete _KvKs;
      if(_random_roughness)
	delete _random_roughness;
      if(_porosity0)
	delete _porosity0;
      if(_kporos)
	delete _kporos;
      if(_psi_ae)
	delete _psi_ae;
      if(_BClambda)
	delete _BClambda;
      if(_theta_rL1)
	delete _theta_rL1;
      if(_theta_rL2)
	delete _theta_rL2;
      if(_theta_rL3)
	delete _theta_rL3;
      if(_soildepth)
	delete _soildepth;
      if(_depth_layer1)
	delete _depth_layer1;
      if(_depth_layer2)
	delete _depth_layer2;
      if(_bedrock_leak)
	delete _bedrock_leak;
      if(_paramWc)
	delete _paramWc;
      if(_paramWp)
	delete _paramWp;
      /* Depth-dependent soil parameters*/
      if(_KsatL1)
	delete _KsatL1;
      if(_KsatL2)
	delete _KsatL2;
      if(_KsatL3)
	delete _KsatL3;
      if(_KvKsL1)
	delete _KvKsL1;
      if(_KvKsL2)
	delete _KvKsL2;
      if(_KvKsL3)
	delete _KvKsL3;
      if(_porosityL1)
	delete _porosityL1;
      if(_porosityL2)
	delete _porosityL2;
      if(_porosityL3)
	delete _porosityL3;
      if(_psi_aeL1)
	delete _psi_aeL1;
      if(_psi_aeL2)
	delete _psi_aeL2;
      if(_psi_aeL3)
	delete _psi_aeL3;
      if(_BClambdaL1)
	delete _BClambdaL1;
      if(_BClambdaL2)
	delete _BClambdaL2;
      if(_BClambdaL3)
	delete _BClambdaL3;
      /* Deep groundwaer */
      if(_DeepGW)
	delete _DeepGW;
      if(_Hydrofrac_DeepGW)
	delete _Hydrofrac_DeepGW;
      if(_chDeepGWparam)
	delete _chDeepGWparam;     
      if(_FluxDeepGWtoChn)
	delete _FluxDeepGWtoChn;
      if(_FluxLattoDeepGW)
	delete _FluxLattoDeepGW;
      if(_FluxDeepGWtoLat)
	delete _FluxDeepGWtoLat;
      if(_AccDeepGWtoChn)
	delete _AccDeepGWtoChn; 
      if(_AccLattoDeepGW)
	delete _AccLattoDeepGW; 
      if(_AccDeepGWtoLat)
	delete _AccDeepGWtoLat; 
      if(_DeepGWupstreamBC)
	delete _DeepGWupstreamBC;
      /* Energy balance */
      if(_meltCoeff)
	delete _meltCoeff;
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
      if(_Temp_L1)
	delete _Temp_L1;
      if(_Temp_L2)
	delete _Temp_L2;
      if(_Temp_L3)
	delete _Temp_L3;      
      /* Channels */
      if(_channelwidth)
	delete _channelwidth;
      if(_channellength)
	delete _channellength;      
      if(_chGWparam)
	delete _chGWparam;      
      if(_Manningn)
	delete _Manningn;
      if(_Temp_w)
	delete _Temp_w;
      if(_FTemp_w)
	delete _FTemp_w;
      if(_chan_roughness)
	delete _chan_roughness;
      /* State-variables*/
      if(_soilmoist1)
	delete _soilmoist1;
      if(_soilmoist2)
	delete _soilmoist2;
      if(_soilmoist3)
	delete _soilmoist3;
      if(_snow)
	delete _snow;
      if(_snow_old)
	delete _snow_old;
      if(_Temp_s_old)
	delete _Temp_s_old;      
      if(_Disch_old)
	delete _Disch_old;
      /*Energy balance and soil prop - state variables*/
      if(_catcharea)
	delete _catcharea;
      if(_fieldcapL1)
	delete _fieldcapL1;
      if(_fieldcapL2)
	delete _fieldcapL2;
      if(_fieldcapL3)
	delete _fieldcapL3;
      if(_infilt_cap)
	delete _infilt_cap;
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
      if(_sensheat_sum)
	delete _sensheat_sum;      
      if(_grndheat)
	delete _grndheat;
      if(_snwheat)
	delete _snwheat;
      if(_Temp_s)
	delete _Temp_s;
      if(_CanopyStorage)
	delete _CanopyStorage;
      /*Channel Storage - state variables*/
      if(_chan_store)
	delete _chan_store;
      if(_chan_store_old)
	delete _chan_store_old;      
      if(_chan_evap)
	delete _chan_evap;
      /*Veg root properties - state variables*/
      if(_Zroot95)
	delete _Zroot95;
      if(_ProotzoneL1)
	delete _ProotzoneL1;
      if(_ProotzoneL2)
	delete _ProotzoneL2;
      if(_ProotzoneL3)
	delete _ProotzoneL3;
      /*Soil storage - state variables*/
      if(_IsSaturated)
	delete _IsSaturated;
      if(_incident_water_depth)
	delete _incident_water_depth;
      if(_soilmoist_av)
	delete _soilmoist_av;
      if(_soilmoist_12)
	delete _soilmoist_12;
      if(_ponding)
	delete _ponding;
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
      if(_GravityWater)
	delete _GravityWater;
      if(_GrndWater)
	delete _GrndWater;
      /* Lateral fluxes - state variables*/
      if(_Layer1upstreamBC)
	delete _Layer1upstreamBC;
      if(_Layer2upstreamBC)
	delete _Layer2upstreamBC;
      if(_GWupstreamBC)
	delete _GWupstreamBC;
      if(_Disch_upstreamBC)
	delete _Disch_upstreamBC;
      //      if(_matGWupstreamBC)
      //	delete _matGWupstreamBC;
      if(_FluxLattoChn)
	delete _FluxLattoChn;
      if(_FluxLattoSrf)
	delete _FluxLattoSrf;
      if(_FluxLattoGW)
	delete _FluxLattoGW;
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
      if(_Layer1Outlet)
	delete _Layer1Outlet;
      if(_Layer2Outlet)
	delete _Layer2Outlet;
      if(_GWOutlet)
	delete _GWOutlet;
      /* Vertical fluxes - state variables*/
      if(_BedrockLeakageFlux)
	delete _BedrockLeakageFlux;
      if(_Evaporation)
	delete _Evaporation;
      if(_EvaporationS_all)
	delete _EvaporationS_all;
      if(_EvaporationI_all)
	delete _EvaporationI_all;
      if(_Transpiration_all)
	delete _Transpiration_all;
      if(_FluxInfilt)
	delete _FluxInfilt;
      if(_FluxExfilt)
	delete _FluxExfilt;
      if(_FluxPercolL2)
	delete _FluxPercolL2;
      if(_FluxPercolL3)
	delete _FluxPercolL3;
      if(_FluxRecharge)
	delete _FluxRecharge;
      if(_ReturnL1)
	delete _ReturnL1;
      if(_ReturnL2)
	delete _ReturnL2;      
      /* Accumulated fluxes - state variables*/
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
      if(_AccEvaporationS)
	delete _AccEvaporationS;
      if(_AccL3toL2)
	delete _AccL3toL2;
      if(_AccLattoChn)
	delete _AccLattoChn;
      if(_AccLattoSrf)
	delete _AccLattoSrf;
      if(_AccLattoGW)
	delete _AccLattoGW;
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
      /* Tracking fluxes - state variables*/
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
      if(_FluxCnptoSrf)
	delete _FluxCnptoSrf;
      if(_FluxCnptoSnow)
	delete _FluxCnptoSnow;
      if(_FluxSnowtoSrf)
	delete _FluxSnowtoSrf;
      if(_FluxSrftoL1)
	delete _FluxSrftoL1;
      if(_FluxL1toL2)
	delete _FluxL1toL2;
      if(_FluxL2toL3)
	delete _FluxL2toL3;
      if(_FluxLeak)
	delete _FluxLeak;
      if(_FluxL1toSrf)
	delete _FluxL1toSrf;
      if(_IncompAlpha)
	delete _IncompAlpha;

      throw;
    }

}

UINT4 rr(UINT4 d){
  if(d<4)
    return 1;
  if(d>6)
    return -1;
  else
    return 0;
}

UINT4 cc(UINT4 d){
  if((d==1) || (d==4) || (d==7))
    return -1;
  if((d==2) || (d==5) || (d==8))
    return 0;
  else
    return 1;
}
