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
 * BasinDestruct.cpp
 *
 *  Created on: Oct 13, 2009
 *      Author: Marco Maneta
 */

#include "Basin.h"

Basin::~Basin(){

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
  //Deep GW
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

}
