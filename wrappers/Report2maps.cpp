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
 * Report2maps.cpp
 *
 *  Created on: Aug 19, 2010
 *      Author: Marco.Maneta
 */

#include "Sativa.h"

int Report2Maps(){

  if(oControl->Rep_Long_Rad_Down)
    WriteMapSeries(oAtmosphere->getIncomingLongWave(), "Ldown", oControl->current_ts_count);
  if(oControl->Rep_Short_Rad_Down)
    WriteMapSeries(oAtmosphere->getIncomingShortWave(), "Sdown", oControl->current_ts_count);
  if(oControl->Rep_Precip)
    WriteMapSeries(oAtmosphere->getPrecipitation(), "Pp", oControl->current_ts_count);
  if(oControl->Rep_Rel_Humidity)
    WriteMapSeries(oAtmosphere->getRelativeHumidty(), "RH", oControl->current_ts_count);
  if(oControl->Rep_Wind_Speed)
    WriteMapSeries(oAtmosphere->getWindSpeed(), "WndSp", oControl->current_ts_count);
  if(oControl->Rep_AvgAir_Temperature)
    WriteMapSeries(oAtmosphere->getTemperature(), "Tp", oControl->current_ts_count);
  if(oControl->Rep_MinAir_Temperature)
    WriteMapSeries(oAtmosphere->getMinTemperature(), "TpMin", oControl->current_ts_count);
  if(oControl->Rep_MaxAir_Temperature)
    WriteMapSeries(oAtmosphere->getMaxTemperature(), "TpMax", oControl->current_ts_count);

  if(oControl->sw_BC){
    if(oControl->Rep_BC_surface)
      WriteMapSeries(oBasin->getBCsurface(), "BCs", oControl->current_ts_count);
    if(oControl->Rep_BC_groundwater)
      WriteMapSeries(oBasin->getBCgroundwater(), "BCgw", oControl->current_ts_count);
  }
  // "Initial" values: stuff that doesn't change over time: only report the first occurence
  // of reportMap_times
  if(oControl->current_t_step == oControl->reportMap_times){

    // Whether each pixel-layer is in the root zone
    // 0 and 1, can be in between when averaging over veg fractions
    if (oControl->Rep_RootZone_in_L1) 
      WriteMap(oBasin->getProotzoneL1(), "RootZone_kL1");
    if (oControl->Rep_RootZone_in_L2)
      WriteMap(oBasin->getProotzoneL2(), "RootZone_kL2");
    if (oControl->Rep_RootZone_in_L3)
      WriteMap(oBasin->getProotzoneL3(), "RootZone_kL3");
    
    // If exponential porosity profile, give field capacity for each.
    // Else, it's the same across the profile
    if(oControl->sw_expPoros){
      if (oControl->Rep_Field_Capacity_L1) 
	WriteMap(oBasin->getFieldCapacityL1(), "FCap1");
      if (oControl->Rep_Field_Capacity_L2)
	WriteMap(oBasin->getFieldCapacityL2(), "FCap2");
      if (oControl->Rep_Field_Capacity_L3)
	WriteMap(oBasin->getFieldCapacityL3(), "FCap3");
    } else
      if (oControl->Rep_Field_Capacity_L1 or oControl->Rep_Field_Capacity_L2 or
	  oControl->Rep_Field_Capacity_L3)
	WriteMap(oBasin->getFieldCapacityL1(), "FCap");
  }
  
  if (oControl->Rep_Canopy_Water_Stor_sum)
    WriteMapSeries(oBasin->getCanopyStorage(), "Cs", oControl->current_ts_count);
  if (oControl->Rep_SWE)
    WriteMapSeries(oBasin->getSnowWaterEquiv(), "SWE", oControl->current_ts_count);
  if (oControl->Rep_Infilt_Cap)
    WriteMapSeries(oBasin->getInfiltCap(), "IfCap", oControl->current_ts_count);
  if (oControl->Rep_Streamflow)
    WriteMapSeries(oBasin->getStreamflow(), "Q", oControl->current_ts_count);
  if (oControl->Rep_Saturation_Area)
    WriteMapSeries(oBasin->getSatArea(), "SatArea_", oControl->current_ts_count);
  if (oControl->Rep_Ponding)
    WriteMapSeries(oBasin->getPondingWater(), "Ponding_", oControl->current_ts_count);
  if (oControl->Rep_Soil_Water_Content_Average)
    WriteMapSeries(oBasin->getSoilMoist_av(), "SWCav", oControl->current_ts_count);
  if (oControl->Rep_Soil_Water_Content_12)
    WriteMapSeries(oBasin->getSoilMoist_12(), "SWCup", oControl->current_ts_count);
  if (oControl->Rep_Soil_Water_Content_L1)
    WriteMapSeries(oBasin->getSoilMoist1(), "SWC1_", oControl->current_ts_count);
  if (oControl->Rep_Soil_Water_Content_L2)
    WriteMapSeries(oBasin->getSoilMoist2(), "SWC2_", oControl->current_ts_count);
  if (oControl->Rep_Soil_Water_Content_L3)
    WriteMapSeries(oBasin->getSoilMoist3(), "SWC3_", oControl->current_ts_count);
  if (oControl->Rep_WaterTableDepth)
    WriteMapSeries(oBasin->getWaterTableDepth(), "WTD_", oControl->current_ts_count);
  if (oControl->Rep_Soil_Sat_Deficit)
    WriteMapSeries(oBasin->getSaturationDeficit(), "SatDef", oControl->current_ts_count);
  if (oControl->Rep_GWater)
    WriteMapSeries(oBasin->getGrndWater(), "GW", oControl->current_ts_count);
  if(oControl->sw_deepGW){
    if (oControl->Rep_DeepGWater)
      WriteMapSeries(oBasin->getDeepGW(), "deepGW", oControl->current_ts_count);
  }
  if (oControl->Rep_Soil_Net_Rad)
    WriteMapSeries(oBasin->getNetRad(), "NRs", oControl->current_ts_count);
  if (oControl->Rep_Net_Rad_sum)
    WriteMapSeries(oBasin->getNetRad_sum(), "NRtot", oControl->current_ts_count);
  if (oControl->Rep_Soil_LE)
    WriteMapSeries(oBasin->getLatheat(), "LEs", oControl->current_ts_count);
  if (oControl->Rep_LE_sum)
    WriteMapSeries(oBasin->getLatheat_sum(), "LEtot", oControl->current_ts_count);
  if (oControl->Rep_Sens_Heat)
    WriteMapSeries(oBasin->getSensHeat(), "SensH", oControl->current_ts_count);
  if (oControl->Rep_Grnd_Heat)
    WriteMapSeries(oBasin->getGrndHeat(), "GrndH", oControl->current_ts_count);
  if (oControl->Rep_Snow_Heat)
    WriteMapSeries(oBasin->getSnwHeat(), "SnowH", oControl->current_ts_count);
  if (oControl->Rep_Soil_Temperature)
    WriteMapSeries(oBasin->getSoilTemp(), "Ts", oControl->current_ts_count);
  if (oControl->Rep_Skin_Temperature)
    WriteMapSeries(oBasin->getSkinTemp(), "Tskin", oControl->current_ts_count);

  if(oControl->toggle_chan_evap > 0){
    if (oControl->Rep_ChanEvap)
      WriteMapSeries(oBasin->getChanEvap(), "EvapC", oControl->current_ts_count);
    if (oControl->toggle_chan_evap == 1){
      if (oControl->Rep_Water_Temperature)
	WriteMapSeries(oBasin->getTemp_w(), "Twater", oControl->current_ts_count);
    }
  }

  if (oControl->Rep_Total_ET)
    WriteMapSeries(oBasin->getEvaporation(), "Evap", oControl->current_ts_count);
  if (oControl->Rep_Transpiration_sum)
    WriteMapSeries(oBasin->getTranspiration_all(), "EvapT", oControl->current_ts_count);
  if (oControl->Rep_Einterception_sum)
    WriteMapSeries(oBasin->getEvaporationI_all(), "EvapI", oControl->current_ts_count);
  if (oControl->Rep_Esoil_sum)
    WriteMapSeries(oBasin->getEvaporationS_all(), "EvapS", oControl->current_ts_count);

  // Internal vertical fluxes
  if (oControl->Rep_GWtoChn)
    WriteMapSeries(oBasin->getFluxGWtoChn(), "GWChn", oControl->current_ts_count);
  if(oControl->sw_deepGW){
    if (oControl->Rep_DeepGWtoChn)
      WriteMapSeries(oBasin->getFluxDeepGWtoChn(), "DeepGWChn", oControl->current_ts_count);
  }
  if (oControl->Rep_SrftoChn)
    WriteMapSeries(oBasin->getFluxSrftoChn(), "SrfChn", oControl->current_ts_count);
  if (oControl->Rep_Infilt)
    WriteMapSeries(oBasin->getFluxInfilt(), "Inf", oControl->current_ts_count);
  if (oControl->Rep_Exfilt)
    WriteMapSeries(oBasin->getFluxExfilt(), "RSrf", oControl->current_ts_count);
  if (oControl->Rep_PercolL2)
    WriteMapSeries(oBasin->getFluxPercolL2(), "PrcL2", oControl->current_ts_count);
  if (oControl->Rep_ReturnL1)
    WriteMapSeries(oBasin->getFluxL2toL1(), "RL1", oControl->current_ts_count);
  if (oControl->Rep_PercolL3)
    WriteMapSeries(oBasin->getFluxPercolL3(), "PrcL3", oControl->current_ts_count);
  if (oControl->Rep_Recharge)
    WriteMapSeries(oBasin->getFluxRecharge(), "Rchg", oControl->current_ts_count);
  if (oControl->Rep_ReturnL2)
    WriteMapSeries(oBasin->getFluxL3toL2(), "RL2", oControl->current_ts_count);
  if (oControl->Rep_Leak)
    WriteMapSeries(oBasin->getBedrockLeakage(),"Leak", oControl->current_ts_count);
  // Internal lateral fluxes
  if (oControl->Rep_LattoSrf)
    WriteMapSeries(oBasin->getFluxLattoSrf(), "LSrfi", oControl->current_ts_count);
  if (oControl->Rep_LattoChn)
    WriteMapSeries(oBasin->getFluxLattoChn(), "LChni", oControl->current_ts_count);
  if (oControl->Rep_LattoGW)
    WriteMapSeries(oBasin->getFluxLattoGW(), "LGWi", oControl->current_ts_count);
  if(oControl->sw_deepGW){
    if (oControl->Rep_LattoDeepGW)
      WriteMapSeries(oBasin->getFluxLattoDeepGW(), "LDeepGWi", oControl->current_ts_count);
  }
  if (oControl->Rep_ChntoLat)
    WriteMapSeries(oBasin->getFluxChntoLat(), "LChno", oControl->current_ts_count);
  if (oControl->Rep_SrftoLat)
    WriteMapSeries(oBasin->getFluxSrftoLat(), "LSrfo", oControl->current_ts_count);
  if (oControl->Rep_GWtoLat)
    WriteMapSeries(oBasin->getFluxGWtoLat(), "LGWo", oControl->current_ts_count);
  if(oControl->sw_deepGW){
    if (oControl->Rep_DeepGWtoLat)
      WriteMapSeries(oBasin->getFluxDeepGWtoLat(), "LDeepGWo", oControl->current_ts_count);
  }
  // Cumulative values
  if (oControl->Rep_GWtoChnacc)
    WriteMapSeries(oBasin->getAccGWtoChn(), "GWChnA", oControl->current_ts_count);
  if(oControl->sw_deepGW){
    if (oControl->Rep_DeepGWtoChnacc)
      WriteMapSeries(oBasin->getAccDeepGWtoChn(), "DeepGWChnA", oControl->current_ts_count);
  }
  if (oControl->Rep_SrftoChnacc)
    WriteMapSeries(oBasin->getAccSrftoChn(), "SrfChnA", oControl->current_ts_count);
  if (oControl->Rep_Infiltacc)
    WriteMapSeries(oBasin->getAccInfilt(), "InfA", oControl->current_ts_count);
  if (oControl->Rep_Exfiltacc)
    WriteMapSeries(oBasin->getAccExfilt(), "RSrfA", oControl->current_ts_count);
  if (oControl->Rep_PercolL2acc)
    WriteMapSeries(oBasin->getAccPercolL2(), "PrcL2A", oControl->current_ts_count);
  if (oControl->Rep_ReturnL1acc)
    WriteMapSeries(oBasin->getAccL2toL1(), "RL1A", oControl->current_ts_count);
  if (oControl->Rep_PercolL3acc)
    WriteMapSeries(oBasin->getAccPercolL3(), "PrcL3A", oControl->current_ts_count);
  if (oControl->Rep_Rechargeacc)
    WriteMapSeries(oBasin->getAccRecharge(), "RchgA", oControl->current_ts_count);
  if (oControl->Rep_ReturnL2acc)
    WriteMapSeries(oBasin->getAccL3toL2(), "RL2A", oControl->current_ts_count);
  if (oControl->Rep_EvaporationSacc)
    WriteMapSeries(oBasin->getAccEvaporationS(), "EsA", oControl->current_ts_count);
  if (oControl->Rep_LattoSrfacc)
    WriteMapSeries(oBasin->getAccLattoSrf(), "LSrfiA", oControl->current_ts_count);
  if (oControl->Rep_LattoChnacc)
    WriteMapSeries(oBasin->getAccLattoChn(), "LChniA", oControl->current_ts_count);
  if (oControl->Rep_LattoGWacc)
    WriteMapSeries(oBasin->getAccLattoGW(), "LGWiA", oControl->current_ts_count);
  if(oControl->sw_deepGW){
    if (oControl->Rep_LattoDeepGWacc)
      WriteMapSeries(oBasin->getAccLattoDeepGW(), "LDeepGWiA", oControl->current_ts_count);
  }
  if (oControl->Rep_SrftoLatacc)
    WriteMapSeries(oBasin->getAccSrftoLat(), "LSrfoA", oControl->current_ts_count);
  if (oControl->Rep_GWtoLatacc)
    WriteMapSeries(oBasin->getAccGWtoLat(), "LGWoA", oControl->current_ts_count);
  if(oControl->sw_deepGW){
    if (oControl->Rep_DeepGWtoLatacc)
      WriteMapSeries(oBasin->getAccDeepGWtoLat(), "LDeepGWoA", oControl->current_ts_count);
  }
  if (oControl->Rep_ChntoLatacc)
    WriteMapSeries(oBasin->getAccChntoLat(), "LChnoA", oControl->current_ts_count);

  if(oControl->sw_trck and oControl->sw_TPD){
    if (oControl->Rep_Moist_MW1)
      WriteMapSeries(oBasin->getMoistureMW1(), "SWCMW1", oControl->current_ts_count);
    if (oControl->Rep_Moist_MW2)
      WriteMapSeries(oBasin->getMoistureMW2(), "SWCMW2", oControl->current_ts_count);
    if (oControl->Rep_Frac_MW1)
      WriteMapSeries(oBasin->getFracMW1(), "fMW1", oControl->current_ts_count);
    if (oControl->Rep_Frac_MW2)
      WriteMapSeries(oBasin->getFracMW2(), "fMW2", oControl->current_ts_count);
    if (oControl->Rep_Frac_MW12)
      WriteMapSeries(oBasin->getFracMW12(), "fMW12", oControl->current_ts_count);
  }

  // Tracking maps
  if(oControl->sw_trck && oControl->sw_2H){
    if (oControl->Rep_d2Hcanopy_sum)
      WriteMapSeries(oTracking->getd2Hcanopy_sum(), "dHcnp", oControl->current_ts_count);
    if (oControl->Rep_d2Hprecip)
      WriteMapSeries(oAtmosphere->getd2Hprecip(), "dHpcp", oControl->current_ts_count);
    if (oControl->Rep_d2Hsnowpack)
      WriteMapSeries(oTracking->getd2Hsnowpack(), "dHsnw", oControl->current_ts_count);
    if (oControl->Rep_d2Hsurface)
      WriteMapSeries(oTracking->getd2Hsurface(), "dHsrf", oControl->current_ts_count);
    if (oControl->Rep_d2Hchan)
      WriteMapSeries(oTracking->getd2Hchan(), "dHchn", oControl->current_ts_count);
    if (oControl->Rep_d2Hsoil1)
      WriteMapSeries(oTracking->getd2Hsoil1(), "dHsL1", oControl->current_ts_count);
    if (oControl->Rep_d2Hsoil2)
      WriteMapSeries(oTracking->getd2Hsoil2(), "dHsL2", oControl->current_ts_count);
    if (oControl->Rep_d2HsoilUp)
      WriteMapSeries(oTracking->getd2Hsoil_12(), "dHsUp", oControl->current_ts_count);
    if (oControl->Rep_d2Hsoil3)
      WriteMapSeries(oTracking->getd2Hsoil3(), "dHsL3", oControl->current_ts_count);
    if (oControl->Rep_d2HsoilAv)
      WriteMapSeries(oTracking->getd2Hsoil_Av(), "dHsAv", oControl->current_ts_count);
    if (oControl->Rep_d2Hgroundwater)
      WriteMapSeries(oTracking->getd2Hgroundwater(), "dHgw", oControl->current_ts_count);
    if(oControl->sw_deepGW){
      if (oControl->Rep_d2HdeepGW)
	WriteMapSeries(oTracking->getd2HdeepGW(), "dHdeepgw", oControl->current_ts_count);
      if (oControl->Rep_d2HdeepGWQ)
	WriteMapSeries(oTracking->getd2HdeepGWQ(), "dHdeepgwQ", oControl->current_ts_count);
    }
    if (oControl->Rep_d2Hleakage)
      WriteMapSeries(oTracking->getd2Hleakage(), "dHlk", oControl->current_ts_count);
    if (oControl->Rep_d2HevapS_sum)
      WriteMapSeries(oTracking->getd2HevapS_sum(), "dHeS", oControl->current_ts_count);
    if (oControl->Rep_d2HevapI_sum)
      WriteMapSeries(oTracking->getd2HevapI_sum(), "dHeI", oControl->current_ts_count);
    if (oControl->Rep_d2HevapT_sum)
      WriteMapSeries(oTracking->getd2HevapT_sum(), "dHeT", oControl->current_ts_count);
    if(oControl->sw_TPD){
      if (oControl->Rep_d2H_MW1)
	WriteMapSeries(oTracking->getd2H_MW1(), "dHMW1", oControl->current_ts_count);
      if (oControl->Rep_d2H_MW2)
	WriteMapSeries(oTracking->getd2H_MW2(), "dHMW2", oControl->current_ts_count);
      if (oControl->Rep_d2H_TB1)
	WriteMapSeries(oTracking->getd2H_TB1(), "dHTB1", oControl->current_ts_count);
      if (oControl->Rep_d2H_TB2)
	WriteMapSeries(oTracking->getd2H_TB2(), "dHTB2", oControl->current_ts_count);
    }
  }

  if(oControl->sw_trck && oControl->sw_18O){
    if (oControl->Rep_d18Ocanopy_sum)
     WriteMapSeries(oTracking->getd18Ocanopy_sum(), "dOcnp", oControl->current_ts_count);
    if (oControl->Rep_d18Oprecip)
      WriteMapSeries(oAtmosphere->getd18Oprecip(), "dOpcp", oControl->current_ts_count);
    if (oControl->Rep_d18Osnowpack)
      WriteMapSeries(oTracking->getd18Osnowpack(), "dOsnw", oControl->current_ts_count);
    if (oControl->Rep_d18Osurface)
      WriteMapSeries(oTracking->getd18Osurface(), "dOsrf", oControl->current_ts_count);
    if (oControl->Rep_d18Ochan)
      WriteMapSeries(oTracking->getd18Ochan(), "dOchn", oControl->current_ts_count);
    if (oControl->Rep_d18Osoil1)
      WriteMapSeries(oTracking->getd18Osoil1(), "dOsL1", oControl->current_ts_count);
    if (oControl->Rep_d18Osoil2)
      WriteMapSeries(oTracking->getd18Osoil2(), "dOsL2", oControl->current_ts_count);
    if (oControl->Rep_d18OsoilUp)
      WriteMapSeries(oTracking->getd18Osoil_12(), "dOsUp", oControl->current_ts_count);
    if (oControl->Rep_d18Osoil3)
      WriteMapSeries(oTracking->getd18Osoil3(), "dOsL3", oControl->current_ts_count);
    if (oControl->Rep_d18OsoilAv)
      WriteMapSeries(oTracking->getd18Osoil_Av(), "dOsAv", oControl->current_ts_count);
    if (oControl->Rep_d18Ogroundwater)
      WriteMapSeries(oTracking->getd18Ogroundwater(), "dOgw", oControl->current_ts_count);
    if(oControl->sw_deepGW){
      if (oControl->Rep_d18OdeepGW)
	WriteMapSeries(oTracking->getd18OdeepGW(), "dOdeepgw", oControl->current_ts_count);
      if (oControl->Rep_d18OdeepGWQ)
	WriteMapSeries(oTracking->getd18OdeepGWQ(), "dOdeepgwQ", oControl->current_ts_count);
    }    
    if (oControl->Rep_d18Oleakage)
      WriteMapSeries(oTracking->getd18Oleakage(), "dOlk", oControl->current_ts_count);
    if (oControl->Rep_d18OevapS_sum)
      WriteMapSeries(oTracking->getd18OevapS_sum(), "dOeS", oControl->current_ts_count);
    if (oControl->Rep_d18OevapI_sum)
      WriteMapSeries(oTracking->getd18OevapI_sum(), "dOeI", oControl->current_ts_count);
    if (oControl->Rep_d18OevapT_sum)
      WriteMapSeries(oTracking->getd18OevapT_sum(), "dOeT", oControl->current_ts_count);
    if(oControl->sw_TPD){
      if (oControl->Rep_d18O_MW1)
	WriteMapSeries(oTracking->getd18O_MW1(), "dOMW1", oControl->current_ts_count);
      if (oControl->Rep_d18O_MW2)
	WriteMapSeries(oTracking->getd18O_MW2(), "dOMW2", oControl->current_ts_count);
      if (oControl->Rep_d18O_TB1)
	WriteMapSeries(oTracking->getd18O_TB1(), "dOTB1", oControl->current_ts_count);
      if (oControl->Rep_d18O_TB2)
	WriteMapSeries(oTracking->getd18O_TB2(), "dOTB2", oControl->current_ts_count);
    }
  }
  if(oControl->sw_trck && oControl->sw_Age){
    if (oControl->Rep_Agecanopy_sum)
     WriteMapSeries(oTracking->getAgecanopy_sum(), "Agecnp", oControl->current_ts_count);
    if (oControl->Rep_Agesnowpack)
      WriteMapSeries(oTracking->getAgesnowpack(), "Agesnw", oControl->current_ts_count);
    if (oControl->Rep_Agesurface)
      WriteMapSeries(oTracking->getAgesurface(), "Agesrf", oControl->current_ts_count);
    if (oControl->Rep_Agechan)
      WriteMapSeries(oTracking->getAgechan(), "Agechn", oControl->current_ts_count);
    if (oControl->Rep_Agesoil1)
      WriteMapSeries(oTracking->getAgesoil1(), "AgesL1", oControl->current_ts_count);
    if (oControl->Rep_Agesoil2)
      WriteMapSeries(oTracking->getAgesoil2(), "AgesL2", oControl->current_ts_count);
    if (oControl->Rep_AgesoilUp)
      WriteMapSeries(oTracking->getAgesoil_12(), "AgesUp", oControl->current_ts_count);
    if (oControl->Rep_Agesoil3)
      WriteMapSeries(oTracking->getAgesoil3(), "AgesL3", oControl->current_ts_count);
    if (oControl->Rep_AgesoilAv)
      WriteMapSeries(oTracking->getAgesoil_Av(), "AgesAv", oControl->current_ts_count);
    if (oControl->Rep_Agegroundwater)
      WriteMapSeries(oTracking->getAgegroundwater(), "Agegw", oControl->current_ts_count);
    if (oControl->Rep_Ageleakage)
      WriteMapSeries(oTracking->getAgeleakage(), "Agelk", oControl->current_ts_count);
    if (oControl->Rep_AgeevapS_sum)
      WriteMapSeries(oTracking->getAgeevapS_sum(), "AgeeS", oControl->current_ts_count);
    if (oControl->Rep_AgeevapI_sum)
      WriteMapSeries(oTracking->getAgeevapI_sum(), "AgeeI", oControl->current_ts_count);
    if (oControl->Rep_AgeevapT_sum)
      WriteMapSeries(oTracking->getAgeevapT_sum(), "AgeeT", oControl->current_ts_count);
    if (oControl->Rep_AgeGWtoChn)
      WriteMapSeries(oTracking->getAgeGWtoChn(), "AgeGWQ", oControl->current_ts_count);
    if(oControl->sw_deepGW){
      if (oControl->Rep_AgedeepGW)
	WriteMapSeries(oTracking->getAgedeepGW(), "Agedeepgw", oControl->current_ts_count);
      if (oControl->Rep_AgedeepGWQ)
	WriteMapSeries(oTracking->getAgedeepGWQ(), "AgedeepgwQ", oControl->current_ts_count);
    }    
    if (oControl->Rep_AgeSrftoChn)
      WriteMapSeries(oTracking->getAgeSrftoChn(), "AgeSrfQ", oControl->current_ts_count);
    if (oControl->Rep_AgeRecharge)
      WriteMapSeries(oTracking->getAgeRecharge(), "AgeR", oControl->current_ts_count);

    if(oControl->sw_TPD){
      if (oControl->Rep_Age_MW1)
	WriteMapSeries(oTracking->getAge_MW1(), "AgeMW1", oControl->current_ts_count);
      if (oControl->Rep_Age_MW2)
	WriteMapSeries(oTracking->getAge_MW2(), "AgeMW2", oControl->current_ts_count);
      if (oControl->Rep_Age_MWUp)
	WriteMapSeries(oTracking->getAge_MW12(), "AgeMW12", oControl->current_ts_count);
      if (oControl->Rep_Age_TB1)
	WriteMapSeries(oTracking->getAge_TB1(), "AgeTB1", oControl->current_ts_count);
      if (oControl->Rep_Age_TB2)
	WriteMapSeries(oTracking->getAge_TB2(), "AgeTB2", oControl->current_ts_count);
      if (oControl->Rep_Age_TBUp)
	WriteMapSeries(oTracking->getAge_TB12(), "AgeTB12", oControl->current_ts_count);
    }
  }

  for (int i = 0; i < oControl->NumSpecs; i++) {
    stringstream name;

    if (oControl->Rep_Veget_frac) {
      name << "p" << i << "_";
      //WriteMapSeries(oBasin->getVegetFrac(i), name.str() , oControl->current_ts_count);
      WriteMap(oBasin->getVegetFrac(i), name.str());
      name.str("");
    }

    if (oControl->Rep_Stem_Density) {
      name << "ntr" << i << "_";
      //WriteMapSeries(oBasin->getStemDensity(i), name.str() , oControl->current_ts_count);
      WriteMap(oBasin->getStemDensity(i), name.str());
      name.str("");
    }

    if (oControl->Rep_RootFrac1Species) {
      name << "L1Ro" << i << "_";
      //WriteMapSeries(oBasin->getRootFrac1(i), name.str() , oControl->current_ts_count);
      WriteMap(oBasin->getRootFrac1(i), name.str());
      name.str("");
    }
    if (oControl->Rep_RootFrac2Species) {
      name << "L2Ro" << i << "_";
      //WriteMapSeries(oBasin->getRootFrac2(i), name.str() , oControl->current_ts_count);
      WriteMap(oBasin->getRootFrac2(i), name.str());
      name.str("");
    }

    if (oControl->Rep_Leaf_Area_Index) {
      name << "lai" << i << "_";
      WriteMapSeries(oBasin->getLAI(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->Rep_Stand_Age) {
      name << "age" << i << "_";
      WriteMapSeries(oBasin->getStandAge(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->Rep_Canopy_Conductance) {
      name << "gc" << i << "_";
      WriteMapSeries(oBasin->getCanopyCond(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->Rep_GPP) {
      name << "gpp" << i << "_";
      WriteMapSeries(oBasin->getGPP(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->Rep_NPP) {
      name << "npp" << i << "_";
      WriteMapSeries(oBasin->getNPP(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->Rep_Basal_Area) {
      name << "bas" << i << "_";
      WriteMapSeries(oBasin->getBasalArea(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->Rep_Tree_Height) {
      name << "hgt" << i << "_";
      //WriteMapSeries(oBasin->getTreeHeight(i), name.str() , oControl->current_ts_count);
      WriteMap(oBasin->getTreeHeight(i), name.str());
      name.str("");
    }

    if (oControl->Rep_Root_Mass) {
      name << "root" << i << "_";
      //WriteMapSeries(oBasin->getRootMass(i), name.str() , oControl->current_ts_count);
      WriteMap(oBasin->getRootMass(i), name.str());
      name.str("");
    }

    if (oControl->Rep_Canopy_Temp) {
      name << "Tc" << i << "_";
      WriteMapSeries(oBasin->getCanopyTemp(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->Rep_Canopy_NetR) {
      name << "NRc" << i << "_";
      WriteMapSeries(oBasin->getCanopyNetRad(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->Rep_Canopy_LE_E) {
      name << "LEEi" << i << "_";
      WriteMapSeries(oBasin->getCanopyLatHeatE(i), name.str() , oControl->current_ts_count);
      name.str("");
    }
    if (oControl->Rep_Canopy_LE_T) {
      name << "LETr" << i << "_";
      WriteMapSeries(oBasin->getCanopyLatHeatT(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->Rep_Canopy_Sens_Heat) {
      name << "Hc" << i << "_";
      WriteMapSeries(oBasin->getCanopySensHeat(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->Rep_Canopy_Water_Stor) {
      name << "Cs" << i << "_";
      WriteMapSeries(oBasin->getCanopyWaterStor(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->Rep_ETspecies) {
      name << "ETc" << i << "_";
      WriteMapSeries(oBasin->getETspecies(i), name.str() , oControl->current_ts_count);
      name.str("");
    }
    if (oControl->Rep_Transpiration) {
      name << "Tr" << i << "_";
      WriteMapSeries(oBasin->getTranspiration(i), name.str() , oControl->current_ts_count);
      name.str("");
    }
    if (oControl->Rep_Einterception) {
      name << "Ei" << i << "_";
      WriteMapSeries(oBasin->getEinterception(i), name.str() , oControl->current_ts_count);
      name.str("");
    }
    if (oControl->Rep_Esoil) {
      name << "Es" << i << "_";
      WriteMapSeries(oBasin->getEsoil(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    // Tracking
    if(oControl->sw_trck){
      if (oControl->sw_2H && oControl->Rep_d2Hcanopy) {
	name << "dHCs" << i << "_";
	WriteMapSeries(oBasin->getd2Hcanopy(i), name.str() , oControl->current_ts_count);
	name.str("");
      }
      
      if (oControl->sw_18O && oControl->Rep_d18Ocanopy) {
	name << "dOCs" << i << "_";
	WriteMapSeries(oBasin->getd18Ocanopy(i), name.str() , oControl->current_ts_count);
	name.str("");
      }
      
      if (oControl->sw_Age && oControl->Rep_Agecanopy) {
	name << "AgeCs" << i << "_";
	WriteMapSeries(oBasin->getAgecanopy(i), name.str() , oControl->current_ts_count);
	name.str("");
      }
      if (oControl->sw_2H && oControl->Rep_d2HevapI) {
	name << "dHeI" << i << "_";
	WriteMapSeries(oBasin->getd2HevapI(i), name.str() , oControl->current_ts_count);
	name.str("");
      }

      if (oControl->sw_18O && oControl->Rep_d18OevapI) {
	name << "dOeI" << i << "_";
	WriteMapSeries(oBasin->getd18OevapI(i), name.str() , oControl->current_ts_count);
	name.str("");
      }

      if (oControl->sw_Age && oControl->Rep_AgeevapI) {
	name << "AgeeI" << i << "_";
	WriteMapSeries(oBasin->getAgeevapI(i), name.str() , oControl->current_ts_count);
	name.str("");
      }
      if (oControl->sw_2H && oControl->Rep_d2HevapT) {
	name << "dHeT" << i << "_";
	WriteMapSeries(oBasin->getd2HevapT(i), name.str() , oControl->current_ts_count);
	name.str("");
      }

      if (oControl->sw_18O && oControl->Rep_d18OevapT) {
	name << "dOeT" << i << "_";
	WriteMapSeries(oBasin->getd18OevapT(i), name.str() , oControl->current_ts_count);
	name.str("");
      }

      if (oControl->sw_Age && oControl->Rep_AgeevapT) {
	name << "AgeeT" << i << "_";
	WriteMapSeries(oBasin->getAgeevapT(i), name.str() , oControl->current_ts_count);
	name.str("");
      }
      if (oControl->sw_2H && oControl->Rep_d2HevapS) {
	name << "dHeS" << i << "_";
	WriteMapSeries(oBasin->getd2HevapS(i), name.str() , oControl->current_ts_count);
	name.str("");
      }

      if (oControl->sw_18O && oControl->Rep_d18OevapS) {
	name << "dOeS" << i << "_";
	WriteMapSeries(oBasin->getd18OevapS(i), name.str() , oControl->current_ts_count);
	name.str("");
      }

      if (oControl->sw_Age && oControl->Rep_AgeevapS) {
	name << "AgeeS" << i << "_";
	WriteMapSeries(oBasin->getAgeevapS(i), name.str() , oControl->current_ts_count);
	name.str("");
      }
    }

  }

  return EXIT_SUCCESS;
}
////////////////////////////////////////Time series output section/////////////////////////////////////////

int Report2Ts(){

  if(oControl->RepTs_OutletDischarge){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "OutletDisch.tab");
    oReport->ReportVectCells(oBasin->getDailyOvlndOutput(),
			     oControl->path_ResultsFolder + "OutletDisch.tab",
			     oControl->current_ts_count);
  }

  if(oControl->RepTs_OutletGW){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "OutletGW.tab");
    oReport->ReportVectCells(oBasin->getDailyGwtrOutput(),
			     oControl->path_ResultsFolder + "OutletGW.tab",
			     oControl->current_ts_count);
  }

  if (oControl->RepTs_Long_Rad_Down){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "Ldown.tab");
    oReport->ReportTimeSeries(oAtmosphere->getIncomingLongWave(),
			      oControl->path_ResultsFolder + "Ldown.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Short_Rad_Down){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "Sdown.tab");
    oReport->ReportTimeSeries(oAtmosphere->getIncomingShortWave(),
			      oControl->path_ResultsFolder + "Sdown.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Precip){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "Precip.tab");
    oReport->ReportTimeSeries(oAtmosphere->getPrecipitation(),
			      oControl->path_ResultsFolder + "Precip.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Rel_Humidity){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "RelHumid.tab");
    oReport->ReportTimeSeries(oAtmosphere->getRelativeHumidty(),
			      oControl->path_ResultsFolder + "RelHumid.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Wind_Speed){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "WindSpeed.tab");
    oReport->ReportTimeSeries(oAtmosphere->getWindSpeed(),
			      oControl->path_ResultsFolder + "WindSpeed.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_AvgAir_Temperature){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "AvgTemp.tab");
    oReport->ReportTimeSeries(oAtmosphere->getTemperature(),
			      oControl->path_ResultsFolder + "AvgTemp.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_MinAir_Temperature){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "MinTemp.tab");
    oReport->ReportTimeSeries(oAtmosphere->getMinTemperature(),
			      oControl->path_ResultsFolder + "MinTemp.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_MaxAir_Temperature){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "MaxTemp.tab");
    oReport->ReportTimeSeries(oAtmosphere->getMaxTemperature(),
			      oControl->path_ResultsFolder + "MaxTemp.tab",
			      oControl->current_ts_count);
  }

  // If exponential porosity profile, give field capacity for each.
  // Else, it's the same across the profile
  if(oControl->sw_expPoros){
    if (oControl->RepTs_Field_Capacity_L1){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "FieldCapL1.tab");
      oReport->ReportTimeSeries(oBasin->getFieldCapacityL1(),
				oControl->path_ResultsFolder + "FieldCapL1.tab",
				oControl->current_ts_count);
    }
    if (oControl->RepTs_Field_Capacity_L2){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "FieldCapL2.tab");
      oReport->ReportTimeSeries(oBasin->getFieldCapacityL2(),
				oControl->path_ResultsFolder + "FieldCapL2.tab",
				oControl->current_ts_count);
    }
    if (oControl->RepTs_Field_Capacity_L3){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "FieldCapL3.tab");
      oReport->ReportTimeSeries(oBasin->getFieldCapacityL3(),
				oControl->path_ResultsFolder + "FieldCapL3.tab",
				oControl->current_ts_count);
    }
  } else
    if (oControl->RepTs_Field_Capacity_L1 or oControl->RepTs_Field_Capacity_L2 or
	oControl->RepTs_Field_Capacity_L3){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "FieldCap.tab");
      oReport->ReportTimeSeries(oBasin->getFieldCapacityL3(),
				oControl->path_ResultsFolder + "FieldCap.tab",
				oControl->current_ts_count);
    }

  
  if (oControl->RepTs_Canopy_Water_Stor_sum){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "CanopyWaterStor.tab");
    oReport->ReportTimeSeries(oBasin->getCanopyStorage(),
			      oControl->path_ResultsFolder + "CanopyWaterStor.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_SWE){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "SWE.tab");
    oReport->ReportTimeSeries(oBasin->getSnowWaterEquiv(),
			      oControl->path_ResultsFolder + "SWE.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Infilt_Cap){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "InfiltCap.tab");
    oReport->ReportTimeSeries(oBasin->getInfiltCap(),
			      oControl->path_ResultsFolder + "InfiltCap.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Streamflow){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "Streamflow.tab");
    oReport->ReportTimeSeries(oBasin->getStreamflow(),
			      oControl->path_ResultsFolder + "Streamflow.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Ponding){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "Ponding.tab");
    oReport->ReportTimeSeries(oBasin->getChanStore(),
			      oControl->path_ResultsFolder + "Ponding.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Soil_Water_Content_Average){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "SoilMoistureAv.tab");
    oReport->ReportTimeSeries(oBasin->getSoilMoist_av(),
			      oControl->path_ResultsFolder + "SoilMoistureAv.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Soil_Water_Content_12){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "SoilMoistureUp.tab");
    oReport->ReportTimeSeries(oBasin->getSoilMoist_12(),
			      oControl->path_ResultsFolder + "SoilMoistureUp.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Soil_Water_Content_L1){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "SoilMoistureL1.tab");
    oReport->ReportTimeSeries(oBasin->getSoilMoist1(),
			      oControl->path_ResultsFolder + "SoilMoistureL1.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Soil_Water_Content_L2){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "SoilMoistureL2.tab");
    oReport->ReportTimeSeries(oBasin->getSoilMoist2(),
			      oControl->path_ResultsFolder + "SoilMoistureL2.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Soil_Water_Content_L3){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "SoilMoistureL3.tab");
    oReport->ReportTimeSeries(oBasin->getSoilMoist3(),
			      oControl->path_ResultsFolder + "SoilMoistureL3.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_WaterTableDepth){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "WaterTableDepth.tab");
    oReport->ReportTimeSeries(oBasin->getWaterTableDepth(),
			      oControl->path_ResultsFolder + "WaterTableDepth.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Soil_Sat_Deficit){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "SoilSatDef.tab");
    oReport->ReportTimeSeries(oBasin->getSaturationDeficit(),
			      oControl->path_ResultsFolder + "SoilSatDef.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_GroundWater) {
    if (oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "GroundWater.tab");
    oReport->ReportTimeSeries(oBasin->getGrndWater(),
			      oControl->path_ResultsFolder + "GroundWater.tab",
			      oControl->current_ts_count);
  }
  if(oControl->sw_deepGW){
    if (oControl->Rep_DeepGWater){
      if (oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "DeepGroundWater.tab");
      oReport->ReportTimeSeries(oBasin->getGrndWater(),
				oControl->path_ResultsFolder + "DeepGroundWater.tab",
				oControl->current_ts_count);
    }
  }
  if (oControl->RepTs_Soil_Net_Rad){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "NetRadS.tab");
    oReport->ReportTimeSeries(oBasin->getNetRad(),
			      oControl->path_ResultsFolder + "NetRadS.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Net_Rad_sum){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "NetRadtot.tab");
    oReport->ReportTimeSeries(oBasin->getNetRad_sum(),
			      oControl->path_ResultsFolder + "NetRadtot.tab",
			      oControl->current_ts_count);
  }

  if (oControl->RepTs_Soil_LE){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "LatHeat.tab");
    oReport->ReportTimeSeries(oBasin->getLatheat(),
			      oControl->path_ResultsFolder + "LatHeat.tab",
			      oControl->current_ts_count);
  }

  if (oControl->RepTs_LE_sum){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "LatHeattot.tab");
    oReport->ReportTimeSeries(oBasin->getLatheat_sum(),
			      oControl->path_ResultsFolder + "LatHeattot.tab",
			      oControl->current_ts_count);
  }

  if (oControl->RepTs_Sens_Heat){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "SensHeat.tab");
    oReport->ReportTimeSeries(oBasin->getSensHeat(),
			      oControl->path_ResultsFolder + "SensHeat.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Grnd_Heat){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "GrndHeat.tab");
    oReport->ReportTimeSeries(oBasin->getGrndHeat(),
			      oControl->path_ResultsFolder + "GrndHeat.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Snow_Heat){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "SnowHeat.tab");
    oReport->ReportTimeSeries(oBasin->getSnwHeat(),
			      oControl->path_ResultsFolder + "SnowHeat.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Soil_Temperature){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "SoilTemp.tab");
    oReport->ReportTimeSeries(oBasin->getSoilTemp(),
			      oControl->path_ResultsFolder + "SoilTemp.tab",
			      oControl->current_ts_count);
  }
	
  if (oControl->RepTs_Skin_Temperature){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "SkinTemp.tab");
    oReport->ReportTimeSeries(oBasin->getSkinTemp(),
			      oControl->path_ResultsFolder + "SkinTemp.tab",
			      oControl->current_ts_count);
  }



  if(oControl->toggle_chan_evap > 0){
    if (oControl->RepTs_ChanEvap){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "EvapC.tab");
      oReport->ReportTimeSeries(oBasin->getChanEvap(),
				oControl->path_ResultsFolder + "EvapC.tab",
				oControl->current_ts_count);
    }
    if(oControl->toggle_chan_evap == 1){
      if (oControl->RepTs_Water_Temperature){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "WaterTemp.tab");
	oReport->ReportTimeSeries(oBasin->getTemp_w(),
				  oControl->path_ResultsFolder + "WaterTemp.tab",
				  oControl->current_ts_count);
      }
    }
  }

  if (oControl->RepTs_Total_ET){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "Evap.tab");
    oReport->ReportTimeSeries(oBasin->getEvaporation(),
			      oControl->path_ResultsFolder + "Evap.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Transpiration_sum){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "EvapT.tab");
    oReport->ReportTimeSeries(oBasin->getTranspiration_all(),
			      oControl->path_ResultsFolder + "EvapT.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Einterception_sum){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "EvapI.tab");
    oReport->ReportTimeSeries(oBasin->getEvaporationI_all(),
			      oControl->path_ResultsFolder + "EvapI.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Esoil_sum){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "EvapS.tab");
    oReport->ReportTimeSeries(oBasin->getEvaporationS_all(),
			      oControl->path_ResultsFolder + "EvapS.tab",
			      oControl->current_ts_count);
  }

  // Some vertical - lateral fluxes
  if (oControl->RepTs_GWtoChn){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "GWtoChn.tab");
    oReport->ReportTimeSeries(oBasin->getFluxGWtoChn(),
			      oControl->path_ResultsFolder + "GWtoChn.tab",
			      oControl->current_ts_count);
  }
  if(oControl->sw_deepGW){
    if (oControl->RepTs_DeepGWtoChn){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "DeepGWtoChn.tab");
      oReport->ReportTimeSeries(oBasin->getFluxDeepGWtoChn(),
				oControl->path_ResultsFolder + "DeepGWtoChn.tab",
				oControl->current_ts_count);
    }
  }
  if (oControl->RepTs_SrftoChn){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "SrftoChn.tab");
    oReport->ReportTimeSeries(oBasin->getFluxSrftoChn(),
			      oControl->path_ResultsFolder + "SrftoChn.tab",
			      oControl->current_ts_count);
  }	
  if (oControl->RepTs_Infilt){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "Infilt.tab");
    oReport->ReportTimeSeries(oBasin->getFluxInfilt(),
			      oControl->path_ResultsFolder + "Infilt.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Exfilt){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "ReturnSrf.tab");
    oReport->ReportTimeSeries(oBasin->getFluxExfilt(),
			      oControl->path_ResultsFolder + "ReturnSrf.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_PercolL2){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "PercolL2.tab");
    oReport->ReportTimeSeries(oBasin->getFluxPercolL2(),
			      oControl->path_ResultsFolder + "PercolL2.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_ReturnL1){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "ReturnL1.tab");
    oReport->ReportTimeSeries(oBasin->getFluxL2toL1(),
			      oControl->path_ResultsFolder + "ReturnL1.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_PercolL3){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "PercolL3.tab");
    oReport->ReportTimeSeries(oBasin->getFluxPercolL3(),
			      oControl->path_ResultsFolder + "PercolL3.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Recharge){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "Recharge.tab");
    oReport->ReportTimeSeries(oBasin->getFluxRecharge(),
			      oControl->path_ResultsFolder + "Recharge.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_ReturnL2){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "ReturnL2.tab");
    oReport->ReportTimeSeries(oBasin->getFluxL3toL2(),
			      oControl->path_ResultsFolder + "ReturnL2.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_Leak){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "Leak.tab");
    oReport->ReportTimeSeries(oBasin->getBedrockLeakage(),
			      oControl->path_ResultsFolder + "Leak.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_LattoSrf){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "SrfLatI.tab");
    oReport->ReportTimeSeries(oBasin->getFluxLattoSrf(),
			      oControl->path_ResultsFolder + "SrfLatI.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_LattoChn){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "ChnLatI.tab");
    oReport->ReportTimeSeries(oBasin->getFluxLattoChn(),
			      oControl->path_ResultsFolder + "ChnLatI.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_LattoGW){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "GWLatI.tab");
    oReport->ReportTimeSeries(oBasin->getFluxLattoGW(),
			      oControl->path_ResultsFolder + "GWLatI.tab",
			      oControl->current_ts_count);
  }
  if(oControl->sw_deepGW){
    if (oControl->RepTs_LattoDeepGW){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "DeepGWLatI.tab");
      oReport->ReportTimeSeries(oBasin->getFluxLattoDeepGW(),
				oControl->path_ResultsFolder + "DeepGWLatI.tab",
				oControl->current_ts_count);
    }
  }
  if (oControl->RepTs_ChntoLat){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "ChnLatO.tab");
    oReport->ReportTimeSeries(oBasin->getFluxChntoLat(),
            oControl->path_ResultsFolder + "ChnLatO.tab",
            oControl->current_ts_count);
  }
  if (oControl->RepTs_SrftoLat){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "SrfLatO.tab");
    oReport->ReportTimeSeries(oBasin->getFluxSrftoLat(),
			      oControl->path_ResultsFolder + "SrfLatO.tab",
			      oControl->current_ts_count);
  }
  if (oControl->RepTs_GWtoLat){
    if(oControl->GetTimeStep() <= oControl->report_times)
      oReport->RenameFile(oControl->path_ResultsFolder + "GWLatO.tab");
    oReport->ReportTimeSeries(oBasin->getFluxGWtoLat(),
			      oControl->path_ResultsFolder + "GWLatO.tab",
			      oControl->current_ts_count);
  }
  if(oControl->sw_deepGW){
    if (oControl->RepTs_DeepGWtoLat){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "DeepGWLatO.tab");
      oReport->ReportTimeSeries(oBasin->getFluxDeepGWtoLat(),
				oControl->path_ResultsFolder + "DeepGWLatO.tab",
				oControl->current_ts_count);
    }
  }
  // Tracking	
  // two-pore domain
  if(oControl->sw_trck and oControl->sw_TPD){
    if (oControl->RepTs_Moist_MW1){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "SWC_MW1.tab");
      oReport->ReportTimeSeries(oBasin->getMoistureMW1(),
				oControl->path_ResultsFolder + "SWC_MW1.tab",
				oControl->current_ts_count);
    }
    if (oControl->RepTs_Moist_MW2){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "SWC_MW2.tab");
      oReport->ReportTimeSeries(oBasin->getMoistureMW2(),
				oControl->path_ResultsFolder + "SWC_MW2.tab",
				oControl->current_ts_count);
    }
    if(oControl->RepTs_Frac_MW1){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "fracMW_L1.tab");
      oReport->ReportTimeSeries(oBasin->getFracMW1(),
				oControl->path_ResultsFolder + "fracMW_L1.tab",
				oControl->current_ts_count);
    }
    if(oControl->RepTs_Frac_MW2){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "fracMW_L2.tab");
      oReport->ReportTimeSeries(oBasin->getFracMW2(),
				oControl->path_ResultsFolder + "fracMW_L2.tab",
				oControl->current_ts_count);
    }
    if(oControl->RepTs_Frac_MW12){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "fracMW_Up.tab");
      oReport->ReportTimeSeries(oBasin->getFracMW12(),
				oControl->path_ResultsFolder + "fracMW_Up.tab",
				oControl->current_ts_count);
    }
  }
  //-----------------------------------------------------------------------------

  if(oControl->sw_trck && oControl->sw_2H){
    if (oControl->RepTs_d2Hprecip){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d2H_precip.tab");
      oReport->ReportTimeSeries(oAtmosphere->getd2Hprecip(),
				oControl->path_ResultsFolder + "d2H_precip.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_d2Hcanopy_sum){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d2H_canopy.tab");
      oReport->ReportTimeSeries(oTracking->getd2Hcanopy_sum(),
				oControl->path_ResultsFolder + "d2H_canopy.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_d2Hsnowpack){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d2H_snowpack.tab");
      oReport->ReportTimeSeries(oTracking->getd2Hsnowpack(),
				oControl->path_ResultsFolder + "d2H_snowpack.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_d2Hsurface){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d2H_surface.tab");
      oReport->ReportTimeSeries(oTracking->getd2Hsurface(),
				oControl->path_ResultsFolder + "d2H_surface.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_d2Hchan){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d2H_chan.tab");
      oReport->ReportTimeSeries(oTracking->getd2Hchan(),
				oControl->path_ResultsFolder + "d2H_chan.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_d2Hsoil1){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d2H_soilL1.tab");
      oReport->ReportTimeSeries(oTracking->getd2Hsoil1(),
				oControl->path_ResultsFolder + "d2H_soilL1.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_d2Hsoil2){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d2H_soilL2.tab");
      oReport->ReportTimeSeries(oTracking->getd2Hsoil2(),
				oControl->path_ResultsFolder + "d2H_soilL2.tab",
				oControl->current_ts_count);
    }
    if (oControl->RepTs_d2HsoilUp){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d2H_soilUp.tab");
      oReport->ReportTimeSeries(oTracking->getd2Hsoil_12(),
				oControl->path_ResultsFolder + "d2H_soilUp.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_d2Hsoil3){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d2H_soilL3.tab");
      oReport->ReportTimeSeries(oTracking->getd2Hsoil3(),
				oControl->path_ResultsFolder + "d2H_soilL3.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_d2HsoilAv){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d2H_soilAv.tab");
      oReport->ReportTimeSeries(oTracking->getd2Hsoil_Av(),
				oControl->path_ResultsFolder + "d2H_soilAv.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_d2Hgroundwater){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d2H_groundwater.tab");
      oReport->ReportTimeSeries(oTracking->getd2Hgroundwater(),
				oControl->path_ResultsFolder + "d2H_groundwater.tab",
				oControl->current_ts_count);
    }
    if(oControl->sw_deepGW){
      if (oControl->RepTs_d2HdeepGW){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "d2H_DeepGW.tab");
	oReport->ReportTimeSeries(oTracking->getd2HdeepGW(),
				  oControl->path_ResultsFolder + "d2H_DeepGW.tab",
				  oControl->current_ts_count);
      }
      if (oControl->RepTs_d2HdeepGWQ){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "d2H_DeepGWQ.tab");
	oReport->ReportTimeSeries(oTracking->getd2HdeepGWQ(),
				  oControl->path_ResultsFolder + "d2H_DeepGWQ.tab",
				  oControl->current_ts_count);
      }
    }
    if (oControl->RepTs_d2Hleakage){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d2H_leakage.tab");
      oReport->ReportTimeSeries(oTracking->getd2Hleakage(),
				oControl->path_ResultsFolder + "d2H_leakage.tab",
				oControl->current_ts_count);
    }
    if (oControl->RepTs_d2HevapS_sum){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d2H_evapS.tab");
      oReport->ReportTimeSeries(oTracking->getd2HevapS_sum(),
				oControl->path_ResultsFolder + "d2H_evapS.tab",
				oControl->current_ts_count);
    }
    if (oControl->RepTs_d2HevapI_sum){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d2H_evapI.tab");
      oReport->ReportTimeSeries(oTracking->getd2HevapI_sum(),
				oControl->path_ResultsFolder + "d2H_evapI.tab",
				oControl->current_ts_count);
    }
    if (oControl->RepTs_d2HevapT_sum){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d2H_evapT.tab");
      oReport->ReportTimeSeries(oTracking->getd2HevapT_sum(),
				oControl->path_ResultsFolder + "d2H_evapT.tab",
				oControl->current_ts_count);
    }
    if(oControl->sw_TPD){
      if (oControl->RepTs_d2H_MW1){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "d2H_MW1.tab");
	oReport->ReportTimeSeries(oTracking->getd2H_MW1(),
				  oControl->path_ResultsFolder + "d2H_MW1.tab",
				  oControl->current_ts_count);
      }
      if (oControl->RepTs_d2H_MW2){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "d2H_MW2.tab");
	oReport->ReportTimeSeries(oTracking->getd2H_MW2(),
				  oControl->path_ResultsFolder + "d2H_MW2.tab",
				  oControl->current_ts_count);
      }
      if (oControl->RepTs_d2H_TB1){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "d2H_TB1.tab");
	oReport->ReportTimeSeries(oTracking->getd2H_TB1(),
				  oControl->path_ResultsFolder + "d2H_TB1.tab",
				  oControl->current_ts_count);
      }
      if (oControl->RepTs_d2H_TB2){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "d2H_TB2.tab");
	oReport->ReportTimeSeries(oTracking->getd2H_TB2(),
				  oControl->path_ResultsFolder + "d2H_TB2.tab",
				  oControl->current_ts_count);
      }
    }

  }

  if(oControl->sw_trck && oControl->sw_18O){
    if (oControl->RepTs_d18Oprecip){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d18O_precip.tab");
      oReport->ReportTimeSeries(oAtmosphere->getd18Oprecip(),
				oControl->path_ResultsFolder + "d18O_precip.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_d18Ocanopy_sum){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d18O_canopy.tab");
      oReport->ReportTimeSeries(oTracking->getd18Ocanopy_sum(),
				oControl->path_ResultsFolder + "d18O_canopy.tab",
				oControl->current_ts_count);
    }
    if (oControl->RepTs_d18Osnowpack){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d18O_snowpack.tab");
      oReport->ReportTimeSeries(oTracking->getd18Osnowpack(),
				oControl->path_ResultsFolder + "d18O_snowpack.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_d18Osurface){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d18O_surface.tab");
      oReport->ReportTimeSeries(oTracking->getd18Osurface(),
				oControl->path_ResultsFolder + "d18O_surface.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_d18Ochan){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d18O_chan.tab");
      oReport->ReportTimeSeries(oTracking->getd18Ochan(),
				oControl->path_ResultsFolder + "d18O_chan.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_d18Osoil1){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d18O_soilL1.tab");
      oReport->ReportTimeSeries(oTracking->getd18Osoil1(),
				oControl->path_ResultsFolder + "d18O_soilL1.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_d18Osoil2){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d18O_soilL2.tab");
      oReport->ReportTimeSeries(oTracking->getd18Osoil2(),
				oControl->path_ResultsFolder + "d18O_soilL2.tab",
				oControl->current_ts_count);
    }
    if (oControl->RepTs_d18OsoilUp){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d18O_soilUp.tab");
      oReport->ReportTimeSeries(oTracking->getd18Osoil_12(),
				oControl->path_ResultsFolder + "d18O_soilUp.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_d18Osoil3){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d18O_soilL3.tab");
      oReport->ReportTimeSeries(oTracking->getd18Osoil3(),
				oControl->path_ResultsFolder + "d18O_soilL3.tab",
				oControl->current_ts_count);
    }
		
    if (oControl->RepTs_d18OsoilAv){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d18O_soilAv.tab");
      oReport->ReportTimeSeries(oTracking->getd18Osoil_Av(),
				oControl->path_ResultsFolder + "d18O_soilAv.tab",
				oControl->current_ts_count);
    }
		
    if (oControl->RepTs_d18Ogroundwater){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d18O_groundwater.tab");
      oReport->ReportTimeSeries(oTracking->getd18Ogroundwater(),
				oControl->path_ResultsFolder + "d18O_groundwater.tab",
				oControl->current_ts_count);
    }
    if(oControl->sw_deepGW){
      if (oControl->RepTs_d18OdeepGW){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "d18O_DeepGW.tab");
	oReport->ReportTimeSeries(oTracking->getd18OdeepGW(),
				  oControl->path_ResultsFolder + "d18O_DeepGW.tab",
				  oControl->current_ts_count);
      }
      if (oControl->RepTs_d18OdeepGWQ){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "d18O_DeepGWQ.tab");
	oReport->ReportTimeSeries(oTracking->getd18OdeepGWQ(),
				  oControl->path_ResultsFolder + "d18O_DeepGWQ.tab",
				  oControl->current_ts_count);
      }
    }
    if (oControl->RepTs_d18Oleakage){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d18O_leakage.tab");
      oReport->ReportTimeSeries(oTracking->getd18Oleakage(),
				oControl->path_ResultsFolder + "d18O_leakage.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_d18OevapS_sum){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d18O_evapS.tab");
      oReport->ReportTimeSeries(oTracking->getd18OevapS_sum(),
				oControl->path_ResultsFolder + "d18O_evapS.tab",
				oControl->current_ts_count);
    }
    if (oControl->RepTs_d18OevapI_sum){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d18O_evapI.tab");
      oReport->ReportTimeSeries(oTracking->getd18OevapI_sum(),
				oControl->path_ResultsFolder + "d18O_evapI.tab",
				oControl->current_ts_count);
    }
    if (oControl->RepTs_d18OevapT_sum){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "d18O_evapT.tab");
      oReport->ReportTimeSeries(oTracking->getd18OevapT_sum(),
				oControl->path_ResultsFolder + "d18O_evapT.tab",
				oControl->current_ts_count);
    }
		
    if(oControl->sw_TPD){
      if (oControl->RepTs_d18O_MW1){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "d18O_MW1.tab");
	oReport->ReportTimeSeries(oTracking->getd18O_MW1(),
				  oControl->path_ResultsFolder + "d18O_MW1.tab",
				  oControl->current_ts_count);
      }
      if (oControl->RepTs_d18O_MW2){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "d18O_MW2.tab");
	oReport->ReportTimeSeries(oTracking->getd18O_MW2(),
				  oControl->path_ResultsFolder + "d18O_MW2.tab",
				  oControl->current_ts_count);
      }
      if (oControl->RepTs_d18O_TB1){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "d18O_TB1.tab");
	oReport->ReportTimeSeries(oTracking->getd18O_TB1(),
				  oControl->path_ResultsFolder + "d18O_TB1.tab",
				  oControl->current_ts_count);
      }
      if (oControl->RepTs_d18O_TB2){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "d18O_TB2.tab");
	oReport->ReportTimeSeries(oTracking->getd18O_TB2(),
				  oControl->path_ResultsFolder + "d18O_TB2.tab",
				  oControl->current_ts_count);
      }
    }
  }

  if(oControl->sw_trck && oControl->sw_Age){

    if (oControl->RepTs_Agecanopy_sum){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "Age_canopy.tab");
      oReport->ReportTimeSeries(oTracking->getAgecanopy_sum(),
				oControl->path_ResultsFolder + "Age_canopy.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_Agesnowpack){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "Age_snowpack.tab");
      oReport->ReportTimeSeries(oTracking->getAgesnowpack(),
				oControl->path_ResultsFolder + "Age_snowpack.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_Agesurface){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "Age_surface.tab");
      oReport->ReportTimeSeries(oTracking->getAgesurface(),
				oControl->path_ResultsFolder + "Age_surface.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_Agechan){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "Age_chan.tab");
      oReport->ReportTimeSeries(oTracking->getAgechan(),
				oControl->path_ResultsFolder + "Age_chan.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_Agesoil1){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "Age_soilL1.tab");
      oReport->ReportTimeSeries(oTracking->getAgesoil1(),
				oControl->path_ResultsFolder + "Age_soilL1.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_Agesoil2){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "Age_soilL2.tab");
      oReport->ReportTimeSeries(oTracking->getAgesoil2(),
				oControl->path_ResultsFolder + "Age_soilL2.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_AgesoilUp){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "Age_soilUp.tab");
      oReport->ReportTimeSeries(oTracking->getAgesoil_12(),
				oControl->path_ResultsFolder + "Age_soilUp.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_Agesoil3){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "Age_soilL3.tab");
      oReport->ReportTimeSeries(oTracking->getAgesoil3(),
				oControl->path_ResultsFolder + "Age_soilL3.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_AgesoilAv){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "Age_soilAv.tab");
      oReport->ReportTimeSeries(oTracking->getAgesoil_Av(),
				oControl->path_ResultsFolder + "Age_soilAv.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_Agegroundwater){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "Age_groundwater.tab");
      oReport->ReportTimeSeries(oTracking->getAgegroundwater(),
				oControl->path_ResultsFolder + "Age_groundwater.tab",
				oControl->current_ts_count);
    }
    if(oControl->sw_deepGW){
      if (oControl->RepTs_AgedeepGW){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "Age_DeepGW.tab");
	oReport->ReportTimeSeries(oTracking->getAgedeepGW(),
				  oControl->path_ResultsFolder + "Age_DeepGW.tab",
				  oControl->current_ts_count);
      }
      if (oControl->RepTs_AgedeepGWQ){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "Age_DeepGWQ.tab");
	oReport->ReportTimeSeries(oTracking->getAgedeepGWQ(),
				  oControl->path_ResultsFolder + "Age_DeepGWQ.tab",
				  oControl->current_ts_count);
      }
    }
    if (oControl->RepTs_Ageleakage){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "Age_leakage.tab");
      oReport->ReportTimeSeries(oTracking->getAgeleakage(),
				oControl->path_ResultsFolder + "Age_leakage.tab",
				oControl->current_ts_count);
    }

    if (oControl->RepTs_AgeevapS_sum){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "Age_evapS.tab");
      oReport->ReportTimeSeries(oTracking->getAgeevapS_sum(),
				oControl->path_ResultsFolder + "Age_evapS.tab",
				oControl->current_ts_count);
    }
    if (oControl->RepTs_AgeevapI_sum){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "Age_evapI.tab");
      oReport->ReportTimeSeries(oTracking->getAgeevapI_sum(),
				oControl->path_ResultsFolder + "Age_evapI.tab",
				oControl->current_ts_count);
    }
    if (oControl->RepTs_AgeevapT_sum){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "Age_evapT.tab");
      oReport->ReportTimeSeries(oTracking->getAgeevapT_sum(),
				oControl->path_ResultsFolder + "Age_evapT.tab",
				oControl->current_ts_count);
    }
    if (oControl->RepTs_AgeGWtoChn){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "Age_GwtoChn.tab");
      oReport->ReportTimeSeries(oTracking->getAgeGWtoChn(),
				oControl->path_ResultsFolder + "Age_GWtoChn.tab",
				oControl->current_ts_count);
    }
    if (oControl->RepTs_AgeSrftoChn){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "Age_SrftoChn.tab");
      oReport->ReportTimeSeries(oTracking->getAgeSrftoChn(),
				oControl->path_ResultsFolder + "Age_SrftoChn.tab",
				oControl->current_ts_count);
    }
    if (oControl->RepTs_AgeRecharge){
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(oControl->path_ResultsFolder + "Age_Recharge.tab");
      oReport->ReportTimeSeries(oTracking->getAgeSrftoChn(),
				oControl->path_ResultsFolder + "Age_Recharge.tab",
				oControl->current_ts_count);
    }


    if(oControl->sw_TPD){
      if (oControl->RepTs_Age_MW1){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "Age_MW1.tab");
	oReport->ReportTimeSeries(oTracking->getAge_MW1(),
				  oControl->path_ResultsFolder + "Age_MW1.tab",
				  oControl->current_ts_count);
      }
      if (oControl->RepTs_Age_MW2){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "Age_MW2.tab");
	oReport->ReportTimeSeries(oTracking->getAge_MW2(),
				  oControl->path_ResultsFolder + "Age_MW2.tab",
				  oControl->current_ts_count);
      }
      if (oControl->RepTs_Age_MWUp){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "Age_MW12.tab");
	oReport->ReportTimeSeries(oTracking->getAge_MW12(),
				  oControl->path_ResultsFolder + "Age_MW12.tab",
				  oControl->current_ts_count);
      }
      if (oControl->RepTs_Age_TB1){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "Age_TB1.tab");
	oReport->ReportTimeSeries(oTracking->getAge_TB1(),
				  oControl->path_ResultsFolder + "Age_TB1.tab",
				  oControl->current_ts_count);
      }
      if (oControl->RepTs_Age_TB2){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "Age_TB2.tab");
	oReport->ReportTimeSeries(oTracking->getAge_TB2(),
				  oControl->path_ResultsFolder + "Age_TB2.tab",
				  oControl->current_ts_count);
      }
      if (oControl->RepTs_Age_TBUp){
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(oControl->path_ResultsFolder + "Age_TB12.tab");
	oReport->ReportTimeSeries(oTracking->getAge_TB12(),
				  oControl->path_ResultsFolder + "Age_TB12.tab",
				  oControl->current_ts_count);
      }

    }

  }

  for (int i = 0; i < oControl->NumSpecs; i++) {
    stringstream name;

    if (oControl->RepTs_Veget_frac) {
      name << oControl->path_ResultsFolder << "p_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getVegetFrac(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->RepTs_Stem_Density) {
      name << oControl->path_ResultsFolder << "num_of_tress_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getStemDensity(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->RepTs_RootFrac1Species) {
      name << oControl->path_ResultsFolder << "RootFrac.L1_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getRootFrac1(i), name.str() , oControl->current_ts_count);
      name.str("");
    }
    if (oControl->RepTs_RootFrac2Species) {
      name << oControl->path_ResultsFolder << "RootFrac.L2_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getRootFrac2(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->RepTs_Leaf_Area_Index) {
      name << oControl->path_ResultsFolder << "lai_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getLAI(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->RepTs_Canopy_Conductance) {
      name << oControl->path_ResultsFolder << "CanopyConduct_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getCanopyCond(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->RepTs_GPP) {
      name << oControl->path_ResultsFolder << "GPP_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getGPP(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->RepTs_NPP) {
      name << oControl->path_ResultsFolder << "NPP_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getNPP(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->RepTs_Basal_Area) {
      name << oControl->path_ResultsFolder << "BasalArea_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getBasalArea(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->RepTs_Tree_Height) {
      name << oControl->path_ResultsFolder << "ThreeHeight_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getTreeHeight(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->RepTs_Root_Mass) {
      name << oControl->path_ResultsFolder << "RootMass_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getRootMass(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->RepTs_Canopy_Temp) {
      name << oControl->path_ResultsFolder << "CanopyTemp_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getCanopyTemp(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->RepTs_Canopy_NetR) {
      name << oControl->path_ResultsFolder << "NetRadC_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getCanopyNetRad(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->RepTs_Canopy_LE_E) {
      name << oControl->path_ResultsFolder << "CanopyLatHeatEi_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getCanopyLatHeatE(i), name.str() , oControl->current_ts_count);
      name.str("");
    }
    if (oControl->RepTs_Canopy_LE_T) {
      name << oControl->path_ResultsFolder << "CanopyLatHeatTr_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getCanopyLatHeatT(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->RepTs_Canopy_Sens_Heat) {
      name << oControl->path_ResultsFolder << "CanopySensHeat_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getCanopySensHeat(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->RepTs_Canopy_Water_Stor) {
      name << oControl->path_ResultsFolder << "CanopyWaterStor_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getCanopyWaterStor(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    if (oControl->RepTs_ETspecies) {
      name << oControl->path_ResultsFolder << "ETc_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getETspecies(i), name.str() , oControl->current_ts_count);
      name.str("");
    }
    if (oControl->RepTs_Transpiration) {
      name << oControl->path_ResultsFolder << "EvapT_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getTranspiration(i), name.str() , oControl->current_ts_count);
      name.str("");
    }
    if (oControl->RepTs_Einterception) {
      name << oControl->path_ResultsFolder << "EvapI_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getEinterception(i), name.str() , oControl->current_ts_count);
      name.str("");
    }
    if (oControl->RepTs_Esoil) {
      name << oControl->path_ResultsFolder << "EvapS_" << i << ".tab";
      if(oControl->GetTimeStep() <= oControl->report_times)
	oReport->RenameFile(name.str());
      oReport->ReportTimeSeries(oBasin->getEsoil(i), name.str() , oControl->current_ts_count);
      name.str("");
    }

    // Tracking
    if(oControl->sw_trck){
      if (oControl->sw_2H && oControl->RepTs_d2Hcanopy) {
	name << oControl->path_ResultsFolder << "d2Hcanopy_" << i << ".tab";
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(name.str());
	oReport->ReportTimeSeries(oBasin->getd2Hcanopy(i), name.str() , oControl->current_ts_count);
	name.str("");
      }
      
      if (oControl->sw_18O && oControl->RepTs_d18Ocanopy) {
	name << oControl->path_ResultsFolder << "d18Ocanopy_" << i << ".tab";
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(name.str());
	oReport->ReportTimeSeries(oBasin->getd18Ocanopy(i), name.str() , oControl->current_ts_count);
	name.str("");
      }
      
      if (oControl->sw_Age && oControl->RepTs_Agecanopy) {
	name << oControl->path_ResultsFolder << "Agecanopy_" << i << ".tab";
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(name.str());
	oReport->ReportTimeSeries(oBasin->getAgecanopy(i), name.str() , oControl->current_ts_count);
	name.str("");
      }

      if (oControl->sw_2H && oControl->RepTs_d2HevapI) {
	name << oControl->path_ResultsFolder << "d2HevapI_" << i << ".tab";
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(name.str());
	oReport->ReportTimeSeries(oBasin->getd2HevapI(i), name.str() , oControl->current_ts_count);
	name.str("");
      }

      if (oControl->sw_18O && oControl->RepTs_d18OevapI) {
	name << oControl->path_ResultsFolder << "d18OevapI_" << i << ".tab";
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(name.str());
	oReport->ReportTimeSeries(oBasin->getd18OevapI(i), name.str() , oControl->current_ts_count);
	name.str("");
      }

      if (oControl->sw_Age && oControl->RepTs_AgeevapI) {
	name << oControl->path_ResultsFolder << "AgeevapI_" << i << ".tab";
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(name.str());
	oReport->ReportTimeSeries(oBasin->getAgeevapI(i), name.str() , oControl->current_ts_count);
	name.str("");
      }
      if (oControl->sw_2H && oControl->RepTs_d2HevapT) {
	name << oControl->path_ResultsFolder << "d2HevapT_" << i << ".tab";
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(name.str());
	oReport->ReportTimeSeries(oBasin->getd2HevapT(i), name.str() , oControl->current_ts_count);
	name.str("");
      }

      if (oControl->sw_18O && oControl->RepTs_d18OevapT) {
	name << oControl->path_ResultsFolder << "d18OevapT_" << i << ".tab";
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(name.str());
	oReport->ReportTimeSeries(oBasin->getd18OevapT(i), name.str() , oControl->current_ts_count);
	name.str("");
      }

      if (oControl->sw_Age && oControl->RepTs_AgeevapT) {
	name << oControl->path_ResultsFolder << "AgeevapT_" << i << ".tab";
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(name.str());
	oReport->ReportTimeSeries(oBasin->getAgeevapT(i), name.str() , oControl->current_ts_count);
	name.str("");
      }
      if (oControl->sw_2H && oControl->RepTs_d2HevapS) {
	name << oControl->path_ResultsFolder << "d2HevapS_" << i << ".tab";
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(name.str());
	oReport->ReportTimeSeries(oBasin->getd2HevapS(i), name.str() , oControl->current_ts_count);
	name.str("");
      }

      if (oControl->sw_18O && oControl->RepTs_d18OevapS) {
	name << oControl->path_ResultsFolder << "d18OevapS_" << i << ".tab";
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(name.str());
	oReport->ReportTimeSeries(oBasin->getd18OevapS(i), name.str() , oControl->current_ts_count);
	name.str("");
      }

      if (oControl->sw_Age && oControl->RepTs_AgeevapS) {
	name << oControl->path_ResultsFolder << "AgeevapS_" << i << ".tab";
	if(oControl->GetTimeStep() <= oControl->report_times)
	  oReport->RenameFile(name.str());
	oReport->ReportTimeSeries(oBasin->getAgeevapS(i), name.str() , oControl->current_ts_count);
	name.str("");
      }
    }
  }

  return EXIT_SUCCESS;
}
