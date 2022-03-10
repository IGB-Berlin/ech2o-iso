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
 * Report2nc.cpp
 *
 *  Created on: Feb 10, 2021
 *      Author: Xiaoqiang Yang
 *      This script stores modeling state variables and fluxes as netcdf4 format file
 */

#include <netcdf>
#include <iostream>
#include <string>
#include "Sativa.h"
using namespace std;


int Report2nc(){
  
  //create the nc file for water fluxes
  if(oControl->current_ts_count == 1){
    oReport->CreatOutputNC(oControl->path_ResultsFolder, "W");
  }  
  if(oControl->Rep_Long_Rad_Down)
    oReport->UpdateOutputNC(oAtmosphere->getIncomingLongWave(), "Ldown","W");
  if(oControl->Rep_Short_Rad_Down)
    oReport->UpdateOutputNC(oAtmosphere->getIncomingShortWave(), "Sdown", "W");
  if(oControl->Rep_Precip)
	oReport->UpdateOutputNC(oAtmosphere->getPrecipitation(), "Pp", "W");
  if(oControl->Rep_Rel_Humidity)
	oReport->UpdateOutputNC(oAtmosphere->getRelativeHumidty(), "RH", "W");
  if(oControl->Rep_Wind_Speed)
	oReport->UpdateOutputNC(oAtmosphere->getWindSpeed(), "WndSp", "W");
  if(oControl->Rep_AvgAir_Temperature)
	oReport->UpdateOutputNC(oAtmosphere->getTemperature(), "Tp", "W");
  if(oControl->Rep_MinAir_Temperature)
	oReport->UpdateOutputNC(oAtmosphere->getMinTemperature(), "TpMin", "W");
  if(oControl->Rep_MaxAir_Temperature)
	oReport->UpdateOutputNC(oAtmosphere->getMaxTemperature(), "TpMax", "W");
  if(oControl->sw_BC){
    if(oControl->Rep_BC_surface)
      oReport->UpdateOutputNC(oBasin->getBCsurface(), "BCs", "W");
    if(oControl->Rep_BC_groundwater)
      oReport->UpdateOutputNC(oBasin->getBCgroundwater(), "BCgw", "W");      
  }

  if (oControl->Rep_Field_Capacity_L1)
	oReport->UpdateOutputNC(oBasin->getFieldCapacityL1(), "FCap1", "W");
  if (oControl->Rep_Field_Capacity_L2)
	oReport->UpdateOutputNC(oBasin->getFieldCapacityL2(), "FCap2", "W");
  if (oControl->Rep_Field_Capacity_L3)
	oReport->UpdateOutputNC(oBasin->getFieldCapacityL3(), "FCap3", "W");

  if (oControl->Rep_Canopy_Water_Stor_sum)
	oReport->UpdateOutputNC(oBasin->getCanopyStorage(), "Cs", "W");
  if (oControl->Rep_SWE)
	oReport->UpdateOutputNC(oBasin->getSnowWaterEquiv(), "SWE", "W");
  if (oControl->Rep_Infilt_Cap)
	oReport->UpdateOutputNC(oBasin->getInfiltCap(), "IfCap",  "W");
  if (oControl->Rep_Streamflow)
	oReport->UpdateOutputNC(oBasin->getStreamflow(), "Q", "W");
  if (oControl->Rep_Saturation_Area)
	oReport->UpdateOutputNC(oBasin->getSatArea(), "SatArea_", "W");
  if (oControl->Rep_Ponding)
	oReport->UpdateOutputNC(oBasin->getPondingWater(), "Ponding_", "W");
  if (oControl->Rep_Soil_Water_Content_Average)
	oReport->UpdateOutputNC(oBasin->getSoilMoist_av(), "SWCav", "W");
  if (oControl->Rep_Soil_Water_Content_12)
	oReport->UpdateOutputNC(oBasin->getSoilMoist_12(), "SWCup", "W");
  if (oControl->Rep_Soil_Water_Content_L1)
	oReport->UpdateOutputNC(oBasin->getSoilMoist1(), "SWC1_", "W");
  if (oControl->Rep_Soil_Water_Content_L2)
	oReport->UpdateOutputNC(oBasin->getSoilMoist2(), "SWC2_", "W");
  if (oControl->Rep_Soil_Water_Content_L3)
	oReport->UpdateOutputNC(oBasin->getSoilMoist3(), "SWC3_", "W");
  if (oControl->Rep_WaterTableDepth)
	oReport->UpdateOutputNC(oBasin->getWaterTableDepth(), "WTD_", "W");
  if (oControl->Rep_Soil_Sat_Deficit)
	oReport->UpdateOutputNC(oBasin->getSaturationDeficit(), "SatDef", "W");
  if (oControl->Rep_GWater)
	oReport->UpdateOutputNC(oBasin->getGrndWater(), "GW", "W");
  if(oControl->sw_deepGW){
    if (oControl->Rep_DeepGWater)
      oReport->UpdateOutputNC(oBasin->getDeepGW(), "deepGW", "W");
  }  
  if (oControl->Rep_Soil_Net_Rad)
	oReport->UpdateOutputNC(oBasin->getNetRad(), "NRs", "W");
  if (oControl->Rep_Net_Rad_sum)
	oReport->UpdateOutputNC(oBasin->getNetRad_sum(), "NRtot", "W");
  if (oControl->Rep_Soil_LE)
	oReport->UpdateOutputNC(oBasin->getLatheat(), "LEs", "W");
  if (oControl->Rep_LE_sum)
	oReport->UpdateOutputNC(oBasin->getLatheat_sum(), "LEtot", "W");  
  if (oControl->Rep_Sens_Heat)
	oReport->UpdateOutputNC(oBasin->getSensHeat(), "SensH", "W");
  if (oControl->Rep_Grnd_Heat)
	oReport->UpdateOutputNC(oBasin->getGrndHeat(), "GrndH", "W");
  if (oControl->Rep_Snow_Heat)
	oReport->UpdateOutputNC(oBasin->getSnwHeat(), "SnowH", "W");
  if (oControl->Rep_Soil_Temperature)
	oReport->UpdateOutputNC(oBasin->getSoilTemp(), "Ts", "W");
  if (oControl->Rep_Skin_Temperature)
	oReport->UpdateOutputNC(oBasin->getSkinTemp(), "Tskin", "W");
  if(oControl->toggle_chan_evap > 0){
    if (oControl->Rep_ChanEvap)
      oReport->UpdateOutputNC(oBasin->getChanEvap(), "EvapC", "W");
    if (oControl->toggle_chan_evap == 1){
      if (oControl->Rep_Water_Temperature)
	oReport->UpdateOutputNC(oBasin->getTemp_w(), "Twater", "W");
    }
  }
  if (oControl->Rep_Total_ET)
	oReport->UpdateOutputNC(oBasin->getEvaporation(), "Evap", "W");
  if (oControl->Rep_Transpiration_sum)
	oReport->UpdateOutputNC(oBasin->getTranspiration_all(), "EvapT", "W");
  if (oControl->Rep_Einterception_sum)
	oReport->UpdateOutputNC(oBasin->getEvaporationI_all(), "EvapI", "W");
  if (oControl->Rep_Esoil_sum)
	oReport->UpdateOutputNC(oBasin->getEvaporationS_all(), "EvapS", "W");

  // Some vertical- lateral fluxes
  if (oControl->Rep_GWtoChn)
    oReport->UpdateOutputNC(oBasin->getFluxGWtoChn(), "GWChn", "W");
  if(oControl->sw_deepGW){
    if (oControl->Rep_DeepGWtoChn)
    oReport->UpdateOutputNC(oBasin->getFluxDeepGWtoChn(), "DeepGWChn", "W");
  }
  if (oControl->Rep_SrftoChn)
    oReport->UpdateOutputNC(oBasin->getFluxSrftoChn(), "SrfChn", "W");
  if (oControl->Rep_Infilt)
    oReport->UpdateOutputNC(oBasin->getFluxInfilt(), "Inf", "W");
  if (oControl->Rep_Exfilt)
    oReport->UpdateOutputNC(oBasin->getFluxExfilt(), "RSrf", "W");
  if (oControl->Rep_PercolL2)
    oReport->UpdateOutputNC(oBasin->getFluxPercolL2(), "PrcL2", "W");
  if (oControl->Rep_ReturnL1)
    oReport->UpdateOutputNC(oBasin->getFluxL2toL1(), "RL1", "W");
  if (oControl->Rep_PercolL3)
    oReport->UpdateOutputNC(oBasin->getFluxPercolL3(), "PrcL3", "W");
  if (oControl->Rep_Recharge)
    oReport->UpdateOutputNC(oBasin->getFluxRecharge(), "Rchg", "W");
  if (oControl->Rep_ReturnL2)
    oReport->UpdateOutputNC(oBasin->getFluxL3toL2(), "RL2", "W");
  if (oControl->Rep_Leak)
    oReport->UpdateOutputNC(oBasin->getBedrockLeakage(),"Leak", "W");
  // Internal lateral fluxes  
  if (oControl->Rep_LattoSrf)
    oReport->UpdateOutputNC(oBasin->getFluxLattoSrf(), "LSrfi", "W");
  if (oControl->Rep_LattoChn)
    oReport->UpdateOutputNC(oBasin->getFluxLattoChn(), "LChni", "W");
  if (oControl->Rep_LattoGW)
    oReport->UpdateOutputNC(oBasin->getFluxLattoGW(), "LGWi", "W");
  if(oControl->sw_deepGW){
    if (oControl->Rep_LattoDeepGW)
    oReport->UpdateOutputNC(oBasin->getFluxLattoDeepGW(), "LDeepGWi", "W");
  }  
  if (oControl->Rep_ChntoLat)
    oReport->UpdateOutputNC(oBasin->getFluxChntoLat(), "LChno", "W");
  if (oControl->Rep_SrftoLat)
    oReport->UpdateOutputNC(oBasin->getFluxSrftoLat(), "LSrfo", "W");
  if (oControl->Rep_GWtoLat)
    oReport->UpdateOutputNC(oBasin->getFluxGWtoLat(), "LGWo", "W");
  if(oControl->sw_deepGW){
    if (oControl->Rep_DeepGWtoLat)
      oReport->UpdateOutputNC(oBasin->getFluxDeepGWtoLat(), "LDeepGWo", "W");
  }  
  // Cumulative values
  if (oControl->Rep_GWtoChnacc)
    oReport->UpdateOutputNC(oBasin->getAccGWtoChn(), "GWChnA", "W");
  if(oControl->sw_deepGW){
    if (oControl->Rep_DeepGWtoChnacc)
    oReport->UpdateOutputNC(oBasin->getAccDeepGWtoChn(), "DeepGWChnA", "W");
  }
  if (oControl->Rep_SrftoChnacc)
    oReport->UpdateOutputNC(oBasin->getAccSrftoChn(), "SrfChnA", "W");
  if (oControl->Rep_Infiltacc)
    oReport->UpdateOutputNC(oBasin->getAccInfilt(), "InfA", "W");
  if (oControl->Rep_Exfiltacc)
    oReport->UpdateOutputNC(oBasin->getAccExfilt(), "RSrfA", "W");
  if (oControl->Rep_PercolL2acc)
    oReport->UpdateOutputNC(oBasin->getAccPercolL2(), "PrcL2A", "W");
  if (oControl->Rep_ReturnL1acc)
    oReport->UpdateOutputNC(oBasin->getAccL2toL1(), "RL1A", "W");
  if (oControl->Rep_PercolL3acc)
    oReport->UpdateOutputNC(oBasin->getAccPercolL3(), "PrcL3A", "W");
  if (oControl->Rep_Rechargeacc)
    oReport->UpdateOutputNC(oBasin->getAccRecharge(), "RchgA", "W");
  if (oControl->Rep_ReturnL2acc)
    oReport->UpdateOutputNC(oBasin->getAccL3toL2(), "RL2A", "W");
  if (oControl->Rep_EvaporationSacc)
    oReport->UpdateOutputNC(oBasin->getAccEvaporationS(), "EsA", "W");
  if (oControl->Rep_LattoSrfacc)
    oReport->UpdateOutputNC(oBasin->getAccLattoSrf(), "LSrfiA", "W");
  if (oControl->Rep_LattoChnacc)
    oReport->UpdateOutputNC(oBasin->getAccLattoChn(), "LChniA", "W");
  if (oControl->Rep_LattoGWacc)
    oReport->UpdateOutputNC(oBasin->getAccLattoGW(), "LGWiA", "W");
  if(oControl->sw_deepGW){
    if (oControl->Rep_LattoDeepGWacc)
      oReport->UpdateOutputNC(oBasin->getAccLattoDeepGW(), "LDeepGWiA", "W");
  }
  if (oControl->Rep_SrftoLatacc)
    oReport->UpdateOutputNC(oBasin->getAccSrftoLat(), "LSrfoA", "W");
  if (oControl->Rep_GWtoLatacc)
    oReport->UpdateOutputNC(oBasin->getAccGWtoLat(), "LGWoA", "W");
  if(oControl->sw_deepGW){
    if (oControl->Rep_DeepGWtoLatacc)
      oReport->UpdateOutputNC(oBasin->getAccDeepGWtoLat(), "LDeepGWoA", "W");
  }
  if (oControl->Rep_ChntoLatacc)
    oReport->UpdateOutputNC(oBasin->getAccChntoLat(), "LChnoA", "W");

  // two pore stuff
  if(oControl->sw_trck and oControl->sw_TPD){
    if (oControl->Rep_Moist_MW1)
	  oReport->UpdateOutputNC(oBasin->getMoistureMW1(), "SWCMW1", "W");
    if (oControl->Rep_Moist_MW2)
	  oReport->UpdateOutputNC(oBasin->getMoistureMW2(), "SWCMW2", "W");
    if (oControl->Rep_Frac_MW1)
	  oReport->UpdateOutputNC(oBasin->getFracMW1(), "fMW1", "W");
    if (oControl->Rep_Frac_MW2)
	  oReport->UpdateOutputNC(oBasin->getFracMW2(), "fMW2", "W");
    if (oControl->Rep_Frac_MW12)
	  oReport->UpdateOutputNC(oBasin->getFracMW12(), "fMW12", "W");
  }

  // Tracking maps
  if(oControl->sw_trck && oControl->sw_2H){
     //create the nc file for tracer 2H
    if(oControl->current_ts_count == 1){
      oReport->CreatOutputNC(oControl->path_ResultsFolder, "TD");
    }
    if (oControl->Rep_d2Hcanopy_sum)
      oReport->UpdateOutputNC(oTracking->getd2Hcanopy_sum(), "dHcnp", "TD");
    if (oControl->Rep_d2Hprecip)
      oReport->UpdateOutputNC(oAtmosphere->getd2Hprecip(), "dHpcp", "TD");
    if (oControl->Rep_d2Hsnowpack)
      oReport->UpdateOutputNC(oTracking->getd2Hsnowpack(), "dHsnw", "TD");
    if (oControl->Rep_d2Hsurface)
      oReport->UpdateOutputNC(oTracking->getd2Hsurface(), "dHsrf", "TD");
    if (oControl->Rep_d2Hchan)
      oReport->UpdateOutputNC(oTracking->getd2Hchan(), "dHchn", "TD");
    if (oControl->Rep_d2Hsoil1)
      oReport->UpdateOutputNC(oTracking->getd2Hsoil1(), "dHsL1", "TD");
    if (oControl->Rep_d2Hsoil2)
      oReport->UpdateOutputNC(oTracking->getd2Hsoil2(), "dHsL2", "TD");
    if (oControl->Rep_d2HsoilUp)
      oReport->UpdateOutputNC(oTracking->getd2Hsoil_12(), "dHsUp", "TD");
    if (oControl->Rep_d2Hsoil3)
      oReport->UpdateOutputNC(oTracking->getd2Hsoil3(), "dHsL3", "TD");
    if (oControl->Rep_d2HsoilAv)
      oReport->UpdateOutputNC(oTracking->getd2Hsoil_Av(), "dHsAv", "TD");
    if (oControl->Rep_d2Hgroundwater)
      oReport->UpdateOutputNC(oTracking->getd2Hgroundwater(), "dHgw", "TD");
    if(oControl->sw_deepGW){
      if (oControl->Rep_d2HdeepGW)
	oReport->UpdateOutputNC(oTracking->getd2HdeepGW(), "dHdeepgw", "TD");
      if (oControl->Rep_d2HdeepGWQ)
	oReport->UpdateOutputNC(oTracking->getd2HdeepGWQ(), "dHdeepgwQ", "TD");
    }    
    if (oControl->Rep_d2HevapS_sum)
      oReport->UpdateOutputNC(oTracking->getd2HevapS_sum(), "dHeS", "TD");
    if (oControl->Rep_d2HevapI_sum)
      oReport->UpdateOutputNC(oTracking->getd2HevapI_sum(), "dHeI", "TD");
    if (oControl->Rep_d2HevapT_sum)
      oReport->UpdateOutputNC(oTracking->getd2HevapT_sum(), "dHeT", "TD");
    if (oControl->Rep_d2Hleakage)
      oReport->UpdateOutputNC(oTracking->getd2Hleakage(), "dHLeak", "TD");		  
    if(oControl->sw_TPD){
      if (oControl->Rep_d2H_MW1)
      oReport->UpdateOutputNC(oTracking->getd2H_MW1(), "dHMW1", "TD");
      if (oControl->Rep_d2H_MW2)
      oReport->UpdateOutputNC(oTracking->getd2H_MW2(), "dHMW2", "TD");
      if (oControl->Rep_d2H_TB1)
      oReport->UpdateOutputNC(oTracking->getd2H_TB1(), "dHTB1", "TD");
      if (oControl->Rep_d2H_TB2)
      oReport->UpdateOutputNC(oTracking->getd2H_TB2(), "dHTB2", "TD");
    }
  }

  if(oControl->sw_trck && oControl->sw_18O){
     //create the nc file for tracer 18O
    if(oControl->current_ts_count == 1){
      oReport->CreatOutputNC(oControl->path_ResultsFolder, "TO");
    }
    if (oControl->Rep_d18Ocanopy_sum)
     oReport->UpdateOutputNC(oTracking->getd18Ocanopy_sum(), "dOcnp", "TO");
    if (oControl->Rep_d18Oprecip)
      oReport->UpdateOutputNC(oAtmosphere->getd18Oprecip(), "dOpcp", "TO");
    if (oControl->Rep_d18Osnowpack)
      oReport->UpdateOutputNC(oTracking->getd18Osnowpack(), "dOsnw", "TO");
    if (oControl->Rep_d18Osurface)
      oReport->UpdateOutputNC(oTracking->getd18Osurface(), "dOsrf", "TO");
    if (oControl->Rep_d18Ochan)
      oReport->UpdateOutputNC(oTracking->getd18Ochan(), "dOchn", "TO");
    if (oControl->Rep_d18Osoil1)
      oReport->UpdateOutputNC(oTracking->getd18Osoil1(), "dOsL1", "TO");
    if (oControl->Rep_d18Osoil2)
      oReport->UpdateOutputNC(oTracking->getd18Osoil2(), "dOsL2", "TO");
    if (oControl->Rep_d18OsoilUp)
      oReport->UpdateOutputNC(oTracking->getd18Osoil_12(), "dOsUp", "TO");
    if (oControl->Rep_d18Osoil3)
      oReport->UpdateOutputNC(oTracking->getd18Osoil3(), "dOsL3", "TO");
    if (oControl->Rep_d18OsoilAv)
      oReport->UpdateOutputNC(oTracking->getd18Osoil_Av(), "dOsAv", "TO");
    if (oControl->Rep_d18Ogroundwater)
      oReport->UpdateOutputNC(oTracking->getd18Ogroundwater(), "dOgw", "TO");
    if(oControl->sw_deepGW){
      if (oControl->Rep_d18OdeepGW)
	oReport->UpdateOutputNC(oTracking->getd18OdeepGW(), "dOdeepgw", "TO");
      if (oControl->Rep_d18OdeepGWQ)
	oReport->UpdateOutputNC(oTracking->getd18OdeepGWQ(), "dOdeepgwQ", "TO");
    }        
    if (oControl->Rep_d18OevapS_sum)
      oReport->UpdateOutputNC(oTracking->getd18OevapS_sum(), "dOeS", "TO");
    if (oControl->Rep_d18OevapI_sum)
      oReport->UpdateOutputNC(oTracking->getd18OevapI_sum(), "dOeI", "TO");
    if (oControl->Rep_d18OevapT_sum)
      oReport->UpdateOutputNC(oTracking->getd18OevapT_sum(), "dOeT", "TO");
    if(oControl->sw_TPD){
      if (oControl->Rep_d18O_MW1)
	oReport->UpdateOutputNC(oTracking->getd18O_MW1(), "dOMW1", "TO");
      if (oControl->Rep_d18O_MW2)
	oReport->UpdateOutputNC(oTracking->getd18O_MW2(), "dOMW2", "TO");
      if (oControl->Rep_d18O_TB1)
	oReport->UpdateOutputNC(oTracking->getd18O_TB1(), "dOTB1", "TO");
      if (oControl->Rep_d18O_TB2)
	oReport->UpdateOutputNC(oTracking->getd18O_TB2(), "dOTB2", "TO");
    }
  }
  if(oControl->sw_trck && oControl->sw_Age){
	//create the nc file for tracer Age
    if(oControl->current_ts_count == 1){
      oReport->CreatOutputNC(oControl->path_ResultsFolder, "TA");
    }
    if (oControl->Rep_Agecanopy_sum)
     oReport->UpdateOutputNC(oTracking->getAgecanopy_sum(), "Agecnp", "TA");
    if (oControl->Rep_Agesnowpack)
      oReport->UpdateOutputNC(oTracking->getAgesnowpack(), "Agesnw", "TA");
    if (oControl->Rep_Agesurface)
      oReport->UpdateOutputNC(oTracking->getAgesurface(), "Agesrf", "TA");
    if (oControl->Rep_Agechan)
      oReport->UpdateOutputNC(oTracking->getAgechan(), "Agechn", "TA");
    if (oControl->Rep_Agesoil1)
      oReport->UpdateOutputNC(oTracking->getAgesoil1(), "AgesL1", "TA");
    if (oControl->Rep_Agesoil2)
      oReport->UpdateOutputNC(oTracking->getAgesoil2(), "AgesL2", "TA");
    if (oControl->Rep_AgesoilUp)
      oReport->UpdateOutputNC(oTracking->getAgesoil_12(), "AgesUp", "TA");
    if (oControl->Rep_Agesoil3)
      oReport->UpdateOutputNC(oTracking->getAgesoil3(), "AgesL3", "TA");
    if (oControl->Rep_AgesoilAv)
      oReport->UpdateOutputNC(oTracking->getAgesoil_Av(), "AgesAv", "TA");
    if (oControl->Rep_Agegroundwater)
      oReport->UpdateOutputNC(oTracking->getAgegroundwater(), "Agegw", "TA");
    if(oControl->sw_deepGW){
      if (oControl->Rep_AgedeepGW)
	oReport->UpdateOutputNC(oTracking->getAgedeepGW(), "Agedeepgw", "TA");
      if (oControl->Rep_AgedeepGWQ)
	oReport->UpdateOutputNC(oTracking->getAgedeepGWQ(), "AgedeepgwQ", "TA");
    }        
    if (oControl->Rep_AgeevapS_sum)
      oReport->UpdateOutputNC(oTracking->getAgeevapS_sum(), "AgeeS", "TA");
    if (oControl->Rep_AgeevapI_sum)
      oReport->UpdateOutputNC(oTracking->getAgeevapI_sum(), "AgeeI", "TA");
    if (oControl->Rep_AgeevapT_sum)
      oReport->UpdateOutputNC(oTracking->getAgeevapT_sum(), "AgeeT", "TA");
    if (oControl->Rep_AgeGWtoChn)
      oReport->UpdateOutputNC(oTracking->getAgeGWtoChn(), "AgeGWQ", "TA");
    if (oControl->Rep_AgeSrftoChn)
      oReport->UpdateOutputNC(oTracking->getAgeSrftoChn(), "AgeSrfQ", "TA");
    if (oControl->Rep_AgeRecharge)
      oReport->UpdateOutputNC(oTracking->getAgeRecharge(), "AgeR", "TA");

    if(oControl->sw_TPD){
      if (oControl->Rep_Age_MW1)
	oReport->UpdateOutputNC(oTracking->getAge_MW1(), "AgeMW1", "TA");
      if (oControl->Rep_Age_MW2)
	oReport->UpdateOutputNC(oTracking->getAge_MW2(), "AgeMW2", "TA");
      if (oControl->Rep_Age_MWUp)
	oReport->UpdateOutputNC(oTracking->getAge_MW12(), "AgeMW12", "TA");
      if (oControl->Rep_Age_TB1)
	oReport->UpdateOutputNC(oTracking->getAge_TB1(), "AgeTB1", "TA");
      if (oControl->Rep_Age_TB2)
	oReport->UpdateOutputNC(oTracking->getAge_TB2(), "AgeTB2", "TA");
      if (oControl->Rep_Age_TBUp)
	oReport->UpdateOutputNC(oTracking->getAge_TB12(), "AgeTB12", "TA");
    }
  }
	
  //create the nc file for vegetation
  if(oControl->current_ts_count == 1){
    oReport->CreatOutputNC(oControl->path_ResultsFolder, "VG");
  }
  for (int i = 0; i < oControl->NumSpecs; i++) {

    stringstream name;

    if (oControl->Rep_Veget_frac) {
      name << "p" << i << "_";
      oReport->UpdateOutputNC(oBasin->getVegetFrac(i), name.str() , "VG");
      name.str("");
    }

    if (oControl->Rep_Stem_Density) {
      name << "ntr" << i << "_";
      oReport->UpdateOutputNC(oBasin->getStemDensity(i), name.str() , "VG");
      name.str("");
    }

    if (oControl->Rep_RootFrac1Species) {
      name << "L1Ro" << i << "_";
      oReport->UpdateOutputNC(oBasin->getRootFrac1(i), name.str() , "VG");
      name.str("");
    }
    if (oControl->Rep_RootFrac2Species) {
      name << "L2Ro" << i << "_";
      oReport->UpdateOutputNC(oBasin->getRootFrac2(i), name.str() , "VG");
      name.str("");
    }

    if (oControl->Rep_Leaf_Area_Index) {
      name << "lai" << i << "_";
      oReport->UpdateOutputNC(oBasin->getLAI(i), name.str() , "VG");
      name.str("");
    }

    if (oControl->Rep_Stand_Age) {
      name << "age" << i << "_";
      oReport->UpdateOutputNC(oBasin->getStandAge(i), name.str() , "VG");
      name.str("");
    }

    if (oControl->Rep_Canopy_Conductance) {
      name << "gc" << i << "_";
      oReport->UpdateOutputNC(oBasin->getCanopyCond(i), name.str() , "VG");
      name.str("");
    }

    if (oControl->Rep_GPP) {
      name << "gpp" << i << "_";
      oReport->UpdateOutputNC(oBasin->getGPP(i), name.str() , "VG");
      name.str("");
    }

    if (oControl->Rep_NPP) {
      name << "npp" << i << "_";
      oReport->UpdateOutputNC(oBasin->getNPP(i), name.str() , "VG");
      name.str("");
    }

    if (oControl->Rep_Basal_Area) {
      name << "bas" << i << "_";
      oReport->UpdateOutputNC(oBasin->getBasalArea(i), name.str() , "VG");
      name.str("");
    }

    if (oControl->Rep_Tree_Height) {
      name << "hgt" << i << "_";
      oReport->UpdateOutputNC(oBasin->getTreeHeight(i), name.str() , "VG");
      name.str("");
    }

    if (oControl->Rep_Root_Mass) {
      name << "root" << i << "_";
      oReport->UpdateOutputNC(oBasin->getRootMass(i), name.str() , "VG");
      name.str("");
    }

    if (oControl->Rep_Canopy_Temp) {
      name << "Tc" << i << "_";
      oReport->UpdateOutputNC(oBasin->getCanopyTemp(i), name.str() , "VG");
      name.str("");
    }

    if (oControl->Rep_Canopy_NetR) {
      name << "NRc" << i << "_";
      oReport->UpdateOutputNC(oBasin->getCanopyNetRad(i), name.str() , "VG");
      name.str("");
    }

    if (oControl->Rep_Canopy_LE_E) {
      name << "LEEi" << i << "_";
      oReport->UpdateOutputNC(oBasin->getCanopyLatHeatE(i), name.str() , "VG");
      name.str("");
    }
    if (oControl->Rep_Canopy_LE_T) {
      name << "LETr" << i << "_";
      oReport->UpdateOutputNC(oBasin->getCanopyLatHeatT(i), name.str() , "VG");
      name.str("");
    }

    if (oControl->Rep_Canopy_Sens_Heat) {
      name << "Hc" << i << "_";
      oReport->UpdateOutputNC(oBasin->getCanopySensHeat(i), name.str() , "VG");
      name.str("");
    }

    if (oControl->Rep_Canopy_Water_Stor) {
      name << "Cs" << i << "_";
      oReport->UpdateOutputNC(oBasin->getCanopyWaterStor(i), name.str() , "VG");
      name.str("");
    }

    if (oControl->Rep_ETspecies) {
      name << "ETc" << i << "_";
      oReport->UpdateOutputNC(oBasin->getETspecies(i), name.str() , "VG");
      name.str("");
    }
    if (oControl->Rep_Transpiration) {
      name << "Tr" << i << "_";
      oReport->UpdateOutputNC(oBasin->getTranspiration(i), name.str() , "VG");
      name.str("");
    }
    if (oControl->Rep_Einterception) {
      name << "Ei" << i << "_";
      oReport->UpdateOutputNC(oBasin->getEinterception(i), name.str() , "VG");
      name.str("");
    }
    if (oControl->Rep_Esoil) {
      name << "Es" << i << "_";
      oReport->UpdateOutputNC(oBasin->getEsoil(i), name.str() , "VG");
      name.str("");
    }

    // Tracking
    if(oControl->sw_trck){
      if (oControl->sw_2H && oControl->Rep_d2Hcanopy) {
	name << "dHCs" << i << "_";
	oReport->UpdateOutputNC(oBasin->getd2Hcanopy(i), name.str() , "VG");
	name.str("");
      }
      
      if (oControl->sw_18O && oControl->Rep_d18Ocanopy) {
	name << "dOCs" << i << "_";
	oReport->UpdateOutputNC(oBasin->getd18Ocanopy(i), name.str() , "VG");
	name.str("");
      }
      
      if (oControl->sw_Age && oControl->Rep_Agecanopy) {
	name << "AgeCs" << i << "_";
	oReport->UpdateOutputNC(oBasin->getAgecanopy(i), name.str() , "VG");
	name.str("");
      }
      if (oControl->sw_2H && oControl->Rep_d2HevapI) {
	name << "dHeI" << i << "_";
	oReport->UpdateOutputNC(oBasin->getd2HevapI(i), name.str() , "VG");
	name.str("");
      }

      if (oControl->sw_18O && oControl->Rep_d18OevapI) {
	name << "dOeI" << i << "_";
	oReport->UpdateOutputNC(oBasin->getd18OevapI(i), name.str() , "VG");
	name.str("");
      }

      if (oControl->sw_Age && oControl->Rep_AgeevapI) {
	name << "AgeeI" << i << "_";
	oReport->UpdateOutputNC(oBasin->getAgeevapI(i), name.str() , "VG");
	name.str("");
      }
      if (oControl->sw_2H && oControl->Rep_d2HevapT) {
	name << "dHeT" << i << "_";
	oReport->UpdateOutputNC(oBasin->getd2HevapT(i), name.str() , "VG");
	name.str("");
      }

      if (oControl->sw_18O && oControl->Rep_d18OevapT) {
	name << "dOeT" << i << "_";
	oReport->UpdateOutputNC(oBasin->getd18OevapT(i), name.str() , "VG");
	name.str("");
      }

      if (oControl->sw_Age && oControl->Rep_AgeevapT) {
	name << "AgeeT" << i << "_";
	oReport->UpdateOutputNC(oBasin->getAgeevapT(i), name.str() , "VG");
	name.str("");
      }
      if (oControl->sw_2H && oControl->Rep_d2HevapS) {
	name << "dHeS" << i << "_";
	oReport->UpdateOutputNC(oBasin->getd2HevapS(i), name.str() , "VG");
	name.str("");
      }

      if (oControl->sw_18O && oControl->Rep_d18OevapS) {
	name << "dOeS" << i << "_";
	oReport->UpdateOutputNC(oBasin->getd18OevapS(i), name.str() , "VG");
	name.str("");
      }

      if (oControl->sw_Age && oControl->Rep_AgeevapS) {
	name << "AgeeS" << i << "_";
	oReport->UpdateOutputNC(oBasin->getAgeevapS(i), name.str() , "VG");
	name.str("");
      }
    }

  }

  return EXIT_SUCCESS;
}




