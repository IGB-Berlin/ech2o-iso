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
 *    Marco Maneta
 *******************************************************************************/
/*
 * InitConf.h
 *
 *  Created on: Oct 13, 2009
 *      Author: Marco Maneta
 */

#ifndef INITCONF_H_
#define INITCONF_H_

#include "ParsingFunctions.h"

struct Control{
  
  /*Folder paths*/
  string path_BasinFolder; //folder where basin property maps are located
  string path_ClimMapsFolder; //folder where weather maps series are located
  string path_ResultsFolder; //folder where results will be placed
  
  /*Time control variables*/
  
  float starttime; //simulation start time (seconds)
  float endtime; //simulation end time (seconds)
  float dt; //simulation time step (seconds)
  float BC_dt; // time step for spatial climatic inputs (seconds)
  float report_times; //times at which report outputs time series
  float reportMap_times; //times at which report outputs maps
  float reportMap_start; // absolute time from which report outputs maps
  
  float current_t_step; //current time step (seconds)
  unsigned int current_ts_count; //current count of time step
  
  
  /*Control switches*/
  string MapType; //indicates if the maps to be read are ASCII (grass) or PCRASTER (csf)
  string ForestStateVarsInputType; //indicates if the forest state variables are input as tables (tables) or maps (maps)
  
  /*Option switches*/
  bool sw_reinfilt; //switch to turn on and off the reinfiltration option
  bool sw_channel; //switch to turn on and off the channel option
  bool sw_expKsat; //switch to turn on and off the exponential profile for hydraulic conductivity
  bool sw_expPoros; //switch to turn on and off the exponential profile for porosity
  bool sw_chan_evap; //switch to turn on and off the channel evaporation processes
  bool sw_ddSoilPar; //swith to turn on and off the import of different soil parameters for each layer

  /*multiple option switches*/
  int toggle_infilt; //switch to 
  //int toggle_soil_water_profile; //toggle between different soil moisture profile calculation
  int toggle_veg_dyn; //switch to turn on and off the dynamic vegetation module (allocation and lai calculation)
  int toggle_ra; //toggle between aerodynamic resistance options
  int toggle_rs; //toggle between different soil resistance option

  // LAI time series binaries (if toggle_veg_dyn ==2)
  // file names (per species) = name below + "_"+ species number (starting at 0) + ".bin"
  string fn_LAI_timeseries;
  
  /*Base maps filenames*/
  string fn_dem; //local base dem filename that forces grid geometry
  string fn_ldd; //local drain direction map filename
  string fn_chwidth; //channel width (m)
  string fn_chgwparam; //channel water transfer parameter
  string fn_chmanningn; //channel roughness parameter
  
  /*Soil properties and parameters*/
  string fn_Ksat0; //top-of-column soil hydraulic conductivity ms-1
  string fn_kKsat; //soil hydraulic conductivity profile coeff m
  string fn_kvkh; //vertical to horizontal ksat anisotropy ratio
  string fn_randrough; //terrain base random roughness to calcualte aerodynamic resistance (m)
  string fn_slope; //surface slope m m-1
  string fn_poros0; //top-of-column porosity (m3.m-3)
  string fn_kporos; // porosity profile coeff (m)
  string fn_psi_ae; //soil air entry pressure m
  string fn_BClambda; //brooks and corey lambda param
  string fn_theta_r; //residual soil moisture
  string fn_soildepth; //soil depth in m
  string fn_depth_layer1; //depth of layer 1 in m
  string fn_depth_layer2;  //depth of layer 2 in m. Layer 3 evaluated from soil depth
  //string fn_root_fraction_lay1; //fraction of roots in soil layer 1
  //string fn_root_fraction_lay2; // fraction of roots in soil layer 2. Soil layer 3 implied
  string fn_Kroot; // coefficient for exponential root profile, in m-1
  string fn_bedrock_leak; //bedrock leakance in s-1
  string fn_paramWc; //empirical parameter in water efficiency function for GPP calculation (see Landsber and Waring, 1997 or TRIPLEX paper
  string fn_paramWp;//empirical parameter in water efficiency function for GPP calculation (see Landsber and Waring, 1997 or TRIPLEX paper
  string fn_snowCf; //empirical parameter that controls the snowmelt rates m s-1 C-1

  string fn_Ksat2, fn_Ksat3; //average soil hydraulic conductivity ms-1
  string fn_kvkh2, fn_kvkh3; //vertical to horizontal hydraulic conductivity
  string fn_poros2, fn_poros3; //soil layer porosity (m3.m-3)
  string fn_psi_ae2, fn_psi_ae3; //soil layer air entry pressure m
  string fn_BClambda2, fn_BClambda3; //soil layer brooks and corey lambda param 
  
  /*Basin state variables*/
  string fn_swe;
  string fn_albedo;
  string fn_emiss;
  string fn_soilheatcap;
  string fn_soilthermcond;
  string fn_dampdepth;
  string fn_tempdamp;
  string fn_streamflow;
  string fn_soilmoist;
  string fn_soilmoist2;
  string fn_soilmoist3;
  string fn_soiltemp;

  string fn_temp_w;    // water temperature (C)
  string fn_chanrough; //channel roughness to calculate aerodyamic resistance (m)
  
  /*Climate zones and climate input files*/
  float snow_rain_temp; //threshold temp for snow rain transition, degC
  string fn_climzones; //base climatic zones map with the grid geometry
  string fn_isohyet; //map with rainfall multipliers
  string fn_Ldown; //Incoming longwave radiation Wm-2
  string fn_Sdown; //Incoming shortwave radiation Wm-2
  string fn_temp; //average air temperature C
  string fn_maxTemp; //max air temp C
  string fn_minTemp; //min air temp C
  string fn_precip; //
  string fn_rel_humid; //relative humidity
  string fn_wind_speed; //wind speed ms-1
  
  /*Forest patches and forest input files*/
  int NumSpecs; //number of tree species in the simulation
  string fn_patches;
  string fn_paramtable;
  string fn_proptable;
  string fn_StemDenstable;
  string fn_LAItable;
  string fn_AGEtable;
  string fn_BasalAreatable;
  string fn_Heighttable;
  string fn_RootMasstable;
  
  /*report flags*/
  bool Rep_Long_Rad_Down;
  bool Rep_Short_Rad_Down;
  bool Rep_Precip;
  bool Rep_Rel_Humidity;
  bool Rep_Wind_Speed;
  bool Rep_AvgAir_Temperature;
  bool Rep_MinAir_Temperature;
  bool Rep_MaxAir_Temperature;
  bool Rep_SWE;
  bool Rep_Infilt_Cap;
  bool Rep_Streamflow;
  bool Rep_Saturation_Area;
  bool Rep_Ponding;
  bool Rep_Soil_Water_Content_Average;
  bool Rep_Soil_Water_Content_12;
  bool Rep_Soil_Water_Content_L1;
  bool Rep_Soil_Water_Content_L2;
  bool Rep_Soil_Water_Content_L3;
  bool Rep_WaterTableDepth;
  bool Rep_RootZone_in_L1;
  bool Rep_RootZone_in_L2;
  bool Rep_RootZone_in_L3;
  bool Rep_Field_Capacity_L1;
  bool Rep_Field_Capacity_L2;
  bool Rep_Field_Capacity_L3;
  bool Rep_Soil_Sat_Deficit;
  bool Rep_GWater;
  bool Rep_Soil_Net_Rad;
  bool Rep_Soil_LE;
  bool Rep_Sens_Heat;
  bool Rep_Grnd_Heat;
  bool Rep_Snow_Heat;
  bool Rep_Soil_Temperature;
  bool Rep_Skin_Temperature;
  
  bool Rep_GWtoChn;
  bool Rep_SrftoChn;
  bool Rep_GWtoChnacc;
  bool Rep_SrftoChnacc;

  bool Rep_Infilt;
  bool Rep_Exfilt;
  bool Rep_PercolL2;
  bool Rep_ReturnL1;
  bool Rep_PercolL3;
  bool Rep_Recharge;
  bool Rep_ReturnL2;
  bool Rep_Leak;  

  bool Rep_LattoSrf;
  bool Rep_LattoChn;
  bool Rep_LattoGW;
  bool Rep_ChntoLat;
  bool Rep_SrftoLat;
  bool Rep_GWtoLat;

  bool Rep_Infiltacc;
  bool Rep_Exfiltacc;
  bool Rep_PercolL2acc;
  bool Rep_ReturnL1acc;
  bool Rep_PercolL3acc;
  bool Rep_Rechargeacc;
  bool Rep_ReturnL2acc;
  bool Rep_EvaporationSacc;
  bool Rep_TranspiL1acc;
  bool Rep_TranspiL2acc;
  bool Rep_TranspiL3acc;
  bool Rep_LattoSrfacc;
  bool Rep_LattoChnacc;
  bool Rep_LattoGWacc;
  bool Rep_ChntoLatacc;
  bool Rep_SrftoLatacc;
  bool Rep_GWtoLatacc;

  bool Rep_Net_Rad_sum;
  bool Rep_LE_sum;
  bool Rep_Total_ET;
  bool Rep_Transpiration_sum;
  bool Rep_Einterception_sum;
  bool Rep_Esoil_sum;
  bool Rep_Canopy_Water_Stor_sum;
  
  bool Rep_Veget_frac;
  bool Rep_Stem_Density;
  bool Rep_RootFrac1Species;
  bool Rep_RootFrac2Species;
  bool Rep_Leaf_Area_Index;
  bool Rep_Stand_Age;
  bool Rep_Canopy_Conductance;
  bool Rep_GPP;
  bool Rep_NPP;
  bool Rep_Basal_Area;
  bool Rep_Tree_Height;
  bool Rep_Root_Mass;
  bool Rep_Canopy_Temp;
  bool Rep_Canopy_NetR;
  bool Rep_Canopy_LE_E;
  bool Rep_Canopy_LE_T;
  bool Rep_Canopy_Sens_Heat;
  bool Rep_Canopy_Water_Stor;
  bool Rep_ETspecies;
  bool Rep_Transpiration;
  bool Rep_Einterception;
  bool Rep_Esoil;

  bool Rep_Water_Temperature;
  bool Rep_ChanEvap;

  bool Rep_Pond_F_Chn;
  
  /*time series reporting input files*/
  string fn_rep_mask;

  bool RepTs_OutletDischarge; //only reported at the outlets
  bool RepTs_OutletGW; //only reported at the outlets
  bool RepTs_Long_Rad_Down;
  bool RepTs_Short_Rad_Down;
  bool RepTs_Precip;
  bool RepTs_Rel_Humidity;
  bool RepTs_Wind_Speed;
  bool RepTs_AvgAir_Temperature;
  bool RepTs_MinAir_Temperature;
  bool RepTs_MaxAir_Temperature;
  bool RepTs_SWE;
  bool RepTs_Infilt_Cap;
  bool RepTs_Streamflow;
  bool RepTs_Ponding;
  bool RepTs_Soil_Water_Content_Average;
  bool RepTs_Soil_Water_Content_12;
  bool RepTs_Soil_Water_Content_L1;
  bool RepTs_Soil_Water_Content_L2;
  bool RepTs_Soil_Water_Content_L3;
  bool RepTs_WaterTableDepth;
  bool RepTs_Field_Capacity_L1;
  bool RepTs_Field_Capacity_L2;
  bool RepTs_Field_Capacity_L3;
  bool RepTs_Soil_Sat_Deficit;
  bool RepTs_GroundWater;
  bool RepTs_Soil_Net_Rad;
  bool RepTs_Soil_LE;
  bool RepTs_Sens_Heat;
  bool RepTs_Grnd_Heat;
  bool RepTs_Snow_Heat;
  bool RepTs_Soil_Temperature;
  bool RepTs_Skin_Temperature;

  bool RepTs_GWtoChn;
  bool RepTs_SrftoChn;

  bool RepTs_Infilt;
  bool RepTs_Exfilt;
  bool RepTs_PercolL2;
  bool RepTs_ReturnL1;
  bool RepTs_PercolL3;
  bool RepTs_Recharge;
  bool RepTs_ReturnL2;
  bool RepTs_Leak;
  
  bool RepTs_LattoSrf;
  bool RepTs_LattoChn;
  bool RepTs_LattoGW;
  bool RepTs_ChntoLat;
  bool RepTs_SrftoLat;
  bool RepTs_GWtoLat;

  bool RepTs_Net_Rad_sum;
  bool RepTs_LE_sum;
  bool RepTs_Total_ET;
  bool RepTs_Transpiration_sum;
  bool RepTs_Einterception_sum;
  bool RepTs_Esoil_sum;
  bool RepTs_Canopy_Water_Stor_sum;
  
  bool RepTs_Veget_frac;
  bool RepTs_Stem_Density;
  bool RepTs_RootFrac1Species;
  bool RepTs_RootFrac2Species;
  bool RepTs_Leaf_Area_Index;
  bool RepTs_Canopy_Conductance;
  bool RepTs_GPP;
  bool RepTs_NPP;
  bool RepTs_Basal_Area;
  bool RepTs_Tree_Height;
  bool RepTs_Root_Mass;
  bool RepTs_Canopy_Temp;
  bool RepTs_Canopy_NetR;
  bool RepTs_Canopy_LE_E;
  bool RepTs_Canopy_LE_T;
  bool RepTs_Canopy_Sens_Heat;
  bool RepTs_Canopy_Water_Stor;
  bool RepTs_ETspecies;
  bool RepTs_Transpiration;
  bool RepTs_Einterception;
  bool RepTs_Esoil;

  bool RepTs_Water_Temperature;
  bool RepTs_ChanEvap;

  bool RepTs_Pond_F_Chn;
  
  // Tracking  -------------------------------------------------
  // Tracking inputs
  string fn_tracking;
  bool sw_trck; //switch to turn on and off the tracking option
  bool sw_2H; //switch to turn on and off the 2H tracking option (if sw_trck = 1)
  bool sw_18O; //switch to turn on and off the 18O tracking option (if sw_trck = 1)
  bool sw_Age; //switch to turn on and off the age tracking option (if sw_trck = 1)
  bool sw_frac; //switch to turn on and off fractionation in soil evap (if sw_trck = 1)
  bool sw_TPD; //switch to turn on the two pore domain option (if sw_trck = 1)

  // Toggle switch for fractionation
  bool sw_chan_frac;  //switch to turn on and off the channel fractionation (if sw_trck = 1)
  int toggle_hs; // toggle to choose which surface relative humidity for fractionation: 0->hs=1, 1->Lee&Pielke 1992, 2->Soderberg 2012
  int toggle_n; // toggle to choose how the turbulent factor is calculated for kinetic fractionation: 0->n=1, 1->follows Mathieu and Bariac 1996
  int toggle_ek; // choose how kinetic fractio. factor is calculated: 0=Merlivat, 1=Vogt, 2=Merlivat & Jouzel
  int toggle_mix; // choose how input+output mixing is done: 
  // 0->immediate mixing without volume change with output, (Vag=V(t), Cag=C(t+1))
  // 1->using average "useful" volume (still without sotrage limit if Fin or Fout >> Vt, Vt+1): 
  //    Cag=(C(t)+C(t+1))/2, Vag=(V(t)+Fin+max(0,V(t)-Fout))/2
 
  string fn_psi_MW; // Transition pressure (in meters of head) from tightly-bound
  // to tightly-bound to mobile water (if sw_TPD =1)

  string fn_IncompMix; // incomplete mixing alpha parameter

  /* input maps for initial values*/
  string fn_d2Hprecip; // deuterium signature in precipitations (2H, per mil)
  //string fn_d2Hcanopy;
  string fn_d2Hsnowpack;
  string fn_d2Hsurface;
  string fn_d2Hsoil1;
  string fn_d2Hsoil2;
  string fn_d2Hsoil3;
  string fn_d2Hgroundwater;
  
  string fn_d18Oprecip; // O eighteen signature in precipitations (18O, per mil)
  //string fn_d2Hcanopy;
  string fn_d18Osnowpack;
  string fn_d18Osurface;
  string fn_d18Osoil1;
  string fn_d18Osoil2;
  string fn_d18Osoil3;
  string fn_d18Ogroundwater;
  
  //string fn_Agecanopy;
  string fn_Agesnowpack;
  string fn_Agesurface;
  string fn_Agesoil1;
  string fn_Agesoil2;
  string fn_Agesoil3;
  string fn_Agegroundwater;
  
  /* maps report */
  bool Rep_Moist_MW1;
  bool Rep_Moist_MW2;
  bool Rep_Frac_MW1;
  bool Rep_Frac_MW2;
  bool Rep_Frac_MW12;

  bool Rep_d2Hprecip;
  bool Rep_d2Hcanopy;
  bool Rep_d2Hcanopy_sum;
  bool Rep_d2Hsnowpack;
  bool Rep_d2Hsurface;
  bool Rep_d2Hchan;
  bool Rep_d2Hsoil1;
  bool Rep_d2Hsoil2;
  bool Rep_d2HsoilUp;
  bool Rep_d2Hsoil3;
  bool Rep_d2HsoilAv;
  bool Rep_d2Hgroundwater;
  bool Rep_d2Hleakage;
  bool Rep_d2HevapS;
  bool Rep_d2HevapS_sum;
  bool Rep_d2HevapI;
  bool Rep_d2HevapI_sum;
  bool Rep_d2HevapT;
  bool Rep_d2HevapT_sum;
  bool Rep_d2H_MW1;
  bool Rep_d2H_MW2;
  bool Rep_d2H_TB1;
  bool Rep_d2H_TB2;
  
  bool Rep_d18Oprecip;
  bool Rep_d18Ocanopy;
  bool Rep_d18Ocanopy_sum;
  bool Rep_d18Osnowpack;
  bool Rep_d18Osurface;
  bool Rep_d18Ochan;
  bool Rep_d18Osoil1;
  bool Rep_d18Osoil2;
  bool Rep_d18OsoilUp;
  bool Rep_d18Osoil3;
  bool Rep_d18OsoilAv;
  bool Rep_d18Ogroundwater;
  bool Rep_d18Oleakage;
  bool Rep_d18OevapS;
  bool Rep_d18OevapS_sum;
  bool Rep_d18OevapI;
  bool Rep_d18OevapI_sum;
  bool Rep_d18OevapT;
  bool Rep_d18OevapT_sum;
  bool Rep_d18O_MW1;
  bool Rep_d18O_MW2;
  bool Rep_d18O_TB1;
  bool Rep_d18O_TB2;
  
  bool Rep_Agecanopy;
  bool Rep_Agecanopy_sum;
  bool Rep_Agesnowpack;
  bool Rep_Agesurface;
  bool Rep_Agechan;
  bool Rep_Agesoil1;
  bool Rep_Agesoil2;
  bool Rep_AgesoilUp;
  bool Rep_Agesoil3;
  bool Rep_AgesoilAv;
  bool Rep_Agegroundwater;
  bool Rep_Ageleakage;
  bool Rep_AgeevapS;
  bool Rep_AgeevapS_sum;
  bool Rep_AgeevapI;
  bool Rep_AgeevapI_sum;
  bool Rep_AgeevapT;
  bool Rep_AgeevapT_sum;
  bool Rep_AgeGWtoChn;
  bool Rep_AgeSrftoChn;
  bool Rep_AgeRecharge;
  bool Rep_Age_MW1;
  bool Rep_Age_MW2;
  bool Rep_Age_MWUp;
  bool Rep_Age_TB1;
  bool Rep_Age_TB2;
  bool Rep_Age_TBUp;
  
  // Time series
  bool RepTs_Moist_MW1;
  bool RepTs_Moist_MW2;
  bool RepTs_Frac_MW1;
  bool RepTs_Frac_MW2;
  bool RepTs_Frac_MW12;

  bool RepTs_d2Hprecip;
  bool RepTs_d2Hcanopy;
  bool RepTs_d2Hcanopy_sum;
  bool RepTs_d2Hsnowpack;
  bool RepTs_d2Hsurface;
  bool RepTs_d2Hchan;
  bool RepTs_d2Hsoil1;
  bool RepTs_d2Hsoil2;
  bool RepTs_d2HsoilUp;
  bool RepTs_d2Hsoil3;
  bool RepTs_d2HsoilAv;
  bool RepTs_d2Hgroundwater;
  bool RepTs_d2Hleakage;
  bool RepTs_d2HevapS;
  bool RepTs_d2HevapS_sum;
  bool RepTs_d2HevapI;
  bool RepTs_d2HevapI_sum;
  bool RepTs_d2HevapT;
  bool RepTs_d2HevapT_sum;
  bool RepTs_d2H_MW1;
  bool RepTs_d2H_MW2;
  bool RepTs_d2H_TB1;
  bool RepTs_d2H_TB2;
  
  bool RepTs_d18Oprecip;
  bool RepTs_d18Ocanopy;
  bool RepTs_d18Ocanopy_sum;
  bool RepTs_d18Osnowpack;
  bool RepTs_d18Osurface;
  bool RepTs_d18Ochan;
  bool RepTs_d18Osoil1;
  bool RepTs_d18Osoil2;
  bool RepTs_d18OsoilUp;
  bool RepTs_d18Osoil3;
  bool RepTs_d18OsoilAv;
  bool RepTs_d18Ogroundwater;
  bool RepTs_d18Oleakage;
  bool RepTs_d18OevapS;
  bool RepTs_d18OevapS_sum;
  bool RepTs_d18OevapI;
  bool RepTs_d18OevapI_sum;
  bool RepTs_d18OevapT;
  bool RepTs_d18OevapT_sum;
  bool RepTs_d18O_MW1;
  bool RepTs_d18O_MW2;
  bool RepTs_d18O_TB1;
  bool RepTs_d18O_TB2;
    
  bool RepTs_Agecanopy;
  bool RepTs_Agecanopy_sum;
  bool RepTs_Agesnowpack;
  bool RepTs_Agesurface;
  bool RepTs_Agechan;
  bool RepTs_Agesoil1;
  bool RepTs_Agesoil2;
  bool RepTs_AgesoilUp;
  bool RepTs_Agesoil3;
  bool RepTs_AgesoilAv;
  bool RepTs_Agegroundwater;
  bool RepTs_Ageleakage;
  bool RepTs_AgeevapS;
  bool RepTs_AgeevapS_sum;
  bool RepTs_AgeevapI;
  bool RepTs_AgeevapI_sum;
  bool RepTs_AgeevapT;
  bool RepTs_AgeevapT_sum;
  bool RepTs_AgeGWtoChn;
  bool RepTs_AgeSrftoChn;
  bool RepTs_AgeRecharge;
  bool RepTs_Age_MW1;
  bool RepTs_Age_MW2;
  bool RepTs_Age_MWUp;
  bool RepTs_Age_TB1;
  bool RepTs_Age_TB2;
  bool RepTs_Age_TBUp;
    // --------------------------
  
  
  Control(){ 	current_ts_count = 1;
    toggle_ra = 0;
  };
  
  int ReadConfigFile(string confilename = "config.ini");
  
  int ReadConfigTrck(string confilename = "configTrck.ini");
  
  void AdvanceTimeStep(){
    current_ts_count++;
    current_t_step = current_ts_count * dt;
  }
  
  float GetTimeStep() {
    return current_t_step;
  }
  
  unsigned int GetTimeStepCount() {
    return current_ts_count;
  }
};

#endif /* INITCONF_H_ */
