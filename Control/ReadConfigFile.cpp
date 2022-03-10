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
 * initConf.cpp
 *
 *  Created on: Oct 13, 2009
 *      Author: Marco Maneta
 *
 *      Reads the configuration file (default name 'config.ini')
 *      and places the file names and control data in the Control object
 */
#include <stdexcept>   // for exception, runtime_error, out_of_range
#include <iostream>
#include <sys/stat.h>
#include <stdlib.h>
#include "ConfigFile.h"
#include "InitConf.h"
using namespace std;

int Control::ReadConfigFile(string confilename /*= "config.ini"*/)
{

  struct stat info; //needed to check the status of the folders

  try{
    ConfigFile Config(confilename);

    Config.readInto(path_BasinFolder, "Maps_Folder");
    //check for slash \ at the end of the path string and if there is none appends it
    if(path_BasinFolder.at(path_BasinFolder.length()-1) != '/')
      path_BasinFolder.append("/");
    if(stat(path_BasinFolder.c_str(), &info)!=0)
      throw std::runtime_error(string("Folder not found: ") + path_BasinFolder);

    Config.readInto(path_ClimMapsFolder, "Clim_Maps_Folder");
    if(path_ClimMapsFolder.at(path_ClimMapsFolder.length()-1) != '/')
      path_ClimMapsFolder.append("/");
    if(stat(path_ClimMapsFolder.c_str(), &info)!=0)
      throw std::runtime_error(string("Folder not found: ") + path_ClimMapsFolder);

    closure_tolerance = Config.read<float>("closure_tolerance");
    
    Config.readInto(path_ResultsFolder, "Output_Folder");
    if(path_ResultsFolder.at(path_ResultsFolder.length()-1) != '/')
      path_ResultsFolder.append("/");
    if(stat(path_ResultsFolder.c_str(), &info)!=0)
      throw std::runtime_error(string("Folder not found: ") + path_ResultsFolder);

    sw_reinfilt = Config.read<bool>("Reinfiltration");
    sw_channel = Config.read<bool>("Channel");
    sw_anthr_heat = Config.read<bool>("Anthropogenic_heat");
    sw_intercept = Config.read<bool>("Interception_proportion");
    sw_deepGW = Config.read<bool>("DeepGW_Storage"); 

    // Vegetation dynamics 
    toggle_veg_dyn = Config.read<int>("Vegetation_dynamics");
    toggle_sm = Config.read<int>("Stomatal_model_opt");    
    toggle_plant_hydr = Config.read<int>("Plant_hydraulics");
    if(toggle_veg_dyn == 2){
      Config.readInto(fn_LAI_timeseries, "TimeSeries_LAI");
      Config.readInto(fn_hgt_timeseries, "TimeSeries_Height");
    }
    
    // Soil dynamics
    toggle_hydrologic_engine  = Config.read<int>("Hydrologic_engine_opt");
    toggle_chan_evap = Config.read<int>("Channel_Evaporation");
    toggle_ra = Config.read<int>("Aerodyn_resist_opt");
    toggle_rs = Config.read<int>("Soil_resistance_opt");

    // Boundary conditions
    sw_BC = Config.read<bool>("Boundary_Condition");
    if(sw_BC){
      Config.readInto(fn_BCsurface, "TimeSeries_BC_Surface");
      Config.readInto(fn_BCgroundwater, "TimeSeries_BC_Groundwater");
      if(sw_deepGW)
	Config.readInto(fn_BCdeepgwtr,"TimeSeries_BC_DeepGW");
    }

    // Soil properties with depth
    toggle_soil_prop = Config.read<int>("Soil_Properties");
    if(toggle_soil_prop == 1){
      sw_expKsat = Config.read<bool>("Hydraulic_Conductivity_profile");
      sw_expPoros = Config.read<bool>("Porosity_profile");
    }

    // Tracking
    sw_trck = Config.read<bool>("Tracking");
    // If there is water tracking, another config file will be read
    if(sw_trck){
      Config.readInto(fn_tracking, "TrackingConfig");
      ReadConfigTrck(fn_tracking);
      cout <<"+ Water tracking is activated: "<< fn_tracking << " read ok... " << "\n";
    }

    starttime = Config.read<float>("Simul_start");
    endtime = Config.read<float>("Simul_end");
    dt = Config.read<float>("Simul_tstep");
    BC_dt = Config.read<float>("Clim_input_tstep");

    report_times = Config.read<float>("Report_interval");
    if(report_times < dt){
      cout << "WARNING: Report time step less than simulation time step. Changing report time step to simulation time step" << endl;
      report_times = dt;
    }

    reportMap_times = Config.read<float>("ReportMap_interval");
    if(reportMap_times < dt){
      cout << "WARNING: ReportMap time step less than simulation time step. Changing reportMap time step to simulation time step" << endl;
      reportMap_times = dt;
    }

    reportMap_start = Config.read<float>("ReportMap_starttime");

    sw_netcdf = Config.read<bool>("NetCDF_output_format");

    Config.readInto(MapType, "MapTypes" );
    Config.readInto(ForestStateVarsInputType, "Species_State_Variable_Input_Method");

    NumSpecs = Config.read<int>("Number_of_Species");

    Config.readInto(fn_dem, "DEM");
    Config.readInto(fn_ttarea, "Fraction_Contributing_Area"); 
    Config.readInto(fn_climzones, "ClimateZones");
    Config.readInto(fn_patches, "ForestPatches");

    Config.readInto(fn_ksat_skin, "Soil_Skin_Infilt_Capacity");
    Config.readInto(fn_Ksat0, "Top-of-profile_Horiz_Hydraulic_Conductivity");
    Config.readInto(fn_kvkh, "Vert_Horz_Anis_ratio");
    Config.readInto(fn_randrough, "Terrain_Random_Roughness");
    Config.readInto(fn_slope, "Slope");
    Config.readInto(fn_fimperv,"Fraction_Impervious_Surface");
    Config.readInto(fn_poros0, "Top-of-profile_Porosity");
    Config.readInto(fn_psi_ae, "Air_entry_pressure");
    Config.readInto(fn_BClambda, "Brooks_Corey_lambda");
    Config.readInto(fn_theta_r, "Residual_soil_moisture");
    Config.readInto(fn_depth_layer1, "Depth_soil_layer_1");
    Config.readInto(fn_depth_layer2, "Depth_soil_layer_2");
    Config.readInto(fn_soildepth,    "Depth_soil_layer_3");
    Config.readInto(fn_paramWc, "Veget_water_use_param1");
    Config.readInto(fn_paramWp, "Veget_water_use_param2");
    Config.readInto(fn_snowCf, "Snow_Melt_Coeff");

    if(toggle_soil_prop == 1){
      if(sw_expKsat)
	Config.readInto(fn_kKsat, "Horiz_Hydraulic_Conductivity_Profile_Coeff");
      if(sw_expPoros)
	Config.readInto(fn_kporos, "Porosity_Profile_Coeff");
    } else if(toggle_soil_prop == 2){
      Config.readInto(fn_Ksat2, "Horiz_Hydraulic_Conductivity_L2");
      Config.readInto(fn_Ksat3, "Horiz_Hydraulic_Conductivity_L3");
      Config.readInto(fn_kvkh2, "Vert_Horz_Anis_ratio_L2");
      Config.readInto(fn_kvkh3, "Vert_Horz_Anis_ratio_L3");
      Config.readInto(fn_poros2, "Porosity_L2");
      Config.readInto(fn_poros3, "Porosity_L3");
      Config.readInto(fn_psi_ae2, "Air_entry_pressure_L2");      
      Config.readInto(fn_psi_ae3, "Air_entry_pressure_L3");
      Config.readInto(fn_BClambda2, "Brooks_Corey_lambda_L2");
      Config.readInto(fn_BClambda3, "Brooks_Corey_lambda_L3");  
    }

    snow_rain_temp = Config.read<float>("Snow_rain_temp_threshold");
    Config.readInto(fn_isohyet, "Isohyet_map");
    Config.readInto(fn_Ldown, "IncomingLongWave");
    Config.readInto(fn_Sdown, "IncomingShortWave");
    Config.readInto(fn_temp, "AirTemperature");
    Config.readInto(fn_maxTemp, "MaxAirTemp");
    Config.readInto(fn_minTemp, "MinAirTemp");
    Config.readInto(fn_precip, "Precipitation");
    Config.readInto(fn_rel_humid, "RelativeHumidity");
    Config.readInto(fn_wind_speed, "WindSpeed");
    Config.readInto(fn_pressure, "Pressure");
    if(sw_anthr_heat)
      Config.readInto(fn_anthrop_heat, "AnthropogenicHeat");

    Config.readInto(fn_ldd, "local_drain_direc");
    Config.readInto(fn_chwidth, "channel_width");
    Config.readInto(fn_chlength, "channel_length");
    Config.readInto(fn_chgwparam, "channel_gw_transfer_param");
    Config.readInto(fn_chmanningn, "mannings_n");

    if(toggle_chan_evap==1){
      Config.readInto(fn_temp_w,"Water_temperature");
      Config.readInto(fn_chanrough,"Channel_roughness");
    }
    /*Basin state variables section*/
    Config.readInto(fn_swe, "snow_water_equivalent");
    Config.readInto(fn_albedo, "Albedo");
    Config.readInto(fn_emiss, "Surface_emissivity");
    Config.readInto(fn_soiltemp, "Soil_temperature");
    Config.readInto(fn_soilheatcap, "Dry_Soil_Heat_Capacity");
    Config.readInto(fn_soilthermcond, "Dry_Soil_Therm_Cond");
    Config.readInto(fn_dampdepth, "Damping_depth");
    Config.readInto(fn_tempdamp, "Temp_at_damp_depth");
    Config.readInto(fn_streamflow, "Streamflow");
    Config.readInto(fn_soilmoist, "Soil_moisture_1");
    Config.readInto(fn_soilmoist2, "Soil_moisture_2");
    Config.readInto(fn_soilmoist3, "Soil_moisture_3");
    Config.readInto(fn_bedrock_leak, "Soil_bedrock_leakance");
    if(sw_deepGW){	
      Config.readInto(fn_chdeepgwparam, "channel_deepgw_transfer_param");
      Config.readInto(fn_deepGW, "Groundwater_DeepStorage");
      Config.readInto(fn_hydro_deepGW, "Fraction_Hydroactive_DeepGW");
    }

    Config.readInto(fn_paramtable, "Species_Parameters");

    if(!ForestStateVarsInputType.compare("tables")){
      Config.readInto(fn_proptable, "Species_Proportion_Table");
      Config.readInto(fn_LAItable, "Species_LAI_Table");
      Config.readInto(fn_AGEtable, "Species_AGE_Table");
      Config.readInto(fn_BasalAreatable, "Species_BasalArea_Table");
      Config.readInto(fn_Heighttable, "Species_Height_table");
      Config.readInto(fn_RootMasstable, "Species_RootMass_table");
      Config.readInto(fn_StemDenstable, "Species_StemDensity_Table");
    }

    current_t_step = current_ts_count * dt;

    //*************Reading the maps*************
    //-------------   Input Maps   -------------
    Rep_Long_Rad_Down = Config.read<bool>("Report_Long_Rad_Down");
    Rep_Short_Rad_Down = Config.read<bool>("Report_Short_Rad_Down");
    Rep_Precip = Config.read<bool>("Report_Precip");
    Rep_Rel_Humidity = Config.read<bool>("Report_Rel_Humidity");
    Rep_Wind_Speed = Config.read<bool>("Report_Wind_Speed");
    Rep_AvgAir_Temperature = Config.read<bool>("Report_AvgAir_Temperature");
    Rep_MinAir_Temperature = Config.read<bool>("Report_MinAir_Temperature");
    Rep_MaxAir_Temperature = Config.read<bool>("Report_MaxAir_Temperature");
    if(sw_anthr_heat)
      Rep_Anthropogenic_Heat = Config.read<bool>("Report_Anthropogenic_Heat");

    if(sw_BC){
      Rep_BC_surface = Config.read<bool>("Report_BC_Surface");
      Rep_BC_groundwater = Config.read<bool>("Report_BC_Groundwater");      
    }    
    //------------- Water Balance  -------------
    Rep_SWE = Config.read<bool>("Report_SWE");
    Rep_Infilt_Cap = Config.read<bool>("Report_Infilt_Cap");
    Rep_Streamflow = Config.read<bool>("Report_Streamflow");
    Rep_Saturation_Area = Config.read<bool>("Report_Saturation_Area");
    Rep_Ponding = Config.read<bool>("Report_Ponding");
    Rep_Soil_Water_Content_Average = Config.read<bool>("Report_Soil_Water_Content_Average");
    Rep_Soil_Water_Content_12 = Config.read<bool>("Report_Soil_Water_Content_Up");
    Rep_Soil_Water_Content_L1 = Config.read<bool>("Report_Soil_Water_Content_L1");
    Rep_Soil_Water_Content_L2 = Config.read<bool>("Report_Soil_Water_Content_L2");
    Rep_Soil_Water_Content_L3 = Config.read<bool>("Report_Soil_Water_Content_L3");
    Rep_WaterTableDepth = Config.read<bool>("Report_WaterTableDepth");
    Rep_Soil_Sat_Deficit = Config.read<bool>("Report_Soil_Sat_Deficit");
    Rep_GWater = Config.read<bool>("Report_Ground_Water");
    if(sw_deepGW)
      Rep_DeepGWater = Config.read<bool>("Report_Deep_Ground_Water");
    Rep_Canopy_Water_Stor_sum = Config.read<bool>("Report_Canopy_Water_Stor_sum");
    // Time-constant outputs: only reported once
    Rep_RootZone_in_L1 = Config.read<bool>("Report_RootZone_in_L1");
    Rep_RootZone_in_L2 = Config.read<bool>("Report_RootZone_in_L2");
    Rep_RootZone_in_L3 = Config.read<bool>("Report_RootZone_in_L3");
    Rep_Field_Capacity_L1 = Config.read<bool>("Report_Field_Capacity_L1");
    Rep_Field_Capacity_L2 = Config.read<bool>("Report_Field_Capacity_L2");
    Rep_Field_Capacity_L3 = Config.read<bool>("Report_Field_Capacity_L3");

    //------------- Energy Balance -------------
    Rep_Soil_Net_Rad = Config.read<bool>("Report_Soil_Net_Rad");
    Rep_Soil_LE = Config.read<bool>("Report_Soil_LE");
    Rep_Sens_Heat = Config.read<bool>("Report_Sens_Heat");
    Rep_Grnd_Heat = Config.read<bool>("Report_Grnd_Heat");
    Rep_Snow_Heat = Config.read<bool>("Report_Snow_Heat");
    Rep_Soil_Temperature = Config.read<bool>("Report_Soil_Temperature");
    Rep_SoilL1_Temperature = Config.read<bool>("Report_SoilL1_Temperature");
    Rep_SoilL2_Temperature = Config.read<bool>("Report_SoilL2_Temperature");
    Rep_SoilL3_Temperature = Config.read<bool>("Report_SoilL3_Temperature");    
    Rep_Skin_Temperature = Config.read<bool>("Report_Skin_Temperature");
    if(toggle_chan_evap == 1)
      Rep_Water_Temperature = Config.read<bool>("Report_Water_Temperature");
    Rep_Net_Rad_sum = Config.read<bool>("Report_Net_Rad_sum");
    Rep_LE_sum = Config.read<bool>("Report_LE_sum");
    Rep_H_sum = Config.read<bool>("Report_H_sum");
    // Maps of species specific energy
    Rep_Canopy_Temp = Config.read<bool>("Report_Canopy_Temp");
    Rep_Canopy_NetR = Config.read<bool>("Report_Canopy_NetR");
    Rep_Canopy_LE_E = Config.read<bool>("Report_Canopy_LE_E");
    Rep_Canopy_LE_T = Config.read<bool>("Report_Canopy_LE_T");
    Rep_Canopy_Sens_Heat = Config.read<bool>("Report_Canopy_Sens_Heat");

    //--------------- ET Outputs -------------
    Rep_Total_ET = Config.read<bool>("Report_Total_ET");
    Rep_Transpiration_sum = Config.read<bool>("Report_Transpiration_sum");
    Rep_Einterception_sum = Config.read<bool>("Report_Einterception_sum");
    Rep_Esoil_sum = Config.read<bool>("Report_Esoil_sum");
    if(toggle_chan_evap > 0)
      Rep_ChanEvap = Config.read<bool>("Report_ChannelE");      
    // Maps of species specific ET
    Rep_ETspecies = Config.read<bool>("Report_species_ET");
    Rep_Transpiration = Config.read<bool>("Report_Transpiration");
    Rep_Einterception = Config.read<bool>("Report_Einterception");
    Rep_Esoil = Config.read<bool>("Report_Esoil");

    //--------------- Vegetation Input Maps -------------
    Rep_Veget_frac = Config.read<bool>("Report_Veget_frac");
    Rep_Stem_Density = Config.read<bool>("Report_Stem_Density");
    Rep_RootFrac1Species = Config.read<bool>("Report_RootFracL1_species");
    Rep_RootFrac2Species = Config.read<bool>("Report_RootFracL2_species");

    //--------------- Vegetation Output Maps -------------
    Rep_Leaf_Area_Index = Config.read<bool>("Report_Leaf_Area_Index");
    Rep_Stand_Age = Config.read<bool>("Report_Stand_Age");
    Rep_Canopy_Conductance = Config.read<bool>("Report_Canopy_Conductance");
    Rep_GPP = Config.read<bool>("Report_GPP");
    Rep_NPP = Config.read<bool>("Report_NPP");
    Rep_Basal_Area = Config.read<bool>("Report_Basal_Area");
    Rep_Tree_Height = Config.read<bool>("Report_Tree_Height");
    Rep_Root_Mass = Config.read<bool>("Report_Root_Mass");
    Rep_Canopy_Water_Stor = Config.read<bool>("Report_Canopy_Water_Stor");
    Rep_RUptake_L1 = Config.read<bool>("Report_RootUptake_Prop_L1");
    Rep_RUptake_L2 = Config.read<bool>("Report_RootUptake_Prop_L2");
    Rep_RUptake_L3 = Config.read<bool>("Report_RootUptake_Prop_L3");    
    Rep_SoilPot = Config.read<bool>("Report_Soil_Water_Potential");
    Rep_VegPot = Config.read<bool>("Report_Leaf_Water_Potential");
    Rep_SapVel = Config.read<bool>("Report_Sapflow_Velocity");    

    //------------- Water Balance Output  -------------
    Rep_GWtoChn = Config.read<bool>("Report_GW_to_Channel");
    if(sw_deepGW)
      Rep_DeepGWtoChn = Config.read<bool>("Report_DeepGW_to_Channel");
    Rep_SrftoChn = Config.read<bool>("Report_Surface_to_Channel");
    Rep_SrftoLat = Config.read<bool>("Report_Overland_Outflow");
    Rep_ChntoLat = Config.read<bool>("Report_Stream_Outflow");
    Rep_GWtoLat = Config.read<bool>("Report_Groundwater_Outflow");
    if(sw_deepGW)
      Rep_DeepGWtoLat = Config.read<bool>("Report_DeepGW_Outflow");
    Rep_Exfilt = Config.read<bool>("Report_Return_Flow_Surface");
    Rep_PercolL2 = Config.read<bool>("Report_Percolation_to_Layer2");
    Rep_PercolL3 = Config.read<bool>("Report_Percolation_to_Layer3");
    Rep_Recharge = Config.read<bool>("Report_Groundwater_Recharge");
    Rep_Leak     = Config.read<bool>("Report_Bedrock_Leakance");
    Rep_LattoSrf = Config.read<bool>("Report_Overland_Inflow");
    Rep_LattoChn = Config.read<bool>("Report_Stream_Inflow");
    Rep_LattoGW = Config.read<bool>("Report_Groundwater_Inflow");
    if(sw_deepGW)
      Rep_LattoDeepGW = Config.read<bool>("Report_DeepGW_Inflow");
    Rep_Infilt = Config.read<bool>("Report_Infiltration");
    Rep_ReturnL1 = Config.read<bool>("Report_Return_Flow_to_Layer1");
    Rep_ReturnL2 = Config.read<bool>("Report_Return_Flow_to_Layer2");

    //------------- Cumulative -------------
    Rep_Infiltacc = Config.read<bool>("Report_Infiltration_acc");
    Rep_Exfiltacc = Config.read<bool>("Report_Return_Flow_Surface_acc");
    Rep_PercolL2acc = Config.read<bool>("Report_Percolation_to_Layer2_acc");
    Rep_ReturnL1acc = Config.read<bool>("Report_Return_Flow_to_Layer1_acc");
    Rep_PercolL3acc = Config.read<bool>("Report_Percolation_to_Layer3_acc");
    Rep_Rechargeacc = Config.read<bool>("Report_Groundwater_Recharge_acc");
    Rep_ReturnL2acc = Config.read<bool>("Report_Return_Flow_to_Layer2_acc");
    Rep_EvaporationSacc = Config.read<bool>("Report_Soil_Evaporation_acc");
    Rep_LattoSrfacc = Config.read<bool>("Report_Overland_Inflow_acc");
    Rep_LattoChnacc = Config.read<bool>("Report_Stream_Inflow_acc");
    Rep_LattoGWacc = Config.read<bool>("Report_Groundwater_Inflow_acc");
    if(sw_deepGW)
      Rep_LattoDeepGWacc = Config.read<bool>("Report_DeepGW_Inflow_acc");
    Rep_SrftoLatacc = Config.read<bool>("Report_Overland_Outflow_acc");
    Rep_ChntoLatacc = Config.read<bool>("Report_Stream_Outflow_acc");
    Rep_GWtoLatacc = Config.read<bool>("Report_Groundwater_Outflow_acc");
    if(sw_deepGW)
      Rep_DeepGWtoLatacc = Config.read<bool>("Report_DeepGW_Outflow_acc");
    Rep_GWtoChnacc = Config.read<bool>("Report_GW_to_Channel_acc");
    if(sw_deepGW)
      Rep_DeepGWtoChnacc = Config.read<bool>("Report_DeepGW_to_Channel_acc");
    Rep_SrftoChnacc = Config.read<bool>("Report_Surface_to_Channel_acc");

    //*************Reading the timeseries*************
    Config.readInto(fn_rep_mask, "TS_mask");
    //-------------   Input Maps   -------------
    RepTs_OutletDischarge = Config.read<bool>("Ts_OutletDischarge");
    RepTs_OutletGW = Config.read<bool>("Ts_OutletGW");
    RepTs_Long_Rad_Down = Config.read<bool>("Ts_Long_Rad_Down");
    RepTs_Short_Rad_Down = Config.read<bool>("Ts_Short_Rad_Down");
    RepTs_Precip = Config.read<bool>("Ts_Precip");
    RepTs_Rel_Humidity = Config.read<bool>("Ts_Rel_Humidity");
    RepTs_Wind_Speed = Config.read<bool>("Ts_Wind_Speed");
    RepTs_AvgAir_Temperature = Config.read<bool>("Ts_AvgAir_Temperature");
    RepTs_MinAir_Temperature = Config.read<bool>("Ts_MinAir_Temperature");
    RepTs_MaxAir_Temperature = Config.read<bool>("Ts_MaxAir_Temperature");
    if(sw_anthr_heat)
      RepTs_Anthropogenic_Heat = Config.read<bool>("Ts_Anthropogenic_Heat");

    //------------- Water Balance  -------------
    RepTs_SWE = Config.read<bool>("Ts_SWE");
    RepTs_Infilt_Cap = Config.read<bool>("Ts_Infilt_Cap");
    RepTs_Streamflow = Config.read<bool>("Ts_Streamflow");
    RepTs_Ponding = Config.read<bool>("Ts_Ponding");
    RepTs_Soil_Water_Content_Average = Config.read<bool>("Ts_Soil_Water_Content_Average");
    RepTs_Soil_Water_Content_12 = Config.read<bool>("Ts_Soil_Water_Content_Up");
    RepTs_Soil_Water_Content_L1 = Config.read<bool>("Ts_Soil_Water_Content_L1");
    RepTs_Soil_Water_Content_L2 = Config.read<bool>("Ts_Soil_Water_Content_L2");
    RepTs_Soil_Water_Content_L3 = Config.read<bool>("Ts_Soil_Water_Content_L3");
    RepTs_WaterTableDepth = Config.read<bool>("Ts_WaterTableDepth");
    RepTs_Soil_Sat_Deficit = Config.read<bool>("Ts_Soil_Sat_Deficit");
    RepTs_GroundWater = Config.read<bool>("Ts_Ground_Water");
    if(sw_deepGW)
      RepTs_DeepGWater = Config.read<bool>("Ts_Deep_Ground_Water");
    RepTs_Canopy_Water_Stor_sum = Config.read<bool>("Ts_Canopy_Water_Stor_sum");
    // Time-constant outputs: only reported once
    RepTs_Field_Capacity_L1 = Config.read<bool>("Ts_Field_Capacity_L1");
    RepTs_Field_Capacity_L2 = Config.read<bool>("Ts_Field_Capacity_L2");
    RepTs_Field_Capacity_L3 = Config.read<bool>("Ts_Field_Capacity_L3");

    //------------- Energy Balance -------------
    RepTs_Soil_Net_Rad = Config.read<bool>("Ts_Soil_Net_Rad");
    RepTs_Soil_LE = Config.read<bool>("Ts_Soil_LE");
    RepTs_Sens_Heat = Config.read<bool>("Ts_Sens_Heat");
    RepTs_Grnd_Heat = Config.read<bool>("Ts_Grnd_Heat");
    RepTs_Snow_Heat = Config.read<bool>("Ts_Snow_Heat");
    RepTs_Soil_Temperature = Config.read<bool>("Ts_Soil_Temperature");
    RepTs_SoilL1_Temperature = Config.read<bool>("Ts_SoilL1_Temperature");
    RepTs_SoilL2_Temperature = Config.read<bool>("Ts_SoilL2_Temperature");
    RepTs_SoilL3_Temperature = Config.read<bool>("Ts_SoilL3_Temperature");        
    RepTs_Skin_Temperature = Config.read<bool>("Ts_Skin_Temperature");
    if(toggle_chan_evap == 1)
      RepTs_Water_Temperature = Config.read<bool>("Ts_Water_Temperature");
    RepTs_Net_Rad_sum = Config.read<bool>("Ts_Net_Rad_sum");
    RepTs_LE_sum = Config.read<bool>("Ts_LE_sum");
    RepTs_H_sum = Config.read<bool>("Ts_H_sum");   
    // Maps of species specific energy
    RepTs_Canopy_Temp = Config.read<bool>("Ts_Canopy_Temp");
    RepTs_Canopy_NetR = Config.read<bool>("Ts_Canopy_NetR");
    RepTs_Canopy_LE_E = Config.read<bool>("Ts_Canopy_LE_E");
    RepTs_Canopy_LE_T = Config.read<bool>("Ts_Canopy_LE_T");
    RepTs_Canopy_Sens_Heat = Config.read<bool>("Ts_Canopy_Sens_Heat");

    //--------------- ET Outputs -------------
    RepTs_Total_ET = Config.read<bool>("Ts_Total_ET");
    RepTs_Transpiration_sum = Config.read<bool>("Ts_Transpiration_sum");
    RepTs_Einterception_sum = Config.read<bool>("Ts_Einterception_sum");
    RepTs_Esoil_sum = Config.read<bool>("Ts_Esoil_sum");
    if(toggle_chan_evap > 0)
      RepTs_ChanEvap = Config.read<bool>("Ts_ChannelE");  
    // Maps of species specific ET
    RepTs_ETspecies = Config.read<bool>("Ts_species_ET");
    RepTs_Transpiration = Config.read<bool>("Ts_Transpiration");
    RepTs_Einterception = Config.read<bool>("Ts_Einterception");
    RepTs_Esoil = Config.read<bool>("Ts_Esoil");

    //--------------- Vegetation Input Maps -------------
    RepTs_Veget_frac = Config.read<bool>("Ts_Veget_frac");
    RepTs_Stem_Density = Config.read<bool>("Ts_Stem_Density");
    RepTs_RootFrac1Species = Config.read<bool>("Ts_RootFracL1_species");
    RepTs_RootFrac2Species = Config.read<bool>("Ts_RootFracL2_species");

    //--------------- Vegetation Output Maps -------------
    RepTs_Leaf_Area_Index = Config.read<bool>("Ts_Leaf_Area_Index");
    RepTs_Stand_Age = Config.read<bool>("Ts_Stand_Age");
    RepTs_Canopy_Conductance = Config.read<bool>("Ts_Canopy_Conductance");
    RepTs_GPP = Config.read<bool>("Ts_GPP");
    RepTs_NPP = Config.read<bool>("Ts_NPP");
    RepTs_Basal_Area = Config.read<bool>("Ts_Basal_Area");
    RepTs_Tree_Height = Config.read<bool>("Ts_Tree_Height");
    RepTs_Root_Mass = Config.read<bool>("Ts_Root_Mass");
    RepTs_Canopy_Water_Stor = Config.read<bool>("Ts_Canopy_Water_Stor");
    RepTs_RUptake_L1 = Config.read<bool>("Ts_RootUptake_Prop_L1");
    RepTs_RUptake_L2 = Config.read<bool>("Ts_RootUptake_Prop_L2");
    RepTs_RUptake_L3 = Config.read<bool>("Ts_RootUptake_Prop_L3");  
    RepTs_SoilPot = Config.read<bool>("Ts_Soil_Water_Potential");
    RepTs_VegPot = Config.read<bool>("Ts_Leaf_Water_Potential");
    RepTs_SapVel = Config.read<bool>("Ts_Sapflow_Velocity");

    //------------- Water Balance Output  -------------
    // Internal vertical fluxes
    RepTs_GWtoChn = Config.read<bool>("Ts_GW_to_Channel");
    if(sw_deepGW)
      RepTs_DeepGWtoChn = Config.read<bool>("Ts_DeepGW_to_Channel");
    RepTs_SrftoChn = Config.read<bool>("Ts_Surface_to_Channel");
    RepTs_SrftoLat = Config.read<bool>("Ts_Overland_Outflow");
    RepTs_ChntoLat = Config.read<bool>("Ts_Stream_Outflow");
    RepTs_GWtoLat = Config.read<bool>("Ts_Groundwater_Outflow");
    if(sw_deepGW)
      RepTs_DeepGWtoLat = Config.read<bool>("Ts_DeepGW_Outflow");
    RepTs_PercolL2 = Config.read<bool>("Ts_Percolation_to_Layer2");
    RepTs_PercolL3 = Config.read<bool>("Ts_Percolation_to_Layer3");
    RepTs_Recharge = Config.read<bool>("Ts_Groundwater_Recharge");
    RepTs_Leak     = Config.read<bool>("Ts_Bedrock_Leakance");
    RepTs_LattoSrf = Config.read<bool>("Ts_Overland_Inflow");
    RepTs_LattoChn = Config.read<bool>("Ts_Stream_Inflow");
    RepTs_LattoGW = Config.read<bool>("Ts_Groundwater_Inflow");
    if(sw_deepGW)
      RepTs_LattoDeepGW = Config.read<bool>("Ts_DeepGW_Inflow");
    RepTs_Infilt = Config.read<bool>("Ts_Infiltration");
    RepTs_Exfilt = Config.read<bool>("Ts_Return_Flow_Surface");
    RepTs_ReturnL1 = Config.read<bool>("Ts_Return_Flow_to_Layer1");
    RepTs_ReturnL2 = Config.read<bool>("Ts_Return_Flow_to_Layer2");

  }
  catch(ConfigFile::file_not_found &fn){
    cout << "File " << fn.filename << " not found\n";
    exit(EXIT_SUCCESS);
  }
  catch(ConfigFile::key_not_found &fn){
    cout << "Key " << fn.key << " not found\n";
    exit(EXIT_SUCCESS);
  }
  catch(std::exception &e){
    cout << e.what();
    exit(EXIT_SUCCESS);
  }




  return 1;
}
