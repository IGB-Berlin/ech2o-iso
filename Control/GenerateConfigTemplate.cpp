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
 * GenerateConfigTemplate.cpp
 *
 *  Created on: Jan 12, 2016
 *      Author: Marco Maneta
 */

#include  <fstream>
#include <unistd.h>

#include "Sativa.h"

void GenerateConfigTemplate(const char *fn){

  ofstream ofOut;

  try{

    if (access(fn, F_OK) != -1) {

      cout << "File exists. Do you want to overwrite? (y, n):  " << endl;
      char c;
      cin.get(c);
      switch (c) {
      case 'y':
	break;
      case 'n':
	exit(EXIT_SUCCESS);
	break;
      default:
	cout << "Not a valid option. Bye" << endl;
	exit(EXIT_SUCCESS);

      }
    }

    ofOut.open(fn);
    if(!ofOut)
      throw std::ios::failure("Error opening file ");


    ofOut << "#ECH2O configuration file v3" << std::endl << std::endl;
    ofOut << "# Please, check Appendix A of Documentation" << std::endl;
    ofOut << "# for units of parameters and variables  " << std::endl;
    ofOut << "# (http://ech2o-iso.readthedocs.io/en/latest/Keywords.html)" << std::endl << std::endl;

    ofOut << "#" << endl << "#Folder section" << endl << "#" << endl << endl;

    ofOut << "Maps_Folder = ./Spatial" << endl;
    ofOut << "Clim_Maps_Folder = ./Climate" << endl;
    ofOut << "Output_Folder = ./Outputs" << endl << endl;

    ofOut << "#" << endl << "#Water tracking (isotopes and/or ages)" << endl;
    ofOut << "Tracking = 1" << endl ;
    ofOut << "TrackingConfig = ./configTrck.ini" << endl << endl; 

    ofOut << "#" << endl << "# Options section" << endl << "#" << endl << endl;
    
    ofOut << "MapTypes = csf" << endl;
    ofOut << "Species_State_Variable_Input_Method = maps # maps or tables" << endl << endl;
    ofOut << "closure_tolerance = 1e-12" << endl << endl;

    ofOut << "#***************************************" << endl;      
    ofOut << "#== Boolean switches" << endl;
    ofOut << "#***************************************" << endl;  
    ofOut << "Reinfiltration = 1" << endl;
    ofOut << "Channel = 1" << endl ;
    ofOut << "Anthropogenic_heat = 0"<< endl;    
    ofOut << "Interception_proportion = 0" <<endl << endl;

    ofOut << "#***************************************" << endl;  
    ofOut << "# TOGGLE SWITCHES:" << endl << endl;
    ofOut << "#***************************************" << endl << endl;    
    ofOut << "#---------------------------------------" << endl;
    ofOut << "#Vegetation growth and dynamics options" << endl;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "#Vegetation dynamics (via allocation): 3 modes" << endl;
    ofOut << "# 0 -> deactivated, no dynamic allocation and constant LAI to initial value" <<endl;
    ofOut << "# 1 -> fully activated" << endl;
    ofOut << "# 2 -> partially activated, (LAI via an input file)" << endl;
    ofOut << "Vegetation_dynamics = 0 " << endl;
    ofOut << "# Used only if Vegetation_dynamics = 2. Files names for each species is" << endl;
    ofOut << "# name below + '_'+ species number (starting at 0) + '.bin'" << endl;
    ofOut << "TimeSeries_LAI = LAI" << endl;
    ofOut << "TimeSeries_Height = Height" << endl;    
    ofOut << "# Stomatal model choices: " << endl;
    ofOut << "# 0 -> Non-linear option (Weibull function) " << endl;
    ofOut << "# 1 -> Linear option " << endl;
    ofOut << "Stomatal_model_opt = 0 " << endl;
    ofOut << "# Vegetation hydraulics (limit of vegetation conductance)" << endl;
    ofOut << "# 0 -> deactivated, only soil controlled" << endl;
    ofOut << "# 1 -> fully activated, soil + xylem controlled" << endl;
    ofOut << "Plant_hydraulics = 1 " << endl << endl;

    ofOut << "#---------------------------------------" << endl;
    ofOut << "#Soil and hydraulics options" << endl;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "# Hydrologic engine choices " << endl;
    ofOut << "# 0 -> Green-Ampt + Gravity flow " << endl;
    //    ofOut << "# 1 -> Richard's equation + instantaneous overland flow" << endl;
    //    ofOut << "# 2 -> Richard's equation + kinematic overland flow" << endl;
    ofOut << "Hydrologic_engine_opt = 0 " << endl;
    ofOut << "# Channel evaporation option (on or off + options)" << endl;
    ofOut << "# 0 -> No channel evaporation "<< endl;
    ofOut << "# 1 -> Energy balance channel evaporation "<< endl;
    ofOut << "# 2 -> Mass transfer channel evaporation "<< endl;        
    ofOut << "Channel_Evaporation = 2 " << endl;
    ofOut << "# Aerodynamic resistance choices: " << endl;
    ofOut << "# 0 = Penman Monteith option " << endl;
    ofOut << "# 1 = Thom and Oliver 1977 " << endl;
    ofOut << "Aerodyn_resist_opt = 0 " << endl;
    ofOut << "# Soil resistance to vapor diffusion choices: " << endl;
    ofOut << "# 0 = No resistance" << endl;
    ofOut << "# 1 = Passerat de Silans et al. 1989 (bare soil)" << endl;
    ofOut << "# 2 = Sellers et al. 1992 (grass)" << endl;
    ofOut << "# 3 = Sakaguchi and Zeng 2009 (CLM 3.5) (mixed vegetation)" << endl;
    ofOut << "# 4 = Dynamic (models 1-3) based on vegetation type" << endl;
    ofOut << "Soil_resistance_opt = 3 " << endl << endl;

    ofOut << "#---------------------------------------" << endl;
    ofOut << "#Additional boundary condition options "<< endl;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "# Boundary conditions options" << endl;
    ofOut << "# 0 -> Boundary conditions not applied" << endl;
    ofOut << "# 1 -> Boundary conditions maps read in" << endl;
    ofOut << "Boundary_Condition = 0 " << endl;
    ofOut << "# Used only if Boundary_Condition = 1. File names for each input: " << endl;
    ofOut << "#   Units: Surface(m3/s) - GW (m2/s) " << endl;    
    ofOut << "TimeSeries_BC_Surface = BCsurface.bin" << endl;
    ofOut << "TimeSeries_BC_Groundwater = BCgroundwater.bin" << endl;
    ofOut << "TimeSeries_BC_DeepGW = BCdeepgroundwater.bin" << endl << endl;
    
    ofOut << "#---------------------------------------" << endl;
    ofOut << "#Soil parameterisation options" << endl;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "# Soil parameterisation for soil layers: " << endl;
    ofOut << "# 0 = Define uniform soil parameters for all layers" << endl;
    ofOut << "# 1 = Defined at the surface (exponential decrease in depth) " << endl;
    ofOut << "# 2 = Define soil parameters for each depth " << endl;
    ofOut << "Soil_Properties = 0 " << endl;
    ofOut << "# Exponential profiles options: (0=vertically-uniform using" << endl;
    ofOut << "# the top-of-profile (see corresponding map inputs below)" << endl;
    ofOut << "Hydraulic_Conductivity_profile = 0" << endl ;
    ofOut << "Porosity_profile = 0" << endl;
    ofOut << "# Activation of a deeper groundwater storage contribution" << endl;
    ofOut << "DeepGW_Storage	= 0" << endl << endl;
    
    ofOut << "#***************************************" << endl;        
    ofOut << "#" << endl << "# Time variables section" << endl << "#" << endl;
    ofOut << "#***************************************" << endl;    
    ofOut << "Simul_start = 0 # Must be less than Simul_end" << endl;
    ofOut << "Simul_end = 31536000 # seconds (365 days)" << endl;
    ofOut << "Simul_tstep = 86400 # seconds (daily)" << endl;
    ofOut << "Clim_input_tstep = 86400 # seconds (daily)" << endl;
    ofOut << "Report_interval = 86400 # seconds (daily)" << endl ;
    ofOut << "ReportMap_interval = 86400 # seconds (daily)" << endl ;
    ofOut << "ReportMap_starttime = 86400 # seconds (from first time step)" << endl;
    ofOut << "NetCDF_output_format = 0" << endl << endl;

    ofOut << "#***************************************" << endl;        
    ofOut << "#" << endl << "# Climate input information" << endl;
    ofOut << "# Maps in this section to be contained in folder pointed by Clim_Maps_Folder" << endl << "#" << endl;
    ofOut << "#***************************************" << endl;        
    ofOut << "ClimateZones = ClimZones.map" << endl;
    ofOut << "Isohyet_map = isohyet.map " << " # Precipitation multiplier map"<< endl;
    ofOut << "Snow_rain_temp_threshold = 2 " << " # Snow to rain temperatures threshold in degC" << endl;
    ofOut << "Precipitation = Precip.bin " << " # Precip rate in meters/second"<< endl;
    ofOut << "AirTemperature = Tavg.bin " << " # Average air temperature in degC" << endl;
    ofOut << "MaxAirTemp = Tmax.bin " << " # Maximum air temperature in degC" << endl;
    ofOut << "MinAirTemp = Tmin.bin " << " # Minimum air temperature in degC"<< endl;
    ofOut << "RelativeHumidity = RH.bin " << " # air relative humidity in kPa/kPa"<< endl;
    ofOut << "WindSpeed = windspeed.bin " << " # Wind speed in meters/second" << endl;
    ofOut << "Pressure = Pressure.bin " << " # air pressure in Pa" << endl;
    ofOut << "IncomingLongWave = Ldown.bin " << " # Downwelling longwave radiation in W/sq.meter" << endl;
    ofOut << "IncomingShortWave = Sdown.bin " << " # Solar radiation in W/sq.meter" << endl << endl;
    ofOut << "AnthropogenicHeat = W.bin" << " # Anthropogenic heat soure/sink in W/sq.meter (only if 'Anthropogenic_heat = 1')" << endl << endl;
   
    ofOut << "#***************************************" << endl;    
    ofOut << "#" << endl << "# Spatial input information" << endl;
    ofOut << "# Maps below this line to be contained in folder pointed by Maps_Folder" << endl;
    ofOut << "#" << endl;
    ofOut << "#***************************************" << endl;    

    ofOut << "#---------------------------------------" << endl;
    ofOut << "#" << endl << "# Drainage network" << endl << "#" << endl;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "DEM = DEM.map" << endl;
    ofOut << "Slope = slope.map " << endl;
    ofOut << "local_drain_direc = ldd.map" << endl;
    ofOut << "Fraction_Contributing_Area = fcontrea.map " << endl;
    ofOut << "Fraction_Impervious_Surface = fImperv.map" << endl;
    ofOut << "channel_width = chanwidth.map" << endl;
    ofOut << "channel_length = chanlength.map" << endl;
    ofOut << "channel_gw_transfer_param = chanparam.map" << endl;
    ofOut << "channel_deepgw_transfer_param = chanDeepparam.map" << endl;
    ofOut << "mannings_n = chanmanningn.map" << endl;
    ofOut << "# Channel evaporation parameters - only if energy balance channel evaporation is on" << endl;
    ofOut << "Water_temperature = water_temp.map" << endl;
    ofOut << "Channel_roughness = chanrough.map" << endl << endl;

    ofOut << "#---------------------------------------" << endl;
    ofOut << "#   " << endl << "# Hydrologic Initial Conditions  " << endl;
    ofOut << "# Forest Initial states are included as maps or tables" << endl << "#   " << endl;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "Streamflow = streamflow.map " << endl;
    ofOut << "snow_water_equivalent = swe.map " << endl;
    ofOut << "Soil_moisture_1 = SWC.L1.map " << endl;
    ofOut << "Soil_moisture_2 = SWC.L2.map " << endl;
    ofOut << "Soil_moisture_3 = SWC.L3.map " << endl;
    ofOut << "Soil_temperature = soiltemp.map " << endl << endl;
    ofOut << "Groundwater_DeepStorage = GW_DeepStorage.map " << endl << endl;

    ofOut << "#---------------------------------------" << endl;
    ofOut << "#   " << endl << "#Soil parameters  " << endl << "#   " << endl;
    ofOut << "#---------------------------------------" << endl;   
    ofOut << "Soil_Skin_Infilt_Capacity = KeffTopSoil.map " << endl;
    ofOut << "Top-of-profile_Horiz_Hydraulic_Conductivity = Keff.map " << endl;
    ofOut << "Horiz_Hydraulic_Conductivity_Profile_Coeff = kKsat.map " << endl;
    ofOut << "Vert_Horz_Anis_ratio = KvKh.map " << endl;
    ofOut << "Terrain_Random_Roughness = randrough.map " << endl;
    ofOut << "Top-of-profile_Porosity = poros.map " << endl;
    ofOut << "Porosity_Profile_Coeff = kporos.map " << endl;
    ofOut << "Air_entry_pressure = psi_ae.map " << endl;
    ofOut << "Brooks_Corey_lambda = BClambda.map " << endl;
    ofOut << "Residual_soil_moisture = theta_r.map " << endl;
    ofOut << "Depth_soil_layer_1 = soildepth.L1.map " << endl;
    ofOut << "Depth_soil_layer_2 = soildepth.L2.map " << endl;
    ofOut << "Depth_soil_layer_3 = soildepth.L3.map " << endl;
    ofOut << "Veget_water_use_param1 = Wc.map " << endl;
    ofOut << "Veget_water_use_param2 = Wp.map " << endl;
    ofOut << "Soil_bedrock_leakance = leakance.map " << endl << endl;
    ofOut << "Fraction_Hydroactive_DeepGW = fActive_DeepGW.map " << endl << endl; 
    ofOut << "Albedo = albedo.map" << endl;
    ofOut << "Surface_emissivity = emissivity.map" << endl;
    ofOut << "Dry_Soil_Heat_Capacity = soilheatcap.map" << endl;
    ofOut << "Dry_Soil_Therm_Cond = soilthermalK.map" << endl;
    ofOut << "Damping_depth = dampdepth.map" << endl;
    ofOut << "Temp_at_damp_depth = temp_damp.map" << endl;
    ofOut << "Snow_Melt_Coeff = snowmeltCoeff.map" << endl << endl;

    ofOut << "#---------------------------------------" << endl;
    ofOut << "#    " << endl << "#Depth-dependent soil parameters (only if activated)" << endl << "#" << endl;
    ofOut << "#---------------------------------------" << endl; 
    ofOut << "Horiz_Hydraulic_Conductivity_L2 = Keff.L2.mat" << endl;
    ofOut << "Horiz_Hydraulic_Conductivity_L3 = Keff.L3.mat" << endl;
    ofOut << "Vert_Horz_Anis_ratio_L2 = KvKh.L2.map " << endl;
    ofOut << "Vert_Horz_Anis_ratio_L3 = KvKh.L3.map " << endl;
    ofOut << "Porosity_L2 = poros.L2.map" << endl;
    ofOut << "Porosity_L3 = poros.L3.map" << endl;
    ofOut << "Air_entry_pressure_L2 = psi_ae.L2.map " << endl;
    ofOut << "Air_entry_pressure_L3 = psi_ae.L3.map " << endl;
    ofOut << "Brooks_Corey_lambda_L2 = BClambda.L2.map " << endl;
    ofOut << "Brooks_Corey_lambda_L3 = BClambda.L3.map " << endl << endl;

    ofOut << "#---------------------------------------" << endl;
    ofOut << "#   " << endl << "#Forest Parameters and initial states " << endl << "#   " << endl;
    ofOut << "#---------------------------------------" << endl; 
    ofOut << "ForestPatches = patches.map" << endl;
    ofOut << "Number_of_Species = 1 " << endl;
    ofOut << "Species_Parameters = SpeciesParams.tab " << endl << endl;
    ofOut << "#Tables below are only needed if Species_State_Variable_Input_Method = tables " << endl;
    ofOut << "Species_Proportion_Table = SpecsProp.tab " << endl;
    ofOut << "Species_StemDensity_Table = SpecsStemDens.tab " << endl;
    ofOut << "Species_LAI_Table = SpecsLAI.tab " << endl;
    ofOut << "Species_AGE_Table = SpecsAge.tab " << endl;
    ofOut << "Species_BasalArea_Table = SpeciesBasalArea.tab " << endl;
    ofOut << "Species_Height_table = SpeciesHeight.tab " << endl;
    ofOut << "Species_RootMass_table = SpecsRootDensity.tab " << endl << endl;

    ofOut << "#***************************************" << endl;        
    ofOut << "#Report map section " << endl;
    ofOut << "#***************************************" << endl << endl;
    ofOut << "#-------------------------------" << endl;
    ofOut << "#Input Maps " << endl;
    ofOut << "#-------------------------------" << endl;
    ofOut << "Report_Long_Rad_Down = 0 " << endl;
    ofOut << "Report_Short_Rad_Down = 0 " << endl;
    ofOut << "Report_Precip = 0 " << endl;
    ofOut << "Report_Rel_Humidity = 0 " << endl;
    ofOut << "Report_Wind_Speed = 0 " << endl;
    ofOut << "Report_AvgAir_Temperature = 0 " << endl;
    ofOut << "Report_MinAir_Temperature = 0 " << endl;
    ofOut << "Report_MaxAir_Temperature = 0 " << endl;
    ofOut << "Report_Anthropogenic_Heat = 0 " << endl << endl;

    ofOut << "Report_BC_Surface = 0" << endl;
    ofOut << "Report_BC_Groundwater = 0" << endl << endl;
    
    ofOut << "#-------------------------------" << endl;
    ofOut << "#Water Balance (Storage) Output Maps " << endl;
    ofOut << "#-------------------------------" << endl;   
    ofOut << "Report_SWE = 0 " << endl;
    ofOut << "Report_Infilt_Cap = 0 " << endl;
    ofOut << "Report_Streamflow = 0 " << endl;
    ofOut << "Report_Saturation_Area = 0 " << endl;
    ofOut << "Report_Ponding = 0 " << endl;
    ofOut << "Report_Soil_Water_Content_Average = 0 " << endl;
    ofOut << "Report_Soil_Water_Content_Up = 0 " << endl;
    ofOut << "Report_Soil_Water_Content_L1 = 0 " << endl;
    ofOut << "Report_Soil_Water_Content_L2 = 0 " << endl;
    ofOut << "Report_Soil_Water_Content_L3 = 0 " << endl;
    ofOut << "Report_WaterTableDepth = 0 " << endl;
    ofOut << "Report_Soil_Sat_Deficit = 0 " << endl;
    ofOut << "Report_Ground_Water = 0 " << endl;
    ofOut << "Report_Deep_Ground_Water = 0 " << endl;
    ofOut << "Report_Canopy_Water_Stor_sum = 0 " << endl;
    ofOut << "# Maps of time-constant variables (only reported once) --" << endl;
    ofOut << "Report_RootZone_in_L1 = 0 " << endl;
    ofOut << "Report_RootZone_in_L2 = 0 " << endl;
    ofOut << "Report_RootZone_in_L3 = 0 " << endl;
    ofOut << "Report_Field_Capacity_L1 = 0 " << endl;
    ofOut << "Report_Field_Capacity_L2 = 0 " << endl;
    ofOut << "Report_Field_Capacity_L3 = 0 " << endl << endl;

    ofOut << "#-------------------------------" << endl;
    ofOut << "#Energy Balance Output Maps " << endl;
    ofOut << "#-------------------------------" << endl;   
    ofOut << "Report_Soil_Net_Rad = 0 " << endl;
    ofOut << "Report_Soil_LE = 0 " << endl;
    ofOut << "Report_Sens_Heat = 0 " << endl;
    ofOut << "Report_Grnd_Heat = 0 " << endl;
    ofOut << "Report_Snow_Heat = 0 " << endl;
    ofOut << "Report_Soil_Temperature = 0 " << endl;
    ofOut << "Report_SoilL1_Temperature = 0 " << endl;
    ofOut << "Report_SoilL2_Temperature = 0 " << endl;
    ofOut << "Report_SoilL3_Temperature = 0 " << endl;  
    ofOut << "Report_Skin_Temperature = 0 " << endl;
    ofOut << "Report_Water_Temperature = 0 " << endl;
    ofOut << "Report_Net_Rad_sum = 0 " << endl ;
    ofOut << "Report_LE_sum = 0 " << endl ;
    ofOut << "Report_H_sum = 0 " << endl ; 
    ofOut << "# Maps of species specific energy " << endl;
    ofOut << "Report_Canopy_Temp = 0 " << endl;
    ofOut << "Report_Canopy_NetR = 0 " << endl;
    ofOut << "Report_Canopy_LE_E = 0 " << endl;
    ofOut << "Report_Canopy_LE_T = 0 " << endl;
    ofOut << "Report_Canopy_Sens_Heat = 0 " << endl << endl;

    ofOut << "#-------------------------------" << endl;
    ofOut << "#Evapotranspiration Output Maps " << endl;
    ofOut << "#-------------------------------" << endl;
    ofOut << "Report_Total_ET = 0 " << endl;
    ofOut << "Report_Transpiration_sum = 0 " << endl;
    ofOut << "Report_Einterception_sum = 0 " << endl;
    ofOut << "Report_Esoil_sum = 0 " << endl;
    ofOut << "Report_ChannelE = 0 " << endl;
    ofOut << "# Maps of species specific evapotranspiration " << endl;
    ofOut << "Report_species_ET = 0 " << endl;
    ofOut << "Report_Transpiration = 0 " << endl;
    ofOut << "Report_Einterception = 0 " << endl;
    ofOut << "Report_Esoil = 0 " << endl << endl;

    ofOut << "#-------------------------------" << endl;
    ofOut << "#Vegetation Species Input Maps " << endl;
    ofOut << "#-------------------------------" << endl;
    ofOut << "Report_Veget_frac = 0 " << endl;
    ofOut << "Report_Stem_Density = 0 " << endl;
    ofOut << "Report_RootFracL1_species = 0 " << endl;
    ofOut << "Report_RootFracL2_species = 0 " << endl << endl;

    ofOut << "#-------------------------------" << endl;
    ofOut << "#Vegetation Species Output Maps " << endl;
    ofOut << "#-------------------------------" << endl; 
    ofOut << "Report_Leaf_Area_Index = 0 " << endl;
    ofOut << "Report_Stand_Age = 0 " << endl;
    ofOut << "Report_Canopy_Conductance = 0 " << endl;
    ofOut << "Report_GPP = 0 " << endl;
    ofOut << "Report_NPP = 0 " << endl;
    ofOut << "Report_Basal_Area = 0 " << endl;
    ofOut << "Report_Tree_Height = 0 " << endl;
    ofOut << "Report_Root_Mass = 0 " << endl;
    ofOut << "Report_Canopy_Water_Stor = 0 " << endl;
    ofOut << "Report_RootUptake_Prop_L1 = 0 " << endl;
    ofOut << "Report_RootUptake_Prop_L2 = 0 " << endl;
    ofOut << "Report_RootUptake_Prop_L3 = 0 " << endl;        
    ofOut << "Report_Soil_Water_Potential = 0 " << endl;
    ofOut << "Report_Leaf_Water_Potential = 0 " << endl;
    ofOut << "Report_Sapflow_Velocity = 0 " << endl << endl;

    ofOut << "#-------------------------------" << endl;
    ofOut << "#Water Balance (Flux) Output Maps " << endl;
    ofOut << "#-------------------------------" << endl; 
    ofOut << "#Maps of outflow from cell/layer " << endl;
    ofOut << "Report_GW_to_Channel = 0 " << endl;
    ofOut << "Report_DeepGW_to_Channel = 0 " << endl;
    ofOut << "Report_Surface_to_Channel = 0 " << endl;
    ofOut << "Report_Overland_Outflow = 0" << endl ;
    ofOut << "Report_Stream_Outflow = 0" << endl;
    ofOut << "Report_Groundwater_Outflow = 0" << endl ;
    ofOut << "Report_DeepGW_Outflow = 0 " << endl;
    ofOut << "Report_Return_Flow_Surface = 0" << endl ;
    ofOut << "Report_Percolation_to_Layer2 = 0" << endl ;
    ofOut << "Report_Percolation_to_Layer3 = 0" << endl ;
    ofOut << "Report_Groundwater_Recharge = 0" << endl ;
    ofOut << "Report_Bedrock_Leakance = 0" << endl;
    ofOut << "#Maps of inflow from cell/layer " << endl;
    ofOut << "Report_Overland_Inflow = 0" << endl ;
    ofOut << "Report_Stream_Inflow = 0" << endl;
    ofOut << "Report_Groundwater_Inflow = 0 " << endl ;
    ofOut << "Report_DeepGW_Inflow = 0 " << endl;
    ofOut << "Report_Infiltration = 0" << endl ;
    ofOut << "Report_Return_Flow_to_Layer1 = 0" << endl ;
    ofOut << "Report_Return_Flow_to_Layer2 = 0" << endl << endl;

    ofOut << "#-------------------------------" << endl;
    ofOut << "#Accumulation Output Maps " << endl;
    ofOut << "#-------------------------------" << endl;    
    ofOut << "Report_Infiltration_acc = 0" << endl ;
    ofOut << "Report_Return_Flow_Surface_acc = 0" << endl ;
    ofOut << "Report_Percolation_to_Layer2_acc = 0" << endl ;
    ofOut << "Report_Return_Flow_to_Layer1_acc = 0" << endl ;
    ofOut << "Report_Percolation_to_Layer3_acc = 0" << endl ;
    ofOut << "Report_Groundwater_Recharge_acc = 0" << endl ;
    ofOut << "Report_Return_Flow_to_Layer2_acc = 0" << endl ;
    ofOut << "Report_Soil_Evaporation_acc = 0" << endl ;
    ofOut << "Report_Transpiration_Layer1_acc = 0" << endl ;
    ofOut << "Report_Transpiration_Layer2_acc = 0" << endl ;
    ofOut << "Report_Transpiration_Layer3_acc = 0" << endl ;
    ofOut << "Report_Overland_Inflow_acc = 0" << endl ;
    ofOut << "Report_Stream_Inflow_acc = 0" << endl ;
    ofOut << "Report_Groundwater_Inflow_acc = 0" << endl ;
    ofOut << "Report_DeepGW_Inflow_acc = 0 " << endl;
    ofOut << "Report_Overland_Outflow_acc = 0" << endl ;
    ofOut << "Report_Stream_Outflow_acc = 0" << endl ;
    ofOut << "Report_Groundwater_Outflow_acc = 0" << endl << endl;
    ofOut << "Report_DeepGW_Outflow_acc = 0 " << endl; 
    ofOut << "Report_GW_to_Channel_acc = 0 " << endl;
    ofOut << "Report_DeepGW_to_Channel_acc = 0 " << endl;
    ofOut << "Report_Surface_to_Channel_acc = 0 " << endl;
    
    ofOut << "#***************************************" << endl;        
    ofOut << "#Report time series section " << endl;
    ofOut << "#***************************************" << endl << endl;
    ofOut << "#-------------------------------" << endl;
    ofOut << "#Input Ts " << endl;
    ofOut << "#-------------------------------" << endl; 
    ofOut << "TS_mask = Tsmask.map " << endl <<"#" << endl;
    ofOut << "Ts_OutletDischarge = 1 " << endl;
    ofOut << "Ts_OutletGW = 1 " << endl;
    ofOut << "Ts_Long_Rad_Down = 0 " << endl;
    ofOut << "Ts_Short_Rad_Down = 0 " << endl;
    ofOut << "Ts_Precip = 0 " << endl;
    ofOut << "Ts_Rel_Humidity = 0 " << endl;
    ofOut << "Ts_Wind_Speed = 0 " << endl;
    ofOut << "Ts_AvgAir_Temperature = 0 " << endl;
    ofOut << "Ts_MinAir_Temperature = 0 " << endl;
    ofOut << "Ts_MaxAir_Temperature = 0 " << endl;
    ofOut << "Ts_Anthropogenic_Heat = 0 " << endl << endl;

    ofOut << "#-------------------------------" << endl;
    ofOut << "#Water Balance (Storage) Output Ts " << endl;
    ofOut << "#-------------------------------" << endl;  
    ofOut << "Ts_SWE = 0 " << endl;
    ofOut << "Ts_Infilt_Cap = 0 " << endl;
    ofOut << "Ts_Streamflow = 0 " << endl;
    ofOut << "Ts_Ponding = 0 " << endl;
    ofOut << "Ts_Soil_Water_Content_Average = 0 " << endl;
    ofOut << "Ts_Soil_Water_Content_Up = 0 " << endl;
    ofOut << "Ts_Soil_Water_Content_L1 = 1 " << endl;
    ofOut << "Ts_Soil_Water_Content_L2 = 1 " << endl;
    ofOut << "Ts_Soil_Water_Content_L3 = 1 " << endl;
    ofOut << "Ts_WaterTableDepth = 0 " << endl;
    ofOut << "Ts_Soil_Sat_Deficit = 0 " << endl;
    ofOut << "Ts_Ground_Water = 0 " << endl;
    ofOut << "Ts_Deep_Ground_Water = 0 " << endl;
    ofOut << "Ts_Canopy_Water_Stor_sum = 0 " << endl;
    ofOut << "# TS of time-constant variables --" << endl; 
    ofOut << "Ts_Field_Capacity_L1 = 0 " << endl;
    ofOut << "Ts_Field_Capacity_L2 = 0 " << endl;
    ofOut << "Ts_Field_Capacity_L3 = 0 " << endl << endl;
    
    ofOut << "#-------------------------------" << endl;
    ofOut << "#Energy Balance Output Ts " << endl;
    ofOut << "#-------------------------------" << endl;
    ofOut << "Ts_Soil_Net_Rad = 0 " << endl;
    ofOut << "Ts_Soil_LE = 0 " << endl;
    ofOut << "Ts_Sens_Heat = 0 " << endl;
    ofOut << "Ts_Grnd_Heat = 0 " << endl;
    ofOut << "Ts_Snow_Heat = 0 " << endl;
    ofOut << "Ts_Soil_Temperature = 0 " << endl;
    ofOut << "Ts_SoilL1_Temperature = 0 " << endl;
    ofOut << "Ts_SoilL2_Temperature = 0 " << endl;
    ofOut << "Ts_SoilL3_Temperature = 0 " << endl; 
    ofOut << "Ts_Skin_Temperature = 0 " << endl;
    ofOut << "Ts_Water_Temperature = 0 " << endl;
    ofOut << "Ts_Net_Rad_sum = 0 " << endl ;
    ofOut << "Ts_LE_sum = 0 " << endl ;
    ofOut << "Ts_H_sum = 0 " << endl ;   
    ofOut << "# Ts of species specific energy " << endl;
    ofOut << "Ts_Canopy_Temp = 0 " << endl;
    ofOut << "Ts_Canopy_NetR = 0 " << endl;
    ofOut << "Ts_Canopy_LE_E = 0 " << endl;
    ofOut << "Ts_Canopy_LE_T = 0 " << endl;
    ofOut << "Ts_Canopy_Sens_Heat = 0 " << endl << endl;

    ofOut << "#-------------------------------" << endl;
    ofOut << "#Evapotranspiration Output Ts " << endl;
    ofOut << "#-------------------------------" << endl;    
    ofOut << "Ts_Total_ET = 1 " << endl;
    ofOut << "Ts_Transpiration_sum = 1 " << endl;
    ofOut << "Ts_Einterception_sum = 1 " << endl;
    ofOut << "Ts_Esoil_sum = 1 " << endl;
    ofOut << "Ts_ChannelE = 0 " << endl;
    ofOut << "# Ts of species specific evapotranspiration " << endl;
    ofOut << "Ts_species_ET = 0 " << endl;
    ofOut << "Ts_Transpiration = 0 " << endl;
    ofOut << "Ts_Einterception = 0 " << endl;
    ofOut << "Ts_Esoil = 0 " << endl << endl;

    ofOut << "#-------------------------------" << endl;
    ofOut << "#Vegetation Species Input Ts " << endl;
    ofOut << "#-------------------------------" << endl;  
    ofOut << "Ts_Veget_frac = 0 " << endl;
    ofOut << "Ts_Stem_Density = 0 " << endl;
    ofOut << "Ts_RootFracL1_species = 0 " << endl;
    ofOut << "Ts_RootFracL2_species = 0 " << endl << endl;

    ofOut << "#-------------------------------" << endl;
    ofOut << "#Vegetation Species Output Ts " << endl;
    ofOut << "#-------------------------------" << endl; 
    ofOut << "Ts_Leaf_Area_Index = 0 " << endl;
    ofOut << "Ts_Stand_Age = 0 " << endl;
    ofOut << "Ts_Canopy_Conductance = 0 " << endl;
    ofOut << "Ts_GPP = 0 " << endl;
    ofOut << "Ts_NPP = 0 " << endl;
    ofOut << "Ts_Basal_Area = 0 " << endl;
    ofOut << "Ts_Tree_Height = 0 " << endl;
    ofOut << "Ts_Root_Mass = 0 " << endl;
    ofOut << "Ts_Canopy_Water_Stor = 0 " << endl;
    ofOut << "Ts_RootUptake_Prop_L1 = 0 "<< endl;
    ofOut << "Ts_RootUptake_Prop_L2 = 0 "<< endl;
    ofOut << "Ts_RootUptake_Prop_L3 = 0 "<< endl;
    ofOut << "Ts_Soil_Water_Potential = 0 " << endl;
    ofOut << "Ts_Leaf_Water_Potential = 0 " << endl;
    ofOut << "Ts_Sapflow_Velocity = 0 " << endl << endl;

    ofOut << "#-------------------------------" << endl;
    ofOut << "#Water Balance (Flux) Output Ts " << endl;
    ofOut << "#-------------------------------" << endl;    
    ofOut << "#Ts of outflow from cell/layer " << endl; 
    ofOut << "Ts_GW_to_Channel = 0 " << endl;
    ofOut << "Ts_DeepGW_to_Channel = 0 " << endl;
    ofOut << "Ts_Surface_to_Channel = 0 " << endl;
    ofOut << "Ts_Overland_Outflow = 0" << endl ;
    ofOut << "Ts_Stream_Outflow = 0" << endl;
    ofOut << "Ts_Groundwater_Outflow = 0" << endl ;
    ofOut << "Ts_DeepGW_Outflow = 0 " << endl;
    ofOut << "Ts_Percolation_to_Layer2 = 0" << endl ;
    ofOut << "Ts_Percolation_to_Layer3 = 0" << endl ;
    ofOut << "Ts_Groundwater_Recharge = 0" << endl ;
    ofOut << "Ts_Bedrock_Leakance = 0" << endl;
    ofOut << "#Ts of inflow from cell/layer " << endl;
    ofOut << "Ts_Overland_Inflow = 0" << endl ;
    ofOut << "Ts_Stream_Inflow = 0" << endl;
    ofOut << "Ts_Groundwater_Inflow = 0 " << endl ;
    ofOut << "Ts_DeepGW_Inflow = 0 " << endl;
    ofOut << "Ts_Infiltration = 0" << endl ;
    ofOut << "Ts_Return_Flow_Surface = 0" << endl ;
    ofOut << "Ts_Return_Flow_to_Layer1 = 0" << endl ;
    ofOut << "Ts_Return_Flow_to_Layer2 = 0" << endl ;

    if (ofOut)
      ofOut.close();
  }
  catch(const std::exception &e){
    cout << "Failure writing configuration template file with  " << e.what() << endl;
    exit(EXIT_FAILURE);
  }
}

