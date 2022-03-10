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
 * GenerateConfigTemplateTrck.cpp
 *
 *  Created on: Jan 30, 2018
 *      Author: Sylvain Kuppel
 */

#include  <fstream>
#include <unistd.h>

#include "Sativa.h"

void GenerateConfigTrckTemplate(const char *fn){
  
  ofstream ofOut;
      
  try {

    ofOut.open(fn);
    if(!ofOut)
      throw std::ios::failure("Error opening file ");
    
    
    ofOut << "# == EcH2O-iso tracking configuration file v3.0 ==" << endl;
    ofOut << "# ----------------                ----------------" << endl;
    ofOut << "# further details can be found at " << endl;
    ofOut << "# http://ech2o-iso.readthedocs.io/en/latest/Keywords.html" << endl << endl;

    ofOut << "# " << endl << "# Units (inputs and outputs): " << endl;
    ofOut << "# - Deuterium ratios (2H):   permille" << endl;
    ofOut << "# - Oxygen 18 ratios (18O): permille" << endl;
    ofOut << "# - Water age:               days" << endl << endl;
    
    ofOut << "#***************************************" << endl;      
    ofOut << "#== Boolean switches" << endl;
    ofOut << "#***************************************" << endl;  
    ofOut << "water_2H = 1 # Deuterium tracking" << endl;
    ofOut << "water_18O = 0 # Oxygen-18 tracking" << endl;
    ofOut << "water_Age = 1 # Age tracking" << endl;
    ofOut << "# Isotopic fractionation from soil evaporation" << endl;
    ofOut << "water_frac = 1" << endl ;
    ofOut << "# Two-pore domain conceptualization in the soil" << endl;
    ofOut << "water_two-pore_domain = 0" << endl;
    ofOut << "# Mixing of precipitation and intercepted water"<< endl;
    ofOut << "interception_mixing = 0" << endl << endl;

    ofOut << "#***************************************" << endl;  
    ofOut << "# TOGGLE SWITCHES:" << endl << endl;
    ofOut << "#***************************************" << endl << endl; 
    ofOut << "# Approach for full mixing" << endl;
    ofOut << "# 0--> mixing computation uses Vres(t) and Cres(t+1)" << endl;
    ofOut << "# 1--> mixing computation uses (Cres(t)+Cres(t+1))/2, and " << endl;
    ofOut << "# 2--> complete mixing computation " << endl;
    ofOut << "# 3--> incomplete mixing computation (< 1 = fast flow & > 1 = slow flow) " << endl;
    ofOut << "#  useful volume = (Vres(t)+Fin+max(0,Vres(t)-Fout))/2 (still ignore storage limitations if Fin or Fout large)" << endl;
    ofOut << "Mixing_mode = 2" << endl << endl ;

    ofOut << "# Toggles switches (only used it water_frac = 1)" << endl;
    ofOut << "# Channel fractionation - option to include channel fractionation if significant" << endl;
    ofOut << "channel_water_frac = 1" << endl;
    ofOut << "# Surface relative humidity - taking into account air space between pores" << endl;
    ofOut << "# 0--> soilRH=1 " << endl; 
    ofOut << "# 1--> soilRH follows Lee and Pielke 1992 (consistent with the evaporation routine)"<< endl;
    ofOut << "# 2--> soilRH follows Sorderberg et al. (2012)" << endl;
    ofOut << "Fractionation_surface_relhum = 1" << endl; 
    ofOut << "# Turbulent factor in kinetic fractionation (n)" << endl;
    ofOut << "# 0--> n=1" <<endl;
    ofOut << "# 1--> n depends on soil water content, following Mathieu and Bariac (1996)" << endl;
    ofOut << "Fractionation_turbulent_factor = 1" << endl;
    ofOut << "# Ratio of isotope diffusivity" << endl;
    ofOut << "# 0--> Di/D = 0.9757 (2H) and 0.9727 (18O), from Merlivat (1965)" << endl;
    ofOut << "# 1--> Di/D = 0.9877 (2H) and 0.9859 (18O), from Vogt (1976)" << endl;
    ofOut << "# 2--> Empirical model by Merlivat and Jouzel (1978)" << endl;
    ofOut << "Fractionation_kinetic_diffusion = 1" << endl << endl ;

    ofOut << "#***************************************" << endl;  
    ofOut << "# -- Inputs files (only necessary if the corresponding switch is =1)" << endl;
    ofOut << "#" << endl;
    ofOut << "#***************************************" << endl << endl;  

    ofOut << "#---------------------------------------" << endl;
    ofOut << "# Climate input for isotopes" << endl;
    ofOut << "# bin files to be contained in folder pointed by Clim_Maps_Folder (see main config file)" << endl;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "d2H_precip = d2H.bin" << endl;
    ofOut << "d18O_precip = d18O.bin" << endl;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "#" << endl << "# Initial states: " << endl;
    ofOut << "# map files to be contained in folder pointed by Maps_Folder (see main config file)" << endl << "#" << endl;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "init_d2H_snowpack = d2H_snowpack.map" << endl;
    ofOut << "init_d2H_surface = d2H_surface.map" << endl;
    ofOut << "init_d2H_soil1 = d2H_soilL1.map" << endl;
    ofOut << "init_d2H_soil2 = d2H_soilL2.map" << endl;
    ofOut << "init_d2H_soil3 = d2H_soilL3.map" << endl;
    ofOut << "init_d2H_groundwater = d2H_groundwater.map" << endl;
    ofOut << "init_d2H_DeepGW = d2H_DeepGW.map" << endl << endl ;

    ofOut << "init_d18O_snowpack = d18O_snowpack.map" << endl;
    ofOut << "init_d18O_surface = d18O_surface.map" << endl;
    ofOut << "init_d18O_soil1 = d18O_soilL1.map" << endl;
    ofOut << "init_d18O_soil2 = d18O_soilL2.map" << endl;
    ofOut << "init_d18O_soil3 = d18O_soilL3.map" << endl;
    ofOut << "init_d18O_groundwater = d18O_groundwater.map" << endl;
    ofOut << "init_d18O_DeepGW = d18O_DeepGW.map" << endl << endl;

    ofOut << "init_Age_snowpack = Age_snowpack.map" << endl;
    ofOut << "init_Age_surface = Age_surface.map" << endl;
    ofOut << "init_Age_soil1 = Age_soilL1.map" << endl;
    ofOut << "init_Age_soil2 = Age_soilL2.map" << endl;
    ofOut << "init_Age_soil3 = Age_soilL3.map" << endl;
    ofOut << "init_Age_groundwater = Age_groundwater.map" << endl;
    ofOut << "init_Age_DeepGW = Age_DeepGW.map" << endl << endl;

    ofOut << "#---------------------------------------" << endl;
    ofOut << "# Two-pores domain:" << endl;
    ofOut << "# if activated, map of pressure head delimiting the two domains" << endl ;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "MobileWater_Transition_Head = TPD_transition_head.map # in meters of head" << endl << endl ;

    ofOut << "#---------------------------------------" << endl;
    ofOut << "# Incomplete mixing: " << endl;
    ofOut << "# if incomplete mixing is activate - layers 1 and 2 incomplete mixing alpha value" << endl;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "Incomplete_Mixing = IncompleteMixing_alpha.map # beta distribution alpha for mixing" << endl << endl;

    ofOut << "#---------------------------------------" << endl;
    ofOut << "# Boundary conditions (only if Boundary_Condition == 1)" << endl;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "TimeSeries_BC_d2H_surface = BCd2Hsurface" << endl;
    ofOut << "TimeSeries_BC_d2H_layer1 = BCd2Hlayer1" << endl;
    ofOut << "TimeSeries_BC_d2H_layer2 = BCd2Hlayer2" << endl;
    ofOut << "TimeSeries_BC_d2H_groundwater = BCd2Hgroundwater" << endl << endl;

    ofOut << "TimeSeries_BC_d18O_surface = BCd18Osurface" << endl;
    ofOut << "TimeSeries_BC_d18O_layer1 = BCd18Olayer1" << endl;
    ofOut << "TimeSeries_BC_d18O_layer2 = BCd18Olayer2" << endl;
    ofOut << "TimeSeries_BC_d18O_groundwater = BCd18Ogroundwater" << endl << endl;

    ofOut << "TimeSeries_BC_Age_surface = BCd2Hsurface" << endl;
    ofOut << "TimeSeries_BC_Age_layer1 = BCd2Hlayer1" << endl;
    ofOut << "TimeSeries_BC_Age_layer2 = BCd2Hlayer2" << endl;
    ofOut << "TimeSeries_BC_Age_groundwater = BCd2Hgroundwater" << endl << endl;

    ofOut << "#***************************************" << endl;  
    ofOut << "#   " << endl << "#Report map section " << endl << "#   " << endl;
    ofOut << "#***************************************" << endl;  
    ofOut << "#---------------------------------------" << endl;
    ofOut << "# If two-pore domain activated  " << endl ;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "Rep_Moisture_MobileWater_L1 = 0" << endl;
    ofOut << "Rep_Moisture_MobileWater_L2 = 0" << endl;
    ofOut << "Rep_Frac_MobileWater_L1 = 0" << endl;
    ofOut << "Rep_Frac_MobileWater_L2 = 0" << endl;
    ofOut << "Rep_Frac_MobileWater_Up = 0" << endl;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "# Storage and Fluxes" << endl;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "Rep_d2Hprecip = 0" << endl;
    ofOut << "Rep_d2Hcanopy = 0" << endl;
    ofOut << "Rep_d2Hcanopy_sum = 0" << endl;
    ofOut << "Rep_d2Hsnowpack = 0" << endl;
    ofOut << "Rep_d2Hsurface = 0" << endl;
    ofOut << "Rep_d2Hchan = 0" << endl;
    ofOut << "Rep_d2Hsoil1 = 0" << endl;
    ofOut << "Rep_d2Hsoil2 = 0" << endl;
    ofOut << "Rep_d2HsoilUp = 0" << endl;
    ofOut << "Rep_d2Hsoil3 = 0" << endl;
    ofOut << "Rep_d2HsoilAv = 0" << endl;
    ofOut << "Rep_d2Hgroundwater = 0" << endl;
    ofOut << "Rep_d2H_DeepGW_Flow = 0" << endl;
    ofOut << "Rep_d2Hleakage = 0" << endl;
    ofOut << "Rep_d2HevapS = 0" << endl;
    ofOut << "Rep_d2HevapS_sum = 0" << endl;
    ofOut << "Rep_d2HevapI = 0" << endl;
    ofOut << "Rep_d2HevapI_sum = 0" << endl;
    ofOut << "Rep_d2HevapT = 0" << endl;
    ofOut << "Rep_d2HevapT_sum = 0" << endl ;
    ofOut << "Rep_d2Hsoil1_MobileWater = 0" << endl ;
    ofOut << "Rep_d2Hsoil2_MobileWater = 0" << endl ;
    ofOut << "Rep_d2Hsoil1_TightlyBound = 0" << endl ;
    ofOut << "Rep_d2Hsoil2_TightlyBound = 0" << endl << endl ;

    ofOut << "Rep_d18Oprecip = 0" << endl;
    ofOut << "Rep_d18Ocanopy = 0" << endl;
    ofOut << "Rep_d18Ocanopy_sum = 0" << endl;
    ofOut << "Rep_d18Osnowpack = 0" << endl;
    ofOut << "Rep_d18Osurface = 0" << endl;
    ofOut << "Rep_d18Ochan = 0" << endl;
    ofOut << "Rep_d18Osoil1 = 0" << endl;
    ofOut << "Rep_d18Osoil2 = 0" << endl;
    ofOut << "Rep_d18OsoilUp = 0" << endl;
    ofOut << "Rep_d18Osoil3 = 0" << endl;
    ofOut << "Rep_d18OsoilAv = 0" << endl;
    ofOut << "Rep_d18Ogroundwater = 0" << endl;
    ofOut << "Rep_d18O_DeepGW_Flow = 0" << endl;
    ofOut << "Rep_d18Oleakage = 0" << endl;
    ofOut << "Rep_d18OevapS = 0" << endl;
    ofOut << "Rep_d18OevapS_sum = 0" << endl;
    ofOut << "Rep_d18OevapI = 0" << endl;
    ofOut << "Rep_d18OevapI_sum = 0" << endl;
    ofOut << "Rep_d18OevapT = 0" << endl;
    ofOut << "Rep_d18OevapT_sum = 0" << endl ;
    ofOut << "Rep_d18Osoil1_MobileWater = 0" << endl ;
    ofOut << "Rep_d18Osoil2_MobileWater = 0" << endl ;
    ofOut << "Rep_d18Osoil1_TightlyBound = 0" << endl ;
    ofOut << "Rep_d18Osoil2_TightlyBound = 0" << endl << endl;

    ofOut << "Rep_Agecanopy = 0" << endl;
    ofOut << "Rep_Agecanopy_sum = 0" << endl;
    ofOut << "Rep_Agesnowpack = 0" << endl;
    ofOut << "Rep_Agesurface = 0" << endl;
    ofOut << "Rep_Agechan = 0" << endl;
    ofOut << "Rep_Agesoil1 = 0" << endl;
    ofOut << "Rep_Agesoil2 = 0" << endl;
    ofOut << "Rep_AgesoilUp = 0" << endl;
    ofOut << "Rep_Agesoil3 = 0" << endl;
    ofOut << "Rep_AgesoilAv = 0" << endl;
    ofOut << "Rep_Agegroundwater = 0" << endl;
    ofOut << "Rep_Age_DeepGW_Flow = 0" << endl;
    ofOut << "Rep_Ageleakage = 0" << endl;
    ofOut << "Rep_AgeevapS = 0" << endl;
    ofOut << "Rep_AgeevapS_sum = 0" << endl;
    ofOut << "Rep_AgeevapI = 0" << endl;
    ofOut << "Rep_AgeevapI_sum = 0" << endl;
    ofOut << "Rep_AgeevapT = 0" << endl;
    ofOut << "Rep_AgeevapT_sum = 0" << endl ;
    ofOut << "Rep_AgeGWtoChn = 0" << endl ;
    ofOut << "Rep_Age_DeepGWtoChn = 0" << endl ;
    ofOut << "Rep_AgeSrftoChn = 0" << endl ;
    ofOut << "Rep_AgeRecharge = 0" << endl ;
    ofOut << "Rep_Agesoil1_MobileWater = 0" << endl ;
    ofOut << "Rep_Agesoil2_MobileWater = 0" << endl ;
    ofOut << "Rep_AgesoilUp_MobileWater = 0" << endl ;
    ofOut << "Rep_Agesoil1_TightlyBound = 0" << endl ;
    ofOut << "Rep_Agesoil2_TightlyBound = 0" << endl ;
    ofOut << "Rep_AgesoilUp_TightlyBound = 0" << endl << endl;

    ofOut << "#***************************************" << endl;  
    ofOut << "#   " << endl << "#Report time series section " << endl;
    ofOut << "#(locations specified in TS_mask map, see main config file)" << endl << "#   " << endl;
    ofOut << "#***************************************" << endl;  
    ofOut << "#---------------------------------------" << endl;
    ofOut << "# -- Report time series" << endl << endl;
    ofOut << "# If two-pore domain activated  " << endl ;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "Ts_Moisture_MobileWater_L1 = 0" << endl;
    ofOut << "Ts_Moisture_MobileWater_L2 = 0" << endl;
    ofOut << "Ts_Frac_MobileWater_L1 = 0" << endl;
    ofOut << "Ts_Frac_MobileWater_L2 = 0" << endl;
    ofOut << "Ts_Frac_MobileWater_Up = 0" << endl << endl;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "# Storage and Fluxes" << endl;
    ofOut << "#---------------------------------------" << endl;
    ofOut << "Ts_d2Hprecip = 0" << endl;
    ofOut << "Ts_d2Hcanopy = 0" << endl;
    ofOut << "Ts_d2Hcanopy_sum = 0" << endl;
    ofOut << "Ts_d2Hsnowpack = 0" << endl;
    ofOut << "Ts_d2Hsurface = 1" << endl;
    ofOut << "Ts_d2Hchan = 1" << endl;
    ofOut << "Ts_d2Hsoil1 = 1" << endl;
    ofOut << "Ts_d2Hsoil2 = 1" << endl;
    ofOut << "Ts_d2HsoilUp = 0" << endl;
    ofOut << "Ts_d2Hsoil3 = 0" << endl;
    ofOut << "Ts_d2HsoilAv = 0" << endl;
    ofOut << "Ts_d2Hgroundwater = 1" << endl;
    ofOut << "Ts_d2H_DeepGW_Flow = 0" << endl;
    ofOut << "Ts_d2Hleakage = 1" << endl;
    ofOut << "Ts_d2HevapS = 0" << endl;
    ofOut << "Ts_d2HevapS_sum = 1" << endl;
    ofOut << "Ts_d2HevapI = 0" << endl;
    ofOut << "Ts_d2HevapI_sum = 0" << endl;
    ofOut << "Ts_d2HevapT = 1" << endl;
    ofOut << "Ts_d2HevapT_sum = 0" << endl ;
    ofOut << "Ts_d2Hsoil1_MobileWater = 0" << endl ;
    ofOut << "Ts_d2Hsoil2_MobileWater = 0" << endl ;
    ofOut << "Ts_d2Hsoil1_TightlyBound = 0" << endl ;
    ofOut << "Ts_d2Hsoil2_TightlyBound = 0" << endl << endl;

    ofOut << "Ts_d18Oprecip = 0" << endl;
    ofOut << "Ts_d18Ocanopy = 0" << endl;
    ofOut << "Ts_d18Ocanopy_sum = 0" << endl;
    ofOut << "Ts_d18Osnowpack = 0 " << endl;
    ofOut << "Ts_d18Osurface = 1" << endl;
    ofOut << "Ts_d18Ochan = 1" << endl;
    ofOut << "Ts_d18Osoil1 = 1" << endl;
    ofOut << "Ts_d18Osoil2 = 1" << endl;
    ofOut << "Ts_d18OsoilUp = 0" << endl;
    ofOut << "Ts_d18Osoil3 = 0" << endl;
    ofOut << "Ts_d18OsoilAv = 0" << endl;
    ofOut << "Ts_d18Ogroundwater = 1" << endl;
    ofOut << "Ts_d18O_DeepGW_Flow = 0" << endl;
    ofOut << "Ts_d18Oleakage = 1" << endl;
    ofOut << "Ts_d18OevapS = 0" << endl;
    ofOut << "Ts_d18OevapS_sum = 1" << endl;
    ofOut << "Ts_d18OevapI = 0" << endl;
    ofOut << "Ts_d18OevapI_sum = 0" << endl;
    ofOut << "Ts_d18OevapT = 1" << endl;
    ofOut << "Ts_d18OevapT_sum = 0" << endl ;
    ofOut << "Ts_d18Osoil1_MobileWater = 0" << endl ;
    ofOut << "Ts_d18Osoil2_MobileWater = 0" << endl ;
    ofOut << "Ts_d18Osoil1_TightlyBound = 0" << endl ;
    ofOut << "Ts_d18Osoil2_TightlyBound = 0" << endl << endl;

    ofOut << "Ts_Agecanopy = 0" << endl;
    ofOut << "Ts_Agecanopy_sum = 0" << endl;
    ofOut << "Ts_Agesnowpack = 1" << endl;
    ofOut << "Ts_Agesurface = 1" << endl;
    ofOut << "Ts_Agechan = 1" << endl;
    ofOut << "Ts_Agesoil1 = 1" << endl;
    ofOut << "Ts_Agesoil2 = 1" << endl;
    ofOut << "Ts_AgesoilUp = 0" << endl;
    ofOut << "Ts_Agesoil3 = 0" << endl;
    ofOut << "Ts_AgesoilAv = 0" << endl;
    ofOut << "Ts_Agegroundwater = 1" << endl;
    ofOut << "Ts_Age_DeepGW_Flow = 0" << endl;
    ofOut << "Ts_Ageleakage = 1" << endl;
    ofOut << "Ts_AgeevapS = 0" << endl;
    ofOut << "Ts_AgeevapS_sum = 0" << endl;
    ofOut << "Ts_AgeevapI = 1" << endl;
    ofOut << "Ts_AgeevapI_sum = 0" << endl;
    ofOut << "Ts_AgeevapT = 0" << endl;
    ofOut << "Ts_AgeevapT_sum = 1" << endl ;
    ofOut << "Ts_AgeGWtoChn = 0" << endl ;
    ofOut << "Ts_Age_DeepGWtoChn = 0" << endl ;
    ofOut << "Ts_AgeSrftoChn = 0" << endl ;
    ofOut << "Ts_AgeRecharge = 0" << endl ;
    ofOut << "Ts_Agesoil1_MobileWater = 0" << endl ;
    ofOut << "Ts_Agesoil2_MobileWater = 0" << endl ;
    ofOut << "Ts_AgesoilUp_MobileWater = 0" << endl ;
    ofOut << "Ts_Agesoil1_TightlyBound = 0" << endl ;
    ofOut << "Ts_Agesoil2_TightlyBound = 0" << endl ;
    ofOut << "Ts_AgesoilUp_TightlyBound = 0" << endl << endl;


   if (ofOut)
	   ofOut.close();
}
catch(const std::exception &e){
	cout << "Failure writing configuration template file with  " << e.what() << endl;
	exit(EXIT_FAILURE);
}
}

