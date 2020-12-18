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
 * ReadConfigTrck.cpp
 *
 *  Created on: Feb 28, 2017
 *      Author: Sylvain Kuppel
 *
 *      Reads the configuration file for tracking (default name 'configTrck.ini')
 *      and places the file names and control data in the Control object
 */
#include <stdexcept>   // for exception, runtime_error, out_of_range
#include <iostream>
#include <sys/stat.h>
#include <stdlib.h>
#include "ConfigFile.h"
#include "InitConf.h"
using namespace std;

int Control::ReadConfigTrck(string confilename /*= "configTrck.ini"*/)
{

  try{
    ConfigFile ConfigTrck(confilename);
    
    sw_2H = ConfigTrck.read<bool>("water_2H");
    sw_18O = ConfigTrck.read<bool>("water_18O");
    sw_Age = ConfigTrck.read<bool>("water_Age");
    sw_frac = ConfigTrck.read<bool>("water_frac");
    sw_TPD = ConfigTrck.read<bool>("water_two-pore_domain");
    
    if(sw_frac){
      toggle_hs = ConfigTrck.read<bool>("Fractionation_surface_relhum");
      toggle_n = ConfigTrck.read<bool>("Fractionation_turbulent_factor");
      toggle_ek = ConfigTrck.read<bool>("Fractionation_kinetic_diffusion");
      sw_chan_frac = ConfigTrck.read<bool>("channel_water_frac");
    }

    if(sw_TPD)
      ConfigTrck.readInto(fn_psi_MW, "MobileWater_Transition_Head");
    //psiMW = ConfigTrck.read<float>("MobileWater_Transition_Head");
		
    toggle_mix = ConfigTrck.read<int>("Mixing_mode");

    if(toggle_mix == 3)
      ConfigTrck.readInto(fn_IncompMix, "Incomplete_Mixing");


    if(sw_TPD){
      Rep_Moist_MW1 = ConfigTrck.read<bool>("Rep_Moisture_MobileWater_L1");
      Rep_Moist_MW2 = ConfigTrck.read<bool>("Rep_Moisture_MobileWater_L2");
      Rep_Frac_MW1 = ConfigTrck.read<bool>("Rep_Frac_MobileWater_L1");
      Rep_Frac_MW2 = ConfigTrck.read<bool>("Rep_Frac_MobileWater_L2");
      Rep_Frac_MW12 = ConfigTrck.read<bool>("Rep_Frac_MobileWater_Up");
      //
      RepTs_Moist_MW1 = ConfigTrck.read<bool>("Ts_Moisture_MobileWater_L1");
      RepTs_Moist_MW2 = ConfigTrck.read<bool>("Ts_Moisture_MobileWater_L2");
      RepTs_Frac_MW1 = ConfigTrck.read<bool>("Ts_Frac_MobileWater_L1");
      RepTs_Frac_MW2 = ConfigTrck.read<bool>("Ts_Frac_MobileWater_L2");
      RepTs_Frac_MW12 = ConfigTrck.read<bool>("Ts_Frac_MobileWater_Up");
      //
    }

    if(sw_2H){
      
      // Precipitation input
      ConfigTrck.readInto(fn_d2Hprecip, "d2H_precip");
      // Initial values
      //ConfigTrck.readInto(fn_d2Hcanopy, "init_d2H_canopy");
      ConfigTrck.readInto(fn_d2Hsnowpack, "init_d2H_snowpack");
      ConfigTrck.readInto(fn_d2Hsurface, "init_d2H_surface");
      ConfigTrck.readInto(fn_d2Hsoil1, "init_d2H_soil1");
      ConfigTrck.readInto(fn_d2Hsoil2, "init_d2H_soil2");
      ConfigTrck.readInto(fn_d2Hsoil3, "init_d2H_soil3");
      ConfigTrck.readInto(fn_d2Hgroundwater, "init_d2H_groundwater");
      //ConfigTrck.readInto(fn_d2Hchannel, "init_d2H_channel");

      // Report maps
      Rep_d2Hprecip = ConfigTrck.read<bool>("Rep_d2Hprecip");
      Rep_d2Hcanopy = ConfigTrck.read<bool>("Rep_d2Hcanopy");
      Rep_d2Hcanopy_sum = ConfigTrck.read<bool>("Rep_d2Hcanopy_sum");
      Rep_d2Hsnowpack = ConfigTrck.read<bool>("Rep_d2Hsnowpack");
      Rep_d2Hsurface = ConfigTrck.read<bool>("Rep_d2Hsurface");
      Rep_d2Hchan = ConfigTrck.read<bool>("Rep_d2Hchan");
      Rep_d2Hsoil1 = ConfigTrck.read<bool>("Rep_d2Hsoil1");
      Rep_d2Hsoil2 = ConfigTrck.read<bool>("Rep_d2Hsoil2");
      Rep_d2HsoilUp = ConfigTrck.read<bool>("Rep_d2HsoilUp");
      Rep_d2Hsoil3 = ConfigTrck.read<bool>("Rep_d2Hsoil3");
      Rep_d2HsoilAv = ConfigTrck.read<bool>("Rep_d2HsoilAv");
      Rep_d2Hgroundwater = ConfigTrck.read<bool>("Rep_d2Hgroundwater");
      Rep_d2Hleakage = ConfigTrck.read<bool>("Rep_d2Hleakage");
      Rep_d2HevapS = ConfigTrck.read<bool>("Rep_d2HevapS");
      Rep_d2HevapS_sum = ConfigTrck.read<bool>("Rep_d2HevapS_sum");
      Rep_d2HevapI = ConfigTrck.read<bool>("Rep_d2HevapI");
      Rep_d2HevapI_sum = ConfigTrck.read<bool>("Rep_d2HevapI_sum");
      Rep_d2HevapT = ConfigTrck.read<bool>("Rep_d2HevapT");
      Rep_d2HevapT_sum = ConfigTrck.read<bool>("Rep_d2HevapT_sum");
      Rep_d2H_MW1 = ConfigTrck.read<bool>("Rep_d2Hsoil1_MobileWater");
      Rep_d2H_MW2 = ConfigTrck.read<bool>("Rep_d2Hsoil2_MobileWater");
      Rep_d2H_TB1 = ConfigTrck.read<bool>("Rep_d2Hsoil1_TightlyBound");
      Rep_d2H_TB2 = ConfigTrck.read<bool>("Rep_d2Hsoil2_TightlyBound");
      // Report time series
      RepTs_d2Hprecip = ConfigTrck.read<bool>("Ts_d2Hprecip");
      RepTs_d2Hcanopy = ConfigTrck.read<bool>("Ts_d2Hcanopy");
      RepTs_d2Hcanopy_sum = ConfigTrck.read<bool>("Ts_d2Hcanopy_sum");
      RepTs_d2Hsnowpack = ConfigTrck.read<bool>("Ts_d2Hsnowpack");
      RepTs_d2Hsurface = ConfigTrck.read<bool>("Ts_d2Hsurface");
      RepTs_d2Hchan = ConfigTrck.read<bool>("Ts_d2Hchan");
      RepTs_d2Hsoil1 = ConfigTrck.read<bool>("Ts_d2Hsoil1");
      RepTs_d2Hsoil2 = ConfigTrck.read<bool>("Ts_d2Hsoil2");
      RepTs_d2HsoilUp = ConfigTrck.read<bool>("Ts_d2HsoilUp");
      RepTs_d2Hsoil3 = ConfigTrck.read<bool>("Ts_d2Hsoil3");
      RepTs_d2HsoilAv = ConfigTrck.read<bool>("Ts_d2HsoilAv");
      RepTs_d2Hgroundwater = ConfigTrck.read<bool>("Ts_d2Hgroundwater");
      RepTs_d2Hleakage = ConfigTrck.read<bool>("Ts_d2Hleakage");
      RepTs_d2HevapS = ConfigTrck.read<bool>("Ts_d2HevapS");
      RepTs_d2HevapS_sum = ConfigTrck.read<bool>("Ts_d2HevapS_sum");
      RepTs_d2HevapI = ConfigTrck.read<bool>("Ts_d2HevapI");
      RepTs_d2HevapI_sum = ConfigTrck.read<bool>("Ts_d2HevapI_sum");
      RepTs_d2HevapT = ConfigTrck.read<bool>("Ts_d2HevapT");
      RepTs_d2HevapT_sum = ConfigTrck.read<bool>("Ts_d2HevapT_sum");
      RepTs_d2H_MW1 = ConfigTrck.read<bool>("Ts_d2Hsoil1_MobileWater");
      RepTs_d2H_MW2 = ConfigTrck.read<bool>("Ts_d2Hsoil2_MobileWater");
      RepTs_d2H_TB1 = ConfigTrck.read<bool>("Ts_d2Hsoil1_TightlyBound");
      RepTs_d2H_TB2 = ConfigTrck.read<bool>("Ts_d2Hsoil2_TightlyBound");
    }
    
    if(sw_18O){
      // Precipitation input
      ConfigTrck.readInto(fn_d18Oprecip, "d18O_precip");
      // Initial values
      //ConfigTrck.readInto(fn_d18Ocanopy, "init_d18O_canopy");
      ConfigTrck.readInto(fn_d18Osnowpack, "init_d18O_snowpack");
      ConfigTrck.readInto(fn_d18Osurface, "init_d18O_surface");
      ConfigTrck.readInto(fn_d18Osoil1, "init_d18O_soil1");
      ConfigTrck.readInto(fn_d18Osoil2, "init_d18O_soil2");
      ConfigTrck.readInto(fn_d18Osoil3, "init_d18O_soil3");
      ConfigTrck.readInto(fn_d18Ogroundwater, "init_d18O_groundwater");
      //ConfigTrck.readInto(fn_d18Ochannel, "init_d18O_channel");
      // Report maps
      Rep_d18Oprecip = ConfigTrck.read<bool>("Rep_d18Oprecip");
      Rep_d18Ocanopy = ConfigTrck.read<bool>("Rep_d18Ocanopy");
      Rep_d18Ocanopy_sum = ConfigTrck.read<bool>("Rep_d18Ocanopy_sum");
      Rep_d18Osnowpack = ConfigTrck.read<bool>("Rep_d18Osnowpack");
      Rep_d18Osurface = ConfigTrck.read<bool>("Rep_d18Osurface");
      Rep_d18Ochan = ConfigTrck.read<bool>("Rep_d18Ochan");
      Rep_d18Osoil1 = ConfigTrck.read<bool>("Rep_d18Osoil1");
      Rep_d18Osoil2 = ConfigTrck.read<bool>("Rep_d18Osoil2");
      Rep_d18OsoilUp = ConfigTrck.read<bool>("Rep_d18OsoilUp");
      Rep_d18Osoil3 = ConfigTrck.read<bool>("Rep_d18Osoil3");
      Rep_d18OsoilAv = ConfigTrck.read<bool>("Rep_d18OsoilAv");
      Rep_d18Ogroundwater = ConfigTrck.read<bool>("Rep_d18Ogroundwater");
      Rep_d18Oleakage = ConfigTrck.read<bool>("Rep_d18Oleakage");
      Rep_d18OevapS = ConfigTrck.read<bool>("Rep_d18OevapS");
      Rep_d18OevapS_sum = ConfigTrck.read<bool>("Rep_d18OevapS_sum");
      Rep_d18OevapI = ConfigTrck.read<bool>("Rep_d18OevapI");
      Rep_d18OevapI_sum = ConfigTrck.read<bool>("Rep_d18OevapI_sum");
      Rep_d18OevapT = ConfigTrck.read<bool>("Rep_d18OevapT");
      Rep_d18OevapT_sum = ConfigTrck.read<bool>("Rep_d18OevapT_sum");
      Rep_d18O_MW1 = ConfigTrck.read<bool>("Rep_d18Osoil1_MobileWater");
      Rep_d18O_MW2 = ConfigTrck.read<bool>("Rep_d18Osoil2_MobileWater");
      Rep_d18O_TB1 = ConfigTrck.read<bool>("Rep_d18Osoil1_TightlyBound");
      Rep_d18O_TB2 = ConfigTrck.read<bool>("Rep_d18Osoil2_TightlyBound");
      // Report time series
      RepTs_d18Oprecip = ConfigTrck.read<bool>("Ts_d18Oprecip");
      RepTs_d18Ocanopy = ConfigTrck.read<bool>("Ts_d18Ocanopy");
      RepTs_d18Ocanopy_sum = ConfigTrck.read<bool>("Ts_d18Ocanopy_sum");
      RepTs_d18Osnowpack = ConfigTrck.read<bool>("Ts_d18Osnowpack");
      RepTs_d18Osurface = ConfigTrck.read<bool>("Ts_d18Osurface");
      RepTs_d18Ochan = ConfigTrck.read<bool>("Ts_d18Ochan");
      RepTs_d18Osoil1 = ConfigTrck.read<bool>("Ts_d18Osoil1");
      RepTs_d18Osoil2 = ConfigTrck.read<bool>("Ts_d18Osoil2");
      RepTs_d18OsoilUp = ConfigTrck.read<bool>("Ts_d18OsoilUp");
      RepTs_d18Osoil3 = ConfigTrck.read<bool>("Ts_d18Osoil3");
      RepTs_d18OsoilAv = ConfigTrck.read<bool>("Ts_d18OsoilAv");
      RepTs_d18Ogroundwater = ConfigTrck.read<bool>("Ts_d18Ogroundwater");
      RepTs_d18Oleakage = ConfigTrck.read<bool>("Ts_d18Oleakage");
      RepTs_d18OevapS = ConfigTrck.read<bool>("Ts_d18OevapS");
      RepTs_d18OevapS_sum = ConfigTrck.read<bool>("Ts_d18OevapS_sum");
      RepTs_d18OevapI = ConfigTrck.read<bool>("Ts_d18OevapI");
      RepTs_d18OevapI_sum = ConfigTrck.read<bool>("Ts_d18OevapI_sum");
      RepTs_d18OevapT = ConfigTrck.read<bool>("Ts_d18OevapT");
      RepTs_d18OevapT_sum = ConfigTrck.read<bool>("Ts_d18OevapT_sum");
      RepTs_d18O_MW1 = ConfigTrck.read<bool>("Ts_d18Osoil1_MobileWater");
      RepTs_d18O_MW2 = ConfigTrck.read<bool>("Ts_d18Osoil2_MobileWater");
      RepTs_d18O_TB1 = ConfigTrck.read<bool>("Ts_d18Osoil1_TightlyBound");
      RepTs_d18O_TB2 = ConfigTrck.read<bool>("Ts_d18Osoil2_TightlyBound");
    }

    if(sw_Age){

      // Initial values
      //ConfigTrck.readInto(fn_Agecanopy, "init_Age_canopy");
      ConfigTrck.readInto(fn_Agesnowpack, "init_Age_snowpack");
      ConfigTrck.readInto(fn_Agesurface, "init_Age_surface");
      ConfigTrck.readInto(fn_Agesoil1, "init_Age_soil1");
      ConfigTrck.readInto(fn_Agesoil2, "init_Age_soil2");
      ConfigTrck.readInto(fn_Agesoil3, "init_Age_soil3");
      ConfigTrck.readInto(fn_Agegroundwater, "init_Age_groundwater");
      //ConfigTrck.readInto(fn_Agechannel, "init_Age_channel");
      // Report maps
      Rep_Agecanopy = ConfigTrck.read<bool>("Rep_Agecanopy");
      Rep_Agecanopy_sum = ConfigTrck.read<bool>("Rep_Agecanopy_sum");
      Rep_Agesnowpack = ConfigTrck.read<bool>("Rep_Agesnowpack");
      Rep_Agesurface = ConfigTrck.read<bool>("Rep_Agesurface");
      Rep_Agechan = ConfigTrck.read<bool>("Rep_Agechan");
      Rep_Agesoil1 = ConfigTrck.read<bool>("Rep_Agesoil1");
      Rep_Agesoil2 = ConfigTrck.read<bool>("Rep_Agesoil2");
      Rep_AgesoilUp = ConfigTrck.read<bool>("Rep_AgesoilUp");
      Rep_Agesoil3 = ConfigTrck.read<bool>("Rep_Agesoil3");
      Rep_AgesoilAv = ConfigTrck.read<bool>("Rep_AgesoilAv");
      Rep_Agegroundwater = ConfigTrck.read<bool>("Rep_Agegroundwater");
      Rep_Ageleakage = ConfigTrck.read<bool>("Rep_Ageleakage");
      Rep_AgeevapS = ConfigTrck.read<bool>("Rep_AgeevapS");
      Rep_AgeevapS_sum = ConfigTrck.read<bool>("Rep_AgeevapS_sum");
      Rep_AgeevapI = ConfigTrck.read<bool>("Rep_AgeevapI");
      Rep_AgeevapI_sum = ConfigTrck.read<bool>("Rep_AgeevapI_sum");
      Rep_AgeevapT = ConfigTrck.read<bool>("Rep_AgeevapT");
      Rep_AgeevapT_sum = ConfigTrck.read<bool>("Rep_AgeevapT_sum");
      Rep_AgeGWtoChn = ConfigTrck.read<bool>("Rep_AgeGWtoChn");
      Rep_AgeSrftoChn = ConfigTrck.read<bool>("Rep_AgeSrftoChn");
      Rep_AgeRecharge = ConfigTrck.read<bool>("Rep_AgeRecharge");
      Rep_Age_MW1 = ConfigTrck.read<bool>("Rep_Agesoil1_MobileWater");
      Rep_Age_MW2 = ConfigTrck.read<bool>("Rep_Agesoil2_MobileWater");
      Rep_Age_MWUp = ConfigTrck.read<bool>("Rep_AgesoilUp_MobileWater");
      //Rep_Age_MWav = ConfigTrck.read<bool>("Rep_AgesoilAv_MobileWater");
      Rep_Age_TB1 = ConfigTrck.read<bool>("Rep_Agesoil1_TightlyBound");
      Rep_Age_TB2 = ConfigTrck.read<bool>("Rep_Agesoil2_TightlyBound");
      Rep_Age_TBUp = ConfigTrck.read<bool>("Rep_AgesoilUp_TightlyBound");
      //Rep_Age_TBav = ConfigTrck.read<bool>("Rep_AgesoilAv_TightlyBound");
      // Report time series
      RepTs_Agecanopy = ConfigTrck.read<bool>("Ts_Agecanopy");
      RepTs_Agecanopy_sum = ConfigTrck.read<bool>("Ts_Agecanopy_sum");
      RepTs_Agesnowpack = ConfigTrck.read<bool>("Ts_Agesnowpack");
      RepTs_Agesurface = ConfigTrck.read<bool>("Ts_Agesurface");
      RepTs_Agechan = ConfigTrck.read<bool>("Ts_Agechan");
      RepTs_Agesoil1 = ConfigTrck.read<bool>("Ts_Agesoil1");
      RepTs_Agesoil2 = ConfigTrck.read<bool>("Ts_Agesoil2");
      RepTs_AgesoilUp = ConfigTrck.read<bool>("Ts_AgesoilUp");
      RepTs_Agesoil3 = ConfigTrck.read<bool>("Ts_Agesoil3");
      RepTs_AgesoilAv = ConfigTrck.read<bool>("Ts_AgesoilAv");
      RepTs_Agegroundwater = ConfigTrck.read<bool>("Ts_Agegroundwater");
      RepTs_Ageleakage = ConfigTrck.read<bool>("Ts_Ageleakage");
      RepTs_AgeevapS = ConfigTrck.read<bool>("Ts_AgeevapS");
      RepTs_AgeevapS_sum = ConfigTrck.read<bool>("Ts_AgeevapS_sum");
      RepTs_AgeevapI = ConfigTrck.read<bool>("Ts_AgeevapI");
      RepTs_AgeevapI_sum = ConfigTrck.read<bool>("Ts_AgeevapI_sum");
      RepTs_AgeevapT = ConfigTrck.read<bool>("Ts_AgeevapT");
      RepTs_AgeevapT_sum = ConfigTrck.read<bool>("Ts_AgeevapT_sum");
      RepTs_AgeGWtoChn = ConfigTrck.read<bool>("Ts_AgeGWtoChn");
      RepTs_AgeSrftoChn = ConfigTrck.read<bool>("Ts_AgeSrftoChn");
      RepTs_AgeRecharge = ConfigTrck.read<bool>("Ts_AgeRecharge");
      RepTs_Age_MW1 = ConfigTrck.read<bool>("Ts_Agesoil1_MobileWater");
      RepTs_Age_MW2 = ConfigTrck.read<bool>("Ts_Agesoil2_MobileWater");
      RepTs_Age_MWUp = ConfigTrck.read<bool>("Ts_AgesoilUp_MobileWater");
      RepTs_Age_TB1 = ConfigTrck.read<bool>("Ts_Agesoil1_TightlyBound");
      RepTs_Age_TB2 = ConfigTrck.read<bool>("Ts_Agesoil2_TightlyBound");
      RepTs_Age_TBUp = ConfigTrck.read<bool>("Ts_AgesoilUp_TightlyBound");
    }
  }
  catch(ConfigFile::file_not_found &fn){
    cout << "Tracking file " << fn.filename << " not found\n";
    exit(EXIT_SUCCESS);
  }
  catch(ConfigFile::key_not_found &fn){
    cout << "Tracking key " << fn.key << " not found\n";
    exit(EXIT_SUCCESS);
  }
  catch(std::exception &e){
    cout << e.what();
    exit(EXIT_SUCCESS);
  }



  return 1;
}
