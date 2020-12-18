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
 * IsotopeConstruct.cpp
 *
 *  Created on: Nov 14, 2016
 *      Author: Sylvain Kuppel
 */

#include "Tracking.h"
#include <errno.h>
#include <stdio.h>
#include <string.h>

Tracking::Tracking(Control &ctrl, Basin &bsn)
{

  // reset the errno value
  errno = 0;
  
  // Construct NULL pointer in case the object is not fully constructed
  // (avoids memory leak)
  _d2Hcanopy_sum = NULL;
  _d2Hsnowpack = NULL;
  _d2Hsnowmelt = NULL;
  _d2Hsurface = NULL;
  _d2Hchan = NULL;
  _d2Hsoil1 = NULL;
  _d2Hsoil2 = NULL;
  _d2Hsoil_12 = NULL;
  _d2Hsoil3 = NULL;
  _d2HsoilAv = NULL;
  _d2Hsoil1_Out = NULL;
  _d2Hsoil2_Out = NULL;
  _d2Hgroundwater = NULL;
  _d2HevapS_sum = NULL;
  _d2HevapC_sum = NULL;
  _d2HevapI_sum = NULL;
  _d2HevapT_sum = NULL;
  _d2HGWtoChn = NULL;
  _d2HSrftoChn = NULL;
  _d2HRecharge = NULL;
  _d2Hleakage = NULL;
  _Fd2HLattoSrf = NULL;
  _Fd2HLattoChn = NULL;
  _Fd2HLattoL1 = NULL;
  _Fd2HLattoL2 = NULL;
  _Fd2HLattoGW = NULL;
  _d2H_MW1 = NULL;
  _d2H_MW2 = NULL;
  _d2H_TB1 = NULL;
  _d2H_TB2 = NULL;

  _d18Ocanopy_sum = NULL;
  _d18Osnowpack = NULL;
  _d18Osnowmelt = NULL;
  _d18Osurface = NULL;
  _d18Ochan= NULL;
  _d18Osoil1 = NULL;
  _d18Osoil2 = NULL;
  _d18Osoil_12 = NULL;
  _d18Osoil3 = NULL;
  _d18OsoilAv = NULL;
  _d18Osoil1_Out = NULL;
  _d18Osoil2_Out = NULL;
  _d18Ogroundwater = NULL;
  _d18OevapS_sum = NULL;
  _d18OevapC_sum = NULL;
  _d18OevapI_sum = NULL;
  _d18OevapT_sum = NULL;
  _d18OGWtoChn = NULL;
  _d18OSrftoChn = NULL;
  _d18ORecharge = NULL;
  _d18Oleakage = NULL;
  _Fd18OLattoSrf = NULL;
  _Fd18OLattoChn = NULL;
  _Fd18OLattoL1 = NULL;
  _Fd18OLattoL2 = NULL;
  _Fd18OLattoGW = NULL;
  _d18O_MW1 = NULL;
  _d18O_MW2 = NULL;
  _d18O_TB1 = NULL;
  _d18O_TB2 = NULL;

  _Agecanopy_sum = NULL;
  _Agesnowpack = NULL;
  _Agesnowmelt = NULL;
  _Agesurface = NULL;
  _Agechan = NULL;
  _Agesoil1 = NULL;
  _Agesoil2 = NULL;
  _Agesoil_12 = NULL;
  _Agesoil3 = NULL;
  _AgesoilAv = NULL;
  _Agesoil1_Out = NULL;
  _Agesoil2_Out = NULL;
  _Agegroundwater = NULL;
  _AgeevapS_sum = NULL;
  _AgeevapC_sum = NULL;
  _AgeevapI_sum = NULL;
  _AgeevapT_sum = NULL;
  _AgeGWtoChn = NULL;
  _AgeSrftoChn = NULL;
  _AgeRecharge = NULL;
  _Ageleakage = NULL;
  _FAgeLattoSrf = NULL;
  _FAgeLattoChn = NULL;
  _FAgeLattoL1 = NULL;
  _FAgeLattoL2 = NULL;
  _FAgeLattoGW = NULL;
  _Age_MW1 = NULL;
  _Age_MW2 = NULL;
  _Age_MW12 = NULL;
  _Age_TB1 = NULL;
  _Age_TB2 = NULL;
  _Age_TB12 = NULL;

  if(!ctrl.sw_trck){
    ctrl.sw_2H = 0;
    ctrl.sw_18O = 0;
    ctrl.sw_Age = 0;
  }

  try{
    if(ctrl.sw_2H){
      /*state variables initialized with the base map*/
      _d2Hcanopy_sum = new grid(*bsn.getDEM());
      _d2Hsnowpack = new grid(ctrl.path_BasinFolder+ctrl.fn_d2Hsnowpack, ctrl.MapType);
      _d2Hsnowmelt = new grid(*bsn.getDEM());
      _d2Hsurface = new grid(ctrl.path_BasinFolder + ctrl.fn_d2Hsurface, ctrl.MapType);
      _d2Hchan = new grid(ctrl.path_BasinFolder + ctrl.fn_d2Hsurface, ctrl.MapType);
      _d2Hsoil1 = new grid(ctrl.path_BasinFolder + ctrl.fn_d2Hsoil1, ctrl.MapType);
      _d2Hsoil2 = new grid(ctrl.path_BasinFolder + ctrl.fn_d2Hsoil2, ctrl.MapType);
      _d2Hsoil_12 = new grid(*bsn.getDEM());
      _d2Hsoil3 = new grid(ctrl.path_BasinFolder + ctrl.fn_d2Hsoil3, ctrl.MapType);
      _d2HsoilAv = new grid(*bsn.getDEM());
      _d2Hsoil1_Out = new grid(*bsn.getDEM());
      _d2Hsoil2_Out = new grid(*bsn.getDEM());
      _d2Hgroundwater = new grid(ctrl.path_BasinFolder + ctrl.fn_d2Hgroundwater, ctrl.MapType);
      _d2HevapS_sum = new grid(*bsn.getDEM());
      _d2HevapC_sum = new grid(*bsn.getDEM());
      _d2HevapI_sum = new grid(*bsn.getDEM());
      _d2HevapT_sum = new grid(*bsn.getDEM());
      _d2HGWtoChn = new grid(*bsn.getDEM());
      _d2HSrftoChn = new grid(*bsn.getDEM());
      _d2HRecharge = new grid(*bsn.getDEM());
      _d2Hleakage = new grid(*bsn.getDEM());
      _Fd2HLattoSrf = new grid(*bsn.getDEM());
      _Fd2HLattoChn = new grid(*bsn.getDEM());
      _Fd2HLattoL1 = new grid(*bsn.getDEM());
      _Fd2HLattoL2 = new grid(*bsn.getDEM());
      _Fd2HLattoGW = new grid(*bsn.getDEM());

      if(ctrl.sw_TPD){
	_d2H_MW1 = new grid(*bsn.getDEM());
	_d2H_MW2 = new grid(*bsn.getDEM());
	_d2H_TB1 = new grid(*bsn.getDEM());
	_d2H_TB2 = new grid(*bsn.getDEM());
      }
    }
	  
    if(ctrl.sw_18O){
      //_d18Oinput = new grid(*bsn.getDEM(), -60);
      _d18Ocanopy_sum = new grid(*bsn.getDEM());
      _d18Osnowpack = new grid(ctrl.path_BasinFolder + ctrl.fn_d18Osnowpack, ctrl.MapType);
      _d18Osnowmelt = new grid(*bsn.getDEM());
      _d18Osurface = new grid(ctrl.path_BasinFolder + ctrl.fn_d18Osurface, ctrl.MapType);
      _d18Ochan = new grid(ctrl.path_BasinFolder + ctrl.fn_d18Osurface, ctrl.MapType);
      _d18Osoil1 = new grid(ctrl.path_BasinFolder + ctrl.fn_d18Osoil1, ctrl.MapType);
      _d18Osoil2 = new grid(ctrl.path_BasinFolder + ctrl.fn_d18Osoil2, ctrl.MapType);
      _d18Osoil_12 = new grid(*bsn.getDEM());
      _d18Osoil3 = new grid(ctrl.path_BasinFolder + ctrl.fn_d18Osoil3, ctrl.MapType);
      _d18OsoilAv = new grid(*bsn.getDEM());
      _d18Osoil1_Out = new grid(*bsn.getDEM());
      _d18Osoil2_Out = new grid(*bsn.getDEM());
      _d18Ogroundwater = new grid(ctrl.path_BasinFolder + ctrl.fn_d18Ogroundwater, ctrl.MapType);
      _d18OevapS_sum = new grid(*bsn.getDEM());
      _d18OevapC_sum = new grid(*bsn.getDEM());
      _d18OevapI_sum = new grid(*bsn.getDEM());
      _d18OevapT_sum = new grid(*bsn.getDEM());
      _d18OGWtoChn = new grid(*bsn.getDEM());
      _d18OSrftoChn = new grid(*bsn.getDEM());
      _d18ORecharge = new grid(*bsn.getDEM());
      _d18Oleakage = new grid(*bsn.getDEM());	    
      _Fd18OLattoSrf = new grid(*bsn.getDEM());
      _Fd18OLattoChn = new grid(*bsn.getDEM());
      _Fd18OLattoL1 = new grid(*bsn.getDEM());
      _Fd18OLattoL2 = new grid(*bsn.getDEM());
      _Fd18OLattoGW = new grid(*bsn.getDEM());

      if(ctrl.sw_TPD){
	_d18O_MW1 = new grid(*bsn.getDEM());	    
	_d18O_MW2 = new grid(*bsn.getDEM());	    
	_d18O_TB1 = new grid(*bsn.getDEM());	    
	_d18O_TB2 = new grid(*bsn.getDEM());	    
      }
    }
	  
    if(ctrl.sw_Age){
      //_Ageinput = new grid(*bsn.getDEM(), -60);
      _Agecanopy_sum = new grid(*bsn.getDEM());
      _Agesnowpack = new grid(ctrl.path_BasinFolder + ctrl.fn_Agesnowpack, ctrl.MapType);
      _Agesnowmelt = new grid(*bsn.getDEM());
      _Agesurface = new grid(ctrl.path_BasinFolder + ctrl.fn_Agesurface, ctrl.MapType);
      _Agechan = new grid(ctrl.path_BasinFolder + ctrl.fn_Agesurface, ctrl.MapType);
      _Agesoil1 = new grid(ctrl.path_BasinFolder + ctrl.fn_Agesoil1, ctrl.MapType);
      _Agesoil2 = new grid(ctrl.path_BasinFolder + ctrl.fn_Agesoil2, ctrl.MapType);
      _Agesoil_12 = new grid(*bsn.getDEM());
      _Agesoil3 = new grid(ctrl.path_BasinFolder + ctrl.fn_Agesoil3, ctrl.MapType);
      _AgesoilAv = new grid(*bsn.getDEM());
      _Agesoil1_Out = new grid(*bsn.getDEM());
      _Agesoil2_Out = new grid(*bsn.getDEM());
      _Agegroundwater = new grid(ctrl.path_BasinFolder + ctrl.fn_Agegroundwater, ctrl.MapType);
      _AgeevapS_sum = new grid(*bsn.getDEM());
      _AgeevapC_sum = new grid(*bsn.getDEM());
      _AgeevapI_sum = new grid(*bsn.getDEM());
      _AgeevapT_sum = new grid(*bsn.getDEM());
      _AgeGWtoChn = new grid(*bsn.getDEM());
      _AgeSrftoChn = new grid(*bsn.getDEM());
      _AgeRecharge = new grid(*bsn.getDEM());
      _Ageleakage = new grid(*bsn.getDEM());
      _FAgeLattoSrf = new grid(*bsn.getDEM());
      _FAgeLattoChn = new grid(*bsn.getDEM());
      _FAgeLattoL1 = new grid(*bsn.getDEM());
      _FAgeLattoL2 = new grid(*bsn.getDEM());
      _FAgeLattoGW = new grid(*bsn.getDEM());
	
      if(ctrl.sw_TPD){
	_Age_MW1 = new grid(*bsn.getDEM());
	_Age_MW2 = new grid(*bsn.getDEM());
	_Age_MW12 = new grid(*bsn.getDEM());
	_Age_TB1 = new grid(*bsn.getDEM());
	_Age_TB2 = new grid(*bsn.getDEM());
	_Age_TB12 = new grid(*bsn.getDEM());
      }
    }

    try{

      //Partial check of maps mainly to make sure no no data is written within the valid domain
      CheckMapsTrck(ctrl, bsn);
      if(errno!=0){
	cout << "Error creating tracking maps: " << endl;
	throw string("  ");
      }
            
      // If Two-pore domain, initialize the maps based on inputs maps
      if(ctrl.sw_TPD){
	CalcInitTPD(bsn, ctrl);
	if(errno!=0){
	  cout << "Error calculating the initial tightly-bound / mobile tracking signatures: " << endl;
	  throw string("  ");
	}
      }
      
    } catch (string e){
      cout << "Check the  " << e << " maps, error " << strerror(errno) << endl;
      throw;
    }

	  
  }catch (std::bad_alloc &)
    { cerr << " Cleaning tracking objects..." << "\n";
	    
      // Ratios
      if(_d2Hcanopy_sum)
	delete _d2Hcanopy_sum;
      if(_d2Hsnowpack)
	delete _d2Hsnowpack;
      if(_d2Hsnowmelt)
	delete _d2Hsnowmelt;
      if(_d2Hsurface)
	delete _d2Hsurface;
      if(_d2Hchan)
	delete _d2Hchan;
      if(_d2Hsoil1)
	delete _d2Hsoil1;
      if(_d2Hsoil2)
	delete _d2Hsoil2;
      if(_d2Hsoil_12)
	delete _d2Hsoil_12;
      if(_d2Hsoil3)
	delete _d2Hsoil3;
      if(_d2HsoilAv)
	delete _d2HsoilAv;
      if(_d2Hsoil1_Out)
	delete _d2Hsoil1_Out;
      if(_d2Hsoil2_Out)
	delete _d2Hsoil2_Out;
      if(_d2Hgroundwater)
	delete _d2Hgroundwater;
      if(_d2HevapS_sum)
	delete _d2HevapS_sum;
      if(_d2HevapI_sum)
	delete _d2HevapI_sum;
      if(_d2HevapT_sum)
	delete _d2HevapT_sum;
      if(_d2HGWtoChn)
	delete _d2HGWtoChn;
      if(_d2HSrftoChn)
	delete _d2HSrftoChn;
      if(_d2HRecharge)
	delete _d2HRecharge;
      if(_d2Hleakage)
	delete _d2Hleakage;
      if(_Fd2HLattoL1)
	delete _Fd2HLattoL1;
      if(_Fd2HLattoL2)
	delete _Fd2HLattoL2;
      if(_Fd2HLattoGW)
	delete _Fd2HLattoGW;
      if(_Fd2HLattoSrf)
	delete _Fd2HLattoSrf;
      if(_Fd2HLattoChn)
	delete _Fd2HLattoChn;

      if(_d2H_MW1)
	delete _d2H_MW1;
      if(_d2H_MW2)
	delete _d2H_MW2;
      if(_d2H_TB1)
	delete _d2H_TB1;
      if(_d2H_TB2)
	delete _d2H_TB2;
	    
      if(_d18Ocanopy_sum)
	delete _d18Ocanopy_sum;
      if(_d18Osnowpack)
	delete _d18Osnowpack;
      if(_d18Osnowmelt)
	delete _d18Osnowmelt;
      if(_d18Osurface)
	delete _d18Osurface;
      if(_d18Ochan)
	delete _d18Ochan;
      if(_d18Osoil1)
	delete _d18Osoil1;
      if(_d18Osoil2)
	delete _d18Osoil2;
      if(_d18Osoil_12)
	delete _d18Osoil_12;
      if(_d18Osoil3)
	delete _d18Osoil3;
      if(_d18OsoilAv)
	delete _d18OsoilAv;
      if(_d18Osoil1_Out)
	delete _d18Osoil1_Out;
      if(_d18Osoil2_Out)
	delete _d18Osoil2_Out;
      if(_d18Ogroundwater)
	delete _d18Ogroundwater;
      if(_d18OevapS_sum)
	delete _d18OevapS_sum;
      if(_d18OevapI_sum)
	delete _d18OevapI_sum;
      if(_d18OevapT_sum)
	delete _d18OevapT_sum;
      if(_d18OGWtoChn)
	delete _d18OGWtoChn;
      if(_d18OSrftoChn)
	delete _d18OSrftoChn;
      if(_d18ORecharge)
	delete _d18ORecharge;
      if(_d18Oleakage)
	delete _d18Oleakage;
      if(_Fd18OLattoL1)
	delete _Fd18OLattoL1;
      if(_Fd18OLattoL2)
	delete _Fd18OLattoL2;
      if(_Fd18OLattoGW)
	delete _Fd18OLattoGW;
      if(_Fd18OLattoSrf)
	delete _Fd18OLattoSrf;
      if(_Fd18OLattoChn)
	delete _Fd18OLattoChn;
      if(_d18O_MW1)
	delete _d18O_MW1;
      if(_d18O_MW2)
	delete _d18O_MW2;
      if(_d18O_TB1)
	delete _d18O_TB1;
      if(_d18O_TB2)
	delete _d18O_TB2;
	    
      if(_Agecanopy_sum)
	delete _Agecanopy_sum;
      if(_Agesnowpack)
	delete _Agesnowpack;
      if(_Agesnowmelt)
	delete _Agesnowmelt;
      if(_Agesurface)
	delete _Agesurface;
      if(_Agechan)
	delete _Agechan;
      if(_Agesoil1)
	delete _Agesoil1;
      if(_Agesoil2)
	delete _Agesoil2;
      if(_Agesoil_12)
	delete _Agesoil_12;
      if(_Agesoil3)
	delete _Agesoil3;
      if(_AgesoilAv)
	delete _AgesoilAv;
      if(_Agesoil1_Out)
	delete _Agesoil1_Out;
      if(_Agesoil2_Out)
	delete _Agesoil2_Out;
      if(_Agegroundwater)
	delete _Agegroundwater;
      if(_AgeevapS_sum)
	delete _AgeevapS_sum;
      if(_AgeevapI_sum)
	delete _AgeevapI_sum;
      if(_AgeevapT_sum)
	delete _AgeevapT_sum;
      if(_AgeGWtoChn)
	delete _AgeGWtoChn;
      if(_AgeSrftoChn)
	delete _AgeSrftoChn;
      if(_AgeRecharge)
	delete _AgeRecharge;
      if(_Ageleakage)
	delete _Ageleakage;
      if(_FAgeLattoL1)
	delete _FAgeLattoL1;
      if(_FAgeLattoL2)
	delete _FAgeLattoL2;
      if(_FAgeLattoGW)
	delete _FAgeLattoGW;
      if(_FAgeLattoSrf)
	delete _FAgeLattoSrf;
      if(_FAgeLattoChn)
	delete _FAgeLattoChn;
      if(_Age_MW1)
	delete _Age_MW1;
      if(_Age_MW2)
	delete _Age_MW2;
      if(_Age_MW12)
	delete _Age_MW12;
      if(_Age_TB1)
	delete _Age_TB1;
      if(_Age_TB2)
	delete _Age_TB2;
      if(_Age_TB12)
	delete _Age_TB12;
	    	    
      throw;
    }
}
