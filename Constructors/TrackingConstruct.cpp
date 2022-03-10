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

Tracking::Tracking(Control &ctrl, Basin &bsn, Atmosphere &atm)
{

  // reset the errno value
  errno = 0;
  //-----------------------------------------------------------------------------
  // If boundary conditions are to be set - read them in
  //-----------------------------------------------------------------------------
  _d2HsurfaceBC     = NULL;
  _d2Hlayer1BC      = NULL;
  _d2Hlayer2BC      = NULL;
  _d2HgroundwaterBC = NULL;
  _d18OsurfaceBC     = NULL;
  _d18Olayer1BC      = NULL;
  _d18Olayer2BC      = NULL;
  _d18OgroundwaterBC = NULL;
  _AgesurfaceBC     = NULL;
  _Agelayer1BC      = NULL;
  _Agelayer2BC      = NULL;
  _AgegroundwaterBC = NULL;
  if(ctrl.sw_BC){
    if(ctrl.sw_2H){
      _d2HsurfaceBC     = new grid(*bsn.getDEM());
      _d2Hlayer1BC      = new grid(*bsn.getDEM());
      _d2Hlayer2BC      = new grid(*bsn.getDEM());
      _d2HgroundwaterBC = new grid(*bsn.getDEM());

      *_d2HsurfaceBC = *bsn.getDEM();
      *_d2Hlayer1BC = *bsn.getDEM();
      *_d2Hlayer2BC = *bsn.getDEM();
      *_d2HgroundwaterBC = *bsn.getDEM();

      try{
	ifd2HsurfaceBC.open((ctrl.path_ClimMapsFolder + ctrl.fn_d2HsurfaceBC).c_str(), ios::binary);
	if(errno!=0) throw ctrl.fn_d2HsurfaceBC;
	ifd2Hlayer1BC.open((ctrl.path_ClimMapsFolder + ctrl.fn_d2Hlayer1BC).c_str(), ios::binary);
	if(errno!=0) throw ctrl.fn_d2Hlayer1BC;
	ifd2Hlayer2BC.open((ctrl.path_ClimMapsFolder + ctrl.fn_d2Hlayer2BC).c_str(), ios::binary);
	if(errno!=0) throw ctrl.fn_d2Hlayer2BC;
	ifd2HgroundwaterBC.open((ctrl.path_ClimMapsFolder + ctrl.fn_d2HgroundwaterBC).c_str(), ios::binary);
	if(errno!=0) throw ctrl.fn_d2HgroundwaterBC;
      } catch (string e){
	cout << "Dang!!: cannot find/read the " << e << " file: error " << strerror(errno) << endl;
	throw;
      }
      // initiate the boundary condition map
      try{
	if(InitiateBCMap_Iso(ifd2HsurfaceBC, *_d2HsurfaceBC, atm) != atm.getSsortedGridTotalCellNumber())
	  throw string("d2H surface BC");
	if(InitiateBCMap_Iso(ifd2Hlayer1BC, *_d2Hlayer1BC, atm) != atm.getSsortedGridTotalCellNumber())
	  throw string("d2H layer 1 BC");
	if(InitiateBCMap_Iso(ifd2Hlayer2BC, *_d2Hlayer2BC, atm) != atm.getSsortedGridTotalCellNumber())
	  throw string("d2H layer 2 BC");
	if(InitiateBCMap_Iso(ifd2HgroundwaterBC, *_d2HgroundwaterBC, atm) != atm.getSsortedGridTotalCellNumber())
	  throw string("d2H groundwater BC");
      } catch (string e) {
	cout << "Error: some sections of the domain were not filled with " << e << "data. " << endl;
	cout << "Please verify that all the boundary zones in the map are presented in the binary boundary data file " << endl;
	cout << "and that the n boundary zones present are the first n zones in the binary boundary data file" << endl;
	throw;
      }
    } //end deuterium

    if(ctrl.sw_18O){
      _d18OsurfaceBC     = new grid(*bsn.getDEM());
      _d18Olayer1BC      = new grid(*bsn.getDEM());
      _d18Olayer2BC      = new grid(*bsn.getDEM());
      _d18OgroundwaterBC = new grid(*bsn.getDEM());

      *_d18OsurfaceBC = *bsn.getDEM();
      *_d18Olayer1BC = *bsn.getDEM();
      *_d18Olayer2BC = *bsn.getDEM();
      *_d18OgroundwaterBC = *bsn.getDEM();

      try{
	ifd18OsurfaceBC.open((ctrl.path_ClimMapsFolder + ctrl.fn_d18OsurfaceBC).c_str(), ios::binary);
	if(errno!=0) throw ctrl.fn_d18OsurfaceBC;
	ifd18Olayer1BC.open((ctrl.path_ClimMapsFolder + ctrl.fn_d18Olayer1BC).c_str(), ios::binary);
	if(errno!=0) throw ctrl.fn_d18Olayer1BC;
	ifd18Olayer2BC.open((ctrl.path_ClimMapsFolder + ctrl.fn_d18Olayer2BC).c_str(), ios::binary);
	if(errno!=0) throw ctrl.fn_d18Olayer2BC;
	ifd18OgroundwaterBC.open((ctrl.path_ClimMapsFolder + ctrl.fn_d18OgroundwaterBC).c_str(), ios::binary);
	if(errno!=0) throw ctrl.fn_d18OgroundwaterBC;
      } catch (string e){
	cout << "Dang!!: cannot find/read the " << e << " file: error " << strerror(errno) << endl;
	throw;
      }
      // initiate the boundary condition map
      try{
	if(InitiateBCMap_Iso(ifd18OsurfaceBC, *_d18OsurfaceBC, atm) != atm.getSsortedGridTotalCellNumber())
	  throw string("d18O surface BC");
	if(InitiateBCMap_Iso(ifd18Olayer1BC, *_d18Olayer1BC, atm) != atm.getSsortedGridTotalCellNumber())
	  throw string("d18O layer 1 BC");
	if(InitiateBCMap_Iso(ifd18Olayer2BC, *_d18Olayer2BC, atm) != atm.getSsortedGridTotalCellNumber())
	  throw string("d18O layer 2 BC");
	if(InitiateBCMap_Iso(ifd18OgroundwaterBC, *_d18OgroundwaterBC, atm) != atm.getSsortedGridTotalCellNumber())
	  throw string("d18O groundwater BC");
      } catch (string e) {
	cout << "Error: some sections of the domain were not filled with " << e << "data. " << endl;
	cout << "Please verify that all the boundary zones in the map are presented in the binary boundary data file " << endl;
	cout << "and that the n boundary zones present are the first n zones in the binary boundary data file" << endl;
	throw;
      }
    } // end oxygen-18

    if(ctrl.sw_Age){
      _AgesurfaceBC     = new grid(*bsn.getDEM());
      _Agelayer1BC      = new grid(*bsn.getDEM());
      _Agelayer2BC      = new grid(*bsn.getDEM());
      _AgegroundwaterBC = new grid(*bsn.getDEM());

      *_AgesurfaceBC = *bsn.getDEM();
      *_Agelayer1BC = *bsn.getDEM();
      *_Agelayer2BC = *bsn.getDEM();
      *_AgegroundwaterBC = *bsn.getDEM();

      try{
	ifAgesurfaceBC.open((ctrl.path_ClimMapsFolder + ctrl.fn_AgesurfaceBC).c_str(), ios::binary);
	if(errno!=0) throw ctrl.fn_AgesurfaceBC;
	ifAgelayer1BC.open((ctrl.path_ClimMapsFolder + ctrl.fn_Agelayer1BC).c_str(), ios::binary);
	if(errno!=0) throw ctrl.fn_Agelayer1BC;
	ifAgelayer2BC.open((ctrl.path_ClimMapsFolder + ctrl.fn_Agelayer2BC).c_str(), ios::binary);
	if(errno!=0) throw ctrl.fn_Agelayer2BC;
	ifAgegroundwaterBC.open((ctrl.path_ClimMapsFolder + ctrl.fn_AgegroundwaterBC).c_str(), ios::binary);
	if(errno!=0) throw ctrl.fn_AgegroundwaterBC;
      } catch (string e){
	cout << "Dang!!: cannot find/read the " << e << " file: error " << strerror(errno) << endl;
	throw;
      }
      // initiate the boundary condition map
      try{
	if(InitiateBCMap_Iso(ifd2HsurfaceBC, *_AgesurfaceBC, atm) != atm.getSsortedGridTotalCellNumber())
	  throw string("Age surface BC");
	if(InitiateBCMap_Iso(ifAgelayer1BC, *_Agelayer1BC, atm) != atm.getSsortedGridTotalCellNumber())
	  throw string("Age layer 1 BC");
	if(InitiateBCMap_Iso(ifAgelayer2BC, *_Agelayer2BC, atm) != atm.getSsortedGridTotalCellNumber())
	  throw string("Age layer 2 BC");
	if(InitiateBCMap_Iso(ifAgegroundwaterBC, *_AgegroundwaterBC, atm) != atm.getSsortedGridTotalCellNumber())
	  throw string("Age groundwater BC");
      } catch (string e) {
	cout << "Error: some sections of the domain were not filled with " << e << "data. " << endl;
	cout << "Please verify that all the boundary zones in the map are presented in the binary boundary data file " << endl;
	cout << "and that the n boundary zones present are the first n zones in the binary boundary data file" << endl;
	throw;
      }
    }
  }
  
  // Construct NULL pointer in case the object is not fully constructed
  // (avoids memory leak)
  _d2Hcanopy_sum = NULL;                              //Total canopy d2H
  _d2Hthroughfall_sum = NULL;                         //Total throughfall d2H
  _d2Hsnowpack = NULL;                                //Total snowpack d2H
  _d2Hsnowmelt = NULL;                                //Total snowmelt d2H
  _d2Hsurface = NULL;                                 //Total surface water d2H
  _d2Hchan = NULL;                                    //Total channel d2H
  _d2Hsoil1 = NULL;                                   //Total soil L1 d2H
  _d2Hsoil2 = NULL;                                   //Total soil L2 d2H
  _d2Hsoil_12 = NULL;                                 //Total soil L1+L2 d2H
  _d2Hsoil3 = NULL;                                   //Total soil L3 d2H
  _d2HsoilAv = NULL;                                  //Total soil 10cm avg d2H
  _d2HdeepGW = NULL;                                  //Total deep GW storage
  _d2Hsoil1Q = NULL;                               //Total soil L1 outflow d2H
  _d2Hsoil2Q = NULL;                               //Total soil L2 outflow d2H
  _d2Hgroundwater = NULL;                             //Total GW outflow d2H
  _d2HdeepGWQ = NULL;                             //Total Deep GW outflow d2H
  _d2HevapS_sum = NULL;                               //Total soil evap d2H
  _d2HevapC_sum = NULL;                               //Total channel evap d2H
  _d2HevapI_sum = NULL;                               //Total interception d2H
  _d2HevapT_sum = NULL;                               //Total transpiration d2H
  _d2HevapI_Vap_sum = NULL;                           //Total interception vapour d2H
  _d2HevapT_Vap_sum = NULL;                           //Total transpiration vapour d2H
  _d2HGWtoChn = NULL;                                 //Total GW to channel d2H
  _d2HDeepGWtoChn = NULL;                             //Total Deep GW to channel d2H
  _d2HSrftoChn = NULL;                                //Total surface to channel d2H
  _d2HRecharge = NULL;                                //Total recharge d2H
  _d2Hleakage = NULL;                                 //Total leakage d2H
  _Fd2HLattoSrf = NULL;                               //Total lateral inflow to surface d2H
  _Fd2HLattoChn = NULL;                               //Total lateral inflow to channel d2H
  _Fd2HLattoGW = NULL;                                //Total lateral inflow to GW d2H
  _Fd2HLattoDeepGW = NULL;                            //Total lateral inflow to DeepGW d2H
  _d2H_MW1 = NULL;                                    //Total mobile water L1 d2H
  _d2H_MW2 = NULL;                                    //Total mobile water L2 d2H
  _d2H_TB1 = NULL;                                    //Total bulk water L1 d2H
  _d2H_TB2 = NULL;                                    //Total bulk water L2 d2H

  _d18Ocanopy_sum = NULL;
  _d18Othroughfall_sum = NULL;  
  _d18Osnowpack = NULL;
  _d18Osnowmelt = NULL;
  _d18Osurface = NULL;
  _d18Ochan= NULL;
  _d18Osoil1 = NULL;
  _d18Osoil2 = NULL;
  _d18Osoil_12 = NULL;
  _d18Osoil3 = NULL;
  _d18OdeepGW = NULL;
  _d18OsoilAv = NULL;
  _d18Osoil1Q = NULL;
  _d18Osoil2Q = NULL;
  _d18Ogroundwater = NULL;
  _d18OdeepGWQ = NULL;                           //Total Deep GW outflow d18O
  _d18OevapS_sum = NULL;
  _d18OevapC_sum = NULL;
  _d18OevapI_sum = NULL;
  _d18OevapT_sum = NULL;
  _d18OevapI_Vap_sum = NULL;
  _d18OevapT_Vap_sum = NULL;  
  _d18OGWtoChn = NULL;
  _d18ODeepGWtoChn = NULL;                           //Total Deep GW to channel d18O
  _d18OSrftoChn = NULL;
  _d18ORecharge = NULL;
  _d18Oleakage = NULL;
  _Fd18OLattoSrf = NULL;
  _Fd18OLattoChn = NULL;
  _Fd18OLattoGW = NULL;
  _Fd18OLattoDeepGW = NULL;                          //Total lateral inflow to DeepGW d18O
  _d18O_MW1 = NULL;
  _d18O_MW2 = NULL;
  _d18O_TB1 = NULL;
  _d18O_TB2 = NULL;

  _Agecanopy_sum = NULL;
  _Agethroughfall_sum = NULL;  
  _Agesnowpack = NULL;
  _Agesnowmelt = NULL;
  _Agesurface = NULL;
  _Agechan = NULL;
  _Agesoil1 = NULL;
  _Agesoil2 = NULL;
  _Agesoil_12 = NULL;
  _Agesoil3 = NULL;
  _AgesoilAv = NULL;
  _AgedeepGW = NULL;
  _Agesoil1Q = NULL;
  _Agesoil2Q = NULL;
  _Agegroundwater = NULL;
  _AgedeepGWQ = NULL;                            //Total Deep GW outflow Age
  _AgeevapS_sum = NULL;
  _AgeevapC_sum = NULL;
  _AgeevapI_sum = NULL;
  _AgeevapT_sum = NULL;
  _AgeGWtoChn = NULL;
  _AgeDeepGWtoChn = NULL;                            //Total Deep GW to channel Age
  _AgeSrftoChn = NULL;
  _AgeRecharge = NULL;
  _Ageleakage = NULL;
  _FAgeLattoSrf = NULL;
  _FAgeLattoChn = NULL;
  _FAgeLattoGW = NULL;
  _FAgeLattoDeepGW = NULL;                           //Total lateral inflow to DeepGW Age

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
      _d2Hthroughfall_sum = new grid(*bsn.getDEM());      
      _d2Hsnowpack = new grid(ctrl.path_BasinFolder+ctrl.fn_d2Hsnowpack, ctrl.MapType);
      _d2Hsnowmelt = new grid(*bsn.getDEM());
      _d2Hsurface = new grid(ctrl.path_BasinFolder + ctrl.fn_d2Hsurface, ctrl.MapType);
      _d2Hchan = new grid(ctrl.path_BasinFolder + ctrl.fn_d2Hsurface, ctrl.MapType);
      _d2Hsoil1 = new grid(ctrl.path_BasinFolder + ctrl.fn_d2Hsoil1, ctrl.MapType);
      _d2Hsoil2 = new grid(ctrl.path_BasinFolder + ctrl.fn_d2Hsoil2, ctrl.MapType);
      _d2Hsoil_12 = new grid(*bsn.getDEM());
      _d2Hsoil3 = new grid(ctrl.path_BasinFolder + ctrl.fn_d2Hsoil3, ctrl.MapType);
      _d2HsoilAv = new grid(*bsn.getDEM());
      _d2Hsoil1Q = new grid(*bsn.getDEM());
      _d2Hsoil2Q = new grid(*bsn.getDEM());
      _d2Hgroundwater = new grid(ctrl.path_BasinFolder + ctrl.fn_d2Hgroundwater, ctrl.MapType);
      _d2HevapS_sum = new grid(*bsn.getDEM());
      _d2HevapC_sum = new grid(*bsn.getDEM());
      _d2HevapI_sum = new grid(*bsn.getDEM());
      _d2HevapT_sum = new grid(*bsn.getDEM());
      _d2HevapI_Vap_sum = new grid(*bsn.getDEM());
      _d2HevapT_Vap_sum = new grid(*bsn.getDEM());            
      _d2HGWtoChn = new grid(*bsn.getDEM());
      _d2HSrftoChn = new grid(*bsn.getDEM());
      _d2HRecharge = new grid(*bsn.getDEM());
      _d2Hleakage = new grid(*bsn.getDEM());
      _Fd2HLattoSrf = new grid(*bsn.getDEM());
      _Fd2HLattoChn = new grid(*bsn.getDEM());
      _Fd2HLattoGW = new grid(*bsn.getDEM());

      if(ctrl.sw_TPD){
	_d2H_MW1 = new grid(*bsn.getDEM());
	_d2H_MW2 = new grid(*bsn.getDEM());
	_d2H_TB1 = new grid(*bsn.getDEM());
	_d2H_TB2 = new grid(*bsn.getDEM());
      }
      if(ctrl.sw_deepGW){
	_d2HdeepGW = new grid(ctrl.path_BasinFolder + ctrl.fn_d2HDeepGW, ctrl.MapType);
	_d2HdeepGWQ = new grid(ctrl.path_BasinFolder + ctrl.fn_d2HDeepGW, ctrl.MapType);
	_Fd2HLattoDeepGW = new grid(*bsn.getDEM());
	_d2HDeepGWtoChn = new grid(*bsn.getDEM());
      }      
    }
	  
    if(ctrl.sw_18O){
      //_d18Oinput = new grid(*bsn.getDEM(), -60);
      _d18Ocanopy_sum = new grid(*bsn.getDEM());
      _d18Othroughfall_sum = new grid(*bsn.getDEM());            
      _d18Osnowpack = new grid(ctrl.path_BasinFolder + ctrl.fn_d18Osnowpack, ctrl.MapType);
      _d18Osnowmelt = new grid(*bsn.getDEM());
      _d18Osurface = new grid(ctrl.path_BasinFolder + ctrl.fn_d18Osurface, ctrl.MapType);
      _d18Ochan = new grid(ctrl.path_BasinFolder + ctrl.fn_d18Osurface, ctrl.MapType);
      _d18Osoil1 = new grid(ctrl.path_BasinFolder + ctrl.fn_d18Osoil1, ctrl.MapType);
      _d18Osoil2 = new grid(ctrl.path_BasinFolder + ctrl.fn_d18Osoil2, ctrl.MapType);
      _d18Osoil_12 = new grid(*bsn.getDEM());
      _d18Osoil3 = new grid(ctrl.path_BasinFolder + ctrl.fn_d18Osoil3, ctrl.MapType);
      _d18OsoilAv = new grid(*bsn.getDEM());
      _d18Osoil1Q = new grid(*bsn.getDEM());
      _d18Osoil2Q = new grid(*bsn.getDEM());
      _d18Ogroundwater = new grid(ctrl.path_BasinFolder + ctrl.fn_d18Ogroundwater, ctrl.MapType);
      _d18OevapS_sum = new grid(*bsn.getDEM());
      _d18OevapC_sum = new grid(*bsn.getDEM());
      _d18OevapI_sum = new grid(*bsn.getDEM());
      _d18OevapT_sum = new grid(*bsn.getDEM());
      _d18OevapI_Vap_sum = new grid(*bsn.getDEM());
      _d18OevapT_Vap_sum = new grid(*bsn.getDEM());            
      _d18OGWtoChn = new grid(*bsn.getDEM());
      _d18OSrftoChn = new grid(*bsn.getDEM());
      _d18ORecharge = new grid(*bsn.getDEM());
      _d18Oleakage = new grid(*bsn.getDEM());	    
      _Fd18OLattoSrf = new grid(*bsn.getDEM());
      _Fd18OLattoChn = new grid(*bsn.getDEM());
      _Fd18OLattoGW = new grid(*bsn.getDEM());

      if(ctrl.sw_TPD){
	_d18O_MW1 = new grid(*bsn.getDEM());	    
	_d18O_MW2 = new grid(*bsn.getDEM());	    
	_d18O_TB1 = new grid(*bsn.getDEM());	    
	_d18O_TB2 = new grid(*bsn.getDEM());	    
      }
      if(ctrl.sw_deepGW){
	_d18OdeepGW = new grid(ctrl.path_BasinFolder + ctrl.fn_d18ODeepGW, ctrl.MapType);
	_d18OdeepGWQ = new grid(ctrl.path_BasinFolder + ctrl.fn_d18ODeepGW, ctrl.MapType);
	_Fd18OLattoDeepGW = new grid(*bsn.getDEM());
	_d18ODeepGWtoChn = new grid(*bsn.getDEM());
      }      
    }
	  
    if(ctrl.sw_Age){
      //_Ageinput = new grid(*bsn.getDEM(), -60);
      _Agecanopy_sum = new grid(*bsn.getDEM());
      _Agethroughfall_sum = new grid(*bsn.getDEM());            
      _Agesnowpack = new grid(ctrl.path_BasinFolder + ctrl.fn_Agesnowpack, ctrl.MapType);
      _Agesnowmelt = new grid(*bsn.getDEM());
      _Agesurface = new grid(ctrl.path_BasinFolder + ctrl.fn_Agesurface, ctrl.MapType);
      _Agechan = new grid(ctrl.path_BasinFolder + ctrl.fn_Agesurface, ctrl.MapType);
      _Agesoil1 = new grid(ctrl.path_BasinFolder + ctrl.fn_Agesoil1, ctrl.MapType);
      _Agesoil2 = new grid(ctrl.path_BasinFolder + ctrl.fn_Agesoil2, ctrl.MapType);
      _Agesoil_12 = new grid(*bsn.getDEM());
      _Agesoil3 = new grid(ctrl.path_BasinFolder + ctrl.fn_Agesoil3, ctrl.MapType);
      _AgesoilAv = new grid(*bsn.getDEM());
      _Agesoil1Q = new grid(*bsn.getDEM());
      _Agesoil2Q = new grid(*bsn.getDEM());
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
      _FAgeLattoGW = new grid(*bsn.getDEM());
	
      if(ctrl.sw_TPD){
	_Age_MW1 = new grid(*bsn.getDEM());
	_Age_MW2 = new grid(*bsn.getDEM());
	_Age_MW12 = new grid(*bsn.getDEM());
	_Age_TB1 = new grid(*bsn.getDEM());
	_Age_TB2 = new grid(*bsn.getDEM());
	_Age_TB12 = new grid(*bsn.getDEM());
      }
      if(ctrl.sw_deepGW){
	_AgedeepGW = new grid(ctrl.path_BasinFolder + ctrl.fn_AgeDeepGW, ctrl.MapType);	
	_AgedeepGWQ = new grid(ctrl.path_BasinFolder + ctrl.fn_AgeDeepGW, ctrl.MapType);
	_FAgeLattoDeepGW = new grid(*bsn.getDEM());
	_AgeDeepGWtoChn = new grid(*bsn.getDEM());
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

      // boundary conditions
      if(_d2HsurfaceBC)
	delete _d2HsurfaceBC;
      if(_d2Hlayer1BC)
	delete _d2Hlayer1BC;
      if(_d2Hlayer2BC)
	delete _d2Hlayer2BC;
      if(_d2HgroundwaterBC)
	delete _d2HgroundwaterBC;	          
      // Ratios
      if(_d2Hcanopy_sum)
	delete _d2Hcanopy_sum;
      if(_d2Hthroughfall_sum)
	delete _d2Hthroughfall_sum;            
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
      if(_d2Hsoil1Q)
	delete _d2Hsoil1Q;
      if(_d2Hsoil2Q)
	delete _d2Hsoil2Q;
      if(_d2HdeepGWQ)
	delete _d2HdeepGWQ;
      if(_d2Hgroundwater)
	delete _d2Hgroundwater;
      if(_d2HevapS_sum)
	delete _d2HevapS_sum;
      if(_d2HevapI_sum)
	delete _d2HevapI_sum;
      if(_d2HevapT_sum)
	delete _d2HevapT_sum;
      if(_d2HevapI_Vap_sum)
	delete _d2HevapI_Vap_sum;
      if(_d2HevapT_Vap_sum)
	delete _d2HevapT_Vap_sum;      
      if(_d2HGWtoChn)
	delete _d2HGWtoChn;
      if(_d2HDeepGWtoChn)
	delete _d2HDeepGWtoChn;      
      if(_d2HSrftoChn)
	delete _d2HSrftoChn;
      if(_d2HRecharge)
	delete _d2HRecharge;
      if(_d2Hleakage)
	delete _d2Hleakage;
      if(_Fd2HLattoGW)
	delete _Fd2HLattoGW;
      if(_Fd2HLattoDeepGW)
	delete _Fd2HLattoDeepGW;
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

      // boundary conditions
      if(_d18OsurfaceBC)
	delete _d18OsurfaceBC;
      if(_d18Olayer1BC)
	delete _d18Olayer1BC;
      if(_d18Olayer2BC)
	delete _d18Olayer2BC;
      if(_d18OgroundwaterBC)
	delete _d18OgroundwaterBC;
      //Ratios
      if(_d18Ocanopy_sum)
	delete _d18Ocanopy_sum;
      if(_d18Othroughfall_sum)
	delete _d18Othroughfall_sum;                  
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
      if(_d18OdeepGW)
	delete _d18OdeepGW;      
      if(_d18Osoil1Q)
	delete _d18Osoil1Q;
      if(_d18Osoil2Q)
	delete _d18Osoil2Q;
      if(_d18Ogroundwater)
	delete _d18Ogroundwater;
      if(_d18OdeepGWQ)
	delete _d18OdeepGWQ;
      if(_d18OevapS_sum)
	delete _d18OevapS_sum;
      if(_d18OevapI_sum)
	delete _d18OevapI_sum;
      if(_d18OevapT_sum)
	delete _d18OevapT_sum;
      if(_d18OevapI_Vap_sum)
	delete _d18OevapI_Vap_sum;
      if(_d18OevapT_Vap_sum)
	delete _d18OevapT_Vap_sum;      
      if(_d18OGWtoChn)
	delete _d18OGWtoChn;
      if(_d18ODeepGWtoChn)
	delete _d18ODeepGWtoChn;      
      if(_d18OSrftoChn)
	delete _d18OSrftoChn;
      if(_d18ORecharge)
	delete _d18ORecharge;
      if(_d18Oleakage)
	delete _d18Oleakage;
      if(_Fd18OLattoGW)
	delete _Fd18OLattoGW;
      if(_Fd18OLattoDeepGW)
	delete _Fd18OLattoDeepGW;
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

      // boundary conditions
      if(_AgesurfaceBC)
	delete _AgesurfaceBC;
      if(_Agelayer1BC)
	delete _Agelayer1BC;
      if(_Agelayer2BC)
	delete _Agelayer2BC;
      if(_AgegroundwaterBC)
	delete _AgegroundwaterBC;	    
      // Ratio	    
      if(_Agecanopy_sum)
	delete _Agecanopy_sum;
      if(_Agethroughfall_sum)
	delete _Agethroughfall_sum;                  
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
      if(_AgedeepGW)
	delete _AgedeepGW;
      if(_Agesoil1Q)
	delete _Agesoil1Q;
      if(_Agesoil2Q)
	delete _Agesoil2Q;
      if(_Agegroundwater)
	delete _Agegroundwater;
      if(_AgedeepGWQ)
	delete _AgedeepGWQ;
      if(_AgeevapS_sum)
	delete _AgeevapS_sum;
      if(_AgeevapI_sum)
	delete _AgeevapI_sum;
      if(_AgeevapT_sum)
	delete _AgeevapT_sum;
      if(_AgeGWtoChn)
	delete _AgeGWtoChn;      
      if(_AgeDeepGWtoChn)
	delete _AgeDeepGWtoChn;
      if(_AgeSrftoChn)
	delete _AgeSrftoChn;
      if(_AgeRecharge)
	delete _AgeRecharge;
      if(_Ageleakage)
	delete _Ageleakage;
      if(_FAgeLattoGW)
	delete _FAgeLattoGW;
      if(_FAgeLattoDeepGW)
	delete _FAgeLattoDeepGW;
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
