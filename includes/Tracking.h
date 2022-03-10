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
 *    Marco Maneta, Sylvain Kuppel, Aaron Smith
 *******************************************************************************/
/*
 * Tracking.h
 *
 *  Created on: Nov 14, 2016
 *      Author: Sylvain Kuppel
 *
 */

#ifndef TRACKING_H_
#define TRACKING_H_

#include "Grid.h"
#include "SortGrid.h"
#include "Basin.h"
#include "Forest.h"
#include "Atmosphere.h"
#include "InitConf.h"

#include <omp.h>
#include <errno.h>
#include <armadillo>

using namespace std;

class Basin;
class Forest;
class Tracking {

  Forest *fForest;

  //Spatial distribution of tracked signatures (d2H, d18O, age) in the compartments
  grid *_d2Hprecip, *_d18Oprecip;
  grid *_d2Hcanopy_sum, *_d18Ocanopy_sum, *_Agecanopy_sum; // Canopy (avg over veg fractions)
  grid *_d2Hthroughfall_sum, *_d18Othroughfall_sum, *_Agethroughfall_sum; // Canopy (avg over veg fractions)  
  grid *_d2Hsnowpack, *_d18Osnowpack, *_Agesnowpack; // Snowpack
  grid *_d2Hsnowmelt, *_d18Osnowmelt, *_Agesnowmelt; // Snowmelt
  grid *_d2Hsurface,  *_d18Osurface, *_Agesurface; // Ponding
  grid *_d2Hchan,  *_d18Ochan, *_Agechan; // channel
  grid *_d2Hsoil1,    *_d18Osoil1, *_Agesoil1; // Vadose layer 1
  grid *_d2Hsoil2,    *_d18Osoil2, *_Agesoil2; // Vadose layer 2
  grid *_d2Hsoil_12,  *_d18Osoil_12, *_Agesoil_12; // Weighted average L1+L2
  grid *_d2Hsoil3,    *_d18Osoil3, *_Agesoil3; // Vadose layer 3
  grid *_d2HsoilAv,   *_d18OsoilAv, *_AgesoilAv; // Vadose weighted average
  grid *_d2HdeepGW,   *_d18OdeepGW, *_AgedeepGW; //Deep GW storage
  grid *_d2H_MW1,    *_d18O_MW1, *_Age_MW1; // Mobile water water 1
  grid *_d2H_MW2,    *_d18O_MW2, *_Age_MW2; // Mobile water water 2
  grid *_d2H_TB1,    *_d18O_TB1, *_Age_TB1; // Tightly-bound water 1
  grid *_d2H_TB2,    *_d18O_TB2, *_Age_TB2; // Tightly-bound water 2
  grid *_Age_MW12, *_Age_TB12;
  grid *_d2Hsoil1Q, *_d18Osoil1Q, *_Agesoil1Q; // Vadose layer 1 outflow
  grid *_d2Hsoil2Q, *_d18Osoil2Q, *_Agesoil2Q; // Vadose layer 2 outflow
  grid *_d2Hgroundwater, *_d18Ogroundwater, *_Agegroundwater; // Groundwater outflow
  grid *_d2HdeepGWQ, *_d18OdeepGWQ, *_AgedeepGWQ; //Deep groundwater outflow
  // Signature of outgoing water
  grid *_d2HevapS_sum, *_d18OevapS_sum, *_AgeevapS_sum;
  grid *_d2HevapC_sum, *_d18OevapC_sum, *_AgeevapC_sum;
  grid *_d2HevapI_sum, *_d18OevapI_sum, *_AgeevapI_sum;
  grid *_d2HevapT_sum, *_d18OevapT_sum, *_AgeevapT_sum;
  grid *_d2Hleakage, *_d18Oleakage, *_Ageleakage; // Groundwater
// Signature of vapour
  grid *_d2HevapI_Vap_sum, *_d18OevapI_Vap_sum;
  grid *_d2HevapT_Vap_sum, *_d18OevapT_Vap_sum;
  // Signature of outgoing water
  grid *_Fd2HLattoSrf, *_Fd18OLattoSrf, *_FAgeLattoSrf;
  grid *_Fd2HLattoChn, *_Fd18OLattoChn, *_FAgeLattoChn;
  grid *_Fd2HLattoGW, *_Fd18OLattoGW, *_FAgeLattoGW;
  grid *_Fd2HLattoDeepGW, *_Fd18OLattoDeepGW, *_FAgeLattoDeepGW; //input CONC. from upstream
  // Internal age contributions
  grid *_d2HGWtoChn, *_d2HDeepGWtoChn, *_d2HSrftoChn, *_d2HRecharge;
  grid *_d18OGWtoChn, *_d18ODeepGWtoChn, *_d18OSrftoChn, *_d18ORecharge;
  grid *_AgeGWtoChn, *_AgeDeepGWtoChn, *_AgeSrftoChn, *_AgeRecharge;
  //-----------------------------------------------------------------------
  // For boundary condition input
  //-----------------------------------------------------------------------
  UINT4 InitiateBCMap_Iso(ifstream & ifHandle, grid & BCMap, Atmosphere &atm);
  int UpdateBCMap_Iso(ifstream & ifHandle, grid & BCMap, Atmosphere &atm);
  ifstream ifd2HsurfaceBC, ifd18OsurfaceBC, ifAgesurfaceBC;
  ifstream ifd2Hlayer1BC,  ifd18Olayer1BC,  ifAgelayer1BC;
  ifstream ifd2Hlayer2BC,  ifd18Olayer2BC,  ifAgelayer2BC;
  ifstream ifd2HgroundwaterBC, ifd18OgroundwaterBC, ifAgegroundwaterBC;
  ifstream ifd2HdeepGWBC,  ifd18OdeepGWBC, ifAgedeepGWBC;
  grid *_d2HsurfaceBC,*_d18OsurfaceBC, *_AgesurfaceBC;
  grid *_d2Hlayer1BC, *_d18Olayer1BC,  *_Agelayer1BC;
  grid *_d2Hlayer2BC, *_d18Olayer2BC,  *_Agelayer2BC;
  grid *_d2HgroundwaterBC, *_d18OgroundwaterBC, *_AgegroundwaterBC;
  grid *_d2HdeepGWBC, *_d18OdeepGWBC, *_AgedeepGWBC;
  
  //vectors containing signature of water output for each cell with no drainage (ldd value of 5). 
  // The vectCell structure contains the row and col  
  //of that cell with no output and the area draining to that cell
  vectCells _d2HOvlndOutput, _d18OOvlndOutput, _AgeOvlndOutput; // surface output
  vectCells _d2HL1wtrOutput, _d18OL1wtrOutput, _AgeL1wtrOutput; // groundwater output
  vectCells _d2HL2wtrOutput, _d18OL2wtrOutput, _AgeL2wtrOutput; // groundwater output
  vectCells _d2HGwtrOutput, _d18OGwtrOutput, _AgeGwtrOutput; // groundwater output
  vectCells _d2HDeepGwtrOutput, _d18ODeepGwtrOutput, _AgeDeepGwtrOutput; // groundwater outlet output
  
  //check maps mainly to make sure no nodata values are in the domain.
  void CheckMapsTrck(Control &ctrl, Basin &bsn);

 public:
  
  //Constructors
  Tracking();
  Tracking(Control &ctrl, Basin &bsn, Atmosphere &atm);
  
  //Destructor
  ~Tracking();
  
  int AdvanceBCMaps_Iso(Atmosphere &atm, Control &ctrl, int iso);

  void BCupstreamMixing(Basin &bsn, Control &ctrl, double BCsrf, double BCGW,
			double BCdeepGW, double dx, double dt, int r, int c);
  
  void MixingV_canopy(Atmosphere &atm, Basin &bsn, Control &ctrl,
		      double DelCanStor, double IntercWater,double evap,
		      double &icanopy, double &ievap, double &ithro,
		      //double &icanopyW, double &ievapVW,
		      double D, double iPrec, double p, int r, int c,
		      int s, int iso);
  
  void MixingV_down(Basin &bsn, Control &ctrl, 
		    double &d1, double &d2, double &d3, double &fc,
		    int r, int c, int step);
  
  void MixingV_latup(Basin &bsn, Control &ctrl, 
		     double &d1, double &d2, double &d3, double &fc,
		     double &Qk1, double &dtdx, double &dx, int r, int c);
  

  void MixingV_richards(Basin &bsn, Control &ctrl, double d1, double d2, double d3, double dt, int reinf, int r, int c);
  
  void MixingV_chan(Basin &bsn, Control &ctrl, double dt, double Qk1, double dtdx, double dx, int r, int c);
  
  void MixingV_evapS(Atmosphere &atm, Basin &bsn, Control &ctrl, 
		     double &d1, double &theta_new,
		     double Ts, double &etp, double &beta,
		     double &d2Hevap, double &d18Oevap, double &Agevap,
		     int r, int c);

  void MixingV_sapflow(Atmosphere &atm, Basin &bsn, Control &ctrl,
		       double pTrp1, double pTrp2, double pTrp3,
		       double transp, double &iTra, double &iTraVap,
		       //double &iTraW, double &iTraVapW,
		       double p, int r, int c, int s, int iso);
  
  void TracerMixing(Basin &bsn, Control &ctrl,
		    double hold, double iold, double &inew, double qin, double iin,
		    double qout, double &iout,double S,int mixmod, int r, int c);

  void MixingV_evapW(Atmosphere &atm, Basin &bsn, Control &ctrl, 
		     double Tw, double evap, double chan_store, int r, int c);
  
  void MixingV_snow(Atmosphere &atm, Basin &bsn, Control &ctrl, double &h, double &h0, double &dh, int r, int c);
  
  void MixingV_through(Atmosphere &atm, Basin &bsn, Control &ctrl, double &rain, double &p, int s, int r, int c);
  
  void FCdownstream(Basin &bsn, Control &ctrl, double &Qk1, double &dtdx, double &dx,
		    int r, int c, int rr, int cc);
  
  void FC_Soildownstream(Basin &bsn, Control &ctrl, int r, int c, int rr, int cc);
  
  void MixingTPD_postET(Basin &bsn, Control &ctrl,
			double &dtheta, double &dtheta2,
			double &kTB_L1, double &kTB_L2,
			double &kMW_L1, double &kMW_L2,
			int r, int c);
  
  // Outlet lateral fluxes' signatures
  void OutletVals(Control &ctrl, int mode, int r, int c);

  // Generic function for both isotopes, using the 'iso' toggle: 0=deuterium, 1=oxygen18
  int Frac_Estorage(Atmosphere &atm, Basin &bsn, Control &ctrl,
		    REAL8 V_old, REAL8 V_new, REAL8 &di_old,
		    REAL8 &di_new, REAL8 &di_evap,REAL8 &beta,
		    REAL8 &Ts, int issoil, int r, int c, int iso);

  // Age increment at end of tiem step
  int IncrementAge(Basin &bsn, Control &ctrl);
  
  // Soil averaged quantities
  int Calcd2Hsoil_12(Basin &bsn);
  int Calcd18Osoil_12(Basin &bsn);
  int CalcAgesoil_12(Basin &bsn);
  int Calcd2Hsoil_Av(Basin &bsn);
  int Calcd18Osoil_Av(Basin &bsn);
  int CalcAgesoil_Av(Basin &bsn);

  // Convert grid from isotopic ratios to isotopic deltas
  void Ratio2DeltaGrid(const Basin &bsn, const grid &m, grid &mO, int iso);
  void Ratio2DeltaGrid_L1L2(const Basin &bsn, const grid &mL1, const grid &mL2, grid &mO, int iso);
  void Ratio2DeltaGrid_SoilAv(const Basin &bsn, const grid &mL1, const grid &mL2,
			      const grid &mL3, const grid &mGW, grid &mO, int iso);
  //int ReadConfigTrck(Control &ctrl, string confilename = "configTrck.ini");
  
  // Conversion from soil to two-pore and vice-versa
  int CalcInitTPD(Basin &bsn, Control &ctrl);
  int CalcTPDtoLayers(Basin &bsn, Control &ctrl);

  //Getters 2H (ratios)
  //grid *getd2Hprecip() const {
  //  return _d2Hprecip;
  //}
  grid *getd2Hthroughfall_sum() const {
    return _d2Hthroughfall_sum;
  }
  grid *getd2Hcanopy_sum() const {
    return _d2Hcanopy_sum;
  }
  grid *getd2Hsnowpack() const {
    return _d2Hsnowpack;
  }
  grid *getd2Hsurface() const {
    return _d2Hsurface;
  }
  grid *getd2HsurfaceBC() const {
    return _d2HsurfaceBC;
  }
  grid *getd2Hchan() const {
    return _d2Hchan;
  }
  grid *getd2Hsoil1() const {
    return _d2Hsoil1;
  }
  grid *getd2Hsoil2() const {
    return _d2Hsoil2;
  }
  grid *getd2Hsoil_12() const {
    return _d2Hsoil_12;
  }
  grid *getd2Hsoil3() const {
    return _d2Hsoil3;
  }
  grid *getd2Hsoil_Av() {
    return _d2HsoilAv;
  }
  grid *getd2Hgroundwater() const {
    return _d2Hgroundwater;
  }
  grid *getd2HgroundwaterBC() const {
    return _d2HgroundwaterBC;
  }
  grid *getd2HdeepGW() const {
    return _d2HdeepGW;
  }
  grid *getd2HdeepGWQ() const {
    return _d2HdeepGWQ;
  }  
  grid *getd2HdeepGWBC() const {
    return _d2HdeepGWBC;
  }  
  grid *getd2HevapS_sum() const {
    return _d2HevapS_sum;
  }
  grid *getd2HevapC_sum() const {
    return _d2HevapC_sum;
  }
  grid *getd2HevapI_sum() const {
    return _d2HevapI_sum;
  }
  grid *getd2HevapT_sum() const {
    return _d2HevapT_sum;
  }
  
  grid *getd2HevapI_Vap_sum() const {
    return _d2HevapI_Vap_sum;
  }
  grid *getd2HevapT_Vap_sum() const {
    return _d2HevapT_Vap_sum;
  }
  
  grid *getd2Hleakage() const {
    return _d2Hleakage;
  }
  grid *getd2HGWtoChn() const {
    return _d2HGWtoChn;
  }
  grid *getd2HDeepGWtoChn() const {
    return _d2HDeepGWtoChn;
  }  
  grid *getd2HSrftoChn() const {
    return _d2HSrftoChn;
  }
  grid *getd2HRecharge() const {
    return _d2HRecharge;
  }
  const vectCells *getd2HOvlndOutput() const {
    return &_d2HOvlndOutput;
  }  
  const vectCells *getd2HL1wtrOutput() const {
    return &_d2HL1wtrOutput;
  }
  const vectCells *getd2HL2wtrOutput() const {
    return &_d2HL2wtrOutput;
  }
  const vectCells *getd2HGwtrOutput() const {
    return &_d2HGwtrOutput;
  }
  const vectCells *getd2HDeepGwtrOutput() const {
    return &_d2HDeepGwtrOutput;
  }  
  grid *getd2H_MW1() const {
    return _d2H_MW1;
  }
  grid *getd2H_MW2() const {
    return _d2H_MW2;
  }
  grid *getd2H_TB1() const {
    return _d2H_TB1;
  }
  grid *getd2H_TB2() const {
    return _d2H_TB2;
  }
  // -------
  
  
  // 18O
  //grid *getd18Oprecip() const {
  //  return _d18Oprecip;
  //}
  grid *getd18Othroughfall_sum() const {
    return _d18Othroughfall_sum;
  }
  grid *getd18Ocanopy_sum() const {
    return _d18Ocanopy_sum;
  }
  grid *getd18Osnowpack() const {
    return _d18Osnowpack;
  }
  grid *getd18Osurface() const {
    return _d18Osurface;
  }
  grid *getd18OsurfaceBC() const {
    return _d18OsurfaceBC;
  }
  grid *getd18Ochan() const {
    return _d18Ochan;
  }
  grid *getd18Osoil1() const {
    return _d18Osoil1;
  }
  grid *getd18Osoil2() const {
    return _d18Osoil2;
  }
  grid *getd18Osoil_12() const {
    return _d18Osoil_12;
  }
  grid *getd18Osoil3() const {
    return _d18Osoil3;
  }
  grid *getd18Osoil_Av() const {
    return _d18OsoilAv;
  }
  grid *getd18Ogroundwater() const {
    return _d18Ogroundwater;
  }
  grid *getd18OgroundwaterBC() const {
    return _d18OgroundwaterBC;
  }
  grid *getd18OdeepGW() const {
    return _d18OdeepGW;
  }
  grid *getd18OdeepGWQ() const {
    return _d18OdeepGWQ;
  }    
  grid *getd18OdeepGWBC() const {
    return _d18OdeepGWBC;
  }  
  grid *getd18OevapS_sum() const {
    return _d18OevapS_sum;
  }
  grid *getd18OevapC_sum() const {
    return _d18OevapC_sum;
  }
  grid *getd18OevapI_sum() const {
    return _d18OevapI_sum;
  }
  grid *getd18OevapT_sum() const {
    return _d18OevapT_sum;
  }

  grid *getd18OevapI_Vap_sum() const {
    return _d18OevapI_Vap_sum;
  }
  grid *getd18OevapT_Vap_sum() const {
    return _d18OevapT_Vap_sum;
  }
  
  grid *getd18Oleakage() const {
    return _d18Oleakage;
  }
  grid *getd18OGWtoChn() const {
    return _d18OGWtoChn;
  }
  grid *getd18ODeepGWtoChn() const {
    return _d18ODeepGWtoChn;
  }  
  grid *getd18OSrftoChn() const {
    return _d18OSrftoChn;
  }
  grid *getd18ORecharge() const {
    return _d18ORecharge;
  }
  const vectCells *getd18OOvlndOutput() const {
    return &_d18OOvlndOutput;
  }  
  const vectCells *getd18OL1wtrOutput() const {
    return &_d18OL1wtrOutput;
  }
  const vectCells *getd18OL2wtrOutput() const {
    return &_d18OL2wtrOutput;
  }
  const vectCells *getd18OGwtrOutput() const {
    return &_d18OGwtrOutput;
  }
  const vectCells *getd18ODeepGwtrOutput() const {
    return &_d18ODeepGwtrOutput;
  }  
  grid *getd18O_MW1() const {
    return _d18O_MW1;
  }
  grid *getd18O_MW2() const {
    return _d18O_MW2;
  }
  grid *getd18O_TB1() const {
    return _d18O_TB1;
  }
  grid *getd18O_TB2() const {
    return _d18O_TB2;
  }
  // -------
    
  // Age
  grid *getAgethroughfall_sum() const {
    return _Agethroughfall_sum;
  }
  grid *getAgecanopy_sum() const {
    return _Agecanopy_sum;
  }
  grid *getAgesnowpack() const {
    return _Agesnowpack;
  }
  grid *getAgesurface() const {
    return _Agesurface;
  }
  grid *getAgesurfaceBC() const {
    return _AgesurfaceBC;
  }  
  grid *getAgechan() const {
    return _Agechan;
  }
  grid *getAgesoil1() const {
    return _Agesoil1;
  }
  grid *getAgesoil2() const {
    return _Agesoil2;
  }
  grid *getAgesoil_12() const {
    return _Agesoil_12;
  }
  grid *getAgesoil3() const {
    return _Agesoil3;
  }
  grid *getAgesoil_Av() const {
    return _AgesoilAv;
  }
  grid *getAgegroundwater() const {
    return _Agegroundwater;
  }
  grid *getAgegroundwaterBC() const {
    return _AgegroundwaterBC;
  }
  grid *getAgedeepGW() const {
    return _AgedeepGW;
  }
  grid *getAgedeepGWQ() const {
    return _AgedeepGWQ;
  }    
  grid *getAgedeepGWBC() const {
    return _AgedeepGWBC;
  }    
  grid *getAgeevapS_sum() const {
    return _AgeevapS_sum;
  }
  grid *getAgeevapC_sum() const {
    return _AgeevapC_sum;
  }
  grid *getAgeevapI_sum() const {
    return _AgeevapI_sum;
  }
  grid *getAgeevapT_sum() const {
    return _AgeevapT_sum;
  }
  grid *getAgeleakage() const {
    return _Ageleakage;
  }
  grid *getAgeGWtoChn() const {
    return _AgeGWtoChn;
  }
  grid *getAgeDeepGWtoChn() const {
    return _AgeDeepGWtoChn;
  }  
  grid *getAgeSrftoChn() const {
    return _AgeSrftoChn;
  }
  grid *getAgeRecharge() const {
    return _AgeRecharge;
  }

  const vectCells *getAgeOvlndOutput() const {
    return &_AgeOvlndOutput;
  }
  const vectCells *getAgeL1wtrOutput() const {
    return &_AgeL1wtrOutput;
  }
  const vectCells *getAgeL2wtrOutput() const {
    return &_AgeL2wtrOutput;
  }  
  const vectCells *getAgeGwtrOutput() const {
    return &_AgeGwtrOutput;
  }
  const vectCells *getAgeDeepGwtrOutput() const {
    return &_AgeDeepGwtrOutput;
  }  
  grid *getAge_MW1() const {
    return _Age_MW1;
  }
  grid *getAge_MW2() const {
    return _Age_MW2;
  }
  grid *getAge_MW12() const {
    return _Age_MW12;
  }
  grid *getAge_TB1() const {
    return _Age_TB1;
  }
  grid *getAge_TB2() const {
    return _Age_TB2;
  }
  grid *getAge_TB12() const {
    return _Age_TB12;
  }
  // ---

  // ---Getters of fForest getters

  grid *getd2Hcanopy(UINT4 n) const;

  grid *getd18Ocanopy(UINT4 n) const;

  grid *getAgecanopy(UINT4 n) const;

  grid *getd2Hthroughfall(UINT4 n) const;
  
  grid *getd18Othroughfall(UINT4 n) const;

  grid *getAgethroughfall(UINT4 n) const;

  // --- reset
  void resetFTracerLat(Control &ctrl) { // reset summed lateral contributions
    if(ctrl.sw_2H){
      _Fd2HLattoGW->reset();
      _Fd2HLattoChn->reset();
      _Fd2HLattoSrf->reset();
      if(ctrl.sw_deepGW)
	_Fd2HLattoDeepGW->reset();
    }
    if(ctrl.sw_18O){
      _Fd18OLattoGW->reset();
      _Fd18OLattoDeepGW->reset();
      _Fd18OLattoChn->reset();
      _Fd18OLattoSrf->reset();
      if(ctrl.sw_deepGW)
	_Fd18OLattoDeepGW->reset();
    }
    if(ctrl.sw_Age){
      _FAgeLattoGW->reset();
      _FAgeLattoDeepGW->reset();
      _FAgeLattoChn->reset();
      _FAgeLattoSrf->reset();
      if(ctrl.sw_deepGW)
	_FAgeLattoDeepGW->reset();
    }
  }  
  
  // --- Setters
  // - 2H
  void setd2Hthroughfall_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d2Hthroughfall_sum->matrix[row][col] = value;
  }  
  void setd2Hcanopy_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d2Hcanopy_sum->matrix[row][col] = value;
  }
  void setd2Hsnowpack(UINT4 row, UINT4 col, REAL8 value) {
    _d2Hsnowpack->matrix[row][col] = value;
  }
  void setd2Hsurface(UINT4 row, UINT4 col, REAL8 value) {
    _d2Hsurface->matrix[row][col] = value;
  }
  void setd2Hchan(UINT4 row, UINT4 col, REAL8 value) {
    _d2Hchan->matrix[row][col] = value;
  }
  void setd2Hsoil1(UINT4 row, UINT4 col, REAL8 value) {
    _d2Hsoil1->matrix[row][col] = value;
  }
  void setd2Hsoil2(UINT4 row, UINT4 col, REAL8 value) {
    _d2Hsoil2->matrix[row][col] = value;
  }
  void setd2Hsoil3(UINT4 row, UINT4 col, REAL8 value) {
    _d2Hsoil3->matrix[row][col] = value;
  }
  void setd2Hgroundwater(UINT4 row, UINT4 col, REAL8 value) {
    _d2Hgroundwater->matrix[row][col] = value;
  }
  void setd2HevapS_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d2HevapS_sum->matrix[row][col] = value;
  }
  void setd2HevapI_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d2HevapI_sum->matrix[row][col] = value;
  }
  void setd2HevapT_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d2HevapT_sum->matrix[row][col] = value;
  }
  void setd2HevapI_Vap_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d2HevapI_Vap_sum->matrix[row][col] = value;
  }
  void setd2HevapT_Vap_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d2HevapT_Vap_sum->matrix[row][col] = value;
  }

  // - 18O
  void setd18Othroughfall_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d18Othroughfall_sum->matrix[row][col] = value;
  }    
  void setd18Ocanopy_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d18Ocanopy_sum->matrix[row][col] = value;
  }  
  void setd18Osnowpack(UINT4 row, UINT4 col, REAL8 value) {
    _d18Osnowpack->matrix[row][col] = value;
  }
  void setd18Osurface(UINT4 row, UINT4 col, REAL8 value) {
    _d18Osurface->matrix[row][col] = value;
  }
  void setd18Ochan(UINT4 row, UINT4 col, REAL8 value) {
    _d18Ochan->matrix[row][col] = value;
  }
  void setd18Osoil1(UINT4 row, UINT4 col, REAL8 value) {
    _d18Osoil1->matrix[row][col] = value;
  }
  void setd18Osoil2(UINT4 row, UINT4 col, REAL8 value) {
    _d18Osoil2->matrix[row][col] = value;
  }
  void setd18Osoil3(UINT4 row, UINT4 col, REAL8 value) {
    _d18Osoil3->matrix[row][col] = value;
  }
  void setd18Ogroundwater(UINT4 row, UINT4 col, REAL8 value) {
    _d18Ogroundwater->matrix[row][col] = value;
  }
  void setd18OevapS_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d18OevapS_sum->matrix[row][col] = value;
  }
  void setd18OevapI_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d18OevapI_sum->matrix[row][col] = value;
  }
  void setd18OevapT_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d18OevapT_sum->matrix[row][col] = value;
  }
  void setd18OevapI_Vap_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d18OevapI_Vap_sum->matrix[row][col] = value;
  }
  void setd18OevapT_Vap_sum(UINT4 row, UINT4 col, REAL8 value) {
    _d18OevapT_Vap_sum->matrix[row][col] = value;
  }  

  // - Age
  void setAgethroughfall_sum(UINT4 row, UINT4 col, REAL8 value) {
    _Agethroughfall_sum->matrix[row][col] = value;
  }    
  void setAgecanopy_sum(UINT4 row, UINT4 col, REAL8 value) {
    _Agecanopy_sum->matrix[row][col] = value;
  }    
  void setAgesnowpack(UINT4 row, UINT4 col, REAL8 value) {
    _Agesnowpack->matrix[row][col] = value;
  }
  void setAgesurface(UINT4 row, UINT4 col, REAL8 value) {
    _Agesurface->matrix[row][col] = value;
  }
  void setAgechan(UINT4 row, UINT4 col, REAL8 value) {
    _Agechan->matrix[row][col] = value;
  }
  void setAgesoil1(UINT4 row, UINT4 col, REAL8 value) {
    _Agesoil1->matrix[row][col] = value;
  }
  void setAgesoil2(UINT4 row, UINT4 col, REAL8 value) {
    _Agesoil2->matrix[row][col] = value;
  }
  void setAgesoil3(UINT4 row, UINT4 col, REAL8 value) {
    _Agesoil3->matrix[row][col] = value;
  }
  void setAgegroundwater(UINT4 row, UINT4 col, REAL8 value) {
    _Agegroundwater->matrix[row][col] = value;
  }
  void setAgeevapS_sum(UINT4 row, UINT4 col, REAL8 value) {
    _AgeevapS_sum->matrix[row][col] = value;
  }
  void setAgeevapI_sum(UINT4 row, UINT4 col, REAL8 value) {
    _AgeevapI_sum->matrix[row][col] = value;
  }
  void setAgeevapT_sum(UINT4 row, UINT4 col, REAL8 value) {
    _AgeevapT_sum->matrix[row][col] = value;
  }

  // ---- Mixing equations ------------------------------------------------------------------
  // When there's only input
  double InputMix(double hold, double iold, double qin, double iin){

    double inew = 0;
    // Assign fluxes
    if(hold + qin > RNDOFFERR)
      inew = (iold*hold + iin*qin)/ (hold + qin);
    else
      inew = iold;
    return inew;
  }
  
  // Mixing when one input one output: two modes to apporximate what aggregated h and i
  // values the DE approximation uses.
  // If mode = 0, h=hold and i=inew. Mixing without change of volume due to output.
  // If mode = 1, h=0.5*(hold+hnew) and i=0.5*(inew+iold). A somewhat linear model of mixing.
  // If mode = 2, h=0.5*(hold+hnew) and i=(hold*inew+hnew*iold)/(hold+hnew). 
  // A more complex, and stable weighted-linear model of mixing.
  // 1000 is added (and then suvstracted) to avoid numerical instabilities with isotopes
  // signatures (it is transparent for age)
  double InOutMix(double hold, double iold, double qin, double iin, double qout, int mode){

    double inew = 0;
    double hsum = 0;

    if(mode==0)
      inew = hold + qin > RNDOFFERR ?
	(iold*hold + iin*qin)/ (hold + qin) : iold;
    else if(mode ==1) {
      hsum = 0.5*(hold+qin+std::max<double>(0,hold-qout));
      inew = 0.5*hsum+qin > RNDOFFERR ? 
	((iold+1000)*(hsum-0.5*qin) + (iin+1000)*qin)/ (hsum + 0.5*qin) -1000 : iold;
    }

    return inew;
  }
  // ------------------------------------------------------------------------------------------
};

#endif /* TRACKING_H_ */
