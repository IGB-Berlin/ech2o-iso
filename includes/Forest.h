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
 * Forest.h
 *
 *  Created on: May 5, 2010
 *      Author: Marco.Maneta
 */

#ifndef FOREST_H_
#define FOREST_H_

#include "Grid.h"
#include "InitConf.h"
#include "ConstAndFuncs.h"
#include "SortGrid.h"
#include "Grove.h"
#include "Atmosphere.h"
#include "Basin.h"
#include "Tracking.h"

using namespace std;
class Basin;
class Tracking;
class Forest {

  UINT4 _NRows;
  UINT4 _NCols;
  REAL8 _dx;
  UINT4 _Nsp; //number of species included in the simulation including clear cut (bare) soil
  grid *_patches; //map with forest patches
  vectCells _vSortedGrid;

  Grove *_species;

  vectCells SortGrid();

  void checkForestDatabase();

  int SetStateVarsTabs(Control &ctrl);
  void SetStateVarsMaps(Control &ctrl);
  int SetSpeciesParameters(Control &ctrl);

  // For LAI input ----------------------------------------------
  UINT4 InitiateLAIMap(ifstream & ifHandle, grid & ClimMap);
  int UpdateLAIMap(ifstream & ifHandle, grid & ClimMap);	  // internal function that updates a LAI map
  ifstream ifLAI; 						  // LAI data file handles
  ifstream ifhgt; 						  // hgt data file handles  
  // -------------------------------------------------------------

  /*State variables*/

  int GrowTrees(UINT4 j, UINT4 r, UINT4 c, REAL8 dt, REAL8 fa, REAL8 ft, REAL8 fw, REAL8 T, REAL8 usablewater);
  int GrowStem(UINT4 spec, UINT4 row, UINT4 col);
  int GrowLAI(UINT4 spec, UINT4 row, UINT4 col, REAL8 T, REAL8 waterusable, REAL8 dt);
  int GrowRoots(UINT4 spec, UINT4 row, UINT4 col, REAL8 dt);

  int GrowGrass(UINT4 j, UINT4 r, UINT4 c, REAL8 dt);
  int GrowGrassLAI(UINT4 spec, UINT4 row, UINT4 col, REAL8 dt);

  double NetRadCanopy(Atmosphere &atm, const double &Ts, REAL8 emiss, REAL8 albedo, REAL8 Kbeers, REAL8 lai, int row, int col);
  double LatHeatCanopy(Basin &bas, Atmosphere &atm, double soilrelhumid, double ra, const double &Ts, int row, int col);
  double SensHeatCanopy(Atmosphere &atm, double ra, const double &Ts, int row, int col);

 public:

  Forest();
  Forest(Control &ctrl);
  ~Forest();

  int CalculateCanopyConduct(const Basin &bas, const Atmosphere &atm, const Control &ctrl, const double &lwp,const double &lwp_max,
			     const double &lwp_min, double &dgsdlwp, UINT4 j, UINT4 r, UINT4 c);

  UINT4 SolveCanopyEnergyBalance(Basin &bas, Atmosphere &atm, Control &ctrl,
				 REAL8 theta_r1, REAL8 theta_r2, REAL8 theta_r3,
				 REAL8 rootdepth,  REAL8 psiae, REAL8 bclambda, REAL8 ra, 
				 REAL8 &DelCanStor, REAL8 &evap_a, REAL8 &transp_a, REAL8 &netR_a, REAL8 &leR_a,
				 REAL8 &h_a, UINT4 s, UINT4 r, UINT4 c);
  int CanopyInterception(Atmosphere &atm, Control &ctrl, REAL8 &DelCanStor, REAL8 &D, UINT4 s, UINT4 r, UINT4 c);
  int GrowForest(Basin &bas, const Atmosphere &atm, const Control &ctrl);
  int SperryModel(Basin &bas, Atmosphere &atm, Control &ctrl, REAL8 rootdepth,REAL8 Sold, REAL8 Keff, REAL8 psiae, REAL8 bclambda, REAL8 airTp,
		  REAL8 airRH, REAL8 rho_a, REAL8 gamma, REAL8 ra, REAL8 poros, REAL8 thetar, REAL8 thetawp, REAL8 evap_a,REAL8 fA, REAL8 fB, 
		  REAL8 fC, REAL8 leavesurfRH, REAL8 leafRH,REAL8 &LET, REAL8 &LE, REAL8 &H,REAL8 &temp0,REAL8 &temp1,REAL8 &temp2, REAL8 &temp3, UINT4 s, 
		  UINT4 r, UINT4 c);

  // Convert grid from isotopic ratios to isotopic deltas
  void Ratio2DeltaGrid(const Basin &bas, const grid &m, grid &mO, int iso);   

  //getters
  UINT4 getNumSpecies() const {
    return _Nsp;
  }

  REAL8 getPropSpecies(UINT4 n, UINT4 row, UINT4 col) const {
    return _species[n]._fraction->matrix[row][col];
  }

  REAL8 getVegSpeciesType(UINT4 n) const {
    return _species[n].vegtype;
  }

  REAL8 getLAISpecies(UINT4 n, UINT4 row, UINT4 col) const {
    return _species[n]._LAI->matrix[row][col];
  }

  REAL8 getMaxCanopyStorage(UINT4 n, UINT4 row, UINT4 col) const {
    Grove *spe = &_species[n]; //ACHTUNG ACHTUNG! MEMORY LEAK??!?!?!
    return spe->MaxCanStorageParamt * spe->_LAI->matrix[row][col];
  }

  REAL8 getThroughfallCoeff(UINT4 n) const {
    return _species[n].Throughfall_coeff;
  }

  REAL8 getCanopyConductance(UINT4 n, UINT4 row, UINT4 col) const {
    Grove *spe = &_species[n]; //ACHTUNG ACHTUNG! MEMORY LEAK??!?!?!
    return spe->_CanopyConductance->matrix[row][col];
  }
  REAL8 getTreeHeight(UINT4 n, UINT4 row, UINT4 col) const {
    Grove *spe = &_species[n]; //ACHTUNG ACHTUNG! MEMORY LEAK??!?!?!
    return spe->_Height->matrix[row][col];
  }

  REAL8 getCanopyTemp(UINT4 n, UINT4 row, UINT4 col) const {
    Grove *spe = &_species[n]; //ACHTUNG ACHTUNG! MEMORY LEAK??!?!?!
    return spe->_Temp_c->matrix[row][col];
  }
  REAL8 getRootUptakeL1(UINT4 n, UINT4 row, UINT4 col) const {
    Grove *spe = &_species[n];
    return spe->_RUptakeL1->matrix[row][col];
  }

  REAL8 getRootUptakeL2(UINT4 n, UINT4 row, UINT4 col) const {
    Grove *spe = &_species[n];
    return spe->_RUptakeL2->matrix[row][col];
  }

  REAL8 getRootUptakeL3(UINT4 n, UINT4 row, UINT4 col) const {
    Grove *spe = &_species[n];
    return spe->_RUptakeL3->matrix[row][col];
  }
  REAL8 getCanopyNetRad(UINT4 n, UINT4 row, UINT4 col) const {
    Grove *spe = &_species[n]; //ACHTUNG ACHTUNG! MEMORY LEAK??!?!?!
    return spe->_NetR_Can->matrix[row][col];
  }

  REAL8 getCanopyLatHeatE(UINT4 n, UINT4 row, UINT4 col) const {
    Grove *spe = &_species[n]; //ACHTUNG ACHTUNG! MEMORY LEAK??!?!?!
    return spe->_LatHeat_CanE->matrix[row][col];
  }
  REAL8 getCanopyLatHeatT(UINT4 n, UINT4 row, UINT4 col) const {
    Grove *spe = &_species[n]; //ACHTUNG ACHTUNG! MEMORY LEAK??!?!?!
    return spe->_LatHeat_CanT->matrix[row][col];
  }

  REAL8 getCanopySensHeat(UINT4 n, UINT4 row, UINT4 col) const {
    Grove *spe = &_species[n]; //ACHTUNG ACHTUNG! MEMORY LEAK??!?!?!
    return spe->_SensHeat_Can->matrix[row][col];
  }

  REAL8 getIntercWater(UINT4 n, UINT4 row, UINT4 col) const {
    return _species[n]._WaterStorage->matrix[row][col];
  }

  REAL8 getBeersCoeff(UINT4 n, UINT4 row, UINT4 col) const {
    return _species[n].KBeers;
  }
    
  REAL8 getKRoot(UINT4 n) const {
    return _species[n].Kroot;
  }

  REAL8 getWiltPoint(UINT4 n) const {
    return _species[n].WiltingPoint;
  }

  REAL8 getRootAspect(UINT4 n) const {
    return _species[n].Aroot;
  }

  REAL8 getSperry_d(UINT4 n, UINT4 row, UINT4 col) const {
    return _species[n].sperry_d;
  }
  REAL8 getSperry_c(UINT4 n, UINT4 row, UINT4 col) const {
    return _species[n].sperry_c;
  }
  REAL8 getSperry_Kp(UINT4 n, UINT4 row, UINT4 col) const {
    return _species[n].sperry_Kp;
  }
  REAL8 getRAI_a(UINT4 n, UINT4 row, UINT4 col) const {
    return _species[n].RAI_a;
  }

  REAL8 getCanopyEmissivity(UINT4 n, UINT4 row, UINT4 col) const {
    return _species[n].emissivity;
  }

  REAL8 getEvapoTransp(UINT4 n, UINT4 row, UINT4 col) const {
    return _species[n]._ET->matrix[row][col];
  }
  REAL8 getEinterception(UINT4 n, UINT4 row, UINT4 col) const {
    return _species[n]._Einterception->matrix[row][col];
  }
  REAL8 getTranspiration(UINT4 n, UINT4 row, UINT4 col) const {
    return _species[n]._Transpiration->matrix[row][col];
  }
  REAL8 getTranspirationFlux(UINT4 n, UINT4 row, UINT4 col) const {
    return _species[n]._TranspirationFlux->matrix[row][col];
  }
  REAL8 getEsoil(UINT4 n, UINT4 row, UINT4 col) const {
    return _species[n]._Esoil->matrix[row][col];
  }
  REAL8 getLeafWaterPotential(UINT4 n, UINT4 row, UINT4 col) const {
    return _species[n]._LeafWatPot->matrix[row][col];
  }
  REAL8 getSoilWaterPotential(UINT4 n, UINT4 row, UINT4 col) const {
    return _species[n]._SoilWatPot->matrix[row][col];
  }  
  grid *getRootFrac1(UINT4 n) const {
    return _species[n]._rootfrac1;
  }
  grid *getRootFrac2(UINT4 n) const {
    return _species[n]._rootfrac2;
  }


  grid *getLAISpeciesMap(UINT4 n) const {
    return _species[n]._LAI;
  }

  grid *getPropSpeciesMap(UINT4 n) const {
    return _species[n]._fraction;
  }

  grid *getStemDensSpeciesMap(UINT4 n) const {
    return _species[n]._StemDensity;
  }

  grid *getAgeSpeciesMap(UINT4 n) const {
    return _species[n]._AGE;
  }

  grid *getCanopyConductSpeciesMap(UINT4 n) const {
    return _species[n]._CanopyConductance;
  }

  grid *getGPPSpeciesMap(UINT4 n) const {
    return _species[n]._GPP;
  }

   grid *getNPPSpeciesMap(UINT4 n) const {
    return _species[n]._NPP;
  }

  grid *getBasalAreaSpeciesMap(UINT4 n) const {
    return _species[n]._BasalArea;
  }

  grid *getTreeHeightSpeciesMap(UINT4 n) const {
    return _species[n]._Height;
  }

  grid *getRootMassSpeciesMap(UINT4 n) const {
    return _species[n]._RootMass;
  }

  grid *getCanopyTempSpeciesMap(UINT4 n) const {
    return _species[n]._Temp_c;
  }

  grid *getRootUptakeL1Map(UINT4 n) const {
    return _species[n]._RUptakeL1;
  }

  grid *getRootUptakeL2Map(UINT4 n) const {
    return _species[n]._RUptakeL2;
  }

  grid *getRootUptakeL3Map(UINT4 n) const {
    return _species[n]._RUptakeL3;
  }
  
  grid *getCanopyNetRadSpeciesMap(UINT4 n) const {
    return _species[n]._NetR_Can;
  }

  grid *getCanopyLatHeatESpeciesMap(UINT4 n) const {
    return _species[n]._LatHeat_CanE;
  }

  grid *getCanopyLatHeatTSpeciesMap(UINT4 n) const {
    return _species[n]._LatHeat_CanT;
  }

  grid *getCanopySensHeatSpeciesMap(UINT4 n) const {
    return _species[n]._SensHeat_Can;
  }

  grid *getCanopyWaterStorSpeciesMap(UINT4 n) const {
    return _species[n]._WaterStorage;
  }
  grid *getEinterceptionSpeciesMap(UINT4 n) const {
    return _species[n]._Einterception;
  }
  grid *getETSpeciesMap(UINT4 n) const {
    return _species[n]._ET;
  }
  grid *getTranspirationSpeciesMap(UINT4 n) const {
    return _species[n]._Transpiration;
  }
  grid *getTranspirationFluxSpeciesMap(UINT4 n) const {
    return _species[n]._TranspirationFlux;
  }
  grid *getEsoilSpeciesMap(UINT4 n) const {
    return _species[n]._Esoil;
  }

  grid *getSoilWaterPotSpeciesMap(UINT4 n) const {
    return _species[n]._SoilWatPot;
  }
  
  grid *getLeafWaterPotSpeciesMap(UINT4 n) const {
    return _species[n]._LeafWatPot;
  }
  grid *getSapVelocitySpeciesMap(UINT4 n) const {
    return _species[n]._SapVelocity;
  }  

  // setters
  void setEsoilSpecies(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._Esoil->matrix[row][col] = value;
  }
  void setETSpecies(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._ET->matrix[row][col] = value;
  }
  void setTranspirationFlow(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._TranspirationFlux->matrix[row][col] = value;
  }
  void setRootFrac1Species(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._rootfrac1->matrix[row][col] = value;
  }
  void setRootFrac2Species(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._rootfrac2->matrix[row][col] = value;
  }

  // -- Tracking
  // getters
  grid *getd2Hcanopy(UINT4 n) const {
    return _species[n]._d2Hcanopy;
  }
  grid *getd2Hthroughfall(UINT4 n) const {
    return _species[n]._d2Hthroughfall;
  }  
  grid *getd2HevapT(UINT4 n) const {
    return _species[n]._d2HevapT;
  }
  grid *getd2HevapI(UINT4 n) const {
    return _species[n]._d2HevapI;
  }
  grid *getd2HevapS(UINT4 n) const {
    return _species[n]._d2HevapS;
  }
  grid *getd2HevapT_Vap(UINT4 n) const {
    return _species[n]._d2HevapT_Vap;
  }
  grid *getd2HevapI_Vap(UINT4 n) const {
    return _species[n]._d2HevapI_Vap;
  }
  
  grid *getd18Ocanopy(UINT4 n) const {
    return _species[n]._d18Ocanopy;
  }
  grid *getd18Othroughfall(UINT4 n) const {
    return _species[n]._d18Othroughfall;
  }  
  grid *getd18OevapT(UINT4 n) const {
    return _species[n]._d18OevapT;
  }
  grid *getd18OevapI(UINT4 n) const {
    return _species[n]._d18OevapI;
  }
  grid *getd18OevapS(UINT4 n) const {
    return _species[n]._d18OevapS;
  }
  grid *getd18OevapT_Vap(UINT4 n) const {
    return _species[n]._d18OevapT_Vap;
  }
  grid *getd18OevapI_Vap(UINT4 n) const {
    return _species[n]._d18OevapI_Vap;
  }
  
  grid *getAgecanopy(UINT4 n) const {
    return _species[n]._Agecanopy;
  }
  grid *getAgethroughfall(UINT4 n) const {
    return _species[n]._Agethroughfall;
  }  
  grid *getAgeevapS(UINT4 n) const {
    return _species[n]._AgeevapS;
  }
  grid *getAgeevapI(UINT4 n) const {
    return _species[n]._AgeevapI;
  }
  grid *getAgeevapT(UINT4 n) const {
    return _species[n]._AgeevapT;
  }

  // setters
  void setd2Hcanopy(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._d2Hcanopy->matrix[row][col] = value;
  }
  void setd2Hthroughfall(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._d2Hthroughfall->matrix[row][col] = value;
  }  
  void setd2HevapS(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._d2HevapS->matrix[row][col] = value;
  }
  void setd2HevapI(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._d2HevapI->matrix[row][col] = value;
  }
  void setd2HevapT(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._d2HevapT->matrix[row][col] = value;
  }
  void setd2HevapI_Vap(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._d2HevapI_Vap->matrix[row][col] = value;
  }
  void setd2HevapT_Vap(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._d2HevapT_Vap->matrix[row][col] = value;
  }  
  
  void setd18Ocanopy(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._d18Ocanopy->matrix[row][col] = value;
  }
  void setd18Othroughfall(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._d18Othroughfall->matrix[row][col] = value;
  }  
  void setd18OevapS(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._d18OevapS->matrix[row][col] = value;
  }
  void setd18OevapI(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._d18OevapI->matrix[row][col] = value;
  } 
  void setd18OevapT(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._d18OevapT->matrix[row][col] = value;
  }
  void setd18OevapI_Vap(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._d18OevapI_Vap->matrix[row][col] = value;
  } 
  void setd18OevapT_Vap(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._d18OevapT_Vap->matrix[row][col] = value;
  }
  
  void setAgecanopy(UINT4 n, UINT4 row, UINT4 col, REAL8 value) { 
    _species[n]._Agecanopy->matrix[row][col] = value;
  }
  void setAgethroughfall(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._Agethroughfall->matrix[row][col] = value;
  }
  void setAgeevapS(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._AgeevapS->matrix[row][col] = value;
  }
  void setAgeevapI(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._AgeevapI->matrix[row][col] = value;
  }
  void setAgeevapT(UINT4 n, UINT4 row, UINT4 col, REAL8 value) {
    _species[n]._AgeevapT->matrix[row][col] = value;
  }

  //external interface that updates all LAI maps by calling UpdateLAIMap
  int AdvanceLAIMaps(); 
  
};

#endif /* FOREST_H_ */
