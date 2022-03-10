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
 * grove.h
 *
 *  Created on: Jun 17, 2010
 *      Author: marco
 */

#ifndef GROVE_H_
#define GROVE_H_

#include <fstream>
#include "Grid.h"

struct Grove {

  /*Constant parameters*/
  REAL8 ID; //species ID
  REAL8 GPP2NPP; //GPP to NPP ratio
  REAL8 alpha; //Canopy quantum efficiency (g C J-1)
  REAL8 gsmax; //maximum stomatal conductance (m s-1)
  REAL8 MaxAge; //maximum age of species
  REAL8 TempOpt; //optimal temperature of species degC
  REAL8 TempMax; //maximum temperature in which specie can grow degC
  REAL8 TempMin; //maximum temperature in which specie can grow degC
  REAL8 Fprn; //coefficient n for allocation of NPP to foliage if vegtype = 0; dead grass leaf turnover rate if vegtype = 1; minimum root allocation factor if vegtype =2
  REAL8 Fpra; //coefficient a for allocation of NPP to foliage if vegtype =0; dead grass leaf turnover adjustment rate if vegtype = 1;  minimum stem allocation factor if vegtype =2
  REAL8 Sprn; //coefficient n for allocation of NPP to stems if vegtype =1; allocation parameter (modulates water and light effect) if vegtype =2
  REAL8 Spra; //coefficient a for allocation of NPP to stems (if vegtype = 0)
  REAL8 gs_light_coeff; //curvature coefficient for light factor in stomatal conductance
  REAL8 gs_vpd_coeff; //curvature coefficient for vapor pressure deficit factor in stomatal conductance
  REAL8 lwp_min; //lwp level for complete closure on stomatal conductance (- MPa)
  REAL8 lwp_max; //lwp level for no hydrologic control on stomatal conductance (- MPa)
  REAL8 WiltingPoint; //volumetric soil moisture at wilting point for the species
  REAL8 SLA; //specific leaf area (m2 g-1)
  REAL8 SRA; //specific root area (m2 g-1)
  REAL8 Crown2Stem; //crown to stem diameter ratio (m m-1)
  REAL8 TreeShapePar; //tree shape parameter from TREEDYN3
  REAL8 WoodDensity; //wood density (gC m-3)
  REAL8 Fhdmin; //minimum tree growth factor [0,1]
  REAL8 Fhdmax; //maximum tree growth factor [0,1]
  REAL8 LeafTurnover; //leaf turnover rate (s-1)
  REAL8 RootTurnover; //root turnover rate (s-1)
  REAL8 MaxCanStorageParamt; //parameter to calcualte maximum canopy storage (after Dickinson 1984 via Liang et al. 1994)
  REAL8 Throughfall_coeff; //parameter to allow direct throughfall - not tied to vegetation proportions  
  REAL8 albedo; //albedo of dry canopy
  REAL8 emissivity; //emissivity of dry canopy
  REAL8 KBeers; //Extinction coefficient in Beers law of light extinction
  REAL8 Kroot; // Decrease coefficient for exponential root profile (m-1)
  REAL8 Aroot; // aspect ratio of radial distance of roots from stem (max distance related to Kroot)
  REAL8 beta; //canopy water efficiency (gC m-1)
  // Sperry parameters
  REAL8 sperry_d; // Sperry model scaling parameter (m)
  REAL8 sperry_c; // Sperry model shape exponent (m)
  REAL8 sperry_Kp; // Sperry model hydraulic conductivity tissue (ms-1) used to calculate conductaance using DBH and plant height
  REAL8 RAI_a; // parameter to scale effective RAI as per Daly et al. (2004). Coupled Dynamics of Photosynthesis,.... J Hydromet.
  UINT4 vegtype; // 0 = evergreen tree . 1 = grass . 2 = deciduous tree
  REAL8 MaxLeafTurnoverWaterStress; //maximum leaf turnover rate due to water stress (s-1)
  REAL8 LeafTurnoverWaterStressShpParam; //Shape parameter leaf turnover rate due to water stress (-)
  REAL8 MaxLeafTurnoverColdStress; //maximum leaf turnover rate due to cold stress (s-1)
  REAL8 ColdStressTemp; //maximum leaf turnover rate due to cold stress (C)
  REAL8 LeafTurnoverColdStressShpParam; //Shape  parameter leaf turnover rate due to cold stress (-)

  /*state variables*/
/*state variables*/
  grid *_fraction;
  grid *_StemDensity; //number of trees per square meter (trees m-2)
  grid *_LAI;
  grid *_grassLAI_g;//Leaf area index of green grass
  grid *_grassLAI_d;//Leaf area index of dry grass
  grid *_AGE; //stand age in years
  grid *_CanopyConductance;//canopy conductance in m s-1
  grid *_GPP; //gross primary production in gCm-2
  grid *_NPP; //net primary production in gCm-2
  grid *_BasalArea; //Average stem cross-sectional area (m2)
  grid *_Height; //tree height (m)
  grid *_RootMass; //density of roots gC m-2
  grid *_Del_FoliageMass; //increment of leaf mass gCm-2
  grid *_Del_StemMass; //increment of stem mass gCm-2
  grid *_Del_RootMass; //increment of root mass gCm-2
  grid *_Temp_c; //canopy temperature C
  grid *_RUptakeL1; //proportion of root uptake from layer 1
  grid *_RUptakeL2; //proportion of root uptake from layer 2
  grid *_RUptakeL3; //proportion of root uptake from layer 3  
  grid *_NetR_Can;//canopy net radiation Wm-2
  grid *_LatHeat_CanE;//Canopy latent heat E interception Wm-2
  grid *_LatHeat_CanT;//Canopy latent heat transpiration Wm-2
  grid *_SensHeat_Can;//Canopy sensible heat Wm-2
  grid *_WaterStorage; //current water stored in canopy m
  grid *_ET; //Actual evapotranspiration ms-1
  grid *_Einterception; //evaporation if interception component ms-1
  grid *_Transpiration; //transpiration component ms-1
  grid *_TranspirationFlux; //transpiration flux component m3.s-1
  grid *_Esoil; // soil evaporation component m.s-1
  grid *_SoilWatPot; // soil water potential (negative m of head)
  grid *_LeafWatPot; // leaf water potential (negative m of head)
  grid *_SapVelocity; //sap water velocity [m.s-1]
  grid *_rootfrac1; // root fraction in first layer
  grid *_rootfrac2; // root fraction in second layer

  ifstream ifLAI; // LAI files handle
  ifstream ifhgt; // hgt files handle
  
  // Tracking
  grid *_d2Hcanopy; // d2H of interception water
  grid *_d18Ocanopy; // d18O of interception water
  grid *_Agecanopy; // Age of interception water

  grid *_d2Hthroughfall; // d2H of throughfall water
  grid *_d18Othroughfall; // d18O of throughfall water
  grid *_Agethroughfall; // Age of throughfall water

  grid *_d2HevapT; // d2H of transpirated water
  grid *_d18OevapT; // d18O of transpirated water
  grid *_AgeevapT; // Age of transpirated water
  grid *_d2HevapT_Vap; // d2H of transpirated water
  grid *_d18OevapT_Vap; // d18O of transpirated water  

  grid *_d2HevapI; // d2H of evaporated interception
  grid *_d18OevapI; // d18O of evaporated interception
  grid *_AgeevapI; // Age of evaporated interception
  grid *_d2HevapI_Vap; // d2H of evaporated interception
  grid *_d18OevapI_Vap; // d18O of evaporated interception  

  grid *_d2HevapS; // d2H of evaporated soil water (below each canopy type)
  grid *_d18OevapS; // d18O of evaporated soil water (below each canopy type)
  grid *_AgeevapS; // Age of evaporated soil water (below each canopy type)
      
  Grove();
  ~Grove();

  int CreateGrids(grid *base);
  int CreateGridsd2H(grid *base);
  int CreateGridsd18O(grid *base);
  int CreateGridsAge(grid *base);

};


#endif /* GROVE_H_ */
