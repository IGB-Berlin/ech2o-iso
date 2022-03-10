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
 * GroveConstruct.cpp
 *
 *  Created on: Jun 17, 2010
 *      Author: marco
 */

#include "Grove.h"

Grove::Grove(){

  _fraction = NULL;                   //Species fraction in pixel
  _StemDensity = NULL;                //stem density []
  _LAI = NULL;                        //leaf area index [m2.m-2]
  _grassLAI_g = NULL;                 //grass LAI green decay []
  _grassLAI_d = NULL;                 //grass LAI dry decay []
  _AGE = NULL;                        //vegetation age [years]
  _CanopyConductance = NULL;          //canopy conductance of vegetation []
  _GPP = NULL;                        //gross primary production []
  _NPP = NULL;                        //net primary production []
  _BasalArea = NULL;                  //tree basal area []
  _Height = NULL;                     //vegetation height [m]
  _RootMass = NULL;                   //vegetation root mass []
  _Del_FoliageMass = NULL;            //increment foliage mass [gC.m-2]
  _Del_StemMass = NULL;               //increment stem mass [gC.m-2]
  _Del_RootMass = NULL;               //increment root mass [gC.m-2]
  _Temp_c = NULL;                     //canopy temperature [C]
  _RUptakeL1 = NULL;                  //Root-uptake from layer 1 [ratio]
  _RUptakeL2 = NULL;                  //Root-uptake from layer 2 [ratio]
  _RUptakeL3 = NULL;                  //Root-uptake from layer 3 [ratio]  
  _NetR_Can = NULL;                   //canopy net radiation [W.m-2]
  _LatHeat_CanE = NULL;               //latent heat of interception evaporation [W.m-2]
  _LatHeat_CanT = NULL;               //latent heat of transpiration [W.m-2]
  _SensHeat_Can = NULL;               //sensible heat of canopy [W.m-2]
  _WaterStorage = NULL;               //water storage in the canopy [m]
  _ET = NULL;                         //total evapotranspiration [m.s-1]
  _Einterception = NULL;              //interception evaporation [m.s-1]
  _Transpiration = NULL;              //transpiration [m.s-1]
  _TranspirationFlux = NULL;          //transpiration flux [m3.s-1]
  _Esoil = NULL;                      //soil evaporation [m.s-1]
  _SoilWatPot = NULL;                 //soil water potential [m]
  _LeafWatPot = NULL;                 //leaf water potential [m]
  _SapVelocity = NULL;                //Sap water velocity [m.s-1]
  _rootfrac1 = NULL;                  //rooting fraction in layer 1 [-]
  _rootfrac2 = NULL;                  //rooting fraction in layer 2 [-]
  
  //------------------------------------------------------------------------------------
  // Tracking
  //------------------------------------------------------------------------------------
  _d2Hcanopy = NULL;                  //deuterium of the canopy water
  _d18Ocanopy = NULL;                 //oxygen-18 of the canopy water
  _Agecanopy = NULL;                  //age of the canopy water
  _d2Hthroughfall = NULL;                 //deuterium of throughfall
  _d18Othroughfall = NULL;                //oxygen-18 of throughfall
  _Agethroughfall = NULL;                 //Age of throughfall
  _d2HevapT = NULL;                   //deuterium of the transpiration
  _d18OevapT = NULL;                  //oxygen-18 of the transpiration
  _AgeevapT = NULL;                   //age of the transpiration
  _d2HevapI = NULL;                   //deuterium of the interception evaporation
  _d18OevapI = NULL;                  //oxygen-18 of the interception evaporation
  _AgeevapI = NULL;                   //age of the interception evaporation
  _d2HevapS = NULL;                   //deuterium of the soil evaporation
  _d18OevapS = NULL;                  //oxygen-18 of the soil evaporation
  _AgeevapS = NULL;                   //age of the soil evaporation

  _d2HevapT_Vap = NULL;               //deuterium of the transpiration vapour
  _d18OevapT_Vap = NULL;              //oxygen-18 of the transpiration vapour
  _d2HevapI_Vap = NULL;               //deuterium of the interception evaporation vapour
  _d18OevapI_Vap = NULL;              //oxygen-18 of the interception evaporation vapour
  
}
