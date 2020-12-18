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
 *    Marco Maneta
 *******************************************************************************/
/*
 * TotalStorage.cpp
 *
 *  Created on: Mar 18, 2010
 *      Author: Marco Maneta
 */

#include "Budget.h"

void Budget::TotalStorage( const grid *Canopy,
			   const grid *Snow,
			   const grid *Ponding,
			   const grid *ChanStore,
			   const grid *SoilL1,
			   const grid *SoilL2,
			   const grid *SoilL3,
			   //const grid *GravWater,
			   const grid *GrndWater,
			   const grid *ProotzoneL1,
			   const grid *ProotzoneL2,
			   const grid *ProotzoneL3,
			   const Basin *b)
{
	canopy = AccountStorages(Canopy, b);
	snowpack = AccountStorages(Snow, b);
	ponding = AccountStorages(Ponding, b);
	chan_store = AccountStorages(ChanStore, b);
	soilL1 = AccountStorages(SoilL1, b);
	soilL2 = AccountStorages(SoilL2, b);
	soilL3 = AccountStorages(SoilL3, b);
	//gravwater = AccountStorages(GravWater, b);
	grndwater = AccountStorages(GrndWater, b);
	vadose = soilL1 + soilL2 + soilL3 + grndwater;
	// Root zone storage: sums contributions of respective layers.
	// These fractions do not represent root fraction, but pondered average
	// of 0 and 1 between species fractions
	rootzone = AccountStorages(SoilL1, ProotzoneL1, b) +
	  AccountStorages(SoilL2, ProotzoneL2, b) +
	  AccountStorages(SoilL3, ProotzoneL3, b) ;
}

void Budget::TotalStorage_d2H( const grid *Canopy, const grid *Canopy_d2H,
			       const grid *Snow, const grid *Snow_d2H,
			       const grid *Ponding, const grid *Ponding_d2H,
			       const grid *ChanStore, const grid *ChanStore_d2H,
			       const grid *SoilL1, const grid *SoilL1_d2H,
			       const grid *ProotzoneL1,
			       const grid *SoilL2, const grid *SoilL2_d2H,
			       const grid *ProotzoneL2,			       
			       const grid *SoilL3, const grid *SoilL3_d2H,
			       const grid *ProotzoneL3,
			       const grid *GWater, const grid *GWater_d2H,
			       const Basin *b)//, const Control *ctrl)
{
  //cout << d2HTot << " "<< d2HsoilL1 << " " << d2HsoilL2 << " " << d2HsoilL3 << " " << 
  //  d2Hgrndwater << endl;
  //cout << d2Hcanopy << " "<< d2Hsnowpack << " " << d2Hponding << " " << endl;

  canopy_d2H = AccountTrckStorages(Canopy, Canopy_d2H, b);
  snowpack_d2H = AccountTrckStorages(Snow, Snow_d2H, b);
  ponding_d2H = AccountTrckStorages(Ponding, Ponding_d2H, b);
  chan_d2H = AccountTrckStorages(ChanStore, ChanStore_d2H, b);
  soilL1_d2H = AccountTrckStorages(SoilL1, SoilL1_d2H, b);
  soilL2_d2H = AccountTrckStorages(SoilL2, SoilL2_d2H, b);
  soilL3_d2H = AccountTrckStorages(SoilL3, SoilL3_d2H, b);
  grndwater_d2H = AccountTrckStorages(GWater, GWater_d2H, b);

  // For Basind2HSummary.txt
  //cout << d2Hsnowpack << endl;
  d2Hcanopy = AccountTrckStorages2(Canopy, Canopy_d2H, b);
  d2Hsnowpack = AccountTrckStorages2(Snow, Snow_d2H, b);
  d2Hponding = AccountTrckStorages2(Ponding, Ponding_d2H, b);
  d2Hchannel = AccountTrckStorages2(ChanStore,ChanStore_d2H,b);
  d2HsoilL1 = AccountTrckStorages2(SoilL1, SoilL1_d2H, b);
  d2HsoilL2 = AccountTrckStorages2(SoilL2, SoilL2_d2H, b);
  d2HsoilL3 = AccountTrckStorages2(SoilL3, SoilL3_d2H, b);
  d2Hgrndwater = AccountTrckStorages2(GWater, GWater_d2H, b);
  //cout << d2Hsnowpack << endl;
  d2Hvadose = AccountTrckVadose(SoilL1, SoilL1_d2H,
				SoilL2, SoilL2_d2H,
				SoilL3, SoilL3_d2H,
				GWater, GWater_d2H, b);

  d2Hrootzone = AccountTrckRootZone(SoilL1, SoilL1_d2H, ProotzoneL1,
				    SoilL2, SoilL2_d2H, ProotzoneL2,
				    SoilL3, SoilL3_d2H, ProotzoneL3,
				    GWater, GWater_d2H, b);

  d2HTot = AccountTrckDomain(Canopy, Canopy_d2H,
			     Snow, Snow_d2H,
			     Ponding, Ponding_d2H,
			     SoilL1, SoilL1_d2H,
			     SoilL2, SoilL2_d2H,
			     SoilL3, SoilL3_d2H,
			     GWater, GWater_d2H, b);
 
 //cout << d2HTot << " "<< d2HsoilL1 << " " << d2HsoilL2 << " " << d2HsoilL3 << " " << 
  //  d2Hgrndwater << endl;
  //cout << d2Hcanopy << " "<< d2Hsnowpack << " " << d2Hponding << " " << endl;

}

void Budget::TotalStorage_d18O( const grid *Canopy, const grid *Canopy_d18O,
				const grid *Snow, const grid *Snow_d18O,
				const grid *Ponding, const grid *Ponding_d18O,
 			        const grid *ChanStore, const grid *ChanStore_d18O,
				const grid *SoilL1, const grid *SoilL1_d18O,
				const grid *ProotzoneL1,
				const grid *SoilL2, const grid *SoilL2_d18O,
				const grid *ProotzoneL2,			       
				const grid *SoilL3, const grid *SoilL3_d18O,
				const grid *ProotzoneL3,
				const grid *GWater, const grid *GWater_d18O,
				const Basin *b)//, const Control *ctrl)
{
  //cout << d18OTot << " "<< d18OsoilL1 << " " << d18OsoilL2 << " " << d18OsoilL3 << " " << 
  //  d18Ogrndwater << endl;
  //cout << d18Ocanopy << " "<< d18Osnowpack << " " << d18Oponding << " " << endl;

  canopy_d18O = AccountTrckStorages(Canopy, Canopy_d18O, b);
  snowpack_d18O = AccountTrckStorages(Snow, Snow_d18O, b);
  ponding_d18O = AccountTrckStorages(Ponding, Ponding_d18O, b);
  chan_d18O = AccountTrckStorages(ChanStore, ChanStore_d18O, b);
  soilL1_d18O = AccountTrckStorages(SoilL1, SoilL1_d18O, b);
  soilL2_d18O = AccountTrckStorages(SoilL2, SoilL2_d18O, b);
  soilL3_d18O = AccountTrckStorages(SoilL3, SoilL3_d18O, b);
  grndwater_d18O = AccountTrckStorages(GWater, GWater_d18O, b);

  // For Basind18OSummary.txt
  //cout << d18Osnowpack << endl;
  d18Ocanopy = AccountTrckStorages2(Canopy, Canopy_d18O, b);
  d18Osnowpack = AccountTrckStorages2(Snow, Snow_d18O, b);
  d18Oponding = AccountTrckStorages2(Ponding, Ponding_d18O, b);
  d18Ochannel = AccountTrckStorages2(ChanStore,ChanStore_d18O,b);
  d18OsoilL1 = AccountTrckStorages2(SoilL1, SoilL1_d18O, b);
  d18OsoilL2 = AccountTrckStorages2(SoilL2, SoilL2_d18O, b);
  d18OsoilL3 = AccountTrckStorages2(SoilL3, SoilL3_d18O, b);
  d18Ogrndwater = AccountTrckStorages2(GWater, GWater_d18O, b);

  //cout << d18Osnowpack << endl;
  d18Ovadose = AccountTrckVadose(SoilL1, SoilL1_d18O,
				SoilL2, SoilL2_d18O,
				SoilL3, SoilL3_d18O,
				GWater, GWater_d18O, b);

  d18Orootzone = AccountTrckRootZone(SoilL1, SoilL1_d18O, ProotzoneL1,
				     SoilL2, SoilL2_d18O, ProotzoneL2,
				     SoilL3, SoilL3_d18O, ProotzoneL3,
				     GWater, GWater_d18O, b);

  d18OTot = AccountTrckDomain(Canopy, Canopy_d18O,
			     Snow, Snow_d18O,
			     Ponding, Ponding_d18O,
			     SoilL1, SoilL1_d18O,
			     SoilL2, SoilL2_d18O,
			     SoilL3, SoilL3_d18O,
			     GWater, GWater_d18O, b);

  //cout << d18OTot << " "<< d18OsoilL1 << " " << d18OsoilL2 << " " << d18OsoilL3 << " " << 
  //  d18Ogrndwater << endl;
  //cout << d18Ocanopy << " "<< d18Osnowpack << " " << d18Oponding << " " << endl;
}

void Budget::TotalStorage_Age( const grid *Canopy, const grid *Canopy_Age,
			       const grid *Snow, const grid *Snow_Age,
			       const grid *Ponding, const grid *Ponding_Age,
			       const grid *ChanStore, const grid *ChanStore_Age,
			       const grid *SoilL1, const grid *SoilL1_Age,
			       const grid *ProotzoneL1,
			       const grid *SoilL2, const grid *SoilL2_Age,
			       const grid *ProotzoneL2,			       
			       const grid *SoilL3, const grid *SoilL3_Age,
			       const grid *ProotzoneL3,
			       const grid *GWater, const grid *GWater_Age,
			       const Basin *b)//, const Control *ctrl)
{

  //cout << AgeTot << " "<< AgesoilL1 << " " << AgesoilL2 << " " << AgesoilL3 << " " << 
  //  Agegrndwater << endl;
  //cout << Agecanopy << " "<< Agesnowpack << " " << Ageponding << " " << endl;

  canopy_Age = AccountTrckStorages(Canopy, Canopy_Age, b);
  snowpack_Age = AccountTrckStorages(Snow, Snow_Age, b);
  ponding_Age = AccountTrckStorages(Ponding, Ponding_Age, b);
  chan_Age = AccountTrckStorages(ChanStore, ChanStore_Age, b);
  soilL1_Age = AccountTrckStorages(SoilL1, SoilL1_Age, b);
  soilL2_Age = AccountTrckStorages(SoilL2, SoilL2_Age, b);
  soilL3_Age = AccountTrckStorages(SoilL3, SoilL3_Age, b);
  grndwater_Age = AccountTrckStorages(GWater, GWater_Age, b);

  // For BasinAgeSummary.txt
  //cout << Agesnowpack << endl;
  Agecanopy = AccountTrckStorages2(Canopy, Canopy_Age, b);
  Agesnowpack = AccountTrckStorages2(Snow, Snow_Age, b);
  Ageponding = AccountTrckStorages2(Ponding, Ponding_Age, b);
  Agechannel = AccountTrckStorages2(ChanStore,ChanStore_Age,b);
  AgesoilL1 = AccountTrckStorages2(SoilL1, SoilL1_Age, b);
  AgesoilL2 = AccountTrckStorages2(SoilL2, SoilL2_Age, b);
  AgesoilL3 = AccountTrckStorages2(SoilL3, SoilL3_Age, b);
  Agegrndwater = AccountTrckStorages2(GWater, GWater_Age, b);

  //cout << Agesnowpack << endl;
  Agevadose = AccountTrckVadose(SoilL1, SoilL1_Age,
				SoilL2, SoilL2_Age,
				SoilL3, SoilL3_Age,
				GWater, GWater_Age, b);

  Agerootzone = AccountTrckRootZone(SoilL1, SoilL1_Age, ProotzoneL1,
				    SoilL2, SoilL2_Age, ProotzoneL2,
				    SoilL3, SoilL3_Age, ProotzoneL3,
				    GWater, GWater_Age, b);

  AgeTot = AccountTrckDomain(Canopy, Canopy_Age,
			     Snow, Snow_Age,
			     Ponding, Ponding_Age,
			     SoilL1, SoilL1_Age,
			     SoilL2, SoilL2_Age,
			     SoilL3, SoilL3_Age,
			     GWater, GWater_Age, b);

  //cout << AgeTot << " "<< AgesoilL1 << " " << AgesoilL2 << " " << AgesoilL3 << " " << 
  //  Agegrndwater << endl;
  //cout << Agecanopy << " "<< Agesnowpack << " " << Ageponding << " " << endl;

}
