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
 * MapGetters.cpp
 *
 *  Created on: Aug 19, 2010
 *      Author: Marco.Maneta
 */

#include "Basin.h"

UINT4 Basin::getNumSpecies() const {
  return fForest->getNumSpecies();
}

grid *Basin::getVegetFrac(UINT4 n) const {
		return fForest->getPropSpeciesMap(n);
}

grid *Basin::getLAI(UINT4 n) const {
		return fForest->getLAISpeciesMap(n);
}

grid *Basin::getStemDensity(UINT4 n) const {
		return fForest->getStemDensSpeciesMap(n);
}

grid *Basin::getStandAge(UINT4 n) const {
	return fForest->getAgeSpeciesMap(n);
}

grid *Basin::getCanopyCond(UINT4 n) const {
	return fForest->getCanopyConductSpeciesMap(n);
}

grid *Basin::getGPP(UINT4 n) const {
	return fForest->getGPPSpeciesMap(n);
}

grid *Basin::getNPP(UINT4 n) const {
	return fForest->getNPPSpeciesMap(n);
}

grid *Basin::getBasalArea(UINT4 n) const {
	return fForest->getBasalAreaSpeciesMap(n);
}

grid *Basin::getTreeHeight(UINT4 n) const {
	return fForest->getTreeHeightSpeciesMap(n);
}

grid *Basin::getRootMass(UINT4 n) const {
	return fForest->getRootMassSpeciesMap(n);
}

grid *Basin::getCanopyTemp(UINT4 n) const {
	return fForest->getCanopyTempSpeciesMap(n);
}

grid *Basin::getCanopyNetRad(UINT4 n) const {
	return fForest->getCanopyNetRadSpeciesMap(n);
}

grid *Basin::getCanopyLatHeatE(UINT4 n) const {
	return fForest->getCanopyLatHeatESpeciesMap(n);
}

grid *Basin::getCanopyLatHeatT(UINT4 n) const {
	return fForest->getCanopyLatHeatTSpeciesMap(n);
}

grid *Basin::getCanopySensHeat(UINT4 n) const {
	return fForest->getCanopySensHeatSpeciesMap(n);
}

grid *Basin::getCanopyWaterStor(UINT4 n) const {
	return fForest->getCanopyWaterStorSpeciesMap(n);
}

grid *Basin::getRootUptakeL1(UINT4 n) const {
  return fForest->getRootUptakeL1Map(n);
}

grid *Basin::getRootUptakeL2(UINT4 n) const {
  return fForest->getRootUptakeL2Map(n);
}

grid *Basin::getRootUptakeL3(UINT4 n) const {
  return fForest->getRootUptakeL3Map(n);
}

grid *Basin::getETspecies(UINT4 n) const {
	return fForest->getETSpeciesMap(n);
}

grid *Basin::getTranspiration(UINT4 n) const {
	return fForest->getTranspirationSpeciesMap(n);
}

grid *Basin::getTranspirationFlux(UINT4 n) const {
	return fForest->getTranspirationFluxSpeciesMap(n);
}

grid *Basin::getEinterception(UINT4 n) const {
	return fForest->getEinterceptionSpeciesMap(n);
}
grid *Basin::getEsoil(UINT4 n) const {
	return fForest->getEsoilSpeciesMap(n);
}

grid *Basin::getSoilWaterPotential(UINT4 n) const {
	return fForest->getSoilWaterPotSpeciesMap(n);
}

grid *Basin::getLeafWaterPotential(UINT4 n) const {
	return fForest->getLeafWaterPotSpeciesMap(n);
}

grid *Basin::getSapVelocity(UINT4 n) const {
	return fForest->getSapVelocitySpeciesMap(n);
}

grid *Basin::getRootFrac1(UINT4 n) const {
  return fForest->getRootFrac1(n);
}

grid *Basin::getRootFrac2(UINT4 n) const {
  return fForest->getRootFrac2(n);
}

// -- Tracking: calls tracking data contained in forest class
// 2H
grid *Basin::getd2Hcanopy(UINT4 n) const {
  return fForest->getd2Hcanopy(n);
}
grid *Basin::getd2Hthroughfall(UINT4 n) const {
  return fForest->getd2Hthroughfall(n);
}
grid *Basin::getd2HevapI(UINT4 n) const {
  return fForest->getd2HevapI(n);
}
grid *Basin::getd2HevapT(UINT4 n) const {
  return fForest->getd2HevapT(n);
}
grid *Basin::getd2HevapS(UINT4 n) const {
  return fForest->getd2HevapS(n);
}
grid *Basin::getd2HevapI_Vap(UINT4 n) const {
  return fForest->getd2HevapI_Vap(n);
}
grid *Basin::getd2HevapT_Vap(UINT4 n) const {
  return fForest->getd2HevapT_Vap(n);
}
// 18O
grid *Basin::getd18Ocanopy(UINT4 n) const {
  return fForest->getd18Ocanopy(n);
}
grid *Basin::getd18Othroughfall(UINT4 n) const {
  return fForest->getd18Othroughfall(n);
}
grid *Basin::getd18OevapI(UINT4 n) const {
  return fForest->getd18OevapI(n);
}
grid *Basin::getd18OevapT(UINT4 n) const {
  return fForest->getd18OevapT(n);
}
grid *Basin::getd18OevapS(UINT4 n) const {
  return fForest->getd18OevapS(n);
}
grid *Basin::getd18OevapI_Vap(UINT4 n) const {
  return fForest->getd18OevapI_Vap(n);
}
grid *Basin::getd18OevapT_Vap(UINT4 n) const {
  return fForest->getd18OevapT_Vap(n);
}
// Age
grid *Basin::getAgecanopy(UINT4 n) const {
	return fForest->getAgecanopy(n);
}
grid *Basin::getAgethroughfall(UINT4 n) const {
  return fForest->getAgethroughfall(n);
}
grid *Basin::getAgeevapI(UINT4 n) const {
	return fForest->getAgeevapI(n);
}
grid *Basin::getAgeevapT(UINT4 n) const {
	return fForest->getAgeevapT(n);
}
grid *Basin::getAgeevapS(UINT4 n) const {
	return fForest->getAgeevapS(n);
}

void Basin::setAgecanopy(UINT4 n, UINT4 r, UINT4 c, REAL8 value) const {
  fForest->setAgecanopy(n, r, c, value);
}
void Basin::setAgeevapS(UINT4 n, UINT4 r, UINT4 c, REAL8 value) const {
  fForest->setAgeevapS(n, r, c, value);
}
void Basin::setAgeevapI(UINT4 n, UINT4 r, UINT4 c, REAL8 value) const {
  fForest->setAgeevapI(n, r, c, value);
}
void Basin::setAgeevapT(UINT4 n, UINT4 r, UINT4 c, REAL8 value) const {
  fForest->setAgeevapT(n, r, c, value);
}
