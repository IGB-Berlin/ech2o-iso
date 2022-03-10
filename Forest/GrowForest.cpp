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
 * CalculateGPP.cpp
 *
 *  Created on: Jun 20, 2010
 *      Author: marco
 */

#include "Forest.h"

int Forest::GrowForest(Basin &bas, const Atmosphere &atm, const Control &ctrl) {

  UINT4 r, c;
  REAL8 dt;
  REAL8 alpha = 0; //canopy quantum efficiency
  REAL8 beta = 0; //canopy water efficiency
  REAL8 par = 0;
  REAL8 forestAge = 0;
  REAL8 airTemp, optTemp, maxTemp, minTemp;
  REAL8 Wc, Wp, UsableTheta, Wr, gsmax;
  REAL8 fa, ft, fw;
  REAL8 theta_wp;
  REAL8 fc1, fc2, fc3;
  REAL8 theta, theta2, theta3;
  REAL8 froot1, froot2, froot3;
  REAL8 E, lai, BeerK ;
  unsigned int j;

  dt = ctrl.dt;

  for (j = 0; j < _Nsp - 1; j++) //grow forest up to Nsp -1 because Nsp is bare soil

#pragma omp parallel for default(none)					\
  private( r, c, alpha, beta, par, E, lai,forestAge,			\
	   airTemp, optTemp, maxTemp, minTemp, Wc, Wp, gsmax,		\
	   UsableTheta, Wr, fa, ft ,fw, fc1, fc2, fc3, BeerK,		\
	   theta_wp, theta, theta2, theta3, froot1, froot2, froot3)	\
  shared(j,bas, atm, ctrl,dt)

    for (unsigned int k = 0; k < _vSortedGrid.cells.size(); k++) {
      r = _vSortedGrid.cells[k].row;
      c = _vSortedGrid.cells[k].col;

      if (_species[j]._fraction->matrix[r][c] < RNDOFFERR) //if there are no trees of teh species there is nothign to grow
	continue;

      alpha = _species[j].alpha;
      beta = _species[j].beta;
      par = atm.getIncomingShortWave()->matrix[r][c] * dt * 0.47; //dt to convert watts to joules. Assume 47% of incoming solar radiation is par
      E = _species[j]._Transpiration->matrix[r][c] * dt; //total amount of transpiration for the time period (m)
      lai = _species[j]._LAI->matrix[r][c];
      forestAge = _species[j]._AGE->matrix[r][c];
      airTemp = atm.getTemperature()->matrix[r][c];
      optTemp = _species[j].TempOpt;
      maxTemp = _species[j].TempMax;
      minTemp = _species[j].TempMin;
      Wc = bas.getParamWc()->matrix[r][c];
      Wp = bas.getParamWp()->matrix[r][c];
      gsmax = _species[j].gsmax;
      theta_wp = _species[j].WiltingPoint;
      theta = bas.getSoilMoist1()->matrix[r][c]; // moisture in soil L1 at time t
      theta2 = bas.getSoilMoist2()->matrix[r][c];
      theta3 = bas.getSoilMoist3()->matrix[r][c];
      froot1 = _species[j]._rootfrac1->matrix[r][c];
      froot2 = _species[j]._rootfrac2->matrix[r][c];  
      froot3 = 1-froot1-froot2;
      BeerK = _species[j].KBeers;
      fc1 = bas.getFieldCapacityL1()->matrix[r][c];
      fc2 = bas.getFieldCapacityL2()->matrix[r][c];
      fc3 = bas.getFieldCapacityL3()->matrix[r][c];

      UsableTheta = max<REAL8>(0,min<REAL8>(1,(theta-theta_wp)/(fc1-theta_wp)))*froot1 +
	max<REAL8>(0,min<REAL8>(1,(theta2-theta_wp)/(fc2-theta_wp)))*froot2 +
	max<REAL8>(0,min<REAL8>(1,(theta3-theta_wp)/(fc3-theta_wp)))*froot3;

      Wr = UsableTheta; //TODO: URGENT: improve competition for water. no competition now*/

      fa = Calculate_fa(_species[j].MaxAge, forestAge);
      ft = Calculate_ft(airTemp, maxTemp, minTemp, optTemp);
      fw = Calculate_fw(_species[j]._CanopyConductance->matrix[r][c],gsmax, Wr, Wc, Wp);
      
      _species[j]._GPP->matrix[r][c] = sqrtl(alpha * par * beta * E) * fa* ft; // * fw;
      _species[j]._NPP->matrix[r][c] = _species[j]._GPP->matrix[r][c] * _species[j].GPP2NPP;
      
      if (_species[j].vegtype == 2){
          
	ft = expl(-BeerK * lai) ;
	fw = Wr;
          
      }

      // Dynamic allocation:
      // - veg_dyn = 0 : none if 
      // - veg_dyn = 1 : calculated dynamically
      // - veg_dyn = 2 : LAI forced from input times series (update from ech2o.cpp)
      if(ctrl.toggle_veg_dyn == 1){
	if (_species[j].vegtype == 1)
	  GrowGrass(j, r, c, dt);
        else	  
	  GrowTrees(j, r, c, dt, fa, ft, fw,atm.getMinTemperature()->matrix[r][c], UsableTheta);
      }
    }

  return EXIT_SUCCESS;

}
