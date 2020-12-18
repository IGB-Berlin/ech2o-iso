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
 * MixingTPD_postET.cpp
 *
 *  Created on: Apr 19, 2018
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"

void Tracking::MixingTPD_postET(Basin &bsn, Control &ctrl,
				double &dth1, double &dth2,
				double &kTB_L1, double &kTB_L2,
				double &kMW_L1, double &kMW_L2,
				int r, int c)
{  

  double theta_MW1 = bsn.getMoistureMW1()->matrix[r][c];
  double theta_MW2 = bsn.getMoistureMW2()->matrix[r][c];
  // Tightly-bound "content" after root water uptake
  double TB1_new = min<double>(bsn.getSoilMoist1()->matrix[r][c],theta_MW1) - kTB_L1*dth1;
  double TB2_new = min<double>(bsn.getSoilMoist2()->matrix[r][c],theta_MW2) - kTB_L2*dth2;
  // Mobile water "content" after root water uptake
  double MW1_new = max<double>(0,bsn.getSoilMoist1()->matrix[r][c] - theta_MW1 - kMW_L1*dth1);
  double MW2_new = max<double>(0,bsn.getSoilMoist2()->matrix[r][c] - theta_MW2 - kMW_L2*dth2);
    
  // Layer 1 ----------------------------------------------------------------------------------
  if(kTB_L1 > RNDOFFERR and MW1_new > RNDOFFERR) {

    if(ctrl.sw_2H)
      _d2H_TB1->matrix[r][c] = InputMix(TB1_new, _d2H_TB1->matrix[r][c],
					min<double>(theta_MW1 - TB1_new,MW1_new), 
					_d2H_MW1->matrix[r][c]) ;
    if(ctrl.sw_18O)
      _d18O_TB1->matrix[r][c] = InputMix(TB1_new, _d18O_TB1->matrix[r][c],
					 min<double>(theta_MW1-TB1_new,MW1_new), 
					 _d18O_MW1->matrix[r][c]) ;

    if(ctrl.sw_Age)
      _Age_TB1->matrix[r][c] = InputMix(TB1_new, _Age_TB1->matrix[r][c],
					min<double>(theta_MW1-TB1_new,MW1_new), 
					_Age_MW1->matrix[r][c]) ;
  }
 
  // Layer 2 ----------------------------------------------------------------------------------
  if(kTB_L2 > RNDOFFERR and MW2_new > RNDOFFERR) {

    if(ctrl.sw_2H)
      _d2H_TB2->matrix[r][c] = InputMix(TB2_new, _d2H_TB2->matrix[r][c],
					min<double>(theta_MW2-TB2_new,MW2_new), 
					_d2H_MW2->matrix[r][c]) ;

    if(ctrl.sw_18O)
      _d18O_TB2->matrix[r][c] = InputMix(TB2_new, _d18O_TB2->matrix[r][c],
					 min<double>(theta_MW2-TB2_new,MW2_new), 
					 _d18O_MW2->matrix[r][c]) ;

    if(ctrl.sw_Age)
      _Age_TB2->matrix[r][c] = InputMix(TB2_new, _Age_TB2->matrix[r][c],
					min<double>(theta_MW2-TB2_new,MW2_new), 
					_Age_MW2->matrix[r][c]) ;
  }
}
  
