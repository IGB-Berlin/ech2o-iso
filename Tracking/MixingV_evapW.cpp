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
 *    Aaron Smith
 *******************************************************************************/
/*
 * MixingV_evapW.cpp
 *
 *  Created on: Jun 21, 2020
 *      Author: Aaron Smith
 */

#include "Basin.h"

void Tracking::MixingV_evapW(Atmosphere &atm, Basin &bsn, Control &ctrl, 
			     double Tw, double evap, double chan_store,
			     int r, int c)
{

  double EVAP   = evap * ctrl.dt; 	  	// p-weighted soil evap
  double V_new = chan_store; 			// useful effective volume
  double V_old = chan_store + EVAP;
  double beta = 0;
  int issoil = 0;
  
  
  //no fractionation if water storage is too low
  if(ctrl.sw_frac and ctrl.sw_chan_frac and V_new/V_old <0.999 and V_new/V_old > 0.5 and V_new > 0.015) {
    if(ctrl.sw_2H)
      Frac_Estorage(atm, bsn, ctrl, V_old, V_new,_d2Hchan->matrix[r][c],
		    _d2Hchan->matrix[r][c],_d2HevapC_sum->matrix[r][c],
		    beta,Tw,issoil, r, c, 0);
     if(ctrl.sw_18O)
      Frac_Estorage(atm, bsn, ctrl, V_old, V_new,_d18Ochan->matrix[r][c],
		    _d18Ochan->matrix[r][c],_d18OevapC_sum->matrix[r][c],
		    beta,Tw,issoil, r, c, 1);

   } else { // If no fractionation (or negligible E), only update E
     if(ctrl.sw_2H)
       _d2HevapC_sum->matrix[r][c] = _d2Hchan->matrix[r][c]; 
     if(ctrl.sw_18O)
       _d18OevapC_sum->matrix[r][c] = _d18Ochan->matrix[r][c]; 
   }

   if(ctrl.sw_Age)
     _AgeevapC_sum->matrix[r][c] = _Agechan->matrix[r][c]; 

// -----------------------------------------------------------------------------------
}

