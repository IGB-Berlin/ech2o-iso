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
 * Created on: Feb 2, 2020
 *     Author: Aaron Smith
 */
#include "Basin.h"
#include "Forest.h"
#include <armadillo>
void Tracking::MixingV_canopy(Atmosphere &atm, Basin &bsn, Control &ctrl,
			      double SR_frac,     // snow/rain fraction
			      double IntercWater, // intercepted water (after precip addition and after evap)
			      double evap,        // interception evaporation [m/s]
			      double &icanopy,    // species canopy tracer
			      double &ievapV,     // species evaporation vapour
			      double &ithrough,   // species throughfall
			      double D,           // throughfall (outgoing)   [m/s]
			      double iPrec,       // incoming precipitation tracer
			      double p,           // species proportion
			      int r,
			      int c,
			      int s,              // vegetation species
			      int iso) 
{
  // -----------------------------------------------
  // If SR_frac == 1 (only snow) 
  //     there is no mixing in the canopy
  //     throughfall  == rain
  //     equal mixing within the canopy
  // Else
  //     mix according to the options
  // -----------------------------------------------
  
  // Initialize the values for mixing
  //  double h_canold = IntercWater - (DelCanStor + evap*ctrl.dt); // initial canopy storage before incoming precip and evap
  double Pp = atm.getPrecipitation()->matrix[r][c];
  double h_canold = IntercWater + (D + evap - Pp)*ctrl.dt;
  double qin;  
  double qout; 
  int mixmod;
  double i_mixcan = 0; //Set the mixed canopy tracer to be updated in TracerMixing
  double beta = 0;
  int issoil = 0;
  // --------------------------------------------------------------------------------------------
  // MIXING OF PRE-CANOPY WATER WITH RAIN/SNOW
  // --------------------------------------------------------------------------------------------
  if(ctrl.sw_int_mix){
    mixmod = 2;
    qin = (Pp*ctrl.dt) < RNDOFFERR ? 0.0 : (Pp*ctrl.dt);     // total incoming water to mix
    qout = (D*ctrl.dt) < RNDOFFERR ? 0.0 : (D*ctrl.dt);
    if(SR_frac == 1){
      TracerMixing(bsn,ctrl,h_canold,icanopy,i_mixcan,qin,iPrec,0.0,ithrough,1.0,mixmod,r,c);
      ithrough = iPrec;
    } else {
      TracerMixing(bsn,ctrl, h_canold, icanopy,i_mixcan,qin,iPrec,qout,ithrough,1.0,mixmod,r,c);
                          // voli    , Tvoli  ,Tvoln   ,Vin,Tin  ,Vout,Tout    ,Sat
    } //end check if snow
  } else {
    qin = (Pp - D) * ctrl.dt < RNDOFFERR ? 0.0 : (Pp - D) * ctrl.dt;
    qout = 0;
    mixmod = 0;

    if(IntercWater > RNDOFFERR && h_canold > RNDOFFERR){ //only if new interception and old interception are above 0m
      TracerMixing(bsn,ctrl,h_canold,icanopy,i_mixcan,qin,iPrec,qout,ithrough,1.0,mixmod,r,c);
      //                    voli    ,Tvoli  , Tvoln  ,Vin,Tin  ,Vout,Tout    ,Sat
      ithrough = iPrec;
      
    } else {
      ithrough = iPrec;
      i_mixcan = iPrec;
    }

  }

  // --------------------------------------------------------------------------------------------

  // --------------------------------------------------------------------------------------------
  // FRACTIONATE WATER IN THE CANOPY IF THERE IS EVAPORATION
  // --------------------------------------------------------------------------------------------  
  double V_new = IntercWater;
  double V_old = V_new + evap*ctrl.dt;
  if(iso < 2){
    if((V_new/V_old)<0.999 and V_new > RNDOFFERR and ctrl.sw_frac){
      Frac_Estorage(atm,bsn,ctrl,
		   V_old,       // Water storage before evaporation
		   V_new,                // Water storage after evaporation
		   i_mixcan,             // tracer in storage before evaporation
		   icanopy,              // tracer in storage after evaporation
		   ievapV,               // tracer in vapour
		   beta,
		   atm.getTemperature()->matrix[r][c],
		   issoil,
		   r,c,iso);
    } else {
      icanopy = i_mixcan;
      ievapV  = i_mixcan;
    }
  }

  if(iso == 2){ //Water ages
    icanopy = i_mixcan;
    ievapV  = i_mixcan;
  }
  
}


