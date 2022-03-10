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
void Tracking::MixingV_sapflow(Atmosphere &atm, Basin &bsn, Control &ctrl,
			       double pTrp1, double pTrp2, double pTrp3, double transp,
			       double &iTra, double &iTraVap,
			       //double &iTraW, double &iTraVapW,    //sapflow and leaf vap
			       double p, int r, int c, int s, int iso) 
{

  double w_iTr = 0; //weighted tracer sapflow
  double w_iTrV= 0; //weighted tracer vapour 
  //  double iLeaf = 0; //Final leaf tracer

  // For two pore domain
  double d1, d2;
  double kMW_L1, kTB_L1, kMW_L2, kTB_L2;
  double dth1, dth2;
  double theta_MW1, theta_MW2;
  double theta1, theta2;
  double theta_r1, theta_r2;
  //---------------------------------------------------------------
  // TWO PORE DOMAIN
  //---------------------------------------------------------------  
  if(ctrl.sw_TPD){  //if two pore domain
    kMW_L1 = kTB_L1 = kMW_L2 = kTB_L2 = 0;
    theta_MW1 = theta_MW2 = 0;
    d1 = bsn.getSoilDepth1()->matrix[r][c];
    d2 = bsn.getSoilDepth2()->matrix[r][c];
    dth1 = transp * p * ctrl.dt * pTrp1 /d1;
    dth2 = transp * p * ctrl.dt * pTrp2 /d2;
    theta1 = bsn.getSoilMoist1()->matrix[r][c];
    theta2 = bsn.getSoilMoist2()->matrix[r][c];
    theta_r1 = bsn.getSoilMoistR1()->matrix[r][c];
    theta_r2 = bsn.getSoilMoistR2()->matrix[r][c];
    
    theta_MW1 = bsn.getMoistureMW1()->matrix[r][c];
    theta_MW2 = bsn.getMoistureMW2()->matrix[r][c];

    if(dth1 > RNDOFFERR){
      kMW_L1 = max<double>(0,(theta1 - theta_MW1) / (theta1 - theta_r1));
      kTB_L1 = max<double>(0,1 - kMW_L1);
    }
    if(dth2 > RNDOFFERR){
      kMW_L2 = max<double>(0,(theta2 - theta_MW2) / (theta2 - theta_r2));
      kTB_L2 = max<double>(0,1 - kMW_L2);
    }

    if (iso == 0){
      w_iTr = pTrp1*(kTB_L1*getd2H_TB1()->matrix[r][c] +
		     kMW_L1*getd2H_MW1()->matrix[r][c]) +
	      pTrp2*(kTB_L2*getd2H_TB2()->matrix[r][c] +
	             kMW_L2*getd2H_MW2()->matrix[r][c]) +
              pTrp3*getd2Hsoil3()->matrix[r][c];
    } else if (iso == 1){
      w_iTr = pTrp1*(kTB_L1*getd18O_TB1()->matrix[r][c] +
		     kMW_L1*getd18O_MW1()->matrix[r][c]) +
	      pTrp2*(kTB_L2*getd18O_TB2()->matrix[r][c] +
	             kMW_L2*getd18O_MW2()->matrix[r][c]) +
              pTrp3*getd18Osoil3()->matrix[r][c];
    } else if (iso == 2){
      w_iTr = pTrp1*(kTB_L1*getAge_TB1()->matrix[r][c] +
		     kMW_L1*getAge_MW1()->matrix[r][c]) +
	      pTrp2*(kTB_L2*getAge_TB2()->matrix[r][c] +
	             kMW_L2*getAge_MW2()->matrix[r][c]) +
              pTrp3*getAgesoil3()->matrix[r][c];
    }
  }
  //---------------------------------------------------------------
  // ONE PORE DOMAIN
  //---------------------------------------------------------------    
  else {          // if only one pore domain
    if (iso == 0){
      w_iTr = pTrp1*getd2Hsoil1()->matrix[r][c] +
	      pTrp2*getd2Hsoil2()->matrix[r][c] +
              pTrp3*getd2Hsoil3()->matrix[r][c];
    } else if (iso == 1){
      w_iTr = pTrp1*getd18Osoil1()->matrix[r][c] +
	      pTrp2*getd18Osoil2()->matrix[r][c] +
              pTrp3*getd18Osoil3()->matrix[r][c];
    } else if (iso == 2){
      w_iTr = pTrp1*getAgesoil1()->matrix[r][c] +
	      pTrp2*getAgesoil2()->matrix[r][c] +
              pTrp3*getAgesoil3()->matrix[r][c];
    }
  }

  /* TODO input ther Craig-Gordon model at the leaf atmosphere interface
  if(ctrl.toggle_plant_hydr == 1 && ctrl.toggle_leaf_frac){ //needed to get the pore pressure in the leaf
    if(iso == 0){ //deuterium
      Frac_Ecanopy(atm,ctrl,transp*dt,0.0,w_iTr,iLeaf,w_iTrV,0,r,c,0);
    }
    if(iso == 1){ //oxygen-18
      Frac_Ecanopy(atm,ctrl,transp*dt,0.0,w_iTr,iLeaf,w_iTrV,0,r,c,1);
    }
  } else { // set vapour to the liquid value
    iTraVap = w_iTr;
  }
  */
  
  w_iTrV = iTraVap;
  
  iTra = w_iTr;
  iTraVap = w_iTrV;
  
}


