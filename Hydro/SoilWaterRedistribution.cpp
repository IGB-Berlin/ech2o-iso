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
 * SoilWaterRedistribution.cpp
 *
 *  Created on: Jul 22, 2015
 *      Author: marco
 */

#define ARMA_NO_DEBUG //disables armadillo bound checks for speed optimization
//#include <armadillo>
#include"Basin.h"

//using namespace arma;

void Basin::SoilWaterRedistribution(Control &ctrl, const double &F, double &theta1,
				    double &theta2, double &theta3, double &pond,
				    double &gw, double &leak,  double dt,
				    int r, int c) {


  //double K1, K12, K2, K23, K3;
  double theta_r1 = _theta_rL1->matrix[r][c];
  double theta_r2 = _theta_rL2->matrix[r][c];
  double theta_r3 = _theta_rL3->matrix[r][c];
  double thetafc1 = _fieldcapL1->matrix[r][c];
  double thetafc2 = _fieldcapL2->matrix[r][c];
  double thetafc3 = _fieldcapL3->matrix[r][c];
  double poros1 = _porosityL1->matrix[r][c];
  double poros2 = _porosityL2->matrix[r][c];
  double poros3 = _porosityL3->matrix[r][c];
  double KvKh1 = _KvKsL1->matrix[r][c];
  double KvKh2 = _KvKsL2->matrix[r][c];
  double KvKh3 = _KvKsL3->matrix[r][c];
  double Ks1 = _KsatL1->matrix[r][c] * KvKh1;
  double Ks2 = _KsatL2->matrix[r][c] * KvKh2;
  double Ks3 = _KsatL3->matrix[r][c] * KvKh3;
  double L = _bedrock_leak->matrix[r][c];

  //depth of soil layers
  double depth = _soildepth->matrix[r][c];
  double d1 = _depth_layer1->matrix[r][c];
  double d2 = _depth_layer2->matrix[r][c];
  double d3 = depth - d1 - d2;

  double  x[3] = {};

  double L1 = theta1*d1;
  double L2 = theta2*d2;
  double L3 = theta3*d3;

  x[0] = L1;
  x[1] = L2;
  x[2] = L3;

  double a1 = dt*Ks1/(poros1-theta_r1);
  double a2 = dt*Ks2/(poros2-theta_r2);
  double a3 = dt*Ks3/(poros3-theta_r3);

  // == Gravitational drainage -------------------------------------------------

  // -- First layer 1 + updating the layer below
  if(x[0]/d1 > thetafc1){
    x[0] =( L1 + a1* theta_r1) / (1 + a1/d1);
    //check if too much drainage
    if(x[0]/d1 < thetafc1)
      x[0] = thetafc1 * d1;

    L2 += L1 - x[0];

    // Tracking
    _FluxPercolL2->matrix[r][c] += max<REAL8>(0,L1 - x[0]);
    if(ctrl.sw_trck)
      _FluxL1toL2->matrix[r][c] += max<REAL8>(0,L1 - x[0]);

    x[1]=L2;
  }

  // -- Second layer + updating the layer below
  if(L2/d2 > thetafc2){
    x[1] =(L2 + a2* theta_r2) / (1 + a2/d2);
    //check if too much drainage
    if(x[1]/d2 < thetafc2)
      x[1] = thetafc2 * d2;

    L3 += L2 - x[1];

    _FluxPercolL3->matrix[r][c] += max<REAL8>(0,L2 - x[1]);
    // Tracking
    if(ctrl.sw_trck)
      _FluxL2toL3->matrix[r][c] += max<REAL8>(0,L2 - x[1]);

    x[2]=L3;
  }

  // -- Third layer 
  if(L3/d3 > thetafc3){
    x[2] =(L3 + L*a3* theta_r3) / (1 + L*a3/d3);
    //check if too much drainage
    if(x[2]/d3 < thetafc3)
      x[2] = thetafc3 * d3;

  }

  theta1 = x[0]/d1;
  theta2 = x[1]/d2;
  theta3 = x[2]/d3;

  // -- Check if it creates over-filling -> cascading adjustment
  //pond -=F;
  // -- L3
  if(theta3 > poros3){
    theta2 += (theta3 - poros3) * d3/d2;

    // Tracking : remove the excess
    _FluxPercolL3->matrix[r][c] -= (theta3 - poros3) * d3 ;
    if(ctrl.sw_trck)
      _FluxL2toL3->matrix[r][c] -= (theta3 - poros3) * d3 ;

    theta3 = poros3;
  }
  // -- L2
  if(theta2 > poros2){
    theta1 += (theta2 - poros2) * d2/d1;

    // Tracking
    _FluxPercolL2->matrix[r][c] -= (theta2 - poros2) * d2;
    if(ctrl.sw_trck)
      _FluxL1toL2->matrix[r][c] -= (theta2 - poros2) * d2;

    theta2 = poros2;
  }
  // -- L1
  if(theta1 > poros1){
    pond += -(poros1 - theta1) * d1;

    // Tracking
    _FluxInfilt->matrix[r][c] -= (theta1 - poros1) * d1;
    if(ctrl.sw_trck)
      _FluxSrftoL1->matrix[r][c] -= (theta1 - poros1) * d1;

    theta1 = poros1;
  }

  // -- Gravitational water (L3) and bedrock leakage
  // (matrices are updated in SolveSurfaceFluxes.cpp)
  _FluxRecharge->matrix[r][c] += max<double>(0,max<double>(0,(theta3 - thetafc3) * d3) - gw);
  gw = max<double>(0,(theta3 - thetafc3) * d3);
  leak = std::max<double>(0,(L3 - x[2])/dt);

  if(ctrl.sw_trck)
    _FluxLeak->matrix[r][c] = (L3 - x[2]);

}

