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
 * CalcSoilMoistureProfile.cpp
 *
 *  Created on: Jul 26, 2011
 *      Author: Marco.Maneta
 *
 *      This function derives a soil moisture profile from the average soil moisture assuming the soil is at hydrostatic equilibrium (except when it is raining).
 *      Hydrostatic equilibrium depends on soil characteristics as per the Brooks and Corey Formula. Average soil moisture is calculated by integrating the hydrostatic
 *      moisture profile and solving for the equivalent depth to saturation (H) that corresponds to the calculated average soil moisture content
 *
 *      Returns by reference the equivalent depth to saturation
 *
 */

#include"Basin.h"

void Basin::CalcSoilMoistureProfile(Atmosphere &atm, Control &ctrl, REAL8 theta,
				    UINT4 row, UINT4 col) {

  REAL8 d = _soildepth->matrix[row][col];

  REAL8 poros0 = _porosity0->matrix[row][col];
  REAL8 kp, poros;
	
  REAL8 theta_r;
  REAL8 psi = _psi_ae->matrix[row][col];
  REAL8 lambda = 1 / _BClambda->matrix[row][col];
  REAL8 psilam = powl(psi, lambda);
  REAL8 H = _WaterTableDepth->matrix[row][col] + 0.1; //equivalent depth to saturation (m). Initial estimate is the old value for H + 0.1 to avoid initial value of 0
  REAL8 H1 = H; //updated equivalent depth to saturation (m)
  REAL8 u = d; // depth of the unsaturated layer needed to integrate the piecewise function
  REAL8 fH, dfH;
  UINT4 k = 0;

  if(ctrl.sw_expPoros){
    kp = _kporos->matrix[row][col];
    poros = kp*poros0/d*(1-expl(-d/kp));
    theta_r = 1/d * (_depth_layer1->matrix[row][col]*_theta_rL1->matrix[row][col] +
		     _depth_layer2->matrix[row][col]*_theta_rL2->matrix[row][col] +
		     (d-_depth_layer1->matrix[row][col]-_depth_layer2->matrix[row][col]) *
		     _theta_rL3->matrix[row][col]);
  } else {
    poros = poros0;
    theta_r = _theta_rL1->matrix[row][col];
  }
	
  //TODO: prepare case of lambda = 1, which results in a different integrated function for soil moisture
  //for the time being, if lambda =1 change it to 0.95 (lambda shouldnt be that high anyway)

  if (lambda == 1)
    lambda = 0.95;

  //if we have a fully saturated soil, then there is no point in calculating any soil moisture profile

  if (((poros - theta) * d) < (psi + H)) {
    REAL8 theta1 = _soilmoist1->matrix[row][col];
    REAL8 theta2 = _soilmoist2->matrix[row][col];
    REAL8 theta3 = _soilmoist3->matrix[row][col];
    REAL8 poros1 = _porosityL1->matrix[row][col];
    REAL8 poros2 = _porosityL2->matrix[row][col];
    REAL8 poros3 = _porosityL3->matrix[row][col];
    REAL8 d1 = _depth_layer1->matrix[row][col];
    REAL8 d2 = _depth_layer2->matrix[row][col];
    REAL8 d3 = _soildepth->matrix[row][col] - d1 - d2;
    REAL8 S1 = 1 - std::max<double>((poros1 - theta1)/poros1,0); 
    REAL8 S2 = 1 - std::max<double>((poros2 - theta2)/poros2,0);
    REAL8 S3 = 1 - std::max<double>((poros3 - theta3)/poros3,0);

    _WaterTableDepth->matrix[row][col] = (poros - (theta1*S1*d1 + theta2*S2*d2 + theta3*S3*d3)/(d1 + d2 + d3))*(d1 + d2 + d3)/poros;
    /* THIS IS THE ORIGINAL -- PROBLEMS WITH 
    _WaterTableDepth->matrix[row][col] = (poros - theta) * d;
    //_soilmoist10cm->matrix[row][col] = theta;
    */
    return;
  }

  do { //uses N-R to calculate the corresponding depth to the equivalent water table for the given average soil moisture

    H = H1;
    if (H == d)
      H += 0.01; //in the unlikely case H equals d, make it 1 cm deeper to avoid division by zero
    if (H > 1000) { //if the depth water table is larger than some large arbitrary depth (m), assume soil moisture at the topsoil equals average soil moisture and return
      //_soilmoist10cm->matrix[row][col] = theta;
      _WaterTableDepth->matrix[row][col] = H1;
      return;
    }
    u = ((H - d) > psi) ? d : H - psi;

    fH = (1 / d)
      * (((((poros - theta_r) * psilam
	    * ((H - u) * pow(H, lambda) - H * pow(H - u, lambda)))
	   / ((lambda - 1) * pow(H - u, lambda) * pow(H, lambda)))
	  + theta_r * u) + poros * (d - u)) - theta;
    dfH = ((H - d) > psi) ?
      (1 / d)
      * -((poros - theta_r)
	  * (pow(psi / (H - u), lambda)
	     - pow(psi / (H), lambda))) :
      -poros / d;

    H1 = H - fH / dfH;

    k++;

  } while (fabs(H - H1) > 0.01 && k < MAX_ITER);

  if (k >= MAX_ITER)
    std::cout
      << "WARNING: non-convergence soil moisture profile at cell row: "
      << row << " col: " << col << " closure err: " << (H1 - H)
      << endl;

  //  cout << "Unsaturated --> " << "poros: " << poros << " |theta: " << theta << " |d: " << d << " |u: " << u << endl;
  //  cout << "                lambda" << lambda << " |psi: " << psi << " |H: " << H1 << " |row: " << row << " |col: " << col << endl;

  //calculates average soil moisture of first 10 cm of soil given teh calculated depth to the water table
  d = 0.1;
  u = ((H1 - d) > psi) ? d : H1 - psi; //

  /*_soilmoist10cm->matrix[row][col] = (1 / d)
   * (((((poros - theta_r) * psilam
   * ((H1 - u) * powl(H1, lambda) - H1 * powl(H1 - u, lambda)))
   / ((lambda - 1) * powl(H1 - u, lambda) * powl(H1, lambda)))
   + theta_r * u) + poros * (u - d)); */

  //stores the calculated depth to saturation
  _WaterTableDepth->matrix[row][col] = H1;

}
