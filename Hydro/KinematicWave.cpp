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
 * KinematicWave.cpp
 *
 *  Created on: Aug 6, 2015
 *      Author: marco
 */


#include "Basin.h"



void Basin::KinematicWave(REAL8 &Qk1, REAL8 &S,  REAL8 &Qij1,  REAL8 &qall,  REAL8 dt, UINT4 r, UINT4 c)
{
  
  REAL8 C = 0; //rhs of kinematic wave equation
  
  
  REAL8 a, n, w,  sqrtS, abQ, Qk,  fQj1i1, dfQj1i1; //kinematic wave factors
  REAL8 dtdx = dt/_dx;
  UINT4 count = 0;
  
  //kinematic wave
  
  Qij1 = _Disch_upstreamBC->matrix[r][c]; //Q at the upstream end of the channel at t+1
  Qk1 = 0;
  
  if(qall+Qij1 > 0){ //if there is water to route
    sqrtS = powl(_slope->matrix[r][c], 0.5);
    
    w = _channelwidth->matrix[r][c];
    n = _Manningn->matrix[r][c];
    a = powl(powl(w,0.67)*n/sqrtS, 0.6); //wetted perimeter is approximated with channel width
    //initial guess through solution of linear kw
    REAL8 avQ = 0.5*(Qij1);
    if (avQ==0) abQ=0;
    else abQ = a*0.6*powl(avQ, 0.6-1);
    Qk1 = ((dtdx*Qij1) + dt*qall)/(dtdx+abQ);
    
    //	if(Qk1<1e-6){
    //		Qk1 = Qij1+qall*_dx  + a * powl(Qi1j,0.6)*_dx/dt;
    
    //}
    {
      
      C =  dtdx * Qij1 + dt*qall;
      
      
      do{
	Qk=Qk1;
	fQj1i1 = dtdx*Qk+a*powl(Qk, 0.6)-C;
	dfQj1i1 = dtdx+a*0.6*powl(Qk, 0.6-1);
	Qk1 = Qk - (fQj1i1/dfQj1i1);
	if (Qk1 <=0){// if NR cannot converge then get some of the available water out and exit the loop
	  Qk1 = 0.61803*((dtdx*Qij1) + dt*qall)/(dtdx+abQ);//Qij1+qall*_dx + a * powl(Qi1j,0.6)*_dx/dt;
	  break;
	}
	count++;
      }while(fabs(fQj1i1)>0.00001 && count < MAX_ITER);
      if(count >=MAX_ITER)
	cout << "Kinematic wave solution did not converge" << endl;
    }
    
  }
  
  // double Salt =  a * powl(Qk1,0.6)*_dx; //best way to backcalcuate channel storage but produces mass balance errors when NR fails to converge
  S = std::max<double>(0.0,(Qij1+qall*_dx  - Qk1)*dt);//a*powl(Qk1, 0.6)*_dx;//(Qij1 +  a * powl(Qi1j,0.6)*_dx/dt + qall*_dx - Qk1)*dt;
  
  //end kinematic wave
  
}


