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
 * FCdownstream.cpp
 *
 *  Created on: Apr 24, 2018
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"

void Tracking::FCdownstream(Basin &bsn, Control &ctrl,double &Qk1, double &dtdx, double &dx, int r, int c, int rr, int cc)
{
  int chan = (ctrl.sw_channel && bsn.getChannelWidth()->matrix[r][c]) ? 1 : 0 ;
  //***************************************************************************************************
  //--- Deuterium -------------------------------------------------------------------------------------
  //***************************************************************************************************
  if(ctrl.sw_2H){
    _Fd2HLattoGW->matrix[rr][cc] += bsn.getFluxGWtoLat()->matrix[r][c]*_d2Hgroundwater->matrix[r][c];
    if(ctrl.sw_deepGW)
      _Fd2HLattoDeepGW->matrix[rr][cc] += bsn.getFluxDeepGWtoLat()->matrix[r][c]*_d2HdeepGWQ->matrix[r][c];
    _Fd2HLattoSrf->matrix[rr][cc] += bsn.getFluxSrftoLat()->matrix[r][c] * _d2Hsurface->matrix[r][c];
    _Fd2HLattoChn->matrix[rr][cc] += Qk1*dtdx / dx * ( (chan == 1) ? _d2Hchan->matrix[r][c] : _d2Hsurface->matrix[r][c] );
  }

  //***************************************************************************************************  
  //--- Oxygen 18 -------------------------------------------------------------------------------------
  //***************************************************************************************************
  if(ctrl.sw_18O){
    _Fd18OLattoGW->matrix[rr][cc] += bsn.getFluxGWtoLat()->matrix[r][c] * _d18Ogroundwater->matrix[r][c];
    if(ctrl.sw_deepGW)
      _Fd18OLattoDeepGW->matrix[rr][cc] += bsn.getFluxDeepGWtoLat()->matrix[r][c] * _d18OdeepGWQ->matrix[r][c];
    _Fd18OLattoSrf->matrix[rr][cc] += bsn.getFluxSrftoLat()->matrix[r][c] * _d18Osurface->matrix[r][c];
    _Fd18OLattoChn->matrix[rr][cc] += Qk1*dtdx / dx * ( (chan == 1) ? _d18Ochan->matrix[r][c] : _d18Osurface->matrix[r][c] );
  }

  //***************************************************************************************************  
  //--- Water age -------------------------------------------------------------------------------------
  //***************************************************************************************************
  if(ctrl.sw_Age){
    _FAgeLattoGW->matrix[rr][cc] += bsn.getFluxGWtoLat()->matrix[r][c] * _Agegroundwater->matrix[r][c];
    if(ctrl.sw_deepGW)
      _FAgeLattoDeepGW->matrix[rr][cc] += bsn.getFluxDeepGWtoLat()->matrix[r][c] * _AgedeepGWQ->matrix[r][c];
    _FAgeLattoSrf->matrix[rr][cc] += bsn.getFluxSrftoLat()->matrix[r][c] * _Agesurface->matrix[r][c];
    _FAgeLattoChn->matrix[rr][cc] += Qk1*dtdx / dx * ( (chan == 1) ? _Agechan->matrix[r][c] : _Agesurface->matrix[r][c] );
  }
}


