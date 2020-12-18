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
 * UpdateSoilMoist.cpp
 *
 *  Created on: Mar 8, 2010
 *      Author: Marco Maneta
 */
#include "Basin.h"

void Basin::ChannelEvaporation(REAL8 LE, //input latent heat
			       REAL8 lambda,
			       REAL8 rs,// input the potential exfiltration capacity
			       REAL8 &evap,
			       REAL8 &chan_store,
			       REAL8 dt, //time step
			       UINT4 r,
			       UINT4 c)
{
  REAL8 EVAP;
  evap = 0;
  EVAP = 0;
  REAL8 le = lambda;

  if(LE>0){
    evap = std::min<REAL8>(1/rs , LE/(rho_w*le));
    EVAP = evap * dt;
    if (chan_store < EVAP){ //makes sure we are not evaporating more that it is available
      EVAP = chan_store;
      evap = EVAP / dt;
    }
  }

  chan_store -= EVAP;

  _chan_store->matrix[r][c] = chan_store;
  _chan_evap->matrix[r][c] = evap;

}
