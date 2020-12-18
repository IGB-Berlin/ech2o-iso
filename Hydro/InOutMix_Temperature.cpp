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

double Basin::InOutMix_Temperature(double hold, //old channel volume
				 double Told, //old channel temperature
				 double qin, //incoming volume
				 double Tin, //incoming temperature
				 double qout, //outgoing volume
				 int mode)
{
  double Tnew = 0;
  double hsum = 0;

  if(mode==0){
    Tnew = hold + qin > RNDOFFERR ? (Told*hold + Tin*qin) / (hold + qin) : Told;

  }else if(mode==1){
    hsum = 0.5 * (hold + qin + std::max<double>(0,hold - qout));
    Tnew = 0.5 * hsum + qin > RNDOFFERR ? (Told*(hsum - 0.5*qin) + Tin*qin) / (hsum + 0.5*qin) : Told;

  }

  return Tnew;
}
