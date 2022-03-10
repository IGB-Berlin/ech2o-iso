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

void Tracking::TracerMixing(Basin &bsn, Control &ctrl, 
			         double hold, //old volume
				 double iold, //old volume tracer
			         double &inew, //new volume tracer
				 double qin,  //incoming volume
				 double iin,  //incoming tracer
				 double qout, //outgoing volume
			         double &iout,//output tracer
			         double S,    // storage saturation
 			         int mixmod,
			         int r,
			         int c) 
{

  double hsum = 0;
  double f = 0;
  double intexp  = 0;
  double intfrac = 0;

  // MIX MOD 0
  if(mixmod==0){
     inew = hold + qin > RNDOFFERR ? (iold*hold + iin*qin)/ (hold + qin) : iold;
     iout = inew;
  }

  // MIX MOD 1  
  if(mixmod ==1) {
     hsum = 0.5*(hold+qin+std::max<double>(0,hold-qout));
     inew = 0.5*hsum+qin > RNDOFFERR ? 
	((iold+1000)*(hsum-0.5*qin) + (iin+1000)*qin)/ (hsum + 0.5*qin) -1000 : iold;
     iout = inew;
  }
  
  // MIX MOD 2
  if(mixmod == 2) {  // Gibson 2002 complete mixing approach -- corrected outflow of each soil layer
    if(qin > RNDOFFERR && hold > RNDOFFERR){// Check in there is inflow and initial volume
      // -------------TIME DEPENDENT APPROACH -------------------------
      if(abs(qin - qout) < RNDOFFERR){ 
	inew = iin - (iin - iold) * exp(-qin/hold);
        intexp = -(hold/qin)*exp(-qin/hold)-(-hold/qin); // definite integral of flux over the time-step
	iout = iin - (iin - iold) * intexp; // actual averaged outflow tracer
      // -------------FRACTION DEPENDENT APPROACH ---------------------
      } else {                              
	f = (hold + qin - qout)/hold;
	inew = iin - (iin - iold) * powl(f,(-1/(1-(qout/qin))));

	if(qout < RNDOFFERR){ // if there is no output then just use the volume tracer
	  iout = inew;
	} else {
	  intfrac = -hold * powl((1+(qin-qout)/hold) , (-qout/(qin - qout)))/qout - (-hold/qout);
	  iout = iin - (iin - iold) * intfrac;
	}
      } // end time/fraction checks

    } else {                                // If no inflow
      if(qin > RNDOFFERR){ // if inflow but no volume
       inew = iin;
       iout = iin;
      } else {             // if volume but no inflow
       inew = iold;
       iout = inew;
      }

    } // end mixing
  }

}
