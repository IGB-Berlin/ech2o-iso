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
 * sativa.cpp
 *
 *  Created on: Feb 19, 2010
 *      Author: Marco Maneta
 */

/*
 * main.cpp
 *
 *  Created on: May 21, 2009
 *      Author: Marco Maneta
 */
#include <time.h>
#include "Sativa.h"

float report_time = 0; //resets to zero when Report_interval time interval passes
float reportMap_time = 0; //resets to zero when ReportMap_interval time interval passes
float advance_climate = 0; // resets to zero when Clim_input_tstep passess

time_t start, theend;
int main(int argc, char* argv[]) {
# ifdef _OPENMP
  printf("Compiled by an OpenMP-compliant implementation.\n");
# endif
  try {
    time(&start);
    Splash(argc, argv);
    CreateWorld(argv);

    if(oControl->starttime > 0){
      cout << "Start time later than climate start time" << endl;
      cout << "   -- Spin up climate maps to start time" << endl;
      // need to initialize nc files here if start time is not 1
      /* if(oControl->sw_netcdf){ //commented out for now
	oReport->CreatOutputNC(oControl->path_ResultsFolder, "W");
	oReport->CreatOutputNC(oControl->path_ResultsFolder, "VG");
	if(oControl->sw_trck){
	  if(oControl->sw_2H)
	    oReport->CreatOutputNC(oControl->path_ResultsFolder, "TD");
	  if(oControl->sw_18O)
	    oReport->CreatOutputNC(oControl->path_ResultsFolder, "TO");
	  if(oControl->sw_Age)
	    oReport->CreatOutputNC(oControl->path_ResultsFolder, "TA");
	}
      }*/ //end net cdf initialization
      // Run loop until correct starting time
      while (oControl->current_t_step <= oControl->starttime) {
	oControl->AdvanceTimeStep();
	advance_climate += oControl->dt;
	if (advance_climate >= oControl->BC_dt) {
	  oAtmosphere->AdvanceClimateMaps(*oControl);

	  if(oControl->toggle_veg_dyn==2){
	    oBasin->AdvanceLAIMaps();
	  }
	  
	  if(oControl->sw_BC){
	    oBasin->AdvanceBCMaps(*oAtmosphere,*oControl);
	    if(oControl->sw_trck){
	      if(oControl->sw_2H)
		 oTracking->AdvanceBCMaps_Iso(*oAtmosphere,*oControl,0);
	      if(oControl->sw_18O)
		 oTracking->AdvanceBCMaps_Iso(*oAtmosphere,*oControl,1);
	      if(oControl->sw_Age)
		 oTracking->AdvanceBCMaps_Iso(*oAtmosphere,*oControl,2);
	    }
	  }
	
	  advance_climate = 0;
	}	
      } // end while loop
    } // end skipping to start time

    while (oControl->current_t_step <= oControl->endtime) {

      SolveTimeStep();

      CalculateBudgets();

      Report2Screen();

      // Report time series
      report_time += oControl->dt;
      if (report_time >= oControl->report_times) { //if report time overdue
	Report2Ts(); //report results
	report_time = 0; //reset the counter
      }

      // Report maps (only from a certain time step, e.g. to avoid spinup map reporting)
      if (oControl->current_t_step >= oControl->reportMap_start) { 
	reportMap_time += oControl->dt;
	if (reportMap_time >= oControl->reportMap_times) { //if report time overdue
	  //if(oControl->sw_netcdf){ //commented out for now
	  //  Report2nc(); //report results in nc
	  //} else {
	    Report2Maps(); //report results
	  //}
	  reportMap_time = 0; //reset the counter
	}
      }

      cout << "\nEnd time step " << oControl->current_ts_count;
      cout << "\nSimulation time " << oControl->current_t_step
	   << " seconds (" << oControl->current_t_step / 86400
	   << " days)\n\n";

      oControl->AdvanceTimeStep();

      advance_climate += oControl->dt;
      if (advance_climate >= oControl->BC_dt) {
	oAtmosphere->AdvanceClimateMaps(*oControl);

	if(oControl->toggle_veg_dyn==2){
	  oBasin->AdvanceLAIMaps();
	}
	
	if(oControl->sw_BC){
	  oBasin->AdvanceBCMaps(*oAtmosphere,*oControl);
	  if(oControl->sw_trck){
	    if(oControl->sw_2H)
	      oTracking->AdvanceBCMaps_Iso(*oAtmosphere,*oControl,0);
	    if(oControl->sw_18O)
	      oTracking->AdvanceBCMaps_Iso(*oAtmosphere,*oControl,1);
	    if(oControl->sw_Age)
	      oTracking->AdvanceBCMaps_Iso(*oAtmosphere,*oControl,2);
	  }
	}
	
	advance_climate = 0;
      }

    }

  } catch (...) {
    cerr
      << "Something bad happened that I cannot really handle until I have a better exception management"
      << endl;
    CrunchWorld();
    return 0;
  }
  CrunchWorld();
  time(&theend);
  int tot_sec = difftime(theend, start);

  int dd = tot_sec / 86400;
  tot_sec = tot_sec % 86400;
  int hh = tot_sec / 3600;
  tot_sec = tot_sec % 3600;
  int mm = tot_sec / 60;
  tot_sec = tot_sec % 60;
  int ss = tot_sec;
  printf("\nTotal run time elapsed:  %i (days) %02i:%02i:%02i (hh:mm:ss)\n", dd,
	 hh, mm, ss);

  return 0;

}
