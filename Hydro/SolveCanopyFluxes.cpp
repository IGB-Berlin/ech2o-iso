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
 * SolveCanopyFluxes.cpp
 *
 *  Created on: Jul 9, 2010
 *      Author: Marco.Maneta
 */

#include"Basin.h"
#include <armadillo>
using namespace arma;

int Basin::SolveCanopyFluxes(Atmosphere &atm, Control &ctrl, Tracking &trck) {
  
  UINT4 r, c;
  REAL8 dt = ctrl.dt;
  
  REAL8 maxTp = 0;
  REAL8 minTp = 0;
  REAL8 snow = 0; //amount of snow reaching the ground ms-1
  REAL8 rain = 0;//amount of rain reaching the ground ms-1
  REAL8 sno_rain_thres = 0; //temperature threshold for snow rain transition, degC

  REAL8 ra; //soil aerodynamic resistance


  REAL8 evap = 0; //evaporation for the tree groves
  REAL8 transp = 0; //transpiratin for the tree groves
  REAL8 netR = 0; //net radiation for the tree groves
  REAL8 leR = 0; //latent heat for the tree groves
  REAL8 haR = 0; //sensible heat for the tree groves
  REAL8 evap_f = 0; //total evaporation for the entire cell
  REAL8 transp_f = 0; //total transpiration for the entire cell

  //canopy storage parameters
  REAL8 D = 0; //canopy trascolation (amount of water that actually reach the ground)
  REAL8 DelCanStor = 0; //Canopy Storage

  //soil parameters
  REAL8 rootdepth;
  REAL8 theta_r1, theta_r2, theta_r3;
  REAL8 psi_ae;
  REAL8 bclambda;

  REAL8 froot1,froot2,froot3; 	// root fractions
  REAL8 d1,d2,d3;  		//soil layer depths

  //aerodynamic resistance parameters
  REAL8 za; //height of wind speed measurements
  REAL8 z0o; // roughness length
  REAL8 zdo; //zero plane displacement
  REAL8 wind; //wind speed
  REAL8 treeheight;

  REAL8 theta,theta2,theta3; 	//soil moisture
  REAL8 theta_wp;
  REAL8 theta_available=0; //water available to roots

  UINT4 nsp;
  REAL8 p; //fraction of species s
  REAL8 veg_p=0; //summed fraction of vegetation

  //unsigned int j;
  UINT4 s;
  int  thre=0;

  // Tracking
  REAL8 pTrp1, pTrp2, pTrp3;
  // Initialize to zero
  _Rn_sum->reset();
  _latheat_sum->reset();
  _sensheat_sum->reset();
  if(ctrl.sw_trck){
    //For tracking
    _FluxCnptoSrf->reset(); // canopy/sky to surface
    _FluxCnptoSnow->reset(); // canopy/sky to snowpack
  }

#pragma omp parallel default(none)					\
  private( s, r,c, p,  treeheight, wind, za, z0o, zdo,			\
	   maxTp, minTp, snow, rain, sno_rain_thres, evap,		\
	   transp, netR, leR,haR,evap_f, transp_f, D, DelCanStor,         \
	   theta, theta2, theta3, theta_wp,theta_available, ra,		\
	   psi_ae, bclambda, rootdepth, froot1, froot2, froot3, d1, d2, d3, \
	   theta_r1, theta_r2, theta_r3, pTrp1, pTrp2, pTrp3, veg_p)		\
  shared(nsp, atm, ctrl, trck, dt, thre, std::cout)

  {
    thre = omp_get_num_threads();
#pragma omp single
    printf("\nnum threads %d: \n", thre);
#pragma omp for nowait
    for (unsigned int j = 0; j < _vSortedGrid.cells.size(); j++) {

      r = _vSortedGrid.cells[j].row;
      c = _vSortedGrid.cells[j].col;

      nsp = fForest->getNumSpecies();
      UINT4 numT = 0;
      
      // *****************************************************************************
      // ---Get the number of tracers to be solved
      // *****************************************************************************
      if(ctrl.sw_trck)
	numT = (ctrl.sw_2H ? 1 : 0) + (ctrl.sw_18O ? 1 : 0) + (ctrl.sw_Age ? 1 : 0);
      colvec tracer(numT);
      tracer = zeros<colvec>(numT,1);
      if(ctrl.sw_trck){
	if(numT == 1)
	  tracer(0) = (ctrl.sw_2H ? 0 : (ctrl.sw_18O ? 1 : 2));
	if(numT == 2){
	  tracer(0) = (ctrl.sw_2H ? 0 : (ctrl.sw_18O ? 1 : 2));
	  tracer(1) = (ctrl.sw_2H ? (ctrl.sw_18O ? 1 : 2) : 2);
	} else {
	  tracer(0) = 0;
	  tracer(1) = 1;
	  tracer(2) = 2;	  
        }
      }

      // *****************************************************************************
      // ---Get the separation of precipitation into rain and snow--------------------
      // *****************************************************************************
      maxTp = atm.getMaxTemperature()->matrix[r][c];
      minTp = atm.getMinTemperature()->matrix[r][c];
      sno_rain_thres = atm.getRainSnowTempThreshold();
      REAL8 SRfrac   = 0; 				//snow/rain fraction (1 is all snow - 0 is no snow)
      REAL8 Throughfall = 0;

      // *****************************************************************************
      // ---Soil and vegetation initialization----------------------------------------
      // *****************************************************************************
      treeheight = 0;
      evap_f = 0;
      transp_f = 0;

      theta_r1 = _theta_rL1->matrix[r][c];
      theta_r2 = _theta_rL2->matrix[r][c];
      theta_r3 = _theta_rL3->matrix[r][c];

      psi_ae = _psi_ae->matrix[r][c];
      bclambda = _BClambda->matrix[r][c];

      d1 = _depth_layer1->matrix[r][c];
      d2 = _depth_layer2->matrix[r][c];
      d3 = _soildepth->matrix[r][c] - d1 - d2;

      // *****************************************************************************
      // ---Tracking initialization---------------------------------------------------
      // *****************************************************************************
      REAL8 Pin = 0; 							//tracer precip in
      mat evapTV(nsp,3), evapTV_f(nsp,3),evapTVp(nsp,3);                //Tra vapour weighting
      mat evapT(nsp,3), evapT_f(nsp,3), evapTp(nsp,3);                  //RU weighting
      mat evapI(nsp,3), evapI_f(nsp,3), evapIp(nsp,3);                  //canopy E weighting      
      mat canopy(nsp,3),canopy_f(nsp,3),canopyp(nsp,3);                 //canopy water weighting
      mat thro(nsp,3),  thro_f(nsp,3),  throp(nsp,3);                   //throughfall weighting

      evapTV = evapTV_f = evapTVp = zeros<mat>(nsp,3);
      evapT  = evapT_f  = evapTp  = zeros<mat>(nsp,3);
      evapI  = evapI_f  = evapIp  = zeros<mat>(nsp,3);
      canopy = canopy_f = canopyp = zeros<mat>(nsp,3);
      thro   = thro_f   = throp   = zeros<mat>(nsp,3);
      // *****************************************************************************
      // --- LOOP ON SPECIES ---------------------------------------------------------
      // *****************************************************************************
      for (s = 0; s < nsp; s++) {
     
	p = fForest->getPropSpecies(s, r, c);
	if (p == 0)
	  continue; //if no species j present, continue
     
	D = DelCanStor = 0; 		// initialize precip and dStorage
	evap = transp = 0; 		// initialize evap and transp
        netR = leR = haR = 0; 		// initialize energy balance terms

	SRfrac = min<REAL8>(1.0,max<REAL8>(0.0, (sno_rain_thres - minTp) / (maxTp - minTp)));
	
	if (s == nsp - 1) { //if this is bare ground set D to precip and skip the tree stuff
	  D = atm.getPrecipitation()->matrix[r][c];
	  thro(s,0) = 1;//ctrl.sw_2H  ? atm.getd2Hprecip()->matrix[r][c]  : 0.0; 	//throughfall d2H
	  thro(s,1) = 1;//ctrl.sw_18O ? atm.getd18Oprecip()->matrix[r][c] : 0.0; 	//throughfall d18O
	  thro(s,2) = 1;//0.0; 							//throughfall age
	  veg_p += 0;
	} else {
	  // **************************************************************************************
          // --- Solve Canopy Energy Balance-------------------------------------------------------
	  // **************************************************************************************      
	  wind = atm.getWindSpeed()->matrix[r][c];
	  veg_p += p;
	  treeheight = max<REAL8>(0.01, fForest->getTreeHeight(s, r, c));
	  za = treeheight + 2;  //*TODO: Tentative rel between forest height and wind velocity profile parameters*/
	  z0o = powl(treeheight, 1.19) * 0.057544; 
	  zdo = powl(treeheight, 0.98) * 0.707946; 

	  theta_wp = fForest->getWiltPoint(s);              
	  theta = _soilmoist1->matrix[r][c]; //soil moisture at time t
	  theta2 = _soilmoist2->matrix[r][c];
	  theta3 = _soilmoist3->matrix[r][c];
	  froot1 = fForest->getRootFrac1(s)->matrix[r][c];
	  froot2 = fForest->getRootFrac2(s)->matrix[r][c];
	  froot3 = 1 - froot1 - froot2;

	  theta_available = std::max(0.0,(theta -theta_wp)) * froot1 + 
	    		    std::max(0.0,(theta2-theta_wp)) * froot2 + 
			    std::max(0.0,(theta3-theta_wp)) * froot3;

	  rootdepth = _Zroot95->matrix[r][c];  	//root depth is the depth of layers that contain 95% of roots

	  ra = CalcAerodynResist(wind, za, 0, 0, z0o, zdo, treeheight,
				 fForest->getLAISpecies(s, r, c),
				 getCanopyTemp(s)->matrix[r][c],
				 atm.getTemperature()->matrix[r][c], ctrl.toggle_ra,
				 false);
       
	  fForest->CanopyInterception(atm, ctrl, DelCanStor, D, s, r, c); //calculates canopy interception and trasc

	  fForest->SolveCanopyEnergyBalance(*this, atm, ctrl, theta_r1, theta_r2, theta_r3,
					    rootdepth, psi_ae, bclambda, ra, DelCanStor, evap, 
					    transp, netR, leR, haR, s, r, c);

	  _CanopyStorage->matrix[r][c] += DelCanStor * p;
          _CanopyStorage->matrix[r][c] = (_CanopyStorage->matrix[r][c]<RNDOFFERR) ? 0.0 : _CanopyStorage->matrix[r][c];

	  _Rn_sum->matrix[r][c] += netR * p ;                         //update net radiation
	  _latheat_sum->matrix[r][c] += leR * p;                      //update LE (int + tra)
	  _sensheat_sum->matrix[r][c] += haR * p;                     //update sensible heat
	  evap_f += evap * p; 					      //evaporation at t=t+1 [m/s]
	  transp_f += transp * p; 				      //transpiration at t=t+1 [m/s]
       
	  pTrp1 = (std::max(0.0,(theta-theta_wp))*froot1) / theta_available;
	  pTrp2 = (std::max(0.0,(theta2-theta_wp))*froot2) / theta_available;
          pTrp3 = (std::max(0.0,(theta3-theta_wp))*froot3) / theta_available;

	  theta  -= transp * p * dt * pTrp1 /d1 ;  		     //soil moisture at t=t+1
	  theta2 -= transp * p * dt * pTrp2 /d2 ; 		     //soil moisture at t=t+1
	  theta3 -= transp * p * dt * pTrp3 /d3 ; 		     //soil moisture at t=t+1

	  fForest->setTranspirationFlow(s,r,c, ((pTrp1 + pTrp2 + pTrp3) * transp * p  * _dx * _dx) ); // m3/s

	  // *********************************************************************************************
	  // ---TRACKING----------------------------------------------------------------------------------
	  // *********************************************************************************************
	  if(ctrl.sw_trck){
	    for(UINT4 n = 0; n<numT ; n++ ) { //loop through all tracers
   	      canopy(s, n) = (tracer(n) == 0) ? fForest->getd2Hcanopy(s)->matrix[r][c] :
                     		((tracer(n) == 1) ? fForest->getd18Ocanopy(s)->matrix[r][c] : 
					   fForest->getAgecanopy(s)->matrix[r][c]); // initialize the canopy
              evapT(s, n)  = (tracer(n) == 0) ? fForest->getd2HevapT(s)->matrix[r][c] :
		                ((tracer(n) == 1) ? fForest->getd18OevapT(s)->matrix[r][c] : 
					   fForest->getAgeevapT(s)->matrix[r][c]);
	      Pin            = (tracer(n) == 0) ? atm.getd2Hprecip()->matrix[r][c] :
				          ((n == 1) ? atm.getd18Oprecip()->matrix[r][c] : 0);

 	      trck.MixingV_canopy(atm,*this,ctrl, SRfrac, fForest->getIntercWater(s,r,c),evap,
				  canopy(s,tracer(n)),evapI(s,tracer(n)),thro(s,tracer(n)),
				  D,Pin,p,r,c,s,tracer(n) );

	      trck.MixingV_sapflow(atm,*this,ctrl,pTrp1,pTrp2,pTrp3, transp,
				   evapT(s,tracer(n)),evapTV(s,tracer(n)),p,r,c,s,tracer(n));
	      
	      canopy_f(s,tracer(n)) = canopy(s,tracer(n)) * p * fForest->getIntercWater(s,r,c);
	    } // close tracer loop
	    canopyp.row(s)  	= canopy.row(s) *  p;          	//Canopy Storage
	    throp.row(s)        = thro.row(s) * (D * p * dt);   //throughfall
	    evapTp.row(s) 	= evapT.row(s) * p; 		//transpiration
	    evapT_f.row(s) 	= evapT.row(s) * (p * transp); 	//transpiration
	    evapTVp.row(s)	= evapTV.row(s) * p; 		//Transp vapour
	    evapTV_f.row(s)     = evapTV.row(s) * (p * transp); //Transp vapour
	    evapIp.row(s) 	= evapI.row(s) * p; 		//Int E vapour
	    evapI_f.row(s) 	= evapI.row(s) * (p * evap);	//Int E vapour
          } // close tracking

	  // *********************************************************************************************
	  // ---Update soil moisture objects--------------------------------------------------------------
	  // *********************************************************************************************
	  _soilmoist1->matrix[r][c] = theta;
	  _soilmoist2->matrix[r][c] = theta2;
	  _soilmoist3->matrix[r][c] = theta3;
	      
	} // end if bare / veg
     
	snow = D * SRfrac;
	rain = D - snow;
	Throughfall += (snow + rain)*p*dt;

        // *********************************************************************************************     
	// ---Water tracking----------------------------------------------------------------------------
        // *********************************************************************************************     
	if(ctrl.sw_trck){
	  // Used later in SolveSurfaceFluxes
	  _FluxCnptoSrf->matrix[r][c] += rain * p * dt;
	  _FluxCnptoSnow->matrix[r][c] += snow * p * dt;
	  trck.MixingV_through(atm, *this, ctrl, rain, p,s, r, c);  // Mix in surface pool from rain throughfall
	  if(ctrl.sw_2H){
	    fForest->setd2Hcanopy(s,r,c,canopy(s,0));
	    fForest->setd2Hthroughfall(s,r,c,thro(s,0));
	    fForest->setd2HevapI(s,r,c,evapI(s,0));
	    fForest->setd2HevapT(s,r,c,evapT(s,0));
	  }
	  if(ctrl.sw_18O){
	    fForest->setd18Ocanopy(s,r,c,canopy(s,1));
	    fForest->setd18Othroughfall(s,r,c,thro(s,1));
	    fForest->setd18OevapI(s,r,c,evapI(s,1));
	    fForest->setd18OevapT(s,r,c,evapT(s,1));
	  }
	  if(ctrl.sw_Age){
	    fForest->setAgecanopy(s,r,c,canopy(s,2));
	    fForest->setAgethroughfall(s,r,c,thro(s,2));
	    fForest->setAgeevapI(s,r,c,evapI(s,2));
	    fForest->setAgeevapT(s,r,c,evapT(s,2));
	  }
	}
     
	_snow->matrix[r][c] +=  snow * dt * p;
	_incident_water_depth->matrix[r][c] += rain *dt *p;
      } //end for over species
   
      _Evaporation->matrix[r][c] = evap_f + transp_f; //total evaporation for the entire cell
      _Transpiration_all->matrix[r][c]  = transp_f ;
      _EvaporationI_all->matrix[r][c] = evap_f ;

      // *********************************************************************************************
      // ---TRACKING----------------------------------------------------------------------------------
      // summed evapT 2H, 18O and Age weighted by flux magnitude and cover ONLY over vegetated fraction
      // *********************************************************************************************
      if(ctrl.sw_trck){
        arma::rowvec thro_all   = arma::sum(throp,0)   / (Throughfall); // set throughfall
	arma::rowvec evapT_all  = arma::sum(evapT_f,0) / transp_f;
	arma::rowvec evapT_all2 = arma::sum(evapTp,0)  / veg_p;
        arma::rowvec evapI_all  = arma::sum(evapI_f,0) / evap_f;
        arma::rowvec evapI_all2 = arma::sum(evapIp,0)  / veg_p;
	arma::rowvec canopy_all = arma::sum(canopy_f,0)/ _CanopyStorage->matrix[r][c];
	arma::rowvec canopy_all2= arma::sum(canopyp,0) / veg_p;
	UINT4 n;
	
	if(ctrl.sw_2H){
	  n = 0;
	  trck.setd2Hthroughfall_sum(r,c,Throughfall > RNDOFFERR ? thro_all[n] : atm.getd2Hprecip()->matrix[r][c]);
	  trck.setd2HevapT_sum(r,c,transp_f > RNDOFFERR ? evapT_all[n] : evapT_all2[n]);
	  trck.setd2HevapI_sum(r,c,evap_f > RNDOFFERR ? evapI_all[n] :evapI_all2[n]);
	  trck.setd2Hcanopy_sum(r,c,_CanopyStorage->matrix[r][c] > RNDOFFERR ? canopy_all[n] : canopy_all2[n]);
	}
	if(ctrl.sw_18O){
	  n = 1;
	  trck.setd18Othroughfall_sum(r,c,Throughfall > RNDOFFERR ? thro_all[n] : atm.getd18Oprecip()->matrix[r][c]);
	  trck.setd18OevapT_sum(r,c,transp_f > RNDOFFERR ? evapT_all[n] : evapT_all2[n]);
	  trck.setd18OevapI_sum(r,c,evap_f > RNDOFFERR ? evapI_all[n] :evapI_all2[n]);
	  trck.setd18Ocanopy_sum(r,c,_CanopyStorage->matrix[r][c] > RNDOFFERR ? canopy_all[n] : canopy_all2[n]);
	}
	if(ctrl.sw_Age){
	  n = 2;
	  trck.setAgethroughfall_sum(r,c,Throughfall > RNDOFFERR ? thro_all[n] : 0);
	  trck.setAgeevapT_sum(r,c,transp_f > RNDOFFERR ? evapT_all[n] : evapT_all2[n]);
	  trck.setAgeevapI_sum(r,c,evap_f > RNDOFFERR ? evapI_all[n] :evapI_all2[n]);
	  trck.setAgecanopy_sum(r,c,_CanopyStorage->matrix[r][c] > RNDOFFERR ? canopy_all[n] : canopy_all2[n]);
	}
      } // End tracking
      veg_p = 0;
   
    }//end for
  }//end omp parallel
   
  return EXIT_SUCCESS;
}
