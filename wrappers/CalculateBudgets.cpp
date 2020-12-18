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
 * CalculateBudgets.cpp
 *
 *  Created on: Aug 2, 2010
 *      Author: Marco.Maneta
 */

#include "Sativa.h"

int CalculateBudgets(){

  //REAL8 output;

  oBudget->TotalPrecipitation(oAtmosphere->getPrecipitation(), oAtmosphere);
  oBudget->TotalEvaporation(oBasin->getEvaporation(), oBasin);
  oBudget->TotalEvaporationS(oBasin->getEvaporationS_all(), oBasin);
  oBudget->TotalEvaporationC(oBasin->getChanEvap(), oBasin);
  oBudget->TotalEvaporationI(oBasin->getEvaporationI_all(), oBasin);
  oBudget->TotalTranspiration(oBasin->getTranspiration_all(), oBasin);
  oBudget->TotalBedrockLeakage(oBasin->getBedrockLeakage(), oBasin);
  oBudget->TotalOvlndFlow(oBasin->getDailyOvlndOutput(), oBasin);
  oBudget->TotalGrndFlow(oBasin->getDailyGwtrOutput(), oBasin);
  oBudget->TotalStorage(oBasin->getCanopyStorage(),
			oBasin->getSnowWaterEquiv(),
			oBasin->getPondingWater(),
			oBasin->getChanStore(),
			oBasin->getSoilWaterDepthL1(),
			oBasin->getSoilWaterDepthL2(),
			oBasin->getSoilWaterDepthL3(),
			//oBasin->getGravityWater(),
			oBasin->getGrndWater(),
			oBasin->getProotzoneL1(),
			oBasin->getProotzoneL2(),
			oBasin->getProotzoneL3(),
			oBasin);
  oBudget->TotalSrftoChn(oBasin->getFluxSrftoChn(), oBasin);
  oBudget->TotalGWtoChn(oBasin->getFluxGWtoChn(), oBasin);
  oBudget->TotalRecharge(oBasin->getFluxRecharge(), oBasin);
  oBudget->TotalSaturationArea(oBasin->getSatArea(), oBasin);
  
  // ---------------------------------------------------------------------------------------------

  // ################################################################################################
  // Tracking ---------------------------------------------------------------------------------------
  if(oControl->sw_trck){
    // === Deuterium ================================================================================
    if (oControl->sw_2H){
      oBudget->TotalPrecipitation_d2H(oAtmosphere->getPrecipitation(), oAtmosphere->getd2Hprecip(), 
				     oAtmosphere);
      oBudget->TotalEvaporationS_d2H(oBasin->getEvaporationS_all(), oTracking->getd2HevapS_sum(), 
				    oBasin);
      oBudget->TotalEvaporationC_d2H(oBasin->getChanEvap(), oTracking->getd2HevapC_sum(), 
				    oBasin);
      oBudget->TotalEvaporationI_d2H(oBasin->getEvaporationI_all(), oTracking->getd2HevapI_sum(), 
				    oBasin);
      oBudget->TotalTranspiration_d2H(oBasin->getTranspiration_all(), oTracking->getd2HevapT_sum(), 
				     oBasin);
      oBudget->TotalBedrockLeakage_d2H(oBasin->getBedrockLeakage(), oTracking->getd2Hleakage(), oBasin);
      oBudget->TotalOvlndFlow_d2H(oBasin->getDailyOvlndOutput(), oTracking->getd2HOvlndOutput());
      oBudget->TotalGrndFlow_d2H(oBasin->getDailyGwtrOutput(), oTracking->getd2HGwtrOutput());
      //      oBudget->TotalStorage_d2H(oBasin->getCanopyStorage(), oTracking->getd2Hcanopy_sum(),
      //			       oBasin->getSnowWaterEquiv(),  oTracking->getd2Hsnowpack(),
      //			       oBasin->getPondingWater(), oTracking->getd2Hsurface(),
      //			       oBasin->getSoilWaterDepthL1(), oTracking->getd2Hsoil1(),
      //				oBasin->getProotzoneL1(),
      //				oBasin->getSoilWaterDepthL2(), oTracking->getd2Hsoil2(),
      //				oBasin->getProotzoneL2(),
      //				oBasin->getSoilWaterDepthL3(), oTracking->getd2Hsoil3(),
      //				oBasin->getProotzoneL3(),
      //				oBasin->getGrndWater(), oTracking->getd2Hgroundwater(),
      //			       oBasin);//, oControl);

      oBudget->TotalStorage_d2H(oBasin->getCanopyStorage(), oTracking->getd2Hcanopy_sum(),
			       oBasin->getSnowWaterEquiv(),  oTracking->getd2Hsnowpack(),
			       oBasin->getPondingWater(), oTracking->getd2Hsurface(),
			       oBasin->getChanStore(),oTracking->getd2Hchan(),
			       oBasin->getSoilWaterDepthL1(), oTracking->getd2Hsoil1(),
				oBasin->getProotzoneL1(),
				oBasin->getSoilWaterDepthL2(), oTracking->getd2Hsoil2(),
				oBasin->getProotzoneL2(),
				oBasin->getSoilWaterDepthL3(), oTracking->getd2Hsoil3(),
				oBasin->getProotzoneL3(),
				oBasin->getGrndWater(), oTracking->getd2Hsoil3(),
			       oBasin);//, oControl);

      // d2H for Basind2HSummary.txt
      oBudget->InstEvaporation_d2H(oBasin->getEvaporationS_all(), oTracking->getd2HevapS_sum(), 
				   oBasin->getEvaporationI_all(), oTracking->getd2HevapI_sum(), 
				   oBasin->getTranspiration_all(), oTracking->getd2HevapT_sum(), 
				   oBasin);
      oBudget->InstEvaporationS_d2H(oBasin->getEvaporationS_all(), oTracking->getd2HevapS_sum(), 
				     oBasin);
      oBudget->InstEvaporationC_d2H(oBasin->getChanEvap(), oTracking->getd2HevapC_sum(), 
				     oBasin);
      oBudget->InstEvaporationI_d2H(oBasin->getEvaporationI_all(), oTracking->getd2HevapI_sum(), 
				     oBasin);
      oBudget->InstTranspiration_d2H(oBasin->getTranspiration_all(), oTracking->getd2HevapT_sum(), 
				      oBasin);
      oBudget->InstBedrockLeakage_d2H(oBasin->getBedrockLeakage(), oTracking->getd2Hleakage(), 
				       oBasin);

      oBudget->InstOvlndFlow_d2H(oBasin->getDailyOvlndOutput(), oTracking->getd2HOvlndOutput());
      oBudget->InstGrndFlow_d2H(oBasin->getDailyGwtrOutput(), oTracking->getd2HGwtrOutput());

      oBudget->InstOut_d2H(oBasin->getEvaporationS_all(), oTracking->getd2HevapS_sum(), 
			   oBasin->getEvaporationI_all(), oTracking->getd2HevapI_sum(), 
			   oBasin->getTranspiration_all(), oTracking->getd2HevapT_sum(), 
			   oBasin->getBedrockLeakage(), oTracking->getd2Hleakage(), 
			   oBasin->getDailyOvlndOutput(), oTracking->getd2HOvlndOutput(),
			   oBasin->getDailyGwtrOutput(), oTracking->getd2HGwtrOutput(),
			   oBasin);

      oBudget->InstSrftoChn_d2H(oBasin->getFluxSrftoChn(), oTracking->getd2HSrftoChn(), oBasin);
      oBudget->InstGWtoChn_d2H(oBasin->getFluxGWtoChn(), oTracking->getd2HGWtoChn(), oBasin);
      oBudget->InstRecharge_d2H(oBasin->getFluxRecharge(), oTracking->getd2HRecharge(), oBasin);     

    }

    // === Oxygen 18 ====================================================================================
    if (oControl->sw_18O){
      oBudget->TotalPrecipitation_d18O(oAtmosphere->getPrecipitation(), oAtmosphere->getd18Oprecip(), 
				       oAtmosphere);
      oBudget->TotalEvaporationS_d18O(oBasin->getEvaporationS_all(), oTracking->getd18OevapS_sum(), 
					oBasin);
      oBudget->TotalEvaporationC_d18O(oBasin->getChanEvap(), oTracking->getd18OevapC_sum(), 
					oBasin);
      oBudget->TotalEvaporationI_d18O(oBasin->getEvaporationI_all(), oTracking->getd18OevapI_sum(), 
				      oBasin);
      oBudget->TotalTranspiration_d18O(oBasin->getTranspiration_all(), oTracking->getd18OevapT_sum(), 
				       oBasin);
      oBudget->TotalBedrockLeakage_d18O(oBasin->getBedrockLeakage(), oTracking->getd18Oleakage(), 
					oBasin);
      oBudget->TotalOvlndFlow_d18O(oBasin->getDailyOvlndOutput(), oTracking->getd18OOvlndOutput());
      oBudget->TotalGrndFlow_d18O(oBasin->getDailyGwtrOutput(), oTracking->getd18OGwtrOutput());
      //      oBudget->TotalStorage_d18O(oBasin->getCanopyStorage(), oTracking->getd18Ocanopy_sum(),
      //				 oBasin->getSnowWaterEquiv(),  oTracking->getd18Osnowpack(),
      //				 oBasin->getPondingWater(), oTracking->getd18Osurface(),
      //				 oBasin->getSoilWaterDepthL1(), oTracking->getd18Osoil1(),
      //				 oBasin->getProotzoneL1(),
      //				 oBasin->getSoilWaterDepthL2(), oTracking->getd18Osoil2(),
      //				 oBasin->getProotzoneL2(),
      //				 oBasin->getSoilWaterDepthL3(), oTracking->getd18Osoil3(),
      //				 oBasin->getProotzoneL3(),
      //				 oBasin->getGrndWater(), oTracking->getd18Ogroundwater(),
      //				 oBasin);//, oControl);

      oBudget->TotalStorage_d18O(oBasin->getCanopyStorage(), oTracking->getd18Ocanopy_sum(),
				 oBasin->getSnowWaterEquiv(),  oTracking->getd18Osnowpack(),
				 oBasin->getPondingWater(), oTracking->getd18Osurface(),
				 oBasin->getChanStore(),oTracking->getd18Ochan(),
				 oBasin->getSoilWaterDepthL1(), oTracking->getd18Osoil1(),
				 oBasin->getProotzoneL1(),
				 oBasin->getSoilWaterDepthL2(), oTracking->getd18Osoil2(),
				 oBasin->getProotzoneL2(),
				 oBasin->getSoilWaterDepthL3(), oTracking->getd18Osoil3(),
				 oBasin->getProotzoneL3(),
				 oBasin->getGrndWater(), oTracking->getd18Osoil3(),
				 oBasin);//, oControl);

      // d18O for Basind18OSummary.txt
      oBudget->InstEvaporation_d18O(oBasin->getEvaporationS_all(), oTracking->getd18OevapS_sum(), 
				   oBasin->getEvaporationI_all(), oTracking->getd18OevapI_sum(), 
				   oBasin->getTranspiration_all(), oTracking->getd18OevapT_sum(), 
				   oBasin);
      oBudget->InstEvaporationS_d18O(oBasin->getEvaporationS_all(), oTracking->getd18OevapS_sum(), 
				     oBasin);
      oBudget->InstEvaporationC_d18O(oBasin->getChanEvap(), oTracking->getd18OevapC_sum(), 
				     oBasin);
      oBudget->InstEvaporationI_d18O(oBasin->getEvaporationI_all(), oTracking->getd18OevapI_sum(), 
				     oBasin);
      oBudget->InstTranspiration_d18O(oBasin->getTranspiration_all(), oTracking->getd18OevapT_sum(), 
				      oBasin);
      oBudget->InstBedrockLeakage_d18O(oBasin->getBedrockLeakage(), oTracking->getd18Oleakage(), 
				       oBasin);

      oBudget->InstOvlndFlow_d18O(oBasin->getDailyOvlndOutput(), oTracking->getd18OOvlndOutput());
      oBudget->InstGrndFlow_d18O(oBasin->getDailyGwtrOutput(), oTracking->getd18OGwtrOutput());

      oBudget->InstOut_d18O(oBasin->getEvaporationS_all(), oTracking->getd18OevapS_sum(), 
			   oBasin->getEvaporationI_all(), oTracking->getd18OevapI_sum(), 
			   oBasin->getTranspiration_all(), oTracking->getd18OevapT_sum(), 
			   oBasin->getBedrockLeakage(), oTracking->getd18Oleakage(), 
			   oBasin->getDailyOvlndOutput(), oTracking->getd18OOvlndOutput(),
			   oBasin->getDailyGwtrOutput(), oTracking->getd18OGwtrOutput(),
			   oBasin);

      oBudget->InstSrftoChn_d18O(oBasin->getFluxSrftoChn(), oTracking->getd18OSrftoChn(), oBasin);
      oBudget->InstGWtoChn_d18O(oBasin->getFluxGWtoChn(), oTracking->getd18OGWtoChn(), oBasin);
      oBudget->InstRecharge_d18O(oBasin->getFluxRecharge(), oTracking->getd18ORecharge(), oBasin);     

    }

    // === Age =====================================================================================
    if (oControl->sw_Age){
     // Age for mass balance calculation
      oBudget->TotalPrecipitation_Age();
      //oAtmosphere->getPrecipitation(), oAtmosphere);//, oControl);
      //oBudget->precipitation_Age += oBudget
      oBudget->TotalEvaporationS_Age(oBasin->getEvaporationS_all(), oTracking->getAgeevapS_sum(), 
				     oBasin);
      oBudget->TotalEvaporationC_Age(oBasin->getChanEvap(), oTracking->getAgeevapC_sum(), 
				     oBasin);
      oBudget->TotalEvaporationI_Age(oBasin->getEvaporationI_all(), oTracking->getAgeevapI_sum(), 
				     oBasin);
      oBudget->TotalTranspiration_Age(oBasin->getTranspiration_all(), oTracking->getAgeevapT_sum(), 
				      oBasin);
      oBudget->TotalBedrockLeakage_Age(oBasin->getBedrockLeakage(), oTracking->getAgeleakage(), 
				       oBasin);
      oBudget->TotalOvlndFlow_Age(oBasin->getDailyOvlndOutput(), oTracking->getAgeOvlndOutput());
      oBudget->TotalGrndFlow_Age(oBasin->getDailyGwtrOutput(), oTracking->getAgeGwtrOutput());

      //output = oTracking->getAgesnowpack()->maxi();
      //std::cout << output << endl;
      //output = oTracking->getAgesnowpack()->mini();
      //std::cout << output << endl;
      //      oBudget->TotalStorage_Age(oBasin->getCanopyStorage(), oTracking->getAgecanopy_sum(),
      //				oBasin->getSnowWaterEquiv(),  oTracking->getAgesnowpack(),
      //				oBasin->getPondingWater(), oTracking->getAgesurface(),
      //				oBasin->getSoilWaterDepthL1(), oTracking->getAgesoil1(),
      //				oBasin->getProotzoneL1(),
      //				oBasin->getSoilWaterDepthL2(), oTracking->getAgesoil2(),
      //				oBasin->getProotzoneL2(),
      //				oBasin->getSoilWaterDepthL3(), oTracking->getAgesoil3(),
      //				oBasin->getProotzoneL3(),
      //				oBasin->getGrndWater(), oTracking->getAgegroundwater(),
      //				oBasin);//, oControl);

      oBudget->TotalStorage_Age(oBasin->getCanopyStorage(), oTracking->getAgecanopy_sum(),
				oBasin->getSnowWaterEquiv(),  oTracking->getAgesnowpack(),
				oBasin->getPondingWater(), oTracking->getAgesurface(),
				oBasin->getChanStore(),oTracking->getAgechan(),
				oBasin->getSoilWaterDepthL1(), oTracking->getAgesoil1(),
				oBasin->getProotzoneL1(),
				oBasin->getSoilWaterDepthL2(), oTracking->getAgesoil2(),
				oBasin->getProotzoneL2(),
				oBasin->getSoilWaterDepthL3(), oTracking->getAgesoil3(),
				oBasin->getProotzoneL3(),
				oBasin->getGrndWater(), oTracking->getAgesoil3(),
				oBasin);//, oControl);

      // Age for BasinAgeSummary.txt
      oBudget->InstEvaporation_Age(oBasin->getEvaporationS_all(), oTracking->getAgeevapS_sum(), 
				   oBasin->getEvaporationI_all(), oTracking->getAgeevapI_sum(), 
				   oBasin->getTranspiration_all(), oTracking->getAgeevapT_sum(), 
				   oBasin);
      oBudget->InstEvaporationS_Age(oBasin->getEvaporationS_all(), oTracking->getAgeevapS_sum(), 
				     oBasin);
      oBudget->InstEvaporationC_Age(oBasin->getChanEvap(), oTracking->getAgeevapC_sum(), 
				     oBasin);
      oBudget->InstEvaporationI_Age(oBasin->getEvaporationI_all(), oTracking->getAgeevapI_sum(), 
				     oBasin);
      oBudget->InstTranspiration_Age(oBasin->getTranspiration_all(), oTracking->getAgeevapT_sum(), 
				      oBasin);
      oBudget->InstBedrockLeakage_Age(oBasin->getBedrockLeakage(), oTracking->getAgeleakage(), 
				       oBasin);

      oBudget->InstOvlndFlow_Age(oBasin->getDailyOvlndOutput(), oTracking->getAgeOvlndOutput());
      oBudget->InstGrndFlow_Age(oBasin->getDailyGwtrOutput(), oTracking->getAgeGwtrOutput());

      oBudget->InstOut_Age(oBasin->getEvaporationS_all(), oTracking->getAgeevapS_sum(), 
			   oBasin->getEvaporationI_all(), oTracking->getAgeevapI_sum(), 
			   oBasin->getTranspiration_all(), oTracking->getAgeevapT_sum(), 
			   oBasin->getBedrockLeakage(), oTracking->getAgeleakage(), 
			   oBasin->getDailyOvlndOutput(), oTracking->getAgeOvlndOutput(),
			   oBasin->getDailyGwtrOutput(), oTracking->getAgeGwtrOutput(),
			   oBasin);

      oBudget->InstSrftoChn_Age(oBasin->getFluxSrftoChn(), oTracking->getAgeSrftoChn(), oBasin);
      oBudget->InstGWtoChn_Age(oBasin->getFluxGWtoChn(), oTracking->getAgeGWtoChn(), oBasin);
      oBudget->InstRecharge_Age(oBasin->getFluxRecharge(), oTracking->getAgeRecharge(), oBasin);     
      
    }

  } // ---------------------------------------------------------------------------------------------
  // Mass balance check
  oBudget->TrckBalanceError(oControl);
  //oBudget->MassBalanceError();

  return EXIT_SUCCESS;
}
