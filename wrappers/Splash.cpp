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
 * Splash.cpp
 *
 *  Created on: Sep 2, 2010
 *      Author: Marco.Maneta, Sylvain Kuppel
 */

#include "Sativa.h"
#include <string.h>

void message(){
	cout << "USAGE: ech2o <config.ini>" << endl;
	cout << "\twhere config.ini is a valid ech2o configuration file" << endl;
	cout << "OR: ech2o -g <output_file>" << endl;
	cout << "\tto write two configuration file templates:" << endl;
	cout << "\t<output_file> + its tracking-specific counterpart." << endl;
}

void Splash(int argc, char* argv[]){

  cout << "|______________________________________________________________________________|" << endl;
  cout << "                                                                              ||" << endl;
  cout << "  ------------                  --       --       -------------               ||" << endl;
  cout << "  888888888888                  88       88       8888888888888               ||" << endl;
  cout << "  88                            88       88       88         88               ||" << endl;
  cout << "  88            -------------   88       88       88         88               ||" << endl;
  cout << "  888888888     8888888888888   88888888888       88         88               ||" << endl;
  cout << "  888888888     88              88888888888       88         88               ||" << endl;
  cout << "  88            88              88       88 8888  88         88               ||" << endl;
  cout << "  88            88              88       88    88 88         88               ||" << endl;
  cout << "  888888888888  8888888888888   88       88   88  8888888888888  === ==== === ||" << endl;
  cout << "  ------------  -------------   --      -- 888    --------------- 8  88   8 8 ||" << endl;
  cout << "                                           8888888                8    88 8 8 ||" << endl;
  cout << "                                                                 *** **** *** ||" << endl;
  cout << "   ECOHYDROLOGICAL MODEL, isotopes and age tracking (EcH2O-iso)               |" << endl;
  cout << "   University of Montana (USA), University of Aberdeen (UK), IGB-Berlin (DE)  |" << endl;
  cout << "   "<< VERSION << "                                                         " << endl;
  cout << " _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ " <<endl << endl;


     if (argc==3 && strcmp(argv[1],"-g")==0 ){
    	GenerateConfigTemplate(argv[2]);
	cout << "Writing configuration file template in " << argv[2] << endl;
	// For now the name for configTrck.ini is fixed
    	GenerateConfigTrckTemplate("configTrck.ini");
	cout << "Writing tracking configuration file template in configTrck.ini" << endl;
	exit(EXIT_SUCCESS);
     }

    if (argc!=2 || strcmp(argv[1],"-g") == 0){
		message();
		exit(EXIT_SUCCESS);
	}




}
