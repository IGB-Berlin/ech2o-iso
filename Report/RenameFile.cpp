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
 * RenameFile.cpp
 *
 *  Created on: Sep 23, 2010
 *      Author: Marco.Maneta
 */

#include<errno.h>
#include "Report.h"

int Report::RenameFile(string oldname){


	remove(oldname.c_str());
	if(errno != 0 && errno != 2)
	{
		cout << errno << " " << oldname;
		perror("error");
		throw (EXIT_FAILURE);
	}
/*	rename(oldname.c_str(), (oldname + ".old").c_str());
	if(errno != 0 && errno != 2)
	{
		cout << errno << " " << oldname;
		perror("Error creating backup .old file");
		throw (EXIT_FAILURE);
	}*/

	return EXIT_SUCCESS;
}
