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
 * Report.h
 *
 *  Created on: Sep 21, 2010
 *      Author: Marco.Maneta
 */

#ifndef REPORT_H_
#define REPORT_H_

#include "SortGrid.h"
#include "InitConf.h"

//Map write defines
#define WriteMapSeries(_var, _s, _tscount)  \
     _var->grid2PCRMap( (oControl->path_ResultsFolder + ParseString(_s, _tscount)).c_str(), CR_REAL4, VS_SCALAR);

#define WriteMap(_var, _s) \
  _var->grid2PCRMap( (oControl->path_ResultsFolder +_s+".map").c_str(), CR_REAL4, VS_SCALAR);


struct Report{

	vectCells mask;

	Report(){};
	Report(Control &ctrl);

	int ReportTimeSeries(const grid *input, string filename, float timestep);
	int ReportVectCells(const vectCells *input, string filename, float timestep);
	int RenameFile(string oldname);
	int UpdateOutputNC(const grid *input, string varname, string outTP);
	int CreatOutputNC(string filepath, string outTP);
};
#endif /* REPORT_H_ */
