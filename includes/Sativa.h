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
 * sativa.h
 *
 *  Created on: Jul 30, 2010
 *      Author: Marco Maneta, Sylvain Kuppel
 */

#ifndef SATIVA_H_
#define SATIVA_H_

#define VERSION "V 3.0"

#include "Basin.h"
#include "Atmosphere.h"
#include "Budget.h"
#include "Report.h"
#include "Tracking.h"

extern Basin *oBasin;
extern Atmosphere *oAtmosphere;
extern Control *oControl;
extern Budget *oBudget;
extern Report *oReport;
extern Tracking *oTracking;

extern ofstream ofSummary;
extern ofstream ofd2HSummary;
extern ofstream ofd18OSummary;
extern ofstream ofAgeSummary;

void Splash(int argc, char* argv[]);
int CreateWorld(char* argv[]);
int SolveTimeStep();
int CalculateBudgets();
int Report2Screen();
int Report2Maps();
int Report2Ts();
int CrunchWorld();

void GenerateConfigTemplate(const char *fn);
void GenerateConfigTrckTemplate(const char *fn);
#endif /* SATIVA_H_ */
